"""Tests for the blabber module."""

import argparse
import gzip
import pathlib
import re
from itertools import product

import pytest

from biomate.blabber.blabber import (
    SURFACE_COUNT,
    SWATH_COUNT,
    TILE_COUNT,
    assign_tiles_to_samples,
    main,
    parse_sequence_mask,
    partition_tiles_by_lane,
    tiles_per_lane,
    validate_args,
)


# ---------------------------------------------------------------------------
# parse_sequence_mask
# ---------------------------------------------------------------------------


class TestParseSequenceMask:
    """Tests for parse_sequence_mask."""

    def test_r1_only(self):
        """A single-segment Y mask sets R1 and leaves all other fields at 0."""
        result = parse_sequence_mask("Y151")
        assert result["R1"] == 151
        assert result["R2"] == 0

    def test_r1_and_r2(self):
        """A four-segment mask sets both R1 and R2 correctly."""
        # A 2-segment mask "Y151;Y151" maps both to index-1 (R1) due to i_substring//2+1.
        # R2 requires 4 segments; use a full 4-part mask to verify R2.
        result = parse_sequence_mask("Y151;I0;I0;Y151")
        assert result["R1"] == 151
        assert result["R2"] == 151
        assert result["I1"] == 0
        assert result["I2"] == 0

    def test_dual_index(self):
        """A four-segment mask with two I segments sets I1, I2, R1 and R2."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert result["R1"] == 151
        assert result["I1"] == 10
        assert result["I2"] == 10
        assert result["R2"] == 151

    def test_single_index_with_trailing_trim(self):
        """An I8N2 segment sets I1=8 and the trailing-trim I1A=2."""
        # I8N2 → I1=8, I1A=2 (the N after setting I1 is a suffix)
        result = parse_sequence_mask("Y151;I8N2;Y151")
        assert result["R1"] == 151
        assert result["I1"] == 8
        assert result["I1A"] == 2
        assert result["R2"] == 151

    def test_prefix_and_suffix_around_index(self):
        """N before I sets I1B; N after I sets I1A."""
        # N1I8N2 → I1B=1, I1=8, I1A=2
        result = parse_sequence_mask("Y151;N1I8N2;Y151")
        assert result["I1B"] == 1
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_umi_read1(self):
        """A UMI prefix (U8) before a read segment is parsed into U1."""
        # U8Y143 → U1=8, R1=143
        result = parse_sequence_mask("U8Y143;I10;I10;Y151")
        assert result["U1"] == 8
        assert result["R1"] == 143
        assert result["I1"] == 10
        assert result["R2"] == 151

    def test_too_many_sections_raises(self):
        """A mask with more than 4 semicolon-separated sections raises ValueError."""
        with pytest.raises(ValueError, match="OverrideCycles"):
            parse_sequence_mask("Y151;I10;I10;Y151;Y50")

    def test_with_index1_matching_length(self):
        """When index1 length matches the I1 mask value, index1 replaces the int."""
        result = parse_sequence_mask("Y151;I10;I10;Y151", index1="GAACTGAGCG")
        assert result["I1"] == "GAACTGAGCG"

    def test_with_index1_length_mismatch_raises(self):
        """A mismatch between index1 length and the I1 mask value raises ValueError."""
        with pytest.raises(ValueError, match="Index 1"):
            parse_sequence_mask("Y151;I10;I10;Y151", index1="ACGT")  # len 4 ≠ 10

    def test_with_index2_matching_length(self):
        """When index2 length matches the I2 mask value, index2 replaces the int."""
        result = parse_sequence_mask("Y151;I10;I10;Y151", index2="CGCTCCACGA")
        assert result["I2"] == "CGCTCCACGA"

    def test_with_index2_length_mismatch_raises(self):
        """A mismatch between index2 length and the I2 mask value raises ValueError."""
        with pytest.raises(ValueError, match="Index 2"):
            parse_sequence_mask("Y151;I10;I10;Y151", index2="ACGT")  # len 4 ≠ 10

    def test_all_zero_fields_default(self):
        """Unused mask fields default to 0 in the returned dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        # Unused fields must default to 0
        assert result["U1"] == 0
        assert result["U2"] == 0
        assert result["R1B"] == 0
        assert result["R1A"] == 0


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    """Tests for validate_args."""

    def _base_args(self, **overrides):
        """Build a valid Namespace, overridable per test."""
        args = argparse.Namespace(
            alphabet="ACGT",
            seq_length=50,
            seq_number=100,
            seq_mask=None,
            index1=None,
            index2=None,
            format="text",
            sample_sheet=None,
            output=None,
        )
        for k, v in overrides.items():
            setattr(args, k, v)
        return args

    def test_valid_args_returned(self):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._base_args()
        assert validate_args(args) is args

    def test_non_alpha_alphabet_raises(self):
        """A non-alphabetic character in --alphabet raises ArgumentTypeError."""
        args = self._base_args(alphabet="ACG1")
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]etters|[Aa]lpha"):
            validate_args(args)

    def test_seq_length_zero_raises(self):
        """--seq-length 0 raises ArgumentTypeError."""
        args = self._base_args(seq_length=0)
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]ength"):
            validate_args(args)

    def test_seq_length_negative_raises(self):
        """A negative --seq-length raises ArgumentTypeError."""
        args = self._base_args(seq_length=-5)
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]ength"):
            validate_args(args)

    def test_seq_number_zero_raises(self):
        """--seq-number 0 raises ArgumentTypeError."""
        args = self._base_args(seq_number=0)
        with pytest.raises(argparse.ArgumentTypeError, match="[Nn]umber"):
            validate_args(args)

    def test_invalid_index1_non_alpha_raises(self):
        """A non-alphabetic character in --index1 raises ArgumentTypeError."""
        args = self._base_args(index1="ACGT1")
        with pytest.raises(argparse.ArgumentTypeError, match="Index 1"):
            validate_args(args)

    def test_invalid_index2_non_alpha_raises(self):
        """A non-alphabetic character in --index2 raises ArgumentTypeError."""
        args = self._base_args(index2="123")
        with pytest.raises(argparse.ArgumentTypeError, match="Index 2"):
            validate_args(args)

    def test_invalid_seq_mask_raises(self):
        """An invalid --seq-mask raises ValueError."""
        # 5 sections is invalid
        args = self._base_args(seq_mask="Y151;I10;I10;Y151;Y50")
        with pytest.raises(ValueError, match="OverrideCycles"):
            validate_args(args)

    def test_fastq_format_with_sample_sheet_requires_output(self, tmp_path):
        """FASTQ format with --sample-sheet but no --output raises ArgumentTypeError."""
        ss = tmp_path / "ss.csv"
        ss.write_text("")
        args = self._base_args(format="fastq", sample_sheet=ss, output=None)
        with pytest.raises(argparse.ArgumentTypeError, match="[Oo]utput"):
            validate_args(args)


# ---------------------------------------------------------------------------
# partition_tiles_by_lane / assign_tiles_to_samples
# ---------------------------------------------------------------------------


def all_tile_ids() -> set:
    """Return the set of all tile identifiers used by blabber."""
    return {
        "".join(x)
        for x in product(
            [str(i + 1) for i in range(SURFACE_COUNT)],
            [str(i + 1) for i in range(SWATH_COUNT)],
            [f"{i + 1:02}" for i in range(TILE_COUNT)],
        )
    }


class TestPartitionTilesByLane:
    """Tests for partition_tiles_by_lane."""

    def test_equal_partition_when_divisible(self):
        """128 tiles across 4 lanes gives exactly 32 tiles per lane."""
        tiles = sorted(all_tile_ids())
        result = partition_tiles_by_lane(tiles, ["L001", "L002", "L003", "L004"])
        assert set(result) == {"L001", "L002", "L003", "L004"}
        assert all(len(pool) == 32 for pool in result.values())

    def test_balanced_partition_when_not_divisible(self):
        """128 tiles across 3 lanes gives pools of 43, 43 and 42 tiles."""
        tiles = sorted(all_tile_ids())
        result = partition_tiles_by_lane(tiles, ["L001", "L002", "L003"])
        assert sorted(len(pool) for pool in result.values()) == [42, 43, 43]

    def test_partition_is_disjoint_and_covers_all_tiles(self):
        """The pools are pairwise disjoint and their union is the full tile set."""
        tiles = sorted(all_tile_ids())
        result = partition_tiles_by_lane(tiles, ["L001", "L002", "L003", "L004"])
        seen = []
        for pool in result.values():
            seen.extend(pool)
        assert len(seen) == len(set(seen)) == len(tiles)
        assert set(seen) == set(tiles)

    def test_single_lane_gets_all_tiles(self):
        """A single lane receives the complete tile set."""
        tiles = sorted(all_tile_ids())
        result = partition_tiles_by_lane(tiles, ["L001"])
        assert result["L001"] == tiles

    def test_more_lanes_than_tiles_raises(self):
        """More lanes than tiles raises ValueError."""
        with pytest.raises(ValueError, match="Cannot partition"):
            partition_tiles_by_lane(["1101"], ["L001", "L002"])


class TestAssignTilesToSamples:
    """Tests for assign_tiles_to_samples."""

    def test_fewer_samples_than_tiles_covers_all_tiles(self):
        """With 4 samples and 16 tiles, every sample gets 4 tiles and all are used."""
        pool = [f"11{i:02d}" for i in range(1, 17)]
        result = assign_tiles_to_samples(pool, 4)
        assert len(result) == 4
        assert all(len(tiles) == 4 for tiles in result)
        assert sorted(t for tiles in result for t in tiles) == pool

    def test_more_samples_than_tiles_covers_all_tiles(self):
        """With 6 samples and 4 tiles, each sample gets 1 tile and all 4 are used."""
        pool = ["1101", "1102", "1103", "1104"]
        result = assign_tiles_to_samples(pool, 6)
        assert len(result) == 6
        assert all(len(tiles) == 1 for tiles in result)
        assert {t for tiles in result for t in tiles} == set(pool)

    def test_equal_samples_and_tiles(self):
        """With 5 samples and 5 tiles, each sample gets exactly its own tile."""
        pool = ["1101", "1102", "1103", "1104", "1105"]
        result = assign_tiles_to_samples(pool, 5)
        assert [tiles[0] for tiles in result] == pool

    def test_single_sample_gets_entire_pool(self):
        """A lane with a single sample assigns the whole pool to it."""
        pool = [f"11{i:02d}" for i in range(1, 9)]
        result = assign_tiles_to_samples(pool, 1)
        assert result[0] == pool

    @pytest.mark.parametrize("pool_size", [1, 5, 16, 64])
    @pytest.mark.parametrize("n_samples", [1, 3, 16, 40])
    def test_every_sample_and_tile_is_used(self, pool_size, n_samples):
        """Every sample gets at least one tile and every tile is assigned."""
        pool = [f"11{i:02d}" for i in range(1, pool_size + 1)]
        result = assign_tiles_to_samples(pool, n_samples)
        assert len(result) == n_samples
        assert all(tiles for tiles in result)
        assert {t for tiles in result for t in tiles} == set(pool)

    def test_zero_samples_raises(self):
        """Zero samples raises ValueError."""
        with pytest.raises(ValueError, match="n_samples"):
            assign_tiles_to_samples(["1101"], 0)

    def test_empty_pool_raises(self):
        """An empty tile pool raises ValueError."""
        with pytest.raises(ValueError, match="tile_pool"):
            assign_tiles_to_samples([], 1)


class TestTilesPerLane:
    """Tests for tiles_per_lane."""

    def test_uses_full_pool_when_all_lanes_can_cover_it(self):
        """No cap applies when every lane has enough sequences for its pool."""
        assert tiles_per_lane(16, [8, 4, 5], 10) == 16

    def test_caps_to_smallest_lane_sequence_count(self):
        """A lane with 1 sample x 10 sequences caps every lane to 10 tiles."""
        assert tiles_per_lane(64, [1, 4, 8], 10) == 10

    def test_pool_size_is_the_upper_bound(self):
        """The cap never exceeds the available pool size."""
        assert tiles_per_lane(8, [4, 4], 100) == 8

    def test_mixed_indexes_sheet_at_seq_10(self):
        """Mirrors assets/SampleSheet_MixedIndexes.csv (8 lanes, 16 tiles each)."""
        counts = [1, 8, 4, 5, 11, 4, 8, 6]
        assert tiles_per_lane(16, counts, 10) == 10

    def test_large_seq_number_uses_full_pool(self):
        """A seq_number large enough to cover the pool leaves it uncapped."""
        counts = [1, 8, 4, 5, 11, 4, 8, 6]
        assert tiles_per_lane(16, counts, 1000) == 16


# ---------------------------------------------------------------------------
# End-to-end tile distribution with a sample sheet
# ---------------------------------------------------------------------------

TWO_LANE_SAMPLE_SHEET = """\
[Header]
FileFormatVersion,2
RunName,DummyRun
InstrumentPlatform,NovaSeqXSeries
[Reads]
Read1Cycles,151
Read2Cycles,151
Index1Cycles,10
Index2Cycles,10
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
1,S1,Sample1,GAACTGAGCG,CGCTCCACGA,ProjA,Y151;I10;I10;Y151
1,S2,Sample2,AGGTCAGATA,TATCTTGTAG,ProjA,Y151;I10;I10;Y151
1,S3,Sample3,TAAACCCTAG,TTCCTATCAG,ProjB,Y151;I10;I10;Y151
1,S4,Sample4,CGTCTCATAT,AGCTACTATA,ProjB,Y151;I10;I10;Y151
2,S5,Sample5,GACGAGATTA,GCGCGGTTAA,ProjA,Y151;I10;I10;Y151
2,S6,Sample6,ATTCCATAAG,CCACCAGGCA,ProjA,Y151;I10;I10;Y151
2,S7,Sample7,CGAGGCTGAC,AGGATCTGAC,ProjB,Y151;I10;I10;Y151
2,S8,Sample8,GCTACGCTAC,TATCTTGTAG,ProjB,Y151;I10;I10;Y151
"""

CONSTRAINED_LANE_SAMPLE_SHEET = """\
[Header]
FileFormatVersion,2
RunName,DummyRun
InstrumentPlatform,NovaSeqXSeries
[Reads]
Read1Cycles,151
Read2Cycles,151
Index1Cycles,10
Index2Cycles,10
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
1,S1,Sample1,GAACTGAGCG,CGCTCCACGA,ProjA,Y151;I10;I10;Y151
2,S2,Sample2,AGGTCAGATA,TATCTTGTAG,ProjA,Y151;I10;I10;Y151
2,S3,Sample3,TAAACCCTAG,TTCCTATCAG,ProjA,Y151;I10;I10;Y151
2,S4,Sample4,CGTCTCATAT,AGCTACTATA,ProjB,Y151;I10;I10;Y151
2,S5,Sample5,GACGAGATTA,GCGCGGTTAA,ProjB,Y151;I10;I10;Y151
"""


def read_lane_tiles(demux_path: pathlib.Path) -> dict[str, set]:
    """Extract the distinct tiles used by each lane from the R1 FASTQ headers.

    The read name has the structure ``prefix:lane:tile:x:y`` optionally
    followed by a UMI suffix (``:umi``) and/or the `` 1:N:0:idx`` extension,
    so the lane and tile are captured with a regular expression that is
    robust against both.
    """
    tiles_by_lane: dict[str, set] = {}
    pattern = re.compile(r":(\d+):(\d{4}):\d{4}:\d{4}(?::\S+)?(?: 1:N:0:\S+)?$")
    for fastq in sorted(demux_path.rglob("*_R1_001.fastq.gz")):
        with gzip.open(fastq, "rt") as fh:
            for line in fh:
                if not line.startswith("@"):
                    continue
                match = pattern.search(line)
                assert match, f"Could not parse lane/tile from read name: {line}"
                lane, tile = match.groups()
                tiles_by_lane.setdefault(lane, set()).add(tile)
    return tiles_by_lane


class TestSampleSheetTileDistribution:
    """Tests for the tile distribution driven by a sample sheet."""

    def test_equal_tiles_per_lane_and_full_coverage(self, tmp_path):
        """Each lane uses an equal number of tiles and all 128 tiles are covered."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(TWO_LANE_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=16))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        tiles_by_lane = read_lane_tiles(demux_path)

        # Both lanes are present and use an equal number of tiles (64 = 128/2)
        assert set(tiles_by_lane) == {"1", "2"}
        assert all(len(tiles) == 64 for tiles in tiles_by_lane.values())

        # The lanes use disjoint tile sets whose union is the full tile set
        lane1, lane2 = tiles_by_lane["1"], tiles_by_lane["2"]
        assert not (lane1 & lane2)
        assert lane1 | lane2 == all_tile_ids()

    def test_every_sample_produces_reads(self, tmp_path):
        """All 8 samples produce R1 files with 16 reads each."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(TWO_LANE_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=16))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        r1_files = sorted(demux_path.rglob("*_R1_001.fastq.gz"))
        assert len(r1_files) == 8
        for fastq in r1_files:
            assert _count_reads(fastq) == 16

    def test_tile_count_capped_when_lane_is_constrained(self, tmp_path):
        """When a lane has fewer sequences than tiles, every lane is capped to
        the same number of tiles and every used tile gets at least one
        sequence."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(CONSTRAINED_LANE_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=10))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        tiles_by_lane = read_lane_tiles(demux_path)

        # Lane 1 has only 10 sequences, so both lanes are capped to 10 tiles
        assert set(tiles_by_lane) == {"1", "2"}
        assert all(len(tiles) == 10 for tiles in tiles_by_lane.values())
        # Lanes stay disjoint and use 20 tiles in total
        lane1, lane2 = tiles_by_lane["1"], tiles_by_lane["2"]
        assert not (lane1 & lane2)
        assert len(lane1 | lane2) == 20

        # Every sample still produces its full number of reads
        r1_files = sorted(demux_path.rglob("*_R1_001.fastq.gz"))
        assert _count_reads(r1_files[0]) == 10
        assert all(_count_reads(f) == 10 for f in r1_files[1:])


def _count_reads(fastq: pathlib.Path) -> int:
    """Count the number of records in a gzip-compressed FASTQ file."""
    with gzip.open(fastq, "rt") as fh:
        return sum(1 for line in fh if line.startswith("@"))


def _blabber_args(
    tmp_path: pathlib.Path, sheet: pathlib.Path, seq_number: int, taint: bool = False
) -> argparse.Namespace:
    """Build the argument namespace used to run blabber.main in the tests."""
    return argparse.Namespace(
        verbose=False,
        quiet=True,
        random_seed=42,
        alphabet="ACGT",
        format="fastq",
        sample_sheet=sheet,
        output=tmp_path,
        flowcell_id="TESTFC",
        seq_number=seq_number,
        seq_length=151,
        seq_mask=None,
        index1=None,
        index2=None,
        taint=taint,
    )


# ---------------------------------------------------------------------------
# Single-end recipes and taint (Undetermined files)
# ---------------------------------------------------------------------------

SINGLE_END_SAMPLE_SHEET = """\
[Header]
FileFormatVersion,2
RunName,DummyRun
InstrumentPlatform,NextSeq2000
[Reads]
Read1Cycles,151
Read2Cycles,151
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
1,S1,Sample1,,,ProjA,Y151;Y151
"""


class TestSingleEndAndTaint:
    """Tests for single-end generation and base-calling taints."""

    def test_single_end_generates_r1_only(self, tmp_path):
        """A single-end recipe produces an R1 file with reads and no R2 file."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(SINGLE_END_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=16))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        r1_files = sorted(demux_path.rglob("*_R1_001.fastq.gz"))
        r2_files = sorted(demux_path.rglob("*_R2_001.fastq.gz"))
        assert len(r1_files) == 1
        assert not r2_files
        assert _count_reads(r1_files[0]) == 16

    def test_single_end_taint_creates_r1_undetermined_only(self, tmp_path):
        """Tainting a single-end lane creates a non-empty R1 Undetermined file only."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(SINGLE_END_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=100, taint=True))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        undetermined = sorted(demux_path.glob("Undetermined*"))
        assert [f.name for f in undetermined] == [
            "Undetermined_S0_L001_R1_001.fastq.gz"
        ]
        assert _count_reads(undetermined[0]) > 0

    def test_paired_end_taint_creates_both_undetermined_files(self, tmp_path):
        """Tainting paired-end lanes creates R1 and R2 Undetermined files with
        the same number of reads in each."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(TWO_LANE_SAMPLE_SHEET)
        main(_blabber_args(tmp_path, sheet, seq_number=100, taint=True))

        demux_path = tmp_path / "TESTFC" / "Demultiplexing"
        for lane in ["L001", "L002"]:
            r1 = demux_path / f"Undetermined_S0_{lane}_R1_001.fastq.gz"
            r2 = demux_path / f"Undetermined_S0_{lane}_R2_001.fastq.gz"
            assert r1.is_file(), f"Missing {r1.name}"
            assert r2.is_file(), f"Missing {r2.name}"
            assert _count_reads(r1) > 0
            assert _count_reads(r1) == _count_reads(r2)
