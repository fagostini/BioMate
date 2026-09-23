"""Tests for the fastrewind module."""

import argparse
from contextlib import contextmanager
import gzip
import pathlib
import struct

import pytest
import biomate.fastrewind.fastrewind as fastrewind_module

from biomate.fastrewind.fastrewind import (
    INTEROP_INDEX_VERSION,
    build_index_records,
    clean_directory,
    parse_fastq_groups,
    parse_sequence_mask,
    preprocess_and_write_bcls,
    validate_args,
    write_index_metrics,
)


# ---------------------------------------------------------------------------
# parse_sequence_mask
# ---------------------------------------------------------------------------


class TestParseSequenceMask:
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
        result = parse_sequence_mask("Y151;I8N2;Y151")
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_prefix_and_suffix_around_index(self):
        """N before I sets I1B; N after I sets I1A."""
        result = parse_sequence_mask("Y151;N1I8N2;Y151")
        assert result["I1B"] == 1
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_umi_first_position(self):
        """A UMI prefix (U8) before a read segment is parsed into U1."""
        result = parse_sequence_mask("U8Y143;I10;I10;Y151")
        assert result["U1"] == 8
        assert result["R1"] == 143

    def test_too_many_sections_raises(self):
        """A mask with more than 4 semicolon-separated sections raises ValueError."""
        with pytest.raises(ValueError, match="OverrideCycles"):
            parse_sequence_mask("Y151;I10;I10;Y151;Y50")

    def test_all_zero_defaults(self):
        """Unused mask fields default to 0 in the returned dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert result["U1"] == 0
        assert result["U2"] == 0
        assert result["R1B"] == 0
        assert result["R1A"] == 0
        assert result["I1B"] == 0
        assert result["I2B"] == 0

    def test_return_type_is_dict(self):
        """parse_sequence_mask returns a dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# clean_directory
# ---------------------------------------------------------------------------


class TestCleanDirectory:
    def test_removes_a_file(self, tmp_path):
        """A plain file is deleted."""
        f = tmp_path / "test.txt"
        f.write_text("data")
        clean_directory(f)
        assert not f.exists()

    def test_removes_empty_directory(self, tmp_path):
        """An empty directory is deleted."""
        d = tmp_path / "emptydir"
        d.mkdir()
        clean_directory(d)
        assert not d.exists()

    def test_removes_directory_with_files(self, tmp_path):
        """A directory containing files is deleted recursively."""
        d = tmp_path / "mydir"
        d.mkdir()
        (d / "a.txt").write_text("hello")
        (d / "b.txt").write_text("world")
        clean_directory(d)
        assert not d.exists()

    def test_removes_nested_directories(self, tmp_path):
        """A deeply nested directory tree is deleted completely."""
        root = tmp_path / "root"
        root.mkdir()
        nested = root / "a" / "b" / "c"
        nested.mkdir(parents=True)
        (nested / "file.bin").write_bytes(b"\x00\x01")
        clean_directory(root)
        assert not root.exists()


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _make_args(self, tmp_path, force=False):
        input_path = tmp_path / "flowcell"
        input_path.mkdir()
        # Add a dummy FASTQ file so the glob check passes
        (input_path / "sample.fastq.gz").write_bytes(b"")
        # Add a sample sheet
        (input_path / "SampleSheet.csv").write_text("[BCLConvert_Data]\n")
        return argparse.Namespace(
            input_path=input_path,
            output_path=tmp_path / "output",
            sample_sheet=None,
            threads=0,
            force=force,
            interop_dir=None,
        )

    def test_valid_args_returned(self, tmp_path):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._make_args(tmp_path)
        result = validate_args(args)
        assert result is args

    def test_nonexistent_input_path_raises(self, tmp_path):
        """A non-existent input path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.input_path = tmp_path / "missing"
        with pytest.raises(argparse.ArgumentTypeError, match="[Vv]alid directory"):
            validate_args(args)

    def test_input_path_is_file_raises(self, tmp_path):
        """A plain file passed as input path raises ArgumentTypeError."""
        f = tmp_path / "not_a_dir.txt"
        f.write_text("")
        args = self._make_args(tmp_path)
        args.input_path = f
        with pytest.raises(argparse.ArgumentTypeError, match="[Vv]alid directory"):
            validate_args(args)

    def test_no_fastq_files_raises(self, tmp_path):
        """A directory with no FASTQ files raises ArgumentTypeError."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        (empty_dir / "SampleSheet.csv").write_text("")
        args = self._make_args(tmp_path)
        args.input_path = empty_dir
        with pytest.raises(argparse.ArgumentTypeError, match="[Nn]o FASTQ"):
            validate_args(args)

    def test_missing_sample_sheet_raises(self, tmp_path):
        """Absence of both SampleSheet.csv and --sample-sheet raises ArgumentTypeError."""
        input_path = tmp_path / "flowcell2"
        input_path.mkdir()
        (input_path / "sample.fastq.gz").write_bytes(b"")
        # No SampleSheet.csv and no --sample-sheet flag
        args = argparse.Namespace(
            input_path=input_path,
            output_path=tmp_path / "output",
            sample_sheet=None,
            threads=0,
            force=False,
            interop_dir=None,
        )
        with pytest.raises(argparse.ArgumentTypeError, match="[Ss]ample[Ss]heet"):
            validate_args(args)

    def test_explicit_sample_sheet_nonexistent_raises(self, tmp_path):
        """A non-existent --sample-sheet path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.sample_sheet = tmp_path / "nonexistent.csv"
        with pytest.raises(argparse.ArgumentTypeError, match="[Ss]ample"):
            validate_args(args)

    def test_force_cleans_default_interop_dir(self, tmp_path):
        """With --force, an existing InterOp directory is cleaned like Data."""
        args = self._make_args(tmp_path, force=True)
        (args.output_path / "Data").mkdir(parents=True)
        (args.output_path / "InterOp").mkdir(parents=True)
        (args.output_path / "InterOp" / "IndexMetricsOut.bin").write_bytes(b"\x02")
        validate_args(args)
        assert not (args.output_path / "Data").exists()
        assert not (args.output_path / "InterOp").exists()

    def test_force_cleans_custom_interop_dir(self, tmp_path):
        """With --force, a custom --interop-dir location is cleaned."""
        args = self._make_args(tmp_path, force=True)
        (args.output_path / "Data").mkdir(parents=True)
        custom = tmp_path / "custom_interop"
        custom.mkdir()
        (custom / "IndexMetricsOut.bin").write_bytes(b"\x02")
        args.interop_dir = custom
        validate_args(args)
        assert not custom.exists()

    def test_no_force_keeps_existing_interop_dir(self, tmp_path):
        """Without --force an existing InterOp directory is left untouched."""
        args = self._make_args(tmp_path, force=False)
        interop_dir = args.output_path / "InterOp"
        interop_dir.mkdir(parents=True)
        (interop_dir / "IndexMetricsOut.bin").write_bytes(b"\x02")
        validate_args(args)
        assert interop_dir.exists()


# ---------------------------------------------------------------------------
# preprocess_and_write_bcls
# ---------------------------------------------------------------------------


def test_preprocess_and_write_bcls_handles_odd_read_count(tmp_path, monkeypatch):
    """An odd number of reads still produces valid cbcl files for each cycle."""

    class FakeBaseCall:
        def __init__(self, base: str, quality: str):
            self.sequence = base
            self.qualities = quality

    class FakeRead:
        def __init__(self, sequence: str, qualities: str):
            self.sequence = sequence
            self.qualities = qualities

        def __len__(self):
            return len(self.sequence)

        def __getitem__(self, idx):
            return FakeBaseCall(self.sequence[idx], self.qualities[idx])

    lane = "L002"
    tiles_root = tmp_path / "tiles"
    lane_dir = tiles_root / lane
    lane_dir.mkdir(parents=True)
    fq_path = lane_dir / "1_1001.fastq.gz"
    fq_path.write_bytes(b"")

    fake_reads = [
        FakeRead("AC", "II"),
        FakeRead("TG", "II"),
        FakeRead("NN", "!!"),
    ]

    @contextmanager
    def fake_dnaio_open(path, open_threads=0):
        assert pathlib.Path(path) == fq_path
        assert open_threads == 0
        yield iter(fake_reads)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    out_dir = tmp_path / "out"
    preprocess_and_write_bcls(
        output_path=out_dir,
        lane=lane,
        tiles_path=tiles_root,
        total_cycles=2,
        threads=0,
        pattern_suffix="*_*.fastq.gz",
    )

    cycle_1 = out_dir / "Data/Intensities/BaseCalls/L002/C1.1/L002_1.cbcl"
    cycle_2 = out_dir / "Data/Intensities/BaseCalls/L002/C2.1/L002_1.cbcl"

    assert cycle_1.exists()
    assert cycle_2.exists()
    assert cycle_1.stat().st_size > 0
    assert cycle_2.stat().st_size > 0


def _cbcl_body(path: pathlib.Path) -> bytes:
    """Decompress the body of the first tile in a cbcl file."""
    data = path.read_bytes()
    hsize = struct.unpack_from("<HLBBI", data, 0)[1]
    ntiles, = struct.unpack_from("<I", data, 44)
    assert ntiles == 1
    _, _, _, comp = struct.unpack_from("<IIII", data, 48)
    return gzip.decompress(data[hsize : hsize + comp])


def test_preprocess_and_write_bcls_packs_clusters_little_endian(tmp_path, monkeypatch):
    """The first cluster of each byte goes in the low nibble and the odd tail
    is zero-padded in the high nibble (the order bcl-convert expects)."""

    class FakeBaseCall:
        def __init__(self, base: str, quality: str):
            self.sequence = base
            self.qualities = quality

    class FakeRead:
        def __init__(self, sequence: str, qualities: str):
            self.sequence = sequence
            self.qualities = qualities

        def __len__(self):
            return len(self.sequence)

        def __getitem__(self, idx):
            return FakeBaseCall(self.sequence[idx], self.qualities[idx])

    lane = "L002"
    tiles_root = tmp_path / "tiles"
    lane_dir = tiles_root / lane
    lane_dir.mkdir(parents=True)
    fq_path = lane_dir / "1_1001.fastq.gz"
    fq_path.write_bytes(b"")

    # "II" is Phred 40 -> quality bits "11", "!!" is Phred 0 -> "00";
    # base map: A="00" C="01" G="10" T="11", N -> "00"
    fake_reads = [
        FakeRead("AC", "II"),
        FakeRead("TG", "II"),
        FakeRead("NN", "!!"),
    ]

    @contextmanager
    def fake_dnaio_open(path, open_threads=0):
        assert pathlib.Path(path) == fq_path
        assert open_threads == 0
        yield iter(fake_reads)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    out_dir = tmp_path / "out"
    preprocess_and_write_bcls(
        output_path=out_dir,
        lane=lane,
        tiles_path=tiles_root,
        total_cycles=2,
        threads=0,
        pattern_suffix="*_*.fastq.gz",
    )

    cycle_1 = out_dir / "Data/Intensities/BaseCalls/L002/C1.1/L002_1.cbcl"
    cycle_2 = out_dir / "Data/Intensities/BaseCalls/L002/C2.1/L002_1.cbcl"

    # Cycle 1: A="1100" (low) | T="1111" (high) -> 0xFC, N tail -> 0x00
    assert _cbcl_body(cycle_1) == b"\xfc\x00"
    # Cycle 2: C="1101" (low) | G="1110" (high) -> 0xED, N tail -> 0x00
    assert _cbcl_body(cycle_2) == b"\xed\x00"


def test_parse_fastq_groups_uses_only_r1_r2_with_extra_files(tmp_path, monkeypatch):
    """Index-read files and duplicate R1/R2 files are ignored during paired-end parsing."""

    class FakeRead:
        def __init__(self, name: str, sequence: str, qualities: str):
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        def write(self, read):
            return None

    class FakeWriterContext:
        def __enter__(self):
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    class FakeReaderContext:
        def __init__(self, records):
            self.records = records

        def __enter__(self):
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    sample_name = "Sample_1_S1_L001_001"
    lane = "L001"
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)

    # Input includes R1/R2 duplicates and index-read files; only one R1 and one R2 must be used.
    filenames = [
        tmp_path / "A" / "Sample_1_S1_L001_R1_001.fastq.gz",
        tmp_path / "A" / "Sample_1_S1_L001_R2_001.fastq.gz",
        tmp_path / "A" / "Sample_1_S1_L001_I1_001.fastq.gz",
        tmp_path / "A" / "Sample_1_S1_L001_I2_001.fastq.gz",
        tmp_path / "B" / "Sample_1_S1_L001_R1_001.fastq.gz",
        tmp_path / "B" / "Sample_1_S1_L001_R2_001.fastq.gz",
    ]

    read_name = "INST:1:FLOWCELL:1:2101:1001:1002 1:N:0:ACGT+TGCA"
    records = [
        (
            FakeRead(read_name, "AC", "II"),
            FakeRead(read_name.replace(" 1:", " 2:"), "GT", "II"),
        )
    ]

    reader_calls = []

    def fake_dnaio_open(*args, **kwargs):
        if "mode" in kwargs:
            return FakeWriterContext()
        reader_calls.append(tuple(pathlib.Path(x).name for x in args))
        return FakeReaderContext(records)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    masks_table = {
        sample_name: {
            "R1B": "",
            "R1": 2,
            "R1A": "",
            "I1": "ACGT",
            "I2": "TGCA",
            "R2B": "",
            "R2": 2,
            "R2A": "",
        }
    }

    result = parse_fastq_groups(
        lane=lane,
        name=sample_name,
        filenames=filenames,
        masks_table=masks_table,
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
    )

    assert len(reader_calls) == 1
    assert len(reader_calls[0]) == 2
    assert "_R1_" in reader_calls[0][0] or "_R1_" in reader_calls[0][1]
    assert "_R2_" in reader_calls[0][0] or "_R2_" in reader_calls[0][1]
    assert all("_I1_" not in name and "_I2_" not in name for name in reader_calls[0])

    assert "2101" in result
    assert (1001, 1002) in result["2101"]


def test_parse_fastq_groups_skips_when_no_r1_r2(tmp_path, monkeypatch):
    """If a sample group has no R1/R2 files, parsing is skipped and no reads are opened."""

    class FakeReaderContext:
        def __enter__(self):
            raise AssertionError(
                "Reader should not be opened when no R1/R2 files exist"
            )

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    sample_name = "Sample_2_S2_L001_001"
    lane = "L001"
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)

    filenames = [
        tmp_path / "A" / "Sample_2_S2_L001_I1_001.fastq.gz",
        tmp_path / "A" / "Sample_2_S2_L001_I2_001.fastq.gz",
    ]

    def fake_dnaio_open(*args, **kwargs):
        return FakeReaderContext()

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    masks_table = {
        sample_name: {
            "R1B": "",
            "R1": 2,
            "R1A": "",
            "I1": "AAAA",
            "I2": "CCCC",
            "R2B": "",
            "R2": 2,
            "R2A": "",
        }
    }

    result = parse_fastq_groups(
        lane=lane,
        name=sample_name,
        filenames=filenames,
        masks_table=masks_table,
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
    )

    assert dict(result) == {}


# ---------------------------------------------------------------------------
# write_index_metrics
# ---------------------------------------------------------------------------


def _unpack_index_metrics(data: bytes, version: int) -> list:
    """Parse the bytes of an IndexMetricsOut.bin file into record tuples."""
    assert data[0] == version
    id_format = "<HHH" if version == 1 else "<HIH"
    count_format = "<I" if version == 1 else "<Q"
    records = []
    offset = 1
    while offset < len(data):
        lane, tile, read = struct.unpack_from(id_format, data, offset)
        offset += struct.calcsize(id_format)
        (index_len,) = struct.unpack_from("<H", data, offset)
        offset += 2
        index = data[offset : offset + index_len].decode("utf-8")
        offset += index_len
        (count,) = struct.unpack_from(count_format, data, offset)
        offset += struct.calcsize(count_format)
        (sample_len,) = struct.unpack_from("<H", data, offset)
        offset += 2
        sample = data[offset : offset + sample_len].decode("utf-8")
        offset += sample_len
        (project_len,) = struct.unpack_from("<H", data, offset)
        offset += 2
        project = data[offset : offset + project_len].decode("utf-8")
        offset += project_len
        records.append((lane, tile, read, index, count, sample, project))
    return records


class TestWriteIndexMetrics:
    RECORDS = [
        (1, 1101, 1, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
        (1, 1101, 2, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
        (2, 2102, 1, "AAAA", 7, "Sample_2", "Proj_B"),
    ]

    @pytest.mark.parametrize("version", [1, 2])
    def test_round_trip(self, tmp_path, version):
        """Written bytes can be parsed back into the original records."""
        path = tmp_path / "InterOp" / "IndexMetricsOut.bin"
        write_index_metrics(path, self.RECORDS, version=version)
        assert _unpack_index_metrics(path.read_bytes(), version) == self.RECORDS

    def test_version_byte(self, tmp_path):
        """The first byte of the file is the format version."""
        path = tmp_path / "IndexMetricsOut.bin"
        write_index_metrics(path, self.RECORDS, version=2)
        assert path.read_bytes()[0] == 2

    def test_default_version_matches_constant(self, tmp_path):
        """The default version is the module INTEROP_INDEX_VERSION constant."""
        path = tmp_path / "IndexMetricsOut.bin"
        write_index_metrics(path, self.RECORDS)
        assert path.read_bytes()[0] == INTEROP_INDEX_VERSION

    def test_invalid_version_raises(self, tmp_path):
        """An unsupported version raises ValueError before anything is written."""
        with pytest.raises(ValueError, match="[Vv]ersion"):
            write_index_metrics(
                tmp_path / "IndexMetricsOut.bin", self.RECORDS, version=3
            )

    def test_empty_records_writes_version_only(self, tmp_path):
        """No records results in a one-byte (version only) file."""
        path = tmp_path / "IndexMetricsOut.bin"
        write_index_metrics(path, [], version=2)
        assert path.read_bytes() == b"\x02"

    def test_empty_strings_default_to_na(self, tmp_path):
        """Empty index/sample/project strings are written as 'NA'."""
        path = tmp_path / "IndexMetricsOut.bin"
        write_index_metrics(path, [(1, 1101, 1, "", 3, "", "")], version=2)
        assert _unpack_index_metrics(path.read_bytes(), 2) == [
            (1, 1101, 1, "NA", 3, "NA", "NA")
        ]

    def test_creates_parent_directories(self, tmp_path):
        """Missing parent directories of the target file are created."""
        path = tmp_path / "a" / "b" / "IndexMetricsOut.bin"
        write_index_metrics(path, self.RECORDS)
        assert path.exists()


# ---------------------------------------------------------------------------
# build_index_records
# ---------------------------------------------------------------------------


class TestBuildIndexRecords:
    MASKS_TABLE = {
        "Sample_1_S1_L001_001": {
            "sample_name": "Sample_1",
            "sample_project": "Proj_A",
            "I1": "ACGT",
            "I2": "TGCA",
        },
        "Sample_2_S2_L001_001": {
            "sample_name": "Sample_2",
            "sample_project": "Proj_B",
            "I1": "AAAA",
            "I2": "CCCC",
        },
        "Sample_3_S3_L002_001": {
            "sample_name": "Sample_3",
            "sample_project": "Proj_C",
            "I1": "GGGG",
            "I2": "",
        },
        "Sample_4_S4_L001_001": {
            "sample_name": "Sample_4",
            "sample_project": "Proj_D",
            "I1": "",
            "I2": "",
        },
    }

    def test_dual_index_expands_to_one_record_per_index_read(self):
        """A dual-index sample yields one record per index read with the
        combined 'I1-I2' name and the same cluster count."""
        counts = {("L001", 1101, "Sample_1_S1_L001_001"): 42}
        records = build_index_records(counts, self.MASKS_TABLE, {"I1": 8, "I2": 8})
        assert records == [
            (1, 1101, 1, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
            (1, 1101, 2, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
        ]

    def test_single_index_yields_one_record_per_run_index_read(self):
        """A single-index sample yields one record per index read of the run,
        all carrying the bare I1 name and the same count."""
        counts = {("L002", 2101, "Sample_3_S3_L002_001"): 3}
        records = build_index_records(counts, self.MASKS_TABLE, {"I1": 8, "I2": 8})
        assert records == [
            (2, 2101, 1, "GGGG", 3, "Sample_3", "Proj_C"),
            (2, 2101, 2, "GGGG", 3, "Sample_3", "Proj_C"),
        ]

    def test_sample_without_indexes_is_skipped(self):
        """A sample without indexes produces no records."""
        counts = {("L001", 1101, "Sample_4_S4_L001_001"): 9}
        records = build_index_records(counts, self.MASKS_TABLE, {"I1": 8, "I2": 8})
        assert records == []

    def test_run_without_index_reads_produces_no_records(self):
        """A run without index cycles produces no records."""
        counts = {("L001", 1101, "Sample_1_S1_L001_001"): 9}
        records = build_index_records(counts, self.MASKS_TABLE, {"R1": 10, "R2": 10})
        assert records == []

    def test_read_numbers_are_index_ordinals(self):
        """Read numbers are the position among the run's index reads (I1 -> 1,
        I2 -> 2), independent of the run's overall read numbering."""
        counts = {
            ("L001", 1101, "Sample_1_S1_L001_001"): 1,
        }
        records = build_index_records(
            counts, self.MASKS_TABLE, {"R1": 151, "I1": 8, "I2": 8, "R2": 151}
        )
        assert [r[2] for r in records] == [1, 2]

    def test_records_ordered_by_lane_sample_tile_read(self):
        """Records are sorted by lane, sample, tile and read like instrument files."""
        counts = {
            ("L002", 1101, "Sample_3_S3_L002_001"): 3,
            ("L001", 1102, "Sample_2_S2_L001_001"): 7,
            ("L001", 1101, "Sample_1_S1_L001_001"): 42,
        }
        records = build_index_records(counts, self.MASKS_TABLE, {"I1": 8, "I2": 8})
        assert records == [
            (1, 1101, 1, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
            (1, 1101, 2, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
            (1, 1102, 1, "AAAA-CCCC", 7, "Sample_2", "Proj_B"),
            (1, 1102, 2, "AAAA-CCCC", 7, "Sample_2", "Proj_B"),
            (2, 1101, 1, "GGGG", 3, "Sample_3", "Proj_C"),
            (2, 1101, 2, "GGGG", 3, "Sample_3", "Proj_C"),
        ]

    def test_unknown_group_and_non_positive_count_are_skipped(self):
        """Groups missing from the masks table or with no clusters are skipped."""
        counts = {
            ("L001", 1101, "Unknown_S9_L001_001"): 5,
            ("L001", 1101, "Sample_1_S1_L001_001"): 0,
        }
        records = build_index_records(counts, self.MASKS_TABLE, {"I1": 8, "I2": 8})
        assert records == []


# ---------------------------------------------------------------------------
# parse_fastq_groups index counts
# ---------------------------------------------------------------------------


def test_parse_fastq_groups_accumulates_index_counts(tmp_path, monkeypatch):
    """Cluster counts are accumulated per (lane, tile, sample group) when an
    index_counts dict is provided."""

    class FakeRead:
        def __init__(self, name: str, sequence: str, qualities: str):
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        def write(self, read):
            return None

    class FakeWriterContext:
        def __enter__(self):
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    class FakeReaderContext:
        def __init__(self, records):
            self.records = records

        def __enter__(self):
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    sample_name = "Sample_1_S1_L001_001"
    lane = "L001"
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)

    def make_read(tile: str, read: str) -> str:
        return f"INST:1:FLOWCELL:1:{tile}:1001:1002 {read}:N:0:ACGT+TGCA"

    records = [
        (
            FakeRead(make_read("1101", "1"), "AC", "II"),
            FakeRead(make_read("1101", "2"), "GT", "II"),
        ),
        (
            FakeRead(make_read("1101", "1"), "AC", "II"),
            FakeRead(make_read("1101", "2"), "GT", "II"),
        ),
        (
            FakeRead(make_read("1102", "1"), "AC", "II"),
            FakeRead(make_read("1102", "2"), "GT", "II"),
        ),
    ]

    def fake_dnaio_open(*args, **kwargs):
        if "mode" in kwargs:
            return FakeWriterContext()
        return FakeReaderContext(records)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    masks_table = {
        sample_name: {
            "R1B": "",
            "R1": 2,
            "R1A": "",
            "I1": "ACGT",
            "I2": "TGCA",
            "R2B": "",
            "R2": 2,
            "R2A": "",
        }
    }

    index_counts = {}
    parse_fastq_groups(
        lane=lane,
        name=sample_name,
        filenames=[
            tmp_path / "A" / "Sample_1_S1_L001_R1_001.fastq.gz",
            tmp_path / "A" / "Sample_1_S1_L001_R2_001.fastq.gz",
        ],
        masks_table=masks_table,
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
        index_counts=index_counts,
    )

    assert index_counts == {
        (lane, 1101, sample_name): 2,
        (lane, 1102, sample_name): 1,
    }


def test_parse_fastq_groups_without_index_counts_unchanged(tmp_path, monkeypatch):
    """Omitting index_counts keeps the previous behaviour (no accumulation)."""

    class FakeRead:
        def __init__(self, name: str, sequence: str, qualities: str):
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        def write(self, read):
            return None

    class FakeWriterContext:
        def __enter__(self):
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    class FakeReaderContext:
        def __init__(self, records):
            self.records = records

        def __enter__(self):
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            return False

    sample_name = "Sample_1_S1_L001_001"
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)
    read_name = "INST:1:FLOWCELL:1:1101:1001:1002 1:N:0:ACGT+TGCA"
    records = [
        (
            FakeRead(read_name, "AC", "II"),
            FakeRead(read_name.replace(" 1:", " 2:"), "GT", "II"),
        )
    ]

    def fake_dnaio_open(*args, **kwargs):
        if "mode" in kwargs:
            return FakeWriterContext()
        return FakeReaderContext(records)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)

    masks_table = {
        sample_name: {
            "R1B": "",
            "R1": 2,
            "R1A": "",
            "I1": "ACGT",
            "I2": "TGCA",
            "R2B": "",
            "R2": 2,
            "R2A": "",
        }
    }

    result = parse_fastq_groups(
        lane="L001",
        name=sample_name,
        filenames=[
            tmp_path / "A" / "Sample_1_S1_L001_R1_001.fastq.gz",
            tmp_path / "A" / "Sample_1_S1_L001_R2_001.fastq.gz",
        ],
        masks_table=masks_table,
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
    )

    assert "1101" in result
    assert (1001, 1002) in result["1101"]


# ---------------------------------------------------------------------------
# InterOp read-back with the official Illumina interop library
# ---------------------------------------------------------------------------


class TestInteropReadBack:
    """Validate the generated IndexMetricsOut.bin with the official Illumina
    interop library (https://github.com/Illumina/interop)."""

    RUN_INFO = """<?xml version="1.0"?>
<RunInfo Version="6">
  <Run Id="260921LH00217_0424_B23NWY7LT4" Number="424">
    <Flowcell>B23NWY7LT4</Flowcell>
    <Instrument>LH00217</Instrument>
    <Date>2026-09-21T13:49:23Z</Date>
    <Reads>
      <Read Number="1" NumCycles="6" IsIndexedRead="N" IsReverseComplement="N"/>
      <Read Number="2" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="N"/>
      <Read Number="3" NumCycles="8" IsIndexedRead="Y" IsReverseComplement="Y"/>
      <Read Number="4" NumCycles="6" IsIndexedRead="N" IsReverseComplement="N"/>
    </Reads>
    <FlowcellLayout LaneCount="1" SurfaceCount="2" SwathCount="4" TileCount="16">
      <TileSet TileNamingConvention="FourDigit">
        <Tiles>
          <Tile>1_1101</Tile>
          <Tile>1_1102</Tile>
        </Tiles>
      </TileSet>
    </FlowcellLayout>
    <ImageDimensions Width="320" Height="160"/>
    <ImageChannels>
      <Name>blue</Name>
      <Name>green</Name>
    </ImageChannels>
  </Run>
</RunInfo>
"""

    def _build_run_folder(self, tmp_path: pathlib.Path, records: list) -> pathlib.Path:
        run_dir = tmp_path / "run"
        (run_dir / "InterOp").mkdir(parents=True)
        # Note: the interop library requires numeric lane prefixes in the tile
        # names (1_1101), unlike the L001_1101 form bcl-convert expects in
        # RunInfo.xml, so a minimal RunInfo.xml is written here directly
        (run_dir / "RunInfo.xml").write_text(self.RUN_INFO)
        write_index_metrics(run_dir / "InterOp" / "IndexMetricsOut.bin", records)
        return run_dir

    def test_index_metrics_are_parsed_by_illumina_interop(self, tmp_path):
        """Records written by write_index_metrics are read back by the
        official interop library with the expected values."""
        interop = pytest.importorskip("interop")
        records = [
            (1, 1101, 1, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
            (1, 1101, 2, "ACGT-TGCA", 42, "Sample_1", "Proj_A"),
            (1, 1102, 1, "AAAA-CCCC", 7, "Sample_2", "Proj_B"),
            (1, 1102, 2, "AAAA-CCCC", 7, "Sample_2", "Proj_B"),
        ]
        run_dir = self._build_run_folder(tmp_path, records)

        table = interop.core.indexing(str(run_dir), per_sample=True)

        parsed = {}
        for row in table:
            if isinstance(row["SampleID"], str):
                parsed[(int(row["Lane"]), int(row["Tile"]))] = (
                    row["Barcode"],
                    row["SampleID"],
                    int(row["Cluster Count"]),
                )
        assert parsed == {
            (1, 1101): ("ACGT-TGCA", "Sample_1", 42),
            (1, 1102): ("AAAA-CCCC", "Sample_2", 7),
        }
