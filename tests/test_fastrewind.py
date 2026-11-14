"""Tests for the fastrewind module."""

import argparse
from contextlib import contextmanager
import gzip
import logging
import pathlib
import re
import struct
import xml.etree.ElementTree as ET

import dnaio
import pytest
import biomate.fastrewind.fastrewind as fastrewind_module

from biomate.blabber.blabber import main as blabber_main
from biomate.fastrewind.fastrewind import (
    INTEROP_INDEX_VERSION,
    build_index_lookup,
    build_index_records,
    clean_directory,
    generate_run_id,
    parse_fastq_groups,
    parse_sequence_mask,
    preprocess_and_write_bcls,
    synthesize_undetermined_recipe,
    validate_args,
    write_index_metrics,
    write_run_info_xml,
    write_run_parameters_xml,
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
        """A two-read mask sets R1 and R2 correctly."""
        result = parse_sequence_mask("Y151;Y151")
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
    """Tests for clean_directory."""

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
    """Tests for validate_args."""

    def _make_args(self, tmp_path, force=False):
        """Create a valid Namespace together with the files validate_args expects."""
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
            instrument="NovaSeqXPlus",
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
            instrument="NovaSeqXPlus",
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
        """With --force, only fastrewind's own file leaves a custom --interop-dir."""
        args = self._make_args(tmp_path, force=True)
        (args.output_path / "Data").mkdir(parents=True)
        custom = tmp_path / "custom_interop"
        custom.mkdir()
        (custom / "IndexMetricsOut.bin").write_bytes(b"\x02")
        (custom / "other.bin").write_bytes(b"\x02")
        args.interop_dir = custom
        validate_args(args)
        assert not (custom / "IndexMetricsOut.bin").exists()
        assert (custom / "other.bin").exists()

    def test_force_does_not_follow_symlinks(self, tmp_path):
        """A symlink inside Data is removed, its target is left alone."""
        args = self._make_args(tmp_path, force=True)
        (args.output_path / "Data").mkdir(parents=True)
        outside = tmp_path / "outside"
        outside.mkdir()
        (outside / "keep.txt").write_text("x")
        (args.output_path / "Data" / "link").symlink_to(outside)
        validate_args(args)
        assert not (args.output_path / "Data").exists()
        assert (outside / "keep.txt").exists()

    def test_unsupported_instrument_rejected_before_cleanup(self, tmp_path):
        """An unsupported instrument fails before --force deletes anything."""
        args = self._make_args(tmp_path, force=True)
        args.instrument = "MiSeq"
        (args.output_path / "Data").mkdir(parents=True)
        with pytest.raises(argparse.ArgumentTypeError):
            validate_args(args)
        assert (args.output_path / "Data").exists()

    def test_no_force_keeps_existing_interop_dir(self, tmp_path):
        """Without --force an existing InterOp directory is left untouched."""
        args = self._make_args(tmp_path, force=False)
        interop_dir = args.output_path / "InterOp"
        interop_dir.mkdir(parents=True)
        (interop_dir / "IndexMetricsOut.bin").write_bytes(b"\x02")
        validate_args(args)
        assert interop_dir.exists()


# ---------------------------------------------------------------------------
# generate_run_id / write_run_info_xml
# ---------------------------------------------------------------------------


class TestGenerateRunId:
    """Tests for generate_run_id."""

    READ_NAME = "LH00217:424:B23NWY7LT4:1:1101:1001:1002 1:N:0:ACGT+TGCA"

    def test_run_id_uses_read_name_fields(self):
        """The run id is built from the date, instrument, run number and flowcell of the read name."""
        result = generate_run_id(self.READ_NAME)
        assert result.endswith("_B23NWY7LT4")
        assert "LH00217" in result
        assert "0424" in result


class TestWriteRunInfoXml:
    """Tests for write_run_info_xml."""

    READ_NAME = "LH00217:424:B23NWY7LT4:1:1101:1001:1002 1:N:0:ACGT+TGCA"
    RUN_ID = "20260310_LM43899_0385_A12GGASZR5"
    FLOWCELL_ID = "A12GGASZR5"  # Last part of RUN_ID
    CYCLES_DICT = {"R1": 150, "I1": 8, "I2": 8, "R2": 150}

    def _write(self, tmp_path, **kwargs) -> str:
        """Write a RunInfo.xml with the given overrides and return its content."""
        xml_path = tmp_path / "RunInfo.xml"
        write_run_info_xml(
            xml_path,
            self.CYCLES_DICT,
            lane_count=1,
            surface_count=2,
            swath_count=4,
            tile_count=16,
            tile_names=["L001_1101", "L001_1102"],
            desc=self.READ_NAME,
            **kwargs,
        )
        return xml_path.read_text()

    def test_run_id_and_flowcell_from_read_name_by_default(self, tmp_path):
        """Without an override, the run id and <Flowcell> are derived from the read name and the date is the current time."""
        xml = self._write(tmp_path)
        assert "<Flowcell>B23NWY7LT4</Flowcell>" in xml
        assert 'B23NWY7LT4" Number="424"' in xml
        assert self.RUN_ID not in xml
        assert re.search(r"<Date>\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z</Date>", xml)

    def test_provided_flowcell_name_used_as_run_id(self, tmp_path):
        """A provided flowcell name (full run folder name) is used as-is as the run id, skipping the generation of a new one, while <Flowcell> holds only its last part and <Date> is pinned to its date segment at midnight."""
        xml = self._write(tmp_path, flowcell_name=self.RUN_ID)
        assert f'<Run Id="{self.RUN_ID}" Number="424">' in xml
        assert f"<Flowcell>{self.FLOWCELL_ID}</Flowcell>" in xml
        assert "<Date>2026-03-10T00:00:00Z</Date>" in xml
        assert "B23NWY7LT4" not in xml


class TestWriteRunParametersXml:
    """Tests for write_run_parameters_xml."""

    READ_NAME = "LH00217:424:B23NWY7LT4:1:1101:1001:1002 1:N:0:ACGT+TGCA"
    RUN_ID = "20260310_LM43899_0385_A12GGASZR5"

    def test_writes_run_context_and_template_structure(self, tmp_path):
        xml_path = tmp_path / "flowcell" / "RunParameters.xml"
        output_path = tmp_path / "flowcell"
        write_run_parameters_xml(
            xml_path,
            output_path,
            {"R1": 150, "I1": 8, "I2": 8, "R2": 150},
            self.READ_NAME,
            "NovaSeqXPlus",
            self.RUN_ID,
        )

        root = ET.parse(xml_path).getroot()
        assert root.findtext("OutputFolder") == str(output_path.resolve())
        assert root.findtext("InstrumentType") == "NovaSeqXPlus"
        assert root.findtext("InstrumentSerialNumber") == "LH00217"
        assert root.findtext("RunId") == self.RUN_ID
        assert root.findtext("RunCounter") == "424"
        assert [
            (read.get("ReadName"), read.get("Cycles"))
            for read in root.findall("./PlannedReads/Read")
        ] == [("Read1", "150"), ("Index1", "8"), ("Index2", "8"), ("Read2", "150")]
        assert [child.tag for child in root] == [
            "Application",
            "SystemSuiteVersion",
            "OutputFolder",
            "CustomPrimerSelections",
            "CloudUploadMode",
            "RunSetupMode",
            "SecondaryAnalysisMode",
            "InstrumentType",
            "InstrumentSerialNumber",
            "RunId",
            "ConsumableInfo",
            "PlannedReads",
            "SecondaryAnalysisInfo",
            "RunCounter",
            "RecipeName",
            "RecipeVersion",
            "ExperimentName",
            "PurgeReagentCartridge",
        ]


def test_parser_flowcell_id_option():
    """--flowcell-id is an optional argument defaulting to None."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    fastrewind_module.init_parser(subparsers)
    args = parser.parse_args(
        ["fastrewind", "--input-path", "in", "--output-path", "out"]
    )
    assert args.flowcell_id is None
    args = parser.parse_args(
        ["fastrewind", "--input-path", "in", "--flowcell-id", "A12GGASZR5"]
    )
    assert args.flowcell_id == "A12GGASZR5"


# ---------------------------------------------------------------------------
# preprocess_and_write_bcls
# ---------------------------------------------------------------------------


def test_preprocess_and_write_bcls_handles_odd_read_count(tmp_path, monkeypatch):
    """An odd number of reads still produces valid cbcl files for each cycle."""

    class FakeBaseCall:
        """Stand-in for a dnaio base call."""

        def __init__(self, base: str, quality: str):
            """Store the base and its quality score."""
            self.sequence = base
            self.qualities = quality

    class FakeRead:
        """Stand-in for a dnaio read (sequence plus qualities)."""

        def __init__(self, sequence: str, qualities: str):
            """Store the sequence and its quality scores."""
            self.sequence = sequence
            self.qualities = qualities

        def __len__(self):
            """Return the read length."""
            return len(self.sequence)

        def __getitem__(self, idx):
            """Return the base call at position idx."""
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
        """Yield the fake reads when the expected FASTQ path is opened."""
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
    (ntiles,) = struct.unpack_from("<I", data, 44)
    assert ntiles == 1
    _, _, _, comp = struct.unpack_from("<IIII", data, 48)
    return gzip.decompress(data[hsize : hsize + comp])


def test_preprocess_and_write_bcls_packs_clusters_little_endian(tmp_path, monkeypatch):
    """The first cluster of each byte goes in the low nibble and the odd tail
    is zero-padded in the high nibble (the order bcl-convert expects)."""

    class FakeBaseCall:
        """Stand-in for a dnaio base call."""

        def __init__(self, base: str, quality: str):
            """Store the base and its quality score."""
            self.sequence = base
            self.qualities = quality

    class FakeRead:
        """Stand-in for a dnaio read (sequence plus qualities)."""

        def __init__(self, sequence: str, qualities: str):
            """Store the sequence and its quality scores."""
            self.sequence = sequence
            self.qualities = qualities

        def __len__(self):
            """Return the read length."""
            return len(self.sequence)

        def __getitem__(self, idx):
            """Return the base call at position idx."""
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
        """Yield the fake reads when the expected FASTQ path is opened."""
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
        """Stand-in for a named dnaio read."""

        def __init__(self, name: str, sequence: str, qualities: str):
            """Store the read name, sequence and qualities."""
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        """Writer that discards everything written to it."""

        def write(self, read):
            """Discard the read."""
            return None

    class FakeWriterContext:
        """Context manager yielding a FakeWriter."""

        def __enter__(self):
            """Return a new FakeWriter."""
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    class FakeReaderContext:
        """Context manager yielding an iterator over the fake records."""

        def __init__(self, records):
            """Store the records to be iterated."""
            self.records = records

        def __enter__(self):
            """Return an iterator over the stored records."""
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
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
        """Return a writer or a reader context, recording the reader calls."""
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
        """Reader context that fails if it is ever entered."""

        def __enter__(self):
            """Fail: no reader may be opened without R1/R2 files."""
            raise AssertionError(
                "Reader should not be opened when no R1/R2 files exist"
            )

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
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
        """Return the reader context that fails on entry."""
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
    """Tests for write_index_metrics."""

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
    """Tests for build_index_records."""

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
        """Stand-in for a named dnaio read."""

        def __init__(self, name: str, sequence: str, qualities: str):
            """Store the read name, sequence and qualities."""
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        """Writer that discards everything written to it."""

        def write(self, read):
            """Discard the read."""
            return None

    class FakeWriterContext:
        """Context manager yielding a FakeWriter."""

        def __enter__(self):
            """Return a new FakeWriter."""
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    class FakeReaderContext:
        """Context manager yielding an iterator over the fake records."""

        def __init__(self, records):
            """Store the records to be iterated."""
            self.records = records

        def __enter__(self):
            """Return an iterator over the stored records."""
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    sample_name = "Sample_1_S1_L001_001"
    lane = "L001"
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)

    def make_read(tile: str, read: str) -> str:
        """Build an instrument-style read name for the given tile and read."""
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
        """Return a writer or a reader context depending on the open mode."""
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
        """Stand-in for a named dnaio read."""

        def __init__(self, name: str, sequence: str, qualities: str):
            """Store the read name, sequence and qualities."""
            self.name = name
            self.sequence = sequence
            self.qualities = qualities

    class FakeWriter:
        """Writer that discards everything written to it."""

        def write(self, read):
            """Discard the read."""
            return None

    class FakeWriterContext:
        """Context manager yielding a FakeWriter."""

        def __enter__(self):
            """Return a new FakeWriter."""
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    class FakeReaderContext:
        """Context manager yielding an iterator over the fake records."""

        def __init__(self, records):
            """Store the records to be iterated."""
            self.records = records

        def __enter__(self):
            """Return an iterator over the stored records."""
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
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
        """Return a writer or a reader context depending on the open mode."""
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
# build_index_lookup
# ---------------------------------------------------------------------------


class TestBuildIndexLookup:
    """Tests for build_index_lookup."""

    SAMPLE_SHEET = """\
[Header]
FileFormatVersion,2
RunName,DummyRun
InstrumentPlatform,NextSeq2000
[Reads]
Read1Cycles,151
Read2Cycles,151
Index1Cycles,10
Index2Cycles,10
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
1,S1,Sample_1,ACGT,TGCA,Proj_A,Y151;I10;I10;Y151
1,S2,Sample_2,AAAA,,Proj_B,Y151;I10;I10;Y151
2,S3,Sample_3,TTTT,,Proj_C,Y151;I10;I10;Y151
"""

    RECIPES = {
        "Sample_1_S1_L001_001": {
            "sample_name": "Sample_1",
            "sample_project": "Proj_A",
            "R1B": "",
            "R1": 151,
            "R1A": "",
            "I1": "ACGTNNNNNN",
            "I2": "TGCA",
            "R2B": "",
            "R2": 151,
            "R2A": "",
        },
        "Sample_2_S2_L001_001": {
            "sample_name": "Sample_2",
            "sample_project": "Proj_B",
            "R1B": "",
            "R1": 151,
            "R1A": "",
            "I1": "AAAA",
            "I2": "",
            "R2B": "",
            "R2": 151,
            "R2A": "",
        },
        "Sample_3_S3_L002_001": {
            "sample_name": "Sample_3",
            "sample_project": "Proj_C",
            "R1B": "",
            "R1": 151,
            "R1A": "",
            "I1": "TTTT",
            "I2": "",
            "R2B": "",
            "R2": 151,
            "R2A": "",
        },
    }

    def _sheet(self, tmp_path, content=SAMPLE_SHEET) -> pathlib.Path:
        """Write the given sheet content and return its path."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(content)
        return sheet

    def test_maps_informative_index_pairs_per_lane(self, tmp_path):
        """Each (lane, index1, index2) informative pair maps to the sample's recipe."""
        lookup = build_index_lookup(self._sheet(tmp_path), self.RECIPES)
        assert lookup[("L001", "ACGT", "TGCA")] is self.RECIPES["Sample_1_S1_L001_001"]
        assert lookup[("L001", "AAAA", "")] is self.RECIPES["Sample_2_S2_L001_001"]
        assert lookup[("L002", "TTTT", "")] is self.RECIPES["Sample_3_S3_L002_001"]
        assert len(lookup) == 3

    def test_rows_without_sample_name_are_skipped(self, tmp_path):
        """A blank data row is skipped but still consumes a sample index number,
        so the following sample keeps the group name of its physical row."""
        sheet = self._sheet(
            tmp_path,
            self.SAMPLE_SHEET
            + ",,,,,,\n2,S4,Sample_4,CCCC,,Proj_D,Y151;I10;I10;Y151\n",
        )
        recipes = dict(
            self.RECIPES,
            # Sample_4 sits on data row 5: the blank row above it is skipped
            # by the lookup but still consumes the index number 4
            Sample_4_S5_L002_001={
                "sample_name": "Sample_4",
                "sample_project": "Proj_D",
                "R1B": "",
                "R1": 151,
                "R1A": "",
                "I1": "CCCC",
                "I2": "",
                "R2B": "",
                "R2": 151,
                "R2A": "",
            },
        )
        lookup = build_index_lookup(sheet, recipes)
        assert lookup[("L002", "CCCC", "")] is recipes["Sample_4_S5_L002_001"]
        assert len(lookup) == 4

    def test_duplicate_informative_index_keeps_first_and_warns(self, tmp_path, caplog):
        """When two samples of the same lane share the same full index, the first
        recipe is kept and a warning is logged."""
        caplog.set_level(logging.WARNING)
        sheet = self._sheet(
            tmp_path,
            self.SAMPLE_SHEET.replace(
                "1,S2,Sample_2,AAAA,,Proj_B", "1,S2,Sample_2,ACGT,TGCA,Proj_B"
            ),
        )
        lookup = build_index_lookup(sheet, self.RECIPES)
        assert lookup[("L001", "ACGT", "TGCA")] is self.RECIPES["Sample_1_S1_L001_001"]
        assert len(lookup) == 2
        assert "Duplicate informative index" in caplog.text

    def test_sheet_without_index_columns_returns_empty(self, tmp_path):
        """A sheet whose data section has no index columns yields an empty lookup."""
        sheet = self._sheet(
            tmp_path,
            """\
[Header]
FileFormatVersion,2
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,Sample_Project,OverrideCycles
1,S1,Sample_1,Proj_A,Y151
""",
        )
        assert build_index_lookup(sheet, self.RECIPES) == {}

    def test_rows_missing_from_masks_table_are_skipped(self, tmp_path):
        """Sheet rows without a matching masks table recipe produce no entries."""
        sheet = self._sheet(tmp_path)
        assert build_index_lookup(sheet, {}) == {}
        partial = {"Sample_1_S1_L001_001": self.RECIPES["Sample_1_S1_L001_001"]}
        lookup = build_index_lookup(sheet, partial)
        assert len(lookup) == 1
        assert lookup[("L001", "ACGT", "TGCA")] is self.RECIPES["Sample_1_S1_L001_001"]


# ---------------------------------------------------------------------------
# synthesize_undetermined_recipe
# ---------------------------------------------------------------------------


class TestSynthesizeUndeterminedRecipe:
    """Tests for synthesize_undetermined_recipe."""

    RECIPE = {
        "sample_name": "Sample_1",
        "sample_project": "Proj_A",
        "R1B": "N" * 3,
        "R1": 148,
        "R1A": "N" * 5,
        "I1": "ACGTNNNNNN",
        "I2": "TGCA",
        "R2B": "N" * 2,
        "R2": 149,
        "R2A": "N",
    }

    def test_preserves_geometry_and_replaces_indexes_with_n(self):
        """The original recipe's geometry is preserved while the informative
        index cycles are replaced by N, one per original index cycle."""
        result = synthesize_undetermined_recipe(
            self.RECIPE, {"R1": 151, "I1": 10, "I2": 10, "R2": 151}
        )
        assert result["sample_name"] == "Undetermined"
        assert result["sample_project"] == "Undetermined"
        assert result["R1B"] == "N" * 3
        assert result["R1"] == 148
        assert result["R1A"] == "N" * 5
        assert result["I1"] == "N" * 10
        assert result["I2"] == "N" * 4
        assert result["R2B"] == "N" * 2
        assert result["R2"] == 149
        assert result["R2A"] == "N"

    def test_fallback_without_recipe_uses_full_length_n_indexes(self):
        """Without a resolved recipe, unpadded reads and full-length N indexes
        from the run's cycle counts are used."""
        result = synthesize_undetermined_recipe(
            None, {"R1": 151, "I1": 10, "I2": 10, "R2": 151}
        )
        assert result["R1B"] == ""
        assert result["R1A"] == ""
        assert result["R2B"] == ""
        assert result["R2A"] == ""
        assert result["R1"] == 0
        assert result["R2"] == 0
        assert result["I1"] == "N" * 10
        assert result["I2"] == "N" * 10

    def test_fallback_without_index_cycles(self):
        """A run without index cycles yields empty index segments."""
        result = synthesize_undetermined_recipe(None, {"R1": 151, "R2": 151})
        assert result["I1"] == ""
        assert result["I2"] == ""


# ---------------------------------------------------------------------------
# parse_fastq_groups undetermined groups
# ---------------------------------------------------------------------------


class _FakeRead:
    """Stand-in for a named dnaio read."""

    def __init__(self, name: str, sequence: str, qualities: str):
        """Store the read name, sequence and qualities."""
        self.name = name
        self.sequence = sequence
        self.qualities = qualities


def _patch_dnaio(monkeypatch, records: list, written: list) -> None:
    """Patch fastrewind_module.dnaio.open with a fake reader over `records`
    and a writer that records every written read into `written`."""

    class FakeWriter:
        """Writer that records every written read."""

        def write(self, read):
            """Record the read's name, sequence and qualities."""
            written.append((read.name, read.sequence, read.qualities))
            return None

    class FakeWriterContext:
        """Context manager yielding a FakeWriter."""

        def __enter__(self):
            """Return a new FakeWriter."""
            return FakeWriter()

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    class FakeReaderContext:
        """Context manager yielding an iterator over the fake records."""

        def __init__(self, records):
            """Store the records to be iterated."""
            self.records = records

        def __enter__(self):
            """Return an iterator over the stored records."""
            return iter(self.records)

        def __exit__(self, exc_type, exc_val, exc_tb):
            """Never suppress exceptions."""
            return False

    def fake_dnaio_open(*args, **kwargs):
        """Return a writer or a reader context depending on the open mode."""
        if "mode" in kwargs:
            return FakeWriterContext()
        return FakeReaderContext(records)

    monkeypatch.setattr(fastrewind_module.dnaio, "open", fake_dnaio_open)


INDEXED_RECIPES = {
    ("L001", "ACGT", "TGCA"): {
        "sample_name": "Sample_1",
        "sample_project": "Proj_A",
        "R1B": "",
        "R1": 2,
        "R1A": "",
        "I1": "ACGTNNNNNN",
        "I2": "TGCA",
        "R2B": "",
        "R2": 2,
        "R2A": "",
    },
    ("L001", "TTTT", "CCCC"): {
        "sample_name": "Sample_2",
        "sample_project": "Proj_B",
        "R1B": "",
        "R1": 2,
        "R1A": "",
        "I1": "TTTT",
        "I2": "CCCC",
        "R2B": "",
        "R2": 2,
        "R2A": "",
    },
}


def test_parse_fastq_groups_undetermined_resolves_recipe_from_index(
    tmp_path, monkeypatch, caplog
):
    """Undetermined paired-end reads are written with the N-padded index of
    their original sample's recipe, resolved per read from the index carried in
    the read name, and are not counted for the InterOp index metrics."""
    caplog.set_level(logging.WARNING)
    written = []
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)
    read_a = "INST:1:FC:1:1101:1001:1002 1:N:0:ACGT+TGCA"
    read_b = "INST:1:FC:1:1102:1003:1004 1:N:0:TTTT+CCCC"
    records = [
        (
            _FakeRead(read_a, "AC", "II"),
            _FakeRead(read_a.replace(" 1:", " 2:"), "TA", "II"),
        ),
        (
            _FakeRead(read_b, "GT", "II"),
            _FakeRead(read_b.replace(" 1:", " 2:"), "CG", "II"),
        ),
    ]
    _patch_dnaio(monkeypatch, records, written)

    index_counts = {}
    result = parse_fastq_groups(
        lane="L001",
        name="Undetermined_S0_L001_001",
        filenames=[
            tmp_path / "Undetermined_S0_L001_R1_001.fastq.gz",
            tmp_path / "Undetermined_S0_L001_R2_001.fastq.gz",
        ],
        masks_table={},
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
        index_counts=index_counts,
        index_lookup=INDEXED_RECIPES,
        cycles_dict={"R1": 2, "I1": 4, "I2": 4, "R2": 2},
    )

    # Undetermined clusters are excluded from the InterOp index metrics
    assert index_counts == {}
    # The N-padded indexes of each original sample's recipe are used:
    # read_a comes from a sample with I1="ACGTNNNNNN" (10 cycles) + I2="TGCA"
    # (4 cycles), read_b from one with I1="TTTT" (4) + I2="CCCC" (4)
    by_name = {name: (seq, qual) for name, seq, qual in written}
    assert by_name[read_a][0] == "AC" + "N" * 10 + "N" * 4 + "TA"
    assert by_name[read_b][0] == "GT" + "N" * 4 + "N" * 4 + "CG"
    # Qualities are all 'I', one per cycle of the synthesised layout
    assert by_name[read_a][1] == "I" * 18
    assert by_name[read_b][1] == "I" * 12
    # Positions are returned like for regular sample groups
    assert (1001, 1002) in result["1101"]
    assert (1003, 1004) in result["1102"]
    assert caplog.text == ""


def test_parse_fastq_groups_undetermined_falls_back_for_unknown_index(
    tmp_path, monkeypatch, caplog
):
    """An undetermined read whose index matches no sample of the lane is
    written with the fallback recipe and a warning is logged."""
    caplog.set_level(logging.WARNING)
    written = []
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)
    known = "INST:1:FC:1:1101:1001:1002 1:N:0:ACGT+TGCA"
    unknown = "INST:1:FC:1:1101:1003:1004 1:N:0:GATC+TGCA"
    records = [
        (
            _FakeRead(known, "AC", "II"),
            _FakeRead(known.replace(" 1:", " 2:"), "TA", "II"),
        ),
        (
            _FakeRead(unknown, "GT", "II"),
            _FakeRead(unknown.replace(" 1:", " 2:"), "CG", "II"),
        ),
    ]
    _patch_dnaio(monkeypatch, records, written)

    parse_fastq_groups(
        lane="L001",
        name="Undetermined_S0_L001_001",
        filenames=[
            tmp_path / "Undetermined_S0_L001_R1_001.fastq.gz",
            tmp_path / "Undetermined_S0_L001_R2_001.fastq.gz",
        ],
        masks_table={},
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
        index_lookup=INDEXED_RECIPES,
        cycles_dict={"R1": 2, "I1": 4, "I2": 4, "R2": 2},
    )

    sequences = {name: seq for name, seq, _ in written}
    # Known index: original recipe geometry (10 + 4 N index cycles)
    assert sequences[known] == "AC" + "N" * 10 + "N" * 4 + "TA"
    # Unknown index: fallback recipe (full-length 4 + 4 N index cycles)
    assert sequences[unknown] == "GT" + "N" * 4 + "N" * 4 + "CG"
    assert "carries index 'GATC+TGCA'" in caplog.text


def test_parse_fastq_groups_undetermined_single_end(tmp_path, monkeypatch, caplog):
    """Undetermined single-end reads get the N-padded I1 of their original
    sample appended to the read, like regular single-end reads."""
    caplog.set_level(logging.WARNING)
    written = []
    tempdir = tmp_path / "lane_tmp"
    tempdir.mkdir(parents=True)
    read_name = "INST:1:FC:1:1101:1001:1002 1:N:0:ACGT"
    records = [_FakeRead(read_name, "AC", "II")]
    _patch_dnaio(monkeypatch, records, written)

    index_counts = {}
    result = parse_fastq_groups(
        lane="L001",
        name="Undetermined_S0_L001_001",
        filenames=[tmp_path / "Undetermined_S0_L001_R1_001.fastq.gz"],
        masks_table={},
        tempdir=tempdir,
        threads=0,
        instrument="NovaSeqXPlus",
        index_counts=index_counts,
        index_lookup={
            ("L001", "ACGT", ""): {
                "sample_name": "Sample_1",
                "sample_project": "Proj_A",
                "R1B": "",
                "R1": 2,
                "R1A": "",
                "I1": "ACGTNNNNNN",
                "I2": "",
                "R2B": "",
                "R2": 0,
                "R2A": "",
            }
        },
        cycles_dict={"R1": 2, "I1": 10},
    )

    assert index_counts == {}
    # Only the N-padded I1 is appended in single-end mode
    assert [seq for _, seq, _ in written] == ["AC" + "N" * 10]
    assert [qual for _, _, qual in written] == ["I" * 12]
    assert (1001, 1002) in result["1101"]
    assert caplog.text == ""


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
        """Create a minimal run folder the interop library can parse."""
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


# ---------------------------------------------------------------------------
# blabber -> fastrewind round trip with taint count
# ---------------------------------------------------------------------------


def _count_reads(fastq: pathlib.Path) -> int:
    """Count the number of records in a gzip-compressed FASTQ file."""
    with gzip.open(fastq, "rt") as fh:
        return sum(1 for line in fh if line.startswith("@"))


def _cbcl_cluster_count(cbcl_path: pathlib.Path) -> int:
    """Sum the per-tile cluster counts stored in a cbcl file header."""
    data = cbcl_path.read_bytes()
    (ntiles,) = struct.unpack_from("<I", data, 44)
    clusters = 0
    offset = 48
    for _ in range(ntiles):
        _, tile_clusters, _, _ = struct.unpack_from("<IIII", data, offset)
        clusters += tile_clusters
        offset += 16
    return clusters


def test_blabber_taint_count_round_trip_through_fastrewind(tmp_path, caplog):
    """Undetermined files produced by blabber --taint-count are processed by
    fastrewind: every undetermined read becomes a real BCL cluster on top of
    the sample clusters, no group is skipped, and no InterOp index metrics
    record is produced for the undetermined clusters."""
    caplog.set_level(logging.WARNING)
    sample_sheet = (
        pathlib.Path(__file__).resolve().parents[1]
        / "assets"
        / "SampleSheet_SingleIndex.csv"
    )
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(sample_sheet.read_text())
    flowcell_id = "20260310_LM43898_0388_D74LKJHGF8"

    blabber_main(
        argparse.Namespace(
            verbose=False,
            quiet=True,
            random_seed=42,
            alphabet="ACGT",
            format="fastq",
            sample_sheet=sheet,
            output=tmp_path,
            flowcell_id=flowcell_id,
            seq_number=4,
            seq_length=151,
            taint=True,
            taint_count=2,
        )
    )

    output_path = tmp_path / "bcl_out"
    fastrewind_module.main(
        argparse.Namespace(
            verbose=False,
            quiet=True,
            input_path=tmp_path,
            output_path=output_path,
            sample_sheet=sheet,
            total_cycles=None,
            interop_dir=None,
            instrument="NovaSeqXPlus",
            flowcell_id=flowcell_id,
            threads=0,
            force=False,
        )
    )

    # No FASTQ group was skipped and every undetermined read found its
    # original sample in the index lookup
    assert "Skipping FASTQ group" not in caplog.text
    assert "matches no sample in lane" not in caplog.text

    demux_path = tmp_path / flowcell_id / "Demultiplexing"
    for lane, sample_count in (("L001", 3), ("L002", 5), ("L003", 4), ("L004", 5)):
        lane_files = sorted(demux_path.rglob(f"*_{lane}_R1_001.fastq.gz"))
        sample_reads = sum(
            _count_reads(f) for f in lane_files if "Undetermined" not in f.name
        )
        undetermined = [f for f in lane_files if "Undetermined" in f.name]
        assert len(undetermined) == 1
        assert _count_reads(undetermined[0]) == 2
        assert sample_reads == sample_count * 4

        clusters = 0
        for surface in ("1", "2"):
            cbcl = (
                output_path
                / f"Data/Intensities/BaseCalls/{lane}/C1.1/{lane}_{surface}.cbcl"
            )
            if cbcl.exists():
                clusters += _cbcl_cluster_count(cbcl)
        # The BCL holds every sample cluster plus the two undetermined ones
        assert clusters == sample_reads + 2

    records = _unpack_index_metrics(
        (output_path / "InterOp" / "IndexMetricsOut.bin").read_bytes(), 2
    )
    assert records
    assert all(record[5] != "Undetermined" for record in records)


# ---------------------------------------------------------------------------
# fit_read_to_recipe / encode_cluster_bits
# ---------------------------------------------------------------------------


def test_fit_read_pads_short_and_rejects_long():
    """A short read is N-padded to the recipe length; a long one raises."""
    from biomate.fastrewind.fastrewind import fit_read_to_recipe

    read = dnaio.SequenceRecord("r", "ACGT", "IIII")
    fit_read_to_recipe(read, 6, "R1")
    assert (read.sequence, read.qualities) == ("ACGTNN", "IIII##")
    fit_read_to_recipe(read, 0, "R1")  # unknown recipe length: untouched
    assert read.sequence == "ACGTNN"
    with pytest.raises(ValueError, match="allows 5 cycles"):
        fit_read_to_recipe(read, 5, "R1")


def test_encode_cluster_bits_bins_quality_and_zeroes_no_calls():
    """Off-bin qualities take the nearest lower bin; N is an all-zero nibble."""
    from biomate.fastrewind.fastrewind import encode_cluster_bits

    qmap = {0: "00", 9: "01", 24: "10", 40: "11"}
    assert encode_cluster_bits(dnaio.SequenceRecord("r", "C", "F"), qmap) == "1001"
    assert encode_cluster_bits(dnaio.SequenceRecord("r", "C", "-"), qmap) == "0101"
    assert encode_cluster_bits(dnaio.SequenceRecord("r", "C", "#"), qmap) == "0001"
    assert encode_cluster_bits(dnaio.SequenceRecord("r", "N", "I"), qmap) == "0000"
