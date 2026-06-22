"""Tests for the dirstruct module."""

import argparse

import pytest

from biomate.dirstruct.dirstruct import validate_args, main


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def test_extract_valid_source_path(self, tmp_path):
        """validate_args accepts an existing directory for the extract command."""
        src = tmp_path / "src"
        src.mkdir()
        args = argparse.Namespace(command="extract", source_path=src, output_file=None)
        assert validate_args(args) is args

    def test_extract_nonexistent_source_raises(self, tmp_path):
        """A non-existent source path for extract raises ArgumentTypeError."""
        args = argparse.Namespace(
            command="extract",
            source_path=tmp_path / "missing",
            output_file=None,
        )
        with pytest.raises(
            argparse.ArgumentTypeError, match="[Dd]oes not exist|not a directory"
        ):
            validate_args(args)

    def test_create_valid_source_file(self, tmp_path):
        """validate_args accepts an existing file for the create command."""
        f = tmp_path / "struct.txt"
        f.write_text("")
        args = argparse.Namespace(command="create", source_file=f, output_path=tmp_path)
        assert validate_args(args) is args

    def test_create_nonexistent_source_file_raises(self, tmp_path):
        """A non-existent source file for create raises ArgumentTypeError."""
        args = argparse.Namespace(
            command="create",
            source_file=tmp_path / "missing.txt",
            output_path=tmp_path,
        )
        with pytest.raises(
            argparse.ArgumentTypeError, match="[Dd]oes not exist|not a file"
        ):
            validate_args(args)


# ---------------------------------------------------------------------------
# main — extract command
# ---------------------------------------------------------------------------


def _extract_args(src, output_file=None, no_tags=False, verbose=False, quiet=False):
    return argparse.Namespace(
        command="extract",
        source_path=src,
        output_file=output_file,
        no_tags=no_tags,
        verbose=verbose,
        quiet=quiet,
    )


class TestExtractCommand:
    def test_extract_writes_to_file(self, tmp_path):
        """Extracted directory tree is written to the output file."""
        src = tmp_path / "src"
        src.mkdir()
        (src / "subdir").mkdir()
        (src / "file.txt").write_text("hello")

        out = tmp_path / "struct.txt"
        main(_extract_args(src, output_file=out, quiet=True))

        content = out.read_text()
        assert "file.txt" in content
        assert "subdir" in content

    def test_extract_writes_tags_by_default(self, tmp_path):
        """Output file contains dir/file type tags by default."""
        src = tmp_path / "src"
        src.mkdir()
        (src / "a.txt").write_text("")

        out = tmp_path / "struct.txt"
        main(_extract_args(src, output_file=out, quiet=True))

        content = out.read_text()
        assert "file" in content
        assert "dir" in content

    def test_extract_no_tags_omits_tags(self, tmp_path):
        """With --no-tags the output contains no tab-separated type labels."""
        src = tmp_path / "src"
        src.mkdir()
        (src / "a.txt").write_text("")

        out = tmp_path / "struct.txt"
        main(_extract_args(src, output_file=out, no_tags=True, quiet=True))

        content = out.read_text()
        assert "\tfile" not in content
        assert "\tdir" not in content

    def test_extract_creates_parent_directory(self, tmp_path):
        """Parent directories of the output file are created if they do not exist."""
        src = tmp_path / "src"
        src.mkdir()
        out = tmp_path / "nested" / "struct.txt"
        main(_extract_args(src, output_file=out, quiet=True))
        assert out.exists()


# ---------------------------------------------------------------------------
# main — create command
# ---------------------------------------------------------------------------


def _create_args(source_file, output_path, verbose=False, quiet=False):
    return argparse.Namespace(
        command="create",
        source_file=source_file,
        output_path=output_path,
        output_file=None,  # required by main() before branching on command
        verbose=verbose,
        quiet=quiet,
    )


class TestCreateCommand:
    def test_create_directories_from_file(self, tmp_path):
        """Directories listed in the struct file are created at the destination."""
        struct_file = tmp_path / "struct.txt"
        struct_file.write_text("mydir\tdir\nmydir/subdir\tdir\n")
        dest = tmp_path / "dest"
        dest.mkdir()
        main(_create_args(struct_file, dest, quiet=True))
        assert (dest / "mydir").is_dir()
        assert (dest / "mydir" / "subdir").is_dir()

    def test_create_files_from_struct(self, tmp_path):
        """Files listed in the struct file are created at the destination."""
        struct_file = tmp_path / "struct.txt"
        struct_file.write_text("mydir\tdir\nmydir/hello.txt\tfile\n")
        dest = tmp_path / "dest"
        dest.mkdir()
        main(_create_args(struct_file, dest, quiet=True))
        assert (dest / "mydir" / "hello.txt").is_file()

    def test_round_trip_extract_then_create(self, tmp_path):
        """Extract a real directory structure, then recreate it elsewhere."""
        src = tmp_path / "original"
        src.mkdir()
        (src / "a").mkdir()
        (src / "a" / "b").mkdir()
        (src / "a" / "file.txt").write_text("hi")

        struct_file = tmp_path / "struct.txt"
        main(_extract_args(src, output_file=struct_file, quiet=True))

        dest = tmp_path / "restored"
        dest.mkdir()
        main(_create_args(struct_file, dest, quiet=True))

        assert (dest / src.name / "a").is_dir()
        assert (dest / src.name / "a" / "b").is_dir()
        assert (dest / src.name / "a" / "file.txt").is_file()
