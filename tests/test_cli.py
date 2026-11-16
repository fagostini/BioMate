"""Tests for the top-level biomate command line."""

import logging

import pytest

from biomate import main


def _run(monkeypatch, *argv):
    """Run the biomate CLI with the given arguments."""
    monkeypatch.setattr("sys.argv", ["biomate", *argv])
    main()


def test_invalid_input_exits_with_message_not_traceback(monkeypatch, capsys, tmp_path):
    """Validation errors are reported as a usage error (exit 2)."""
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, "dirstruct", "create", "--source-file", str(tmp_path / "x"))
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "error: The source file" in err
    assert "Traceback" not in err


@pytest.mark.parametrize(
    ("flag", "level"), [("--verbose", logging.DEBUG), ("--quiet", logging.ERROR)]
)
def test_verbosity_flags_after_subcommand(monkeypatch, tmp_path, flag, level):
    """--verbose/--quiet are accepted after a nested sub-command."""
    struct = tmp_path / "s.txt"
    struct.write_text("d\tdir\n")
    _run(
        monkeypatch,
        "dirstruct",
        "create",
        "--source-file",
        str(struct),
        "--output-path",
        str(tmp_path / "out"),
        flag,
    )
    assert (tmp_path / "out" / "d").is_dir()
    assert logging.getLogger().level == level
