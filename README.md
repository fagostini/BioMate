# BioMate

A collection of scripts and utility tools for everyday processing of biological data.

For detailed tool documentation, see [docs/README.md](docs/README.md).

## Installation

### Requirements
- Python 3.12 or later
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

> _**Note:** Running the demultiplexing test recipes additionally requires the Illumina [BCL-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert/downloads.html) executable, which is not bundled with this repository because its redistribution is prohibited. Download it from the [Illumina website](https://support.illumina.com/sequencing/sequencing_software/bcl-convert/downloads.html) and place it in the `assets/` directory. See [docs/README.md](docs/README.md#external-dependencies) for details._

### Install as a command-line tool

To only use `biomate` (no development setup), install it as an isolated tool:

```bash
uv tool install git+https://github.com/fagostini/BioMate

# With InterOp support (nspector --run-folder)
uv tool install "biomate[interop] @ git+https://github.com/fagostini/BioMate"

biomate --help
```

### Development setup with uv

```bash
# Clone repository
git clone https://github.com/fagostini/BioMate.git
cd BioMate

# Install all dependencies + dev tools
uv sync --all-extras --all-groups

# Run commands with isolated environment
uv run biomate --help
uv run pytest tests/
```

### Alternative: Using pip

```bash
# Generate requirements.txt from lock file
uv export --frozen --output-file=requirements.txt

# Install dependencies
pip install -r requirements.txt

# Install package in development mode
pip install -e .

# Run commands
biomate --help
pytest tests/
```

### Rust engine (optional)

`nspector` and `index --engine rust` can use `biomate_kernels`, a compiled Rust extension in `rust/biomate_kernels/`. It is not on PyPI and is not built by `uv sync` or `pip install`. The Docker and Singularity images build and include it. Without it, `nspector --engine auto` (the default) uses the Python engine, which gives identical results but is about 20× slower on large runs (one hour against three minutes for 277 M reads).

Requirements:
- a Rust toolchain (`cargo`, `rustc`), e.g. from [rustup](https://rustup.rs/)
- [maturin](https://www.maturin.rs/) `>=1.5,<2.0`, installed in the same environment as `biomate`
- a clone of this repository (the kernel is built from source)

In a development setup (`make sync` done, `cargo` on `PATH`), build it into `.venv` and check it:

```bash
make kernels
```

To build it into another environment that runs `biomate` by hand:

```bash
uv pip install "maturin>=1.5,<2.0"
cd rust/biomate_kernels
VIRTUAL_ENV=../../.venv maturin develop --release
```

- Always pass `--release`. A debug build is many times slower.
- In a conda shell, maturin refuses to build when both `CONDA_PREFIX` and `VIRTUAL_ENV` are set. Prefix the command with `env -u CONDA_PREFIX`.
- `uv sync` removes packages the lock file does not list, including this one. Use `uv sync --inexact` to keep it, or rebuild after each sync.
- Rebuild after pulling changes to `rust/biomate_kernels/`. When the installed kernel is from another BioMate version, `nspector --engine auto` warns and uses Python, and `--engine rust` stops with an error.

Check the install:

```bash
uv run python -c "import biomate_kernels; print(biomate_kernels.SCAN_FORMAT)"
```

## Dependency Management

BioMate uses a **single source of truth** for dependency management:

### Architecture

- **`pyproject.toml`** - Source of truth
  - Defines all dependencies with version constraints
  - Developers edit this file to change dependencies
  - What gets published to PyPI

- **`uv.lock`** - Reproducible builds
  - Machine-readable lock file with exact resolved versions
  - Committed to repository
  - Ensures reproducible builds across machines

- **`requirements.txt`** - Generated artifact
  - Auto-generated from `uv.lock` for pip-only environments
  - NOT committed (regenerated as needed)
  - Contains all transitive dependencies with hashes

### Updating Dependencies

To update dependencies:

1. Edit `pyproject.toml` to change version constraints
2. Run `uv lock` to update lock file
3. Regenerate `requirements.txt` if needed:
   ```bash
   uv export --frozen --output-file=requirements.txt
   ```
4. Commit: `pyproject.toml` and `uv.lock`

### For CI/CD

Use `uv` in CI pipelines (recommended):

```yaml
- run: uv sync --all-extras --all-groups
- run: uv run pytest tests/
```

Or use pip with generated requirements:

```yaml
- run: uv export --frozen --output-file=requirements.txt
- run: pip install -r requirements.txt
- run: pytest tests/
```

## Development

### Run Tests

```bash
uv run pytest tests/ -v
```

### Run Linting/Formatting

```bash
uv run ruff check src/
uv run ruff format src/
```

### Build Documentation

```bash
uv run mkdocs serve
```

## Quick Start

```bash
# Install everything
uv sync --all-extras --all-groups

# See available tools
uv run biomate --help

# Example: Generate random sequences
uv run biomate blabber --seq-number 100 --seq-length 50

# Run tests
uv run pytest tests/ -v
```

## Project Structure

```
BioMate/
├── src/biomate/           # Main package source code
│   ├── blabber/          # Sequence generator
│   ├── dirstruct/        # Directory structure tool
│   ├── fastrewind/       # FASTQ file conversion
│   ├── index/            # Indexing tool
│   ├── nspector/         # FASTQ inspector
│   ├── strainer/         # Index analysis
│   ├── triage/           # Undetermined reads triage
│   ├── waif/             # Always-undetermined barcodes
│   ├── design/, static/  # Web UI templates and assets
│   ├── headers.py        # Read header index counting
│   ├── samplesheet.py    # Sample sheet parsing
│   └── web_interface.py  # Tornado web UI
├── rust/biomate_kernels/ # Optional Rust kernels
├── tests/                # Test suite
├── docs/                 # Documentation
├── pyproject.toml        # Dependency definitions
├── uv.lock              # Frozen dependencies (committed)
└── README.md            # This file
```

## License

See [LICENSE](LICENSE) file for details.

## Contributing

See existing documentation for code standards and testing requirements.

---

**Note:** This project requires Python 3.12+. For development, using `uv` is strongly recommended for faster, more reliable dependency management.
