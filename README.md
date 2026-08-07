# BioMate

A collection of scripts and utility tools for everyday processing of biological data.

For detailed tool documentation, see [docs/README.md](docs/README.md).

## Installation

### Requirements
- Python 3.12 or later
- [uv](https://docs.astral.sh/uv/) (recommended) or pip

### Recommended: Using uv

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
│   └── web_interface.py  # Tornado web UI
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
