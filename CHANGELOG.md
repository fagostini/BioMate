## 0.3.0 (2026-05-11)

### Feat

- Add `release` recipe to package the test case

### Fix

- Fix bash script to exclude index files from checks

### Refactor

- Refactor `parse_sequence_mask`

## 0.2.0 (2026-04-01)

### BREAKING CHANGE

- Temporarily commented the `taint` procedure

### Feat

- Add `commitizen` as dev dependency
- **Makefile**: Add mixed indexes sample sheet and recipes for it
- Add individual FastQ comparison of original and demultiplexed files
- **fastrewind**: Add support for 'N' characters in cycles masks
- **fastrewind**: Implement `masks_table` to track and correctly handle each sample
- Add `temp` directory to `.gitignore`
- Add test recipes to `Makefile`
- **blabber**: Add `--flowcell-id` option to avoid random names
- **fastrewind**: Add copy of SampleSheet.csv to destintion BCL directory
- **blabber**: Improve SampleSheet v2 parsing and organisation
- **blabber**: Remove complete _randomness_ in tile positions generation
- **fastrewind**: Implement purge of previous data when `--force` is provided
- **fastrewind**: Add parsing of sequence masks
- Implemented support for design masks and samplesheet v2 in `blabber` module
- Add `Makefile` and basic recipes
- Add support and config file for `pre-commit`
- Add `pre-commit` config file
- Implement BCL content generation and dump to file
- Include `FastRewind` module in `Biomate`
- Implement working version of `FastRewind` module (still missing generation of some files)
- Implement patterns from file processing in `index` module
- Add tags to `dirstruct.py` output
- Initialise `mkdocs` folders and files
- Document `blabber`
- Add `--taint` flag in `blabber` module
- Implement 'tainting' of Undetermined files
- Implement `blabber` module
- Add module as `BioMate` submodule
- Initial commit
- Add `dirstruct` and `index` modules
- Initialise repository with `uv`
- Add `index` dependencies

### Fix

- **fastrewind**: Fix regex
- Remove `set -o pipefail` from the comparison script to view intermediate outputs
- **fastrewind**: Fix labelling when the sample sheet contains more than 10 samples
- **fastrewind**: Improve handling of missing indexes
- **blabber**: Small fix to index parser
- **fastrewind**: Improve parsing of total number of cycles
- **fastrewind**: Fix handling of lanes with different number of samples
- **fastrewind**: Fix data type and value in BCL header
- **fastrewind**: Fix order of little endian bits in clusters encoding
- **fastrewind**: Fix tile pattern to match placeholder values
- Add _test_ folder to _.gitiignore_
- **fastrewind**: Fix generation of `*.filter` files
- Remove debug line
- Fix README path in `pyproject.toml`
- Add `notebooks` folder to `.gitignore`
- Add outputs extensions to `.gitignore`
- Add verbosity control to `dirstruct` module
- Fix arguments of `index.py` module
- Change default parameters of `index` module
- Fix `validate_args` return value in `dirstruct.py`
- Fix `setup_logging` arguments in `dirstruct.py`
- Fix `output_file` check in `dirstruct.py`
- Fix logging level setting
- Remove redundant line in `indexing` submodule
- Add `.DS_store` to `.gitignore`

### Refactor

- **fastrewind**: Drop support for SampleSheet version 1
- **fastrewind**: Improve parsing of total cycles from SampleSheet
- **fastrewind**: Refactor `parse_sequence_mask` function so that it can handle correctly mixed masks
- **blabber**: Refactor sample sheet parsing and processing functions
- **blabber**: Streamline `parse_sequence_mask` function
- **fastrewind**: Streamline `parse_sequence_mask` function
- **fastrewind**:  Code cleanup
- Code clean-up
- Re-implement logic and procedure for input parsing
- Rename `indexing` module to `index`
- Rename project from `BioFun` to `BioMate` and fix files accordingly

### Perf

- Refactor code to improve performance
