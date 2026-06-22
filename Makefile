# Makefile for the project
include Makefile.help

# Makefile containing the project's variables
PROJECT_NAME := BioMate
PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
UV ?= uv
BCLCONVERT ?= bcl-convert
SHELL=/bin/bash

# Test assets and variables
NISS = assets/SampleSheet_NoIndex.csv
SISS = assets/SampleSheet_SingleIndex.csv
DISS = assets/SampleSheet_DualIndex.csv
MISS = assets/SampleSheet_MixedIndexes.csv

TEMP_DIR = temp
FLOWCELL_ID = 20260310_LM43899_0385_A12GGASZR5
SEQ_NUM = 10

RELEASE_DATA = DemuxCasesData.tar.gz
TEST_PATH ?= tests
PYTEST_ARGS ?=

.PHONY: check-uv
check-uv:
	@command -v $(UV) >/dev/null 2>&1 || { \
		echo "Error: '$(UV)' not found in PATH."; \
		echo "Install uv (https://docs.astral.sh/uv/) or set UV=/full/path/to/uv."; \
		exit 1; \
	}

.PHONY: install-uv
## Install uv using the official installer script
install-uv:
	@curl -LsSf https://astral.sh/uv/install.sh | sh

.PHONY: check-bclconvert
check-bclconvert:
	@command -v $(BCLCONVERT) >/dev/null 2>&1 || { \
		echo "Error: '$(BCLCONVERT)' not found in PATH."; \
		echo "Install bcl-convert or set BCLCONVERT=/full/path/to/bcl-convert."; \
		exit 1; \
	}


.PHONY: sync
## Sync the project's dependencies with the environment
sync: check-uv
	@$(UV) sync --all-extras --all-groups


uv.lock: check-uv
	@$(UV) lock


.PHONY: lock
## Create a lockfile for the project's dependencies
lock: uv.lock


.PHONY: upgrade
## Upgrade the project's dependencies
upgrade: check-uv
	@$(UV) lock --upgrade


.PHONY: build
## Build the project into distribution archives
build: check-uv
	@$(UV) build


.venv/bin/activate: check-uv
	@$(UV) venv


.PHONY: venv
## Create a new virtual environment
venv: .venv/bin/activate

.PHONY: singularity
## Create the Singularity image for the project
singularity: remove_singularity singularity.sif

.PHONY: remove_singularity

singularity.sif: singularity.def
	@singularity build --fakeroot --fix-perms --force --bind $(PROJECT_DIR):/mnt $@ $< > singularity_build.out 2> singularity_build.err && echo "Singularity image built successfully!" || echo "Error building Singularity image! Check singularity_build.err for details."

remove_singularity:
	@if [[ -f singularity.sif ]]; then rm singularity.sif; fi


requirements.txt: check-uv uv.lock
	@$(UV) export --output-file requirements.txt


.PHONY: export
## Export the project's dependencies to a requirements.txt file
export: requirements.txt


.PHONY: interrogate
## Run interrogate in verbose mode to check for code quality
interrogate: check-uv
	@$(UV)x interrogate --verbose


.PHONY: test
## Run the pytest test suite
test: check-uv
	@$(UV) run pytest $(TEST_PATH) $(PYTEST_ARGS)


.PHONY: test_verbose
## Run the pytest suite in verbose mode
test_verbose: check-uv
	@$(UV) run pytest -v --tb=short $(TEST_PATH) $(PYTEST_ARGS)


.PHONY: test_failed
## Re-run only tests that failed in the previous pytest run
test_failed: check-uv
	@$(UV) run pytest --lf -v --tb=short $(TEST_PATH) $(PYTEST_ARGS)

.PHONY: deploy
## Deploy documentation on GitHub
deploy: check-uv
	@$(UV) run mkdocs gh-deploy


.PHONY: release
## Create a new release
release: $(RELEASE_DATA)


$(RELEASE_DATA): test_mix
	mv $(TEMP_DIR) $(basename $(basename $(RELEASE_DATA)))
	tar -czf $(RELEASE_DATA) $(basename $(basename $(RELEASE_DATA)))
	@rm -r $(basename $(basename $(RELEASE_DATA)))

$(TEMP_DIR):
	@mkdir -p $(TEMP_DIR)

copy_none: $(TEMP_DIR) $(NISS) clean_samplesheet
	@cp $(NISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (No Index) to working folder." || { echo "Error copying SampleSheet to working folder!"; exit 1; }

copy_single: $(TEMP_DIR) $(SISS) clean_samplesheet
	@cp $(SISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Single Index) to working folder." || { echo "Error copying SampleSheet to working folder!"; exit 1; }

copy_dual: $(TEMP_DIR) $(DISS) clean_samplesheet
	@cp $(DISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Dual Indexes) to working folder." || { echo "Error copying SampleSheet to working folder!"; exit 1; }

copy_mix: $(TEMP_DIR) $(MISS) clean_samplesheet
	@cp $(MISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Mixed Indexes) to working folder." || { echo "Error copying SampleSheet to working folder!"; exit 1; }

$(TEMP_DIR)/$(FLOWCELL_ID): check-uv $(TEMP_DIR) clean_flowcell
	@$(UV) run biomate --verbose blabber --format fastq --sample-sheet $(TEMP_DIR)/SampleSheet.csv --seq-number $(SEQ_NUM) --output $(TEMP_DIR) --flowcell-id $(FLOWCELL_ID) > $(TEMP_DIR)/blabber.out 2> $(TEMP_DIR)/blabber.err && echo "Blabber module executed successfully!" || { echo "Error executing Blabber module!"; exit 1; }

$(TEMP_DIR)/Data $(TEMP_DIR)/RunInfo.xml: check-uv $(TEMP_DIR)/$(FLOWCELL_ID) clean_data
	@$(UV) run biomate --verbose fastrewind --input-path $(TEMP_DIR) --output-path $(TEMP_DIR) --threads 8 > $(TEMP_DIR)/fastrewind.out 2> $(TEMP_DIR)/fastrewind.err  && echo "Fastrewind module executed successfully!" || { echo "Error executing Fastrewind module!"; exit 1; }

validate_samplesheet: check-bclconvert $(TEMP_DIR)/RunInfo.xml
	@$(BCLCONVERT) --output-directory $(TEMP_DIR)/Demultiplexing --bcl-input-directory $(TEMP_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true > $(TEMP_DIR)/bcl-validate.out 2> $(TEMP_DIR)/bcl-validate.err && echo "BCL-convert SampleSheet validation was successful!" || { echo "Error executing BCL-convert SampleSheet validation!"; exit 1; }

$(TEMP_DIR)/Demultiplexing: check-bclconvert $(TEMP_DIR)/RunInfo.xml clean_demux
	@$(BCLCONVERT) --output-directory $(TEMP_DIR)/Demultiplexing --bcl-input-directory $(TEMP_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true > $(TEMP_DIR)/bcl-convert.out 2> $(TEMP_DIR)/bcl-convert.err && echo "BCL-convert executed successfully!" || { echo "Error executing BCL-convert!"; exit 1; }


.PHONY: report_results
report_results: $(TEMP_DIR)/Demultiplexing results_message
	@find $(TEMP_DIR)/Demultiplexing -name "*.fastq.gz" | grep -v -e "Undetermined" -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true
	@find $(TEMP_DIR)/Demultiplexing -name "*.fastq.gz" | grep "Undetermined" | grep -v -e "_I1_" -e "_I2_" | xargs -r zgrep -c ^@ || true


compare_results: $(TEMP_DIR)/Demultiplexing comparison_message
	@bash assets/compare_results.sh && echo "All files pairwise comparisons were successful!" || echo "WARNING: Some pairwise comparisons yield different results!"


test_none: deepclean run_message copy_none validate_samplesheet report_results compare_results

test_single: deepclean run_message copy_single validate_samplesheet report_results compare_results

test_dual: deepclean run_message copy_dual validate_samplesheet report_results compare_results

test_mix: deepclean run_message copy_mix validate_samplesheet report_results compare_results

.PHONY: copy_none copy_single copy_dual copy_mix
.PHONY: validate_samplesheet compare_results
.PHONY: test_none test_single test_dual test_mix

.PHONY: cleanup_message
cleanup_message:
	@echo "-------------------------------------"
	@echo "- Past runs cleanup -----------------"
	@echo "-------------------------------------"

.PHONY: run_message
run_message:
	@echo "-------------------------------------"
	@echo "- Executing tests -------------------"
	@echo "-------------------------------------"

.PHONY: results_message
results_message:
	@echo "-------------------------------------"
	@echo "- Test results ($(SEQ_NUM)) -----------------"
	@echo "-------------------------------------"

.PHONY: comparison_message
comparison_message:
	@echo "-------------------------------------"
	@echo "- Original vs artificial comparison -"
	@echo "-------------------------------------"


.PHONY: deepclean
deepclean: cleanup_message clean_samplesheet clean_flowcell clean_data
	if [[ -d $(TEMP_DIR) ]]; then rm -r $(TEMP_DIR); fi

.PHONY: clean_samplesheet
clean_samplesheet:
	if [[ -f $(TEMP_DIR)/SampleSheet.csv ]]; then rm $(TEMP_DIR)/SampleSheet.csv; fi

.PHONY: clean_flowcell
clean_flowcell:
	if [[ -d $(TEMP_DIR)/$(FLOWCELL_ID) ]]; then rm -r $(TEMP_DIR)/$(FLOWCELL_ID); fi

.PHONY: clean_data
clean_data:
	if [[ -d $(TEMP_DIR)/Data ]]; then rm -r $(TEMP_DIR)/Data; fi

.PHONY: clean_demux
clean_demux:
	@if [[ -d $(TEMP_DIR)/Demultiplexing ]]; then rm -r $(TEMP_DIR)/Demultiplexing; fi
