# Makefile for the project
include Makefile.help

# Makefile containing the project's variables
PROJECT_NAME := BioMate
PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
UV=$(shell which uv)
BCLCONVERT=$(shell which bcl-convert)
SHELL=/bin/bash

# Test assets and variables
NISS = assets/SampleSheet_NoIndex.csv
SISS = assets/SampleSheet_SingleIndex.csv
DISS = assets/SampleSheet_DualIndex.csv
MISS = assets/SampleSheet_MixedIndexes.csv

TEMP_DIR = temp
FLOWCELL_ID = 20260310_LM43899_0385_A12GGASZR5
SEQ_NUM = 10

# Download and install uv
$(UV):
	curl -LsSf https://astral.sh/uv/install.sh | sh


.PHONY: sync
## Sync the project's dependencies with the environment
sync: $(UV)
	@$(UV) sync --all-extras --all-groups


uv.lock: $(UV)
	@$(UV) lock


.PHONY: lock
## Create a lockfile for the project's dependencies
lock: uv.lock


.PHONY: upgrade
## Upgrade the project's dependencies
upgrade: $(UV)
	@$(UV) lock --upgrade


.PHONY: build
## Build the project into distribution archives
build: $(UV)
	@$(UV) build


.venv/bin/activate: $(UV)
	@$(UV) venv


.PHONY: venv
## Create a new virtual environment
venv: .venv/bin/activate


requirements.txt: $(UV) uv.lock
	@$(UV) export --output-file requirements.txt


.PHONY: export
## Export the project's dependencies to a requirements.txt file
export: requirements.txt


.PHONY: interrogate
## Run interrogate in verbose mode to check for code quality
interrogate: $(UV)
	@$(UV)x interrogate --verbose

.PHONY: deploy
## Deploy documentation on GitHub
deploy:
	@$(UV) run mkdocs gh-deploy

$(TEMP_DIR):
	@mkdir -p $(TEMP_DIR)

copy_none: $(TEMP_DIR) $(NISS) clean_samplesheet
	@cp $(NISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (No Index) to working folder." || "Error copying SampleSheet to working folder!"

copy_single: $(TEMP_DIR) $(SISS) clean_samplesheet
	@cp $(SISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Single Index) to working folder." || "Error copying SampleSheet to working folder!"

copy_dual: $(TEMP_DIR) $(DISS) clean_samplesheet
	@cp $(DISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Dual Indexes) to working folder." || "Error copying SampleSheet to working folder!"

copy_mix: $(TEMP_DIR) $(MISS) clean_samplesheet
	@cp $(MISS) $(TEMP_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Mixed Indexes) to working folder." || "Error copying SampleSheet to working folder!"

$(TEMP_DIR)/$(FLOWCELL_ID): $(UV) $(TEMP_DIR) clean_flowcell
	@$(UV) run biomate --verbose blabber --format fastq --sample-sheet $(TEMP_DIR)/SampleSheet.csv --seq-number $(SEQ_NUM) --output $(TEMP_DIR) --flowcell-id $(FLOWCELL_ID) > $(TEMP_DIR)/blabber.out 2> $(TEMP_DIR)/blabber.err && echo "Blabber module executed successfully!" || "Error executing Blabber module!"

$(TEMP_DIR)/Data $(TEMP_DIR)/RunInfo.xml: $(UV) $(TEMP_DIR)/$(FLOWCELL_ID) clean_data
	@$(UV) run biomate --verbose fastrewind --input-path $(TEMP_DIR) --output-path $(TEMP_DIR) > $(TEMP_DIR)/fastrewind.out 2> $(TEMP_DIR)/fastrewind.err  && echo "Fastrewind module executed successfully!" || "Error executing Fastrewind module!"

validate_samplesheet: $(BCLCONVERT) $(TEMP_DIR)/RunInfo.xml
	@$(BCLCONVERT) --output-directory $(TEMP_DIR)/Demultiplexing --bcl-input-directory $(TEMP_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true > $(TEMP_DIR)/bcl-validate.out 2> $(TEMP_DIR)/bcl-validate.err && echo "BCL-convert SampleSheet validation was successful!" || "Error executing BCL-convert SampleSheet validation!"

$(TEMP_DIR)/Demultiplexing: $(BCLCONVERT) $(TEMP_DIR)/RunInfo.xml clean_demux
	@$(BCLCONVERT) --output-directory $(TEMP_DIR)/Demultiplexing --bcl-input-directory $(TEMP_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true > $(TEMP_DIR)/bcl-convert.out 2> $(TEMP_DIR)/bcl-convert.err && echo "BCL-convert executed successfully!" || "Error executing BCL-convert!"


.PHONY: report_results
report_results: $(TEMP_DIR)/Demultiplexing results_message
	@find $(TEMP_DIR)/Demultiplexing -name "*.fastq.gz" | grep -v "Undetermined" | xargs zgrep -c ^@ || true
	@find $(TEMP_DIR)/Demultiplexing -name "*.fastq.gz" | grep "Undetermined" | xargs zgrep -c ^@ || true


compare_results: $(TEMP_DIR)/Demultiplexing comparison_message
	@bash assets/compare_results.sh && echo "All files pairwise comparisons were successful!" || echo "WARNING: Some pairwise comparisons yield different results!"


test_none: deepclean run_message copy_none validate_samplesheet report_results compare_results

test_single: deepclean run_message copy_single validate_samplesheet report_results compare_results

test_dual: deepclean run_message copy_dual validate_samplesheet report_results compare_results

test_mix: deepclean run_message copy_mix validate_samplesheet report_results compare_results

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
