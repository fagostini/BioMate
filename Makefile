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

OUTPUT_DIR = temp
FLOWCELL_ID = 20260310_LM43899_0385_A12GGASZR5

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

$(OUTPUT_DIR):
	@mkdir -p $(OUTPUT_DIR)

copy_none: $(OUTPUT_DIR) $(NISS) clean_samplesheet
	@cp $(NISS) $(OUTPUT_DIR)/SampleSheet.csv && echo "Copied SampleSheet (No Index) to working folder." || "Error copying SampleSheet to working folder!"

copy_single: $(OUTPUT_DIR) $(SISS) clean_samplesheet
	@cp $(SISS) $(OUTPUT_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Single Index) to working folder." || "Error copying SampleSheet to working folder!"

copy_dual: $(OUTPUT_DIR) $(DISS) clean_samplesheet
	@cp $(DISS) $(OUTPUT_DIR)/SampleSheet.csv && echo "Copied SampleSheet (Dual Indexes) to working folder." || "Error copying SampleSheet to working folder!"

$(OUTPUT_DIR)/$(FLOWCELL_ID): $(UV) $(OUTPUT_DIR) clean_flowcell
	@$(UV) run biomate --verbose blabber --format fastq --sample-sheet $(OUTPUT_DIR)/SampleSheet.csv --seq-number 10 --output $(OUTPUT_DIR) --flowcell-id $(FLOWCELL_ID) > $(OUTPUT_DIR)/blabber.out 2> $(OUTPUT_DIR)/blabber.err && echo "Blabber module executed successfully!" || "Error executing Blabber module!"


$(OUTPUT_DIR)/Data $(OUTPUT_DIR)/RunInfo.xml: $(UV) $(OUTPUT_DIR)/$(FLOWCELL_ID) clean_data
	@$(UV) run biomate --verbose fastrewind --input-path $(OUTPUT_DIR) --output-path $(OUTPUT_DIR) > $(OUTPUT_DIR)/fastrewind.out 2> $(OUTPUT_DIR)/fastrewind.err  && echo "Fastrewind module executed successfully!" || "Error executing Fastrewind module!"

validate_samplesheet: $(BCLCONVERT) $(OUTPUT_DIR)/RunInfo.xml
	@$(BCLCONVERT) --output-directory $(OUTPUT_DIR)/Demultiplexed --bcl-input-directory $(OUTPUT_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true --bcl-validate-sample-sheet-only true > $(OUTPUT_DIR)/bcl-validate.out 2> $(OUTPUT_DIR)/bcl-validate.err && echo "BCL-convert SampleSheet validation was successful!" || "Error executing BCL-convert SampleSheet validation!"

$(OUTPUT_DIR)/Demultiplexed: $(BCLCONVERT) $(OUTPUT_DIR)/RunInfo.xml clean_demux
	@$(BCLCONVERT) --output-directory $(OUTPUT_DIR)/Demultiplexed --bcl-input-directory $(OUTPUT_DIR) --strict-mode true --bcl-sampleproject-subdirectories true --sample-name-column-enabled true > $(OUTPUT_DIR)/bcl-convert.out 2> $(OUTPUT_DIR)/bcl-convert.err && echo "BCL-convert executed successfully!" || "Error executing BCL-convert!"


.PHONY: report_results
report_results: $(OUTPUT_DIR)/Demultiplexed results_message
	@find $(OUTPUT_DIR)/Demultiplexed -name "*.fastq.gz" | grep -v "Undetermined" | xargs zgrep -c ^@ || true
	@find $(OUTPUT_DIR)/Demultiplexed -name "*.fastq.gz" | grep "Undetermined" | xargs zgrep -c ^@ || true
	@echo "--------------------------------"


test_none: deepclean run_message copy_none validate_samplesheet report_results

test_single: deepclean run_message copy_single validate_samplesheet report_results

test_dual: deepclean run_message copy_dual validate_samplesheet report_results

.PHONY: cleanup_message
cleanup_message:
	@echo "--------------------------------"
	@echo "- Past runs cleanup ------------"
	@echo "--------------------------------"

.PHONY: run_message
run_message:
	@echo "--------------------------------"
	@echo "- Executing tests --------------"
	@echo "--------------------------------"

.PHONY: results_message
results_message:
	@echo "--------------------------------"
	@echo "- Test results -----------------"
	@echo "--------------------------------"


.PHONY: deepclean
deepclean: cleanup_message clean_samplesheet clean_flowcell clean_data
	if [[ -d $(OUTPUT_DIR) ]]; then rm -r $(OUTPUT_DIR); fi

.PHONY: clean_samplesheet
clean_samplesheet:
	if [[ -f $(OUTPUT_DIR)/SampleSheet.csv ]]; then rm $(OUTPUT_DIR)/SampleSheet.csv; fi

.PHONY: clean_flowcell
clean_flowcell:
	if [[ -d $(OUTPUT_DIR)/$(FLOWCELL_ID) ]]; then rm -r $(OUTPUT_DIR)/$(FLOWCELL_ID); fi

.PHONY: clean_data
clean_data:
	if [[ -d $(OUTPUT_DIR)/Data ]]; then rm -r $(OUTPUT_DIR)/Data; fi

.PHONY: clean_demux
clean_demux:
	@if [[ -d $(OUTPUT_DIR)/Demultiplexed ]]; then rm -r $(OUTPUT_DIR)/Demultiplexed; fi
