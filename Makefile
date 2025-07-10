# Makefile for the project
include Makefile.help

# Makefile containing the project's variables
PROJECT_NAME := BioMate
PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
UV=$(shell which uv)
SHELL=/bin/bash

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
