SHELL := /bin/bash
.PHONY: venv clean-venv

# Create a Python virtual environment and install requirements.
# Usage:
#   make venv
#   make venv VENV_DIR=.venv PYTHON=python3

VENV_DIR ?= .venv
PYTHON ?= python3

venv:
	@echo "Setting up virtualenv at $(VENV_DIR) using $(PYTHON)"
	@tools/setup_venv.sh || \
	{ echo "tools/setup_venv.sh failed; trying with VENV_DIR=$(VENV_DIR) PYTHON=$(PYTHON)"; \
	  VENV_DIR=$(VENV_DIR) PYTHON=$(PYTHON) bash tools/setup_venv.sh; }

clean-venv:
	@if [ -d "$(VENV_DIR)" ]; then \
		rm -rf "$(VENV_DIR)"; echo "Removed $(VENV_DIR)"; else \
		echo "No venv at $(VENV_DIR)"; fi
