SHELL := /bin/bash
.PHONY: venv clean-venv plot plot3d

# Create a Python virtual environment and install requirements.
# Usage:
#   make venv
#   make venv VENV_DIR=.venv PYTHON=python3

VENV_DIR ?= .venv
PYTHON ?= python3

# Plotting targets
PLOT_SCRIPT := tools/plot_estimates.py
ARTIFACTS_DIR := artifacts

venv:
	@echo "Setting up virtualenv at $(VENV_DIR) using $(PYTHON)"
	@tools/setup_venv.sh || \
	{ echo "tools/setup_venv.sh failed; trying with VENV_DIR=$(VENV_DIR) PYTHON=$(PYTHON)"; \
	  VENV_DIR=$(VENV_DIR) PYTHON=$(PYTHON) bash tools/setup_venv.sh; }

clean-venv:
	@if [ -d "$(VENV_DIR)" ]; then \
		rm -rf "$(VENV_DIR)"; echo "Removed $(VENV_DIR)"; else \
		echo "No venv at $(VENV_DIR)"; fi

# Plot 2D results (after running tests)
plot:
	@echo "Generating 2D plot..."
	$(PYTHON) $(PLOT_SCRIPT) \
		--sensors $(ARTIFACTS_DIR)/sensors.csv \
		--targets $(ARTIFACTS_DIR)/targets.csv \
		--estimates $(ARTIFACTS_DIR)/estimates.csv \
		--cov $(ARTIFACTS_DIR)/covariances.npy \
		--outfile $(ARTIFACTS_DIR)/plot.png
	@echo "Plot saved to $(ARTIFACTS_DIR)/plot.png"

# Plot 3D results (requires 3D data)
# Usage: make plot3d
# Or with custom view angle: make plot3d ELEV=20 AZIM=45
ELEV ?= 30
AZIM ?= -60
plot3d:
	@echo "Generating 3D plot..."
	$(PYTHON) $(PLOT_SCRIPT) \
		--sensors $(ARTIFACTS_DIR)/sensors_3d.csv \
		--targets $(ARTIFACTS_DIR)/targets_3d.csv \
		--estimates $(ARTIFACTS_DIR)/estimates_3d.csv \
		--cov $(ARTIFACTS_DIR)/covariances_3d.npy \
		--3d \
		--elev $(ELEV) \
		--azim $(AZIM) \
		--outfile $(ARTIFACTS_DIR)/plot_3d.png
	@echo "3D plot saved to $(ARTIFACTS_DIR)/plot_3d.png"
