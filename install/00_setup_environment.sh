#!/usr/bin/env bash
# =============================================================================
# 00_setup_environment.sh
# Creates a conda environment with all Python dependencies for the pipeline.
#
# Tested on: Rocky Linux 8 / BlueCrystal BC4 (University of Bristol)
# Adjust CONDA_ROOT and ENV_NAME if needed.
#
# Usage (on HPC login node):
#   bash install/00_setup_environment.sh
# =============================================================================

set -euo pipefail

ENV_NAME="pep_docking"
PYTHON_VER="3.10"

echo "=== Creating conda environment: $ENV_NAME ==="

# ── Load conda (adjust module name for your HPC) ─────────────────────────────
# BC4 / BlueCrystal:
#   module load languages/anaconda3/2021-3.8.8
# ARCHER2 / others:
#   module load miniconda3
# If conda is already in PATH, comment out the module load line.
module load languages/anaconda3/2021-3.8.8 2>/dev/null || \
    module load miniconda3 2>/dev/null || \
    echo "WARNING: could not load conda module – assuming conda is already in PATH"

# ── Create environment ────────────────────────────────────────────────────────
conda create -y -n "$ENV_NAME" python="$PYTHON_VER"

# ── Activate ──────────────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "Python: $(python --version)"

# ── Core scientific stack (conda-forge for reliable HPC builds) ───────────────
conda install -y -c conda-forge \
    numpy \
    scipy \
    pandas \
    matplotlib \
    mdanalysis \
    rdkit \
    freesasa \
    biopython \
    tqdm

# ── pip-only packages ─────────────────────────────────────────────────────────
pip install --upgrade pip
pip install \
    haddock3 \
    pdb-tools

echo ""
echo "=== Environment '$ENV_NAME' ready ==="
echo "Activate with:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Next: bash install/01_install_haddock3.sh"
