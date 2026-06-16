#!/usr/bin/env bash
# =============================================================================
# 01_install_haddock3.sh
# Installs HADDOCK3 from source into the pep_docking conda environment.
#
# HADDOCK3 is a modular Python-based docking engine.
# GitHub: https://github.com/haddocking/haddock3
# Docs:   https://www.bonvinlab.org/haddock3/
#
# Usage:
#   conda activate pep_docking
#   bash install/01_install_haddock3.sh
# =============================================================================

set -euo pipefail

INSTALL_DIR="$HOME/software/haddock3"

echo "=== Installing HADDOCK3 ==="

# ── Check conda env ───────────────────────────────────────────────────────────
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: activate the pep_docking conda environment first."
    echo "  conda activate pep_docking"
    exit 1
fi

# ── Clone HADDOCK3 ────────────────────────────────────────────────────────────
mkdir -p "$(dirname "$INSTALL_DIR")"
if [[ -d "$INSTALL_DIR" ]]; then
    echo "Updating existing HADDOCK3 clone …"
    git -C "$INSTALL_DIR" pull
else
    git clone --depth 1 https://github.com/haddocking/haddock3.git "$INSTALL_DIR"
fi

cd "$INSTALL_DIR"

# ── Install into active conda environment ────────────────────────────────────
pip install -e ".[all]"

# ── Verify installation ───────────────────────────────────────────────────────
echo ""
haddock3 --version && echo "HADDOCK3 installed successfully." || \
    echo "ERROR: haddock3 command not found – check pip install output above."

# ── CNS (optional) ───────────────────────────────────────────────────────────
# HADDOCK3 does NOT require CNS for the modules used in this pipeline
# (topoaa, rigidbody, flexref, emref, clustfcc).
# If you need CNS-based modules, contact bonvinlab@science.uu.nl for access.

echo ""
echo "=== HADDOCK3 installation complete ==="
echo "Test: haddock3 --help"
echo ""
echo "Next: bash install/02_install_rosetta.sh"
