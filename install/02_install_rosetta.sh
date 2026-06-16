#!/usr/bin/env bash
# =============================================================================
# 02_install_rosetta.sh
# Compiles Rosetta (FlexPepDocking application) on the HPC.
#
# IMPORTANT: Rosetta requires a FREE academic license before downloading.
#   Register at: https://els2.comotion.uw.edu/product/rosetta
#   After registration you will receive a download link for the tarball.
#
# Also check first whether Rosetta is already installed on your HPC:
#   module avail rosetta
#   which FlexPepDocking.linuxgccrelease
#
# Usage (on a compute node or in an interactive session):
#   bash install/02_install_rosetta.sh
#
# Approximate compile time: 30–90 min on 16 cores.
# =============================================================================

set -euo pipefail

# ── Settings (adjust paths) ───────────────────────────────────────────────────
ROSETTA_TARBALL="$HOME/downloads/rosetta_src_*.tgz"  # path to downloaded tarball
ROSETTA_INSTALL="$HOME/software/rosetta"
NCORES=16   # compile threads — use a full node for speed

# ── Step 1: check if already installed ───────────────────────────────────────
echo "=== Checking for existing Rosetta installation ==="
if command -v FlexPepDocking.linuxgccrelease &>/dev/null; then
    echo "FlexPepDocking.linuxgccrelease already in PATH – skipping compilation."
    FlexPepDocking.linuxgccrelease --version 2>/dev/null || true
    exit 0
fi

# Try module system
for mod in rosetta Rosetta; do
    if module load "$mod" 2>/dev/null; then
        echo "Loaded module: $mod"
        if command -v FlexPepDocking.linuxgccrelease &>/dev/null; then
            echo "Rosetta available via module – no compilation needed."
            echo "Add 'module load $mod' to your ~/.bashrc"
            exit 0
        fi
        module unload "$mod"
    fi
done

# ── Step 2: Load required compilers ──────────────────────────────────────────
echo ""
echo "=== Loading compiler modules ==="
# Adjust module names for your HPC:
module load tools/cmake/3.20.0   2>/dev/null || true
module load languages/gcc/11.3.0 2>/dev/null || true
module load mpi/openmpi/4.1.1    2>/dev/null || true

# ── Step 3: Extract tarball ───────────────────────────────────────────────────
echo ""
echo "=== Extracting Rosetta tarball ==="
TARBALL=$(ls $ROSETTA_TARBALL 2>/dev/null | head -1)
if [[ -z "$TARBALL" ]]; then
    echo "ERROR: No Rosetta tarball found at $ROSETTA_TARBALL"
    echo ""
    echo "Download Rosetta from: https://els2.comotion.uw.edu/product/rosetta"
    echo "Place the .tgz file in $HOME/downloads/ then re-run this script."
    exit 1
fi

mkdir -p "$ROSETTA_INSTALL"
tar -xzf "$TARBALL" -C "$ROSETTA_INSTALL" --strip-components=1

# ── Step 4: Compile (MPI + static binary for HPC portability) ────────────────
echo ""
echo "=== Compiling Rosetta (this takes 30–90 min) ==="
cd "$ROSETTA_INSTALL/main/source"

# Compile FlexPepDocking only (much faster than full Rosetta build)
./scons.py \
    -j "$NCORES" \
    mode=release \
    extras=mpi \
    bin/FlexPepDocking

echo ""
echo "Compile finished."

# ── Step 5: Set up PATH ───────────────────────────────────────────────────────
BIN_DIR="$ROSETTA_INSTALL/main/source/bin"
echo ""
echo "=== Adding Rosetta to PATH ==="
echo "export ROSETTA_BIN=$BIN_DIR" >> "$HOME/.bashrc"
echo "export PATH=\$ROSETTA_BIN:\$PATH"   >> "$HOME/.bashrc"
echo "Added to ~/.bashrc:  export ROSETTA_BIN=$BIN_DIR"

# Verify
source "$HOME/.bashrc"
ls "$BIN_DIR"/FlexPepDocking* && echo "Rosetta FlexPepDocking binary found." || \
    echo "ERROR: binary not found – check compilation log above."

echo ""
echo "=== Rosetta installation complete ==="
echo "Test: FlexPepDocking.mpi.linuxgccrelease --help"
echo ""
echo "Next: run the pipeline with  bash slurm/submit_pipeline.sh"
