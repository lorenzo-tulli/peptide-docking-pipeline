#!/bin/bash
#SBATCH --job-name=haddock3
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16          # matches NCORES in generate_haddock3_configs.py
#SBATCH --mem=64G
#SBATCH --time=24:00:00             # 12h per representative (conservative)
#SBATCH --array=1-15                # one task per representative
#SBATCH --account=chem021482
#SBATCH --output=logs/03_haddock_%A_%a.out
#SBATCH --error=logs/03_haddock_%A_%a.err

# =============================================================================
# Step 3: HADDOCK3 docking – runs all 10 representatives as a SLURM array.
# SLURM_ARRAY_TASK_ID (1..10) selects the TOML config from toml_list.txt.
#
# Prerequisites:
#   cd pipeline/03_haddock
#   python prepare_haddock_inputs.py
#   python generate_haddock3_configs.py
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

TOML_LIST="03_haddock/haddock3_runs/toml_list.txt"

if [[ ! -f "$TOML_LIST" ]]; then
    echo "ERROR: $TOML_LIST not found."
    echo "Run: cd 03_haddock && python prepare_haddock_inputs.py && python generate_haddock3_configs.py"
    exit 1
fi

# Pick the TOML for this array task (1-indexed)
TOML=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TOML_LIST")

if [[ -z "$TOML" ]]; then
    echo "ERROR: No TOML found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "=== HADDOCK3 – representative $SLURM_ARRAY_TASK_ID ==="
echo "Config: $TOML"
echo "Started: $(date)"
echo "Node:    $(hostname)"
echo "Cores:   $SLURM_CPUS_PER_TASK"

haddock3 "$TOML"

echo ""
echo "Finished: $(date)"
