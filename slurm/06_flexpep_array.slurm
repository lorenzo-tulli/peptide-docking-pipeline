#!/bin/bash
#SBATCH --job-name=flexpep_pyrosetta
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1           # PyRosetta is single-threaded
#SBATCH --mem=16G                   # PyRosetta needs ~8-12 GB per process
#SBATCH --time=48:00:00             # 200 models × ~8 min each ≈ 27 h; 48 h is safe
#SBATCH --account=chem021482
#SBATCH --array=1-10                 # one task per top HADDOCK pose
#SBATCH --output=logs/06_flexpep_%A_%a.out
#SBATCH --error=logs/06_flexpep_%A_%a.err

# =============================================================================
# Step 6: FlexPepDock refinement via PyRosetta
#
# Each array task refines one of the 5 top HADDOCK poses.
# 200 models are generated sequentially within each task.
# All 5 tasks run in parallel (subject to cluster availability).
#
# Prerequisites:
#   conda activate pep_docking
#   python test_pyrosetta.py          ← must pass 4/4 before submitting
#   cd pipeline/06_flexpep
#   python prepare_flexpep.py         ← creates flexpep_inputs/fpd_pose_XX.pdb
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

POSE_ID=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")

echo "=== FlexPepDock (PyRosetta) – pose $POSE_ID ==="
echo "Node:    $(hostname)"
echo "Started: $(date)"

python 06_flexpep/run_flexpep_pyrosetta.py \
    --pose_id   "$POSE_ID" \
    --n_models  200 \
    --input_dir  "06_flexpep/flexpep_inputs" \
    --output_dir "06_flexpep/flexpep_outputs"

echo ""
echo "Finished: $(date)"
SCORE_CSV="06_flexpep/flexpep_outputs/pose_${POSE_ID}/scores.csv"
if [[ -f "$SCORE_CSV" ]]; then
    N=$(tail -n +2 "$SCORE_CSV" | wc -l)
    echo "Models scored: $N"
fi
