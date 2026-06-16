#!/bin/bash
#SBATCH --job-name=score_bsa
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=logs/05_score_%j.out
#SBATCH --error=logs/05_score_%j.err

# =============================================================================
# Step 5: Score cluster representatives (HADDOCK score + BSA),
#         write top 5 candidates for FlexPepDock.
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

echo "=== Step 5: Score and select top 5 ==="
echo "Started: $(date)"

cd 05_scoring
python rank_poses.py

echo "Top 5: $(ls top5/top_pose_*.pdb 2>/dev/null | wc -l) poses"
echo ""
echo "IMPORTANT: Visually inspect the top5/ structures in PyMOL before"
echo "           running FlexPepDock. Check Tyr-4 is near the ATP γ-phosphate."
echo ""
echo "Finished: $(date)"
