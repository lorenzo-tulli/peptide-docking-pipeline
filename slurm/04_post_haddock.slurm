#!/bin/bash
#SBATCH --job-name=merge_recluster
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=logs/04_post_haddock_%j.out
#SBATCH --error=logs/04_post_haddock_%j.err

# =============================================================================
# Step 4: Collect all HADDOCK3 outputs, merge pool, re-cluster.
# Runs only after all 10 HADDOCK array tasks complete (enforced by submit_pipeline.sh).
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

echo "=== Step 4: Merge and re-cluster HADDOCK outputs ==="
echo "Started: $(date)"

# Point collect_and_recluster.py at the haddock3_runs directory
# (the script already looks in HADDOCK_BASE_DIR = haddock_runs,
#  but HADDOCK3 writes to haddock3_runs – patch it here via env var)
cd 04_post_haddock

# Symlink so the script finds the runs under the expected name
if [[ ! -e haddock_runs ]]; then
    ln -s ../03_haddock/haddock3_runs haddock_runs
fi

python collect_and_recluster.py

echo "Merged poses:  $(ls merged_pool/pose_*.pdb 2>/dev/null | wc -l)"
echo "Cluster reps:  $(ls clustered_poses/cluster_rep_*.pdb 2>/dev/null | wc -l)"
echo "Finished: $(date)"
