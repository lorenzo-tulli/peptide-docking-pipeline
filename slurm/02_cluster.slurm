#!/bin/bash
#SBATCH --job-name=pep_cluster
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=logs/02_cluster_%j.out
#SBATCH --error=logs/02_cluster_%j.err

# =============================================================================
# Step 2: RMSD clustering of 500 conformers → 10 representatives
# The pairwise RMSD matrix (500×500) is the bottleneck; expect ~15 min.
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

echo "=== Step 2: RMSD clustering ==="
echo "Started: $(date)"

cd 02_clustering
python rmsd_cluster.py

echo "Representatives: $(ls representatives/rep_*.pdb 2>/dev/null | wc -l)"
echo "Finished: $(date)"
