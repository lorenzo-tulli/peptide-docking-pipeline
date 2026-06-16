#!/bin/bash
#SBATCH --job-name=pep_conformers
#SBATCH --partition=compute               # adjust: short / cpu / compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=logs/01_conformers_%j.out
#SBATCH --error=logs/01_conformers_%j.err

# =============================================================================
# Step 1: Extract peptide/receptor and generate 500 RDKit conformers
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

echo "=== Step 1: Conformer generation ==="
echo "Started: $(date)"
echo "Node:    $(hostname)"
echo "Cores:   $SLURM_CPUS_PER_TASK"

# 1a – split PDB into peptide and receptor
cd 01_conformers
python extract_peptide.py

# 1b – generate 500 RDKit conformers (uses all available cores via RDKit)
python generate_rdkit_conformers.py

echo ""
echo "Finished: $(date)"
echo "Conformers: $(ls conformers/conf_*.pdb 2>/dev/null | wc -l)"
