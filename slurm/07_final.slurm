#!/bin/bash
#SBATCH --job-name=final_ranking
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=logs/07_final_%j.out
#SBATCH --error=logs/07_final_%j.err

# =============================================================================
# Step 7: Parse Rosetta scores, rank by total_score, extract 5 MD-ready PDBs.
# =============================================================================

set -euo pipefail
PIPELINE_DIR="$(pwd)"
mkdir -p logs

echo "=== Step 7: Final ranking and PDB extraction ==="
echo "Started: $(date)"

cd 07_final
python final_ranking.py

echo ""
echo "=== PIPELINE COMPLETE ==="
echo "MD starting structures: $(ls final_top5/md_start_*.pdb 2>/dev/null | wc -l)"
echo ""
echo "Next steps:"
echo "  1. Add Mg²⁺ and ATP back (they were stripped for FlexPepDock)"
echo "  2. pdb4amber -i md_start_XX.pdb -o md_start_XX_amber.pdb"
echo "  3. tleap → solvate → AMBER MD"
echo ""
echo "Finished: $(date)"
