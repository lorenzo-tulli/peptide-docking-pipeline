#!/usr/bin/env bash
# =============================================================================
# submit_pipeline.sh
# Submits the full docking pipeline as a chain of SLURM jobs with
# afterok dependencies so each step only starts after the previous succeeds.
#
# Usage:
#   cd pipeline
#   bash slurm/submit_pipeline.sh
#
# To re-run from a specific step, comment out the earlier sbatch calls
# and replace the first dependency with the last completed job ID.
# =============================================================================

set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PIPELINE_DIR"
mkdir -p logs

if [[ ! -f "01_conformers/peptide_reference.pdb" ]] && \
   [[ ! -f "03_haddock/haddock_inputs/ambig.tbl" ]]; then
    echo ""
    echo "NOTE: Running setup scripts before submitting the pipeline …"
    echo ""
fi

jobid() { echo "$1" | grep -oP '(?<=job )\d+'; }

echo "======================================================"
echo " Submitting pipeline from: $PIPELINE_DIR"
echo "======================================================"

echo ""
echo "[1/7] Conformers …"
J1=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    slurm/01_conformers.slurm)
JID1=$(jobid "$J1")
echo "  Submitted job $JID1"

echo "[2/7] Clustering …"
J2=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID1" \
    slurm/02_cluster.slurm)
JID2=$(jobid "$J2")
echo "  Submitted job $JID2  (depends on $JID1)"

echo "[2.5] Preparing HADDOCK inputs (login node) …"
echo "  Will execute when job $JID2 finishes."
echo "  Registering a wrapper job …"

PREP_SCRIPT=$(mktemp /tmp/prep_haddock_XXXX.sh)
cat > "$PREP_SCRIPT" << EOFPREP
#!/bin/bash
#SBATCH --job-name=prep_haddock
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=$PIPELINE_DIR/logs/prep_haddock_%j.out
#SBATCH --error=$PIPELINE_DIR/logs/prep_haddock_%j.err
set -euo pipefail
cd $PIPELINE_DIR/03_haddock
python prepare_haddock_inputs.py
python generate_haddock3_configs.py
EOFPREP

J_PREP=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID2" \
    "$PREP_SCRIPT")
JID_PREP=$(jobid "$J_PREP")
echo "  Submitted job $JID_PREP  (depends on $JID2)"

echo "[3/7] HADDOCK3 docking array (10 reps) …"
J3=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID_PREP" \
    slurm/03_haddock_array.slurm)
JID3=$(jobid "$J3")
echo "  Submitted array job $JID3  (depends on $JID_PREP)"
echo "  (10 tasks: ${JID3}_1 … ${JID3}_10)"

echo "[4/7] Merge and re-cluster …"
J4=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID3" \
    slurm/04_post_haddock.slurm)
JID4=$(jobid "$J4")
echo "  Submitted job $JID4  (depends on all HADDOCK tasks)"

echo "[5/7] Score + BSA ranking …"
J5=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID4" \
    slurm/05_score.slurm)
JID5=$(jobid "$J5")
echo "  Submitted job $JID5  (depends on $JID4)"

PREP_FPD=$(mktemp /tmp/prep_flexpep_XXXX.sh)
cat > "$PREP_FPD" << EOFPFPD
#!/bin/bash
#SBATCH --job-name=prep_flexpep
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --account=chem021482
#SBATCH --output=$PIPELINE_DIR/logs/prep_flexpep_%j.out
#SBATCH --error=$PIPELINE_DIR/logs/prep_flexpep_%j.err
set -euo pipefail
cd $PIPELINE_DIR/06_flexpep
python prepare_flexpep.py
EOFPFPD

J_PFD=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID5" \
    "$PREP_FPD")
JID_PFD=$(jobid "$J_PFD")
echo "  Submitted FlexPepDock prep job $JID_PFD  (depends on $JID5)"

echo "[6/7] FlexPepDock array (5 poses) …"
J6=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID_PFD" \
    slurm/06_flexpep_array.slurm)
JID6=$(jobid "$J6")
echo "  Submitted array job $JID6  (depends on $JID_PFD)"
echo "  (5 tasks: ${JID6}_1 … ${JID6}_5)"

echo "[7/7] Final ranking …"
J7=$(sbatch \
    --chdir="$PIPELINE_DIR" \
    --dependency=afterok:"$JID6" \
    slurm/07_final.slurm)
JID7=$(jobid "$J7")
echo "  Submitted job $JID7  (depends on all FlexPepDock tasks)"

echo ""
echo "======================================================"
echo " All jobs submitted. Dependency chain:"
echo ""
echo "  $JID1 (conformers)"
echo "    → $JID2 (cluster)"
echo "      → $JID_PREP (prepare haddock)"
echo "        → ${JID3}[1-10] (HADDOCK3 array)"
echo "          → $JID4 (merge+recluster)"
echo "            → $JID5 (score+BSA)"
echo "              → $JID_PFD (prepare flexpep)"
echo "                → ${JID6}[1-5] (FlexPepDock array)"
echo "                  → $JID7 (final ranking)"
echo ""
echo " Monitor: squeue -u \$USER"
echo " Logs:    pipeline/logs/"
echo ""
echo " IMPORTANT: After step 5 ($JID5) completes, visually inspect"
echo "  pipeline/05_scoring/top5/*.pdb in PyMOL before FlexPepDock starts."
echo "  To pause: scancel $JID_PFD $JID6 $JID7"
echo "======================================================"
