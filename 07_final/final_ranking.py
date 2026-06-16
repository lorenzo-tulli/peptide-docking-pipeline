#!/usr/bin/env python3
"""
Step 7 Final ranking of PyRosetta FlexPepDock models by total_score.

Reads scores.csv from each of the FlexPepDock output directories,
ranks by total_score, picks the best model per pose, and copies the
PDB files to final_top5/ ready for AMBER MD.

Dependencies:
  pip install pandas numpy

Input:
  ../06_flexpep/flexpep_outputs/pose_XX/scores.csv
  ../06_flexpep/flexpep_outputs/pose_XX/model_XXXX.pdb

Output:
  final_top5/md_start_XX.pdb
  final_ranking.csv
  all_flexpep_scores.csv

Usage:
  cd pipeline/07_final
  python final_ranking.py
"""

import os, sys, glob, shutil
import pandas as pd

FLEXPEP_OUT = "../06_flexpep/flexpep_outputs"
OUT_DIR     = "final_top5"
SCORE_COL   = "total_score"

os.makedirs(OUT_DIR, exist_ok=True)

csv_files = sorted(glob.glob(os.path.join(FLEXPEP_OUT, "pose_*/scores.csv")))

if not csv_files:
    sys.exit(
        f"No scores.csv files found under {FLEXPEP_OUT}/\n"
        "Ensure all FlexPepDock array tasks have completed:\n"
        "  squeue -u $USER\n"
        "  ls 06_flexpep/flexpep_outputs/pose_*/scores.csv"
    )

print(f"Found {len(csv_files)} score files.\n")

dfs = []
for csv in csv_files:
    pose_dir  = os.path.basename(os.path.dirname(csv))  # pose_01, pose_02 …
    df        = pd.read_csv(csv)
    df["run"] = pose_dir
    dfs.append(df)

    n_models = len(df)
    best     = df[SCORE_COL].min()
    best_model = df.loc[df[SCORE_COL].idxmin(), "model"]
    print(f"  {pose_dir}: {n_models} models  "
          f"best {SCORE_COL} = {best:.2f}  (model {best_model:04d})")

combined = pd.concat(dfs, ignore_index=True)
combined  = combined.sort_values(SCORE_COL).reset_index(drop=True)
combined.to_csv("all_flexpep_scores.csv", index=False)
print(f"\nAll {len(combined)} models → all_flexpep_scores.csv")

best_idx   = combined.groupby("run")[SCORE_COL].idxmin()
best_models = combined.loc[best_idx].sort_values(SCORE_COL).reset_index(drop=True)

print(f"\nTop 5 starting structures for AMBER MD:")
print(f"{'Rank':<5} {'Pose':<10} {'Model':<8} {SCORE_COL:<14} {'interface_dG':<14}")
print("-" * 55)

for rank, row in best_models.head(10).iterrows():
    out_pdb   = os.path.join(OUT_DIR, f"md_start_{rank+1:02d}.pdb")
    src_pdb   = str(row["pdb"])

    i_dg = f"{row['interface_dG']:.2f}" if pd.notna(row.get("interface_dG")) else "N/A"
    print(f"  {rank+1:<3}  {row['run']:<10}  {int(row['model']):04d}    "
          f"{row[SCORE_COL]:<14.2f}  {i_dg}")

    if os.path.isfile(src_pdb):
        shutil.copy(src_pdb, out_pdb)
    else:
        # Path stored in CSV may be absolute from a different working dir;
        # reconstruct relative to FLEXPEP_OUT
        fallback = os.path.join(
            FLEXPEP_OUT,
            row["run"],
            f"model_{int(row['model']):04d}.pdb"
        )
        if os.path.isfile(fallback):
            shutil.copy(fallback, out_pdb)
        else:
            print(f"     WARNING: PDB not found at {src_pdb} or {fallback}")
            continue

print(f"\nMD starting structures → {OUT_DIR}/md_start_XX.pdb")

best_models.head(10).assign(
    md_start_pdb=[f"md_start_{i+1:02d}.pdb" for i in range(min(5, len(best_models)))]
).to_csv("final_ranking.csv", index=False)
print(f"Final ranking table   → final_ranking.csv")
print(
    "\nNext steps for AMBER MD:"
    "\n  1. Re-add Mg²⁺ and ATP from 01_conformers/receptor.pdb"
    "\n     (they were stripped for PyRosetta — add them back in PyMOL)"
    "\n  2. pdb4amber -i md_start_XX.pdb -o md_start_XX_amber.pdb"
    "\n  3. tleap → solvate → minimise → equilibrate → production MD"
)
