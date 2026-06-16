#!/usr/bin/env python3
"""
Step 5 Rank HADDOCK cluster representatives by HADDOCK score + BSA,
then write the top 5 candidates for FlexPepDock.

BSA (buried surface area) is computed with freeSASA:
  BSA = (SASA_receptor + SASA_peptide SASA_complex) / 2
  HADDOCK score alone measures the internal energy of the complex (van der Waals, electrostatics, AIR violations) 
  but doesn't directly tell you how well the two molecules bury their interface. A pose could have a good HADDOCK 
  score while barely touching the receptor. BSA compensates for this by measuring how much surface area becomes 
  hidden upon binding a large BSA means the peptide is genuinely embedded in the binding site, not just loosely associated.

Dependencies:
  pip install freesasa MDAnalysis pandas numpy

Input:  ../04_post_haddock/clustered_poses/cluster_rep_XX.pdb
        ../04_post_haddock/cluster_ranking.csv
Output: top5/top_pose_XX.pdb
        pose_rankings.csv

Usage:
  cd pipeline/05_scoring
  python rank_poses.py
"""

import os, sys, glob, shutil, re, tempfile
import numpy as np
import pandas as pd
import MDAnalysis as mda

try:
    import freesasa
    HAS_FREESASA = True
except ImportError:
    print("WARNING: freesasa not installed – BSA not computed.\n"
          "Install: pip install freesasa")
    HAS_FREESASA = False

CLUSTER_DIR  = "../04_post_haddock/clustered_poses"
SCORE_CSV    = "../04_post_haddock/cluster_ranking.csv"
TOP_N        = 10
WEIGHT_SCORE = 0.7    # HADDOCK score contribution to composite rank
WEIGHT_BSA   = 0.3    # BSA contribution (higher BSA → better)
PEPTIDE_SEL  = "resid 1:7 or chainid P or segid P"
RECEPTOR_SEL = "not (resid 1:7 or chainid P or segid P or resname MG ATP)"
OUT_DIR      = "top5"

os.makedirs(OUT_DIR, exist_ok=True)

pdbs = sorted(glob.glob(os.path.join(CLUSTER_DIR, "cluster_rep_*.pdb")))
if not pdbs:
    sys.exit(f"No cluster_rep_*.pdb in {CLUSTER_DIR}")

try:
    df_rank = pd.read_csv(SCORE_CSV)
except FileNotFoundError:
    df_rank = pd.DataFrame()

def calc_bsa(pdb_path):
    if not HAS_FREESASA:
        return None
    try:
        u   = mda.Universe(pdb_path)
        pep = u.select_atoms(PEPTIDE_SEL)
        rec = u.select_atoms(RECEPTOR_SEL)

        if len(pep) == 0 or len(rec) == 0:
            return None

        with tempfile.TemporaryDirectory() as tmp:
            pep_pdb = os.path.join(tmp, "pep.pdb")
            rec_pdb = os.path.join(tmp, "rec.pdb")
            pep.write(pep_pdb)
            rec.write(rec_pdb)

            sasa_complex = freesasa.calc(freesasa.Structure(pdb_path)).totalArea()
            sasa_pep     = freesasa.calc(freesasa.Structure(pep_pdb)).totalArea()
            sasa_rec     = freesasa.calc(freesasa.Structure(rec_pdb)).totalArea()

        return (sasa_rec + sasa_pep - sasa_complex) / 2.0
    except Exception as e:
        print(f"  BSA error ({os.path.basename(pdb_path)}): {e}")
        return None

print(f"Scoring {len(pdbs)} cluster representatives …\n")
rows = []
for pdb in pdbs:
    name = os.path.basename(pdb)
    bsa  = calc_bsa(pdb)

    # Try to get HADDOCK score from CSV
    hs = None
    if not df_rank.empty:
        idx = df_rank["medoid"].str.endswith(name) if "medoid" in df_rank.columns else pd.Series(dtype=bool)
        if idx.any():
            hs = df_rank.loc[idx.idxmax(), "haddock_score"]

    bsa_str = f"{bsa:.1f}" if bsa else "N/A"
    hs_str  = f"{hs:.2f}" if hs is not None else "N/A"
    print(f"  {name}  HADDOCK={hs_str}  BSA={bsa_str} Å²")
    rows.append({"pdb": pdb, "name": name, "haddock_score": hs, "bsa_A2": bsa})

df = pd.DataFrame(rows)

def normalise(series, lower_better=True):
    s = pd.to_numeric(series, errors="coerce")
    mn, mx = s.min(), s.max()
    if mn == mx:
        return pd.Series(0.5, index=series.index)
    norm = (s - mn) / (mx - mn)
    return norm if lower_better else (1.0 - norm)

df["score_norm"] = normalise(df["haddock_score"], lower_better=True)
df["bsa_norm"]   = normalise(df["bsa_A2"],         lower_better=False)  # higher BSA better
df["composite"]  = WEIGHT_SCORE * df["score_norm"] + WEIGHT_BSA * df["bsa_norm"]

df_sorted = df.sort_values("composite").reset_index(drop=True)
df_sorted.to_csv("pose_rankings.csv", index=False)
print(f"\nRankings saved → pose_rankings.csv")

print(f"\nTop {TOP_N} poses for FlexPepDock:")
for i, row in df_sorted.head(TOP_N).iterrows():
    dst    = os.path.join(OUT_DIR, f"top_pose_{i+1:02d}.pdb")
    shutil.copy(row["pdb"], dst)
    bsa_s  = f"{row['bsa_A2']:.1f} Å²" if pd.notna(row["bsa_A2"]) else "N/A"
    hs_s   = f"{row['haddock_score']:.2f}" if pd.notna(row["haddock_score"]) else "N/A"
    print(f"  {i+1}. {row['name']}  score={hs_s}  BSA={bsa_s}  → top_pose_{i+1:02d}.pdb")

print(f"\nVisual check recommended in PyMOL before FlexPepDock.")
print(f"Next: run 06_flexpep/prepare_flexpep.py")
