#!/usr/bin/env python3
"""
Step 4 Collect all HADDOCK outputs, merge into one pool, re-cluster,
and parse HADDOCK scores.

Handles both HADDOCK2.4 (water/ directory) and HADDOCK3 (emref/ directory)
output layouts.

Dependencies:
  pip install MDAnalysis numpy scipy pandas

Input layout (adjust HADDOCK_BASE_DIR):
  haddock_runs/
    run_rep_01/
      water/    (HADDOCK2.4)   OR   emref/  (HADDOCK3)
        *.pdb
    run_rep_02/ ...

Output:
  merged_pool/pose_XXXX.pdb         all collected structures
  merged_pool/haddock_scores.csv    parsed scores
  clustered_poses/cluster_rep_XX.pdb   re-clustered representatives

Usage:
  cd pipeline/04_post_haddock
  python collect_and_recluster.py
"""

import os, sys, glob, shutil, re, gzip
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import MDAnalysis as mda

HADDOCK_BASE_DIR = "haddock_runs"    # directory containing run_rep_01 … run_rep_10
N_RECLUST        = 40                # intermediate re-clustering before top-poses
PEPTIDE_RESIDS   = "resid 1:7"      # selection string for peptide in output PDBs
BB_SEL           = "backbone"
OUT_MERGED       = "merged_pool"
OUT_CLUSTERED    = "clustered_poses"

# Geometry filter applied to HADDOCK output poses (direct distances, no superposition)
FILTER_TYR4_OH_ASP171_CG = 6.0   # Å: TYR4:OH  → ASP171:CG
FILTER_GLU3_CD_LYS213_NZ = 4.0   # Å: GLU3:CD  → LYS213:NZ

os.makedirs(OUT_MERGED,    exist_ok=True)
os.makedirs(OUT_CLUSTERED, exist_ok=True)

# HADDOCK3 compresses output after the run (clean_steps); the final docked
# structures live in 6_seletopclusts/ as .pdb.gz files.
search_patterns = [
    f"{HADDOCK_BASE_DIR}/*/water/*.pdb",                # HADDOCK 2.4
    f"{HADDOCK_BASE_DIR}/*/6_seletopclusts/*.pdb",      # HADDOCK 3 (uncompressed)
    f"{HADDOCK_BASE_DIR}/*/6_seletopclusts/*.pdb.gz",   # HADDOCK 3 (compressed)
    f"{HADDOCK_BASE_DIR}/*/emref/*.pdb",                # HADDOCK 3 emref fallback
    f"{HADDOCK_BASE_DIR}/*/emref/*.pdb.gz",             # HADDOCK 3 emref fallback (compressed)
]

found = set()
for pat in search_patterns:
    found.update(glob.glob(pat, recursive=True))

if not found:
    sys.exit(
        f"No PDB files found under {HADDOCK_BASE_DIR}/\n"
        "Check HADDOCK_BASE_DIR path and ensure HADDOCK runs completed.\n"
        "Expected: run_rep_XX/6_seletopclusts/*.pdb.gz"
    )
found = sorted(found)
print(f"Found {len(found)} PDB files across all HADDOCK runs.")

def parse_haddock_score(pdb_path):
    """Returns HADDOCK score (float) or None if not found."""
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("REMARK"):
                continue
            # HADDOCK2.4 format: REMARK HADDOCK-SCORE: -123.45
            m = re.search(r"HADDOCK.SCORE[:\s]+([-\d.]+)", line, re.I)
            if m:
                return float(m.group(1))
            # HADDOCK3 format: REMARK energies: ...  score= -123.45
            m = re.search(r"score=\s*([-\d.]+)", line)
            if m:
                return float(m.group(1))
    return None

rows = []
for idx, pdb in enumerate(found):
    dst = os.path.join(OUT_MERGED, f"pose_{idx:04d}.pdb")
    if pdb.endswith(".gz"):
        with gzip.open(pdb, "rb") as f_in, open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(pdb, dst)
    score = parse_haddock_score(dst)
    rows.append({"pose_id": f"pose_{idx:04d}", "source": pdb,
                 "haddock_score": score})

df_scores = pd.DataFrame(rows)
df_scores.to_csv(os.path.join(OUT_MERGED, "haddock_scores.csv"), index=False)
n_scored = df_scores["haddock_score"].notna().sum()
print(f"  {n_scored}/{len(rows)} poses have HADDOCK scores.")
print(f"Merged pool → {OUT_MERGED}/")

def passes_geometry_filter(pdb_path):
    """Direct distance check on the docked complex (no superposition needed)."""
    try:
        u = mda.Universe(pdb_path)
        tyr4_oh  = u.select_atoms("resid 4 and name OH")
        glu3_cd  = u.select_atoms("resid 3 and name CD")
        asp171_cg = u.select_atoms("resid 171 and name CG")
        lys213_nz = u.select_atoms("resid 213 and name NZ")
        if any(len(ag) == 0 for ag in [tyr4_oh, glu3_cd, asp171_cg, lys213_nz]):
            return False
        d1 = float(np.linalg.norm(tyr4_oh.positions[0]  - asp171_cg.positions[0]))
        d2 = float(np.linalg.norm(glu3_cd.positions[0]  - lys213_nz.positions[0]))
        return d1 <= FILTER_TYR4_OH_ASP171_CG and d2 <= FILTER_GLU3_CD_LYS213_NZ
    except Exception:
        return False

print(f"\nApplying geometry filter to {len(all_pdbs_merged := sorted(glob.glob(os.path.join(OUT_MERGED, 'pose_*.pdb'))))} merged poses …")
print(f"  TYR4:OH  → ASP171:CG ≤ {FILTER_TYR4_OH_ASP171_CG} Å")
print(f"  GLU3:CD  → LYS213:NZ ≤ {FILTER_GLU3_CD_LYS213_NZ} Å")
geo_pass = [p for p in all_pdbs_merged if passes_geometry_filter(p)]
print(f"  {len(geo_pass)}/{len(all_pdbs_merged)} poses pass geometry filter.")
if len(geo_pass) == 0:
    sys.exit("No poses passed the geometry filter. Consider relaxing the thresholds.")

print("Loading peptide coordinates for RMSD re-clustering")
coords_list = []
valid_pdbs  = []

for pdb in geo_pass:
    try:
        u  = mda.Universe(pdb)
        pep = u.select_atoms(PEPTIDE_RESIDS)
        if len(pep) == 0:
            # Try chain-based selection
            pep = u.select_atoms("chainid P or segid P")
        bb = pep.select_atoms(BB_SEL)
        if len(bb) == 0:
            bb = pep.select_atoms("not name H*")
        if len(bb) == 0:
            continue
        coords_list.append(bb.positions.copy())
        valid_pdbs.append(pdb)
    except Exception as e:
        print(f"  WARNING – skipping {pdb}: {e}")

n       = len(valid_pdbs)
n_atoms = min(c.shape[0] for c in coords_list)
coords  = np.array([c[:n_atoms] for c in coords_list])
coords -= coords.mean(axis=1, keepdims=True)   # centre

if n < N_RECLUST:
    N_RECLUST = max(n // 2, 1)

print(f"  {n} poses loaded. Computing RMSD matrix")
rmsd_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        diff = coords[i] - coords[j]
        v    = np.sqrt((diff**2).sum(axis=1).mean())
        rmsd_matrix[i, j] = v
        rmsd_matrix[j, i] = v

Z      = linkage(squareform(rmsd_matrix, checks=False), method="average")
labels = fcluster(Z, t=N_RECLUST, criterion="maxclust")

print(f"\nRe-clustering into {N_RECLUST} clusters:")

cluster_info = []
for c in range(1, N_RECLUST + 1):
    members = np.where(labels == c)[0]
    if len(members) == 0:
        continue
    sub       = rmsd_matrix[np.ix_(members, members)]
    medoid_i  = members[np.argmin(sub.sum(axis=1))]
    rep_pdb   = valid_pdbs[medoid_i]

    # Get HADDOCK score for medoid
    pose_id   = os.path.basename(rep_pdb)
    hs_row    = df_scores[df_scores["pose_id"] == pose_id.replace(".pdb", "")]
    hs        = hs_row["haddock_score"].values[0] if len(hs_row) else None

    cluster_info.append({
        "cluster": c, "size": len(members),
        "medoid": rep_pdb, "haddock_score": hs
    })
    print(f"  Cluster {c:2d}: {len(members):3d} members  score={hs}")

# Sort clusters by HADDOCK score (lower is better), then by size
df_cl = pd.DataFrame(cluster_info)
df_cl_sorted = df_cl.sort_values(
    ["haddock_score", "size"], ascending=[True, False]
).reset_index(drop=True)

print(f"\nTop clusters by HADDOCK score:")
for rank, row in df_cl_sorted.head(20).iterrows():
    dst = os.path.join(OUT_CLUSTERED, f"cluster_rep_{rank+1:02d}.pdb")
    shutil.copy(row["medoid"], dst)
    print(f"  {rank+1:2d}. cluster {int(row['cluster'])}  "
          f"size={int(row['size'])}  score={row['haddock_score']}  → {dst}")

df_cl_sorted.to_csv("cluster_ranking.csv", index=False)
print(f"\nFull cluster ranking → cluster_ranking.csv")
print(f"Candidate structures → {OUT_CLUSTERED}/cluster_rep_XX.pdb")
print(
    "\nNEXT: visually inspect cluster_rep_01..N_CLUSTER in PyMOL/VMD, "
    "then run 05_scoring/rank_poses.py for BSA-weighted ranking."
)
