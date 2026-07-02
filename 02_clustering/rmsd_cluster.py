#!/usr/bin/env python3
"""
Step 2 RMSD clustering of peptide conformers to obtain representative structures.

Algorithm: 

1) For each conformer from RDKit, extract the N, CA, C, O atoms (backbone).
2) Each conformer is superimposed onto the reference peptide backbone using 
the Kabsch algorithm (finds the optimal rotation/translation that minimises RMSD 
between two sets of points) and filtering only conformers where:
  -TYR4:OH is within 0.6nm ASP171:CG
  -GLU3:CD is within 0.4nm LYS213:NZ
3) For every pair of surviving conformers, compute the RMSD between their backbone coordinates.
4) The RMSD matrix is fed into scipy's hierarchical clustering using average linkage: the distance 
between two clusters is the average of all pairwise RMSDs between their members. 
The tree is cut to produce exactly 10 clusters, then the medoid is selected.

Dependencies:
  pip install MDAnalysis numpy scipy pandas

Input:  ../01_conformers/conformers/conf_XXXX.pdb
Output: representatives/rep_XX.pdb   (10 files, HADDOCK-ready)
        cluster_summary.csv

Usage:
  cd pipeline/02_clustering
  python rmsd_cluster.py
"""

import os, sys, glob, shutil
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import MDAnalysis as mda

CONF_DIR    = "../01_conformers/conformers"
N_CLUSTERS  = 15
BB_SEL      = "backbone"
OUT_DIR     = "representatives"

# Reference complex for geometry pre-filter (must contain receptor + peptide)
REFERENCE_COMPLEX        = "/group/chem/oliveira2t/vi24769/EGFR_STRING/system-building/equilibration_protocol/1Mg/peptide_docking_pipeline/WT-DIMER_peptide_substrate_ATP-Mg.pdb"
FILTER_TYR4_OH_ASP171_CG = 6.0
FILTER_GLU3_CD_LYS213_NZ = 4.0

os.makedirs(OUT_DIR, exist_ok=True)

print("Loading reference complex for geometry pre-filter")
u_ref = mda.Universe(REFERENCE_COMPLEX)

ref_pep_bb = u_ref.select_atoms("resid 1:7 and backbone")
if len(ref_pep_bb) == 0:
    sys.exit("ERROR: no backbone atoms found for resid 1:7 in REFERENCE_COMPLEX.\n"
             "Check the path and peptide residue numbering.")

asp171_cg = u_ref.select_atoms("resid 171 and name CG")
lys213_nz = u_ref.select_atoms("resid 213 and name NZ")
if len(asp171_cg) == 0:
    sys.exit("ERROR: ASP171:CG not found in reference complex.\n"
             "Verify residue number and chain in REFERENCE_COMPLEX.")
if len(lys213_nz) == 0:
    sys.exit("ERROR: LYS213:NZ not found in reference complex.\n"
             "Verify residue number and chain in REFERENCE_COMPLEX.")

REF_BB_COORDS   = ref_pep_bb.positions.copy()
ASP171_CG_POS   = asp171_cg.positions[0].copy()
LYS213_NZ_POS   = lys213_nz.positions[0].copy()
N_REF_BB_ATOMS  = len(REF_BB_COORDS)

print(f"  Peptide backbone atoms in reference : {N_REF_BB_ATOMS}")
print(f"  ASP171:CG  {ASP171_CG_POS}")
print(f"  LYS213:NZ  {LYS213_NZ_POS}")


def _kabsch(mob, ref):
    """Return (R, t) minimising ||mob @ R.T + t – ref||.
    mob and ref are (N,3) arrays of corresponding atom positions."""
    mob_c = mob - mob.mean(axis=0)
    ref_c = ref - ref.mean(axis=0)
    H     = mob_c.T @ ref_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = ref.mean(axis=0) - mob.mean(axis=0) @ R.T
    return R, t


def passes_geometry_filter(pdb_path):
    """Superimpose conformer backbone onto reference peptide (Kabsch), then
    check TYR4:OH-ASP171:CG and GLU3:CD-LYS213:NZ distances."""
    try:
        u      = mda.Universe(pdb_path)
        mob_bb = u.select_atoms(BB_SEL)

        if len(mob_bb) != N_REF_BB_ATOMS:
            return False   # atom-count mismatch skip

        R, t = _kabsch(mob_bb.positions, REF_BB_COORDS)

        tyr4_oh = u.select_atoms("resid 4 and name OH")
        glu3_cd = u.select_atoms("resid 3 and name CD")
        if len(tyr4_oh) == 0 or len(glu3_cd) == 0:
            return False

        tyr4_pos = tyr4_oh.positions[0] @ R.T + t
        glu3_pos = glu3_cd.positions[0] @ R.T + t

        d1 = float(np.linalg.norm(tyr4_pos - ASP171_CG_POS))
        d2 = float(np.linalg.norm(glu3_pos - LYS213_NZ_POS))
        return d1 <= FILTER_TYR4_OH_ASP171_CG and d2 <= FILTER_GLU3_CD_LYS213_NZ
    except Exception:
        return False


pdbs = sorted(glob.glob(os.path.join(CONF_DIR, "conf_*.pdb")))
if not pdbs:
    sys.exit(f"No conf_*.pdb files found in {CONF_DIR}.\nRun 01_conformers first.")
print(f"\nFound {len(pdbs)} conformers.")

print("Loading backbone coordinates")
coords_list = []
valid_pdbs  = []

for pdb in pdbs:
    try:
        u  = mda.Universe(pdb)
        ag = u.select_atoms(BB_SEL)
        if len(ag) == 0:
            ag = u.select_atoms("not name H*")
        coords_list.append(ag.positions.copy())
        valid_pdbs.append(pdb)
    except Exception as e:
        print(f"  WARNING skipping {pdb}: {e}")

print(f"  → {len(valid_pdbs)} conformers loaded.")

print(f"\nApplying geometry pre-filter (superimpose onto reference peptide)")
print(f"  TYR4:OH  → ASP171:CG  ≤ {FILTER_TYR4_OH_ASP171_CG} Å")
print(f"  GLU3:CD  → LYS213:NZ  ≤ {FILTER_GLU3_CD_LYS213_NZ} Å")

keep       = [passes_geometry_filter(p) for p in valid_pdbs]
valid_pdbs  = [p for p, k in zip(valid_pdbs,  keep) if k]
coords_list = [c for c, k in zip(coords_list, keep) if k]

n_pass = len(valid_pdbs)
print(f"  {n_pass}/{len(pdbs)} conformers pass.")

if n_pass == 0:
    sys.exit("No conformers passed the geometry filter.\n"
             "Check REFERENCE_COMPLEX path and consider relaxing the thresholds.")
if n_pass < N_CLUSTERS:
    print(f"  WARNING: only {n_pass} conformers pass – reducing N_CLUSTERS to {n_pass}.")
    N_CLUSTERS = n_pass

n       = len(valid_pdbs)
n_atoms = min(c.shape[0] for c in coords_list)
coords  = np.array([c[:n_atoms] for c in coords_list])
coords -= coords.mean(axis=1, keepdims=True)   # centre each conformer

print("\nComputing pairwise RMSD matrix", flush=True)
rmsd_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        diff = coords[i] - coords[j]
        v    = np.sqrt((diff**2).sum(axis=1).mean())
        rmsd_matrix[i, j] = v
        rmsd_matrix[j, i] = v
    if i % 100 == 0:
        print(f"  {i}/{n} …", flush=True)

print(f"RMSD range: 0 – {rmsd_matrix.max():.2f} Å")

Z      = linkage(squareform(rmsd_matrix, checks=False), method="average")
labels = fcluster(Z, t=N_CLUSTERS, criterion="maxclust")

print(f"\nCluster composition:")
rows = []
for c in range(1, N_CLUSTERS + 1):
    members = np.where(labels == c)[0]
    sub       = rmsd_matrix[np.ix_(members, members)]
    medoid_i  = members[np.argmin(sub.sum(axis=1))]
    rep_path  = os.path.join(OUT_DIR, f"rep_{c:02d}.pdb")
    shutil.copy(valid_pdbs[medoid_i], rep_path)
    avg_rmsd  = sub[sub > 0].mean() if sub[sub > 0].size > 0 else 0.0
    print(f"  Cluster {c:2d} : {len(members):3d} members  "
          f"avg RMSD {avg_rmsd:.2f} Å  medoid → rep_{c:02d}.pdb")
    rows.append({
        "cluster": c, "size": len(members),
        "avg_intra_rmsd": round(avg_rmsd, 3),
        "medoid_pdb": os.path.basename(valid_pdbs[medoid_i]),
        "rep_file": f"rep_{c:02d}.pdb"
    })

pd.DataFrame(rows).to_csv("cluster_summary.csv", index=False)
print(f"\nSummary → cluster_summary.csv")
print(f"Representatives → {OUT_DIR}/rep_XX.pdb")
print(f"\nNext: run 03_haddock/prepare_haddock_inputs.py")
