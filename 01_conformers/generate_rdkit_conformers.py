#!/usr/bin/env python3
"""
Step 1b Generate N_CONF conformers of the free peptide with RDKit ETKDG + MMFF94s.

MMFF94s = Merck Molecular Force Field 94, static variant. It's used here instead of 
the plain MMFF94 because the "s" variant is more accurate for reproducing experimental 
geometries of organic molecules, and peptides in this context are treated as small 
molecules by RDKit

Dependencies:
  pip install rdkit numpy

Input:  peptide_reference.pdb   (from extract_peptide.py, chain P removed use peptide)
Output: conformers/conf_XXXX.pdb   (heavy atoms only)
        conformers/all_conformers.sdf

Usage:
  cd pipeline/01_conformers
  python generate_rdkit_conformers.py

Notes:
  RDKit ETKDG treats the peptide as a small molecule for generating diverse
  starting geometries before docking.
  If generation yields < 200 conformers, consider increasing NUM_ATTEMPTS or using
  the PEP-FOLD3 server instead (see pepfold3_instructions.md).
"""

import os, sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

IN_PDB      = "peptide_reference.pdb"
OUT_DIR     = "conformers"
N_CONF      = 5000
SEED        = int(np.random.randint(0, 2**31 - 1))
MMFF_ITER   = 2000
NUM_ATTEMPTS = 10000   # embedding attempts; raise if < N_CONF succeed

print(f"Random seed: {SEED}")
os.makedirs(OUT_DIR, exist_ok=True)

# Load
mol = Chem.MolFromPDBFile(IN_PDB, removeHs=False, sanitize=True)
if mol is None:
    sys.exit(
        "ERROR: RDKit could not parse peptide_reference.pdb.\n"
        "Try sanitize=False or check for non-standard residues."
    )

mol = Chem.AddHs(mol)
print(f"Loaded peptide: {mol.GetNumAtoms()} atoms (with H)")

# ETKDG v3
params = AllChem.ETKDGv3()
params.randomSeed           = SEED
params.numThreads           = 0        # use all available cores
params.useSmallRingTorsions = True
params.pruneRmsThresh       = 0.1      # prune conformers closer than 0.1 Å RMSD

print(f"Embedding {N_CONF} conformers (max attempts={NUM_ATTEMPTS})")
cids = AllChem.EmbedMultipleConfs(mol, N_CONF, params)
print(f"{len(cids)} conformers embedded")

if len(cids) == 0:
    sys.exit(
        "ERROR: No conformers generated.\n"
        "Try: mol = Chem.MolFromPDBFile(IN_PDB, removeHs=True, sanitize=False)"
    )

# MMFF94s minimisation
print("Minimising with MMFF94s")
results = AllChem.MMFFOptimizeMoleculeConfs(
    mol, mmffVariant="MMFF94s", maxIters=MMFF_ITER, numThreads=0
)
converged = sum(1 for r in results if r[0] == 0)
print(f"{converged}/{len(cids)} conformers converged")

# Write SDF (all conformers, no H)
mol_noh  = Chem.RemoveHs(mol)
sdf_path = os.path.join(OUT_DIR, "all_conformers.sdf")
writer   = Chem.SDWriter(sdf_path)
for cid in cids:
    writer.write(mol_noh, confId=cid)
writer.close()
print(f"SDF{sdf_path}")

# Write individual PDBs
for i, cid in enumerate(cids):
    Chem.MolToPDBFile(mol_noh, os.path.join(OUT_DIR, f"conf_{i:04d}.pdb"), confId=cid)

print(f"Individual PDBs {OUT_DIR}/conf_XXXX.pdb")
print(f"\nNext: run 02_clustering/rmsd_cluster.py")
