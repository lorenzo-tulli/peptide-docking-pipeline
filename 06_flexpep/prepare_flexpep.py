#!/usr/bin/env python3
"""
Step 6a Prepare PDB files for FlexPepDock refinement via PyRosetta.

FlexPepDock conventions:
  Chain A = receptor (protein)
  Chain B = peptide
  Peptide must already be in approximately the correct binding site
    (provided by HADDOCK output)
  Receptor and peptide must be in the same PDB; no ligands or waters

This script:
  1. Reads each top_pose_XX.pdb from 05_scoring/top5/
  2. Remaps chains: longest chain(s) A, shortest B
  3. Removes cofactors (MG, ATP) PyRosetta needs a clean structure
     (add them back for MD after FlexPepDock)
  4. Writes flexpep_inputs/fpd_pose_XX.pdb

The actual refinement is done by run_flexpep_pyrosetta.py (called by SLURM).

Dependencies:
  pip install MDAnalysis numpy

Usage:
  cd pipeline/06_flexpep
  python prepare_flexpep.py
"""

import os, sys, glob, shutil
import numpy as np
import MDAnalysis as mda

TOP5_DIR  = "../05_scoring/top5"
OUT_DIR   = "flexpep_inputs"
SKIP_RES  = {"MG", "ATP", "HOH", "WAT", "SOL"}   # residues to strip

os.makedirs(OUT_DIR, exist_ok=True)

def write_pdb_manual(atoms_rec, atoms_pep, out_path):
    """
    Write a clean PDB with receptor (chain A) and peptide (chain B).
    Uses MDAnalysis atom data; workaround for write-chain limitation.
    """
    lines = []

    def atom_line(atom, chain, serial):
        name = f"{atom.name:<4}"
        res  = atom.resname
        rid  = atom.resid
        x, y, z = atom.position
        elem = getattr(atom, "element", atom.name[0])
        return (f"ATOM  {serial:5d} {name} {res:3s} {chain}{rid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2s}\n")

    serial = 1
    for atom in atoms_rec:
        lines.append(atom_line(atom, "A", serial))
        serial += 1
    lines.append("TER\n")
    for atom in atoms_pep:
        lines.append(atom_line(atom, "B", serial))
        serial += 1
    lines.append("TER\nEND\n")

    with open(out_path, "w") as fh:
        fh.writelines(lines)

top5 = sorted(glob.glob(os.path.join(TOP5_DIR, "top_pose_*.pdb")))
if not top5:
    sys.exit(f"No top_pose_*.pdb in {TOP5_DIR}. Run 05_scoring first.")

fpd_pdbs = []

for i, pdb in enumerate(top5, 1):
    u = mda.Universe(pdb)

    # Strip cofactors / solvent
    clean = u.select_atoms(
        "not (resname " + " ".join(SKIP_RES) + ")"
    )

    # Identify peptide vs receptor:
    # Group by segment; assume the segment with fewest residues is the peptide
    segments = list(clean.segments)
    seg_sizes = [(len(s.residues), s) for s in segments]
    seg_sizes.sort(key=lambda x: x[0])

    if len(seg_sizes) >= 2:
        pep_seg = seg_sizes[0][1]
        rec_segs = [s for _, s in seg_sizes[1:]]
        pep_atoms = pep_seg.atoms
        rec_indices = np.concatenate([s.atoms.indices for s in rec_segs])
        rec_atoms = u.atoms[rec_indices]
    else:
        # Fallback: residues 1-7 = peptide
        pep_atoms = clean.select_atoms("resid 1:7")
        rec_atoms = clean.select_atoms("not resid 1:7")

    out_pdb = os.path.join(OUT_DIR, f"fpd_pose_{i:02d}.pdb")
    write_pdb_manual(rec_atoms, pep_atoms, out_pdb)
    fpd_pdbs.append(out_pdb)
    print(f"Pose {i}: {out_pdb}")

print(f"\n{len(fpd_pdbs)} PDBs written to {OUT_DIR}/")
print("Next: sbatch slurm/06_flexpep_array.slurm")
print("      (calls run_flexpep_pyrosetta.py for each pose)")
