#!/usr/bin/env python3
"""
Step 3a Prepare HADDOCK-compatible input files.

For each of the N_CLUSTER representative peptide conformers:
  Aigns chain P and writes to haddock_inputs/
For the receptor:
  Copies receptor.pdb (chains A, B, C) to haddock_inputs/

Also:
  Identifies interface residues from the reference complex
  Writes haddock_inputs/ambig.tbl   (HADDOCK AIR restraints, CNS format)
  Writes haddock3_run_TEMPLATE.toml (HADDOCK3 config template)

Dependencies:
  pip install MDAnalysis numpy

Usage:
  cd pipeline/03_haddock
  python prepare_haddock_inputs.py
"""

import os, sys, glob, shutil
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances as mda_dist

REP_DIR      = "../02_clustering/representatives"
PEPTIDE_REF  = "../01_conformers/peptide_reference.pdb"
RECEPTOR_PDB = "../01_conformers/receptor.pdb"
REFERENCE_COMPLEX = "/group/chem/oliveira2t/vi24769/EGFR_STRING/system-building/equilibration_protocol/1Mg/peptide_docking_pipeline/WT-DIMER_peptide_substrate_ATP-Mg.pdb"   # original complex for interface ID
OUT_DIR      = "haddock_inputs"
CONTACT_DIST = 5.0   # Å contact threshold for active residue detection

os.makedirs(OUT_DIR, exist_ok=True)

STRIP_RESNAMES = {"MG", "ATP", "ACE", "NME", "NHE", "NH2"}
rec_out = os.path.join(OUT_DIR, "receptor.pdb")
with open(RECEPTOR_PDB) as fi, open(rec_out, "w") as fo:
    for line in fi:
        if line.startswith(("ATOM", "HETATM")):
            chain   = line[21]
            resname = line[17:20].strip()
            if chain in ("A", "B") and resname not in STRIP_RESNAMES:
                fo.write(line[:21] + "A" + line[22:])  # relabel B → A
    fo.write("TER\nEND\n")
print(f"receptor.pdb → {rec_out}  (single chain A; MG/ATP/ACE/NME stripped)")

rep_pdbs = sorted(glob.glob(os.path.join(REP_DIR, "rep_*.pdb")))
if not rep_pdbs:
    sys.exit(f"No rep_*.pdb files in {REP_DIR}. Run 02_clustering first.")

for rp in rep_pdbs:
    name    = os.path.basename(rp).replace("rep_", "peptide_rep_")
    out_pdb = os.path.join(OUT_DIR, name)
    with open(rp) as fi, open(out_pdb, "w") as fo:
        for line in fi:
            if line.startswith(("ATOM", "HETATM")):
                fo.write(line[:21] + "P" + line[22:])
            else:
                fo.write(line)
    print(f"  {name} written")

print(f"\nIdentifying interface residues (d ≤ {CONTACT_DIST} Å) …")
u = mda.Universe(REFERENCE_COMPLEX)

pep_atoms = u.select_atoms("resid 1:7")
rec_atoms = u.select_atoms("resid 8:9999 and not (resname MG ATP)")

dist_mat  = mda_dist.distance_array(pep_atoms.positions, rec_atoms.positions)

pep_active = sorted({pep_atoms[i].resid for i in range(len(pep_atoms))
                     if dist_mat[i].min() <= CONTACT_DIST})
rec_active = sorted({rec_atoms[j].resid for j in range(len(rec_atoms))
                     if dist_mat[:, j].min() <= CONTACT_DIST})

# All receptor residues now live in segid A (chains A+B merged → single chain A)
print(f"  Peptide active residues (segid P) : {pep_active}")
print(f"  Receptor active residues (segid A): {sorted(rec_active)}")

def fmt_group(residues, segid):
    parts = [f"(resi {r} and segid {segid})" for r in residues]
    return "\n        or\n        ".join(parts)

tbl_out = os.path.join(OUT_DIR, "ambig.tbl")
with open(tbl_out, "w") as f:
    f.write("! Ambiguous Interaction Restraints for HADDOCK\n")
    f.write("! Peptide (segid P) ↔ Receptor (segid A, chains A+B merged)\n")
    f.write("! Contact threshold used: 5.0 Å from reference complex\n\n")

    rec_block = fmt_group(sorted(rec_active), "A")

    for p_res in pep_active:
        f.write(f"assign (resi {p_res} and segid P)\n")
        f.write( "       (\n")
        f.write(f"        {rec_block}\n")
        f.write( "       ) 2.0 2.0 0.0\n\n")

    # Reciprocal: one restraint per receptor active residue → all peptide active
    pep_block = fmt_group(pep_active, "P")
    for r_res in sorted(rec_active):
        f.write(f"assign (resi {r_res} and segid A)\n")
        f.write( "       (\n")
        f.write(f"        {pep_block}\n")
        f.write( "       ) 2.0 2.0 0.0\n\n")

    # ── TYR4 ring-atom specific restraints ───────────────────────────────────
    # Per-atom restraints for the TYR4 aromatic ring to guide side chain
    # orientation. Each ring atom is restrained to whichever receptor residues
    # are within CONTACT_DIST of that atom in the reference complex.
    TYR4_RING_ATOMS = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
    f.write("! TYR4 (resi 4) ring-atom restraints — orientation guidance\n\n")
    for atom_name in TYR4_RING_ATOMS:
        ring_sel = u.select_atoms(f"resid 4 and name {atom_name}")
        if len(ring_sel) == 0:
            continue
        d = mda_dist.distance_array(ring_sel.positions, rec_atoms.positions)
        near_resids = sorted({rec_atoms[j].resid
                              for j in range(len(rec_atoms))
                              if d[0, j] <= CONTACT_DIST})
        if not near_resids:
            continue
        near_block = fmt_group(near_resids, "A")
        f.write(f"assign (resi 4 and segid P and name {atom_name})\n")
        f.write( "       (\n")
        f.write(f"        {near_block}\n")
        f.write( "       ) 2.0 2.0 0.0\n\n")
    found_ring_atoms = [a for a in TYR4_RING_ATOMS
                        if len(u.select_atoms(f"resid 4 and name {a}")) > 0]
    print(f"  TYR4 ring atoms restrained: {found_ring_atoms}")

print(f"\nAIR restraints → {tbl_out}")

toml_out = "haddock3_run_TEMPLATE.toml"
with open(toml_out, "w") as f:
    f.write(f"""\
# HADDOCK3 run configuration – peptide-protein docking
# Usage (for each representative, rep 01..N_CLUSTER):
#   haddock3 haddock3_run_repXX.toml
#
# Replace PEPTIDE_FILE with the appropriate peptide_rep_XX.pdb

run_dir  = "run_rep_01"    # change per representative
ncores   = 4               # adjust for your system

molecules = [
    "haddock_inputs/receptor.pdb",
    "haddock_inputs/peptide_rep_01.pdb"   # <<< change per representative
]

[topoaa]
# Topology and coordinates

[rigidbody]
sampling       = 2000
ambig_fname    = "haddock_inputs/ambig.tbl"
w_elec         = 0.2
w_vdw          = 1.0

[seletop]
select = 300

[flexref]
ambig_fname = "haddock_inputs/ambig.tbl"
cool2_steps = 500
cool3_steps = 500

[emref]
# Final energy minimisation

[clustfcc]
fraction_cutoff = 0.60
min_population  = 4

[seletopclusts]
top_models = 4

[caprieval]
# optional, requires reference structure
""")

print(f"HADDOCK3 template → {toml_out}")
print(
    "\nFor N_CLUSTER representatives, create N_CLUSTER copies of this TOML "
    "and change run_dir + molecules accordingly.\n"
    "See also haddock_webserver_guide.md for the web-server workflow."
)
