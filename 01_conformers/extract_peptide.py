#!/usr/bin/env python3
"""
Step 1a Split INH_PROTEIN.pdb into:
  peptide_reference.pdb  : VPEYINQ (residues 1-7), chain P
  receptor.pdb           : protein dimer (chains A, B) + Mg-ATP (chain C)

The original PDB has no chain IDs (AMBER output). Chain assignment:
  residues  1-7   peptide  (chain P)
  residues  8-352 chain A  (protein monomer 1)
  residues 353-697 chain B  (protein monomer 2)
  residues 698-699 chain C  (MG, ATP cofactors)

Usage:
  cd pipeline/01_conformers
  python extract_peptide.py
"""

INPUT_PDB = "/group/chem/oliveira2t/vi24769/EGFR_STRING/system-building/equilibration_protocol/1Mg/peptide_docking_pipeline/WT-DIMER_peptide_substrate_ATP-Mg.pdb"

PEPT_RANGE   = range(1, 8)        # residues 1-7
CHAIN_A_RANGE = range(8, 353)     # monomer 1
CHAIN_B_RANGE = range(353, 698)   # monomer 2
COFACTOR_RESIDS = {698, 699}      # MG, ATP

def chain_for(resnum):
    if resnum in PEPT_RANGE:
        return "P"
    if resnum in CHAIN_A_RANGE:
        return "A"
    if resnum in CHAIN_B_RANGE:
        return "B"
    if resnum in COFACTOR_RESIDS:
        return "C"
    return "X"

peptide_lines  = []
receptor_lines = []

with open(INPUT_PDB) as fh:
    for line in fh:
        rec = line[:6].strip()
        if rec not in ("ATOM", "HETATM"):
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        ch = chain_for(resnum)
        new = line[:21] + ch + line[22:]
        if ch == "P":
            peptide_lines.append(new)
        else:
            receptor_lines.append(new)

with open("peptide_reference.pdb", "w") as fh:
    fh.writelines(peptide_lines)
    fh.write("TER\nEND\n")

with open("receptor.pdb", "w") as fh:
    fh.writelines(receptor_lines)
    fh.write("TER\nEND\n")

print(f"peptide_reference.pdb : {len(peptide_lines):4d} atoms  (chain P, residues 1-7)")
print(f"receptor.pdb          : {len(receptor_lines):4d} atoms  (chains A B C)")
