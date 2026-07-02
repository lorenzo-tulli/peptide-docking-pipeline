# Peptide Docking & Refinement Pipeline

A fully automated SLURM pipeline for peptide–protein docking and refinement on HPC clusters.

---

## Overview

The pipeline docks a short peptide (VPEYINQ, 7 residues) against a kinase dimer receptor
using a geometry-guided conformer selection → HADDOCK3 docking → FlexPepDock refinement
workflow, producing ranked structures ready for AMBER MD simulation.

## Pipeline steps

```
Step 1  generate_rdkit_conformers.py   ETKDG conformers + MMFF94s minimisation
Step 2  rmsd_cluster.py                geometry pre-filter + backbone RMSD clustering → representative structures
Step 3  prepare_haddock_inputs.py      receptor/peptide PDBs + AIR restraints (ambig.tbl)
        generate_haddock3_configs.py   one TOML config per representative
        03_haddock_array.slurm         HADDOCK3 docking
Step 4  collect_and_recluster.py       collect HADDOCK outputs, geometry filter, re-cluster → top poses
Step 5  rank_poses.py                  composite rank: HADDOCK score (60%) + BSA (40%) → top poses
Step 6  prepare_flexpep.py             prepare FlexPepDock inputs
        06_flexpep_array.slurm         PyRosetta FlexPepDock refinement
Step 7  final_ranking.py               rank by total score + interface ΔG
```

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥ 3.10 | all scripts |
| RDKit | ≥ 2023 | conformer generation (ETKDG v3 + MMFF94s) |
| MDAnalysis | ≥ 2.7 | structure parsing, RMSD, geometry checks |
| NumPy / SciPy / pandas | latest | numerics, clustering |
| HADDOCK3 | v2026.5.0 | rigid-body + flexible docking |
| CNS | 1.3 | HADDOCK3 internal engine |
| PyRosetta | ≥ 2024 | FlexPepDock refinement (ref2015) |
| freeSASA | ≥ 2.1 | buried surface area calculation |

---

## Installation

### Conda environment

```bash
conda create -n peptide_docking python=3.10 -c conda-forge -y
conda activate peptide_docking
conda install -c conda-forge rdkit mdanalysis scipy pandas freesasa -y
pip install pyrosetta  # requires PyRosetta wheel from the Rosetta Commons
```

### HADDOCK3

```bash
conda activate peptide_docking
git clone https://github.com/haddocking/haddock3.git
cd haddock3
pip install -e .
```

CNS binaries must be separately obtained (free for non-commercial use):
https://cns-online.org/v1.3/

---

## Usage

### Prepare inputs

Place the following files in `pipeline/01_conformers/`:
- `peptide_reference.pdb` — peptide structure (chain P, residues 1–7)
- `receptor.pdb` — receptor structure (chains A, B, C)

Update the `REFERENCE_COMPLEX` path in `02_clustering/rmsd_cluster.py` and
`03_haddock/prepare_haddock_inputs.py` to point to your full reference complex PDB.

### Run the pipeline

Each step is submitted as a SLURM job. Steps must be run in order:

```bash
cd pipeline

# Step 1 – conformer generation (~2 h)
sbatch slurm/01_conformers.slurm

# Step 2 – clustering (~1 h)
sbatch slurm/02_cluster.slurm

# Step 3 – HADDOCK3 docking (~12 h, 15 parallel jobs)
cd 03_haddock
python prepare_haddock_inputs.py
python generate_haddock3_configs.py
cd ..
sbatch slurm/03_haddock_array.slurm

# Step 4 – collect and re-cluster (~2 h)
sbatch slurm/04_post_haddock.slurm

# Step 5 – scoring (~1 h)
sbatch slurm/05_score.slurm

# --- recommended: inspect top10/ in PyMOL before continuing ---

# Step 6 – FlexPepDock refinement (~24 h, 10 parallel jobs)
cd 06_flexpep && python prepare_flexpep.py && cd ..
sbatch slurm/06_flexpep_array.slurm

# Step 7 – final ranking (~30 min)
sbatch slurm/07_final.slurm
```

## Outputs

| Path | Contents |
|------|----------|
| `01_conformers/conformers/` | 5000 individual conformer PDBs |
| `02_clustering/representatives/` | 15 cluster medoids (HADDOCK inputs) |
| `03_haddock/haddock3_runs/` | HADDOCK3 run directories |
| `04_post_haddock/clustered_poses/` | top 20 geometry-filtered re-clustered poses |
| `05_scoring/top10/` | top 10 poses for FlexPepDock |
| `06_flexpep/flexpep_outputs/` | 200 refined models per pose |
| `07_final/final_top10/` | 10 best structures ready for AMBER MD |

---

## Geometry filters

Two distance-based filters enforce the expected binding mode:

| Contact | Threshold | Applied at |
|---------|-----------|------------|
| TYR4:OH → ASP171:CG | ≤ 7 Å | Step 2 (free conformers, after Kabsch superposition) |
| GLU3:CD → LYS213:NZ | ≤ 4 Å | Step 2 (free conformers, after Kabsch superposition) |
| TYR4:OH → ASP171:CG | ≤ 7 Å | Step 4 (docked poses, direct measurement) |
| GLU3:CD → LYS213:NZ | ≤ 4 Å | Step 4 (docked poses, direct measurement) |

---

## Known limitations

### HIP histidine gets dropped by HADDOCK/CNS

`HID`/`HIE` (neutral histidine tautomers) dock fine — CNS folds them into its own
`HIS` topology. `HIP` (doubly-protonated, positively charged histidine) has no CNS
topology at all, and CNS silently drops any such residue instead of erroring out,
leaving a one-residue gap in the receptor.

`01_conformers/extract_peptide.py` now renames any `HIP` residue to `HIS` before
docking, and prints the residue number(s) it renamed. **If your structure has a
`HIP` residue, you must rename it back from `HIS` to `HIP` at that same residue
number before running `pdb4amber`/`tleap`**, so the final AMBER MD system keeps the
intended doubly-protonated, positively charged histidine instead of a generic
re-protonated one.

---
