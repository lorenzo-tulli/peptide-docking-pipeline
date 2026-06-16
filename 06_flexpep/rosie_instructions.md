# FlexPepDock via ROSIE Web Server
## Refining top 5 HADDOCK poses (200 models per pose)

**Server URL:** https://rosie.rosettacommons.org/flexpep_docking

Free academic access — requires ROSIE account registration.

---

## Input file requirements

Each input PDB must:
- Contain **exactly two chains**: A (receptor) and B (peptide)
- Be free of cofactors and solvent (MG, ATP removed — add back for MD)
- Have correct peptide placement in the binding site (from HADDOCK)

The script `prepare_flexpep.py` produces ready-to-upload PDBs in `flexpep_inputs/`.

---

## Submission protocol (repeat for each of the 5 poses)

**1. Upload the PDB**
- Go to https://rosie.rosettacommons.org/flexpep_docking
- Under **"Input"**, upload `flexpep_inputs/fpd_pose_01.pdb`

**2. Set run parameters**

| Parameter | Recommended value | Notes |
|-----------|-------------------|-------|
| Number of output structures | **200** | full refinement |
| Receptor chain | **A** | protein dimer |
| Peptide chain | **B** | ACE-VPEYINQ |
| Refinement mode | **pep_refine** | full-atom refinement |
| Pre-optimize (low-res) | **ON** | |
| Extra chi sampling | **ex1, ex2aro** | |
| Score function | **ref2015** | |

**3. Submit and record the job ID**

**4. Download results when complete (~30 min per run)**
- Download `score.sc` (score file)
- Download `silent.out` or individual PDB models

**5. Repeat for all 5 poses**

---

## Expected output structure (per pose)

```
rosie_job_XXXXXX/
├── score.sc           ← Rosetta score file (key output)
├── S_0001.pdb ... S_0200.pdb   ← refined models
└── logs/
```

Organise downloads as:
```
06_flexpep/
  rosie_outputs/
    pose_01/
      score.sc
      S_0001.pdb ... S_0200.pdb
    pose_02/ ...
```

---

## Key score terms in score.sc

| Column | Meaning |
|--------|---------|
| `total_score` | Overall Rosetta energy (lower = better) |
| `pep_score` | Peptide-specific score |
| `I_sc` | Interface score |
| `reweighted_sc` | Reweighted score (recommended for final ranking) |

---

## Next step

```bash
cd ../../07_final
python final_ranking.py
```

---

## Alternative: Local Rosetta installation

If you have Rosetta installed locally:

```bash
# Load flags file (from prepare_flexpep.py)
cd pipeline/06_flexpep/flexpep_inputs
bash run_all_flexpep.sh
```

Rosetta download requires academic license registration at:
https://els2.comotion.uw.edu/product/rosetta

For HPC clusters (SLURM):
```bash
sbatch --ntasks=4 --time=4:00:00 \
    --wrap="FlexPepDocking.linuxgccrelease @fpd_pose_01.flags"
```
