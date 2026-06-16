# HADDOCK Web Server Guide
## Peptide–protein docking: ACE-VPEYINQ × kinase dimer

---

## Which server to use?

| Server | URL | Limit | Best for |
|--------|-----|-------|----------|
| **HADDOCK 2.4** (recommended) | https://wenmr.science.uu.nl/haddock2.4/ | 2 molecules, 50 models/cluster | Academic free access |
| **HADDOCK 3** (local/HPC) | https://github.com/haddocking/haddock3 | unlimited | Large-scale / automated |

Both require a **free ORCA/WeNMR account**: https://wenmr.science.uu.nl/

---

## HADDOCK 2.4 – Step-by-step (repeat for each of the 10 representatives)

### Preparation (already done by prepare_haddock_inputs.py)

```
haddock_inputs/
├── receptor.pdb          ← protein dimer (chains A, B) + Mg/ATP (chain C)
├── peptide_rep_01.pdb    ← representative 1 (chain P)
├── ...
├── peptide_rep_10.pdb    ← representative 10 (chain P)
└── ambig.tbl             ← AIR restraints
```

### Submission (one run per representative)

**1. Login and start a new docking run**
- Go to https://wenmr.science.uu.nl/haddock2.4/
- Click **"Start a new docking run"**

**2. Molecule 1 – Receptor**
- Upload `haddock_inputs/receptor.pdb`
- Segid: `A` (chain A) and `B` (chain B)
- Active residues: paste the chain A list from prepare_haddock_inputs.py output
- Passive residues: leave blank (or let HADDOCK auto-define surface neighbours)
- **Check** "Define passive residues automatically around the active residues"

**3. Molecule 2 – Peptide**
- Upload `haddock_inputs/peptide_rep_01.pdb`  (change for each run)
- Segid: `P`
- Active residues: paste the peptide list from prepare_haddock_inputs.py output
- Passive residues: leave blank

**4. AIR restraints**
- Under **"Distance restraints"** → **"Ambiguous restraints"**
- Upload `haddock_inputs/ambig.tbl`

**5. Docking parameters** (under "Advanced parameters")

| Parameter | Value | Reason |
|-----------|-------|--------|
| Structures for rigid-body | 1000 | broad sampling |
| Structures for semi-flexible | 200 | |
| Structures for water refinement | 200 | |
| Clustering method | RMSD | |
| RMSD cutoff | 2.0 Å | tight for 7-mer |
| Number of structures to analyse | 200 | |
| Peptide flexibility | **ON** | peptide_flexibility = true |

**6. Submit and wait (~30–60 min per run)**

**7. Download results**
- Download the full archive (`.tgz`)
- Extract to `../04_post_haddock/haddock_runs/run_rep_01/` … `run_rep_10/`

### Repeat for all 10 representatives

You can submit all 10 simultaneously (server allows multiple queued jobs).

---

## HADDOCK 3 – Local/HPC run (batch all 10 at once)

```bash
# Install (once)
pip install haddock3

# Run all 10 representatives in sequence
for i in $(seq -w 1 10); do
    cp haddock3_run_TEMPLATE.toml haddock3_run_rep${i}.toml
    # Edit run_dir and molecules in the TOML
    sed -i "s/run_rep_01/run_rep_${i}/g" haddock3_run_rep${i}.toml
    sed -i "s/peptide_rep_01/peptide_rep_${i}/g" haddock3_run_rep${i}.toml
    haddock3 haddock3_run_rep${i}.toml &
done
wait
echo "All HADDOCK3 runs complete"
```

---

## What to expect in the output

Each run produces:
```
run_rep_XX/
├── it0/              ← rigid-body docking (1000 structures)
├── it1/              ← semi-flexible refinement (200)
├── water/            ← explicit water refinement (200)
└── analysis/
    ├── cluster_*.pdb      ← cluster representatives
    └── haddock_score.dat  ← HADDOCK score breakdown
```

The **HADDOCK score** (lower is better) is:
```
HADDOCK score = 1.0×Evdw + 0.2×Eelec + 0.1×Edesol − 0.01×BSA + 1.0×Eair
```

---

## Next step

After downloading all 10 run archives:
```bash
cd ../04_post_haddock
python collect_and_recluster.py
```
