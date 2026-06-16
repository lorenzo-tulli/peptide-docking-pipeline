# PEP-FOLD3 Web Server Instructions
## Generating 500 peptide conformers for ACE-VPEYINQ

PEP-FOLD3 is a fast _ab initio_ peptide structure prediction server.
Because ETKDG can fail on ACE-capped peptides, PEP-FOLD3 is the recommended
alternative for diverse conformer generation.

**Server URL:** https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/

---

### Peptide information

| Property | Value |
|----------|-------|
| Sequence (1-letter) | VPEYINQ |
| Length | 7 residues |
| N-terminal cap | ACE (acetyl) |
| C-terminal | free (NH₂ in PEP-FOLD, but OK for docking) |

> Note: PEP-FOLD3 uses sequence only – it does not accept the ACE cap.
> Enter the bare 7-residue sequence and re-attach the cap after docking.

---

### Submission protocol (10 × 50 models = 500 total)

PEP-FOLD3 generates up to 50 models per submission. Submit 10 times with
different random seeds to obtain 500 diverse structures.

**Step-by-step:**

1. Navigate to https://bioserv.rpbs.univ-paris-diderot.fr/services/PEP-FOLD3/
2. Paste the sequence: `VPEYINQ`
3. Set **Number of simulations** = 50
4. Change the **Seed** to a unique value for each submission (e.g., 1, 2, …, 10)
5. Click **Submit** and wait (~2-5 min per run)
6. Download all models (ZIP archive)
7. Rename each batch of 50 to `batch_01/` … `batch_10/`
8. Collect all PDB files into `01_conformers/conformers/`
   and rename as `conf_XXXX.pdb` (sequential numbering)

```bash
# Collect and rename after downloading all 10 batches
mkdir -p conformers
n=0
for batch in batch_0*/; do
    for pdb in "$batch"*.pdb; do
        cp "$pdb" conformers/conf_$(printf "%04d" $n).pdb
        n=$((n+1))
    done
done
echo "Total conformers: $n"
```

---

### Alternative: AlphaFold + ESMFold for diverse seeds

For even more conformational diversity, generate models with both:
- **PEP-FOLD3** (structure prediction)
- **ESMFold** (language-model based, https://esmatlas.com/resources?action=fold)

Mix all structures into the `conformers/` directory before clustering.

---

### Next step

Once you have ≥ 100 PDBs in `01_conformers/conformers/`:

```bash
cd ../02_clustering
python rmsd_cluster.py
```
