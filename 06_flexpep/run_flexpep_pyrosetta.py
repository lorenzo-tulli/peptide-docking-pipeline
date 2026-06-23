#!/usr/bin/env python3
"""
Step 6b FlexPepDock refinement via PyRosetta.

Runs N_MODELS independent refinement trajectories for one HADDOCK pose,
writes each model as a PDB, and records scores in a CSV.

FlexPepDock is a Rosetta protocol specifically designed for peptide-protein docking and refinement. 
Unlike HADDOCK which uses CNS/PARALLHDG, Rosetta uses its own Rosetta energy function (ref2015) which is a knowledge-based + physics-based hybrid.
It starts from the HADDOCK pose, then:
 -Low-resolution preoptimisation to sample the optimal backbone conformation
 -High-resolution Monte Carlo full atom representation to sample both backbone and sidechain conformation
 -Accept/reject by Metropolis criterion

Called by slurm/06_flexpep_array.slurm (one SLURM array task per pose).

Dependencies:
  pyrosetta (installed via conda), pandas

Usage:
  python run_flexpep_pyrosetta.py --pose_id 01 --n_models 200

Input:  flexpep_inputs/fpd_pose_01.pdb   (receptor chain A, peptide chain B)
Output: flexpep_outputs/pose_01/model_XXXX.pdb
        flexpep_outputs/pose_01/scores.csv
"""

import os, sys, time, argparse, json
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--pose_id",    required=True,     help="Zero-padded pose ID, e.g. 01")
ap.add_argument("--n_models",   type=int, default=200)
ap.add_argument("--input_dir",  default="flexpep_inputs")
ap.add_argument("--output_dir", default="flexpep_outputs")
args = ap.parse_args()

POSE_ID   = args.pose_id
N_MODELS  = args.n_models
INPUT_PDB = os.path.abspath(os.path.join(args.input_dir, f"fpd_pose_{POSE_ID}.pdb"))
OUT_DIR   = os.path.abspath(os.path.join(args.output_dir, f"pose_{POSE_ID}"))

if not os.path.isfile(INPUT_PDB):
    sys.exit(f"ERROR: input PDB not found: {INPUT_PDB}\n"
             "Run prepare_flexpep.py first.")

os.makedirs(OUT_DIR, exist_ok=True)

# ── initialise PyRosetta (must happen before any other pyrosetta calls) ───────
print(f"Initialising PyRosetta …", flush=True)
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import FlatHarmonicFunc
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core import scoring as rosetta_scoring

init(extra_options=" ".join([
    "-ex1", "-ex2aro",
    "-use_input_sc",
    "-detect_disulf false",
    "-flexPepDocking:receptor_chain A",
    "-flexPepDocking:peptide_chain B",
    "-flexPepDocking:pep_refine true",
    "-lowres_preoptimize true",
    "-score:weights ref2015",
    # suppress routine Rosetta verbosity; remove -mute if you need debug output
    "-mute protocols.moves core.pack core.scoring core.optimization",
]), silent=True)

sfxn = ScoreFunctionFactory.create_score_function("ref2015")
sfxn.set_weight(rosetta_scoring.atom_pair_constraint, 2.0)
fpd  = FlexPepDockingProtocol()

_TYR4_RING_ATOMS = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]

def add_tyr4_ring_constraints(pose):
    """Flat-harmonic distance constraints from each TYR4 ring atom to its
    closest receptor heavy atom in the HADDOCK input structure.
    Flat region ±1 Å, harmonic sd=0.5 Å outside — keeps the ring in the
    HADDOCK-predicted orientation throughout FlexPepDock sampling."""
    tyr4_rosnum = pose.pdb_info().pdb2pose("B", 4)
    if tyr4_rosnum == 0:
        print("  WARNING: TYR4 (chain B resid 4) not found – no ring constraints")
        return 0
    rec_end = pose.chain_end(1)
    n = 0
    for atom_name in _TYR4_RING_ATOMS:
        res = pose.residue(tyr4_rosnum)
        if not res.has(atom_name):
            continue
        pep_id  = AtomID(res.atom_index(atom_name), tyr4_rosnum)
        pep_xyz = res.atom(atom_name).xyz()
        min_d, best_id = 999.0, None
        for rn in range(1, rec_end + 1):
            rres = pose.residue(rn)
            for ai in range(1, rres.natoms() + 1):
                if rres.atom_is_hydrogen(ai):
                    continue
                d = pep_xyz.distance(rres.atom(ai).xyz())
                if d < min_d:
                    min_d, best_id = d, AtomID(ai, rn)
        if best_id and min_d < 8.0:
            pose.add_constraint(
                AtomPairConstraint(pep_id, best_id,
                                   FlatHarmonicFunc(min_d, 0.5, 1.0))
            )
            n += 1
    return n


def add_glu3_lys213_constraint(pose):
    """Flat-harmonic constraint between GLU3:CD (peptide) and LYS213:NZ (receptor).
    Holds the salt-bridge distance observed in the HADDOCK input pose."""
    glu3_rosnum  = pose.pdb_info().pdb2pose("B", 3)
    lys213_rosnum = pose.pdb_info().pdb2pose("A", 213)
    if glu3_rosnum == 0:
        print("  WARNING: GLU3 (chain B resid 3) not found – skipping constraint")
        return 0
    if lys213_rosnum == 0:
        print("  WARNING: LYS213 (chain A resid 213) not found – skipping constraint")
        return 0
    glu3_res  = pose.residue(glu3_rosnum)
    lys213_res = pose.residue(lys213_rosnum)
    if not glu3_res.has("CD"):
        print("  WARNING: GLU3 has no CD atom – skipping constraint")
        return 0
    if not lys213_res.has("NZ"):
        print("  WARNING: LYS213 has no NZ atom – skipping constraint")
        return 0
    glu3_id  = AtomID(glu3_res.atom_index("CD"),   glu3_rosnum)
    lys213_id = AtomID(lys213_res.atom_index("NZ"), lys213_rosnum)
    d = glu3_res.atom("CD").xyz().distance(lys213_res.atom("NZ").xyz())
    pose.add_constraint(
        AtomPairConstraint(glu3_id, lys213_id, FlatHarmonicFunc(d, 0.5, 1.0))
    )
    return 1


# Optional: interface analyser for I_dG (graceful fallback if unavailable)
try:
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    _iam_available = True
except ImportError:
    _iam_available = False

def interface_dG(pose):
    if not _iam_available:
        return None
    try:
        iam = InterfaceAnalyzerMover()
        iam.set_interface("A_B")
        iam.set_pack_separated(True)
        iam.set_scorefunction(sfxn)
        iam.apply(pose)
        return round(iam.get_interface_dG(), 3)
    except Exception:
        return None

print(f"=== FlexPepDock – pose {POSE_ID} ===")
print(f"Input:   {INPUT_PDB}")
print(f"Output:  {OUT_DIR}")
print(f"Models:  {N_MODELS}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n", flush=True)

rows = []
for i in range(1, N_MODELS + 1):
    t0 = time.time()

    # Fresh pose each model so MC starts from the same input geometry
    pose = pose_from_pdb(INPUT_PDB)
    n_tyr4 = add_tyr4_ring_constraints(pose)
    n_salt = add_glu3_lys213_constraint(pose)
    if i == 1:
        print(f"  TYR4 ring constraints added:      {n_tyr4}", flush=True)
        print(f"  GLU3:CD–LYS213:NZ constraint added: {n_salt}", flush=True)

    fpd.apply(pose)

    total  = round(sfxn(pose), 3)
    i_dg   = interface_dG(pose)

    out_pdb = os.path.join(OUT_DIR, f"model_{i:04d}.pdb")
    pose.dump_pdb(out_pdb)

    elapsed = time.time() - t0
    rows.append({
        "pose_id":     POSE_ID,
        "model":       i,
        "pdb":         out_pdb,
        "total_score": total,
        "interface_dG": i_dg,
        "time_s":      round(elapsed, 1),
    })

    print(f"  [{i:3d}/{N_MODELS}]  total={total:9.2f}  "
          f"I_dG={str(i_dg):>9}  ({elapsed:.0f}s)", flush=True)

df = pd.DataFrame(rows)
csv_path = os.path.join(OUT_DIR, "scores.csv")
df.to_csv(csv_path, index=False)

print(f"\nScores → {csv_path}")
print(f"Best total_score : {df['total_score'].min():.2f}  "
      f"(model {df.loc[df['total_score'].idxmin(), 'model']:04d})")
if df["interface_dG"].notna().any():
    print(f"Best interface_dG: {df['interface_dG'].min():.2f}")
print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
