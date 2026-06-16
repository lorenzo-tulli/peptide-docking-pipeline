#!/usr/bin/env python3
"""
PyRosetta installation test.
Run this before submitting any SLURM jobs to confirm everything works.

Usage:
  conda activate pep_docking
  python pipeline/06_flexpep/test_pyrosetta.py

Expected output: 4/4 tests passed
"""

import sys

PASS = []
FAIL = []

def ok(msg):
    print(f"  PASS  {msg}")
    PASS.append(msg)

def fail(msg, err):
    print(f"  FAIL  {msg}")
    print(f"        → {err}")
    FAIL.append(msg)

# ── Test 1: basic import ──────────────────────────────────────────────────────
print("\n[1/4] Import pyrosetta …")
try:
    import pyrosetta
    ok(f"pyrosetta imported (version: {pyrosetta.version()})")
except Exception as e:
    fail("import pyrosetta", e)
    print("\nCannot continue without a working import. Check your installation.")
    sys.exit(1)

# ── Test 2: initialise (loads the Rosetta database) ──────────────────────────
print("\n[2/4] Initialise PyRosetta (loads database, ~10 s) …")
try:
    pyrosetta.init(
        extra_options="-mute all",
        silent=True
    )
    ok("init() completed successfully")
except Exception as e:
    fail("pyrosetta.init()", e)

# ── Test 3: build a pose from sequence and score it ───────────────────────────
print("\n[3/4] Build peptide pose from sequence and score …")
try:
    from pyrosetta import pose_from_sequence
    from pyrosetta.rosetta.core.scoring import ScoreFunctionFactory

    pose = pose_from_sequence("VPEYINQ", "fa_standard")
    sfxn = ScoreFunctionFactory.create_score_function("ref2015")
    score = sfxn(pose)
    ok(f"Pose built ({pose.total_residue()} residues), total_score = {score:.2f}")
except Exception as e:
    fail("pose_from_sequence / scoring", e)

# ── Test 4: FlexPepDockingProtocol can be imported ───────────────────────────
print("\n[4/4] Import FlexPepDockingProtocol …")
try:
    from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
    fpd = FlexPepDockingProtocol()
    ok("FlexPepDockingProtocol instantiated")
except Exception as e:
    fail("FlexPepDockingProtocol", e)
    print("\n  If this fails, try reinstalling PyRosetta with:")
    print("  conda install pyrosetta -c <your-channel>")

# ── Summary ───────────────────────────────────────────────────────────────────
total = len(PASS) + len(FAIL)
print(f"\n{'='*50}")
print(f"  {len(PASS)}/{total} tests passed")

if FAIL:
    print(f"\n  Failed tests:")
    for f in FAIL:
        print(f"    - {f}")
    print("\n  PyRosetta is NOT ready. Fix the errors above before running the pipeline.")
    sys.exit(1)
else:
    print("\n  PyRosetta is ready. Run the pipeline with:")
    print("  bash pipeline/slurm/submit_pipeline.sh")
    print(f"{'='*50}\n")
