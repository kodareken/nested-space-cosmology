#!/usr/bin/env python3
"""Inspect authenticated demonstrations; use --recompute to execute them again."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

from reproduce_public_results import load_manifest, validate_checkout

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ("results/nsc-3-boundary-response.json", "scripts/check_nsc_boundary_response.py",
     "Geometry to boundary response",
     "Finite spatial maps; physical Lorentzian matching and the stationary scale remain open."),
    ("results/nsc-6-vacuum-work.json", "scripts/check_nsc_vacuum_work.py",
     "Geometric work to Dirac excitations",
     "The radius pulse is prescribed; its self-consistent gravitational source remains open."),
    ("results/nsc-7-observable-bridge.json", "scripts/check_nsc_observable_bridge.py",
     "A scalar sheet link and the Dirac mass equation",
     "The actual throat coupling, physical sector and particle identification remain open."),
    ("results/nsc-10-influence.json", "scripts/check_nsc_influence.py",
     "Normalized histories and force fluctuations",
     "Finite Gaussian-state construction; absolute metric-source matching remains open."),
    ("results/nsc-11-response-matching.json", "scripts/check_nsc_response_matching.py",
     "Static and causal geometric response",
     "Relative ultrastatic radius response; the complete finite-cutoff metric theory remains open."),
    ("results/development/compact-interaction.json", "scripts/check_nsc_compact_interaction.py",
     "Interactions between compact modes",
     "A candidate bulk contact; overall coupling, boundary action and quantum state remain open."),
    ("results/development/charged-sector.json", "scripts/check_nsc_charged_sector.py",
     "A consistent charged compact sector",
     "Candidate field content and parity; gauge representation, link mass and self-sourcing remain open."),
)


def headlines(relative: str, record: dict) -> list[str]:
    if relative.endswith("nsc-3-boundary-response.json"):
        return [f"Joined and eliminated responses agree: {record['gate']['joined_schur_identity_holds']}",
                f"Independent half-domain maps resolved: {record['gate']['discrete_continuum_maps_resolved']}"]
    if relative.endswith("nsc-6-vacuum-work.json"):
        clock_pairs = record["pure_clock_control"]["pair_number_coefficient"]
        return [record["work_identity"],
                f"Uniform-clock control pair coefficient: {clock_pairs:.3g} (recorded numerical residual)."]
    if relative.endswith("nsc-7-observable-bridge.json"):
        sector = record["spinor_algebra"]["candidate_restriction"]
        return [f"Invariant projector: {sector['projector']}",
                f"Sector rank: {sector['rank']}; complementary sector rank: {sector['complement_rank']}"]
    if relative.endswith("nsc-10-influence.json"):
        equal = record["equal_history"]["amplitude"]
        smeared = record["smeared_generator"]
        return [f"Equal-history amplitude: {equal['real']:.8g} + ({equal['imag']:.3g})i",
                f"Smeared force variance: {smeared['variance']:.8g}, width {smeared['width']:g} (declared units)."]
    if relative.endswith("nsc-11-response-matching.json"):
        gate = record["gate"]
        return [f"Relative Euclidean/retarded limits agree: {gate['relative_euclidean_and_Abel_limit_match']}",
                f"State-dependent coordinate contact retained: {gate['state_dependent_coordinate_contact_retained']}"]
    if relative.endswith("compact-interaction.json"):
        overlap = next(row for row in record["overlaps"]
                       if row["levels"] == [3, 1, 1, 1] and row["chiralities"] == "LLLL")
        return [f"Geometric 3111/LLLL overlap: {overlap['warped']:.8g}",
                "The full fermionic vertex retains this cross-level interaction; one massive level is not an exact closed sector."]
    if relative.endswith("charged-sector.json"):
        allowed = [row for row in record["parity_and_anomaly_rows"]
                   if row["anomaly_free_zero_sector"] and row["neutral_even_scalar_link_allowed"]]
        pairs = [(row["parent_parity_sign"], row["child_parity_sign"]) for row in allowed]
        multiplicity = record["multiplicity"]
        return [f"Allowed opposite-parity pairs in this minimal candidate: {pairs}",
                f"Bulk Dirac copies: {multiplicity['bulk_complex_Dirac_fields']}; four-dimensional Dirac zero fields before link mass: {multiplicity['massless_4D_Dirac_zero_fields_before_link_mass']}"]
    raise ValueError(f"No demonstration summary for {relative}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute", action="store_true",
                        help="rerun the seven selected generators with their all-field checks")
    args = parser.parse_args()
    manifest = load_manifest()
    validate_checkout(manifest)
    declared = {step["output"] for step in manifest["steps"]}
    missing = [path for path, *_ in CASES if path not in declared]
    if missing:
        raise RuntimeError(f"Demonstrations are missing from the release manifest: {missing}")
    env = os.environ.copy()
    env.update({key: "1" for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
                                     "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
    env["PYTHONHASHSEED"] = "0"
    print("Authenticated research demonstrations", flush=True)
    if not args.recompute:
        print("Displaying recorded results; no scientific generators are being rerun.", flush=True)
    for number, (relative, generator, title, scope) in enumerate(CASES, start=1):
        if args.recompute:
            subprocess.run([sys.executable, str(ROOT / generator), "--check"],
                           cwd=ROOT, env=env, check=True)
        record = json.loads((ROOT / relative).read_text())
        print(f"\n{number}. {title}")
        for line in headlines(relative, record):
            print(f"   {line}")
        print(f"   Scope: {scope}")
        print(f"   Record: {relative}", flush=True)
    action = "recomputed and compared" if args.recompute else "displayed from authenticated records"
    print(f"\n{len(CASES)} demonstrations {action} within the {manifest['result_count']}-record release.")
    if not args.recompute:
        print("To execute these checks again: python3 scripts/demonstrate.py --recompute")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
