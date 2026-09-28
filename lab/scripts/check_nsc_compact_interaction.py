#!/usr/bin/env python3
"""Reproduce the compact-mode contact map, including its fermionic vertex."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_compact_interaction import (
    WARP_COEFFICIENT, clifford_channels, fermionic_vertex,
    independent_overlap, profile_overlap, quadrature, vertex,
)

OUTPUT = ROOT / "results/development/compact-interaction.json"
SOURCES = (
    "src/recursive_horizons/nsc_compact_interaction.py",
    "src/recursive_horizons/nsc_spinor_bridge.py",
    "src/recursive_horizons/nsc_boundary.py",
    "scripts/check_nsc_compact_interaction.py",
    "docs/nsc-compact-interaction.md",
)
INPUTS = ("docs/nsc-compact-mass-map.md",)
ATOL = RTOL = 3e-12


def hashes(paths):
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def compare(expected, actual, path="$"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for k in expected:
            compare(expected[k], actual[k], f"{path}/{k}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare(a, b, f"{path}/{i}")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (float, int)) or not np.isclose(expected, actual, atol=ATOL, rtol=RTOL):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def calculate():
    cases = (
        ((0, 0, 0, 0), "LLLL", 1.),
        ((1, 1, 1, 1), "LLLL", 1.5),
        ((1, 1, 1, 1), "RRRR", 1.5),
        ((1, 1, 1, 1), "LLRR", .5),
        ((0, 0, 1, 1), "LLLL", 1.),
        ((0, 0, 1, 1), "LLRR", 1.),
        ((3, 1, 1, 1), "LLLL", .5),
        ((3, 1, 1, 1), "RRRR", -.5),
    )
    overlaps = []
    for levels, chiralities, exact_flat in cases:
        flat = profile_overlap(levels, chiralities, warp=0.)
        fixed = profile_overlap(levels, chiralities)
        adaptive = independent_overlap(levels, chiralities)
        if abs(flat - exact_flat) > 2e-12 or abs(fixed - adaptive) > 2e-12:
            raise AssertionError("overlap integration mismatch")
        overlaps.append({"levels": list(levels), "chiralities": chiralities,
                         "flat_exact": exact_flat, "flat_quadrature": flat,
                         "warped": fixed, "adaptive": adaptive,
                         "fixed_minus_adaptive": fixed - adaptive})
    data = clifford_channels()
    triple_matrix = sum((c * np.kron(g, g) for c, g in data["triples"]), np.zeros((16, 16), complex))
    decomposed_matrix = sum((c * np.kron(g, g) for c, g in data["axial"] + data["tensor"]), np.zeros((16, 16), complex))
    if not np.array_equal(triple_matrix, decomposed_matrix):
        raise AssertionError("published Clifford identity in project signature")
    witnesses = []
    ll_index = (2, 0, 3, 1)
    tensor_index = (0, 0, 1, 1)
    for levels in ((0, 0, 0, 0), (1, 1, 1, 1), (3, 1, 1, 1), (0, 0, 1, 1), (0, 1, 0, 1), (0, 1, 1, 1)):
        raw = vertex(levels)
        split = vertex(levels, channels="decomposed")
        anti = fermionic_vertex(levels)
        tensor = fermionic_vertex(levels, channels="tensor")
        error = float(np.max(np.abs(raw - split)))
        if error > 2e-12 or np.max(np.abs(anti.imag)) > 2e-12:
            raise AssertionError("dressed vertex mismatch")
        ll = float(anti[ll_index].real)
        if levels == (3, 1, 1, 1) and abs(ll) < 1:
            raise AssertionError("missing cross-level interaction")
        if levels == (0, 1, 1, 1) and np.max(np.abs(anti)) > 2e-12:
            raise AssertionError("symmetric-interval parity selection")
        if levels == (0, 0, 0, 0) and np.max(np.abs(tensor)) != 0:
            raise AssertionError("pure Weyl zero mode has no tensor bilinear")
        witnesses.append({"levels": list(levels), "direct_minus_decomposed_max": error,
                          "fermionic_left_current_component": ll,
                          "fermionic_tensor_component": float(tensor[tensor_index].real),
                          "fermionic_max_absolute": float(np.max(np.abs(anti)))})
    convergence = [{"points": n,
                    "overlap_3111_LLLL": profile_overlap((3, 1, 1, 1), "LLLL", points=n),
                    "overlap_1111_LLRR": profile_overlap((1, 1, 1, 1), "LLRR", points=n)}
                   for n in (32, 64, 96)]
    return {
        "schema": "NSC-COMPACT-INTERACTION-MAP-v1",
        "status": "development calculation; not part of the v0.3.0 release collection",
        "sources": ["https://arxiv.org/html/1405.0397v1", "https://arxiv.org/html/1512.06074v1"],
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "domain": {"coordinate_interval": [-1., 1.], "dimensional_coordinate": "Y=L_star*y",
                   "sigma_coefficient": WARP_COEFFICIENT, "sigma": "coefficient*y^2",
                   "boundary": "P_R chi=0 at both endpoints; chosen fundamental interval",
                   "state": "not specified; this is an action vertex"},
        "normalization": {"kappa4_bulk_squared": "kappa5_squared/I3", "I3": "integral exp(3*sigma) dY",
                          "I3_over_L_star": quadrature()[2],
                          "vertex_prefactor_omitted": "kappa4_bulk_squared/(32*s_T)",
                          "s_T": None, "s_T_meaning": "five-dimensional torsion stiffness relative to Einstein-Cartan; EC is1; NSC value not selected",
                          "overlap": "I3*integral exp(-3*sigma)*f_i*f_j*f_k*f_l dY",
                          "overlap_units": "dimensionless",
                          "fermionic_projection": "minus one quarter times antisymmetrization in barred and unbarred mode/spin legs"},
        "overlaps": overlaps, "quadrature_convergence": convergence,
        "left_current_witness_axes_a_b_c_d": list(ll_index),
        "tensor_witness_axes_a_b_c_d": list(tensor_index),
        "vertex_witnesses": witnesses,
        "conclusions": {"generic_single_massive_mode_exact_truncation": False,
                        "odd_bilinear_implies_zero_contact_square": False,
                        "mode_sum_odd_vertex_vanishes_for_this_symmetric_bulk": True,
                        "current_sign_determines_stress_or_cosmic_acceleration": False},
        "unresolved": ["five-dimensional common-functional stiffness and derivative terms",
                       "compact boundary action and total gravitational normalization",
                       "physical compact size and species", "state and contour propagators",
                       "renormalized absolute stress and metric backreaction"],
        "comparison": {"fields": "all", "atol": ATOL, "rtol": RTOL,
                       "exact": "keys, types, source/input hashes; authenticate dependencies before reproduction",
                       "exceptions": []},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite a recorded result")
    expected = None
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        compare(expected["source_hashes"], hashes(SOURCES), "$/source_hashes")
        compare(expected["input_hashes"], hashes(INPUTS), "$/input_hashes")
    result = calculate()
    if args.check:
        compare(expected, result)
        print("compact interaction: authenticated dependencies, all fields reproduced")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
