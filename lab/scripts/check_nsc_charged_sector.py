#!/usr/bin/env python3
"""Check the minimal charged compact sector and its UV gauge normalization.

This adds a candidate U(1) connection to the two-copy free carrier. It does
not derive the gauge group, the scalar link, a state or a self-sourced metric.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_spinor_bridge import pauli, weyl_matrices
from recursive_horizons.nsc_covariant_identities import _matrices

OUTPUT = ROOT / "results/development/charged-sector.json"
SOURCES = ("scripts/check_nsc_charged_sector.py", "src/recursive_horizons/nsc_spinor_bridge.py",
           "src/recursive_horizons/nsc_covariant_identities.py", "docs/nsc-charged-self-sourcing-route.md")


def hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCES}


def calculate():
    w = weyl_matrices()
    beta, chirality = w["gamma"][0], w["gamma5"]
    tau1, _, tau3 = pauli()
    kp = sp.kronecker_product
    full_chirality = kp(sp.eye(2), chirality)
    alpha_y = kp(sp.eye(2), sp.I * beta * chirality)
    neutral_even_link = kp(tau1, beta)
    rows = []
    for parent, child in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        sheet_parity = sp.diag(parent, child)
        reflection = -kp(sheet_parity, chirality)
        even = (sp.eye(8) + reflection) / 2
        if sp.simplify(reflection * alpha_y + alpha_y * reflection) != sp.zeros(8):
            raise AssertionError("compact reflection of first-order derivative")
        if even * alpha_y * even != sp.zeros(8):
            raise AssertionError("compact current boundary form")
        left = int(sp.trace(even * kp(sp.eye(2), w["left"])) / 2)
        right = int(sp.trace(even * kp(sp.eye(2), w["right"])) / 2)
        # A right-handed charge-one field counts as left-handed charge minus one.
        linear_anomaly = -sp.trace(even * full_chirality) / 2
        cubic_anomaly = linear_anomaly  # both physical fields have unit charge
        mass_allowed = reflection * neutral_even_link == neutral_even_link * reflection
        anomaly_free = linear_anomaly == cubic_anomaly == 0
        if mass_allowed != (parent != child) or anomaly_free != mass_allowed:
            raise AssertionError("minimal charged-domain conditions")
        if parent != child:
            target = (sp.eye(8) - parent * kp(tau3, chirality)) / 2
            if even != target:
                raise AssertionError("previous invariant-sector projector")
        rows.append({"parent_parity_sign": parent, "child_parity_sign": child,
                     "left_Weyl_zero_modes": left, "right_Weyl_zero_modes": right,
                     "mixed_gravitational_U1_anomaly": str(linear_anomaly),
                     "cubic_U1_anomaly": str(cubic_anomaly),
                     "neutral_even_scalar_link_allowed": mass_allowed,
                     "anomaly_free_zero_sector": anomaly_free})
    # Universal a4 formula, evaluated only to match charge/rank conventions.
    # One unit F12 component has F_ab F_ab=2 in the Euclidean frame.
    gamma = _matrices()[0]
    e_linear = sp.I * gamma[0] * gamma[1]
    e_part = sp.trace(e_linear * e_linear) / 2
    omega_part = 2 * sp.trace((-sp.I * sp.eye(4)) ** 2) / 12
    per_dirac = sp.simplify((e_part + omega_part) / 2)
    if per_dirac != sp.Rational(2, 3):
        raise AssertionError("unit-charge Dirac heat coefficient")
    return {
        "schema": "NSC-CHARGED-SECTOR-v1",
        "status": "candidate operator/domain matching; not a new wormhole construction",
        "source_hashes": hashes(),
        "sources": ["https://arxiv.org/html/1807.04726v3", "https://arxiv.org/html/hep-ph/9912408v3",
                    "https://arxiv.org/html/hep-th/0306138v3"],
        "assumptions": ["two five-dimensional Dirac copies on the same four-dimensional development background",
                        "unit U(1) charge on both copies; neutral off-diagonal scalar coupling",
                        "room-diagonal compact reflection with signs plus or minus one",
                        "even canonical scalar-link profile for the parity test",
                        "same declared length-2 compact interval; its global extension remains open"],
        "parity_and_anomaly_rows": rows,
        "allowed_zero_projectors": ["(I8-tau3 tensor gamma5)/2", "(I8+tau3 tensor gamma5)/2"],
        "multiplicity": {"bulk_complex_Dirac_fields": 2, "bulk_spinor_rank": 8,
                         "massless_4D_Dirac_zero_fields_before_link_mass": 1,
                         "massive_4D_Dirac_fields_per_nonzero_compact_level": 2},
        "gauge_heat_coefficient": {"per_unit_charge_Dirac_a4_F_squared": str(per_dirac),
                                   "two_copy_a4_F_squared": str(2 * per_dirac),
                                   "heat_prefactor_excluded": "(4*pi*t)^(-5/2)",
                                   "proper_time_window_weight": "Lambda-nu_match",
                                   "two_copy_bulk_Maxwell_coefficient": "4*(Lambda-nu_match)/(3*(4*pi)^(5/2))",
                                   "total_renormalized_gauge_coupling_fixed": False},
        "unfulfilled": ["derivation of the U(1) internal/gauge representation from the full NSC algebra",
                        "actual scalar-link profile and mass from the throat/operator equations",
                        "five-dimensional determinant phase and boundary anomaly/domain matching",
                        "same-action renormalized gravitational, gauge and vacuum coefficients",
                        "magnetic flux geometry, complete fermion loop and state matching",
                        "NSC self-sourcing, child formation and cosmological predictions"],
        "comparison": "all exact fields and authenticated hashes; imported Einstein and Casimir solutions are not recomputed",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite evidence")
    expected = None
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        if expected["source_hashes"] != hashes():
            raise AssertionError("source authentication failed")
    result = calculate()
    if args.check:
        if json.dumps(expected, sort_keys=True) != json.dumps(result, sort_keys=True):
            raise AssertionError("exact fields differ")
        print("charged compact sector: parity, anomaly and gauge normalization reproduced")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
