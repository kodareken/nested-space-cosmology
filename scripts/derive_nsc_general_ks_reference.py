#!/usr/bin/env python3
"""Verify the supplied-history general-KS fourth-order reference owner."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_general_ks_reference import (  # noqa: E402
    FIELDS,
    GeneralKSFourthOrderReferenceHistory,
)
from recursive_horizons.nsc_landau_cauchy_isometry import KSCauchyHistory  # noqa: E402
from recursive_horizons.nsc_mode_resolved_cauchy_state import (  # noqa: E402
    ModeResolvedCauchyState,
)


MANIFEST = ROOT / "results/development/nsc-mode-resolved-cauchy-state.json"
OUTPUT = ROOT / "results/development/nsc-general-ks-reference.json"
SOURCE_PATHS = (
    "src/recursive_horizons/nsc_general_ks_reference.py",
    "scripts/derive_nsc_general_ks_reference.py",
    "docs/nsc-general-ks-reference.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def diagnostic_history(points: int, seed_a: float) -> KSCauchyHistory:
    """Manufactured local stencil control; it is not a physical duration."""
    x = np.linspace(-0.25, 0.25, points)
    return KSCauchyHistory(
        time=x,
        a_parallel=seed_a*np.exp(0.07*x-0.03*x*x+0.01*x**3),
        radius=np.exp(-0.04*x+0.025*x*x-0.008*x**3),
        lapse=np.exp(0.03*x+0.012*x*x),
        shift=0.02*x-0.005*x*x,
        label="manufactured local derivative stencil; not a metric solution",
    )


def calculate():
    manifest = json.loads(MANIFEST.read_text())
    payload = ROOT / manifest["payload"]["path"]
    if sha256(payload) != manifest["payload"]["sha256"]:
        raise RuntimeError("mode-state payload hash mismatch")
    state = ModeResolvedCauchyState.load(payload, manifest["channels"])
    seed_a = float(manifest["cauchy_surfaces"]["seed"]["a_parallel"])
    seed_r = float(manifest["cauchy_surfaces"]["seed"]["r"])
    owner = GeneralKSFourthOrderReferenceHistory(
        state.arrays, state.channels, seed_a_parallel=seed_a, seed_radius=seed_r,
    )
    coarse = owner.evaluate(diagnostic_history(65, seed_a))
    fine = owner.evaluate(diagnostic_history(129, seed_a))

    # The grids share every second node.  Ignore only the one-sided boundary
    # stencils when measuring refinement of the four nodal force channels.
    shared = slice(4, -4)
    refinement = {}
    for name in FIELDS:
        left = coarse.forces[name][shared]
        right = fine.forces[name][::2][shared]
        scale = max(1.0, float(np.max(abs(right))))
        refinement[name] = float(np.max(abs(left-right))/scale)
    maximum_refinement = max(refinement.values())
    center = 64
    tolerances = {
        "seed_hamiltonian_axis": 3e-15,
        "bloch_recursion": 3e-12,
        "order_normalization": 3e-12,
        "projector_trace": 3e-15,
        "parity_completed_beta": 3e-15,
        "force_refinement_relative": 3e-8,
    }
    residuals = fine.residuals
    checks = {
        "seed_hamiltonian_axes_match_serialized_state": (
            residuals["seed_hamiltonian_axis_residual"]
            <= tolerances["seed_hamiltonian_axis"]
        ),
        "bloch_orders_zero_through_four_close": (
            residuals["maximum_bloch_recursion_residual"]
            <= tolerances["bloch_recursion"]
        ),
        "order_by_order_normalization_closes": (
            residuals["maximum_order_normalization_residual"]
            <= tolerances["order_normalization"]
        ),
        "reference_projector_trace_is_one": (
            residuals["maximum_projector_trace_residual"]
            <= tolerances["projector_trace"]
        ),
        "parity_completed_momentum_subtraction_is_zero": (
            residuals["maximum_parity_completed_beta_projection"]
            <= tolerances["parity_completed_beta"]
        ),
        "four_force_channels_refine": (
            maximum_refinement <= tolerances["force_refinement_relative"]
        ),
    }
    if not all(np.isfinite(fine.forces[name]).all() for name in FIELDS):
        raise RuntimeError("nonfinite reference action force")

    return {
        "schema": "NSC-GENERAL-KS-FOURTH-ORDER-REFERENCE-HISTORY-v1",
        "owner": "GeneralKSFourthOrderReferenceHistory",
        "status": (
            "PASS: supplied-history order-0..4 reference and all four ADM "
            "variation channels are executable; extended composition remains OPEN"
        ),
        "domain": {
            "in": (
                "the 33 authenticated ModeResolvedCauchyState channels on a supplied "
                "smooth homogeneous KS history, block-diagonal in axial frequency"
            ),
            "out": [
                "selection of a physical history or duration",
                "frequency mixing and the transmitting tilted interface",
                "local induced action forces",
                "a finite physical stress, null signs, updated constraints, or metric evolution",
            ],
            "hamiltonian_axes": "h=(-m_j,lambda_n/r,k/a_parallel)",
            "normal_derivative": "D=N^-1 d/dtau on the supplied history nodes",
            "orders": [0, 1, 2, 3, 4],
        },
        "construction": {
            "state_manifest": str(MANIFEST.relative_to(ROOT)),
            "state_payload": str(payload.relative_to(ROOT)),
            "state_channels": len(state.channels),
            "state_blocks": len(state.arrays["hx_seed"]),
            "action_measure": "copy_count*degeneracy*dk/pi after explicit +/-k parity average",
            "reference_force_sign": "+Tr(C_ref delta H), composing with -Tr(C delta H)",
            "beta_policy": (
                "pair k and -k before summation; never use the one-sided "
                "positive-frequency payload as standalone vacuum momentum"
            ),
        },
        "diagnostic": {
            "kind": "manufactured local differentiation stencil only",
            "physical_duration_or_history_selected": False,
            "coarse_nodes": 65,
            "fine_nodes": 129,
            "center_reference_action_forces": {
                name: float(fine.forces[name][center]) for name in FIELDS
            },
            "force_refinement_relative": refinement,
            "maximum_force_refinement_relative": maximum_refinement,
        },
        "residuals": {name: float(value) for name, value in residuals.items()},
        "tolerances": tolerances,
        "checks": checks,
        "gate": {
            "owner_residuals_pass": bool(all(checks.values())),
            "all_four_ADM_variations_exposed": True,
            "parity_completed_momentum_subtraction": True,
            "physical_history_selected": False,
            "extended_same_action_gate_ready": False,
            "remaining_owners": [
                "GeneralKSLocalInducedHistory",
                "transmitting tilted-Landau interface with frequency-mixing Cauchy map",
            ],
            "coupled_evolution_reopened": False,
        },
        "locked_inputs": {
            "A": manifest["locked_inputs"]["A"],
            "magnetic_flux": manifest["locked_inputs"]["magnetic_flux"],
            "Omega": manifest["locked_inputs"]["Omega"],
            "zeta": manifest["locked_inputs"]["zeta"],
            "V_full": manifest["locked_inputs"]["V_full"],
            "changed": False,
        },
        "scope": {
            "new_physical_term_or_counterflow": False,
            "finite_stress_or_nulls_fabricated": False,
            "metric_timestep_started": False,
            "old_generator_rerun": False,
        },
        "input_hashes": {
            str(MANIFEST.relative_to(ROOT)): sha256(MANIFEST),
            str(payload.relative_to(ROOT)): sha256(payload),
        },
        "source_hashes": {
            relative: sha256(ROOT / relative) for relative in SOURCE_PATHS
        },
        "comparison": {
            "float_atol": 3e-13,
            "float_rtol": 3e-13,
            "exact": "strings, booleans, integers, keys, input and source hashes",
        },
    }


def compare(expected, actual, path="root"):
    if isinstance(expected, dict):
        if set(expected) != set(actual):
            raise AssertionError(f"{path}: keys differ")
        for key in expected:
            compare(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if len(expected) != len(actual):
            raise AssertionError(f"{path}: list lengths differ")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}[{index}]")
    elif isinstance(expected, float):
        if not np.isclose(expected, actual, atol=3e-13, rtol=3e-13):
            raise AssertionError(f"{path}: {expected!r} != {actual!r}")
    elif expected != actual:
        raise AssertionError(f"{path}: {expected!r} != {actual!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("general-KS fourth-order reference record reproduced")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
