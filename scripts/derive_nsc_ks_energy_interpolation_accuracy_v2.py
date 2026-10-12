#!/usr/bin/env python3
"""Ideal energy interpolation bound for a declared nonzero w/U history."""
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

import numpy as np
from flint import arb, ctx, fmpq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_current_history_bounds import (
    radius_bounds, rational_record, value_integral_bounds)
from recursive_horizons.nsc_ks_energy_interpolation_bound import (
    energy_derivative_majorants, weighted_A_up_frobenius, field_remainders_from_basis)
from recursive_horizons.nsc_ks_energy_node_bound import (
    actual_node_remainder_factor, multiply_derivative_bound)
from recursive_horizons.nsc_ks_energy_propagator import chebyshev_energy_nodes
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_finite_matter_error import source_norm_upper, finite_matter_error
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from derive_nsc_ks_coupled_retained_control import interpolation_interval

HISTORY = "results/development/nsc-ks-gate-history-lm-broyden.json"
BACKGROUND = "results/development/nsc-ks-radius-bound.json"
OWNERS = (
    "scripts/derive_nsc_ks_energy_interpolation_accuracy_v2.py",
    "src/recursive_horizons/nsc_ks_current_history_bounds.py",
    "src/recursive_horizons/nsc_ks_energy_interpolation_bound.py",
    "src/recursive_horizons/nsc_ks_energy_node_bound.py",
    "src/recursive_horizons/nsc_ks_finite_matter_error.py",
    "src/recursive_horizons/nsc_ks_retained_upstream_archive.py",
    "src/recursive_horizons/nsc_local_incoming_family.py",
    "src/recursive_horizons/nsc_ks_energy_propagator.py",
    "src/recursive_horizons/nsc_local_incoming_constraints.py",
    "scripts/derive_nsc_ks_coupled_retained_control.py",
    "src/recursive_horizons/evidence_io.py",
    "requirements-validation.txt",
)
BITS = 160


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def ball(value):
    value = Q(value)
    return arb(fmpq(value.numerator, value.denominator))


def compute(history_path, degree):
    if degree not in (32, 48, 64):
        raise ValueError("declared operator degree 32, 48 or 64 required")
    history_path = Path(history_path)
    if not history_path.is_absolute():
        history_path = ROOT / history_path
    if not history_path.resolve().is_relative_to(ROOT):
        raise ValueError("history must lie inside the laboratory")
    history = json.loads(history_path.read_text())
    family = LocalIncomingFamily(np.array(history["history"]["coefficients"], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != history["profile_identity"]:
        raise ValueError("current history identity changed")
    background = json.loads((ROOT / BACKGROUND).read_text())
    archive = RetainedUpstreamArchive(ROOT)
    records = []
    total = [arb(0), arb(0)]
    geometry = None
    with ctx.workprec(BITS):
        for key in archive.family_keys:
            entries = archive.family_entries(key)
            first = next(batch for batch, _ in entries if batch.energy_sign > 0)
            bounds = radius_bounds(family, float(first.rho_up))
            for name in ("rho_upper", "axial_lower", "reference_radius_lower"):
                if bounds[name] != Q(background["bounds"][name]["exact_rational"]):
                    raise ValueError("imported background lemma differs")
            geometry = rational_record(bounds)
            integrals = value_integral_bounds(bounds, float(first.mass), float(first.angular))
            K0, D, K1 = (ball(integrals[name]) for name in ("K0", "D", "K1"))
            interval = interpolation_interval(archive._families[key])
            nodes = chebyshev_energy_nodes(interval, degree)
            factor = actual_node_remainder_factor(nodes, interval, bits=BITS)
            # Only W and W_z are consumed. The zero J inputs do not establish
            # a bound for the actual retarded history derivatives.
            majorants = energy_derivative_majorants(K0, D, K1, 0, 0, len(nodes), bits=BITS)
            ew = multiply_derivative_bound(factor, majorants["dE_n_W_upper"], bits=BITS)
            ewz = multiply_derivative_bound(factor, majorants["dE_n_Wz_upper"], bits=BITS)
            norm2 = arb(0)
            gamma = {1: arb(0), -1: arb(0)}
            max_energy = arb(0)
            multiplicity = None
            for batch, channel in entries:
                normalized = _channel_record("_", {"_": channel})
                if multiplicity is not None and normalized["multiplicity"] != multiplicity:
                    raise ValueError("signed family multiplicity differs")
                multiplicity = normalized["multiplicity"]
                # C_src is block diagonal in the retained energy labels.
                source = batch.source
                if len(source.energies) % 3:
                    raise ValueError("owned three-fibre preparation required")
                for offset in range(0, len(source.energies), 3):
                    block = source.covariance[offset:offset+3, offset:offset+3]
                    outside = np.array(source.covariance[offset:offset+3], copy=True)
                    outside[:, offset:offset+3] = 0
                    if np.count_nonzero(outside):
                        raise ValueError("source has cross-energy coherence outside this bound")
                    gamma[batch.energy_sign] = arb.max(
                        gamma[batch.energy_sign], source_norm_upper(block, bits=BITS))
                if batch.energy_sign > 0:
                    norm = weighted_A_up_frobenius(
                        batch.initial_columns, source.column_weights, bits=BITS)
                    norm2 += norm**2
                    max_energy = arb.max(max_energy, arb(float(np.max(source.energies))))
            norm = norm2.sqrt().upper()
            remainder = field_remainders_from_basis(ew, ewz, norm, max_energy, bits=BITS)
            fnorm = arb(2).sqrt()*K0.exp()*norm
            fznorm = arb(2).sqrt()*K0.exp()*(K1+max_energy)*norm
            signed = {}
            for energy_sign in (1, -1):
                error = finite_matter_error(
                    fnorm, fznorm,
                    restored_upper(remainder["weighted_F_remainder_upper"]),
                    restored_upper(remainder["weighted_F_z_remainder_upper"]),
                    gamma[energy_sign], mass=ball(float(first.mass)),
                    absolute_angular=ball(abs(float(first.angular))),
                    axial_lower=ball(bounds["axial_lower"]),
                    radius_lower=ball(bounds["radius_lower"]),
                    multiplicity=ball(float(multiplicity)), bits=BITS)
                signed[str(energy_sign)] = error
                for i, name in enumerate(("N", "beta")):
                    total[i] += restored_upper(error[name])
            records.append({
                "positive_family": list(key), "interval": list(interval),
                "actual_nodes": nodes.tolist(),
                "integral_bounds": rational_record(integrals),
                "actual_node_factor": factor,
                "W_interpolation_remainder": ew, "Wz_interpolation_remainder": ewz,
                "weighted_A_up_norm_upper": exact_upper(norm),
                "signed_covariance_norm_uppers": {str(k): exact_upper(v) for k, v in gamma.items()},
                "separate_signed_matter_error_uppers": signed,
            })
        return {
            "schema": "NSC-KS-CURRENT-ENERGY-INTERPOLATION-ACCURACY-v2",
            "status": "COMPONENT: ideal interpolation enclosure; other gate errors OPEN",
            "profile_identity": identity, "degree": degree, "precision_bits": BITS,
            "geometry_bounds": geometry, "families": records,
            "covered_positive_families": len(records), "covered_signed_families": 2*len(records),
            "constraint_order": ["N", "beta"],
            "total_N_beta_interpolation_error_upper": [exact_upper(value) for value in total],
            "within_provisional_allocation": [bool(value.upper() <= ball(Q(1, 10**12))) for value in total],
            "scope": {
                "actual_binary_nodes_enclosed": True, "U_profile_included": True,
                "original_source_weights_and_coherences": True,
                "node_evolution_error_bound": None, "interpolation_arithmetic_error_bound": None,
                "source_preparation_error_bound": None, "source_quadrature_error_bound": None,
                "changed_history_UV_bound": None, "history_tangent_interpolation_bound": None,
                "between_node_residual_bound": None, "field_runs": 0, "source_runs": 0,
                "physical_EXISTENCE_certificate": False, "physical_NONEXISTENCE_certificate": False,
            },
            "source_hashes": {name: digest(ROOT/name) for name in OWNERS},
            "input_hashes": {**archive.input_hashes,
                str(history_path.relative_to(ROOT)): digest(history_path), BACKGROUND: digest(ROOT/BACKGROUND)},
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--degree", type=int, default=48)
    parser.add_argument("--history", default=HISTORY)
    args = parser.parse_args()
    output = ROOT / "results/development" / f"nsc-ks-current-energy-interpolation-d{args.degree}.json"
    started = time.process_time()
    result = compute(args.history, args.degree)
    if args.record:
        raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()
        publish_exclusive_file(ROOT, str(output.relative_to(ROOT)), raw)
    elif result != json.loads(output.read_text()):
        raise ValueError("current-history interpolation replay differs")
    print(json.dumps({
        "status": result["status"], "degree": result["degree"],
        "families": len(result["families"]),
        "N_beta_interpolation_upper": [float(restored_upper(v)) for v in result["total_N_beta_interpolation_error_upper"]],
        "within_allocation": result["within_provisional_allocation"],
        "CPU_seconds": time.process_time()-started,
    }, indent=2))
