#!/usr/bin/env python3
"""One-cell local-Fourier residual enclosure for the corrected current history.

This successor replaces only the unusable global derivative-L1 profile tail
in v3.  It reuses the immutable saved DOP853 cell, directly integrates the
retained profile coefficients, and encloses the omitted band on local panels.
It remains a one-cell method pilot, not a whole-field certificate.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import validate_nsc_ks_current_field_pilot_v3 as V3
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_operator import (
    AnalyticRadiusFamily, time_operator_enclosure)
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_difference_residual_polynomial import (
    DifferenceResidualPolynomial, ball_split, difference_operator_remainder_bounds,
    difference_polynomial_bounds, monomial_keys)
from recursive_horizons.nsc_ks_fourier_residual_bound import (
    polynomial_residual_bounds, profile_sup_bounds)
from recursive_horizons.nsc_ks_local_profile_fourier_bound import (
    enclose_profile_fourier_local, serialize_local_profile_fourier)


OUTPUT = ROOT / "results/development/nsc-ks-current-field-pilot-v4.json"
BITS = 90
SETTINGS = {
    "retained_index": 256,
    "taylor_order": 7,
    "exterior_panels": 512,
    "transition_panels": 256,
    "interior_panels": 512,
    "max_derivative": 2,
    "flat_denominator": 256,
    "relative_tolerance_bits": 45,
    "absolute_tolerance_bits": 65,
    "max_depth": 8,
    "bits": BITS,
}
SOURCE_PATHS = (
    "scripts/validate_nsc_ks_current_field_pilot_v4.py",
    "scripts/validate_nsc_ks_current_field_pilot_v3.py",
    "src/recursive_horizons/nsc_ks_local_profile_fourier_bound.py",
    "src/recursive_horizons/nsc_ks_difference_residual_polynomial.py",
    "src/recursive_horizons/nsc_ks_residual_polynomial.py",
    "src/recursive_horizons/nsc_ks_difference_error.py",
    "tests/test_nsc_ks_local_profile_fourier_bound.py",
    "docs/nsc-ks-current-field-pilot-v4.md",
)


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT / path).read_bytes()).hexdigest()


def _arb_fraction(value):
    return arb(value.numerator) / arb(value.denominator)


def _stable(value):
    return {key: item for key, item in value.items() if key != "runtime"}


def compute():
    before = {path: digest(path) for path in SOURCE_PATHS}
    capture, arrays = V3.load_capture()
    previous = json.loads(V3.OUTPUT.read_text())
    family = V3.load_family(capture)
    bounds, _rejected, radius_tail = V3.geometry_bounds(family, capture)
    segment, grid, origin, period, period_q, weights, energies = V3.segment_from_capture(
        capture, arrays)
    keys = monomial_keys(V3.RECIPROCAL_ORDER)
    started_cpu, started_wall = time.process_time(), time.monotonic()
    with ctx.workprec(BITS):
        enclosed = enclose_profile_fourier_local(
            family, keys, origin, period, **SETTINGS)
        profiles = {key: {
            "coefficients": list(item.coefficients),
            "tail": list(item.uniform_tail_bounds),
        } for key, item in enclosed.items()}
        field_x, field_a, field_d = ball_split(
            segment, weights, len(grid), period, bits=BITS)
        polynomials, time_errors = time_operator_enclosure(
            AnalyticRadiusFamily(family), field_x.rho_start, field_x.rho_end,
            capture["angular"], degree=V3.TIME_DEGREE,
            reciprocal_order=V3.RECIPROCAL_ORDER, bits=BITS)
        potential = {key: value for key, value in polynomials.items()
                     if isinstance(key, tuple)}
        if set(potential) != set(keys):
            raise ValueError("time-series monomials differ from the reciprocal keys")
        residual = DifferenceResidualPolynomial(
            field_x, field_a, polynomials["inv_a2"], polynomials["inv_a"],
            polynomials["inv_ar"], potential, capture["mass"],
            capture["angular"], energies)
        polynomial = difference_polynomial_bounds(residual, profiles)
        finite_profiles = {key: {**row, "tail": (arb(0), arb(0))}
                           for key, row in profiles.items()}
        finite = polynomial_residual_bounds(residual, finite_profiles)
        radius_remainder = [_arb_fraction(value) for value in
                            radius_tail["potential_derivative_tail_bounds"][:2]]
        remainder = difference_operator_remainder_bounds(
            field_d, field_a, time_errors,
            profile_sup_bounds(profiles, field_x.length), radius_remainder,
            capture["mass"], capture["angular"], energies)
        total = tuple((left + right).upper()
                      for left, right in zip(polynomial, remainder))
    after = {path: digest(path) for path in SOURCE_PATHS}
    if before != after:
        raise ValueError("field-v4 source changed during evaluation")
    inventory = {
        f"{p}_{q}": {
            "powers": [p, q],
            "retained_index": item.retained_index,
            "tail_upper": [exact_upper(value) for value in item.uniform_tail_bounds],
            "coefficient_flat_strip_error_upper": exact_upper(
                item.coefficient_flat_strip_error),
            "local_cells": item.local_cells,
            "maximum_subdivision_depth": item.maximum_subdivision_depth,
        } for (p, q), item in enclosed.items()
    }
    largest = max(inventory, key=lambda key: float(
        arb(inventory[key]["tail_upper"][0]["mantissa"])
        * arb(2) ** inventory[key]["tail_upper"][0]["exponent"]))
    previous_total = [restored_upper(value) for value in
                      previous["bounds"]["total_continuous_normalized_residual"]]
    improvement = [float(old / new) for old, new in zip(previous_total, total)]
    return {
        "schema": "NSC-KS-CURRENT-FIELD-PILOT-v4",
        "status": (
            "OPEN: one-cell local-Fourier residual enclosure; whole backward "
            "cone, all families and source error remain missing"),
        "profile_identity": capture["profile_identity"],
        "family": capture["family"],
        "cell_index": capture["cell_index"],
        "method": (
            "direct acb Fourier coefficients plus local midpoint-Taylor "
            "enclosure of analytic-profile minus retained Fourier polynomial"),
        "settings": SETTINGS,
        "effective_period_rational": [period_q.numerator, period_q.denominator],
        "profile_enclosure": inventory,
        "profile_coefficient_payload": serialize_local_profile_fourier(enclosed),
        "largest_order0_tail_key_display": largest,
        "bounds": {
            "polynomial_finite_band": [exact_upper(value) for value in finite],
            "polynomial_and_local_profile_tail": [exact_upper(value) for value in polynomial],
            "time_coefficient_and_radius_remainder": [exact_upper(value) for value in remainder],
            "total_continuous_normalized_residual": [exact_upper(value) for value in total],
        },
        "v3_total_continuous_normalized_residual": [
            exact_upper(value) for value in previous_total],
        "v3_to_v4_improvement_factor_display": improvement,
        "coverage": {
            "whole_time_cell": True,
            "whole_spatial_period": True,
            "all_history_cells": False,
            "all_source_families": False,
            "source_columns": len(weights),
        },
        "certificate_use": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "next_required_step": (
            "capture and sum all current-history cells for representative "
            "families, propagate the residual through nsc_ks_difference_error, "
            "then test whole-source scaling"),
        "source_hashes": before,
        "input_hashes": {
            str(V3.CAPTURE.relative_to(ROOT)): digest(V3.CAPTURE),
            capture["payload"]["path"]: capture["payload"]["sha256"],
            capture["history_path"]: digest(capture["history_path"]),
            str(V3.OUTPUT.relative_to(ROOT)): digest(V3.OUTPUT),
        },
        "runtime": {
            "CPU_seconds": time.process_time() - started_cpu,
            "wall_seconds": time.monotonic() - started_wall,
            "python": platform.python_version(),
        },
    }


def display(value):
    print(json.dumps({
        "status": value["status"],
        "bounds": value["bounds"],
        "largest_order0_tail_key_display": value["largest_order0_tail_key_display"],
        "runtime": value["runtime"],
    }, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    modes.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    if args.check and not OUTPUT.exists():
        raise FileNotFoundError("field-v4 record is missing")
    if args.check:
        recorded = json.loads(OUTPUT.read_text())
        if recorded.get("schema") != "NSC-KS-CURRENT-FIELD-PILOT-v4":
            raise ValueError("unexpected field-v4 schema")
        for path, expected in {**recorded["source_hashes"],
                               **recorded["input_hashes"]}.items():
            if digest(path) != expected:
                raise ValueError("field-v4 dependency changed: " + path)
        if recorded["certificate_use"] or recorded["coverage"]["all_history_cells"]:
            raise ValueError("one-cell field pilot cannot be certificate evidence")
        display(recorded)
        return
    result = compute()
    if args.record:
        raw = (json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)), raw)
    elif args.replay:
        recorded = json.loads(OUTPUT.read_text())
        if _stable(result) != _stable(recorded):
            raise ValueError("field-v4 replay differs")
    display(result)


if __name__ == "__main__":
    main()
