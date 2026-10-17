#!/usr/bin/env python3
"""Read the saved three-grid mass split. No evolution and no replay.

The completed comparison and its three compact receipts are authenticated
by the same byte bindings as the register and the endpoint figure. ``compare``
is recomputed. Each saved common cut then reports the physical actual and
matched mass residuals, the conformal flux, ``alpha`` times the normalized
trace and balance terms, the explicit time remainder, the clock-scaled
remainder identity, the wall projection gap, and the strong/weak trace gaps.
Algebra residuals are scales and exact-zero-safe ratios. They are not a
pass/fail rule, a dominance test, or a closure classification.

The all-stage envelope is the runner's absolute maximum. No maximum state
is stored. The nine cuts are not every stage. Their distance from the
envelope can guide interpretation, but does not locate the maximizing stage.

    python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_gaps.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json --output results/summaries/nsc-atlas-fixed-endpoint-gap-diagnostic-v1.json

    python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_endpoint_gaps.py -q
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from analyze_nsc_atlas_fixed_endpoint import GAPS, compare
from check_prediction_register import (
    assert_comparison_matches, endpoint_receipts, finite_amplitude_block,
    load_comparison, local_file, read_json,
)
from register_nsc_atlas_fixed_endpoint import endpoint_rows

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-endpoint-gap-decomposition-v1"
ENVELOPE_KEYS = {
    "physical_mass_actual": "temporal.physical_temporal_mass_actual",
    "physical_mass_matched": "temporal.physical_temporal_mass_matched",
    "conformal_wall_flux": "boundary.conformal_wall_flux",
}
CONNECTION_KEYS = {
    "physical_mass_actual": GAPS["mass_trace"],
    "physical_mass_matched": GAPS["mass_matched"],
    "conformal_wall_flux": GAPS["wall_flux"],
}
NOT_INFERRED = (
    "No maximum state is stored. The all-stage envelope is an absolute maximum over stage calls.",
    "The nine common cuts are saved samples. They are not every substage.",
    "A sample whose absolute value equals the envelope is still not a stored maximum state.",
    "The sample-to-envelope ratio describes coverage; exact equality is not a tolerance verdict or identification of the maximizing substage.",
    "Identity residuals are scales and ratios. An exact-zero denominator is reported and is not divided.",
    "No numeric tolerance classifies closure, term dominance, or pass/fail.",
    "Agreement of the mass-rate identities does not establish the physical mass or lapse closure.",
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def real(value, label):
    require(type(value) in (int, float) and math.isfinite(value), label + " must be a finite real")
    return float(value)


def mapping(value, label):
    require(isinstance(value, dict), label + " must be an object")
    return value


def safe_quotient(numerator, denominator):
    """Signed ratio; division by zero, including 0/0, has no numeric value."""
    numerator = float(numerator)
    denominator = float(denominator)
    exact_zero = denominator == 0.0
    if exact_zero:
        ratio = None
    else:
        ratio = numerator / denominator
        require(math.isfinite(ratio), "quotient is not representable")
    return dict(numerator=numerator, denominator=denominator, ratio=ratio,
                exact_zero_denominator=exact_zero)


def comparison_binding(root, path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    resolved = candidate.resolve()
    require(resolved.is_file() and resolved.is_relative_to(root.resolve()),
            "comparison must be an existing repository file")
    relative = resolved.relative_to(root.resolve()).as_posix()
    raw = local_file(root, relative).read_bytes()
    return dict(path=relative, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())


def prepare_output(root, path):
    candidate = Path(path)
    if not candidate.is_absolute():
        candidate = root / candidate
    require(candidate.suffix.lower() == ".json", "output must be a JSON path")
    require(candidate.parent.is_dir(), "output parent directory must already exist")
    resolved = candidate.resolve()
    require(resolved.is_relative_to(root.resolve()), "output must stay inside the repository")
    require(not candidate.exists() and not resolved.exists(), "output already exists")
    return resolved


def authenticate(root, comparison_path):
    """Bind the comparison bytes, its three receipts, and a fresh ``compare``."""
    binding = comparison_binding(root, comparison_path)
    saved = read_json(local_file(root, binding["path"]).read_bytes())
    loaded, found = load_comparison(root, binding)
    require(loaded == saved, "comparison bytes changed while loading")
    rows = endpoint_rows(found)
    records = endpoint_receipts(found, rows)
    recomputed = compare(records)
    assert_comparison_matches(saved, recomputed)
    finite_amplitude_block(saved, binding, rows)
    require(saved["physical_evolutions"] == 0, "gap diagnostic cannot record a new evolution")
    require(len(records) == 3 and [record["level"] for record in records] == [0, 1, 2],
            "three levels 0, 1, 2 required")
    return binding, saved, records, recomputed


def _cut_numbers(sample, alpha, label):
    temporal = mapping(sample.get("temporal_residual"), label + " temporal residual")
    decomposition = mapping(sample.get("mass_gap_decomposition"), label + " mass decomposition")
    background = mapping(sample.get("background_projection"), label + " background projection")
    mass_actual = real(temporal.get("mass_actual"), label + " mass actual")
    mass_matched = real(temporal.get("mass_matched"), label + " mass matched")
    flux = real(temporal.get("flux"), label + " flux")
    clock = real(decomposition.get("clock_rate"), label + " clock rate")
    trace = real(decomposition.get("normalized_trace_term"), label + " normalized trace")
    balance = real(decomposition.get("normalized_balance_term"), label + " normalized balance")
    stored_reconstructed = real(decomposition.get("reconstructed_physical_actual"),
                                label + " stored reconstructed actual")
    stored_identity = real(decomposition.get("identity_defect"), label + " stored identity defect")
    remainder = real(background.get("explicit_time_remainder"), label + " explicit time remainder")
    remainder_identity = real(background.get("remainder_identity_defect"),
                              label + " remainder identity defect")
    wall_gap = real(background.get("wall_time_projection_gap"), label + " wall projection gap")
    momentum = real(decomposition.get("p_trace_minus_weak"), label + " momentum gap")
    normal = real(decomposition.get("h_trace_minus_weak"), label + " normal gap")
    physical_trace = alpha * trace
    physical_balance = alpha * balance
    reconstructed = physical_trace + physical_balance
    identity = mass_actual - reconstructed
    scaled = safe_quotient(remainder_identity, clock)
    balance_minus_remainder = balance - remainder
    if scaled["exact_zero_denominator"]:
        closed_residual = None
    else:
        closed_residual = balance_minus_remainder - scaled["ratio"]
    return dict(
        mass_actual=mass_actual, mass_matched=mass_matched, flux=flux,
        physical_trace_term=physical_trace, physical_balance_term=physical_balance,
        normalized_trace_term=trace, normalized_balance_term=balance,
        explicit_time_remainder=remainder, clock_rate=clock,
        clock_scaled_remainder_identity_defect=scaled,
        wall_time_projection_gap=wall_gap,
        p_trace_minus_weak=momentum, h_trace_minus_weak=normal,
        reconstruction_identity_defect=identity,
        algebra_residuals=dict(
            stored_identity_defect=stored_identity,
            stored_minus_recomputed_identity=stored_identity - identity,
            stored_reconstructed_physical_actual=stored_reconstructed,
            stored_minus_recomputed_physical_actual=stored_reconstructed - reconstructed,
            balance_minus_remainder=balance_minus_remainder,
            balance_minus_remainder_and_scaled_identity=closed_residual,
            absolute_scales=dict(
                reconstruction_identity_defect=abs(identity),
                stored_minus_recomputed_identity=abs(stored_identity - identity),
                stored_minus_recomputed_physical_actual=abs(stored_reconstructed - reconstructed),
                balance_minus_remainder=abs(balance_minus_remainder),
                balance_minus_remainder_and_scaled_identity=(
                    None if closed_residual is None else abs(closed_residual)),
            ),
            ratios=dict(
                physical_trace_over_mass_actual=safe_quotient(physical_trace, mass_actual),
                physical_balance_over_mass_actual=safe_quotient(physical_balance, mass_actual),
                remainder_over_normalized_balance=safe_quotient(remainder, balance),
            ),
        ),
    )


def cut_report(sample, alpha, index):
    label = "sample " + str(index)
    station = real(sample.get("x"), label + " x")
    return dict(x=station, role="common_cut", **_cut_numbers(sample, alpha, label))


def endpoint_numbers(record, alpha):
    endpoint = mapping(record.get("endpoint"), "endpoint")
    station = real(record["execution"]["last_accepted_x"], "endpoint station")
    return station, _cut_numbers(endpoint, alpha, "endpoint")


def signed_difference(left, right):
    skip = {"clock_scaled_remainder_identity_defect", "algebra_residuals"}
    return {key: left[key] - right[key] for key in left if key not in skip and key in right}


def envelope_coverage(samples, maximum_absolute):
    maxima = mapping(maximum_absolute, "maximum absolute")
    coverage = {}
    for name, temporal_key in ENVELOPE_KEYS.items():
        require(temporal_key in maxima, "missing all-stage envelope " + temporal_key)
        require(CONNECTION_KEYS[name] in maxima, "missing connection envelope " + CONNECTION_KEYS[name])
        envelope = real(maxima[temporal_key], temporal_key)
        connection = real(maxima[CONNECTION_KEYS[name]], CONNECTION_KEYS[name])
        require(envelope >= 0 and connection >= 0, name + " envelope must be an absolute maximum")
        field = {"physical_mass_actual": "mass_actual",
                 "physical_mass_matched": "mass_matched",
                 "conformal_wall_flux": "flux"}[name]
        sample_abs = [abs(sample[field]) for sample in samples]
        saved = max(sample_abs) if sample_abs else 0.0
        coverage[name] = dict(
            saved_sample_absolute_maximum=saved,
            all_stage_envelope=envelope,
            connection_envelope=connection,
            connection_minus_temporal_envelope=connection - envelope,
            sample_over_envelope=safe_quotient(saved, envelope),
            saved_sample_absolute_equals_envelope=saved == envelope,
        )
    return coverage


def level_report(record):
    alpha = real(record["input"]["alpha"], "alpha")
    samples = record.get("samples")
    require(isinstance(samples, list) and len(samples) == 9, "nine common samples required")
    cuts = [cut_report(sample, alpha, index) for index, sample in enumerate(samples)]
    station, endpoint = endpoint_numbers(record, alpha)
    matches = [cut for cut in cuts if cut["x"] == station]
    if matches:
        saved = matches[-1]
        difference = signed_difference(endpoint, saved)
        represented = True
        extra = None
    else:
        difference = None
        represented = False
        extra = dict(x=station, role="endpoint", **endpoint)
    rows = cuts if extra is None else [*cuts, extra]
    return dict(
        level=record["level"],
        source_commit=record["source_commit"],
        cap_cells=record["geometry"]["cap_cells"],
        exterior_cells=record["geometry"]["exterior_cells"],
        alpha=alpha,
        endpoint_station=station,
        endpoint_already_a_common_cut=represented,
        endpoint_minus_saved_cut=difference,
        cuts=rows,
        envelope_coverage=envelope_coverage(cuts, record.get("maximum_absolute")),
    )


def gap_report(records, recomputed):
    levels = [level_report(record) for record in records]
    observable = recomputed["physical_observable"]
    return dict(
        schema=SCHEMA,
        status="completed_saved_receipt_gap_decomposition",
        physical_evolutions=0,
        replays=0,
        closure_classification=None,
        dominance_classification=None,
        units=dict(
            convention="Dimensionless evolution variables; physical mass normalization includes alpha. Only the radius comparison is converted to SI.",
            mass_actual="physical: alpha*(temporal - A^2*exp(lapse)*p*h/2)",
            mass_matched="physical: alpha*(temporal - A^2*exp(lapse)*R^2*Pi*Phi/2)",
            flux="conformal wall flux, not multiplied by alpha",
            physical_trace_term="physical: alpha*normalized_trace_term",
            physical_balance_term="physical: alpha*normalized_balance_term",
            normalized_trace_term="normalized strong/weak product gap",
            normalized_balance_term="normalized weak_balance/clock_rate",
            explicit_time_remainder="normalized background time-projection remainder",
            clock_scaled_remainder_identity_defect="normalized remainder_identity_defect/clock_rate",
            wall_time_projection_gap="k_0(Z)-e·a_h, not multiplied by alpha",
            p_trace_minus_weak="strong minus weak momentum trace",
            h_trace_minus_weak="strong minus weak normal trace",
            reconstruction_identity_defect="physical mass_actual - alpha*(trace + balance)",
            envelopes="absolute maxima",
        ),
        endpoint_comparison=dict(
            status=recomputed["status"],
            prediction_status=recomputed["prediction_status"],
            requested_mass_and_flux_contraction=recomputed["requested_mass_and_flux_contraction"],
            endpoint_radius_pair_difference_contracts=recomputed[
                "endpoint_radius_pair_difference_contracts"],
            endpoint_radius_pair_differences=recomputed["endpoint_radius_pair_differences"],
            shared_alpha=recomputed["shared_alpha"],
            sample_cuts=recomputed["sample_cuts"],
            jordan_radius_metres=observable["value_metres"],
            resolution_half_width_metres=observable["resolution_half_width_metres"],
            resolution_range_metres=observable["resolution_range_metres"],
            difference_from_quiet_reference_metres=observable[
                "difference_from_quiet_reference_metres"],
        ),
        levels=levels,
        not_inferred=list(NOT_INFERRED),
    )


def assert_compact(value, path="report"):
    if isinstance(value, dict):
        for key, item in value.items():
            require(isinstance(key, str), "non-text key in " + path)
            assert_compact(item, path + "." + key)
        return
    if isinstance(value, list):
        require(not value or not isinstance(value[0], (int, float)) or len(value) <= 16,
                "raw numeric series in " + path)
        for index, item in enumerate(value):
            assert_compact(item, path + "[" + str(index) + "]")
        return
    require(value is None or type(value) in (str, bool, int, float), "unsupported value in " + path)
    if type(value) is float:
        require(math.isfinite(value), "nonfinite value in " + path)


def build_report(root, comparison_path):
    binding, saved, records, recomputed = authenticate(root, comparison_path)
    report = gap_report(records, recomputed)
    report["comparison"] = binding
    report["inputs"] = list(saved["inputs"])
    report["comparison_provenance"] = dict(
        analysis_source_commit=saved["analysis_source_commit"],
        analysis_source_sha256=saved["analysis_source_sha256"],
    )
    report["diagnostic_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    assert_compact(report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    output = prepare_output(ROOT, args.output)
    report = build_report(ROOT, args.comparison)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(text)
    print(json.dumps(dict(status=report["status"], physical_evolutions=0, replays=0,
                          levels=len(report["levels"]))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
