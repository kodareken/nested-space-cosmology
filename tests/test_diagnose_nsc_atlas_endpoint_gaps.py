"""Saved-receipt mass split. Synthetic receipts only; no evolution."""
import copy
import hashlib
import json
from pathlib import Path

import pytest

import diagnose_nsc_atlas_endpoint_gaps as tool
from analyze_nsc_atlas_fixed_endpoint import GAPS, compare
from recursive_horizons.nsc_atlas_endpoint_comparison import ENDPOINT_SOURCE_X

ALPHA = 0.00012983127142551443
LEVELS = [
    "results/summaries/synthetic-gap-level0.json",
    "results/summaries/synthetic-gap-level1.json",
    "results/summaries/synthetic-gap-level2.json",
]
COMPARISON = "results/summaries/synthetic-gap-comparison.json"


def dump(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def piece(*, alpha, trace, balance, mass_actual, mass_matched, flux, clock, remainder,
          remainder_identity, wall_gap, momentum, normal, stored_identity=None,
          stored_reconstructed=None):
    physical = alpha * (trace + balance)
    if stored_identity is None:
        stored_identity = mass_actual - physical
    if stored_reconstructed is None:
        stored_reconstructed = physical
    return dict(
        temporal_residual=dict(mass_actual=mass_actual, mass_matched=mass_matched, flux=flux),
        mass_gap_decomposition=dict(
            clock_rate=clock, p_trace_minus_weak=momentum, h_trace_minus_weak=normal,
            normalized_trace_term=trace, normalized_balance_term=balance,
            reconstructed_physical_actual=stored_reconstructed, identity_defect=stored_identity),
        background_projection=dict(
            explicit_time_remainder=remainder, remainder_identity_defect=remainder_identity,
            wall_time_projection_gap=wall_gap),
    )


def consistent_piece(alpha, sign, clock, trace, balance_shift, flux, momentum, normal, wall_gap,
                     identity_offset=0.0):
    remainder_identity = sign * 4e-12
    scaled = remainder_identity / clock
    remainder = sign * -0.01
    balance = remainder + scaled + balance_shift
    mass_actual = alpha * trace + alpha * balance + identity_offset
    return piece(
        alpha=alpha, trace=trace, balance=balance, mass_actual=mass_actual,
        mass_matched=sign * -0.002, flux=flux, clock=clock, remainder=remainder,
        remainder_identity=remainder_identity, wall_gap=wall_gap, momentum=momentum,
        normal=normal, stored_identity=identity_offset,
        stored_reconstructed=alpha * (trace + balance))


def receipt(level, samples, endpoint, envelopes, alpha=ALPHA, radius=.54):
    maximum = {key: (1.0, 0.5, 0.25)[level] for key in GAPS.values()}
    maximum.update(envelopes)
    return dict(
        schema="nsc-atlas-fixed-endpoint-v1", level=level, status="completed_numerical_endpoint",
        source_commit=str(level) * 40,
        radial_accuracy=dict(rtol=2.5e-14, atol=1e-15, adjoint_order=32,
                             gradient_absolute_tolerance=2e-8),
        execution=dict(last_accepted_x=ENDPOINT_SOURCE_X, raw_history_saved=False,
                       final_state_saved=False),
        geometry=dict(x_left=.4, x_right=ENDPOINT_SOURCE_X, reference_radius_endpoint=.55,
                      cap_cells=[2, 4, 8][level], exterior_cells=[32, 64, 128][level]),
        input=dict(sha256="a" * 64, alpha=alpha, source_rows=257),
        preparation=dict(input_sha256="b" * 64, original_input_sha256="c" * 64),
        endpoint=dict(einstein_radius=radius, jordan_radius=radius, c=0.,
                      frame_formula_needs_check=False, **endpoint),
        samples=samples,
        maximum_absolute=maximum,
        transport_comparison={key: (1.0, 0.5, 0.25)[level] for key in (
            "mass_absolute_gap", "scalar_absolute_gap",
            "max_abs_running_mass_gap", "max_abs_running_scalar_gap")},
    )


def signed_samples(alpha=ALPHA):
    samples = []
    for index in range(9):
        sign = -1.0 if index % 2 else 1.0
        clock = 2.0 + 0.25 * index
        trace = sign * (0.02 + 0.001 * index)
        flux = sign * (0.2 + 0.01 * index)
        samples.append(dict(
            x=.4 + (ENDPOINT_SOURCE_X - .4) * index / 8,
            **consistent_piece(alpha, sign, clock, trace, 0.0, flux,
                               sign * 1e-5 * (index + 1), -sign * 0.014, sign * (0.5 + index),
                               identity_offset=(1e-6 if index == 4 else 0.0))))
    return samples


def write_tree(root, rows):
    for relative, row in zip(LEVELS, rows):
        dump(root / relative, row)
    saved = compare(copy.deepcopy(rows))
    saved["inputs"] = []
    for relative in LEVELS:
        raw = (root / relative).read_bytes()
        saved["inputs"].append(dict(path=relative, bytes=len(raw),
                                    sha256=hashlib.sha256(raw).hexdigest()))
    saved["analysis_source_commit"] = "ab" * 20
    saved["analysis_source_sha256"] = "cd" * 32
    dump(root / COMPARISON, saved)
    return saved


def run(root, rows):
    write_tree(root, rows)
    output = root / "results/summaries/synthetic-gap-decomposition.json"
    tool.ROOT = root
    try:
        assert tool.main(["--comparison", COMPARISON, "--output", str(output)]) == 0
    finally:
        tool.ROOT = Path(__file__).resolve().parents[1]
    return json.loads(output.read_text(encoding="utf-8"))


def manufactured_rows():
    rows = []
    for level in range(3):
        samples = signed_samples()
        endpoint = {key: copy.deepcopy(samples[-1][key]) for key in (
            "temporal_residual", "mass_gap_decomposition", "background_projection")}
        actual_envelope = max(abs(sample["temporal_residual"]["mass_actual"]) for sample in samples)
        flux_envelope = max(abs(sample["temporal_residual"]["flux"]) for sample in samples) + 1.0
        rows.append(receipt(
            level, samples, endpoint,
            {tool.ENVELOPE_KEYS["physical_mass_actual"]: actual_envelope,
             tool.ENVELOPE_KEYS["physical_mass_matched"]: 0.002,
             tool.ENVELOPE_KEYS["conformal_wall_flux"]: flux_envelope},
            radius=(.54, .545, .546)[level]))
    return rows


def test_manufactured_signs_clock_alpha_and_endpoint_identity(tmp_path):
    report = run(tmp_path, manufactured_rows())
    assert report["schema"] == tool.SCHEMA
    assert report["physical_evolutions"] == 0 and report["replays"] == 0
    assert report["closure_classification"] is None
    assert report["dominance_classification"] is None
    assert report["endpoint_comparison"]["requested_mass_and_flux_contraction"] is True
    assert report["endpoint_comparison"]["jordan_radius_metres"] == pytest.approx(5460.)
    assert report["endpoint_comparison"]["resolution_half_width_metres"] == pytest.approx(60.)
    assert report["comparison"]["sha256"] == hashlib.sha256(
        (tmp_path / COMPARISON).read_bytes()).hexdigest()
    assert [row["sha256"] for row in report["inputs"]] == [
        hashlib.sha256((tmp_path / relative).read_bytes()).hexdigest() for relative in LEVELS]
    assert report["diagnostic_sha256"] == hashlib.sha256(Path(tool.__file__).read_bytes()).hexdigest()
    assert "coefficient" not in json.dumps(report)
    level = report["levels"][1]
    assert level["cap_cells"] == 4 and level["exterior_cells"] == 64
    assert len(level["cuts"]) == 9
    assert level["endpoint_already_a_common_cut"] is True
    assert level["cuts"][-1]["role"] == "common_cut"
    assert level["endpoint_minus_saved_cut"]["mass_actual"] == 0.0
    cut = level["cuts"][5]
    sign = -1.0
    clock = 2.0 + 0.25 * 5
    trace = sign * (0.02 + 0.001 * 5)
    remainder_identity = sign * 4e-12
    scaled = remainder_identity / clock
    remainder = sign * -0.01
    balance = remainder + scaled
    assert cut["mass_actual"] == ALPHA * trace + ALPHA * balance
    assert cut["mass_actual"] < 0 and cut["flux"] < 0 and cut["normalized_trace_term"] < 0
    assert cut["mass_matched"] == sign * -0.002
    assert cut["flux"] == sign * (0.2 + 0.01 * 5)
    assert cut["physical_trace_term"] == ALPHA * trace
    assert cut["physical_balance_term"] == ALPHA * balance
    assert cut["normalized_trace_term"] == trace
    assert cut["explicit_time_remainder"] == remainder
    assert cut["clock_scaled_remainder_identity_defect"]["ratio"] == scaled
    assert cut["clock_scaled_remainder_identity_defect"]["exact_zero_denominator"] is False
    assert cut["wall_time_projection_gap"] == sign * (0.5 + 5)
    assert cut["p_trace_minus_weak"] == sign * 1e-5 * 6
    assert cut["h_trace_minus_weak"] == -sign * 0.014
    assert cut["reconstruction_identity_defect"] == 0.0
    offset = level["cuts"][4]
    assert offset["reconstruction_identity_defect"] == pytest.approx(1e-6, abs=1e-18)
    assert offset["algebra_residuals"]["stored_minus_recomputed_identity"] == pytest.approx(0.0, abs=1e-18)
    assert offset["algebra_residuals"]["balance_minus_remainder_and_scaled_identity"] == pytest.approx(0.0, abs=1e-18)
    assert offset["algebra_residuals"]["absolute_scales"]["reconstruction_identity_defect"] == pytest.approx(1e-6, abs=1e-18)
    assert cut["algebra_residuals"]["balance_minus_remainder_and_scaled_identity"] == pytest.approx(0.0, abs=1e-18)
    assert report["not_inferred"] == list(tool.NOT_INFERRED)
    coverage = level["envelope_coverage"]["physical_mass_actual"]
    assert coverage["saved_sample_absolute_equals_envelope"] is True
    flux_coverage = level["envelope_coverage"]["conformal_wall_flux"]
    assert flux_coverage["saved_sample_absolute_equals_envelope"] is False
    assert flux_coverage["sample_over_envelope"]["exact_zero_denominator"] is False
    assert flux_coverage["sample_over_envelope"]["ratio"] < 1


def test_zero_clock_zero_alpha_and_zero_envelope_stay_unclassified(tmp_path):
    rows = []
    for level in range(3):
        samples = []
        for index in range(9):
            clock = 0.0
            remainder_identity = 5.0 if index == 3 else 0.0
            flux = 0.25 if index == 2 and level == 0 else 0.0
            samples.append(dict(
                x=.4 + (ENDPOINT_SOURCE_X - .4) * index / 8,
                **piece(alpha=0.0, trace=0.0, balance=0.0, mass_actual=0.0, mass_matched=0.0,
                        flux=flux, clock=clock, remainder=0.0,
                        remainder_identity=remainder_identity, wall_gap=0.0, momentum=0.0,
                        normal=0.0)))
        endpoint = {key: copy.deepcopy(samples[-1][key]) for key in (
            "temporal_residual", "mass_gap_decomposition", "background_projection")}
        rows.append(receipt(
            level, samples, endpoint,
            {tool.ENVELOPE_KEYS["physical_mass_actual"]: 0.0,
             tool.ENVELOPE_KEYS["physical_mass_matched"]: 0.0,
             tool.ENVELOPE_KEYS["conformal_wall_flux"]: 0.0 if level else 0.0},
            alpha=0.0, radius=(.54, .545, .546)[level]))
    report = run(tmp_path, rows)
    cut = report["levels"][0]["cuts"][3]
    assert cut["physical_trace_term"] == 0.0
    assert cut["physical_balance_term"] == 0.0
    assert cut["reconstruction_identity_defect"] == 0.0
    scaled = cut["clock_scaled_remainder_identity_defect"]
    assert scaled["exact_zero_denominator"] is True
    assert scaled["numerator"] == 5.0 and scaled["ratio"] is None
    assert cut["algebra_residuals"]["balance_minus_remainder_and_scaled_identity"] is None
    zero = report["levels"][0]["cuts"][0]["clock_scaled_remainder_identity_defect"]
    assert zero["ratio"] is None and zero["exact_zero_denominator"] is True
    actual = report["levels"][0]["envelope_coverage"]["physical_mass_actual"]
    assert actual["sample_over_envelope"]["ratio"] is None
    assert actual["saved_sample_absolute_equals_envelope"] is True
    flux = report["levels"][0]["envelope_coverage"]["conformal_wall_flux"]
    assert flux["saved_sample_absolute_maximum"] == 0.25
    assert flux["all_stage_envelope"] == 0.0
    assert flux["sample_over_envelope"]["ratio"] is None
    assert flux["sample_over_envelope"]["exact_zero_denominator"] is True
    assert flux["saved_sample_absolute_equals_envelope"] is False
    assert report["closure_classification"] is None
    assert report["not_inferred"] == list(tool.NOT_INFERRED)


def test_endpoint_difference_stays_on_the_ninth_cut(tmp_path):
    rows = manufactured_rows()
    rows[0]["endpoint"]["temporal_residual"]["mass_actual"] += 0.25
    report = run(tmp_path, rows)
    level = report["levels"][0]
    assert len(level["cuts"]) == 9
    assert level["endpoint_minus_saved_cut"]["mass_actual"] == pytest.approx(0.25)
    assert level["endpoint_minus_saved_cut"]["explicit_time_remainder"] == 0.0
    assert level["cuts"][-1]["x"] == ENDPOINT_SOURCE_X


def test_endpoint_off_the_common_cuts_is_one_added_row(tmp_path):
    rows = manufactured_rows()
    for row in rows:
        for index, sample in enumerate(row["samples"]):
            sample["x"] = 0.41 + 0.001 * index
    report = run(tmp_path, rows)
    level = report["levels"][2]
    assert level["endpoint_already_a_common_cut"] is False
    assert len(level["cuts"]) == 10
    assert level["cuts"][-1]["role"] == "endpoint"
    assert level["cuts"][-1]["x"] == ENDPOINT_SOURCE_X
    assert level["endpoint_minus_saved_cut"] is None


def test_tampered_receipt_and_malformed_sample_write_nothing(tmp_path):
    rows = manufactured_rows()
    write_tree(tmp_path, rows)
    output = tmp_path / "results/summaries/synthetic-gap-decomposition.json"
    before = {relative: (tmp_path / relative).read_bytes() for relative in (*LEVELS, COMPARISON)}
    raw = bytearray((tmp_path / LEVELS[0]).read_bytes())
    raw[-2] ^= 1
    (tmp_path / LEVELS[0]).write_bytes(raw)
    tool.ROOT = tmp_path
    try:
        with pytest.raises(ValueError, match="bytes changed"):
            tool.main(["--comparison", COMPARISON, "--output", str(output)])
        assert not output.exists()
        (tmp_path / LEVELS[0]).write_bytes(before[LEVELS[0]])
        broken = json.loads(before[LEVELS[1]])
        del broken["samples"][2]["background_projection"]
        dump(tmp_path / LEVELS[1], broken)
        # Rebuild only the comparison binding so authentication reaches the sample.
        saved = json.loads(before[COMPARISON])
        fresh = (tmp_path / LEVELS[1]).read_bytes()
        saved["inputs"][1]["bytes"] = len(fresh)
        saved["inputs"][1]["sha256"] = hashlib.sha256(fresh).hexdigest()
        dump(tmp_path / COMPARISON, saved)
        with pytest.raises(ValueError, match="background projection"):
            tool.main(["--comparison", COMPARISON, "--output", str(output)])
        assert not output.exists()
        assert (tmp_path / LEVELS[0]).read_bytes() == before[LEVELS[0]]
        assert (tmp_path / LEVELS[2]).read_bytes() == before[LEVELS[2]]
    finally:
        tool.ROOT = Path(__file__).resolve().parents[1]


def test_existing_output_and_refusal_do_not_mutate_inputs(tmp_path):
    rows = manufactured_rows()
    write_tree(tmp_path, rows)
    output = tmp_path / "results/summaries/synthetic-gap-decomposition.json"
    output.write_text("keep\n", encoding="utf-8", newline="\n")
    before = {relative: (tmp_path / relative).read_bytes() for relative in (*LEVELS, COMPARISON)}
    tool.ROOT = tmp_path
    try:
        with pytest.raises(ValueError, match="already exists"):
            tool.main(["--comparison", COMPARISON, "--output", str(output)])
        assert output.read_text(encoding="utf-8") == "keep\n"
        outside = tmp_path.parent / "nsc-gap-outside-not-in-repo.json"
        with pytest.raises(ValueError, match="inside the repository"):
            tool.main(["--comparison", COMPARISON, "--output", str(outside)])
        assert not outside.exists()
        with pytest.raises(ValueError, match="JSON"):
            tool.main(["--comparison", COMPARISON,
                       "--output", str(tmp_path / "results/summaries/gap.txt")])
        with pytest.raises(ValueError):
            tool.main(["--comparison", "results/summaries/missing-comparison.json",
                       "--output", str(tmp_path / "results/summaries/other-gap.json")])
    finally:
        tool.ROOT = Path(__file__).resolve().parents[1]
    assert {relative: (tmp_path / relative).read_bytes() for relative in (*LEVELS, COMPARISON)} == before
    assert not (tmp_path / "results/summaries/other-gap.json").exists()
