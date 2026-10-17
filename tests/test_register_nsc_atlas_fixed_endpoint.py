"""NSC-P001 keeps its perturbative record when a finite-amplitude extension is added."""
import copy
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

import pytest

import check_prediction_register as checker
import register_nsc_atlas_fixed_endpoint as registrar
from analyze_nsc_atlas_fixed_endpoint import GAPS, compare
from recursive_horizons.nsc_atlas_endpoint_comparison import ENDPOINT_SOURCE_X

ROOT = Path(__file__).resolve().parents[1]
RELATIVE_LEVELS = [
    "results/summaries/synthetic-fixed-endpoint-level0.json",
    "results/summaries/synthetic-fixed-endpoint-level1.json",
    "results/summaries/synthetic-fixed-endpoint-level2.json",
]
COMPARISON = "results/summaries/synthetic-fixed-endpoint-comparison.json"


def dump(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def receipts(radii=(.54, .545, .546), gaps=(1., .5, .25)):
    rows = []
    for level, (radius, gap) in enumerate(zip(radii, gaps)):
        rows.append(dict(
            schema="nsc-atlas-fixed-endpoint-v1", level=level, status="completed_numerical_endpoint",
            source_commit=str(level) * 40,
            radial_accuracy=dict(rtol=2.5e-14, atol=1e-15, adjoint_order=32,
                                 gradient_absolute_tolerance=2e-8),
            execution=dict(last_accepted_x=ENDPOINT_SOURCE_X, raw_history_saved=False,
                           final_state_saved=False),
            geometry=dict(x_left=.4, x_right=ENDPOINT_SOURCE_X, reference_radius_endpoint=.55,
                          cap_cells=[2, 4, 8][level], exterior_cells=[32, 64, 128][level]),
            input=dict(sha256="a" * 64, alpha=.0001, source_rows=257),
            preparation=dict(input_sha256="b" * 64, original_input_sha256="c" * 64),
            endpoint=dict(einstein_radius=radius, jordan_radius=radius, c=0.,
                          frame_formula_needs_check=False),
            samples=[dict(x=.4 + (ENDPOINT_SOURCE_X - .4) * cut / 8) for cut in range(9)],
            maximum_absolute={key: gap for key in GAPS.values()},
            transport_comparison={key: gap for key in (
                "mass_absolute_gap", "scalar_absolute_gap",
                "max_abs_running_mass_gap", "max_abs_running_scalar_gap")},
        ))
    return rows


def historical_root(tmp_path, source=ROOT):
    """Copy one repository and restore the explicit pre-extension P001 view.

    Only the optional finite-amplitude object is removed, and its SI flag is
    reset. A source that already carries a completed extension still yields
    the historical register.
    """
    source = Path(source)
    for name in (checker.REGISTER, checker.DOCUMENT, *checker.EVIDENCE.values()):
        destination = tmp_path / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source / name, destination)
    path = tmp_path / checker.REGISTER
    document = json.loads(path.read_text(encoding="utf-8"))
    saved = document["entries"][0]
    changed = "finite_amplitude" in saved
    saved.pop("finite_amplitude", None)
    if saved["observational_test"].get("SI_scale_selected") is not False:
        saved["observational_test"]["SI_scale_selected"] = False
        changed = True
    if changed:
        path.write_text(json.dumps(document, indent=2, allow_nan=False) + "\n",
                        encoding="utf-8", newline="\n")
    return tmp_path


def write_comparison(root, rows, input_order=(0, 1, 2)):
    for relative, row in zip(RELATIVE_LEVELS, rows):
        dump(root / relative, row)
    comparison = compare(copy.deepcopy(rows))
    comparison["inputs"] = []
    for index in input_order:
        raw = (root / RELATIVE_LEVELS[index]).read_bytes()
        comparison["inputs"].append(dict(
            path=RELATIVE_LEVELS[index], bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    comparison["analysis_source_commit"] = "ab" * 20
    comparison["analysis_source_sha256"] = "cd" * 32
    dump(root / COMPARISON, comparison)
    return comparison


def register_bytes(root):
    return (root / checker.REGISTER).read_bytes()


def entry(root):
    return json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))["entries"][0]


def test_current_register_checker_accepts_either_valid_state():
    saved = json.loads((ROOT / checker.REGISTER).read_text(encoding="utf-8"))["entries"][0]
    result = checker.check(ROOT)
    assert result["predictions"] == 1
    assert result["evidence_files"] == 7
    assert result["coefficient_projections"] == 4
    assert result["exact_identity_projections"] == 2
    assert result["scientific_replays"] == 0
    assert saved["observational_test"]["status"] == "pending"
    assert saved["observational_test"]["calibration"] == "none"
    assert saved["observational_test"]["physical_source_selected"] is False
    assert saved["observational_test"]["SI_scale_selected"] is ("finite_amplitude" in saved)
    if "finite_amplitude" in saved:
        limits = saved["finite_amplitude"]["limits"]
        assert limits["calibration"] == "none"
        assert limits["physical_source_selected"] is False
        assert limits["physical_emitter_selected"] is False
        assert limits["empirical_confirmation"] is False
        assert result["finite_amplitude_endpoints"] == 3
    else:
        assert "finite_amplitude_endpoints" not in result


def test_summary_coverage_cannot_be_promoted(tmp_path):
    saved=write_comparison(tmp_path, receipts())
    recomputed=compare(receipts())
    checker.assert_comparison_matches(saved,recomputed)
    saved['summary_coverage'][0]['original_receipt_status']='stopped'
    with pytest.raises(ValueError,match='summary_coverage'):
        checker.assert_comparison_matches(saved,recomputed)


def test_historical_view_drops_an_installed_extension(tmp_path):
    source = historical_root(tmp_path / "source")
    write_comparison(source, receipts())
    registrar.register(source, COMPARISON)
    installed = entry(source)
    assert installed["finite_amplitude"]["limits"]["calibration"] == "none"
    assert installed["finite_amplitude"]["limits"]["physical_source_selected"] is False
    assert installed["finite_amplitude"]["limits"]["physical_emitter_selected"] is False
    assert installed["observational_test"]["SI_scale_selected"] is True
    view = historical_root(tmp_path / "view", source)
    saved = entry(view)
    assert "finite_amplitude" not in saved
    assert saved["observational_test"]["SI_scale_selected"] is False
    assert saved["observational_test"]["calibration"] == "none"
    assert saved["observational_test"]["physical_source_selected"] is False
    assert saved["observational_test"]["status"] == "pending"
    for key in ("coefficients", "identities", "uncertainty", "open_limits", "evidence", "status",
                "rejection", "parameters"):
        assert saved[key] == installed[key]
    assert checker.check(view) == dict(
        predictions=1, evidence_files=7, coefficient_projections=4,
        exact_identity_projections=2, scientific_replays=0)
    assert checker.check(source)["finite_amplitude_endpoints"] == 3


def test_synthetic_extension_preserves_historical_coefficients(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("registration must not run a generator")
    monkeypatch.setattr(subprocess, "run", forbidden)
    root = historical_root(tmp_path)
    before = entry(root)
    write_comparison(root, receipts(), input_order=(2, 0, 1))
    result = registrar.register(root, COMPARISON)
    saved = entry(root)
    assert result["prediction_status"] == checker.ENVELOPE_STATUS
    assert result["physical_evolutions"] == 0
    assert checker.check(root)["finite_amplitude_endpoints"] == 3
    for key in ("coefficients", "identities", "uncertainty", "open_limits", "evidence", "status",
                "registration", "parameters", "domain", "expansion", "prediction_point", "rejection"):
        assert saved[key] == before[key]
    assert saved["observational_test"]["status"] == "pending"
    assert saved["observational_test"]["calibration"] == "none"
    assert saved["observational_test"]["physical_source_selected"] is False
    assert saved["observational_test"]["SI_scale_selected"] is True
    block = saved["finite_amplitude"]
    assert block["schema"] == "nsc-p001-finite-amplitude-v1"
    assert block["prediction_status"] == checker.ENVELOPE_STATUS
    assert block["event"]["distinct_from_perturbative_rho"] is True
    assert block["event"]["kind"] == "fixed_incoming_label"
    assert block["limits"]["historical_coefficients_predict_this_event"] is False
    assert block["limits"]["physical_emitter_selected"] is False
    assert block["limits"]["empirical_confirmation"] is False
    assert block["limits"]["SI_scale_role"] == "declared_fiducial_input"
    assert block["limits"]["measured_astrophysical_scale"] is False
    assert block["physical_observable"]["uncertainty"]["confidence_level"] is None
    assert block["physical_observable"]["uncertainty"]["rigorous_error_enclosure"] is False
    assert [row["level"] for row in block["evidence"]["endpoints"]] == [0, 1, 2]
    assert saved["coefficients"]["rho"]["value"] == pytest.approx(-5.2868390416906355)
    again = registrar.register(root, COMPARISON)
    assert again["prediction_status"] == result["prediction_status"]
    assert entry(root)["coefficients"] == before["coefficients"]
    assert len(json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))["entries"]) == 1


def test_unconverged_comparison_keeps_its_own_status(tmp_path):
    root = historical_root(tmp_path)
    rows = receipts()
    rows[2]["maximum_absolute"][GAPS["wall_flux"]] = .7
    write_comparison(root, rows)
    result = registrar.register(root, COMPARISON)
    block = entry(root)["finite_amplitude"]
    assert result["prediction_status"] == checker.UNCONVERGED_STATUS
    assert block["prediction_status"] == checker.UNCONVERGED_STATUS
    assert block["requested_mass_and_flux_contraction"] is False
    assert block["comparison_status"] == "completed_saved_receipt_comparison"
    assert checker.check(root)["finite_amplitude_endpoints"] == 3
    document = json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))
    document["entries"][0]["finite_amplitude"]["prediction_status"] = checker.ENVELOPE_STATUS
    dump(root / checker.REGISTER, document)
    with pytest.raises(ValueError, match="finite-amplitude status"):
        checker.check(root)


@pytest.mark.parametrize("mutation", [
    lambda rows: rows[2].update(status="stopped"),
    lambda rows: rows[2]["execution"].update(last_accepted_x=.45),
    lambda rows: rows[2]["input"].update(alpha=.0002),
    lambda rows: rows[1]["preparation"].update(input_sha256="d" * 64),
    lambda rows: rows[2]["radial_accuracy"].update(rtol=2e-13),
    lambda rows: rows[2]["geometry"].update(cap_cells=3),
])
def test_incomplete_or_changed_inputs_do_not_change_the_register(tmp_path, mutation):
    root = historical_root(tmp_path)
    rows = receipts()
    mutation(rows)
    before = register_bytes(root)
    for relative, row in zip(RELATIVE_LEVELS, rows):
        dump(root / relative, row)
    bindings = []
    for relative in RELATIVE_LEVELS:
        raw = (root / relative).read_bytes()
        bindings.append(dict(path=relative, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    dump(root / COMPARISON, dict(
        schema="nsc-atlas-fixed-endpoint-comparison-v1", status="completed_saved_receipt_comparison",
        inputs=bindings, analysis_source_commit="ab" * 20, analysis_source_sha256="cd" * 32,
        prediction_status=checker.ENVELOPE_STATUS, shared_alpha=.0001,
    ))
    with pytest.raises(ValueError):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before


def test_tampered_bytes_and_promoted_status_preserve_the_register(tmp_path):
    root = historical_root(tmp_path)
    write_comparison(root, receipts())
    before = register_bytes(root)
    receipt = root / RELATIVE_LEVELS[0]
    receipt.write_bytes(receipt.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="evidence bytes changed"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before
    write_comparison(root, receipts())
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    document["prediction_status"] = checker.UNCONVERGED_STATUS
    document["requested_mass_and_flux_contraction"] = False
    dump(root / COMPARISON, document)
    with pytest.raises(ValueError, match="prediction_status"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before
    document = json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))
    document["entries"][0]["coefficients"]["rho"]["value"] = -5.0
    dump(root / checker.REGISTER, document)
    tampered = register_bytes(root)
    with pytest.raises(ValueError, match="rho saved-value projection"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == tampered


def test_recompute_tolerance_accepts_only_a_tight_float_band(tmp_path):
    root = historical_root(tmp_path)
    rows = receipts()
    fresh = compare(rows)["physical_observable"]["value_metres"]
    write_comparison(root, rows)
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    inside = fresh * checker.RECOMPUTE_REL_TOL * 1e-2
    document["physical_observable"]["value_metres"] = fresh + inside
    dump(root / COMPARISON, document)
    registrar.register(root, COMPARISON)
    copied = entry(root)["finite_amplitude"]["physical_observable"]["value_metres"]
    assert copied == pytest.approx(fresh + inside, rel=0, abs=0)
    assert copied != fresh
    assert checker.check(root)["scientific_replays"] == 0

    root = historical_root(tmp_path)
    write_comparison(root, rows)
    before = register_bytes(root)
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    outside = fresh * checker.RECOMPUTE_REL_TOL * 1e2
    document["physical_observable"]["value_metres"] = fresh + outside
    dump(root / COMPARISON, document)
    with pytest.raises(ValueError, match="recompute tolerance"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before

    document["physical_observable"]["value_metres"] = fresh
    document["physical_observable"]["event"]["affine_source_x"] = .48
    dump(root / COMPARISON, document)
    with pytest.raises(ValueError, match="incoming-label event|physical_observable.event"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before

    write_comparison(root, rows)
    before = register_bytes(root)
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    document["inputs"][0]["path"] = "results/summaries/other-endpoint.json"
    dump(root / COMPARISON, document)
    with pytest.raises(ValueError, match="path|evidence"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before

    write_comparison(root, rows)
    before = register_bytes(root)
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    document["inputs"][0]["sha256"] = "e" * 64
    dump(root / COMPARISON, document)
    with pytest.raises(ValueError, match="evidence bytes changed"):
        registrar.register(root, COMPARISON)
    assert register_bytes(root) == before

    write_comparison(root, rows)
    registrar.register(root, COMPARISON)
    document = json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))
    document["entries"][0]["finite_amplitude"]["physical_observable"]["value_metres"] = fresh + inside
    dump(root / checker.REGISTER, document)
    with pytest.raises(ValueError, match="physical observable copy"):
        checker.check(root)


def test_one_ulp_radius_subtraction_is_allowed_and_material_error_is_not(tmp_path):
    root = historical_root(tmp_path)
    rows = receipts(radii=(.55, .55, .55))
    write_comparison(root, rows)
    document = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
    observable = document["physical_observable"]
    radii = [value / observable["scale"]["ell_metres"] for value in observable["grid_values_metres"]]
    pair_tol = checker.subtraction_abs_tol((radii[0], radii[1]))
    half_tol = checker.subtraction_abs_tol(observable["grid_values_metres"])
    quiet_tol = checker.subtraction_abs_tol(
        (observable["value_metres"], observable["quiet_reference_metres"]))
    assert document["endpoint_radius_pair_differences"] == [0.0, 0.0]
    assert observable["resolution_half_width_metres"] == 0.0
    assert observable["difference_from_quiet_reference_metres"] == 0.0
    assert math.ulp(radii[0]) < pair_tol
    assert math.ulp(observable["value_metres"]) < quiet_tol
    document["endpoint_radius_pair_differences"][0] = math.ulp(radii[0])
    observable["resolution_half_width_metres"] = math.ulp(observable["grid_values_metres"][-1])
    observable["difference_from_quiet_reference_metres"] = math.ulp(observable["value_metres"])
    dump(root / COMPARISON, document)
    registrar.register(root, COMPARISON)
    block = entry(root)["finite_amplitude"]
    assert block["endpoint_radius_pair_difference_contracts"] is False
    assert block["prediction_status"] == checker.UNCONVERGED_STATUS
    assert block["requested_mass_and_flux_contraction"] is True
    assert block["limits"]["calibration"] == "none"
    assert block["limits"]["physical_source_selected"] is False
    assert block["limits"]["physical_emitter_selected"] is False
    assert block["physical_observable"]["resolution_half_width_metres"] == (
        math.ulp(observable["grid_values_metres"][-1]))

    def reject_saved_change(mutate, message):
        historical_root(root)
        write_comparison(root, rows)
        before = register_bytes(root)
        changed = json.loads((root / COMPARISON).read_text(encoding="utf-8"))
        mutate(changed)
        dump(root / COMPARISON, changed)
        with pytest.raises(ValueError, match=message):
            registrar.register(root, COMPARISON)
        assert register_bytes(root) == before

    reject_saved_change(
        lambda saved: saved["endpoint_radius_pair_differences"].__setitem__(0, pair_tol * 100),
        "recompute tolerance")
    reject_saved_change(
        lambda saved: saved["physical_observable"].__setitem__(
            "resolution_half_width_metres", half_tol * 100),
        "recompute tolerance")
    reject_saved_change(
        lambda saved: saved["physical_observable"].__setitem__(
            "difference_from_quiet_reference_metres", quiet_tol * 100),
        "recompute tolerance")
    reject_saved_change(
        lambda saved: saved.update(endpoint_radius_pair_difference_contracts=True),
        "endpoint_radius_pair_difference_contracts")
    reject_saved_change(
        lambda saved: saved["residual_maxima"]["wall_flux"]["values"].__setitem__(
            0, saved["residual_maxima"]["wall_flux"]["values"][0] + math.ulp(1.0)),
        "residual_maxima")
    reject_saved_change(
        lambda saved: saved.__setitem__(
            "shared_alpha", saved["shared_alpha"] + math.ulp(saved["shared_alpha"])),
        "shared_alpha")
    reject_saved_change(
        lambda saved: saved["radial_accuracy"].__setitem__(
            "rtol", saved["radial_accuracy"]["rtol"] + math.ulp(saved["radial_accuracy"]["rtol"])),
        "radial_accuracy")


def test_scale_event_and_uncertainty_cannot_be_promoted(tmp_path):
    root = historical_root(tmp_path)
    write_comparison(root, receipts())
    registrar.register(root, COMPARISON)
    original = json.loads((root / checker.REGISTER).read_text(encoding="utf-8"))

    def reject(mutate, message):
        document = copy.deepcopy(original)
        mutate(document["entries"][0])
        dump(root / checker.REGISTER, document)
        with pytest.raises(ValueError, match=message):
            checker.check(root)

    reject(lambda item: item["finite_amplitude"]["physical_observable"]["scale"].update(
        measured_astrophysical_scale=True), "physical observable copy")
    reject(lambda item: item["finite_amplitude"]["physical_observable"]["scale"].update(ell_metres=1.),
           "physical observable copy")
    reject(lambda item: item["finite_amplitude"]["physical_observable"]["uncertainty"].update(
        rigorous_error_enclosure=True), "physical observable copy")
    reject(lambda item: item["finite_amplitude"]["physical_observable"]["uncertainty"].update(
        confidence_level=.95), "physical observable copy")
    reject(lambda item: item["finite_amplitude"]["limits"].update(empirical_confirmation=True),
           "finite-amplitude limits")
    reject(lambda item: item["finite_amplitude"]["limits"].update(
        historical_coefficients_predict_this_event=True), "finite-amplitude limits")
    reject(lambda item: item["finite_amplitude"]["event"].update(distinct_from_perturbative_rho=False),
           "finite-amplitude event")
    reject(lambda item: item["observational_test"].update(SI_scale_selected=False), "observational status")
    reject(lambda item: item["observational_test"].update(status="confirmed"), "observational status")
    reject(lambda item: item.update(id="NSC-P002"), "prediction identity")
    reject(lambda item: item["finite_amplitude"].update(id="NSC-P002"), "finite-amplitude fields")


def test_registration_command_reports_failure_without_writing(tmp_path, capsys):
    root = historical_root(tmp_path)
    before = register_bytes(root)
    missing = root / "results/summaries/missing-comparison.json"
    assert registrar.main(["--root", str(root), "--comparison", str(missing)]) == 1
    assert register_bytes(root) == before
    assert "prediction register update failed" in capsys.readouterr().out
