"""Resumable direct-vacuum window. At most two original rows are solved."""
from hashlib import sha256
from io import BytesIO, StringIO
import json
from pathlib import Path
import runpy
import shutil
from contextlib import redirect_stdout

import numpy as np
import pytest
from flint import arb, ctx

from recursive_horizons.evidence_io import PrepublicationError
from recursive_horizons.nsc_direct_source_window import (
    EQUATION, EXPECTED_POSITIVE_ROWS, LOW16_PANEL, LOW32_PANEL, MAX_STEP, SCHEMA,
    START, aggregate, checkpoint_name, cover, direct_catalogue, encode_record,
    in_direct_window, publish_checkpoint,
)
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
import recursive_horizons.nsc_direct_source_window as window

ROOT = Path(__file__).resolve().parents[1]
PRODUCTION = (
    ROOT/"results"/"development"/"nsc-direct-source-window-v1",
    ROOT/"results"/"development"/"nsc-middle-source-coverage-v1",
)
PUBLISHED_COVERAGE = ROOT/"results"/"development"/"nsc-ks-source-operator-majorant-v3.json"
SOLVES = {"count": 0}


def _snapshot(path):
    path = Path(path)
    if not path.exists():
        return None
    digest = sha256()
    if path.is_file():
        digest.update(path.read_bytes())
        return digest.hexdigest()
    for item in sorted(path.rglob("*"), key=lambda candidate: candidate.relative_to(path).as_posix()):
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b"\0")
        if item.is_symlink():
            digest.update(b"link")
        elif item.is_file():
            digest.update(item.read_bytes())
        else:
            digest.update(b"dir")
        digest.update(b"\0")
    return digest.hexdigest()


@pytest.fixture(scope="module", autouse=True)
def production_directories_stay_byte_identical():
    before = {str(path): _snapshot(path) for path in PRODUCTION}
    yield
    after = {str(path): _snapshot(path) for path in PRODUCTION}
    assert after == before


@pytest.fixture(scope="module")
def archive():
    return RetainedUpstreamArchive(ROOT)


def test_window_is_the_original_low_energy_census(archive):
    source = (ROOT/"src/recursive_horizons/nsc_direct_source_window.py").read_text()
    assert "VacuumCorrection" not in source
    assert "capture_correction" not in source
    assert "nsc_middle_source_coverage" not in source
    script = (ROOT/"scripts/derive_nsc_direct_source_window.py").read_text()
    assert "nsc-middle-source-coverage" not in script
    assert "nsc-direct-source-window-v1" in script

    rows = direct_catalogue(archive)
    assert len(rows) == EXPECTED_POSITIVE_ROWS == 87
    low16 = sorted(row.row for row in rows if row.panel == LOW16_PANEL)
    low32 = sorted(row.row for row in rows if row.panel == LOW32_PANEL)
    assert low16 == list(range(16, 48)) and len(low16) == 32
    assert low32 == list(range(41, 96)) and len(low32) == 55
    assert {row.panel for row in rows} == {LOW16_PANEL, LOW32_PANEL}
    energies = [row.energy for row in rows]
    assert energies == sorted(energies) and len(set(energies)) == 87
    assert rows[0].energy == 2.071557525252559
    assert rows[-1].energy == 15.9788018699833
    assert rows[0].panel == LOW32_PANEL and rows[0].row == 41
    assert rows[-1].panel == LOW16_PANEL and rows[-1].row == 47
    assert in_direct_window(2.0) and not in_direct_window(16.0)
    assert all(in_direct_window(energy) for energy in energies)
    assert all(row.positive_columns.shape == (2, 3) and row.negative_columns.shape == (2, 3)
               for row in rows)
    assert all(row.positive_covariance.shape == (3, 3) and row.negative_covariance.shape == (3, 3)
               for row in rows)
    assert all(row.negative_energy_fiber == tuple(-energy for energy in row.energy_fiber)
               for row in rows)
    assert all(len(row.energy_fiber) == 3 and len(set(row.energy_fiber)) == 1 for row in rows)
    assert all(row.positive_angular_sign == 1 and row.negative_angular_sign == -1 for row in rows)
    assert all(not np.array_equal(row.positive_columns, row.negative_columns) for row in rows)
    assert all(row.weight > 0 and row.column_weight > 0 for row in rows)
    assert all(not (row.panel == LOW16_PANEL and row.row in (0, 15)) for row in rows)
    with pytest.raises(ValueError):
        in_direct_window(float("nan"))
    seen = set()
    window._remember(seen, rows[0].panel, rows[0].row, float(rows[0].energy).hex())
    with pytest.raises(ValueError, match="duplicate"):
        window._remember(seen, rows[0].panel, rows[0].row, float(rows[1].energy).hex())

    published = json.loads(PUBLISHED_COVERAGE.read_text())
    covered = {(item["panel"], item["row"], item["energy_sign"])
               for item in published["source_coverage"]["rows"]}
    assert len(covered) == 868
    window_keys = {(row.panel, row.row, sign) for row in rows for sign in (1, -1)}
    assert len(window_keys) == 174 and window_keys.isdisjoint(covered)


def test_finished_finite_window_stays_open(archive):
    catalogue = direct_catalogue(archive)
    upper = exact_upper(arb(2)**-10)
    completed = []
    for spec in catalogue:
        positive = {"energy_sign": 1, "angular_sign": 1, "signed_total_upper": upper}
        negative = {"energy_sign": -1, "angular_sign": -1, "signed_total_upper": upper}
        completed.append({
            "panel": spec.panel, "row": spec.row,
            "energy_fiber_hex": [float(spec.energy).hex()]*3,
            "negative_energy_fiber_hex": [float(spec.negative_energy).hex()]*3,
            "payload": {"path": f"rows/{checkpoint_name(spec.panel, spec.row)}/witness.npz"},
            "cells": 1, "nfev": 1, "equation": EQUATION,
            "vacuum_bloch_error_upper": upper,
            "occupation": {"operator_distance_upper": upper},
            "signed_comparisons": [positive, negative],
            "weight_hex": float(spec.weight).hex(),
            "column_weight_hex": float(spec.column_weight).hex(),
        })
    common = dict(mode="check", row_budget=None, cpu_budget=None, cpu_seconds=0.0,
                  capture_cpu_seconds=0.0, new_rows_captured=0)
    partial = aggregate(catalogue, completed[:-1], **common)
    assert partial["status"].startswith("OPEN")
    assert "ENCLOSED" not in partial["status"]
    assert partial["coverage_complete"] is False
    assert partial["missing_positive_rows"] == 1
    assert partial["completed_signed_rows"] == 172
    assert set(partial["missing"][0]) == {"panel", "row", "energy_hex"}
    assert partial["maximum_completed_signed_total_upper"] == upper
    assert partial["maximum_completed_signed_total_upper"] != exact_upper(arb(0))
    assert partial["physical_upstream_budget_component"] is None
    assert partial["physical_local_gate"] == "OPEN"
    assert partial["all_source_families"] is False
    assert partial["weight_applied_to_error"] is False
    assert partial["correction_forcing_used"] is False
    assert partial["finite_window_closes_gate"] is False
    finished = aggregate(catalogue, completed, **common)
    assert finished["coverage_complete"] is True
    assert finished["missing"] == []
    assert finished["status"].startswith("OPEN")
    assert "ENCLOSED" not in finished["status"]
    assert finished["physical_local_gate"] == "OPEN"
    assert finished["physical_upstream_budget_component"] is None
    assert finished["finite_window_closes_gate"] is False
    assert finished["schema"] == SCHEMA
    broken = dict(completed[-1])
    broken["signed_comparisons"] = [broken["signed_comparisons"][0]]
    with pytest.raises(ValueError, match="missing bound"):
        aggregate(catalogue, completed[:-1]+[broken], **common)


def test_publication_budgets_and_frozen_campaign_do_not_solve(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    payload = b"witness-bytes"
    name = checkpoint_name(LOW32_PANEL, 41)
    record = {"panel": LOW32_PANEL, "row": 41, "payload": {
        "path": f"rows/{name}/witness.npz", "sha256": sha256(payload).hexdigest(),
        "bytes": len(payload)}}
    publish_checkpoint(root, record, payload)
    with pytest.raises(PrepublicationError, match="destination exists"):
        publish_checkpoint(root, record, payload)
    assert (root/"rows"/name/"witness.npz").read_bytes() == payload
    before_middle = _snapshot(PRODUCTION[1])
    with pytest.raises(ValueError, match="frozen middle-source"):
        cover(ROOT, PRODUCTION[1], mode="check")
    assert _snapshot(PRODUCTION[1]) == before_middle

    def explode(*_args, **_kwargs):
        raise AssertionError("budgeted resume solved an original row")

    monkeypatch.setattr(window, "capture_direct_vacuum", explode)
    absent = tmp_path/"absent-check"
    checked = cover(ROOT, absent, mode="check")
    assert not absent.exists()
    assert checked["status"].startswith("OPEN")
    assert checked["completed_positive_rows"] == 0
    assert checked["new_rows_captured"] == 0
    assert checked["physical_upstream_budget_component"] is None
    assert checked["finite_window_closes_gate"] is False
    assert checked["missing_positive_rows"] == 87
    assert checked["maximum_completed_signed_total_upper"] is None
    assert set(checked["missing"][0]) == {"panel", "row", "energy_hex"}
    paused = cover(ROOT, tmp_path/"paused", mode="resume", row_budget=0, cpu_budget=1e-9)
    assert paused["new_rows_captured"] == 0
    assert paused["coverage_complete"] is False
    assert paused["status"].startswith("OPEN")
    assert not (tmp_path/"paused"/"rows").exists()
    with pytest.raises(ValueError, match="row budget"):
        cover(ROOT, tmp_path/"bad-budget", mode="resume", row_budget=-1, cpu_budget=1)
    with pytest.raises(ValueError, match="CPU budget"):
        cover(ROOT, tmp_path/"bad-cpu", mode="resume", row_budget=1, cpu_budget=0)
    with pytest.raises(ValueError, match="solve budget"):
        cover(ROOT, tmp_path/"bad-check", mode="check", row_budget=1, cpu_budget=1)

    driver = runpy.run_path(str(ROOT/"scripts/derive_nsc_direct_source_window.py"))
    assert driver["DEFAULT_OUTPUT"] == "results/development/nsc-direct-source-window-v1"
    cli_absent = tmp_path/"cli-absent"
    with redirect_stdout(StringIO()):
        assert driver["main"](["--check", "--output", str(cli_absent)]) == 0
    assert not cli_absent.exists()
    with pytest.raises(SystemExit) as caught:
        driver["main"](["--resume", "--output", str(tmp_path/"cli-resume")])
    assert caught.value.code == 2
    with pytest.raises(SystemExit) as caught:
        driver["main"](["--check", "--row-budget", "1", "--cpu-budget", "1",
                        "--output", str(tmp_path/"cli-check")])
    assert caught.value.code == 2

    def timed_out(*_args, **kwargs):
        assert 0 < kwargs["cpu_limit"] <= 120
        raise TimeoutError("direct vacuum Bloch CPU budget exhausted")

    monkeypatch.setattr(window, "capture_direct_vacuum", timed_out)
    stopped = cover(ROOT, tmp_path/"timed-out", mode="resume", row_budget=1, cpu_budget=120)
    assert stopped["new_rows_captured"] == 0
    assert stopped["status"].startswith("OPEN")
    assert stopped["maximum_completed_signed_total_upper"] is None
    assert stopped["physical_upstream_budget_component"] is None
    assert not (tmp_path/"timed-out"/"rows").exists()


def _install_solver(monkeypatch, events, real_capture, real_validate, *, allow_capture):
    def capture(*args, **kwargs):
        events.append("capture")
        if not allow_capture:
            raise AssertionError("saved replay or resume captured a row")
        SOLVES["count"] += 1
        if SOLVES["count"] > 2:
            raise AssertionError("more than two original rows were solved")
        assert kwargs["cpu_limit"] > 0
        assert kwargs["max_step"] == MAX_STEP
        assert args[1] == START
        assert args[0].bits == 192 and args[0].metric_terms == 48
        return real_capture(*args, **kwargs)

    def validate(*args, **kwargs):
        events.append("validate")
        assert kwargs.get("degree") == 12
        return real_validate(*args, **kwargs)

    monkeypatch.setattr(window, "capture_direct_vacuum", capture)
    monkeypatch.setattr(window, "validate_direct_vacuum", validate)


def _copy_row(root, destination, panel, row):
    source = root/"rows"/checkpoint_name(panel, row)
    target = destination/"rows"/checkpoint_name(panel, row)
    shutil.copytree(source, target)
    return target


def _row_bytes(root):
    rows = root/"rows"
    return {path.name: ((path/"record.json").read_bytes(), (path/"witness.npz").read_bytes())
            for path in sorted(rows.iterdir(), key=lambda item: item.name)}


def _expect_refusal(destination, match):
    with pytest.raises(ValueError, match=match):
        cover(ROOT, destination, mode="check")
    with pytest.raises(ValueError, match=match):
        cover(ROOT, destination, mode="resume", row_budget=1, cpu_budget=30)


def _totals_match(record):
    with ctx.workprec(192):
        occupation = restored_upper(record["occupation"]["operator_distance_upper"])
        assert occupation > 0
        for comparison in record["signed_comparisons"]:
            archived = restored_upper(comparison["archived_distance_upper"])
            total = restored_upper(comparison["signed_total_upper"])
            assert exact_upper((archived+occupation).upper()) == comparison["signed_total_upper"]
            assert total > archived
            assert comparison["signed_total_lower"] == comparison["archived_distance_lower"]


def test_two_original_rows_then_replay_resume_and_refusals(tmp_path, monkeypatch):
    assert SOLVES["count"] == 0
    root = tmp_path.resolve()/"window"
    events = []
    real_capture = window.capture_direct_vacuum
    real_validate = window.validate_direct_vacuum
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=True)
    first = cover(ROOT, root, mode="resume", row_budget=1, cpu_budget=120)
    assert events == ["capture", "validate"]
    assert SOLVES["count"] == 1
    assert first["new_rows_captured"] == 1
    assert first["completed_positive_rows"] == 1
    assert first["completed_negative_partners"] == 1
    assert first["completed_signed_rows"] == 2
    assert first["missing_positive_rows"] == 86
    assert first["status"].startswith("OPEN")
    assert first["coverage_complete"] is False
    assert first["finite_window_closes_gate"] is False
    assert first["physical_upstream_budget_component"] is None
    assert first["physical_local_gate"] == "OPEN"
    assert first["correction_forcing_used"] is False
    assert first["weight_applied_to_error"] is False
    assert set(first["missing"][0]) == {"panel", "row", "energy_hex"}
    row = first["rows"][0]
    assert row["panel"] == LOW32_PANEL and row["row"] == 41
    assert row["equation"] == EQUATION
    assert row["energy_fiber_hex"] == [float(2.071557525252559).hex()]*3
    assert row["negative_energy_fiber_hex"] == [float(-2.071557525252559).hex()]*3
    assert [item["energy_sign"] for item in row["signed_comparisons"]] == [1, -1]
    assert [item["angular_sign"] for item in row["signed_comparisons"]] == [1, -1]
    assert row["occupation"]["incoming_gap_used"] is True
    assert row["occupation"]["source_occupation_law_changed"] is False
    _totals_match(row)
    assert (row["signed_comparisons"][0]["archived_distance_upper"]
            != row["signed_comparisons"][1]["archived_distance_upper"])
    saved = json.loads((root/"rows"/checkpoint_name(LOW32_PANEL, 41)/"record.json").read_text())
    assert saved["schema"] == "NSC-DIRECT-SOURCE-ROW-v1"
    assert saved["coverage_schema"] == SCHEMA
    assert saved["correction_forcing_used"] is False
    assert saved["replay_uses_saved_trace"] is True
    assert saved["positive_scalar_bound_copied_to_negative"] is False
    assert saved["negative_comparison_uses_complemented_endpoint"] is True
    assert saved["archived_source_replaced"] is False
    assert saved["initial_archived_arrays_changed"] is False
    assert saved["weight_applied_to_error"] is False
    assert saved["covariance_errors_unweighted"] is True
    assert saved["frame_order"] == 16 and saved["bits"] == 192
    assert saved["metric_terms"] == 48 and saved["defect_degree"] == 12
    assert saved["max_step_hex"] == float(MAX_STEP).hex()
    assert saved["start_log_delta_hex"] == float(START).hex()
    assert saved["rtol_hex"] == float(2e-13).hex() and saved["atol_hex"] == float(2e-15).hex()
    assert saved["analytic_radius"] == "0.1"
    assert len(saved["energy_fiber_hex"]) == 3 and len(saved["negative_energy_fiber_hex"]) == 3
    assert saved["positive_source_digest"] != saved["negative_source_digest"]
    assert saved["positive_preparation_digest"] != saved["negative_preparation_digest"]
    assert saved["group"] == 14 and saved["mass_hex"] and saved["angular_hex"]
    assert saved["rho_up_hex"] and saved["horizon_rho_hex"] and saved["kappa_hex"] and saved["omega_hex"]
    assert saved["archive_digests"] and saved["source_hashes"]
    assert saved["weight_hex"] and saved["column_weight_hex"]
    with np.load(root/"rows"/checkpoint_name(LOW32_PANEL, 41)/"witness.npz") as payload:
        assert not np.array_equal(payload["positive_columns"], payload["negative_columns"])
        assert not np.array_equal(payload["positive_covariance"], payload["negative_covariance"])
        assert payload["trace"].shape[1] == 26 and len(payload["trace"]) == saved["cells"]

    events.clear()
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=False)
    replayed = cover(ROOT, root, mode="check")
    assert events == ["validate"]
    assert replayed["new_rows_captured"] == 0
    assert replayed["completed_positive_rows"] == 1
    events.clear()
    resumed = cover(ROOT, root, mode="resume", row_budget=0, cpu_budget=120)
    assert events == ["validate"]
    assert resumed["new_rows_captured"] == 0
    assert SOLVES["count"] == 1

    events.clear()
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=True)
    second = cover(ROOT, root, mode="resume", row_budget=1, cpu_budget=120)
    assert events == ["validate", "capture", "validate"]
    assert SOLVES["count"] == 2
    assert second["new_rows_captured"] == 1
    assert second["completed_positive_rows"] == 2
    assert second["completed_signed_rows"] == 4
    assert second["missing_positive_rows"] == 85
    assert second["status"].startswith("OPEN")
    assert second["physical_upstream_budget_component"] is None
    assert [item["row"] for item in second["rows"]] == [41, 42]
    assert [item["panel"] for item in second["rows"]] == [LOW32_PANEL, LOW32_PANEL]
    assert all(item["occupation"]["incoming_gap_used"] is True for item in second["rows"])
    assert all(item["signed_comparisons"][0]["archived_distance_upper"]
               != item["signed_comparisons"][1]["archived_distance_upper"]
               for item in second["rows"])
    for item in second["rows"]:
        _totals_match(item)
    uppers = [restored_upper(comparison["signed_total_upper"])
              for item in second["rows"] for comparison in item["signed_comparisons"]]
    assert second["maximum_completed_signed_total_upper"] == exact_upper(max(uppers))

    events.clear()
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=False)
    again = cover(ROOT, root, mode="resume", row_budget=0, cpu_budget=120)
    assert events == ["validate", "validate"]
    assert again["new_rows_captured"] == 0
    assert SOLVES["count"] == 2
    events.clear()
    checked = cover(ROOT, root, mode="check")
    assert events == ["validate", "validate"]
    assert checked["new_rows_captured"] == 0
    assert checked["completed_positive_rows"] == 2

    saved_bytes = _row_bytes(root)
    events.clear()

    def timed_out(*_args, **_kwargs):
        events.append("capture")
        raise TimeoutError("direct vacuum Bloch CPU budget exhausted")

    monkeypatch.setattr(window, "capture_direct_vacuum", timed_out)
    stopped = cover(ROOT, root, mode="resume", row_budget=1, cpu_budget=120)
    assert events == ["validate", "validate", "capture"]
    assert stopped["new_rows_captured"] == 0
    assert stopped["completed_positive_rows"] == 2
    assert _row_bytes(root) == saved_bytes
    assert SOLVES["count"] == 2

    events.clear()
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=False)

    def refuse(name, mutate, match):
        destination = tmp_path.resolve()/name
        target = _copy_row(root, destination, LOW32_PANEL, 41)
        mutate(target)
        before = list(events)
        events.clear()
        _expect_refusal(destination, match)
        assert events == []
        events.extend(before)

    def drop_payload(target):
        (target/"witness.npz").unlink()

    def change_payload(target):
        payload = bytearray((target/"witness.npz").read_bytes())
        payload[-1] ^= 0x01
        (target/"witness.npz").write_bytes(bytes(payload))

    def change_settings(target):
        record = json.loads((target/"record.json").read_text())
        record["rtol_hex"] = float(1e-8).hex()
        (target/"record.json").write_bytes(encode_record(record))

    def change_source(target):
        raw = (target/"witness.npz").read_bytes()
        with np.load(BytesIO(raw), allow_pickle=False) as data:
            arrays = {key: data[key].copy() for key in data.files}
        arrays["positive_columns"][0, 0] += 0.01
        payload = deterministic_npz_bytes(arrays)
        record = json.loads((target/"record.json").read_text())
        record["payload"]["sha256"] = sha256(payload).hexdigest()
        record["payload"]["bytes"] = len(payload)
        (target/"witness.npz").write_bytes(payload)
        (target/"record.json").write_bytes(encode_record(record))

    def drop_bound(target):
        record = json.loads((target/"record.json").read_text())
        del record["positive_preparation_digest"]
        (target/"record.json").write_bytes(encode_record(record))

    def change_bound(target):
        record = json.loads((target/"record.json").read_text())
        key = next(iter(record["archive_digests"]))
        digest = record["archive_digests"][key]
        record["archive_digests"][key] = digest[:-1]+("0" if digest[-1] != "0" else "1")
        (target/"record.json").write_bytes(encode_record(record))

    def corrupt_record(target):
        (target/"record.json").write_bytes(b"{\n")

    def extra_file(target):
        (target/"note.txt").write_bytes(b"x")

    refuse("missing-payload", drop_payload, "payload")
    refuse("changed-payload", change_payload, "payload")
    refuse("changed-settings", change_settings, "settings")
    refuse("changed-source", change_source, "source")
    refuse("missing-bound", drop_bound, "missing bound")
    refuse("changed-bound", change_bound, "changed bound")
    refuse("corrupt-record", corrupt_record, "corrupt")
    refuse("corrupt-extra", extra_file, "corrupt")

    duplicate = tmp_path.resolve()/"duplicate"
    _copy_row(root, duplicate, LOW32_PANEL, 41)
    copied = _copy_row(root, duplicate, LOW32_PANEL, 42)
    record = json.loads((copied/"record.json").read_text())
    record["row"] = 41
    record["energy_fiber_hex"] = saved["energy_fiber_hex"]
    (copied/"record.json").write_bytes(encode_record(record))
    events.clear()
    _expect_refusal(duplicate, "duplicate")
    assert events == []
    assert SOLVES["count"] == 2
    print({
        "first_capture_cpu": first["capture_cpu_seconds"],
        "first_total_cpu": first["cpu_seconds"],
        "second_capture_cpu": second["capture_cpu_seconds"],
        "second_total_cpu": second["cpu_seconds"],
        "two_row_replay_cpu": checked["cpu_seconds"],
        "rows": [{
            "panel": item["panel"], "row": item["row"],
            "energy": float.fromhex(item["energy_fiber_hex"][0]),
            "cells": item["cells"], "nfev": item["nfev"],
            "vacuum_error": float(restored_upper(item["vacuum_bloch_error_upper"])),
            "occupation": float(restored_upper(item["occupation"]["operator_distance_upper"])),
            "positive_archived": float(restored_upper(
                item["signed_comparisons"][0]["archived_distance_upper"])),
            "negative_archived": float(restored_upper(
                item["signed_comparisons"][1]["archived_distance_upper"])),
            "positive_total": float(restored_upper(
                item["signed_comparisons"][0]["signed_total_upper"])),
            "negative_total": float(restored_upper(
                item["signed_comparisons"][1]["signed_total_upper"])),
            "weight": float.fromhex(item["weight_hex"]),
        } for item in second["rows"]],
    })
