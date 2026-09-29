"""Resumable middle-window coverage. At most two original rows are solved."""
from hashlib import sha256
from io import BytesIO, StringIO
import json
from pathlib import Path
import runpy
import shutil
from contextlib import redirect_stdout

import numpy as np
import pytest

from recursive_horizons.evidence_io import PrepublicationError
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from recursive_horizons.nsc_middle_source_coverage import (
    ENERGY_LOWER, ENERGY_UPPER, EXPECTED_POSITIVE_ROWS, PANEL, RTOL, START,
    checkpoint_name, cover, encode_record, in_middle_window, middle_catalogue,
    publish_checkpoint,
)
import recursive_horizons.nsc_middle_source_coverage as coverage

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT/"results"/"development"/"nsc-middle-source-coverage-v1"
SOLVES = {"count": 0}


@pytest.fixture(scope="module")
def archive():
    return RetainedUpstreamArchive(ROOT)


def test_window_is_the_original_middle_census(archive):
    rows = middle_catalogue(archive)
    assert len(rows) == EXPECTED_POSITIVE_ROWS == 48
    assert [row.row for row in rows] == list(range(48))
    assert {row.panel for row in rows} == {PANEL}
    assert rows[0].positive_angular_sign == 1 and rows[0].negative_angular_sign == -1
    energies = [row.energy for row in rows]
    assert energies == sorted(energies) and len(set(energies)) == 48
    assert in_middle_window(ENERGY_LOWER) and not in_middle_window(ENERGY_UPPER)
    assert all(ENERGY_LOWER <= energy < ENERGY_UPPER for energy in energies)
    assert all(row.positive_columns.shape == (2, 3) and row.negative_columns.shape == (2, 3)
               for row in rows)
    assert all(row.positive_covariance.shape == (3, 3) and row.negative_covariance.shape == (3, 3)
               for row in rows)
    assert all(row.negative_energy == -row.energy for row in rows)
    low = high = 0
    for batch, _channel in archive.family_entries((14, 1)):
        if batch.energy_sign <= 0:
            continue
        for energy in batch.source.energies[::3]:
            if float(energy) < 16:
                low += 1
            elif float(energy) >= 32:
                high += 1
    assert low == 152 and high == 384
    published = json.loads((ROOT/"results/development/nsc-vacuum-source-correction-v1.json").read_text())
    assert rows[0].row == published["row"] == 0
    assert float(rows[0].energy).hex() == published["energy_hex"]
    assert float(START).hex() == published["start_log_delta_hex"]
    assert float(RTOL).hex() == published["rtol_hex"]
    with pytest.raises(ValueError):
        in_middle_window(float("nan"))
    seen = set()
    coverage._remember(seen, PANEL, 0, float(rows[0].energy).hex())
    with pytest.raises(ValueError, match="duplicate"):
        coverage._remember(seen, PANEL, 0, float(rows[1].energy).hex())


def test_publication_and_budgets_do_not_solve(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    payload = b"witness-bytes"
    name = checkpoint_name(PANEL, 3)
    record = {"panel": PANEL, "row": 3, "payload": {
        "path": f"rows/{name}/witness.npz", "sha256": sha256(payload).hexdigest(),
        "bytes": len(payload)}}
    publish_checkpoint(root, record, payload)
    with pytest.raises(PrepublicationError, match="destination exists"):
        publish_checkpoint(root, record, payload)
    assert (root/"rows"/name/"witness.npz").read_bytes() == payload

    def explode(*_args, **_kwargs):
        raise AssertionError("budgeted resume solved an original row")

    monkeypatch.setattr(coverage, "capture_correction", explode)
    absent = tmp_path/"absent-check"
    checked = cover(ROOT, absent, mode="check")
    assert not absent.exists()
    assert checked["status"].startswith("OPEN")
    assert checked["completed_positive_rows"] == 0
    assert checked["new_rows_captured"] == 0
    assert checked["physical_upstream_budget_component"] is None
    assert checked["missing_positive_rows"] == 48
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

    driver = runpy.run_path(str(ROOT/"scripts/derive_nsc_middle_source_coverage.py"))
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
        assert 0 < kwargs["cpu_limit"] <= 30
        raise TimeoutError("vacuum-correction proof-witness CPU budget exhausted")

    monkeypatch.setattr(coverage, "capture_correction", timed_out)
    stopped = cover(ROOT, tmp_path/"timed-out", mode="resume", row_budget=1, cpu_budget=30)
    assert stopped["new_rows_captured"] == 0
    assert stopped["status"].startswith("OPEN")
    assert stopped["maximum_completed_source_operator_error_upper"] is None
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
        return real_capture(*args, **kwargs)

    def validate(*args, **kwargs):
        events.append("validate")
        return real_validate(*args, **kwargs)

    monkeypatch.setattr(coverage, "capture_correction", capture)
    monkeypatch.setattr(coverage, "validate_correction", validate)


def _copy_row(root, destination, row):
    source = root/"rows"/checkpoint_name(PANEL, row)
    target = destination/"rows"/checkpoint_name(PANEL, row)
    shutil.copytree(source, target)
    return target


def _expect_refusal(destination, match):
    with pytest.raises(ValueError, match=match):
        cover(ROOT, destination, mode="check")
    with pytest.raises(ValueError, match=match):
        cover(ROOT, destination, mode="resume", row_budget=1, cpu_budget=30)


def test_two_original_rows_then_replay_resume_and_refusals(tmp_path, monkeypatch):
    assert SOLVES["count"] == 0
    root = tmp_path.resolve()/"coverage"
    events = []
    real_capture = coverage.capture_correction
    real_validate = coverage.validate_correction
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=True)
    first = cover(ROOT, root, mode="resume", row_budget=1, cpu_budget=120)
    assert events == ["capture", "validate"]
    assert SOLVES["count"] == 1
    assert first["new_rows_captured"] == 1
    assert first["completed_positive_rows"] == 1
    assert first["completed_negative_partners"] == 1
    assert first["completed_signed_rows"] == 2
    assert first["missing_positive_rows"] == 47
    assert first["status"].startswith("OPEN")
    assert first["coverage_complete"] is False
    assert first["physical_upstream_budget_component"] is None
    assert first["physical_local_gate"] == "OPEN"
    assert set(first["missing"][0]) == {"panel", "row", "energy_hex"}
    row = first["rows"][0]
    assert row["equation"] == coverage.EQUATION
    assert restored_upper(row["horizon_initial_remainder_upper"]) > 0
    assert restored_upper(row["thermal_coherent_operator_error_upper"]) > 0
    assert row["source_comparisons"][0]["energy_sign"] == 1
    assert row["source_comparisons"][1]["energy_sign"] == -1
    assert row["source_comparisons"][1]["angular_sign"] == -1
    assert (row["source_comparisons"][0]["physical_source_operator_error_upper"]
            != row["source_comparisons"][1]["physical_source_operator_error_upper"])
    published = json.loads((ROOT/"results/development/nsc-vacuum-source-correction-v1.json").read_text())
    saved = json.loads((root/"rows"/checkpoint_name(PANEL, 0)/"record.json").read_text())
    for field in ("vacuum_bloch_error_upper", "initial_bloch_error_upper",
                  "thermal_coherent_operator_error_upper", "source_comparisons",
                  "cells", "nfev", "energy_hex", "panel", "row"):
        assert saved[field] == published[field], field
    assert saved["horizon_initial_remainder_nonzero"] is True
    assert saved["positive_scalar_bound_copied_to_negative"] is False
    assert saved["negative_comparison_uses_complemented_endpoint"] is True
    assert saved["archived_source_replaced"] is False
    assert saved["source_occupations_changed"] is False

    events.clear()
    _install_solver(monkeypatch, events, real_capture, real_validate, allow_capture=False)
    replayed = cover(ROOT, root, mode="check")
    assert events == ["validate"]
    assert replayed["new_rows_captured"] == 0
    assert replayed["completed_positive_rows"] == 1
    events.clear()
    # A saved row is validated, not solved again. Budget zero asks for no new row.
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
    assert second["missing_positive_rows"] == 46
    assert second["status"].startswith("OPEN")
    assert second["physical_upstream_budget_component"] is None
    assert [item["row"] for item in second["rows"]] == [0, 1]
    assert all(restored_upper(item["horizon_initial_remainder_upper"]) > 0 for item in second["rows"])
    assert all(restored_upper(item["thermal_coherent_operator_error_upper"]).is_finite()
               for item in second["rows"])
    assert all(item["source_comparisons"][0]["physical_source_operator_error_upper"]
               != item["source_comparisons"][1]["physical_source_operator_error_upper"]
               for item in second["rows"])
    uppers = [restored_upper(comparison["physical_source_operator_error_upper"])
              for item in second["rows"] for comparison in item["source_comparisons"]]
    assert second["maximum_completed_source_operator_error_upper"] == exact_upper(max(uppers))

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

    def refuse(name, mutate, match):
        destination = tmp_path.resolve()/name
        target = _copy_row(root, destination, 0)
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

    def change_dependency(target):
        record = json.loads((target/"record.json").read_text())
        key = "src/recursive_horizons/nsc_vacuum_source_correction.py"
        digest = record["source_hashes"][key]
        record["source_hashes"][key] = digest[:-1]+("0" if digest[-1] != "0" else "1")
        (target/"record.json").write_bytes(encode_record(record))

    refuse("missing-payload", drop_payload, "payload")
    refuse("changed-payload", change_payload, "payload")
    refuse("changed-settings", change_settings, "settings")
    refuse("changed-source", change_source, "source")
    refuse("changed-dependency", change_dependency, "dependency")

    duplicate = tmp_path.resolve()/"duplicate"
    _copy_row(root, duplicate, 0)
    copied = _copy_row(root, duplicate, 1)
    record = json.loads((copied/"record.json").read_text())
    record["row"] = 0
    record["energy_hex"] = saved["energy_hex"]
    (copied/"record.json").write_bytes(encode_record(record))
    events.clear()
    _expect_refusal(duplicate, "duplicate")
    assert events == []
    assert SOLVES["count"] == 2
    assert not CAMPAIGN.exists()
