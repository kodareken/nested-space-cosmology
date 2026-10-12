"""Actual18 provenance, signed bounds, and checkpoint refusal. No nine-row solve."""
from hashlib import sha256
from io import BytesIO
import json
import math
from pathlib import Path
import runpy
import sys

import numpy as np
import pytest
from flint import arb

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from nsc_threshold18_source import (  # noqa: E402
    ADJUDICATION, BASE_COMMIT, EQUATION, EXPECTED_POSITIVE_ROWS, MAX_ENERGY_HEX,
    MAX_STEP, MIN_ENERGY_HEX, PANEL, PUBLISHED_SIGNED_ROWS, ROWS, SCHEMA, START,
    aggregate, binding_record, checkpoint_name, cover, deterministic_npz_bytes,
    encode_record, in_actual18, in_fixed_direct_window, load_context,
)
from recursive_horizons.evidence_io import PrepublicationError  # noqa: E402
from recursive_horizons.nsc_direct_vacuum_source import archived_vacuum_distance  # noqa: E402
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper  # noqa: E402
from recursive_horizons.nsc_ks_signed_state import (  # noqa: E402
    s3_conjugate, source_complement_residual,
)
import nsc_threshold18_source as owner  # noqa: E402

PARENT = ROOT.parent / "lab"
PIN = ROOT / "pin" / "lab"
WINDOW = PARENT / "src/recursive_horizons/nsc_direct_source_window.py"
PUBLISHED = PIN / "results/development/nsc-ks-source-operator-majorant-v3.json"


@pytest.fixture(scope="module")
def context():
    return load_context(PIN)


def test_recorder_reuses_the_pinned_solver_and_leaves_the_direct_window_fixed():
    source = (ROOT / "src/nsc_threshold18_source.py").read_text()
    assert "capture_direct_vacuum" in source
    assert "solve_ivp" not in source
    assert "VacuumCorrection" not in source
    assert "nsc_direct_source_window" not in source
    assert "nsc_middle_source_coverage" not in source
    frozen = WINDOW.read_text()
    assert "ENERGY_LOWER = 2.0" in frozen
    assert "ENERGY_UPPER = 16.0" in frozen
    assert "LOW32_ROWS = tuple(range(41, 96))" in frozen
    assert "range(32, 41)" not in frozen
    vacuum = (PIN / "src/recursive_horizons/nsc_direct_vacuum_source.py").read_bytes()
    occupation = (PIN / "src/recursive_horizons/nsc_source_occupation_enclosure.py").read_bytes()
    assert vacuum == (PARENT / "src/recursive_horizons/nsc_direct_vacuum_source.py").read_bytes()
    assert occupation == (PARENT / "src/recursive_horizons/nsc_source_occupation_enclosure.py").read_bytes()
    head = (ROOT / "pin" / ".git" / "HEAD").read_text().strip()
    assert BASE_COMMIT.startswith("eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9"[:12])
    assert head == "eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9" or head.startswith("ref:")


def test_original_signed_bounds_and_provenance(context):
    _lab, catalogue, dependencies = context
    assert len(catalogue) == EXPECTED_POSITIVE_ROWS == 9
    assert [row.row for row in catalogue] == list(ROWS) == list(range(32, 41))
    assert {row.panel for row in catalogue} == {PANEL}
    assert owner._hex(catalogue[0].energy) == MIN_ENERGY_HEX
    assert owner._hex(catalogue[-1].energy) == MAX_ENERGY_HEX
    assert [row.energy for row in catalogue] == sorted(row.energy for row in catalogue)
    assert all(in_actual18(row.energy) and not in_fixed_direct_window(row.energy) for row in catalogue)
    assert in_actual18(math.pi / 2) and not in_actual18(2.0)
    assert in_fixed_direct_window(2.0) and not in_fixed_direct_window(math.nextafter(2.0, 0.0))
    for row in catalogue:
        assert row.negative_energy_fiber == tuple(-energy for energy in row.energy_fiber)
        assert len(set(row.energy_fiber)) == 1
        assert row.positive_angular_sign == 1 and row.negative_angular_sign == -1
        assert np.array_equal(row.negative_columns, s3_conjugate(row.positive_columns))
        assert source_complement_residual(row.positive_covariance, row.negative_covariance) == 0.0
        assert row.weight > 0 and row.column_weight > 0
        assert row.column_weight == float(np.sqrt(row.weight / (2 * np.pi)))
        assert row.positive_source_digest != row.negative_source_digest
        assert not np.array_equal(row.positive_columns, row.negative_columns)
    closure = dependencies
    assert closure["base_commit"] == BASE_COMMIT
    assert closure["source_hashes"]["src/recursive_horizons/nsc_direct_vacuum_source.py"] == sha256(
        (PIN / "src/recursive_horizons/nsc_direct_vacuum_source.py").read_bytes()).hexdigest()
    assert closure["successor_hashes"]["src/nsc_threshold18_source.py"] == sha256(
        (ROOT / "src/nsc_threshold18_source.py").read_bytes()).hexdigest()
    family_payload = "results/development/artifacts/nsc-ks-retained-response-sum/family-14-+1.npz"
    assert family_payload in closure["archive_digests"]
    assert closure["baseline_family_source_hashes"]
    published = json.loads(PUBLISHED.read_text())
    covered = {(item["panel"], item["row"], item["energy_sign"])
               for item in published["source_coverage"]["rows"]}
    assert len(covered) == published["source_coverage"]["covered_signed_rows"] == PUBLISHED_SIGNED_ROWS
    signed = {(row.panel, row.row, sign) for row in catalogue for sign in (1, -1)}
    assert len(signed) == 18 and signed.isdisjoint(covered)
    shared = set(closure["archive_digests"]) & set(published["input_hashes"])
    assert shared
    assert all(closure["archive_digests"][name] == published["input_hashes"][name] for name in shared)
    assert closure["myrsa_used"] if False else closure["base_repository"].startswith("https://github.com/")


def test_endpoint_law_is_nx_minus_ny_minus_nz_and_occupation_is_added_once():
    # Q has nonzero nx, ny, and nz, so (nx,-ny,-nz) is distinct from n and from -n.
    columns = np.array([[1 + 0.25j, 0j, 0j], [0.5 - 0.125j, 0j, 0j]], complex)
    covariance = np.eye(3, dtype=complex)
    endpoint = (arb("0.2"), arb("-0.3"), arb("0.4"))
    validation = {"endpoint": endpoint, "bits": 192, "bloch_error": arb(0)}
    flipped = dict(validation)
    flipped["endpoint"] = (endpoint[0], -endpoint[1], -endpoint[2])
    negated = dict(validation)
    negated["endpoint"] = tuple(-component for component in endpoint)
    complemented = archived_vacuum_distance(columns, covariance, validation, complement=True)
    manual = archived_vacuum_distance(columns, covariance, flipped, complement=False)
    full = archived_vacuum_distance(columns, covariance, negated, complement=False)
    plain = archived_vacuum_distance(columns, covariance, validation, complement=False)
    assert exact_upper(complemented["lower"]) == exact_upper(manual["lower"])
    assert exact_upper(complemented["upper"]) == exact_upper(manual["upper"])
    assert exact_upper(complemented["upper"]) != exact_upper(full["upper"])
    assert exact_upper(complemented["upper"]) != exact_upper(plain["upper"])
    archived = {"lower": arb("0.1"), "upper": arb("0.2")}
    occupation = exact_upper(arb("0.05"))
    compared = owner._comparison(1, 1, (1.6, 1.6, 1.6), "pos", "prep", archived, occupation)
    assert compared["signed_total_lower"] == compared["archived_distance_lower"]
    added = (restored_upper(exact_upper(arb("0.2"))) + restored_upper(occupation)).upper()
    assert exact_upper(added) == compared["signed_total_upper"]
    assert restored_upper(compared["signed_total_upper"]) > restored_upper(compared["archived_distance_upper"])
    assert compared["weight_applied_to_error"] is False
    with pytest.raises(ValueError, match="missing bound"):
        owner._comparison(1, 1, (1.6, 1.6, 1.6), "pos", "prep",
                          {"lower": None, "upper": arb("0.2")}, occupation)


def test_absent_check_keeps_missing_rows_null(context, tmp_path):
    catalogue = context[1]
    absent = tmp_path / "absent-check"
    checked = cover(absent, mode="check", lab_root=PIN)
    assert not absent.exists()
    assert checked["status"].startswith("OPEN")
    assert "ENCLOSED" not in checked["status"]
    assert checked["completed_positive_rows"] == 0
    assert checked["completed_signed_rows"] == 0
    assert checked["missing_positive_rows"] == 9
    assert checked["new_rows_captured"] == 0
    assert checked["maximum_completed_signed_total_upper"] is None
    assert checked["physical_upstream_budget_component"] is None
    assert checked["physical_local_gate"] == "OPEN"
    assert checked["finite_window_closes_gate"] is False
    assert checked["fixed_direct_window_claims_these_rows"] is False
    assert checked["subgap108_batch_included"] is False
    assert checked["negative_endpoint_adjudication"] == ADJUDICATION
    assert checked["negative_endpoint_adjudication_flipped"] is False
    assert checked["weight_applied_to_error"] is False
    assert [item["row"] for item in checked["missing"]] == list(range(32, 41))
    assert all(item["signed_total_upper"] is None and item["archived_distance_upper"] is None
               for item in checked["missing"])
    assert checked["closure"]["base_commit"] == BASE_COMMIT
    assert checked["closure"]["stale_nuc_checkout_used"] is False
    assert checked["schema"] == SCHEMA
    upper = exact_upper(arb(2)**-10)
    completed = []
    for spec in catalogue:
        positive = {"energy_sign": 1, "angular_sign": 1, "signed_total_upper": upper,
                    "weight_applied_to_error": False}
        negative = {"energy_sign": -1, "angular_sign": -1, "signed_total_upper": upper,
                    "weight_applied_to_error": False}
        completed.append({
            "panel": spec.panel, "row": spec.row,
            "energy_fiber_hex": [float(spec.energy).hex()] * 3,
            "negative_energy_fiber_hex": [float(spec.negative_energy).hex()] * 3,
            "payload": {"path": f"rows/{checkpoint_name(spec.panel, spec.row)}/witness.npz"},
            "cells": 1, "nfev": 1, "equation": EQUATION,
            "vacuum_bloch_error_upper": upper,
            "occupation": {"operator_distance_upper": upper},
            "signed_comparisons": [positive, negative],
            "weight_hex": float(spec.weight).hex(),
            "column_weight_hex": float(spec.column_weight).hex(),
            "negative_density_endpoint": "(nx,-ny,-nz)",
        })
    common = dict(mode="check", row_budget=None, cpu_budget=None, cpu_seconds=0.0,
                  capture_cpu_seconds=0.0, new_rows_captured=0)
    finished = aggregate(catalogue, completed, context[2], **common)
    assert finished["coverage_complete"] is True
    assert finished["missing"] == []
    assert finished["status"].startswith("OPEN")
    assert finished["completed_signed_rows"] == 18
    assert finished["physical_local_gate"] == "OPEN"
    assert finished["physical_upstream_budget_component"] is None
    assert finished["maximum_completed_signed_total_upper"] == upper
    assert finished["maximum_completed_signed_total_upper"] != exact_upper(arb(0))
    broken = dict(completed[-1])
    broken["signed_comparisons"] = [broken["signed_comparisons"][0]]
    with pytest.raises(ValueError, match="missing bound"):
        aggregate(catalogue, completed[:-1] + [broken], context[2], **common)


def test_budgets_frozen_targets_and_checkpoint_mutations(context, tmp_path):
    catalogue, dependencies = context[1], context[2]
    spec = catalogue[0]
    with pytest.raises(ValueError, match="frozen source campaign"):
        cover(tmp_path / "results" / "development" / "nsc-direct-source-window-v1",
              mode="check", lab_root=PIN)
    with pytest.raises(ValueError, match="frozen source campaign"):
        cover(tmp_path / "results" / "development" / "nsc-middle-source-coverage-v1",
              mode="resume", row_budget=0, cpu_budget=1, lab_root=PIN)
    with pytest.raises(ValueError, match="parent snapshot"):
        cover(PARENT / "results" / "development" / "not-a-campaign", mode="check", lab_root=PIN)
    with pytest.raises(ValueError, match="pinned source"):
        cover(PIN / "results" / "development" / "not-a-campaign", mode="check", lab_root=PIN)
    with pytest.raises(ValueError, match="row budget"):
        cover(tmp_path / "bad-budget", mode="resume", row_budget=-1, cpu_budget=1, lab_root=PIN)
    with pytest.raises(ValueError, match="CPU budget"):
        cover(tmp_path / "bad-cpu", mode="resume", row_budget=1, cpu_budget=0, lab_root=PIN)
    with pytest.raises(ValueError, match="solve budget"):
        cover(tmp_path / "bad-check", mode="check", row_budget=1, cpu_budget=1, lab_root=PIN)

    paused = cover(tmp_path / "paused", mode="resume", row_budget=0, cpu_budget=1e-9, lab_root=PIN)
    assert paused["new_rows_captured"] == 0
    assert paused["missing_positive_rows"] == 9
    assert paused["maximum_completed_signed_total_upper"] is None
    assert not (tmp_path / "paused" / "rows").exists()
    starved = cover(tmp_path / "starved", mode="resume", row_budget=9, cpu_budget=1e-6, lab_root=PIN)
    assert starved["new_rows_captured"] == 0
    assert starved["completed_signed_rows"] == 0
    assert starved["maximum_completed_signed_total_upper"] is None
    assert all(item["signed_total_upper"] is None for item in starved["missing"])
    assert not (tmp_path / "starved" / "rows").exists()

    driver = runpy.run_path(str(ROOT / "scripts/derive_nsc_threshold18_source.py"))
    assert driver["DEFAULT_OUTPUT"] == "results/development/nsc-threshold18-actual18-v1"
    with pytest.raises(SystemExit) as caught:
        driver["main"](["--resume", "--output", str(tmp_path / "cli-resume")])
    assert caught.value.code == 2
    with pytest.raises(SystemExit) as caught:
        driver["main"](["--check", "--row-budget", "1", "--cpu-budget", "1",
                        "--output", str(tmp_path / "cli-check")])
    assert caught.value.code == 2

    def write_row(destination, record, payload):
        directory = destination / "rows" / checkpoint_name(spec.panel, spec.row)
        directory.mkdir(parents=True)
        (directory / "record.json").write_bytes(encode_record(record))
        (directory / "witness.npz").write_bytes(payload)

    def expect(destination, match):
        with pytest.raises(ValueError, match=match):
            cover(destination, mode="check", lab_root=PIN)
        with pytest.raises(ValueError, match=match):
            cover(destination, mode="resume", row_budget=1, cpu_budget=30, lab_root=PIN)

    base = binding_record(spec, dependencies)
    payload = b"not-a-witness"
    record = dict(base)
    record["rtol_hex"] = float(1e-8).hex()
    changed = tmp_path / "changed-settings"
    write_row(changed, record, payload)
    expect(changed, "settings")

    missing_bound = dict(base)
    del missing_bound["positive_preparation_digest"]
    dropped = tmp_path / "missing-bound"
    write_row(dropped, missing_bound, payload)
    expect(dropped, "missing bound")

    changed_bound = dict(base)
    digests = dict(changed_bound["archive_digests"])
    key = next(iter(digests))
    digest = digests[key]
    digests[key] = digest[:-1] + ("0" if digest[-1] != "0" else "1")
    changed_bound["archive_digests"] = digests
    bound_dir = tmp_path / "changed-bound"
    write_row(bound_dir, changed_bound, payload)
    expect(bound_dir, "changed bound")

    arrays = {
        "trace": np.zeros((1, 26), float),
        "positive_columns": spec.positive_columns.copy(),
        "negative_columns": spec.negative_columns.copy(),
        "positive_covariance": spec.positive_covariance.copy(),
        "negative_covariance": spec.negative_covariance.copy(),
    }
    arrays["positive_columns"][0, 0] += 0.01
    tampered = deterministic_npz_bytes(arrays)
    source_record = dict(base)
    source_record["payload"] = {
        "path": f"rows/{checkpoint_name(spec.panel, spec.row)}/witness.npz",
        "sha256": sha256(tampered).hexdigest(), "bytes": len(tampered),
    }
    source_dir = tmp_path / "changed-source"
    write_row(source_dir, source_record, tampered)
    expect(source_dir, "source")

    missing_payload = tmp_path / "missing-payload" / "rows" / checkpoint_name(spec.panel, spec.row)
    missing_payload.mkdir(parents=True)
    (missing_payload / "record.json").write_bytes(encode_record(base))
    expect(tmp_path / "missing-payload", "payload")

    corrupt = tmp_path / "corrupt" / "rows" / checkpoint_name(spec.panel, spec.row)
    corrupt.mkdir(parents=True)
    (corrupt / "record.json").write_bytes(b"{\n")
    (corrupt / "witness.npz").write_bytes(payload)
    expect(tmp_path / "corrupt", "corrupt")

    extra = tmp_path / "extra"
    write_row(extra, base, payload)
    (extra / "rows" / checkpoint_name(spec.panel, spec.row) / "note.txt").write_bytes(b"x")
    expect(extra, "corrupt")

    duplicate = tmp_path / "duplicate"
    write_row(duplicate, base, payload)
    other = duplicate / "rows" / checkpoint_name(spec.panel, spec.row + 1)
    other.mkdir()
    twin = dict(base)
    (other / "record.json").write_bytes(encode_record(twin))
    (other / "witness.npz").write_bytes(payload)
    expect(duplicate, "duplicate")

    witness = deterministic_npz_bytes({
        "trace": np.zeros((1, 26), float),
        "positive_columns": spec.positive_columns,
        "negative_columns": spec.negative_columns,
        "positive_covariance": spec.positive_covariance,
        "negative_covariance": spec.negative_covariance,
    })
    published = dict(base)
    published["payload"] = {
        "path": f"rows/{checkpoint_name(spec.panel, spec.row)}/witness.npz",
        "sha256": sha256(witness).hexdigest(), "bytes": len(witness),
    }
    exclusive = tmp_path / "exclusive"
    owner.publish_checkpoint(exclusive, published, witness)
    with pytest.raises(PrepublicationError, match="destination exists"):
        owner.publish_checkpoint(exclusive, published, witness)
    assert (exclusive / "rows" / checkpoint_name(spec.panel, spec.row) / "witness.npz").read_bytes() == witness
    with np.load(BytesIO(witness), allow_pickle=False) as loaded:
        assert np.array_equal(loaded["negative_columns"], s3_conjugate(loaded["positive_columns"]))
    assert restored_upper(exact_upper(arb("1.5"))) > 0
