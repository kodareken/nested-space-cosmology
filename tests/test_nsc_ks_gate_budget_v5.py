"""Gate Budget v5: PASS-only fill, frozen allocations and OPEN nulls."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from recursive_horizons.nsc_ks_continuous_constraint_enclosure import (
    PROFILE_IDENTITY,
    STATE_LAW,
)
from recursive_horizons.nsc_ks_assembly_arithmetic import (
    DECLARED_POSITIVE_FAMILIES,
    DECLARED_SIGNED_FAMILIES,
)
from recursive_horizons.nsc_ks_gate_budget_v5 import (
    ALLOCATION,
    ALLOCATION_TOTAL,
    ERROR_COMPONENTS,
    RESIDUAL_RESERVE,
    SCHEMA,
    admit_budget,
    compose_open_budget,
    fill_component,
    validate_open_budget,
)
from recursive_horizons.nsc_local_gate_certificate_v2 import ERROR_COMPONENTS as CERT_COMPONENTS
import derive_nsc_ks_gate_budget_v2 as B


def digest(raw):
    return sha256(raw).hexdigest()


def pass_enclosure(bound=(1e-13, 2e-13)):
    path = "results/development/nsc-ks-gate-budget-v4.json"
    return {
        "status": "PASS: directed enclosure",
        "bound": list(bound),
        "certificate_use": True,
        "numerical_indicators_used_as_bounds": False,
        "profile_identity": PROFILE_IDENTITY,
        "state_law": STATE_LAW,
        "node_count": 257,
        "evidence": {"path": path, "sha256": digest((ROOT / path).read_bytes())},
        "scope": "unit-test PASS enclosure",
    }


def test_nine_named_components_match_the_certificate_and_v2_allocation():
    assert tuple(ERROR_COMPONENTS) == CERT_COMPONENTS
    assert B.ALLOCATION == ALLOCATION
    assert sum(ALLOCATION.values()) == pytest.approx(ALLOCATION_TOTAL)
    assert RESIDUAL_RESERVE == 1e-11
    record = compose_open_budget()
    assert set(record["components"]) == set(ERROR_COMPONENTS)
    assert record["missing_components"] == list(ERROR_COMPONENTS)
    assert record["full_error_sum"] is None
    assert record["schema"] == SCHEMA
    assert record["v4_bytes_preserved"]
    assert record["declared_positive_families"] == DECLARED_POSITIVE_FAMILIES
    assert record["declared_signed_families"] == DECLARED_SIGNED_FAMILIES
    assert record["covered_positive_families"] is None
    assert record["covered_signed_families"] is None
    assert record["source_inventory_authenticated"] is False
    validate_open_budget(record)


def test_fill_rejects_enclosed_open_samples_and_profile_mismatch():
    with pytest.raises(ValueError, match="PASS"):
        fill_component("between_node", {
            "status": "ENCLOSED: predecessor", "bound": [1e-13, 1e-13],
            "certificate_use": True})
    with pytest.raises(ValueError, match="PASS"):
        fill_component("arithmetic", {
            "status": "OPEN: missing", "bound": [1e-13, 1e-13],
            "certificate_use": True})
    with pytest.raises(ValueError, match="diagnostic"):
        fill_component("between_node", {
            "status": "PASS: control", "bound": [1e-13, 1e-13],
            "certificate_use": False})
    with pytest.raises(ValueError, match="samples"):
        fill_component("between_node", {
            "status": "PASS: samples", "bound": [1e-13, 1e-13],
            "certificate_use": True, "samples_used_as_bound": True})
    with pytest.raises(ValueError, match="profile"):
        fill_component("arithmetic", {
            **pass_enclosure(), "profile_identity": digest(b"other profile")})
    with pytest.raises(ValueError, match="state law"):
        fill_component("arithmetic", {**pass_enclosure(), "state_law": "frozen C0"})
    with pytest.raises(ValueError, match="profile identity required"):
        fill_component("arithmetic", {
            **pass_enclosure(), "profile_identity": None})
    with pytest.raises(ValueError, match="state law required"):
        fill_component("arithmetic", {**pass_enclosure(), "state_law": None})
    with pytest.raises(ValueError, match="node_count required"):
        fill_component("between_node", {
            k: v for k, v in pass_enclosure().items() if k != "node_count"})
    row = fill_component("between_node", pass_enclosure())
    assert row["bound"] == [1e-13, 2e-13]
    assert row["status"].startswith("PASS")
    assert row["node_count"] == 257


def test_predecessor_enclosed_rows_do_not_silently_fill_v5():
    predecessor = json.loads(
        (ROOT / "results/development/nsc-ks-gate-budget-v4.json").read_text())
    record = compose_open_budget(
        predecessor=predecessor, predecessor_sha256=digest(b"v4"))
    assert record["predecessor_enclosed_components"] == [
        "energy_interpolation", "covered_regions", "phase_value"]
    assert record["predecessor_bounds_are_not_v5_pass_fills"]
    for name in ERROR_COMPONENTS:
        assert record["components"][name]["bound"] is None
        assert record["components"][name]["allocation_per_component"] == ALLOCATION[name]
    assert record["filled_from_pass"] == []
    with pytest.raises(ValueError, match="missing"):
        admit_budget(record)


def test_fill_rejects_n129_phase_and_fabricated_digest():
    predecessor = json.loads(
        (ROOT / "results/development/nsc-ks-gate-budget-v4.json").read_text())
    phase = dict(predecessor["components"]["phase_value"])
    phase["status"] = "PASS: reused n129"
    phase["certificate_use"] = True
    with pytest.raises(ValueError, match="profile identity required"):
        fill_component("phase_value", phase)
    phase["profile_identity"] = PROFILE_IDENTITY
    phase["state_law"] = STATE_LAW
    phase["node_count"] = 129
    with pytest.raises(ValueError, match="129-node"):
        fill_component("phase_value", phase)
    phase["node_count"] = 257
    with pytest.raises(ValueError, match="129-node"):
        fill_component("phase_value", phase)
    fake = {
        **pass_enclosure(),
        "evidence": {
            "path": "results/development/nsc-ks-current-phase-accuracy-v3-n129.json",
            "sha256": "ab" * 32,
        },
    }
    with pytest.raises(ValueError, match="does not match file bytes"):
        fill_component("phase_value", fake)
    missing = dict(pass_enclosure())
    missing.pop("evidence")
    with pytest.raises(ValueError, match="evidence"):
        fill_component("between_node", missing)


def test_admission_requires_all_nine_finite_total_and_reserve():
    enclosures = {name: pass_enclosure((1e-13, 1e-13)) for name in ERROR_COMPONENTS}
    record = compose_open_budget(pass_enclosures=enclosures)
    record["all_declared_source_paths_closed"] = True
    record["status"] = "PASS: nine components enclosed"
    record["missing_components"] = []
    record["covered_positive_families"] = DECLARED_POSITIVE_FAMILIES
    record["covered_signed_families"] = DECLARED_SIGNED_FAMILIES
    record["source_inventory_authenticated"] = True
    with pytest.raises(ValueError, match="authenticated 60/120"):
        admit_budget(record)
    moved = deepcopy(record)
    moved["components"]["between_node"]["allocation_per_component"] = 2e-12
    moved["allocation_changed"] = True
    with pytest.raises(ValueError, match="moved"):
        admit_budget(moved)
    huge = deepcopy(record)
    huge["components"]["field_space_time"]["bound"] = [3e-11, 0.0]
    huge["full_error_sum"] = [4e-11, 1e-12]
    with pytest.raises(ValueError, match="2e-11"):
        admit_budget(huge)
    existence = deepcopy(record)
    existence["physical_EXISTENCE_certificate"] = True
    with pytest.raises(ValueError, match="EXISTENCE"):
        admit_budget(existence)


def test_v4_bytes_are_preserved_and_open_record_stays_null():
    v4_path = ROOT / "results/development/nsc-ks-gate-budget-v4.json"
    before = v4_path.read_bytes()
    predecessor = json.loads(before)
    record = compose_open_budget(predecessor=predecessor)
    validate_open_budget(record)
    assert v4_path.read_bytes() == before
    assert predecessor["schema"] == "NSC-KS-GATE-ERROR-BUDGET-v4"
    forged = deepcopy(record)
    forged["status"] = "PASS: nulls"
    with pytest.raises(ValueError, match="OPEN"):
        validate_open_budget(forged)
