"""Research prerequisites stay distinct from physical connections and evidence."""
from copy import deepcopy
from pathlib import Path

import pytest

from recursive_horizons import nsc_cosmological_atlas as model


def programme_atlas():
    atlas = model.build_atlas()
    atlas["programme"] = model.programme_record()
    return atlas


def annotation(kind="preparation"):
    row = {
        "id": "open-connection", "kind": kind, "status": "unresolved",
        "domain": "Finite early wall-to-centre observer interval in the fixed reference geometry",
        "assumptions": ["fixed geometry and field law"],
        "missing": ["physical exterior preparation or transition mechanism"],
        "owner_refs": ["docs/nsc-atlas-field-transfer.md"], "evidence_refs": [],
        "region_refs": ["child_a"], "link_refs": [],
        "prerequisite_ids": ["parent-preparation"], "unlock_ids": ["observer-response"],
    }
    if kind == "preparation":
        row.update(origin="prescribed_at_boundary", remote_preparation_established=False)
    elif kind == "transition_inheritance":
        row.update(map_ref=None, transferred_quantities=[], physical_inheritance_established=False)
    else:
        row.update(observer_region="child_a", observables=["local scalar response"],
                   alternatives=["exterior A", "exterior B"], equivalence_established=False,
                   discriminator="A second independent observable after a justified preparation")
    return row


def test_original_atlas_bytes_and_causal_results_are_unchanged():
    atlas = model.build_atlas()
    assert model.fingerprint(atlas) == "ceb8d7c8cadd0e6bb97ff78318de4eb47f3c04e77a659119dd88cddeaef648bd"
    original = deepcopy(atlas)
    atlas["programme"] = model.programme_record()
    atlas["scenarios"][0]["connection_records"] = [annotation()]
    assert model.validate_atlas(atlas)
    assert atlas["scenarios"][0]["links"] == original["scenarios"][0]["links"]
    assert model.earliest_arrival(atlas["scenarios"][0], "child_a", "child_b", 0, clock="tau") is None
    assert "programme" not in model.build_atlas()


def test_eight_questions_have_existing_owners_evidence_and_reading_limits():
    record = model.programme_record()
    assert {q["id"] for q in record["questions"]} == set(model.PROGRAMME_THEMES)
    root = Path(__file__).resolve().parents[1]
    for row in record["questions"]:
        assert row["reading_boundary"] and row["missing"]
        for ref in row["owner_refs"] + row["evidence_refs"]:
            assert (root / ref.split("#", 1)[0]).is_file(), ref
    assert record["new_compute_admitted"] is False
    assert record["physical_recurrence_established"] is False
    record["questions"][0]["missing"].clear()
    assert model.programme_record()["questions"][0]["missing"]


def test_research_dag_distinguishes_reused_equations_from_missing_mechanism():
    record = model.programme_record()
    nodes = {row["id"]: row for row in record["calculation_graph"]}
    assert nodes["wall-transmission"]["gap_kind"] == "reusable_result"
    assert nodes["parent-preparation"]["gap_kind"] == "known_equations_missing_data"
    assert nodes["inheritance-map"]["gap_kind"] == "missing_physical_mechanism"
    assert set(nodes["wall-trace"]["prerequisite_ids"]) == {"parent-preparation", "wall-transmission"}
    assert nodes["observer-response"]["unlock_ids"] == ["observer-inference", "clock-entropy"]
    nodes["wall-trace"]["prerequisite_ids"].append("observer-response")
    nodes["observer-response"]["unlock_ids"].append("wall-trace")
    with pytest.raises(ValueError, match="cycle"):
        model.validate_programme(record)


def test_unknown_dependencies_and_inconsistent_unlocks_are_rejected():
    record = model.programme_record()
    record["questions"][0]["unlock_ids"].append("unrecorded-computation")
    with pytest.raises(ValueError, match="unknown research"):
        model.validate_programme(record)
    record = model.programme_record()
    record["calculation_graph"][0]["unlock_ids"] = []
    with pytest.raises(ValueError, match="agree with prerequisite"):
        model.validate_programme(record)


@pytest.mark.parametrize("kind,flag", [
    ("preparation", "remote_preparation_established"),
    ("transition_inheritance", "physical_inheritance_established"),
    ("observer_equivalence", "equivalence_established"),
])
def test_missing_or_unsupported_claims_cannot_be_promoted(kind, flag):
    atlas = programme_atlas()
    row = annotation(kind)
    atlas["scenarios"][0]["connection_records"] = [row]
    assert model.validate_atlas(atlas)
    del row[flag]
    with pytest.raises(ValueError, match="must be explicit"):
        model.validate_atlas(atlas)
    row[flag] = True
    with pytest.raises(ValueError, match="requires evidence and independent"):
        model.validate_atlas(atlas)


def test_prescribed_input_cannot_become_remote_preparation_even_with_citations():
    atlas = programme_atlas()
    row = annotation()
    row.update(status="established_in_domain", remote_preparation_established=True,
               evidence_refs=["results/development/nsc-cosmological-atlas-v6/field-transfer.json"],
               independent_check_refs=["tests/test_nsc_atlas_field_transfer.py"])
    atlas["scenarios"][0]["connection_records"] = [row]
    with pytest.raises(ValueError, match="prescribed boundary data"):
        model.validate_atlas(atlas)


def test_unknown_physical_and_research_references_are_rejected():
    for field in ("region_refs", "link_refs", "prerequisite_ids", "unlock_ids"):
        atlas = programme_atlas()
        row = annotation()
        row[field] = ["unknown"]
        atlas["scenarios"][0]["connection_records"] = [row]
        with pytest.raises(ValueError, match="unknown"):
            model.validate_atlas(atlas)
    atlas = programme_atlas()
    atlas["scenarios"][0]["connection_records"] = [None]
    with pytest.raises(ValueError, match="list of records"):
        model.validate_atlas(atlas)


def test_optional_transition_law_cannot_retune_a_fixed_theory_comparison():
    original = programme_atlas()["scenarios"][0]
    original["connection_records"] = [annotation("transition_inheritance")]
    changed = deepcopy(original)
    changed["regions"][1]["state"] = "changed input state"
    assert model.require_fixed_theory(original, changed)
    changed["connection_records"][0]["map_ref"] = "docs/nsc-nested-qualities.md"
    with pytest.raises(ValueError, match="transition law changed"):
        model.require_fixed_theory(original, changed)


def test_incomplete_observer_equivalence_does_not_remove_a_candidate():
    atlas = programme_atlas()
    row = annotation("observer_equivalence")
    atlas["scenarios"][0]["connection_records"] = [row]
    assert model.validate_atlas(atlas)
    assert len(atlas["scenarios"]) == 4
    row["alternatives"] = ["exterior A"]
    with pytest.raises(ValueError, match="at least two"):
        model.validate_atlas(atlas)
