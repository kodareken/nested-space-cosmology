"""Tests for the narrow-domain failure.

Narrowing the outer radius was tried as a cheap route to more resolution
near the inner edge. It produced gains two to two and a half times larger and
order-one boundary curvature, and it is unstable: about 29% of all node-steps
reverse. These tests pin the instability, the withdrawal, and the reason the
narrow variant must not be read as a stronger result.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
NARROW_SWEEP = DEVELOPMENT / "nsc-rn-edge-sweep-narrow-v2.json"
NARROW_SWEEP_V1 = DEVELOPMENT / "nsc-rn-edge-sweep-narrow-v1.json"
NARROW_STAB = DEVELOPMENT / "nsc-rn-edge-stability-narrow-v1.json"
WIDE_STAB = DEVELOPMENT / "nsc-rn-edge-stability-v2.json"
NARROW_CURV = DEVELOPMENT / "nsc-rn-curvature-n6-exc080-v1.json"
NARROW_SHELL = DEVELOPMENT / "nsc-rn-inner-packet-shell-n6-exc080-v1.json"
WIDE_SWEEP = DEVELOPMENT / "nsc-rn-edge-sweep-v2.json"
NOTE = LAB / "docs/nsc-rn-edge-saturation.md"
FAILURE_NOTE = LAB / "docs/nsc-rn-narrow-failure.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


class TestTheNarrowRunsAreUnstable:
    def test_the_stability_audit_failed(self):
        assert _load(NARROW_STAB)["all_runs_stable"] is False

    def test_every_narrow_run_failed(self):
        for run in _load(NARROW_STAB)["runs"]:
            assert not run["stable"]

    def test_the_reversals_are_pervasive(self):
        for run in _load(NARROW_STAB)["runs"]:
            assert run["mass_decreases_below_limit"] > 10000

    def test_the_wide_runs_have_none(self):
        """The contrast is the whole point of the comparison."""
        for run in _load(WIDE_STAB)["runs"]:
            assert run["mass_decreases_below_limit"] == 0

    def test_all_the_fields_were_still_finite(self):
        """So finiteness alone would not have caught it."""
        for run in _load(NARROW_STAB)["runs"]:
            assert run["all_fields_finite"]
            assert run["lapse_always_positive"]

    def test_the_audit_would_have_caught_it(self):
        """Monotonicity, not finiteness, is the discriminating check."""
        assert "finiteness checks pass" in FAILURE_NOTE.read_text(encoding="utf-8")


class TestTheNarrowNumbersAreWithdrawn:
    def test_the_note_withdraws_them_by_name(self):
        note = FAILURE_NOTE.read_text(encoding="utf-8")
        assert "**Withdrawn**" in note
        assert "nsc-rn-edge-sweep-narrow-v1.json" in note

    def test_it_says_they_are_what_an_unstable_march_produces(self):
        assert "what an unstable march produces" in FAILURE_NOTE.read_text(encoding="utf-8")

    def test_the_records_are_retained_not_deleted(self):
        for path in (NARROW_SWEEP, NARROW_STAB, NARROW_CURV, NARROW_SHELL):
            assert path.is_file()

    def test_the_wide_saturation_is_unchanged(self):
        block = _load(WIDE_SWEEP)["sweep"]
        assert block["amplitude_gain"][-1] == pytest.approx(1.5849, abs=1e-3)


class TestTheCurvatureClaimIsNotTakenAtFaceValue:
    def test_the_narrow_curvature_verdict_says_growth(self):
        record = _load(NARROW_CURV)
        assert record["verdict"] == "edge_curvature_grows_beyond_any_settled_scale"
        assert record["R4_edge_max_over_initial"] > 100.0

    def test_the_narrow_curvature_is_order_one(self):
        assert _load(NARROW_CURV)["R4_edge_max_abs"] > 0.1

    def test_the_edge_layer_is_real_not_a_single_cell(self):
        """So it is not dismissible as a spike, only as unstable."""
        record = _load(NARROW_CURV)
        assert record["R4_edge_over_neighbour"][-1] > 5.0


class TestWhyTheWideDomainIsLoadBearing:
    def test_the_note_says_the_far_field_absorbs_the_interaction(self):
        note = FAILURE_NOTE.read_text(encoding="utf-8")
        assert "carries the room the march needs" in note
        assert "load-bearing in numerics" in note

    def test_it_connects_to_the_earlier_far_field_finding(self):
        assert "negligible in mass and apparently load-bearing" in \
            FAILURE_NOTE.read_text(encoding="utf-8")


class TestScope:
    def test_it_says_it_is_not_the_universe(self):
        assert "does not measure the\nactual universe" in FAILURE_NOTE.read_text(encoding="utf-8")

    def test_it_says_the_outer_radius_matters(self):
        assert "the outer radius is now known to\nmatter" in \
            FAILURE_NOTE.read_text(encoding="utf-8")

    def test_it_names_the_open_question_as_a_resources_question(self):
        note = FAILURE_NOTE.read_text(encoding="utf-8")
        assert "resources question" in note
        assert "treated as one rather than worked around again" in note

    def test_the_question_stays_open(self):
        assert "Still open" in FAILURE_NOTE.read_text(encoding="utf-8")


class TestTheRouteIsClosed:
    def test_every_narrow_run_used_the_outer_override(self):
        for run in _load(NARROW_SWEEP)["runs"]:
            assert run["r_outer"] == pytest.approx(6.0)

    def test_the_narrow_runs_really_are_better_resolved(self):
        """The plan worked; only the march failed."""
        for run in _load(NARROW_SWEEP)["runs"]:
            assert run["gap_in_cells"] > 7.0

    def test_the_first_narrow_sweep_had_the_wrong_cells(self):
        first = _load(NARROW_SWEEP_V1)
        assert min(first["sweep"]["gap_in_cells"]) < 2.0
        assert first["sweep"]["gap_in_cells"] ==             _load(WIDE_SWEEP)["sweep"]["gap_in_cells"][:3]

    def test_the_cells_are_now_read_from_each_run(self):
        wide = _load(WIDE_SWEEP)["sweep"]["gap_in_cells"]
        narrow = _load(NARROW_SWEEP)["sweep"]["gap_in_cells"]
        assert narrow != wide

    def test_the_wide_runs_used_the_default_outer(self):
        for run in _load(WIDE_SWEEP)["runs"]:
            assert run["r_outer"] == pytest.approx(32.0)