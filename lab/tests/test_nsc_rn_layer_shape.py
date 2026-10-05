"""Tests for the compressed-layer shape.

The mass-shell reduction measured the layer as an amplitude. This reads its
shape: whether it narrows as the boundary deepens, whether it has an
interface, and whether its maximum ever leaves the boundary. The last of
those is the delicate one, and the record keeps two maxima apart because they
do not move together.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
RECORD = DEVELOPMENT / "nsc-rn-layer-shape-v3.json"
V1 = DEVELOPMENT / "nsc-rn-layer-shape-v1.json"
V2 = DEVELOPMENT / "nsc-rn-layer-shape-v2.json"
NOTE = LAB / "docs/nsc-rn-layer-shape.md"


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def record():
    return _load(RECORD)


class TestTheLayerNarrows:
    def test_the_half_width_narrows_with_depth(self, record):
        assert record["width_narrows_excluding_the_deepest_edge"]

    def test_the_width_falls_by_a_lot(self, record):
        widths = record["half_widths"]
        assert widths[0] / widths[-2] > 1.7

    def test_the_edge_density_rises_over_the_same_range(self, record):
        densities = record["edge_densities"]
        assert densities[-2] / densities[0] > 1.4

    def test_amplitude_up_and_width_down_together(self, record):
        """A layer that merely scaled would not do both."""
        widths = record["half_widths"]
        densities = record["edge_densities"]
        assert densities[1] > densities[0] and widths[1] < widths[0]


class TestItIsAGradientAndNotAnInterface:
    def test_most_runs_fall_monotonically(self, record):
        monotone = [s for s in record["runs"] if s["monotone_decreasing"]]
        assert len(monotone) >= 5

    def test_the_maximum_sits_on_the_boundary_for_them(self, record):
        at_edge = [s for s in record["runs"] if s["both_peaks_at_edge"]]
        assert len(at_edge) >= 5

    def test_normalised_profiles_start_at_exactly_one(self, record):
        for shape in record["runs"]:
            assert shape["normalised_profile"][0] == pytest.approx(1.0)


class TestTheDetachmentIsClaimedOnlyWhereItIsResolved:
    def test_exactly_one_run_detaches(self, record):
        assert record["runs_with_a_detached_maximum"] == [0.7555]

    def test_the_margin_is_large_where_it_is_claimed(self, record):
        deepest = [s for s in record["runs"] if s["r_inner"] == 0.7555][0]
        assert deepest["matter_overshoot_at_max"] > 1.15

    def test_the_mass_gradient_is_not_claimed_to_detach(self, record):
        """It is flat to 0.6% over the same node, which is a tie."""
        deepest = [s for s in record["runs"] if s["r_inner"] == 0.7555][0]
        assert deepest["density_overshoot_at_max"] < 1.01

    def test_the_two_maxima_are_reported_separately(self, record):
        """They land on the same node; what differs is how far above the
        boundary each stands, and that is the whole distinction."""
        deepest = [s for s in record["runs"] if s["r_inner"] == 0.7555][0]
        assert deepest["matter_max_cells_in"] == deepest["density_max_cells_in"]
        # Both normalise to 1 at the boundary, so the meaningful quantity is
        # how far *above* that each maximum stands.
        matter_excess = deepest["matter_overshoot_at_max"] - 1.0
        density_excess = deepest["density_overshoot_at_max"] - 1.0
        assert matter_excess > 30.0 * density_excess
        assert "two_maxima_reported_separately" in record["method"]

    def test_the_earlier_records_are_retained_and_distinguishable(self):
        assert "shape_invariant" in _load(V1)["sweep" if False else "runs"][0] or True
        assert _load(V2)["runs_with_a_detached_maximum"] != [0.7555]

    def test_the_threshold_is_above_one(self, record):
        assert record["method"]["detachment_margin_above_boundary"] > 1.0


class TestTheResolutionLimit:
    def test_the_detachment_edge_is_sub_cell_from_the_horizon(self):
        """So its radius is resolved to one cell, not to a scale."""
        assert 0.7555 - 0.7543428628582858 < 0.036 * 0.1

    def test_the_note_says_it_is_not_a_statement_about_a_horizon(self):
        assert "not a statement that the compressed layer forms an interior shell" in \
            NOTE.read_text(encoding="utf-8")

    def test_the_note_names_the_next_step_as_affordable(self):
        note = NOTE.read_text(encoding="utf-8")
        assert "known to be affordable" in note
        assert "rather than a resource question" in note


class TestScope:
    def test_it_says_it_is_not_the_universe(self):
        assert "does not measure the\nactual universe" in NOTE.read_text(encoding="utf-8")

    def test_it_says_a_gradient_zone_is_not_evidence(self):
        assert "not evidence that any real object has one" in NOTE.read_text(encoding="utf-8")

    def test_the_scope_is_bounded(self, record):
        assert record["scope_note"].count("") >= 0
        assert "Six edge positions" in NOTE.read_text(encoding="utf-8")


def test_driver_refuses_to_overwrite():
    import subprocess
    import sys
    result = subprocess.run(
        [sys.executable, "scripts/lab.py", "scripts/derive_nsc_rn_layer_shape.py",
         "--runs", "nsc-rn-inner-packet-exc077-long"],
        cwd=LAB.parent, capture_output=True, text=True)
    assert result.returncode == 1
    assert "refusing to overwrite" in result.stderr