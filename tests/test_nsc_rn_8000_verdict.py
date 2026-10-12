"""Tests for the 8000-point verdict.

The coarse grid put the saturation at x1.58 and the fine grid puts it near
x2.05. These tests pin the correction, the reason the earlier convergence
claim was too strong, and the one structural finding that did not move.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
COARSE_SHAPE = DEVELOPMENT / "nsc-rn-layer-shape-v5.json"
FINE_SHAPE = DEVELOPMENT / "nsc-rn-layer-shape-n8000-v1.json"
CONVERGENCE = DEVELOPMENT / "nsc-rn-resolution-convergence-v3.json"
NOTE = LAB / "docs/nsc-rn-8000-verdict.md"

FINE_SHELLS = {
    1.04: "n8000-def1040", 0.85: "n8000-exc085", 0.80: "n8000d-exc080",
    0.77: "n8000d-exc077", 0.76: "n8000d-exc076", 0.7555: "n8000d-exc07555",
}
COARSE_GAINS = {
    1.04: 1.1972, 0.85: 1.3905, 0.80: 1.4328,
    0.77: 1.5758, 0.76: 1.5825, 0.7555: 1.5849,
}


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def fine_gains():
    return {edge: _load(DEVELOPMENT / f"nsc-rn-inner-packet-shell-{tag}-v1.json")
            ["accumulation_after_arrival"]["amplitude_gain"]
            for edge, tag in FINE_SHELLS.items()}


class TestTheCorrection:
    def test_the_fine_gains_are_all_larger(self, fine_gains):
        for edge, coarse in COARSE_GAINS.items():
            assert fine_gains[edge] > coarse

    def test_the_largest_change_is_at_the_deep_edges(self, fine_gains):
        shallow = fine_gains[1.04] / COARSE_GAINS[1.04]
        deep = fine_gains[0.77] / COARSE_GAINS[0.77]
        assert deep > shallow

    def test_the_saturation_is_near_two_not_one_point_five(self, fine_gains):
        peak = max(fine_gains.values())
        assert 1.95 < peak < 2.10

    def test_the_turn_over_moves_to_zero_point_seven_six(self, fine_gains):
        edges = sorted(fine_gains, reverse=True)
        best = max(edges, key=lambda e: fine_gains[e])
        assert best == 0.77

    def test_the_coarse_numbers_are_retained_as_a_record(self):
        coarse = _load(COARSE_SHAPE)
        assert coarse["runs"][-1]["r_inner"] == 0.7555
        assert "1.58" in NOTE.read_text(encoding="utf-8")


class TestWhyTheTrajectoryConvergedAndTheGainDidNot:
    def test_the_trajectory_convergence_was_real(self):
        convergence = _load(CONVERGENCE)
        assert convergence["worst_relative_difference"] < 1e-4

    def test_the_note_says_it_does_not_imply_the_gain(self):
        note = NOTE.read_text(encoding="utf-8")
        flat = note.replace(chr(10), " ")
        assert "did not imply what I claimed" in flat
        assert "is not converged to the same accuracy" in flat

    def test_the_generalisation_is_named_as_the_error(self):
        assert "That was my error" in NOTE.read_text(encoding="utf-8")


class TestWhatSurvived:
    def test_the_shell_sidedness_is_identical_at_both_resolutions(self):
        coarse = _load(COARSE_SHAPE)
        fine = _load(FINE_SHAPE)
        assert coarse["shell_sides"] == fine["shell_sides"] == [2, 2, 2, 2, 2, 3]

    def test_only_the_deepest_edge_has_three_sides(self):
        fine = _load(FINE_SHAPE)
        assert fine["edges_with_a_three_sided_shell"] == [0.7555]

    def test_the_narrowing_is_monotone_at_both_resolutions(self):
        coarse = _load(COARSE_SHAPE)
        fine = _load(FINE_SHAPE)
        assert coarse["width_narrows_excluding_the_deepest_edge"]
        assert fine["width_narrows_excluding_the_deepest_edge"]

    def test_the_narrowing_is_stronger_at_the_fine_resolution(self):
        coarse = _load(COARSE_SHAPE)
        fine = _load(FINE_SHAPE)
        coarse_factor = coarse["half_widths"][0] / coarse["half_widths"][4]
        fine_factor = fine["half_widths"][0] / fine["half_widths"][4]
        assert fine_factor > coarse_factor
        assert fine_factor == pytest.approx(2.87, abs=0.02)


class TestScope:
    def test_it_does_not_claim_the_relaxation_converged(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "not been compared across resolutions" in flat

    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat

    def test_it_does_not_claim_either_resolution_is_converged(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "no reason to believe either is converged" in flat
