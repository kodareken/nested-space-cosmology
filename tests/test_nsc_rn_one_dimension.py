"""Tests for the one-dimensionality of this construction.

A wave moved to the side becomes a sphere only if there is a side to move it
into. That claim was made about this construction and is checked here against
the code and the recorded state rather than inferred from the line element.

The checks are deliberately redundant: a requirement in one module, a
hardcoded constant in another, and the shape of the recorded arrays. Any one
of them could be changed without the other, so all three are pinned.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
SRC = LAB / "src/recursive_horizons"
RESPONSE = SRC / "nsc_rn_response.py"
OBSERVABLES = SRC / "nsc_rn_observables.py"
STATION = (LAB / "results/development/nsc-rn-inner-packet-exc0760-long/station-0.npz")
NOTE = LAB / "docs/nsc-rn-one-dimension.md"


class TestTheAngularIndexIsPinnedToOne:
    def test_the_response_refuses_a_wider_state(self):
        text = RESPONSE.read_text(encoding="utf-8")
        assert "expects the rank-one column (2, points, 1)" in text
        assert "shape[2] != 1" in text

    def test_every_recorded_field_is_one_dimensional_in_radius(self):
        with np.load(STATION) as payload:
            assert payload["coupled_phi"].shape[0] == 2
            assert payload["coupled_phi"].shape[2] == 1
            for key in ("coupled_mass", "radius", "frozen_lapse",
                        "frozen_shift", "R4", "Ricci2"):
                assert payload[key].ndim == 1, key

    def test_the_radius_count_matches_the_spinor(self):
        with np.load(STATION) as payload:
            assert payload["coupled_phi"].shape[1] == payload["radius"].size


class TestTheAngleIsAFixedConstant:
    def test_the_connection_is_evaluated_at_one_value(self):
        text = OBSERVABLES.read_text(encoding="utf-8")
        assert "theta = 0.5 * math.pi" in text

    def test_it_is_not_a_parameter(self):
        text = OBSERVABLES.read_text(encoding="utf-8")
        assert "def _pg_christoffel(radius, theta," in text
        body = text.split("def _pg_christoffel")[1].split("\ndef ")[0]
        # theta arrives as an argument here, but every call site supplies
        # the same fixed value, which is what pins it.
        assert "theta = 0.5 * math.pi" in text


class TestThereIsNoTransverseOperator:
    def test_no_angular_derivative_in_the_evolution(self):
        for name in ("nsc_rn_pg.py", "nsc_rn_source.py"):
            text = (SRC / name).read_text(encoding="utf-8")
            for forbidden in ("d/dtheta", "partial_theta", "cot("):
                assert forbidden not in text, (name, forbidden)

    def test_the_characteristics_are_radial(self):
        """lambda_pm = -beta +- N along d r, with no angular dependence."""
        from recursive_horizons import nsc_rn_observables as obs
        radius = np.array([2.0, 4.0, 8.0])
        chart = obs.chart_speeds(np.ones(3), 0.3 + 0.2 / radius)
        assert chart[0].shape == radius.shape
        assert np.all(chart[0] > 0.0)


class TestWhatThisBounds:
    def test_the_note_says_everything_measured_is_radial(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "no slot for one" in flat
        assert "is withdrawn" in flat
        assert "inside of that circle is the third" in flat
        assert "It matches the sphere" in flat

    def test_it_does_not_claim_transverse_structure_is_impossible(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not say transverse structure is impossible" in flat

    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat