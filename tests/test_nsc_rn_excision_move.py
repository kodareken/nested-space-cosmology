"""Tests for moving the grid's inner edge inboard of the horizon midpoint.

Every earlier inner-packet result was measured at the default excision, which
is the horizon midpoint. Moving it is the first step toward containing the
accumulation, and it changed two headline numbers. These tests pin the
default path's invariance, the override's validity, and the finding itself.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_pg as pg  # noqa: E402

DEVELOPMENT = LAB / "results/development"
DEEP_RUN = DEVELOPMENT / "nsc-rn-inner-packet-exc085"
SHALLOW = DEVELOPMENT / "nsc-rn-inner-packet-shell-v2.json"
DEEP_SHELL = DEVELOPMENT / "nsc-rn-inner-packet-shell-exc085-v1.json"
DEEP_CURV = DEVELOPMENT / "nsc-rn-curvature-exc085-v1.json"


def _geometry(mu=1.04):
    import derive_nsc_rn_parent as driver
    return driver.case_geometry(mu)


class TestTheDefaultPathIsUnchanged:
    def test_the_default_is_still_the_horizon_midpoint(self):
        geometry = _geometry()
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        r_minus = float(model.horizons["r_minus"])
        r_plus = float(model.horizons["r_plus"])
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"])
        assert grid["excision"] == pytest.approx(0.5 * (r_minus + r_plus))
        assert grid["radius"][0] == pytest.approx(0.5 * (r_minus + r_plus))
        assert grid["excision_is_default"]

    def test_two_default_grids_are_identical(self):
        geometry = _geometry()
        a = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"])
        b = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"])
        assert np.array_equal(a["radius"], b["radius"])
        assert np.array_equal(a["derivative"].toarray() if hasattr(a["derivative"], "toarray")
                              else a["derivative"], b["derivative"].toarray()
                              if hasattr(b["derivative"], "toarray") else b["derivative"])

    def test_the_grid_records_which_edge_it_used(self):
        geometry = _geometry()
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"], excision=0.85)
        assert grid["excision"] == pytest.approx(0.85)
        assert not grid["excision_is_default"]


class TestTheOverrideIsValid:
    @pytest.mark.parametrize("edge", [0.80, 0.85, 0.90, 1.00])
    def test_edges_inside_the_horizons_are_accepted(self, edge):
        geometry = _geometry()
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"], excision=edge)
        assert grid["radius"][0] == pytest.approx(edge)

    @pytest.mark.parametrize("edge", [0.70, 0.7542, 1.40, 2.0])
    def test_edges_outside_the_horizons_are_refused(self, edge):
        geometry = _geometry()
        with pytest.raises(pg.ChartAdmissionError, match="strictly between"):
            pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"], excision=edge)

    def test_a_non_finite_edge_is_refused(self):
        geometry = _geometry()
        with pytest.raises(pg.ChartAdmissionError, match="finite"):
            pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"], excision=float("nan"))

    def test_the_edge_stays_above_the_inner_horizon(self):
        assert DEEP_RUN.is_dir()
        summary = json.loads((DEEP_RUN / "summary.json").read_text(encoding="utf-8"))
        assert summary["r_minus"] < summary["r_inner"]
        assert summary["excision_default"] is False
        assert summary["excision_requested"] == pytest.approx(summary["r_inner"])
        assert summary["r_inner"] < summary["horizon_midpoint"]


class TestTheChartAdmitsTheRegion:
    def test_the_lapse_is_unity_across_the_new_region(self):
        """The objection that the chart goes imaginary below r_plus is wrong."""
        geometry = _geometry()
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        r_minus = float(model.horizons["r_minus"])
        radii = np.linspace(r_minus + 1e-3, 1.3, 200)
        mass = float(model.mass)
        charge2 = float(model.magnetic_r2)
        f = 1.0 - 2.0 * mass / radii + charge2 / radii ** 2
        beta2 = 2.0 * mass / radii - charge2 / radii ** 2
        assert np.allclose(f + beta2, 1.0, atol=1e-12)
        assert np.all(beta2 > 0.0)


class TestTheFinding:
    @pytest.fixture(scope="class")
    def deep(self):
        return json.loads(DEEP_SHELL.read_text(encoding="utf-8"))

    @pytest.fixture(scope="class")
    def shallow(self):
        return json.loads(SHALLOW.read_text(encoding="utf-8"))

    def test_the_compression_is_larger_at_the_deeper_edge(self, deep, shallow):
        assert deep["accumulation_after_arrival"]["amplitude_gain"] > \
            shallow["accumulation_after_arrival"]["amplitude_gain"]

    def test_the_shallow_gain_was_a_lower_bound(self, deep, shallow):
        """The earlier number was truncation, not the effect."""
        assert deep["accumulation_after_arrival"]["amplitude_gain"] > 1.25
        assert shallow["accumulation_after_arrival"]["amplitude_gain"] < 1.25

    def test_more_mass_is_retained_with_the_deeper_edge(self, deep, shallow):
        assert deep["mass_budget"]["fraction_below_inner_edge"] > \
            shallow["mass_budget"]["fraction_below_inner_edge"]

    def test_the_relaxation_is_also_truncated_at_the_shallow_edge(self):
        deep = json.loads(DEEP_CURV.read_text(encoding="utf-8"))
        shallow = json.loads((DEVELOPMENT / "nsc-rn-curvature-v2.json").read_text(encoding="utf-8"))
        assert deep["R4_edge_final_over_initial"] > shallow["R4_edge_final_over_initial"]
        assert shallow["R4_edge_final_over_initial"] < 0.05

    def test_curvature_stays_bounded_at_the_deeper_edge(self):
        deep = json.loads(DEEP_CURV.read_text(encoding="utf-8"))
        assert deep["verdict"] == "edge_curvature_bounded_while_mass_compresses"
        assert deep["R4_edge_max_over_initial"] < 10.0

    def test_the_note_records_that_the_edge_is_still_a_boundary(self):
        note = (LAB / "docs/nsc-rn-excision-move.md").read_text(encoding="utf-8")
        assert "still a **boundary**" in note
        assert "is not the horizon" in note
        assert "r_minus = 0.754" in note