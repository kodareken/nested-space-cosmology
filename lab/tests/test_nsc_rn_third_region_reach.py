"""Tests for the default-off guards that open the third region.

Four functions refused to work below the inner horizon. Each refusal was a
choice in the function rather than a property of the geometry, since the
shift is real down to Q^2/(2M). Each is now overridable, and each override
is off by default.

Two things have to hold at once, and the second is the one that matters
most: the override must work, and nothing that does not ask for it may
change. These tests pin both, and they pin that the guards are still there
so a later change cannot quietly remove them.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

SOURCE = LAB / "src/recursive_horizons/nsc_rn_source.py"
PG = LAB / "src/recursive_horizons/nsc_rn_pg.py"
NOTE = LAB / "docs/nsc-rn-third-region-reach.md"
MU = 1.04
REAL_LIMIT = 1.0 / (2.0 * MU)


@pytest.fixture(scope="module")
def geometry():
    return driver.case_geometry(MU)


class TestTheDefaultsAreUnchanged:
    def test_the_default_excision_is_still_the_midpoint(self, geometry):
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"])
        assert grid["below_inner_horizon"] is False
        assert grid["allow_below_inner_horizon"] is False

    def test_two_default_grids_are_identical(self, geometry):
        args = (geometry["mass"],)
        kwargs = dict(r_m=geometry["r_m"], points=201,
                      r_out=geometry["route"]["outer"])
        a = pg.build_grid(*args, **kwargs)
        b = pg.build_grid(*args, **kwargs)
        assert np.array_equal(a["radius"], b["radius"])

    def test_the_chart_still_refuses_by_default(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        below = float(model.horizons["r_minus"]) * 0.9
        with pytest.raises(ValueError, match="above the inner horizon"):
            source.vacuum_chart(model, np.array([below]))

    def test_the_grid_still_refuses_by_default(self, geometry):
        with pytest.raises(pg.ChartAdmissionError, match="inner and outer horizons"):
            pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"], excision=0.70)


class TestTheOverrideOpensIt:
    def test_the_chart_accepts_below_r_minus_with_the_flag(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        below = float(model.horizons["r_minus"]) * 0.9
        chart = source.vacuum_chart(model, np.array([below]),
                                    allow_below_inner_horizon=True)
        assert np.all(np.isfinite(chart["lapse"]))
        assert np.all(np.isfinite(chart["shift"]))

    def test_the_grid_accepts_below_r_minus_with_the_flag(self, geometry):
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"], excision=0.70,
                             allow_below_inner_horizon=True)
        assert grid["radius"][0] == pytest.approx(0.70)
        assert grid["below_inner_horizon"] is True

    def test_it_still_refuses_below_the_shift_real_limit(self, geometry):
        """The flag moves the bound to Q^2/(2M), not to zero."""
        with pytest.raises(pg.ChartAdmissionError):
            pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                          r_out=geometry["route"]["outer"],
                          excision=REAL_LIMIT * 0.5,
                          allow_below_inner_horizon=True)

    def test_the_grid_reports_where_it_sits(self, geometry):
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"], excision=0.70,
                             allow_below_inner_horizon=True)
        assert grid["below_inner_horizon"]
        assert grid["allow_below_inner_horizon"]


class TestTheRealObstructionIsNumerical:
    """The guards open; the integrator stops at the horizon for another reason."""

    def _packet(self, geometry, lower):
        grid = pg.build_grid(geometry["mass"], r_m=geometry["r_m"], points=201,
                             r_out=geometry["route"]["outer"], excision=lower,
                             allow_below_inner_horizon=True)
        return source.prepare_radial_packet(
            mass_over_rm=MU, energy_target=1.0e-3, rank=driver.PILOT_RANK,
            model=grid["model"], radii=grid["radius"], weights=grid["weights"],
            derivative=grid["derivative"], center=3.0 * geometry["r_m"],
            allow_below_inner_horizon=True)

    @pytest.mark.parametrize("lower", [0.95, 0.80, 0.755])
    def test_it_integrates_down_to_the_inner_horizon(self, geometry, lower):
        self._packet(geometry, lower)

    @pytest.mark.parametrize("lower", [0.70, 0.60])
    def test_it_fails_below_it_on_stiffness(self, geometry, lower):
        with pytest.raises(RuntimeError, match="step size"):
            self._packet(geometry, lower)

    def test_the_boundary_is_the_horizon_itself(self, geometry):
        r_minus = float(pg.action_model(geometry["mass"],
                                         geometry["r_m"]).horizons["r_minus"])
        self._packet(geometry, r_minus + 0.0007)
        with pytest.raises(RuntimeError):
            self._packet(geometry, r_minus - 0.05)

    def test_the_note_says_this_one_is_real(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "This one is real" in flat
        assert "It is a stiff boundary, not a choice" in flat


class TestTheHorizonJetIsGeneralised:
    """The clearing construction applies at both horizons; only r+ was built."""

    def test_it_still_builds_at_the_default_horizon(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        jet = source._horizon_jet(model, 1.0e-3)
        assert jet["radius"] == pytest.approx(float(model.horizons["r_plus"]))

    def test_it_can_build_at_the_inner_horizon(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        jet = source._horizon_jet(model, 1.0e-3, inner)
        assert jet["radius"] == pytest.approx(inner)
        assert np.isfinite(jet["plus"]) and np.isfinite(jet["minus"])

    def test_the_split_is_not_half_applied(self):
        """A reverted change must leave no residue in the solver."""
        text = SOURCE.read_text(encoding="utf-8")
        body = text.split("def scattering_mode")[1]
        assert "inner_jet" not in body


class TestTheGuardsRemainInTheSource:
    """A later change must not quietly remove them."""

    def test_the_chart_text_is_still_there(self):
        assert "above the inner horizon" in SOURCE.read_text(encoding="utf-8")

    def test_the_grid_guard_is_still_there(self):
        assert "child and parent cuts are not ordered" in PG.read_text(encoding="utf-8")

    def test_the_override_is_defaulted_off_in_both_modules(self):
        assert "allow_below_inner_horizon=False" in SOURCE.read_text(encoding="utf-8")
        assert "allow_below_inner_horizon=False" in PG.read_text(encoding="utf-8")