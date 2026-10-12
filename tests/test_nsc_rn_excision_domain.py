"""Tests for where the sourced march is admissible.

The first four gates on the way below the inner horizon were conventions and
were opened. The fifth is different: the excision boundary condition requires
a pure outflow, and below the inner horizon the outgoing speed is positive,
so the boundary would emit. That is a property of the geometry rather than
of the code, and these tests pin it as such.
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

NOTE = LAB / "docs/nsc-rn-excision-domain.md"
RUNNER = LAB / "scripts/run_nsc_rn_inner_packet.py"
MU = 1.04


@pytest.fixture(scope="module")
def geometry():
    return driver.case_geometry(MU)


def _speeds(model, edge):
    chart = source._point_chart(model, edge, allow_below_inner_horizon=True)
    return chart["speed_plus"], chart["speed_minus"]


class TestTheSignReversesAtTheInnerHorizon:
    def test_outgoing_is_negative_above_it(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        for edge in (1.04, 0.90, 0.80, inner + 1e-3):
            outgoing, ingoing = _speeds(model, edge)
            assert outgoing < 0.0 and ingoing < 0.0

    def test_outgoing_is_positive_below_it(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        for edge in (inner - 1e-2, 0.70, 0.60, 0.50):
            outgoing, _ = _speeds(model, edge)
            assert outgoing > 0.0

    def test_the_crossing_is_at_the_horizon(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        above, _ = _speeds(model, inner + 1e-3)
        below, _ = _speeds(model, inner - 1e-3)
        assert above < 0.0 < below


class TestTheChartRefusesBelow:
    def test_admit_chart_accepts_above_the_inner_horizon(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        grid = {"spacing": 0.02}
        chart = source.vacuum_chart(model, np.array([inner + 0.01, 8.0]),
                                    allow_below_inner_horizon=True)
        admitted = pg.admit_chart(grid, chart["lapse"], chart["shift"])
        assert admitted["admitted"]

    def test_it_refuses_below_the_inner_horizon(self, geometry):
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        inner = float(model.horizons["r_minus"])
        grid = {"spacing": 0.02}
        chart = source.vacuum_chart(model, np.array([inner - 0.05, 8.0]),
                                    allow_below_inner_horizon=True)
        with pytest.raises(pg.ChartAdmissionError, match="pure outflow"):
            pg.admit_chart(grid, chart["lapse"], chart["shift"])


class TestThisIsNotAFifthConvention:
    def test_the_note_says_the_domain_is_the_trapped_region(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "maps exactly onto the trapped region" in flat

    def test_the_note_says_it_is_not_a_guard_to_flip(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "not a matter of flipping a guard" in flat

    def test_the_note_leaves_the_region_open(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "remains open, and it is now open for the right reason" in flat

    def test_the_runner_flag_is_still_required_below_the_horizon(self):
        text = RUNNER.read_text(encoding="utf-8")
        assert "--below-horizon to go there deliberately" in text

    def test_the_runner_threads_the_flag_to_both_calls(self):
        text = RUNNER.read_text(encoding="utf-8")
        assert text.count("allow_below_inner_horizon=args.below_horizon") == 2


class TestScope:
    def test_it_says_it_is_not_the_universe(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "does not measure the actual universe" in flat

    def test_it_uses_no_station_payloads(self):
        flat = NOTE.read_text(encoding="utf-8").replace(chr(10), " ")
        assert "No station payloads" in flat