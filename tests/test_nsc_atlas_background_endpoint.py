"""Manufactured finite-profile endpoint controls; no prepared scalar solve."""
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_background_evolution import (
    background_stage, evolve_background, quiet_departures,
)
from test_nsc_atlas_reference_wall_gauge import GEOMETRY


TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}


@pytest.fixture(scope="module")
def context():
    budget = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    gauge = ReferenceWallGauge(GEOMETRY, budget=budget)
    yield budget, gauge
    budget.check()
    print("background endpoint controls:", budget.evaluations, "points;",
          time.process_time()-budget.started, "process CPU seconds")


def source(x):
    return {"r": x, "log_C": 0., "c": 0., "r_x": 1., "log_C_x": 0.,
            "c_x": 0., "source_constraint": 0.}


def finite_profile(last):
    return CharacteristicProfile([0., last], [0., 0.], [0., 0.], [0., 0.])


def test_exact_terminal_stage_has_no_spurious_coverage_requirement(context):
    budget, gauge = context
    x = .42
    reference = gauge.sample(x)
    last = float(reference["T"]+np.arctanh(reference["r"]))
    eta, y = np.linspace(0., 1., 5), np.linspace(0., 1., 4)
    state = quiet_departures(exterior_cells=4, cap_cells=3)
    result = background_stage(state, x=x, eta=eta, y_faces=y,
        source=source(x), gauge=reference, alpha=0., tube=TUBE, budget=budget,
        profile=finite_profile(last))
    assert result["physical_wall"]["T"]+np.arctanh(result["physical_wall"]["r"]) == last
    np.testing.assert_array_equal(result["physical_cap"]["Iq"], 0.)
    # A genuine floating-coordinate excursion is rejected rather than hidden
    # by extrapolating the finite profile. This is a domain-policy control.
    shifted = quiet_departures(exterior_cells=4, cap_cells=3)
    shifted["wall"]["r"] = float(np.nextafter(reference["r"], 1.)-reference["r"])
    shifted_label = reference["T"]+np.arctanh(reference["r"]+shifted["wall"]["r"])
    assert shifted_label > last
    with pytest.raises(ValueError, match="finite labels"):
        background_stage(shifted, x=x, eta=eta, y_faces=y,
            source=source(x), gauge=reference, alpha=0., tube=TUBE, budget=budget,
            profile=finite_profile(last))


def test_quiet_evolution_with_profile_ending_at_reference_final_label(context):
    budget, gauge = context
    x = np.linspace(.4, .42, 5)
    end = gauge.sample(x[-1])
    last = float(end["T"]+np.arctanh(end["r"]))
    eta, y = np.linspace(0., 1., 5), np.linspace(0., 1., 4)
    row = evolve_background(quiet_departures(exterior_cells=4, cap_cells=3), x,
        eta=eta, y_faces=y, source=source, gauge=gauge.sample, alpha=0.,
        tube=TUBE, budget=budget, profile=finite_profile(last))
    assert row["steps"] == 4 and row["stages"] == 9
    final = row["snapshots"][-1]
    assert final["physical"]["wall"]["T"]+np.arctanh(final["physical"]["wall"]["r"]) <= last
    assert max(abs(float(v)) for v in final["state"]["wall"].values()) < 4e-14
    np.testing.assert_array_equal(final["physical"]["cap"]["Iq"], 0.)
