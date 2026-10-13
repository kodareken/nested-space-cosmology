"""Replay arithmetic cannot alter scientific bounds or hide unrelated drift."""
import copy
import json
from pathlib import Path

import pytest

from recursive_horizons import nsc_normalization_flow as flow
from recursive_horizons.nsc_normalization_arithmetic import (
    euler_replay_roundoff_envelope, validate_euler_summary,
)


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("index,points", [(0, 32), (1, 64), (2, 128)])
def test_diagnostic_arithmetic_envelope_has_a_negative_control(index, points):
    path = f"$/discrete_euler_refinement/grids/{index}/discrete_E4_el_maximum"
    envelope = euler_replay_roundoff_envelope(points)
    assert envelope > 0
    flow.compare_record(0., envelope / 2, path)
    with pytest.raises(RuntimeError):
        flow.compare_record(0., envelope * 2, path)
    with pytest.raises(RuntimeError):
        flow.compare_record(0., -envelope / 2, path)
    # Adjacent scientific quantities do not receive this diagnostic policy.
    with pytest.raises(RuntimeError):
        flow.compare_record(0., 1e-10, path.replace("discrete_E4_el_maximum", "c2_only_gradient_residual"))


def test_envelope_is_independent_of_the_physics_tolerance(monkeypatch):
    before = [euler_replay_roundoff_envelope(n) for n in (32, 64, 128)]
    monkeypatch.setattr(flow, "ABS_TOL", 1.)
    euler_replay_roundoff_envelope.cache_clear()
    after = [euler_replay_roundoff_envelope(n) for n in (32, 64, 128)]
    assert before == after


def test_diagnostic_mutation_inside_loose_envelope_still_fails_summary_identity():
    saved = json.loads((ROOT / "results/nsc-12-normalization-flow.json").read_text())
    altered = copy.deepcopy(saved)
    row = altered["discrete_euler_refinement"]["grids"][2]
    row["discrete_E4_el_maximum"] += euler_replay_roundoff_envelope(128) / 2
    with pytest.raises(ValueError, match="inconsistent"):
        validate_euler_summary(row)
    with pytest.raises(RuntimeError, match="summary contract"):
        flow.compare_record(saved, altered)


def test_original_bounds_scope_and_source_identities_remain_exact():
    saved = json.loads((ROOT / "results/nsc-12-normalization-flow.json").read_text())
    assert saved["reproduction_tolerance"]["float_atol"] == flow.ABS_TOL == 2e-12
    assert saved["reproduction_tolerance"]["float_rtol"] == flow.REL_TOL == 2e-12
    for change in ("source", "gate", "criterion"):
        altered = copy.deepcopy(saved)
        if change == "source":
            key = next(iter(altered["source_hashes"]))
            altered["source_hashes"][key] = "0" * 64
        elif change == "gate":
            altered["gate"]["discrete_euler_residual_vanishes_under_refinement"] = False
        else:
            altered["discrete_euler_refinement"]["grids"][2]["c2_only_gradient_residual"] += 1e-10
        with pytest.raises(RuntimeError):
            flow.compare_record(saved, altered)


def test_undeclared_grid_has_no_portability_exemption():
    with pytest.raises(ValueError, match="declared three grids"):
        euler_replay_roundoff_envelope(256)
