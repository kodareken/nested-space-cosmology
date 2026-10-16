"""Test accounting itself without spending the declared CPU allowance."""
import pytest

import control_budget


def test_unrelated_cpu_is_excluded_but_control_cpu_and_points_accumulate(monkeypatch):
    now = [100.]
    monkeypatch.setattr(control_budget.time, "process_time", lambda: now[0])
    ledger = control_budget.ControlBudget(maximum_evaluations=10, cpu_seconds=3.)
    now[0] = 1000.  # collection/import time does not start the allowance
    with ledger.active() as first:
        first.check(4)
        now[0] += 1.
    now[0] += 500.  # unrelated tests do not consume this control's allowance
    with ledger.active() as second:
        assert second is first and second.evaluations == 4
        assert now[0]-second.started == 1.
        second.check(3)
        now[0] += 1.
    assert ledger.active_cpu_seconds == 2. and first.evaluations == 7
    with pytest.raises(RuntimeError, match="budget exhausted"):
        with ledger.active():
            now[0] += 2.


def test_nested_imported_fixture_uses_one_clock_and_one_point_ceiling(monkeypatch):
    now = [100.]
    monkeypatch.setattr(control_budget.time, "process_time", lambda: now[0])
    ledger = control_budget.ControlBudget(maximum_evaluations=5, cpu_seconds=3.)
    with ledger.active() as first:
        started = first.started
        now[0] += .5
        with ledger.active() as second:
            assert second is first and second.started == started
            second.check(5)
        now[0] += .5
    assert ledger.active_cpu_seconds == 1.
    with pytest.raises(RuntimeError, match="budget exhausted"):
        with ledger.active() as third:
            third.check(1)
