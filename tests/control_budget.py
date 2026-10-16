"""Cumulative control allowances charged only while their tests are active.

Imported pytest fixtures can share this ledger without charging collection or
unrelated tests. During an activation the ordinary scientific Budget (and any
children sharing its start) still enforces both limits; neither the accumulated
point count nor the CPU time already used by these controls is reset.
"""
from contextlib import contextmanager
import time

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


class ControlBudget:
    def __init__(self, *, maximum_evaluations, cpu_seconds):
        self.maximum_evaluations = maximum_evaluations
        self.cpu_seconds = cpu_seconds
        self.budget = None
        self.active_cpu_seconds = 0.
        self._depth = 0

    @contextmanager
    def active(self):
        if self._depth == 0:
            now = time.process_time()
            if self.budget is None:
                self.budget = Budget(maximum_evaluations=self.maximum_evaluations,
                                     cpu_seconds=self.cpu_seconds)
            self.budget.started = now-self.active_cpu_seconds
            self.budget.check()
        self._depth += 1
        try:
            yield self.budget
        finally:
            self._depth -= 1
            if self._depth == 0:
                self.active_cpu_seconds = time.process_time()-self.budget.started
                self.budget.check()
