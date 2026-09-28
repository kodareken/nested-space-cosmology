from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import sys
import tomllib
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration_v7 import Proto7RunMember  # noqa: E402
from scripts.run_fgc_gr0_calibration_v8 import (  # noqa: E402
    Proto8RunMember,
    _event_assessment,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)


def _member(method: str, order: int, count: int, time: float = 0.0):
    grid = UniformRadialGrid(0.0, 128.0, count)
    u = np.zeros((count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    state = EvolutionState(u, p, q)
    key = f"{method}-{count}"
    return key, SimpleNamespace(
        key=key,
        amplitude="5/2",
        method_label=method,
        spatial_order=order,
        point_count=count,
        state=state,
        initial=SimpleNamespace(grid=grid),
        time=time,
        step_index=0,
        transaction_serial=0,
        CFL_retry_count=0,
        source_retry_count=0,
        transaction=SimpleNamespace(
            state=SimpleNamespace(accepted_stage_count=0),
            causal_state=SimpleNamespace(accumulated_characteristic_distance=0.0),
            boundary_geometry=SimpleNamespace(
                outer_radius=128.0,
                measurement_radius=24.0,
                sbp_stencil_reach=order * grid.spacing,
                minimum_causal_buffer=16.0,
            ),
        ),
        tracers=SimpleNamespace(
            positions=np.arange(0.5, 24.5, 0.5),
            proper_times=np.zeros(48),
            event_proper_times=[],
            event_fields=[],
        ),
    )


def _members(time: float = 0.0):
    return dict(
        _member(method, order, count, time)
        for method, order in (("RK4", 4), ("SSPRK3", 2))
        for count in (1025, 2049, 4097)
    )


class FGCGR0CampaignRunnerV8Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (REPOSITORY / "configs/fgc/fgc-1-cal5-run1.toml").open("rb") as handle:
            cls.plan = tomllib.load(handle)

    def test_evolution_member_is_exactly_the_inherited_PROTO7_member(self) -> None:
        self.assertIs(Proto8RunMember, Proto7RunMember)
        self.assertIs(Proto8RunMember.advance_to, Proto7RunMember.advance_to)

    def test_common_event_uses_PROTO8_admission_without_advancing(self) -> None:
        members = _members()
        before = {
            key: (item.state, item.time, item.step_index)
            for key, item in members.items()
        }
        assessed = _event_assessment(self.plan, members, 0.0)
        primary = assessed["primary_common_event"]
        comparator = assessed["comparator_common_event"]
        self.assertEqual(primary.accepted_stage_counts, (0, 0, 0))
        self.assertEqual(comparator.accepted_stage_counts, (0, 0, 0))
        self.assertEqual(
            primary.admission_passed,
            primary.constraint_admission.admission_passed
            and primary.spatial_spectral_admission.admission_passed,
        )
        self.assertEqual(
            comparator.admission_passed,
            comparator.constraint_admission.admission_passed
            and comparator.spatial_spectral_admission.admission_passed,
        )
        self.assertTrue(primary.spatial_spectral_admission.admission_passed)
        self.assertTrue(comparator.spatial_spectral_admission.admission_passed)
        self.assertFalse(assessed["temporal_spectral_admission"]["available"])
        self.assertFalse(assessed["qualified_trapped_common_event"])
        self.assertEqual(
            primary.constraint_admission.excluded_outer_rows,
            (7, 7, 7),
        )
        self.assertEqual(
            comparator.constraint_admission.excluded_outer_rows,
            (5, 5, 5),
        )
        for key, item in members.items():
            self.assertIs(item.state, before[key][0])
            self.assertEqual((item.time, item.step_index), before[key][1:])

    def test_bitwise_common_time_misalignment_fails_before_diagnostics(self) -> None:
        members = _members()
        members["RK4-1025"].time = np.nextafter(0.0, 1.0)
        with self.assertRaisesRegex(ValueError, "bitwise aligned"):
            _event_assessment(self.plan, members, 0.0)


if __name__ == "__main__":
    unittest.main()
