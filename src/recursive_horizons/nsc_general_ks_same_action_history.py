"""Executable boundary variation for the homogeneous NSC KS history class.

The lapse/shift variation is evaluated before an interior solve.  In a
homogeneous Kantowski--Sachs history every already-declared local metric and
magnetic term has zero normal--axial momentum.  A continuous CTP state with a
nonzero seed T01 therefore excludes a stationary history in this class.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite


@dataclass(frozen=True)
class GeneralKSSameActionHistoryFunctional:
    einstein_coefficient: float
    magnetic_flux: int
    omega: float
    zeta: float
    vacuum_coefficient: float

    def __post_init__(self):
        if (not isfinite(self.einstein_coefficient)
                or self.einstein_coefficient <= 0):
            raise ValueError("positive finite Einstein coefficient required")
        if isinstance(self.magnetic_flux, bool) or not isinstance(self.magnetic_flux, int):
            raise ValueError("integer magnetic flux required")
        if not all(isfinite(value) for value in (
            self.omega, self.zeta, self.vacuum_coefficient
        )):
            raise ValueError("finite locked scale data required")

    def homogeneous_shift_variation(
        self,
        *,
        state_T01: float,
        compact_local_T01: float = 0.0,
        magnetic_T01: float = 0.0,
        vacuum_T01: float = 0.0,
        curvature_reference_T01: float = 0.0,
        boundary_shell_T01: float = 0.0,
    ) -> dict[str, float]:
        values = (
            state_T01, compact_local_T01, magnetic_T01, vacuum_T01,
            curvature_reference_T01, boundary_shell_T01,
        )
        if not all(isfinite(value) for value in values):
            raise ValueError("finite same-action momentum ledger required")
        total = sum(values)
        return {
            "state_T01": state_T01,
            "compact_local_T01": compact_local_T01,
            "magnetic_T01": magnetic_T01,
            "vacuum_T01": vacuum_T01,
            "curvature_reference_T01": curvature_reference_T01,
            "boundary_shell_T01": boundary_shell_T01,
            "total_T01": total,
            "normalized_shift_variation": -total/(2*self.einstein_coefficient),
            "counterflow_required": -total,
        }

    def homogeneous_nonexistence_certificate(self, *, state_T01: float, tolerance: float):
        variation = self.homogeneous_shift_variation(state_T01=state_T01)
        exists = abs(variation["normalized_shift_variation"]) <= tolerance
        return {
            "stationary_history_exists_in_declared_class": exists,
            "variation": variation,
            "continuity_argument": (
                "unitary finite-dimensional C(t) and smooth metric vertices make "
                "T01(t) continuous; nonzero T01 at the seed persists on a "
                "neighborhood, while the homogeneous shift constraint requires "
                "T01(t)=0 on every slice"
            ),
            "tolerance": tolerance,
        }
