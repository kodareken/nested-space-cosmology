"""Joint spectral/charged-throat scale binding.

The light charged Dirac sector is retained explicitly, so the Einstein and
gauge coefficients supplied to the imported charged-throat relation are the
Wilsonian compact complement above the magnetic scale.  The homogeneous
vacuum coefficient is removed only by the already selected relational
zero-tadpole law.  No MMP metric equation is solved here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from math import pi, sqrt
from typing import Any

import numpy as np
from scipy.optimize import brentq

from .nsc_compact_matching import local_coefficients, spectral_moments
from .nsc_mmp_embedding import CommonCoefficients
from .nsc_warped_source import WarpedCompactWeight


@dataclass(frozen=True)
class ScaleBindingConventions:
    """Dimensionless throat-frame conventions for one fixed flux sector."""

    magnetic_flux: int
    interval_over_throat_radius: float = 2.0
    throat_radius: float = 1.0
    warp_amplitude: float = 18.0 / 1015.0
    warp_path: float = 1.0
    compact_modes: int = 34
    quadrature_points: int = 112
    momentum_extent_floor: float = 14.0
    quadrature_tolerance: float = 1e-9

    def __post_init__(self) -> None:
        if (isinstance(self.magnetic_flux, bool)
                or not isinstance(self.magnetic_flux, int)
                or self.magnetic_flux == 0):
            raise ValueError("a nonzero integer magnetic flux is required")
        numeric = (
            self.interval_over_throat_radius,
            self.throat_radius,
            self.warp_amplitude,
            self.warp_path,
            self.momentum_extent_floor,
            self.quadrature_tolerance,
        )
        if not all(np.isfinite(numeric)) or min(numeric[:2]) <= 0:
            raise ValueError("finite positive throat and interval scales required")
        if self.warp_amplitude < 0 or not 0 <= self.warp_path <= 1:
            raise ValueError("invalid compact warp")
        if self.compact_modes < 2 or self.quadrature_points < 2 * self.compact_modes + 4:
            raise ValueError("compact Galerkin representation is unresolved")
        if self.momentum_extent_floor <= 0 or self.quadrature_tolerance <= 0:
            raise ValueError("positive integration controls required")

    @property
    def magnetic_matching_scale(self) -> float:
        """MMP magnetic scale sqrt(|q|)/r_e in throat units."""
        return sqrt(abs(self.magnetic_flux)) / self.throat_radius


@dataclass(frozen=True)
class ScaleBindingBranch:
    magnetic_flux: int
    omega: float
    zeta: float
    coefficients: dict[str, float]
    mmp_map: dict[str, float]
    radius_residual: float
    cutoff_over_magnetic_scale: float
    child_H_over_cutoff: float
    conditional_gap_squared_over_cutoff: float
    conditional_gap_times_throat: float | None
    recursive_endpoint_product: float | None
    numerical: dict[str, Any]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


class ScaleBindingSolver:
    """Solve one fixed-q MMP radius equation with the shared spectral weight."""

    def __init__(self, conventions: ScaleBindingConventions):
        self.conventions = conventions

    @lru_cache(maxsize=64)
    def coefficients(self, omega: float) -> tuple[dict[str, float], dict[str, Any]]:
        c = self.conventions
        if not np.isfinite(omega) or omega <= c.magnetic_matching_scale:
            raise ValueError("the Wilsonian cutoff must exceed the magnetic matching scale")
        cutoff = omega / c.throat_radius
        weight = WarpedCompactWeight(
            cutoff=cutoff,
            interval=c.interval_over_throat_radius * c.throat_radius,
            modes=c.compact_modes,
            points=c.quadrature_points,
            amplitude=c.warp_amplitude,
            path=c.warp_path,
        )
        extent = max(c.momentum_extent_floor / c.throat_radius, 6.0 * cutoff)
        moments = spectral_moments(
            weight, extent=extent, tolerance=c.quadrature_tolerance
        )
        values = local_coefficients(
            weight, moments, c.magnetic_matching_scale
        )
        return values, moments

    def radius_residual(self, omega: float) -> float:
        values, _ = self.coefficients(float(omega))
        q = abs(self.conventions.magnetic_flux)
        radius_squared = q * q * values["C_gauge_Dirac"] / (
            4.0 * values["A_Dirac"]
        )
        return float(radius_squared / self.conventions.throat_radius**2 - 1.0)

    def solve(
        self,
        bracket: tuple[float, float],
        *,
        child_H_times_throat: float = sqrt(3.0 * pi),
        recursive_transfer: float | None = None,
        root_tolerance: float = 2e-11,
    ) -> ScaleBindingBranch:
        low, high = map(float, bracket)
        if not self.conventions.magnetic_matching_scale < low < high:
            raise ValueError("ordered bracket above the magnetic matching scale required")
        left, right = self.radius_residual(low), self.radius_residual(high)
        if left * right >= 0:
            raise ValueError("bracket does not isolate a radius root")
        root = float(brentq(
            self.radius_residual, low, high,
            xtol=root_tolerance, rtol=root_tolerance,
        ))
        values, moments = self.coefficients(root)
        coefficients = CommonCoefficients(
            einstein=values["A_Dirac"],
            gauge=values["C_gauge_Dirac"],
            vacuum=0.0,
        )
        mapped = coefficients.mmp_map(self.conventions.magnetic_flux)
        zeta = (root * self.conventions.throat_radius) ** 2
        gap_squared = 1.0 - 3.0 * pi / (2.0 * zeta)
        gap = sqrt(zeta - 3.0 * pi / 2.0) if gap_squared >= 0 else None
        endpoint = None if recursive_transfer is None else root * recursive_transfer
        return ScaleBindingBranch(
            magnetic_flux=self.conventions.magnetic_flux,
            omega=root,
            zeta=zeta,
            coefficients={
                "A_Wilsonian_at_magnetic_scale": float(values["A_Dirac"]),
                "C_Wilsonian_at_magnetic_scale": float(values["C_gauge_Dirac"]),
                "V_full_relational": 0.0,
                "C_Weyl": float(values["C_Weyl_Dirac"]),
                "C_Euler": float(values["C_Euler_Dirac"]),
                "C_boxR": float(values["C_boxR_Dirac"]),
            },
            mmp_map={key: float(value) for key, value in mapped.items()},
            radius_residual=self.radius_residual(root),
            cutoff_over_magnetic_scale=root / self.conventions.magnetic_matching_scale,
            child_H_over_cutoff=child_H_times_throat / root,
            conditional_gap_squared_over_cutoff=gap_squared,
            conditional_gap_times_throat=gap,
            recursive_endpoint_product=endpoint,
            numerical={
                "compact_modes": self.conventions.compact_modes,
                "quadrature_points": self.conventions.quadrature_points,
                "momentum_extent": max(
                    self.conventions.momentum_extent_floor,
                    6.0 * root * self.conventions.throat_radius,
                ),
                "quadrature_tolerance": self.conventions.quadrature_tolerance,
                "Q1_scaling_residual": float(moments["Q1_scaling_residual"]),
                "Q2_scaling_residual": float(moments["Q2_scaling_residual"]),
                "bracket": [low, high],
            },
        )
