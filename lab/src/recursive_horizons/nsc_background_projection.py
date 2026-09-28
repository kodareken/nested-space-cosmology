"""Proper-volume and child-clock projection of a completed neck tensor.

The stress tensor is an input.  This module does not construct a state, solve
the metric equations, or extend the tensor away from the neck.  It applies the
standard homogeneous Kantowski--Sachs conservation ledger on the already
stored Bronnikov child branch.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import atan, isfinite, pi, sqrt


@dataclass(frozen=True)
class ChildFrameTensor:
    """Angularly homogeneous orthonormal stress at one child-time slice."""

    rho: float
    T01: float
    p_parallel: float
    p_sphere: float
    parent_Killing_power: float

    def __post_init__(self) -> None:
        if not all(isfinite(value) for value in asdict(self).values()):
            raise ValueError("finite child-frame tensor required")

    def plus(self, other: "ChildFrameTensor") -> "ChildFrameTensor":
        return ChildFrameTensor(**{
            name: getattr(self, name) + getattr(other, name)
            for name in asdict(self)
        })


def bronnikov_child_slice(rho_coordinate: float) -> dict[str, float]:
    """Stored child geometry in proper time, for rho_coordinate <= 0.

    The imported metric is
      A=1+3 rho+3(1+rho^2)(atan(rho)-pi/2),
      a_parallel=sqrt(-A), r=sqrt(1+rho^2), dT=-d rho/a_parallel.
    """
    if not isfinite(rho_coordinate) or rho_coordinate > 0:
        raise ValueError("finite expanding-child coordinate rho<=0 required")
    rho = float(rho_coordinate)
    metric_A = 1.0 + 3.0*rho + 3.0*(1.0+rho*rho)*(atan(rho)-pi/2.0)
    if metric_A >= 0:
        raise ValueError("child proper-time chart requires A<0")
    metric_A_prime = 6.0*(1.0-rho*(pi/2.0-atan(rho)))
    a_parallel = sqrt(-metric_A)
    sphere_radius = sqrt(1.0+rho*rho)
    H_parallel = metric_A_prime/(2.0*a_parallel)
    H_sphere = -a_parallel*rho/(1.0+rho*rho)
    volume = 4.0*pi*a_parallel*sphere_radius*sphere_radius
    expansion = H_parallel+2.0*H_sphere
    return {
        "rho_coordinate": rho,
        "A": metric_A,
        "A_prime": metric_A_prime,
        "d_rho_d_child_time": -a_parallel,
        "a_parallel": a_parallel,
        "sphere_radius": sphere_radius,
        "H_parallel": H_parallel,
        "H_sphere": H_sphere,
        "volume_expansion": expansion,
        "proper_volume_per_coordinate_dz": volume,
        "proper_volume_rate_per_coordinate_dz": volume*expansion,
        "volume_scale_H": expansion/3.0,
    }


def conserved_local_jet(
    tensor: ChildFrameTensor,
    geometry: dict[str, float],
    *,
    energy_transfer: float = 0.0,
    momentum_transfer: float = 0.0,
) -> dict[str, float]:
    """Apply the imported homogeneous conservation equations at one slice.

    energy_transfer and momentum_transfer are proper-volume source densities.
    They are zero for the locked free, block-diagonal Gaussian realization.
    """
    if not isfinite(energy_transfer) or not isfinite(momentum_transfer):
        raise ValueError("finite transfer rates required")
    hp = geometry["H_parallel"]
    hs = geometry["H_sphere"]
    theta = geometry["volume_expansion"]
    volume = geometry["proper_volume_per_coordinate_dz"]
    a_parallel = geometry["a_parallel"]

    dilution = -theta*tensor.rho
    pressure_work_density = -(hp*tensor.p_parallel+2.0*hs*tensor.p_sphere)
    rho_rate = dilution+pressure_work_density+energy_transfer
    momentum_rate = -2.0*(hp+hs)*tensor.T01+momentum_transfer
    cell_energy = volume*tensor.rho
    cell_energy_rate = volume*rho_rate + geometry["proper_volume_rate_per_coordinate_dz"]*tensor.rho
    pressure_work_rate = volume*pressure_work_density
    killing_momentum = a_parallel*volume*tensor.T01
    killing_momentum_rate = a_parallel*volume*momentum_transfer
    oriented_parent_power = -killing_momentum

    energy_residual = (
        rho_rate + hp*(tensor.rho+tensor.p_parallel)
        + 2.0*hs*(tensor.rho+tensor.p_sphere) - energy_transfer
    )
    integrated_residual = (
        cell_energy_rate
        + volume*(hp*tensor.p_parallel+2.0*hs*tensor.p_sphere-energy_transfer)
    )
    momentum_residual = (
        momentum_rate+2.0*(hp+hs)*tensor.T01-momentum_transfer
    )
    return {
        "energy_transfer_Q": float(energy_transfer),
        "momentum_transfer": float(momentum_transfer),
        "volume_dilution_density_rate": float(dilution),
        "directional_pressure_work_density_rate": float(pressure_work_density),
        "density_rate_child_time": float(rho_rate),
        "momentum_density_rate_child_time": float(momentum_rate),
        "proper_energy_per_coordinate_dz": float(cell_energy),
        "proper_energy_rate_per_coordinate_dz": float(cell_energy_rate),
        "directional_pressure_work_rate_per_coordinate_dz": float(pressure_work_rate),
        "density_derivative_per_log_volume": float(rho_rate/theta),
        "density_derivative_per_log_volume_scale": float(3.0*rho_rate/theta),
        "child_spatial_Killing_charge_per_coordinate_dz": float(killing_momentum),
        "child_spatial_Killing_charge_rate": float(killing_momentum_rate),
        "oriented_parent_power_from_child_momentum": float(oriented_parent_power),
        "energy_conservation_residual": float(energy_residual),
        "integrated_energy_ledger_residual": float(integrated_residual),
        "momentum_conservation_residual": float(momentum_residual),
        "stored_Killing_power_map_residual": float(
            oriented_parent_power-tensor.parent_Killing_power
        ),
    }
