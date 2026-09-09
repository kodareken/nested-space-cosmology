"""Projected magnetic Dirac/gauge interaction on the recorded exterior.

Reuse QED2 bosonization: the charged collective mode has mu²=q*g4²/(4*pi²*r²).
The background and coupling are inputs. This is not a full 4D reduction or
a self-sourced geometry. No fermion rest mass or added scalar is introduced.
"""
from dataclasses import dataclass
from math import atan, asinh, cosh, exp, log, pi, sqrt

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_lorentzian import geometry


@dataclass(frozen=True)
class Exterior:
    horizon_rho: float
    surface_gravity: float

    def A(self, rho):
        """Same recorded A, stably evaluated near the root and far infinity."""
        delta = rho-self.horizon_rho
        if delta <= 0:
            raise ValueError("the scattering chart is outside the recorded horizon")
        if delta < 1e-4:
            h = self.horizon_rho
            second = -6*atan(1/h)+6*h/(1+h*h)
            third = 12/(1+h*h)**2
            fourth = -48*h/(1+h*h)**3
            return 2*self.surface_gravity*delta+second*delta**2/2+third*delta**3/6+fourth*delta**4/24
        if rho > 20:
            return 1-sum(6*(-1)**n/((2*n+1)*(2*n+3)*rho**(2*n+1)) for n in range(6))
        return float(geometry(rho)[2])


def strength(gauge_coupling, flux=1):
    if not np.isfinite(gauge_coupling) or gauge_coupling < 0:
        raise ValueError("finite nonnegative gauge coupling required")
    if not isinstance(flux, int) or isinstance(flux, bool) or flux == 0:
        raise ValueError("nonzero integer flux required for an LLL sector")
    return abs(flux)*gauge_coupling**2/(4*pi*pi)


def transmission(background, frequency, gauge_coupling, flux=1,
                 phase_extent=60., horizon_offset=1e-8, rtol=2e-9,
                 method="amplitudes"):
    """Current-normalized charge-mode transmission with analytic tail bounds.

    The potential is V=A*strength/r² and dx*/d rho=1/A. The main solver
    uses incoming/outgoing amplitudes; the independent control evolves
    the scalar and its tortoise derivative directly. No horizon is re-solved.
    """
    if not np.isfinite(frequency) or frequency <= 0 or phase_extent <= 0:
        raise ValueError("positive frequency and exterior phase extent required")
    if horizon_offset <= 0 or rtol <= 0 or method not in ("amplitudes", "field"):
        raise ValueError("positive offset/tolerance and a declared solver required")
    c = strength(gauge_coupling, flux)
    if c == 0:
        return {"transmission": 1., "reflection": 0., "raw_transmission": 1.,
                "relative_current_residual": 0., "tail_rapidity_bound": 0.,
                "transmission_lower": 1., "transmission_upper": 1.}
    h = background.horizon_rho
    end = max(50., phase_extent/frequency)
    begin, finish = log(horizon_offset), log(end-h)

    def amplitude_rhs(u, state):
        delta = exp(u);rho = h+delta
        a, b, phase = state
        angle = 2*phase.real
        factor = c*delta/(2*frequency*(1+rho*rho))
        return np.array([1j*factor*(a+b*np.exp(1j*angle)),
                         -1j*factor*(a*np.exp(-1j*angle)+b),
                         frequency*delta/background.A(rho)], dtype=complex)

    def field_rhs(u, state):
        delta = exp(u);rho = h+delta;A = background.A(rho)
        psi, derivative = state
        return np.array([delta*derivative/A,
                         delta*(c/(1+rho*rho)-frequency**2/A)*psi], dtype=complex)

    if method == "amplitudes":
        initial = [1.+0j, 0j, 0j]
        rhs = amplitude_rhs
    else:
        initial = [1.+0j, -1j*frequency]
        rhs = field_rhs
    solution = solve_ivp(rhs, (begin, finish), initial, method="DOP853",
                         rtol=rtol, atol=rtol*.02, max_step=.25)
    if not solution.success:
        raise RuntimeError(solution.message)
    if method == "amplitudes":
        a, b, _ = solution.y[:, -1]
    else:
        psi, derivative = solution.y[:, -1]
        a, b = (psi+1j*derivative/frequency)/2, (psi-1j*derivative/frequency)/2
    a2, b2 = abs(a)**2, abs(b)**2
    # The exact Wronskian is |a|²-|b|²=1. This stable formula avoids T>1
    # from roundoff when reflection is tiny; raw extraction stays visible.
    value = 1/(1+b2)
    theta = asinh(abs(b))
    area_left = atan(h+horizon_offset)-atan(h)
    area_right = atan(1/end)
    tail = c*(area_left+area_right)/(2*frequency)
    return {"transmission": float(value), "reflection": float(b2/(1+b2)),
            "raw_transmission": float(1/a2),
            "relative_current_residual": float(abs(a2-b2-1)/(1+a2+b2)),
            "tail_rapidity_bound": tail,
            "transmission_lower": 1/cosh(theta+tail)**2,
            "transmission_upper": 1/cosh(max(0.,theta-tail))**2}


def outgoing_power(background, gauge_coupling, flux=1, nodes=64,
                   minimum_u=1e-4, maximum_u=24., phase_extent=60.,
                   horizon_offset=1e-8, rtol=2e-9):
    """Conditional Unruh power for the projected interacting sector.

    u=omega/T_H. The q-1 neutral units of central charge are massless;
    the one charged mode is scattered by the calculated variable potential.
    Bounds cover omitted frequency/domain tails, not the ODE or quadrature
    errors, which are checked separately. This is not cosmological Q.
    """
    strength(gauge_coupling, flux)
    if not isinstance(nodes, int) or nodes < 8 or not 0 < minimum_u < maximum_u:
        raise ValueError("resolved positive frequency interval required")
    temperature = background.surface_gravity/(2*pi)
    single_free = background.surface_gravity**2/(48*pi)
    if gauge_coupling == 0:
        # Reuse the exact free CFT result instead of numerically truncating it.
        return {"charged_power_over_one_free_channel": 1.,
                "charged_ratio_lower_with_tail_bounds": 1.,
                "charged_ratio_upper_with_tail_bounds": 1.,
                "missing_IR_ratio_bound": 0., "missing_UV_ratio_bound": 0.,
                "maximum_relative_current_residual": 0.,
                "neutral_central_charge": abs(flux)-1,
                "parent_Killing_power": abs(flux)*single_free,
                "total_power_over_free_LLL": 1.}
    x, weights = leggauss(nodes)
    lo, hi = log(minimum_u), log(maximum_u)
    locations = np.exp(lo+(x+1)*(hi-lo)/2)
    weights = weights*(hi-lo)/2
    integrals = np.zeros(3)
    maximum_current = 0.
    for u, weight in zip(locations, weights):
        response = transmission(background, u*temperature, gauge_coupling, flux,
                                phase_extent, horizon_offset, rtol)
        factor = 6/(pi*pi)*weight*u*u/np.expm1(u)
        integrals += factor*np.array([response["transmission"],response["transmission_lower"],response["transmission_upper"]])
        maximum_current = max(maximum_current, response["relative_current_residual"])
    missing_IR = 6*minimum_u/(pi*pi)
    missing_UV = 6/(pi*pi)*(maximum_u+1)*exp(-maximum_u)/(1-exp(-maximum_u))
    neutral = abs(flux)-1
    return {"charged_power_over_one_free_channel": float(integrals[0]),
            "charged_ratio_lower_with_tail_bounds": float(integrals[1]),
            "charged_ratio_upper_with_tail_bounds": float(integrals[2]+missing_IR+missing_UV),
            "missing_IR_ratio_bound": missing_IR, "missing_UV_ratio_bound": missing_UV,
            "maximum_relative_current_residual": maximum_current,
            "neutral_central_charge": neutral,
            "parent_Killing_power": float(single_free*(neutral+integrals[0])),
            "total_power_over_free_LLL": float((neutral+integrals[0])/abs(flux))}
