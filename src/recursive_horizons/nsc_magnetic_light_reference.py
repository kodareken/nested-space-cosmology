"""Charged light-field restoration in the existing mu=1 normalization.

The only physical charge evaluated here is the retained |q|=4 sector; q=0
is an exact neutral normalization control. The angular spectrum is the owned
monopole spectrum. No mode, scattering, metric or state evolution is performed.

The optional homogeneous map restores NONZERO angular modes only. It uses
the existing neutral four-dimensional conformal owner plus the charged gauge
cocycle, and subtracts the two-dimensional LLL cocycle so that the separately
owned PHYSICAL-metric LLL is counted once by the final source assembly.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import pi

import mpmath as mp
import numpy as np

from .nsc_angular_stress import profile_jets, curvature_squared_tensor, physical_source


def _charge(value):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or abs(value) not in (0, 4):
        raise ValueError("only the retained |q|=4 light field and q=0 neutral control are owned")
    return abs(int(value))


def _mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


def derivative_series_tail_bound(charge, first_omitted):
    """Exact rational upper bound on the omitted part of Z_q'(-1).

    For a=b+1, r=(b/a)^2 and J>=3, use
      zeta(2j-3,a) <= a^(3-2j) [1+a/(2j-4)]
    and bound the remaining geometric series. This bounds spectral-series
    truncation only, not floating-point evaluation of the Hurwitz constants.
    """
    q = _charge(charge)
    if isinstance(first_omitted, bool) or not isinstance(first_omitted, (int, np.integer)) or first_omitted < 3:
        raise ValueError("first omitted Hurwitz-series index must be an integer >=3")
    if q == 0:
        return Fraction(0)
    b = Fraction(q, 2); a = b+1; ratio = (b/a)**2; J = int(first_omitted)
    return (4*a**3*ratio**J/(J*(J-1)*(1-ratio))
            *(1+a/Fraction(2*J-4)))


@dataclass(frozen=True)
class MagneticLightSpectrum:
    charge: int
    z_minus_one: Fraction
    z_nonzero_zero: Fraction
    harmonic_number: Fraction
    zprime_minus_one: object
    positive_series_partial: object
    derivative_tail_bound: Fraction
    series_terms: int
    last_series_index: int
    decimal_precision: int

    @property
    def finite_part_one(self):
        """FP Z_q(1)=4 gamma_E-2 H_q, with residue two."""
        with mp.workdps(self.decimal_precision):
            return +(4*mp.euler-2*_mp_fraction(self.harmonic_number))

    @property
    def fourth_order_harmonic(self):
        """ln(mu)+FP/8 at the already fixed mu=1 and barred sphere R=1."""
        with mp.workdps(self.decimal_precision):
            return +(self.finite_part_one/8)

    def cylinder_density(self, radius=1):
        """Same proper-time finite normalization as the neutral cylinder.

        R is the reference-cylinder sphere radius; the homogeneous physical
        map below uses barred R=1 and handles its conformal factor separately.
        The subtraction normalization is fixed at mu=1, never fitted.
        """
        with mp.workdps(self.decimal_precision):
            R = mp.mpf(str(radius))
            if not mp.isfinite(R) or R <= 0:
                raise ValueError("positive finite reference-cylinder radius required")
            return +(-(self.zprime_minus_one
                       +(mp.log(R*R)+1-mp.euler)*_mp_fraction(self.z_minus_one))
                     /(16*mp.pi**2*R**4))

    def cylinder_series_error_bound(self, radius=1):
        with mp.workdps(self.decimal_precision):
            R = mp.mpf(str(radius))
            if not mp.isfinite(R) or R <= 0:
                raise ValueError("positive finite reference-cylinder radius required")
            return +(_mp_fraction(self.derivative_tail_bound)/(16*mp.pi**2*R**4))


@lru_cache(maxsize=8)
def magnetic_light_spectrum(charge=4, *, series_tolerance="1e-32", decimal_precision=70, max_terms=256):
    """Evaluate the new finite charged spectral sum, with a rational tail bound.

    Z_q(s)=4 sum_{l=b+1}^infinity l(l^2-b^2)^(-s), b=|q|/2.
    The j=2 pole-times-zero contributes -b^4 to Z(-1), and
    b^4(1/4+psi(b+1)/2) inside Z'(-1)/4. It must not be discarded.
    All special values and the harmonic finite part retain exact rationals
    apart from the explicitly evaluated Euler/Hurwitz constants.
    """
    q = _charge(charge)
    if (isinstance(decimal_precision, bool) or not isinstance(decimal_precision, int)
            or decimal_precision < 40 or isinstance(max_terms, bool)
            or not isinstance(max_terms, int) or max_terms < 1):
        raise ValueError("at least 40 decimal digits and a positive finite series budget required")
    tolerance = Fraction(str(series_tolerance))
    if tolerance <= 0:
        raise ValueError("positive spectral-series truncation tolerance required")
    harmonic = sum((Fraction(1, n) for n in range(1, q+1)), Fraction(0))
    with mp.workdps(decimal_precision):
        b = mp.mpf(q)/2; a = b+1
        constants = (2*mp.diff(lambda s: mp.zeta(s, a), -3)
                     +b*b*mp.zeta(-1, a)
                     -2*b*b*mp.diff(lambda s: mp.zeta(s, a), -1)
                     +b**4*(mp.mpf(1)/4+mp.digamma(a)/2))
        partial = mp.mpf(0); last = 2; bound = Fraction(0)
        if q:
            for j in range(3, max_terms+3):
                partial += b**(2*j)*mp.zeta(2*j-3, a)/(j*(j-1))
                last = j; bound = derivative_series_tail_bound(q, j+1)
                if bound <= tolerance:
                    break
            else:
                raise ArithmeticError("new charged spectral sum did not reach its declared tail tolerance")
        derivative = +(4*(constants-partial))
        return MagneticLightSpectrum(
            q, Fraction(1, 30)-Fraction(q*q, 6), -Fraction(q)-Fraction(1, 3),
            harmonic, derivative, +partial, bound, max(0, last-2), last, decimal_precision)


def homogeneous_magnetic_restoration(theta, spectrum):
    """Nonzero-angular local restoration on the owned homogeneous geometry.

    This is an allocation, not the state-minus-reference integral or full
    tensor. The caller separately adds the physical-g2 LLL state/geometric
    source once and the already locked compact complement once. The gauge
    term is the known light Dirac a4 gauge coefficient (2/3)F^2, not a new
    Maxwell coupling. No compact C_F, C_W or Einstein term is included here.
    """
    if not isinstance(spectrum, MagneticLightSpectrum):
        raise ValueError("evaluated q=4 magnetic spectrum or neutral control required")
    q = _charge(spectrum.charge)
    if np.ndim(theta) != 0 or np.iscomplexobj(theta) or not np.isfinite(theta) or not 0 < theta < pi:
        raise ValueError("a real point on the owned homogeneous interior is required")
    theta = float(theta)
    jets = profile_jets(theta); W, W1, W2, _, _ = jets
    if W <= 0:
        raise ValueError("timelike homogeneous interior requires W>0")
    sin, cos = np.sin(theta), np.cos(theta)
    radius = 1/sin; sigma = np.log(radius)
    sigma1, sigma2 = -cos/sin, 1/sin**2
    energy4, parallel4 = curvature_squared_tensor(jets)
    cylinder = float(spectrum.cylinder_density())
    harmonic = float(spectrum.fourth_order_harmonic)
    bar_rho = cylinder+harmonic*energy4/(480*pi*pi)
    bar_parallel = -cylinder-harmonic*parallel4/(480*pi*pi)
    # Reuse only this owner's neutral 4D WZ map and neutral gravitational
    # trace. Its returned sphere pressure is replaced by the charged split.
    neutral_map = physical_source({
        'q': theta, 'W_jets': jets, 'rho_coordinate': -cos/sin,
        'sphere_radius': radius, 'bar_rho': bar_rho, 'bar_p_parallel': bar_parallel,
    })
    gauge_trace = q*q/(48*pi*pi*radius**4)
    gauge_rho = sigma*gauge_trace
    # These are q units of the 2D Weyl cocycle, lifted by 4*pi*r^2.
    # Subtract them here because the inherited LLL uses physical g2.
    lll_weyl_rho = -q*(W1*sigma1+W*sigma1*sigma1)/(96*pi*pi*radius**4)
    lll_weyl_parallel = q*(2*W*sigma2+W1*sigma1-W*sigma1*sigma1)/(96*pi*pi*radius**4)
    # Same A'' curvature as the existing physical radial metric, expressed
    # through its owned q-profile: d/d(rho)=sin(theta)^2 d/d(theta).
    radial_R2 = -(W2+2*W1*sigma1+2*W*sigma2)/radius**2
    physical_lll_trace = q*radial_R2/(96*pi*pi*radius**2)
    rho = neutral_map['rho']+gauge_rho-lll_weyl_rho
    parallel = neutral_map['p_parallel']-gauge_rho-lll_weyl_parallel
    trace = neutral_map['trace']+gauge_trace-physical_lll_trace
    sphere = (rho-parallel-trace)/2
    return {
        'rho': float(rho), 'p_parallel': float(parallel), 'T01': 0.,
        'p_sphere': float(sphere), 'trace': float(trace),
        'theta': theta, 'rho_coordinate': float(-cos/sin), 'radius': float(radius),
        'physical_radial_R2': float(radial_R2),
        'neutral_gravitational_trace': neutral_map['trace'],
        'gauge_trace': float(gauge_trace), 'excluded_physical_LLL_trace': float(physical_lll_trace),
        'bar_static_density': cylinder, 'bar_fourth_order_harmonic': harmonic,
        'bar_rho': float(bar_rho), 'bar_p_parallel': float(bar_parallel),
        'neutral_4D_WZ_rho': neutral_map['local_anomaly_rho'],
        'gauge_WZ': {'rho': float(gauge_rho), 'p_parallel': float(-gauge_rho)},
        'subtracted_LLL_Weyl_cocycle': {'rho': float(lll_weyl_rho), 'p_parallel': float(lll_weyl_parallel)},
        'static_density_series_error_bound': float(spectrum.cylinder_series_error_bound()),
        'normalization_mu': 1., 'magnetic_charge_abs': q,
        'scope': {'nonzero_angular_restoration_only': True, 'LLL_state_or_geometric_source_included': False,
                  'compact_complement_included': False, 'state_integral_included': False,
                  'new_field_or_fitted_coefficient': False, 'full_stress_or_constraints_claimed': False},
    }
