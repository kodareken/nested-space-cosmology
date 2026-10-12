"""Same-action charged light restoration in an explicit canonical KS foliation.

This includes the LLL GEOMETRIC allocation once, but no LLL occupation/state
term and no compact induced action. Curvature is supplied entirely by the
existing spherical_invariants owner. Three-point quadrature integrates the
formal Weyl interpolation parameter, never a physical time or new field.

The LLL normal-order allocation depends on the canonical spatial foliation.
Input coordinates must therefore be KS (T,z), or a reparameterization of KS
time with the same slices. The PG fields of OwnedCompactKSHistory cannot be
passed silently to this API when the source is canonical in KS coordinates.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_magnetic_light_reference import MagneticLightSpectrum
from .nsc_spherical_local_history import spherical_invariants


CHANNELS = ("cylinder", "bar_radial_R2_squared", "LLL_geometry_bar",
            "WZ_Weyl", "WZ_Euler", "WZ_boxR", "WZ_gauge")


def _validated_data(grid, g, radius, jets):
    shape = (len(grid.time), len(grid.radius_coordinate))
    g, r = np.asarray(g), np.asarray(radius)
    if (g.shape != (*shape, 2, 2) or r.shape != shape or jets is None
            or len(jets) != 4):
        raise ValueError("metric, radius and analytic first/second jets on the canonical grid are required")
    dg, ddg, dr, ddr = map(np.asarray, jets)
    expected = ((*shape, 2, 2, 2), (*shape, 2, 2, 2, 2), (*shape, 2), (*shape, 2, 2))
    if any(v.shape != e or not np.isfinite(v).all() for v, e in zip((dg, ddg, dr, ddr), expected)):
        raise ValueError("finite analytic metric/radius derivative arrays through order two required")
    if (not np.isfinite(g).all() or not np.isfinite(r).all()
            or np.any(r.real <= 0) or np.any(g[..., 1, 1].real >= 0)):
        raise ValueError("positive radius and spacelike canonical spatial slices are required")
    errors = (g-g.swapaxes(-1, -2), dg-dg.swapaxes(-1, -2),
              ddg-ddg.swapaxes(-1, -2), ddg-ddg.swapaxes(-4, -3), ddr-ddr.swapaxes(-1, -2))
    if max(float(np.max(abs(e))) for e in errors) > 3e-11:
        raise ValueError("symmetric metric jets with commuting coordinate derivatives required")
    return g, r, (dg, ddg, dr, ddr)


def radius_power_jets(radius, dr, ddr, exponent):
    """Analytic scalar power through second order; no coordinate differences."""
    r = np.asarray(radius)
    value = r**exponent
    first = exponent*(r**(exponent-1))[..., None]*dr
    second = (exponent*(r**(exponent-1))[..., None, None]*ddr
              +exponent*(exponent-1)*(r**(exponent-2))[..., None, None]
              *dr[..., :, None]*dr[..., None, :])
    return value, first, second


def conformal_metric_jets(g, radius, jets, parameter):
    """g_s2=r^(2s-2) g2 and r_s=r^s, for the formal Weyl path.

    s=0 gives bar g2 and unit sphere; s=1 gives the supplied physical metric.
    Its value is an interpolation parameter, not a history coordinate.
    """
    if np.ndim(parameter) != 0 or not np.isfinite(parameter) or np.iscomplexobj(parameter):
        raise ValueError("finite real formal Weyl parameter required")
    dg, ddg, dr, ddr = jets
    factor, df, ddf = radius_power_jets(radius, dr, ddr, 2*parameter-2)
    G = factor[..., None, None]*g
    dG = df[..., :, None, None]*g[..., None, :, :]+factor[..., None, None, None]*dg
    ddG = (ddf[..., :, :, None, None]*g[..., None, None, :, :]
           +df[..., :, None, None, None]*dg[..., None, :, :, :]
           +df[..., None, :, None, None]*dg[..., :, None, :, :]
           +factor[..., None, None, None, None]*ddg)
    rs, drs, ddrs = radius_power_jets(radius, dr, ddr, parameter)
    return G, rs, (dG, ddG, drs, ddrs)


def weyl_euler_average(grid, g, radius, jets, *, quadrature_order=3):
    """Integrate r^(4s) E4[g_s] over formal s in [0,1].

    The weighted four-dimensional Euler density is a polynomial of degree
    at most four in s; three Gauss points suffice. A different numerical
    order is available only for an independent quadrature identity check.
    """
    if isinstance(quadrature_order, bool) or not isinstance(quadrature_order, int) or quadrature_order < 3:
        raise ValueError("at least three formal Weyl quadrature nodes required")
    x, weights = leggauss(quadrature_order)
    result = np.zeros_like(radius, dtype=np.result_type(g, radius, float))
    for s, w in zip((x+1)/2, weights/2):
        G, rs, sjets = conformal_metric_jets(g, radius, jets, s)
        invariants = spherical_invariants(grid, G, rs, sjets)
        result += w*radius**(4*s)*invariants['E4']
    return result


def _barred_lll_density(g, dg, charge):
    """Existing ADM geometric normal-order action, in the input foliation."""
    a = np.sqrt(-g[..., 1, 1])
    beta = -g[..., 0, 1]/a**2
    N = np.sqrt(g[..., 0, 0]+a*a*beta*beta)
    da = -dg[..., :, 1, 1]/(2*a[..., None])
    dbeta = -dg[..., :, 0, 1]/a[..., None]**2-2*beta[..., None]*da/a[..., None]
    dN = (dg[..., :, 0, 0]+2*(a*beta*beta)[..., None]*da
          +2*(a*a*beta)[..., None]*dbeta)/(2*N[..., None])
    H = (da[..., 0]-beta*da[..., 1]-a*dbeta[..., 1])/(N*a)
    Lz = da[..., 1]/a
    return charge/(24*np.pi)*(-N*a*H*H+2*dN[..., 1]*Lz/a-N*Lz*Lz/a)


@dataclass(frozen=True)
class LightRestorationAction:
    """The retained charged light restoration; coefficients are never fitted."""
    spectrum: MagneticLightSpectrum

    def __post_init__(self):
        if not isinstance(self.spectrum, MagneticLightSpectrum) or self.spectrum.charge not in (0, 4):
            raise ValueError("the owned q=4 spectrum or exact q=0 normalization control is required")

    def actions(self, grid, g, radius, jets, *, canonical_coordinates):
        """Return channel actions and diagnostics from supplied canonical jets.

        ``canonical_coordinates`` must explicitly be 'KS' or
        'KS_time_reparametrization'. This is the caller's declaration that
        the input uses the SAME spatial foliation as its canonical C0.
        PG-to-KS conversion must be done by the existing chart owner first.
        Complex increments are permitted only for analytic differentiation.
        """
        if canonical_coordinates not in ("KS", "KS_time_reparametrization"):
            raise ValueError("same canonical KS spatial foliation required; PG input is not silently converted")
        g, r, jets = _validated_data(grid, g, radius, jets)
        physical = spherical_invariants(grid, g, r, jets)
        bar_g, bar_r, bar_jets = conformal_metric_jets(g, r, jets, 0.)
        barred = spherical_invariants(grid, bar_g, bar_r, bar_jets)
        sqrtbar = barred['measure']  # Unit barred sphere, so this is sqrt(-det bar g2).
        sigma = np.log(r)
        charge = self.spectrum.charge
        Cq = float(self.spectrum.cylinder_density())
        hq = float(self.spectrum.fourth_order_harmonic)
        alpha, beta, gamma = -1/20, 11/360, -1/30
        euler_average = weyl_euler_average(grid, g, r, jets)
        prefactor = -4*np.pi/(16*np.pi**2)
        density = {
            'cylinder': -4*np.pi*Cq*sqrtbar,
            'bar_radial_R2_squared': hq/(240*np.pi)*sqrtbar*barred['R2']**2,
            'LLL_geometry_bar': _barred_lll_density(bar_g, bar_jets[0], charge),
            'WZ_Weyl': prefactor*sqrtbar*alpha*sigma*barred['C2'],
            'WZ_Euler': prefactor*sqrtbar*beta*sigma*euler_average,
            'WZ_boxR': prefactor*sqrtbar*gamma*(barred['R']**2-r**4*physical['R']**2)/12,
            'WZ_gauge': prefactor*sqrtbar*(2/3)*sigma*(charge**2/2),
        }
        actions = {name: grid.integral(density[name]) for name in CHANNELS}
        return actions, {
            'total_action': sum(actions.values()), 'densities': density,
            'physical_invariants': physical, 'barred_invariants': barred,
            'formal_weyl_euler_average': euler_average,
            'canonical_coordinates': canonical_coordinates,
            'formal_weyl_quadrature_order': 3, 'normalization_mu': 1.,
            'scope': {'LLL_geometry_included_once': True, 'LLL_state_included': False,
                      'nonzero_angular_state_integral_included': False, 'compact_induced_action_included': False,
                      'same_C0_spatial_foliation_required': True, 'physical_history_selected': False,
                      'Cauchy_data_or_constraints_solved': False,
                      'boundary': 'Euler derivatives require the declared metric-variation boundary terms; compact checks fix them'},
        }
