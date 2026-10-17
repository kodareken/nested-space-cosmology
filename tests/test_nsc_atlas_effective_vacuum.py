"""Check the Step-7 trace independently against the implemented bulk equations."""
import time

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_effective_vacuum import (
    classical_vacuum_calculation,
    conformal_anomaly_trace_projection,
    conformal_anomaly_scale_si,
)
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator


def test_exact_jordan_and_einstein_action_calculation():
    result = classical_vacuum_calculation()
    assert all(value == 0 for value in result['checks'].values())
    assert result['induced_lambda_Lref_squared'] == 0
    assert result['omitted_improvement_control'] != 0


@pytest.mark.parametrize('region,expected', [('cap', 12.), ('exterior', 0.)])
@pytest.mark.parametrize('alpha', [0., 0.00012983127142551443, .02])
def test_geometric_ricci_from_actual_mixed_operator(region, expected, alpha):
    # Independent point jets, not an evolved state or a PDE existence claim.
    r, C, c = .49, 1.3, .21
    rU, rV, cU, cV = -.37, .42, -.53, .19
    row = bulk_operator(dict(r=r, C=C, c=c),
                        dict(r_U=rU, r_V=rV, c_U=cU, c_V=cV),
                        alpha=alpha, region=region)
    rUV, cUV = row['r_UV'], row['c_UV']
    # Ricci of -C dU dV+r^2 dOmega^2, computed from its metric jets.
    R_E = 4*row['log_C_UV']/C+2/r**2+8*rU*rV/(C*r**2)+16*rUV/(C*r)
    u = c/r
    uU, uV = (cU-u*rU)/r, (cV-u*rV)/r
    uUV = cUV/r-(cU*rV+cV*rU+c*rUV)/r**2+2*c*rU*rV/r**3
    D = 1+alpha*u*u/6
    lU, lV = -alpha*u*uU/(3*D), -alpha*u*uV/(3*D)
    lUV = -alpha*(uU*uV+u*uUV)/(3*D)+alpha**2*u*u*uU*uV/(9*D**2)
    box_l = -4*(lUV+(rU*lV+rV*lU)/r)/C
    grad_l_squared = -4*lU*lV/C
    R_J = (R_E+3*box_l-1.5*grad_l_squared)/D
    assert np.isfinite(R_J)
    assert R_J == pytest.approx(expected, rel=0, abs=2e-14)


def test_conformal_anomaly_trace_coefficients():
    started = time.process_time()
    result = conformal_anomaly_trace_projection()
    assert time.process_time()-started < 60
    assert all(value == 0 for value in result['checks'].values())
    assert result['minimal_scalar_R2_coefficient'] != 0
    assert result['exterior_radial_variation'] != 0
    h, m, radius, zeta = sp.symbols('h m r zeta', positive=True)
    assert result['cap_delta_lambda_ell_squared'] == zeta*h**4/(960*sp.pi**2)
    assert result['exterior_delta_lambda_ell_squared'] == -zeta*m**2/(240*sp.pi**2*radius**6)
    assert result['cap_reference_coefficient'] == 1/(960*sp.pi**2)
    assert result['exterior_reference_coefficient'] == -sp.Integer(128000)/(5314683*sp.pi**2)
    assert result['weyl_c'] == sp.Rational(1, 120)
    assert result['euler_a'] == sp.Rational(1, 360)


def test_zero_curvature_limit_of_trace_projection():
    result = conformal_anomaly_trace_projection()
    h, m = sp.symbols('h m', positive=True)
    assert result['cap_coefficient'].subs(h, 0) == 0
    assert sp.factor(result['exterior_coefficient'].subs(m, 0)) == 0


def _radial_invariants(f_expr):
    """Independent sectional reduction, separate from the Christoffel sum."""
    radius = sp.symbols('r', positive=True)
    f = sp.sympify(f_expr)
    df, ddf = sp.diff(f, radius), sp.diff(f, radius, 2)
    radial, mixed, sphere = ddf/2, df/(2*radius), (1-f)/radius**2
    return {
        'ricci_scalar': sp.factor(-ddf-4*df/radius+2*(1-f)/radius**2),
        'ricci_square': sp.factor(2*(radial+2*mixed)**2+2*(sphere-2*mixed)**2),
        'riemann_square': sp.factor(4*radial**2+16*mixed**2+4*sphere**2),
    }


def test_curvature_contraction_matches_radial_reduction():
    result = conformal_anomaly_trace_projection()
    h, m, radius = sp.symbols('h m r', positive=True)
    for f_expr, invariants in (
            (1-h**2*radius**2, result['cap_invariants']),
            (1-2*m/radius, result['exterior_invariants'])):
        reduced = _radial_invariants(f_expr)
        for name in ('ricci_scalar', 'ricci_square', 'riemann_square'):
            assert sp.factor(invariants[name]-reduced[name]) == 0
    assert result['cap_invariants']['riemann_square'] == 24*h**4
    assert result['cap_invariants']['ricci_square'] == 36*h**4
    assert result['cap_invariants']['ricci_scalar'] == 12*h**2
    assert result['exterior_invariants']['riemann_square'] == 48*m**2/radius**6
    assert result['exterior_invariants']['ricci_square'] == 0
    assert result['exterior_invariants']['ricci_scalar'] == 0


def test_si_scale_from_energy_density_and_reference_length():
    result = conformal_anomaly_scale_si()
    ell, G, h, c = (result[key] for key in
                    ('ell_m', 'G_SI', 'planck_h_SI_exact', 'c_SI_exact'))
    # Independent SI trace: -hbar*c*H^4/(240*pi^2); project with -kappa_SI/4.
    trace = -h*c/(480*sp.pi**3*ell**4)
    projected = -(8*sp.pi*G/c**4)*trace/4
    cap = result['regions']['cap']
    assert sp.simplify(cap['delta_lambda_per_m2']-projected) == 0
    assert sp.simplify(cap['standard_uncertainty_from_G']
                       -cap['delta_lambda_ell_squared']*result['G_standard_uncertainty_SI']/G) == 0
    rescaled = conformal_anomaly_scale_si(20000)
    for region in ('cap', 'exterior'):
        for name, factor in [('delta_lambda_ell_squared', 4), ('delta_lambda_per_m2', 16)]:
            assert sp.simplify(factor*rescaled['regions'][region][name]
                               -result['regions'][region][name]) == 0
