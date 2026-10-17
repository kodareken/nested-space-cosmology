"""Check the Step-7 trace independently against the implemented bulk equations."""
import time

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_effective_vacuum import (
    classical_vacuum_calculation,
    conformal_anomaly_trace_projection,
    conformal_anomaly_scale_si,
    conformal_interface_heat_coefficients,
    constant_curvature_trace_condition,
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


def test_unsourced_cap_curvature_from_independent_si_trace():
    result = constant_curvature_trace_condition()
    assert all(value == 0 for value in result['checks'].values())
    G, hbar, c, H2, ell = sp.symbols('G hbar c H2 ell', positive=True)
    # Work directly in dimensional SI: R=12H^2, T=-hbar*c*H^4/(240*pi^2).
    trace_equation = 12*H2-(8*sp.pi*G/c**4)*hbar*c*H2**2/(240*sp.pi**2)
    si_curved_root = sp.solve(trace_equation, H2)[0]  # positive symbol excludes zero
    assert sp.simplify(si_curved_root*G*hbar/c**3-360*sp.pi) == 0
    zeta = sp.symbols('zeta', positive=True)
    geometric_curved_root = next(root for root in result['unsourced_X_roots'] if root != 0)
    assert sp.simplify(geometric_curved_root.subs(zeta, 8*sp.pi*G*hbar/(c**3*ell**2))
                       -si_curved_root*ell**2) == 0


def test_conformal_interface_matches_gauss_and_reference_value():
    started = time.process_time()
    result = conformal_interface_heat_coefficients()
    assert time.process_time()-started < 60
    assert all(value == 0 for value in result['checks'].values())
    k, H = sp.symbols('k H', positive=True)
    b, h, m, x = sp.symbols('b h m x', positive=True)
    d1, d2, d3 = sp.symbols('d1 d2 d3', real=True)
    assert result['normal_dictionary'] == {
        'U': k, 'L_sum': 3*k, 'S': k, 'D_over_mean_K': -2}
    assert result['a2_integrand'] == 0
    assert result['surface_polynomial'] == -4*k*(k**2+3*H**2)
    assert result['surface_density'] == -k*(k**2+3*H**2)/(1440*sp.pi**2)
    assert result['normalized_density_ell3'] == -b*(b**2+3*h**2)/(1440*sp.pi**2)
    assert result['reference_normalized_density'] == -sp.Rational(7, 720)/sp.pi**2
    # Generic diagonal cancellation, independent of the spherical eigenvalues.
    trace_D, trace_D2 = d1+d2+d3, d1**2+d2**2+d3**2
    generic_A1 = (-trace_D2*(3*k)-(k*trace_D)*trace_D+2*k*trace_D2
                  +k*(trace_D**2+trace_D2))
    assert sp.factor(generic_A1) == 0
    # Actual quiet shell versus Q_plus-Q_minus, with 2M/R^3 on the time axis.
    beta_plus = m/(b*x**2)-(b**2+h**2)*x/(2*b)
    beta_minus = m/(b*x**2)+(b**2-h**2)*x/(2*b)
    D_tau = -sp.diff(beta_plus+beta_minus, x)
    D_angle = -(beta_plus+beta_minus)/x
    trace_wall = D_tau+2*D_angle
    assert sp.factor((b/2)*(trace_wall+D_tau)-(2*m/x**3+2*h**2)) == 0
    assert sp.factor((b/2)*(trace_wall+D_angle)-(-m/x**3+2*h**2)) == 0


def test_interface_coefficient_from_euler_boundary_term():
    # Independent Euler boundary polynomial, Davis hep-th/0208205v3 (6)-(7),
    # converted from inward to outward K. No transmittal heat coefficients enter it.
    result = conformal_interface_heat_coefficients()
    b, h, m, r = sp.symbols('b h m r', positive=True)
    beta_minus = m/(b*r**2)+(b**2-h**2)*r/(2*b)
    beta_plus = m/(b*r**2)-(b**2+h**2)*r/(2*b)
    Kminus = sp.diag(sp.diff(beta_minus, r), beta_minus/r, beta_minus/r)
    Kplus = sp.diag(sp.diff(beta_plus, r), beta_plus/r, beta_plus/r)

    def J(K):
        return (sp.trace(K)*sp.trace(K*K)-sp.Rational(2, 3)*sp.trace(K*K*K)
                -sp.trace(K)**3/3)

    # Unit flat four-ball: outward K=I, induced Einstein tensor=-I.
    ball_density = 4*(J(sp.eye(3))+2*sp.trace(sp.eye(3)))
    assert ball_density*(2*sp.pi**2) == 32*sp.pi**2
    intrinsic_R = 6*h**2+sp.trace(Kminus)**2-sp.trace(Kminus*Kminus)
    euler_surface = sp.factor(4*(J(Kminus)-J(Kplus)+b*intrinsic_R))
    scalar_euler_anomaly = -euler_surface/(360*(4*sp.pi)**2)
    assert sp.factor(scalar_euler_anomaly-result['normalized_density_ell3']) == 0
