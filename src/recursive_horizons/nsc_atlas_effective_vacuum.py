"""Exact Step-7 trace calculations for the owned conformal scalar candidate.

The classical function projects the smooth classical stress. The anomaly
function adds one real conformal scalar's curvature-trace projection on each
quiet smooth background. The cap vacuum in the classical calculation is an
input. The anomaly coefficient is not that input and is not a backreaction.
"""

import sympy as sp


def classical_vacuum_calculation():
    """Derive the Jordan trace and check it through the Einstein-frame action.

    X is (nabla phi)^2, B is box(phi), and R is the Jordan Ricci scalar.
    In the independent frame check, q=Psi/sqrt(6), Y=(nabla_E q)^2,
    and Lambda0=kappa*V0. All identities hold in a smooth bulk region.
    """
    phi, X, B, R, Lambda0 = sp.symbols('phi X B R Lambda0', real=True)
    kappa = sp.Symbol('kappa', positive=True)
    Lref = sp.Symbol('Lref', positive=True)
    box_square = 2*phi*B+2*X
    # Contract all four terms of the improved Hilbert stress, including G phi^2.
    scalar_trace = sp.expand(X-2*X+(-R*phi**2+3*box_square)/6)
    wave_residual = B-R*phi/6
    trace_identity = sp.expand(scalar_trace-phi*wave_residual)
    einstein_trace = -R+4*Lambda0-kappa*scalar_trace
    R_on_shell = sp.solve(einstein_trace.subs(B, R*phi/6), R)[0]

    q, Y = sp.symbols('q Y', real=True)
    F = 1/sp.cosh(q)**2
    potential = Lambda0/F**2  # kappa V_E; inherited from the same fixed V0.
    log_F = sp.log(F)
    box_q = sp.diff(potential, q)/6
    box_log_F = sp.diff(log_F, q)*box_q+sp.diff(log_F, q, 2)*Y
    grad_log_F_squared = sp.diff(log_F, q)**2*Y
    R_E = 6*Y+4*potential
    R_J_from_E = sp.simplify(F*(R_E+3*box_log_F-sp.Rational(3, 2)*grad_log_F_squared))
    checks = {
        'improved_trace_is_field_times_equation': trace_identity,
        'jordan_trace_fixes_only_supplied_vacuum': sp.simplify(R_on_shell-4*Lambda0),
        'einstein_frame_gives_same_jordan_curvature': sp.simplify(R_J_from_E-R_on_shell),
    }
    # Deliberately omit the improvement: a kinetic scalar is not trace-free.
    omitted_improvement = -X
    if any(value != 0 for value in checks.values()) or omitted_improvement == 0:
        raise ArithmeticError('conformal vacuum projection identity failed')
    return {
        'scalar_trace': scalar_trace,
        'wave_residual': wave_residual,
        'jordan_ricci': R_on_shell,
        'jordan_ricci_from_einstein_frame': R_J_from_E,
        'induced_lambda_Lref_squared': sp.simplify(Lref**2*(R_on_shell-4*Lambda0)/4),
        'checks': checks,
        'omitted_improvement_control': omitted_improvement,
    }


def _coordinate_spherical_invariants(f_expr):
    """Contract Christoffel symbols of ds^2=-f dt^2+dr^2/f+r^2 dOmega^2.

    Coordinates are dimensionless in the supplied reference length ell. The
    calculation is local on a smooth chart and excludes wall deltas.
    """
    r, theta = sp.symbols('r theta', positive=True)
    coords = (sp.Symbol('t', real=True), r, theta, sp.Symbol('phi', real=True))
    f = sp.sympify(f_expr)
    diagonal = (-f, 1/f, r**2, r**2*sp.sin(theta)**2)
    inverse = tuple(sp.together(1/component) for component in diagonal)
    gamma = [[[0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for rho in range(4):
        for mu in range(4):
            for nu in range(mu, 4):
                bracket = 0
                if rho == mu:
                    bracket += sp.diff(diagonal[rho], coords[nu])
                if rho == nu:
                    bracket += sp.diff(diagonal[rho], coords[mu])
                if mu == nu:
                    bracket -= sp.diff(diagonal[mu], coords[rho])
                if bracket == 0:
                    continue
                value = sp.together(inverse[rho]*bracket/2)
                gamma[rho][mu][nu] = value
                gamma[rho][nu][mu] = value

    riemann = {}
    for rho in range(4):
        for sigma in range(4):
            for mu in range(4):
                for nu in range(mu+1, 4):
                    curvature = (sp.diff(gamma[rho][nu][sigma], coords[mu])
                                 - sp.diff(gamma[rho][mu][sigma], coords[nu]))
                    for lam in range(4):
                        curvature += (gamma[rho][mu][lam]*gamma[lam][nu][sigma]
                                      - gamma[rho][nu][lam]*gamma[lam][mu][sigma])
                    curvature = sp.together(curvature)
                    if curvature == 0:
                        continue
                    riemann[(rho, sigma, mu, nu)] = curvature
                    riemann[(rho, sigma, nu, mu)] = -curvature

    def riemann_up(rho, sigma, mu, nu):
        if mu == nu:
            return 0
        if mu > nu:
            return -riemann.get((rho, sigma, nu, mu), 0)
        return riemann.get((rho, sigma, mu, nu), 0)

    ricci = [[0 for _ in range(4)] for _ in range(4)]
    for sigma in range(4):
        for nu in range(sigma, 4):
            value = sp.together(sum(riemann_up(rho, sigma, rho, nu) for rho in range(4)))
            ricci[sigma][nu] = value
            ricci[nu][sigma] = value
    ricci_scalar = sp.together(sum(inverse[i]*ricci[i][i] for i in range(4)))
    ricci_square = sp.together(sum(
        inverse[mu]*inverse[nu]*ricci[mu][nu]**2
        for mu in range(4) for nu in range(4)))
    riemann_square = sp.together(sum(
        diagonal[rho]*inverse[sigma]*inverse[mu]*inverse[nu]*curvature**2
        for (rho, sigma, mu, nu), curvature in riemann.items()))
    return {
        'ricci_scalar': sp.factor(ricci_scalar),
        'ricci_square': sp.factor(ricci_square),
        'riemann_square': sp.factor(riemann_square),
    }


def _analytic_spherical_invariants(f_expr):
    """Reduced invariants of the same spherical chart.

    Orthonormal sectional components of a smooth f(r) chart are
    R_0101=f''/2, R_0202=R_0303=f'/(2r), R_1212=R_1313=-f'/(2r) and
    R_2323=(1-f)/r^2. Their contractions are the independent check of the
    coordinate Christoffel sum. No third metric is expanded.
    """
    r = sp.symbols('r', positive=True)
    f = sp.sympify(f_expr)
    df, ddf = sp.diff(f, r), sp.diff(f, r, 2)
    radial, mixed, sphere = ddf/2, df/(2*r), (1-f)/r**2
    time_eigenvalue = radial+2*mixed
    angle_eigenvalue = sphere-2*mixed
    return {
        'ricci_scalar': sp.factor(-ddf-4*df/r+2*(1-f)/r**2),
        'ricci_square': sp.factor(2*time_eigenvalue**2+2*angle_eigenvalue**2),
        'riemann_square': sp.factor(4*radial**2+16*mixed**2+4*sphere**2),
    }


def conformal_anomaly_trace_projection():
    """Trace projection of one real conformal scalar on each quiet smooth background.

    Free-field quantum theory and local covariant renormalization are imported.
    This is not a derived quantum mechanics, not a joint state, not the full
    vacuum stress, and not an H induced by nesting. The classical improved
    scalar remains trace-free. The general Laplace-type coefficient below is
    imported, not proved here; only its R^2 cancellation at E=-R/6 is derived.
    Finite-volume and Newton counterterms stay supplied renormalized
    coefficients. Lref is the existing ell, h, m and zeta are supplied, and
    there is no self-consistent backreaction.
    """
    R, ricci2, riemann2, box_R, box_E, E = sp.symbols(
        'R Ricci2 Riemann2 BoxR BoxE E')
    # Imported scalar Seeley-DeWitt coefficient, before the (4*pi)^-2 factor.
    a4_numerator = (5*R**2-2*ricci2+2*riemann2+60*R*E+180*E**2
                    +12*box_R+60*box_E)
    conformal = {E: -R/6, box_E: -box_R/6}
    a4_conformal = sp.expand(sp.together(a4_numerator.subs(conformal))/360)
    duff_numerator = riemann2-ricci2+box_R
    heat_identity = sp.expand(a4_conformal-duff_numerator/180)
    r_square_coefficient = a4_conformal.coeff(R**2)
    minimal_r_square = sp.expand(a4_numerator.subs({E: 0, box_E: 0})/360).coeff(R**2)

    C2 = riemann2-2*ricci2+sp.Rational(1, 3)*R**2
    euler = riemann2-4*ricci2+R**2
    c_coefficient, a_coefficient = sp.Rational(1, 120), sp.Rational(1, 360)
    weyl_euler = sp.expand(c_coefficient*C2-a_coefficient*euler-(riemann2-ricci2)/180)
    normalization_4pi = sp.expand(180*(4*sp.pi)**2-2880*sp.pi**2)

    h, m, radius = sp.symbols('h m r', positive=True)
    zeta = sp.symbols('zeta', positive=True)
    quiet = {
        'cap': _coordinate_spherical_invariants(1-h**2*radius**2),
        'exterior': _coordinate_spherical_invariants(1-2*m/radius),
    }
    analytic = {
        'cap': _analytic_spherical_invariants(1-h**2*radius**2),
        'exterior': _analytic_spherical_invariants(1-2*m/radius),
    }
    expected = {
        'cap': {'ricci_scalar': 12*h**2, 'ricci_square': 36*h**4,
                'riemann_square': 24*h**4},
        'exterior': {'ricci_scalar': 0, 'ricci_square': 0,
                     'riemann_square': 48*m**2/radius**6},
    }
    checks = {
        'conformal_a4_cancels_R2': heat_identity,
        'conformal_R2_coefficient': r_square_coefficient,
        'weyl_euler_matches_riemann_minus_ricci': weyl_euler,
        'duff_denominator_is_2880_pi2': normalization_4pi,
    }
    for region in ('cap', 'exterior'):
        for name, target in expected[region].items():
            checks[region+'_'+name] = sp.factor(quiet[region][name]-target)
            checks[region+'_analytic_'+name] = sp.factor(analytic[region][name]-target)
        checks[region+'_constant_R'] = sp.diff(quiet[region]['ricci_scalar'], radius)

    def coefficient(region):
        curvature = quiet[region]['riemann_square']-quiet[region]['ricci_square']
        # Box R vanishes because R is constant on each quiet smooth region.
        return sp.factor(-curvature/(4*2880*sp.pi**2))

    cap_coefficient = coefficient('cap')
    exterior_coefficient = coefficient('exterior')
    cap_reference = sp.factor(cap_coefficient.subs(h, 1))
    exterior_reference = sp.factor(exterior_coefficient.subs(
        {m: sp.Rational(2, 5), radius: sp.Rational(11, 20)}))
    checks['cap_coefficient'] = sp.factor(cap_coefficient-h**4/(960*sp.pi**2))
    checks['exterior_coefficient'] = sp.factor(
        exterior_coefficient+m**2/(240*sp.pi**2*radius**6))
    checks['zero_curvature_cap'] = cap_coefficient.subs(h, 0)
    checks['zero_curvature_exterior'] = sp.factor(exterior_coefficient.subs(m, 0))
    if any(value != 0 for value in checks.values()):
        failed = [name for name, value in checks.items() if value != 0]
        raise ArithmeticError('conformal anomaly trace projection failed: '+', '.join(failed))
    exterior_radial_variation = sp.factor(sp.diff(exterior_coefficient, radius))
    if exterior_radial_variation == 0 or minimal_r_square == 0:
        raise ArithmeticError('anomaly control vanished')
    return {
        'imported_a4': a4_numerator/360,
        'conformal_endomorphism': -R/6,
        'conformal_a4': a4_conformal,
        'weyl_c': c_coefficient,
        'euler_a': a_coefficient,
        'cap_invariants': quiet['cap'],
        'exterior_invariants': quiet['exterior'],
        'cap_coefficient': cap_coefficient,
        'exterior_coefficient': exterior_coefficient,
        'cap_delta_lambda_ell_squared': sp.factor(zeta*cap_coefficient),
        'exterior_delta_lambda_ell_squared': sp.factor(zeta*exterior_coefficient),
        'cap_reference_coefficient': cap_reference,
        'exterior_reference_coefficient': exterior_reference,
        'checks': checks,
        'minimal_scalar_R2_coefficient': minimal_r_square,
        'exterior_radial_variation': exterior_radial_variation,
    }


def conformal_anomaly_scale_si(ell_m=10000):
    """Evaluate the quiet trace at P001's supplied shape and a fiducial SI scale.

    CODATA 2022 inputs: G=6.67430(15)e-11, h=6.62607015e-34 exactly,
    c=299792458 exactly. Source: https://physics.nist.gov/cuu/pdf/wall_2022.pdf
    h is Planck's constant here, distinct from the shape parameter H*ell=1.
    The returned uncertainty is solely the standard input uncertainty in G;
    it is not a quantum-model, wall, backreaction or observational error bar.
    """
    length = sp.Rational(str(ell_m))
    if length <= 0:
        raise ValueError('ell_m must be a positive finite length in metres')
    G, u_G = sp.Rational('6.67430e-11'), sp.Rational('0.00015e-11')
    planck_h, c = sp.Rational('6.62607015e-34'), sp.Integer(299792458)
    hbar = planck_h/(2*sp.pi)
    # kappa_SI=8*pi*G/c^4 multiplies trace energy density hbar*c*curvature^2.
    zeta = sp.factor(8*sp.pi*G*hbar/(c**3*length**2))
    symbolic = conformal_anomaly_trace_projection()
    regions = {}
    for region in ('cap', 'exterior'):
        value = sp.factor(zeta*symbolic[region+'_reference_coefficient'])
        regions[region] = {
            'delta_lambda_ell_squared': value,
            'delta_lambda_per_m2': value/length**2,
            'standard_uncertainty_from_G': sp.Abs(value)*u_G/G,
        }
    return {
        'ell_m': length,
        'G_SI': G,
        'G_standard_uncertainty_SI': u_G,
        'planck_h_SI_exact': planck_h,
        'c_SI_exact': c,
        'zeta': zeta,
        'relative_standard_uncertainty_from_G': u_G/G,
        'regions': regions,
    }


if __name__ == '__main__':
    calculation = classical_vacuum_calculation()
    for name, value in calculation.items():
        print(f'{name}: {value}')
    print('Scope: smooth classical bulk of the current conformal scalar candidate.')
    print('Supplied baseline at L_ref=ell: cap Lambda0*ell^2=3; exterior=0.')
    print('The zero is the added smooth trace projection, not a cosmological fit.')
    anomaly = conformal_anomaly_trace_projection()
    print('cap DeltaLambda_trace*ell^2:', anomaly['cap_delta_lambda_ell_squared'])
    print('exterior DeltaLambda_trace*ell^2:', anomaly['exterior_delta_lambda_ell_squared'])
    print('cap coefficient at h=1:', anomaly['cap_reference_coefficient'])
    print('exterior coefficient at m=2/5, r=11/20:', anomaly['exterior_reference_coefficient'])
    print('cap numeric coefficient:', anomaly['cap_reference_coefficient'].evalf(20))
    print('exterior numeric coefficient:', anomaly['exterior_reference_coefficient'].evalf(20))
    print('Scope: one real conformal scalar; trace projection on each smooth quiet chart.')
    print('Box R is the scheme-dependent piece and vanishes because R is constant.')
    print('The exterior projection depends on r. The cap trace does not fix the full stress.')
    scaled = conformal_anomaly_scale_si()
    print('Fiducial ell in metres:', scaled['ell_m'])
    print('zeta from CODATA 2022 inputs:', scaled['zeta'].evalf(12))
    for region, row in scaled['regions'].items():
        print(region, 'DeltaLambda_trace*ell^2:', row['delta_lambda_ell_squared'].evalf(12),
              '+/-', row['standard_uncertainty_from_G'].evalf(5), '(G input only)')
    print('This scale illustration is conditional on ell=10000 m; it is not a cosmological prediction.')
