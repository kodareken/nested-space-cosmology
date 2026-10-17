"""Exact first Step-7 calculation for the owned conformal scalar action.

This projects the smooth classical stress, not an averaged cosmology or a
quantum effective action. The cap vacuum is an input to this candidate.
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


if __name__ == '__main__':
    calculation = classical_vacuum_calculation()
    for name, value in calculation.items():
        print(f'{name}: {value}')
    print('Scope: smooth classical bulk of the current conformal scalar candidate.')
    print('Supplied baseline at L_ref=ell: cap Lambda0*ell^2=3; exterior=0.')
    print('The zero is the added smooth trace projection, not a cosmological fit.')
