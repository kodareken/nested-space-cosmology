"""Independent number-state perturbation check of the two-current exchange."""
import sympy as sp


def test_number_state_routes_recover_the_connected_moment_and_kernel():
    a, b, c, d, y = sp.symbols("a b c d y", positive=True)
    L, R = sp.symbols("L R", real=True)
    EL, ER, Et = a+b+y, c+d+y, a+b+c+d
    product = a*b*c*d
    # H_L|0> = i A_L |1,1,0,0,1>; likewise for H_R.
    AL = L*(a-b)/(2*sp.sqrt(2*a*b*y))
    AR = R*(c-d)/(2*sp.sqrt(2*c*d*y))
    # O=Q1Q2Q3Q4 connects vacuum to |1,1,1,1,0> and the two first-order states.
    Oelement = 1/(4*sp.sqrt(product))
    vacuum_second_order = -2*AL*AR*(1/EL+1/ER)*Oelement/Et
    first_order_cross = 2*AL*AR*Oelement/(EL*ER)
    mixed_moment = sp.simplify(vacuum_second_order+first_order_cross)
    Gamma = -L*R*(a-b)*(c-d)/(Et*EL*ER)
    assert sp.simplify(mixed_moment-Gamma/(8*product)) == 0
    # Divide by sqrt(product of the four free variances), each 1/(2 omega).
    normalized = mixed_moment/(1/(4*sp.sqrt(product)))
    assert sp.simplify(normalized-Gamma/(2*sp.sqrt(product))) == 0
