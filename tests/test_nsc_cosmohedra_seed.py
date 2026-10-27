"""Formal low-point cubic seed, not an identification of the NSC spacetime."""
import sympy as sp


def test_mixed_connected_seed_directly_satisfies_the_schrodinger_coefficient():
    q = sp.symbols("q1 q2 q3 q4 s")
    w = sp.symbols("w1 w2 w3 w4 y", positive=True)
    gL, gR = sp.symbols("gL gR")
    EL, ER, Et = w[0]+w[1]+w[4], w[2]+w[3]+w[4], sum(w[:4])
    S0 = sum(wi*qi**2/2 for wi, qi in zip(w, q))
    SL = gL*q[0]*q[1]*q[4]/EL
    SR = gR*q[2]*q[3]*q[4]/ER
    SLR = -gL*gR*q[0]*q[1]*q[2]*q[3]/(Et*EL*ER)
    S = S0+SL+SR+SLR
    V = sum(wi**2*qi**2/2 for wi, qi in zip(w, q))
    V += gL*q[0]*q[1]*q[4]+gR*q[2]*q[3]*q[4]
    residual = V+sum((sp.diff(S, qi, 2)-sp.diff(S, qi)**2)/2 for qi in q)-sum(w)/2
    mixed = sp.diff(residual, gL, gR).subs({gL: 0, gR: 0})
    assert sp.simplify(mixed) == 0


def test_seed_and_formal_hamiltonian_have_the_declared_normalized_scaling():
    w = sp.symbols("w1 w2 w3 w4 y", positive=True)
    Omega = sp.symbols("Omega", positive=True)
    g, q1, q2, q3 = sp.symbols("g q1 q2 q3")
    K = 1/(sum(w[:4])*(w[0]+w[1]+w[4])*(w[2]+w[3]+w[4]))
    shifted = K.subs({wi: Omega*wi for wi in w}, simultaneous=True)
    assert sp.simplify(Omega**3*shifted-K) == 0
    # Canonical q,p scaling multiplies kinetic/quadratic H by Omega.
    p, frequency = sp.symbols("p frequency", positive=True)
    quadratic_shifted = (Omega*p**2+(Omega*frequency)**2*q1**2/Omega)/2
    assert sp.simplify(quadratic_shifted/Omega-(p**2+frequency**2*q1**2)/2) == 0
    cubic_shifted = Omega**sp.Rational(5, 2)*g*q1*q2*q3/Omega**sp.Rational(3, 2)
    assert sp.simplify(cubic_shifted/Omega-g*q1*q2*q3) == 0
    # Mixed log-wavefunction term includes two such couplings and four q's.
    assert sp.simplify(Omega**5*shifted/Omega**2-K) == 0
    assert K.subs(dict(zip(w, [1, 2, 3, 4, 2]))) == sp.Rational(1, 450)
