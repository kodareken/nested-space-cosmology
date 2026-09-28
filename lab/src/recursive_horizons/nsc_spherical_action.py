"""First-order owner for the local spherical metric/source interface.

This reuses generalized 2D dilaton gravity. X=8*pi*A4*r² is area, not a
new scalar or the common resolution compensator. A4,V4,C_F denote full
effective coefficients; the recorded Dirac pieces do not fix their totals.
Higher-curvature/nonlocal terms are outside this two-derivative owner.
"""
import sympy as sp


def potentials(X, A4, V4, C_F, magnetic_flux):
    """Potentials in GKV's conventions; physical bulk action is -L_dil.

    L_lab=sqrt|g2|[-X R_L/2+U(X)(grad X)²/2-V(X)]. The chapter-7
    canonical component convention is used for the constraints below.
    """
    U = -1/(2*X)
    V = -8*sp.pi*A4+V4*X/(2*A4)+16*sp.pi**2*A4*C_F*magnetic_flux**2/X
    return {"U": U, "V": V,
            "mu_squared": A4*sp.Abs(magnetic_flux)/(2*sp.pi*C_F*X)}


def poisson_tensor(X, X_plus, X_minus, U, V):
    """Sign convention: G_i=partial_1 p_i+P^{ij}q_j, p=(X,X+,X-)."""
    total = U*X_plus*X_minus+V
    return sp.Matrix([[0, -X_plus, X_minus],
                      [X_plus, 0, -total], [-X_minus, total, 0]])


def canonical_constraints(p, q, dp, U, V, field, field_prime,
                          field_momentum, mu_squared):
    """Coupled local constraints, including the charge-mode area source.

    q=(omega_1,e^-_1,e^+_1), p=(X,X+,X-). The canonical scalar has
    positive Lorentzian kinetic energy. q2*q3<0 gives spacelike slices.
    Its X-dependent potential is retained in both geometry constraints.
    """
    geometry = sp.Matrix(dp)+poisson_tensor(*p,U,V)*sp.Matrix(q)
    energy = mu_squared*field**2/2
    matter = sp.Matrix([0,
        (field_prime-field_momentum)**2/(4*q[1])-q[2]*energy,
        -(field_prime+field_momentum)**2/(4*q[2])+q[1]*energy])
    return geometry+matter


def primitive(X, A4, V4, C_F, magnetic_flux):
    """w'=exp(Q)V, Q=-log(X)/2; integration constant belongs to Casimir."""
    return (-16*sp.pi*A4*sp.sqrt(X)+V4*X**sp.Rational(3,2)/(3*A4)
            -32*sp.pi**2*A4*C_F*magnetic_flux**2/sp.sqrt(X))


def casimir(X, X_plus, X_minus, A4, V4, C_F, magnetic_flux):
    return X_plus*X_minus/sp.sqrt(X)+primitive(X,A4,V4,C_F,magnetic_flux)


def pg_null_coframe(beta):
    """Rows (e^-,e^+), columns (d tau,d rho); determinant is +1."""
    return sp.Matrix([[1-beta,-1],[1+beta,1]])/sp.sqrt(2)
