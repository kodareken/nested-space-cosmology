"""First-order four-dimensional curvature EFT map for the current source.

This is a perturbative change of variables, not a fundamental fourth-order
gravity theory. Physical observables, boundaries and the quantum Jacobian
must follow the same map. Full coefficients remain symbolic inputs.
"""
import sympy as sp

from .nsc_spherical_action import canonical_constraints, potentials


def tensor_trace(covariant, inverse_metric):
    return sp.trace(inverse_metric*covariant)


def tensor_square(covariant, inverse_metric):
    return sp.trace(inverse_metric*covariant*inverse_metric*covariant)


def inverse_metric_shift(ricci, stress, inverse_metric, A4, cW, cR):
    """old inverse metric = new inverse metric + this O(cW,cR) term.

    S0 variation is -A4*(G-T/(2A4))*delta g^mu nu; T includes vacuum.
    The original a4 density is -cW*C²-cR*R² in the laboratory convention.
    """
    R=tensor_trace(ricci,inverse_metric)
    T=tensor_trace(stress,inverse_metric)
    ricci_up=inverse_metric*ricci*inverse_metric
    stress_up=inverse_metric*stress*inverse_metric
    return (-2*cW/A4*(ricci_up-R*inverse_metric/6
             +(stress_up-T*inverse_metric/3)/(2*A4))
            +cR/A4*(R-T/(2*A4))*inverse_metric)


def stress_contact(stress, inverse_metric, A4, cW, cR):
    """Bulk contact after order reduction; Euler/boundary terms kept separately."""
    T=tensor_trace(stress,inverse_metric)
    return (-cW/(2*A4**2)*(tensor_square(stress,inverse_metric)-T*T/3)
            -cR*T*T/(4*A4**2))


def spherical_contact(X, field, kinetic, A4, V4, C_F, magnetic_flux, cW, cR):
    """Contact density per sqrt|g2|, in the leading spherical source variables.

    kinetic=(grad field)². Includes the eliminated electric Maxwell energy,
    background magnetic energy and vacuum contribution. Substitution of the
    leading electric solution is valid at this perturbative order only.
    """
    area=X/(2*A4)
    radius_squared=X/(8*sp.pi*A4)
    mu2=potentials(X,A4,V4,C_F,magnetic_flux)["mu_squared"]
    magnetic=C_F*magnetic_flux**2/(2*radius_squared**2)
    electric=mu2*field**2/(2*area)
    b=-cW/(2*A4**2)
    bR=-cR/(4*A4**2)
    return b*(kinetic**2/(2*area)+4*area*(magnetic+electric)**2-4*area*V4**2/3)+16*bR*area*V4**2


def corrected_constraints(p, q, dp, field, field_prime, momentum,
                          A4, V4, C_F, magnetic_flux, cW, cR, order):
    """Canonical constraints through first order in the bookkeeping parameter.

    Terms quadratic in order are not an accepted extension of the EFT.
    The leading velocity is used in the perturbative Legendre transform.
    """
    pars=potentials(p[0],A4,V4,C_F,magnetic_flux)
    base=canonical_constraints(p,q,dp,pars["U"],pars["V"],field,field_prime,momentum,pars["mu_squared"])
    kinetic0=(field_prime**2-momentum**2)/(2*q[1]*q[2])
    correction=spherical_contact(p[0],field,kinetic0,A4,V4,C_F,magnetic_flux,cW,cR)
    return base+order*sp.Matrix([0,q[2]*correction,-q[1]*correction])
