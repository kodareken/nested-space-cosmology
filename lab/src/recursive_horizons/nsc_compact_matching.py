"""Match the computed compact modulus to one canonical light Dirac field.

These are the finite Dirac contributions in the stated proper-time scheme.
They do not supply independent finite terms, a quantum state or a complete
compensating action. The matching cutoff nu is not a normalization mass.
"""
from math import pi, sqrt

import numpy as np
from scipy.integrate import quad_vec
from scipy.special import erf, erfc, exp1

from .nsc_warped_source import integrated_exp1


def proper_ratio(weight):
    """J=ell_proper/ell for the existing Gaussian warp."""
    a = weight.path*weight.amplitude
    return 1. if a == 0 else sqrt(pi)*erf(sqrt(a))/(2*sqrt(a))


def zero_limit(weight, matching_cutoff, maximum=24):
    """Finite limit of h_Lambda(y)-E1(y/nu²) as y decreases to zero.

    At lambda=0 the proper compact coordinate gives the exact nonzero
    singular values p*pi/ell_proper. The light singular value is
    epsilon0(y)=y/J²+O(y²), though its canonical Dirac mass stays zero.
    This limit is not a prescription to discard actual spacetime zero modes.
    """
    if not np.isfinite(matching_cutoff) or matching_cutoff <= 0:
        raise ValueError("positive matching cutoff required")
    if not isinstance(maximum, int) or maximum < 1:
        raise ValueError("at least one compact level required")
    ratio = proper_ratio(weight)
    length = weight.interval*ratio
    factor = (pi/(length*weight.cutoff))**2
    massive = 2*float(np.sum(exp1(factor*np.arange(1, maximum+1)**2)))
    value = 2*np.log(weight.cutoff*ratio/matching_cutoff)+massive
    # E1(c*p²)<=exp(-c*p²)/(c*N²), p>N. Integral comparison of the
    # decreasing Gaussian sum gives this bound, evaluated in floating point.
    tail = sqrt(pi)*erfc(sqrt(factor)*maximum)/(factor**1.5*maximum**2)
    return {"proper_ratio": ratio, "proper_length": length,
            "light_singular_slope": 1/ratio**2,
            "matched_weight_at_zero": float(value),
            "massive_zero_argument_weight": massive,
            "zero_limit_tail_bound": float(tail)}


def matched_weight(weight, squared, matching_cutoff):
    if not np.isfinite(squared) or squared < 0:
        raise ValueError("nonnegative squared spectral argument required")
    if squared == 0:
        return zero_limit(weight, matching_cutoff)["matched_weight_at_zero"]
    if not np.isfinite(matching_cutoff) or matching_cutoff <= 0:
        raise ValueError("positive matching cutoff required")
    return weight.response(squared, derivatives=False)["weight"]-float(exp1(squared/matching_cutoff**2))


def spectral_moments(weight, extent=12., tolerance=1e-9):
    """Q1[h]=int h(y)dy and Q2[h]=int y*h(y)dy, with first variations.

    Compact eigenvalues come from the existing owner; no angular or radial
    discretization is introduced. Omitted momentum bounds exclude retained
    Galerkin errors and are evaluated with ordinary floating-point arithmetic.
    """
    if extent < weight.cutoff or tolerance <= 0:
        raise ValueError("resolved positive integration extent/tolerance required")
    names = ("weight", "d_interval", "d_path", "d_cutoff")

    def integrand(k):
        response = weight.response(k*k)
        v = np.array([response[name] for name in names])
        return np.concatenate((2*k*v, 2*k**3*v))

    values, error = quad_vec(integrand, 0, extent, epsabs=tolerance,
                            epsrel=tolerance, norm="max")
    q1, q2 = dict(zip(names, map(float, values[:4]))), dict(zip(names, map(float, values[4:])))
    heat_bound = 1+weight.interval*weight.cutoff/sqrt(pi)
    tail = weight.cutoff**2*heat_bound*np.exp(-(extent/weight.cutoff)**2)
    return {"Q1": q1, "Q2": q2, "quadrature_error_estimate": float(error),
            "Q1_momentum_tail_bound": float(tail),
            "Q2_momentum_tail_bound": float(tail*(extent**2+weight.cutoff**2)),
            "Q1_scaling_residual": weight.interval*q1["d_interval"]-weight.cutoff*q1["d_cutoff"]+2*q1["weight"],
            "Q2_scaling_residual": weight.interval*q2["d_interval"]-weight.cutoff*q2["d_cutoff"]+4*q2["weight"]}


def local_coefficients(weight, moments, matching_cutoff):
    """Known 4D Dirac heat coefficients weighted by the calculated complement.

    Euclidean density: V-A*R_E+C_F*F²+C_W*C_E²+C_Euler*E4+C_box*box_E R_E.
    No R_E² term arises in this free Dirac a4 combination. An independent
    finite R² coefficient in the completed theory is not thereby fixed.
    """
    zero = zero_limit(weight, matching_cutoff)
    q1 = moments["Q1"]["weight"]-matching_cutoff**2
    q2 = moments["Q2"]["weight"]-matching_cutoff**4/2
    q0 = zero["matched_weight_at_zero"]
    norm = (4*pi)**2
    return {**zero, "matching_cutoff": matching_cutoff,
            "Q1_complement": q1, "Q2_complement": q2,
            "V_Dirac": 2*q2/norm, "A_Dirac": q1/(6*norm),
            "C_gauge_Dirac": q0/(3*norm),
            "C_Weyl_Dirac": -q0/(40*norm),
            "C_Euler_Dirac": 11*q0/(720*norm),
            "C_boxR_Dirac": -q0/(60*norm), "C_R2_Dirac": 0.,
            "d_nu_V": -4*matching_cutoff**3/norm,
            "d_nu_A": -matching_cutoff/(3*norm),
            "d_nu_C_gauge": -2/(3*norm*matching_cutoff)}


def flat_moments(weight):
    """Independent analytic momentum integral of the imported unwarped tower."""
    b2 = (np.arange(weight.modes+1)*pi/weight.interval)**2
    copies = np.array([1]+[2]*weight.modes)
    first = integrated_exp1(b2, weight.cutoff)
    second = (weight.cutoff**4*np.exp(-b2/weight.cutoff**2)-b2*first)/2
    return float(copies @ first), float(copies @ second)


def static_basis_coefficients(coefficients, basis_mass=1.):
    """Map to the existing finite_response energy convention, without rerunning it.

    Its basis is (M4,M2 R_L,C2,R_L2,E4,box_L R_L), with R_E=-R_L.
    Static squared invariants and box_E R_E=box_L R_L agree. This reference
    basis mass is neither the matching cutoff nor a determinant normalization.
    """
    if not np.isfinite(basis_mass) or basis_mass <= 0:
        raise ValueError("positive basis reference mass required")
    return [coefficients["V_Dirac"]/basis_mass**4,
            coefficients["A_Dirac"]/basis_mass**2,
            coefficients["C_Weyl_Dirac"], coefficients["C_R2_Dirac"],
            coefficients["C_Euler_Dirac"], coefficients["C_boxR_Dirac"]]
