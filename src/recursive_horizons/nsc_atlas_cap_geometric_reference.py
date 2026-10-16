"""One regular-centre, fixed-radius mass-gradient reference at alpha=0.

Continuous linear chi and u in z=atanh(r) give exact element integrals.
The incoming wall condition is a weak boundary flux; its discrepancy from
the last element's derivative is reported. This is not a time integrator,
a finite-backreaction cap, or a replacement for the background cap producer.
"""
import numpy as np
from scipy.linalg import solve_banded

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_polar_cap import _scalar, _vector


def _positive_gap(z):
    """1/z - 2/sinh(2z), without small-z subtraction.

    For 2z<=1 the positive sinh remainder is summed through x**21/21!.
    Its first omitted term, relative to the leading x**3/6, is at most
    6/23! divided by 1-1/(24*25), below 2.4e-22. This bounds series
    truncation in real arithmetic, not floating-point function evaluation.
    """
    x = 2*z
    if x <= 1:
        term, tail = 1., 1.
        for k in range(1, 10):
            term *= x*x/((2*k+2)*(2*k+3))
            tail += term
        return (2*z/3)*tail/(np.sinh(x)/x)
    return 1/z-2/np.sinh(x)


def geometric_cap_reference(z_nodes, chi, velocity, *, T_x, incoming_w_plus, budget):
    """Evaluate a single weak RHS, mass gradient and physical flux.

    ``chi`` and ``velocity=chi_T`` are continuous linear nodal fields with
    zero centre values. Radius is fixed; alpha and the centre lapse are zero.
    The caller supplies one positive clock rate and w_plus=Pi+Phi. Coefficient
    work is charged conservatively before the tridiagonal solve. Returned
    arrays are owned by this result, never aliases of the caller's inputs.
    """
    budget = _budget(budget)
    z = _vector(z_nodes, "z_nodes").copy()
    chi = _vector(chi, "chi").copy()
    velocity = _vector(velocity, "velocity").copy()
    tx = _scalar(T_x, "T_x", positive=True)
    incoming = _scalar(incoming_w_plus, "incoming_w_plus")
    if (not 2 <= len(z) <= 129 or z[0] != 0 or np.any(np.diff(z) <= 0)
            or chi.shape != z.shape or velocity.shape != z.shape):
        raise ValueError("require 2..129 ordered nodes from zero and matching field arrays")
    if chi[0] != 0 or velocity[0] != 0:
        raise ValueError("regular centre requires chi(0)=velocity(0)=0")
    before = budget.evaluations
    budget.check(40*len(z)+32)
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            R = float(np.tanh(z[-1]))
            if not 0 < R < 1:
                raise ValueError("wall radius must be representable strictly between zero and one")
            A = float(1/np.cosh(z[-1])**2)
            gap = float(_positive_gap(z[-1]))
            beta = gap*z[-1]
            if not np.isfinite(gap) or gap <= 0 or beta <= 0:
                raise ValueError("positive geometric trace gap is not resolved")
            h = np.diff(z)
            slope = np.diff(chi)/h
            mean = chi[-1]/z[-1]
            shifted = slope-mean
            # Free-node stiffness action for K-(A/R)e_W e_W^T. This
            # form retains the small positive gap for a linear chi.
            stiffness = np.r_[shifted[:-1]-shifted[1:], shifted[-1]+gap*chi[-1]]
            diagonal = np.r_[(h[:-1]+h[1:])/3, h[-1]/3]
            off = h[1:]/6
            band = np.zeros((3, len(z)-1))
            band[1] = diagonal
            band[0, 1:] = off
            band[2, :-1] = off
            mass_velocity = diagonal*velocity[1:]
            mass_velocity[:-1] += off*velocity[2:]
            mass_velocity[1:] += off*velocity[1:-1]
            drive = A*R*incoming
            rhs = -stiffness.copy()
            rhs[-1] += drive-velocity[-1]
            acceleration = solve_banded((1, 1), band, rhs, check_finite=True)
            applied = diagonal*acceleration
            applied[:-1] += off*acceleration[1:]
            applied[1:] += off*acceleration[:-1]
            chi_rate = tx*velocity
            velocity_rate = np.r_[0., tx*acceleration]
            gradient_chi, gradient_velocity = stiffness/2, mass_velocity/2
            kinetic = np.sum(h*((velocity[:-1]+velocity[1:])**2
                                +velocity[:-1]**2+velocity[1:]**2)/6)
            shape = np.sum(h*shifted*shifted)
            potential = shape+gap*chi[-1]**2
            mass = (kinetic+potential)/4
            if mass == 0 and (np.any(chi != 0) or np.any(velocity != 0)):
                raise ValueError("positive mass underflows for the supplied nonzero fields")
            mass_rate = float(gradient_chi@chi_rate[1:]
                              +gradient_velocity@velocity_rate[1:])
            flux = velocity[-1]*(drive-velocity[-1])/2
            weak_v = A*(R*incoming+chi[-1]/R)-velocity[-1]
            trace_incoming = (velocity[-1]+slope[-1])/(A*R)-chi[-1]/R**2
            result = dict(
                alpha=0., radius=R, radius_rate=0., centre_lapse=0., T_x=tx,
                z_nodes=z, chi=chi, velocity=velocity,
                rates=dict(chi=chi_rate, velocity=velocity_rate),
                gradient=dict(chi=gradient_chi, velocity=gradient_velocity),
                mass=float(mass), kinetic_mass=float(kinetic/4),
                potential_mass=float(potential/4), positive_trace_gap=gap,
                coercivity_factor=float(beta),
                mass_rate=mass_rate, physical_mass_flux=float(flux),
                balance_defect=float(mass_rate-tx*flux),
                incoming_flux_bound=float(tx*drive*drive/8),
                linear_solve_defect=float(np.max(abs(applied-rhs))),
                centre=dict(f=float(chi[1]/z[1]), Pi=float(velocity[1]/z[1]), Phi=0.),
                wall=dict(A=A, f=float(chi[-1]/R), u=float(velocity[-1]),
                          weak_v=float(weak_v), reconstructed_v=float(slope[-1]),
                          derivative_discrepancy=float(slope[-1]-weak_v),
                          incoming_w_plus=incoming,
                          reconstructed_incoming=float(trace_incoming),
                          incoming_discrepancy=float(trace_incoming-incoming)),
                coefficient_evaluations=budget.evaluations-before,
                domain="zero backreaction, fixed radius, one instantaneous weak RHS",
            )
    except (FloatingPointError, np.linalg.LinAlgError) as error:
        raise ValueError("unresolved geometric reference arithmetic") from error
    def finite(value):
        if isinstance(value, dict):
            return all(finite(v) for v in value.values())
        return isinstance(value, str) or bool(np.all(np.isfinite(value)))
    if not finite(result):
        raise ValueError("nonfinite geometric reference result")
    budget.check()
    return result
