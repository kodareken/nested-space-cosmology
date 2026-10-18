"""One-loop massless Dirac trace on the actual quiet one-sided backgrounds.

Free four-component Dirac quantization and the Laplace-type heat coefficient
are imported. No singular wall curvature is squared. A trace projection is
not the full stress, a state, gravitational slip or self-consistent Lambda.
"""
from __future__ import annotations

import math

from .nsc_atlas_dirac_probe import quiet_reference_curvatures


def dirac_bulk_trace(*, R, Ricci_squared, Riemann_squared):
    """Return <T^mu_mu>/(hbar*c) for a smooth constant-R region, in length^-4.

    For L=D^2=-(nabla^2+E), E=-R/4 and tr I=4,
    tr Omega_mn Omega^mn=-Riemann^2/2. The traced heat numerator is
    5R^2-8Ricci^2-7Riemann^2. The fermionic determinant reverses its sign.
    The total-derivative term is omitted only because Box R=0 here.
    """
    R, Ric, Riem = (float(v) for v in (R, Ricci_squared, Riemann_squared))
    if not all(math.isfinite(v) for v in (R, Ric, Riem)):
        raise ValueError('smooth bulk curvature invariants must be finite')
    return (7*Riem+8*Ric-5*R*R)/(360*(4*math.pi)**2)


def reference_bulk_trace_report(radius=.45, *, ell_m=10000.):
    r, ell = float(radius), float(ell_m)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError('fiducial physical length must be positive and finite')
    geometry = quiet_reference_curvatures(r)
    coefficients = {}
    for name, row in geometry.items():
        shape = dirac_bulk_trace(R=row['R'], Ricci_squared=row['Ricci_squared'],
                                Riemann_squared=row['Riemann_squared'])
        coefficients[name] = dict(trace_times_ell4_over_hbar_c=shape,
            delta_Lambda_trace_shape=-2*math.pi*shape)
    # Same supplied physical constants and 10 km fiducial used by the scalar
    # trace/P001 scale calculations. h and c exact; G is supplied measured input.
    h_SI, c_SI, G_SI = 6.62607015e-34, 299792458., 6.67430e-11
    lP2 = (h_SI/(2*math.pi))*G_SI/c_SI**3
    zeta = lP2/ell**2
    for row in coefficients.values():
        row['delta_Lambda_trace_times_ell2'] = row['delta_Lambda_trace_shape']*zeta
        row['delta_Lambda_trace_per_m2'] = row['delta_Lambda_trace_times_ell2']/ell**2
    return dict(inputs=dict(Rwall_over_ell=r, H_ell=1., M_over_ell=.4,
            k_wall_ell=2., Dirac_mass=0., four_component_Dirac_fields=1,
            ell_m=ell, ell_is_predicted=False, G_SI=G_SI, h_SI=h_SI, c_SI=c_SI,
            Planck_length_squared_over_ell_squared=zeta,
            independent_new_quantum_state_is_selected=False),
        one_sided=coefficients,
        checks=dict(cap_shape_error=abs(coefficients['cap']['delta_Lambda_trace_shape']-11/(120*math.pi)),
            exterior_shape_error=abs(coefficients['exterior']['delta_Lambda_trace_shape']
                                    +7*.4**2/(60*math.pi*r**6))),
        counts=dict(algebraic_wall_points=1, local_one_loop_bulk_coefficients=2,
            nonlocal_loop_polarizations=0, P001_time_steps=0, Dirac_time_steps=0,
            new_evolution_grids=0),
        nested_Psi_over_Phi=None,
        scope='One smooth-bulk trace term of the declared Dirac probe on the original quiet metric. No interface distribution products, full vacuum stress, state, nonlocal polarization, backreaction or derived cosmological Lambda.',
        uncertainty_scope='Algebraic curvature coefficient conditional on imported free Dirac theory; supplied G and fiducial length, not a total physical uncertainty or a converged full quantum model.')
