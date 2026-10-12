"""Conditional continuous error bounds for the owned reference/difference split.

No sampled residual, convergence estimate, source error or tangent error is
promoted to a certificate here. Every input must be supplied by its proof
owner on the same backward cone and with the same source weights.
"""
from flint import arb, ctx

from .nsc_ks_ball_trajectory import exact_upper
from .nsc_ks_characteristic_error import _upper


def propagate_difference_error(reference_initial, reference_residual,
                               reference_offdiagonal_integral,
                               difference_initial, difference_residual_integrals,
                               offdiagonal_integral, Bz_integral,
                               M_integral, Mz_integral,
                               maximum_absolute_energy, *, bits=160):
    """Majorize e_A and e_D without treating A+D as two independent fields.

    A is spatially homogeneous in envelope coordinates. Thus e_A,z = 0.
    The error equations are e_A'=G_ref e_A+R_A and
    e_D'=L_g e_D+M e_A+R_D. After z differentiation the additional
    terms are B_z e_D+M_z e_A. Integrating factors give

      a = exp(K_ref) (a0 + RA),
      q = d0 + RD + Mint*a,
      d <= exp(K) q,
      d_z <= exp(K) (dz0 + RDz + K1*q + Mzint*a).

    All integrals bound absolute residuals or induced row norms continuously;
    using their whole-interval values also bounds every preceding subinterval.
    The diagonal imaginary source phases give no row-norm growth. Row norms
    are l2 in weighted source columns, L-infinity in z, for each spin row.
    """
    if (not isinstance(difference_initial, (tuple, list))
            or not isinstance(difference_residual_integrals, (tuple, list))
            or len(difference_initial) != 2 or len(difference_residual_integrals) != 2):
        raise ValueError("difference field and axial residual pairs required")
    with ctx.workprec(bits):
        a0, ra, kr, k, kz, mint, mzint, energy = (
            _upper(value, name) for value, name in (
                (reference_initial, "reference initial error"),
                (reference_residual, "reference residual integral"),
                (reference_offdiagonal_integral, "reference offdiagonal integral"),
                (offdiagonal_integral, "offdiagonal integral"),
                (Bz_integral, "B_z integral"),
                (M_integral, "M integral"), (Mz_integral, "M_z integral"),
                (maximum_absolute_energy, "source energy")))
        d0, dz0 = (_upper(v, "difference initial error") for v in difference_initial)
        rd, rdz = (_upper(v, "difference residual integral")
                   for v in difference_residual_integrals)
        a = kr.exp() * (a0 + ra)
        q = d0 + rd + mint*a
        d = k.exp()*q
        dz = k.exp()*(dz0 + rdz + kz*q + mzint*a)
        factor = arb(2).sqrt()
        return {
            "reference_row_error_upper": exact_upper(a),
            "difference_row_error_upper": exact_upper(d),
            "difference_axial_envelope_error_upper": exact_upper(dz),
            "weighted_reference_F_error_upper": exact_upper(factor*a),
            "weighted_reference_F_z_error_upper": exact_upper(factor*energy*a),
            "weighted_difference_F_error_upper": exact_upper(factor*d),
            "weighted_difference_F_z_error_upper": exact_upper(factor*(dz+energy*d)),
            "reference_envelope_spatial_derivative": "identically zero",
            "requires_continuous_cone_bounds": True,
            "source_accuracy_included": False,
            "tangent_accuracy_included": False,
            "physical_local_gate": "OPEN: input residual enclosures required",
        }


def difference_matter_error(reference_norm, difference_norm,
                            reference_axial_norm, difference_axial_norm,
                            reference_error, difference_error,
                            reference_axial_error, difference_axial_error,
                            source_norm, *, mass, absolute_angular,
                            axial_lower, radius_lower, multiplicity, bits=160):
    """Bound the cross/quadratic contraction at the same incoming geometry.

    Norms refer to computed weighted F_ref, delta F, F_ref,z and delta F_z.
    The covariance, geometry and weights are fixed. The density difference
    is D C A* + A C D* + D C D*. The unsymmetrized momentum difference
    is Dz C A* + Az C D* + Dz C D*. Expanding perturbations in these
    expressions gives the nonnegative polynomials below. Hermitian and
    anti-Hermitian symmetrization do not enlarge the corresponding bound.
    The same vertex operator norms as finite_matter_error are used.

    Errors in the covariance, weights, vertices and contraction arithmetic
    are separate obligations. The reference error used in the contraction
    may differ from the reference error that forced the nodal D solve;
    the latter belongs in propagate_difference_error's inputs.
    """
    with ctx.workprec(bits):
        values = [reference_norm, difference_norm, reference_axial_norm,
                  difference_axial_norm, reference_error, difference_error,
                  reference_axial_error, difference_axial_error, source_norm,
                  mass, absolute_angular, axial_lower, radius_lower, multiplicity]
        a, d, az, dz, ea, ed, eaz, edz, cov, m, ell, scale, radius, mu = (
            _upper(v, "difference contraction bound") for v in values)
        if not scale > 0 or not radius > 0 or not mu > 0:
            raise ValueError("positive geometry and multiplicity required")
        density = cov*(2*(a*ed + d*ea + ea*ed + d*ed) + ed**2)
        current = cov*(
            dz*ea + a*edz + edz*ea
            + az*ed + d*eaz + eaz*ed
            + dz*ed + d*edz + edz*ed)
        potential = (m*m + (ell/radius)**2).sqrt()
        return {
            "N": exact_upper(mu*(potential*density + current/scale)),
            "beta": exact_upper(mu*current),
            "density_error_upper": exact_upper(density),
            "momentum_error_upper": exact_upper(current),
            "correlated_reference_and_difference": True,
            "source_accuracy_included": False,
            "arithmetic_accuracy_included": False,
            "physical_local_gate": "OPEN: continuous field error inputs required",
        }
