"""Local characteristic error majorant for the canonical KS envelope.

The inputs must bound residual integrals on the backward cone of I. This
consumer does not certify sampled residual values or source preparation.
"""
from fractions import Fraction

from flint import arb, fmpq, ctx

from .nsc_ks_ball_trajectory import exact_upper, restored_upper


def _upper(value, name):
    if value is None or isinstance(value, bool):
        raise ValueError('explicit nonnegative bound required: ' + name)
    if isinstance(value, arb):
        result = value
    else:
        try:
            q = Fraction(value)
            result = arb(fmpq(q.numerator, q.denominator))
        except (ValueError, TypeError, OverflowError) as error:
            raise ValueError('finite bound required: ' + name) from error
    if not result.is_finite() or not result >= 0:
        raise ValueError('finite nonnegative bound required: ' + name)
    return result.upper()


def propagate_local_characteristics(initial_rows, residual_integrals,
                                     offdiagonal_integral, Bz_integral,
                                     maximum_absolute_energy, *, bits=120):
    """Return directed local weighted-field errors, conditional on input bounds.

    Row norms are l2 in the original weighted source columns and L-infinity
    in z on the shrinking backward cone. Each input pair lists orders 0,1.
    """
    if initial_rows is None or residual_integrals is None or len(initial_rows) != 2 or len(residual_integrals) != 2:
        raise ValueError('field and first-derivative bounds required')
    with ctx.workprec(bits):
        e0, e1 = (_upper(v, 'initial row error') for v in initial_rows)
        R0, R1 = (_upper(v, 'continuous cone residual integral') for v in residual_integrals)
        K0, K1, energy = (_upper(v, name) for v, name in (
            (offdiagonal_integral, 'offdiagonal integral'), (Bz_integral, 'B_z integral'),
            (maximum_absolute_energy, 'maximum source energy')))
        growth = K0.exp()
        field = growth*(e0+R0)
        derivative = growth*(e1+R1+K1*(e0+R0))
        factor = arb(2).sqrt()
        return {'row_field_error_upper': exact_upper(field),
                'row_axial_envelope_error_upper': exact_upper(derivative),
                'weighted_F_error_upper': exact_upper(factor*field),
                'weighted_F_z_error_upper': exact_upper(factor*(derivative+energy*field)),
                'growth_upper': exact_upper(growth),
                'requires_continuous_cone_bounds': True,
                'physical_local_gate': 'OPEN: residual and source provenance required'}
