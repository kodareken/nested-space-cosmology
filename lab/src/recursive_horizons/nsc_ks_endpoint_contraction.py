"""Continuous endpoint norms and numerical insertion errors for A+D.

The spatial interpolant is the finite trigonometric polynomial defined by
the captured endpoint. Arb DFT encloses its coefficients; the Fourier
triangle estimate bounds every point of the period, not only target nodes.
The source/preparation and assembly-arithmetic errors remain separate.
"""
from copy import deepcopy
from fractions import Fraction

import numpy as np
from flint import acb, arb, ctx, fmpq

from .nsc_ks_ball_trajectory import complex_ball, exact_upper, restored_upper
from .nsc_ks_finite_matter_error import source_norm_upper
from .nsc_ks_reference_error import homogeneous_reference_norms


def endpoint_difference_norms(difference, weights, energies, period_length, *, bits=120):
    """Uniform weighted Frobenius bounds on delta F and its actual z derivative.

    For each spin row, triangle inequality is applied across Fourier modes,
    while the original source columns retain their weighted l2 norm. The
    carrier derivative is i*(k-E) on each coefficient, including Nyquist's
    owned negative sign. No evolved momentum k=-E substitution is made.
    """
    d = np.asarray(difference, complex)
    weights, energies = np.asarray(weights, float), np.asarray(energies, float)
    if (d.ndim != 3 or d.shape[0] != 2 or d.shape[1] != len(weights)
            or weights.ndim != 1 or energies.shape != weights.shape
            or not len(weights) or d.shape[2] < 8):
        raise ValueError('difference (2, source, grid), weights and energy labels required')
    if (not np.isfinite(d).all() or not np.isfinite(weights).all()
            or not np.isfinite(energies).all() or np.any(weights <= 0)):
        raise ValueError('finite endpoint and original positive weights required')
    n = d.shape[2]
    with ctx.workprec(bits):
        period = Fraction(period_length)
        length = arb(fmpq(period.numerator,period.denominator))
        if not length.is_finite() or not length > 0:
            raise ValueError('positive finite period length required')
        wave = [2*arb.pi()*(k if k < (n+1)//2 else k-n)/length for k in range(n)]
        row, axial = [], []
        for spin in range(2):
            spectra = [[v*arb(float(w))/n for v in acb.dft([complex_ball(x) for x in values])]
                       for w, values in zip(weights, d[spin])]
            norms, derivatives = [], []
            for k in range(n):
                norms.append(sum((c[k].abs_upper()**2 for c in spectra), arb(0)).sqrt())
                derivatives.append(sum(((acb(0,wave[k]-arb(float(e)))*c[k]).abs_upper()**2
                                         for e,c in zip(energies,spectra)),arb(0)).sqrt())
            row.append(sum(norms,arb(0)).upper())
            axial.append(sum(derivatives,arb(0)).upper())
        return {
            'difference_norm': exact_upper(sum((v*v for v in row),arb(0)).sqrt()),
            'difference_axial_norm': exact_upper(sum((v*v for v in axial),arb(0)).sqrt()),
            'whole_spatial_period': True,
            'actual_carrier_derivative_retained': True,
            'nodal_maximum_used_as_bound': False,
        }


def complete_endpoint_inputs(inputs, reference, difference, weights, energies,
                             period_length, propagation, *, bits=120):
    """Fill contraction inputs from the same endpoint and propagated errors.

    Caller authenticates endpoint identity against its stream. This function
    checks source layout and preserves absent covariance/multiplicity slots.
    """
    reference = np.asarray(reference,complex)
    if reference.shape != (2,len(weights)) or not np.isfinite(reference).all():
        raise ValueError('finite reference with original source layout required')
    from .nsc_ks_current_field_cone import MATTER_INPUT_NAMES
    result = deepcopy(inputs)
    norms = endpoint_difference_norms(difference,weights,energies,period_length,bits=bits)
    with ctx.workprec(bits):
        ref, axial = homogeneous_reference_norms(reference,weights,energies,bits=bits)
        values = {'reference_norm':exact_upper(ref),'reference_axial_norm':exact_upper(axial),
                  **{name:norms[name] for name in ('difference_norm','difference_axial_norm')}}
        for name, value in values.items():
            result['matter'][name] = {'value':value,'owner':'authenticated endpoint Fourier norm',
                                     'reason':None}
        for name,key in (
            ('reference_error','weighted_reference_F_error_upper'),
            ('reference_axial_error','weighted_reference_F_z_error_upper'),
            ('difference_error','weighted_difference_F_error_upper'),
            ('difference_axial_error','weighted_difference_F_z_error_upper')):
            result['matter'][name] = {'value':propagation[key],
                'owner':'propagate_difference_error over complete authenticated history',
                'reason':None}
        # A denominator needs a LOWER endpoint. The older generic slot
        # packer makes uppers; reuse the exact rational geometry proof here.
        for name in ('axial_lower','radius_lower'):
            q=Fraction(result['geometry'][name]['exact_rational'])
            lower=arb(fmpq(q.numerator,q.denominator)).lower()
            if not lower > 0:raise ValueError('strict positive geometry lower required')
            result['matter'][name]={'value':exact_upper(lower),
                'owner':'lower dyadic endpoint of owned exact rational geometry bound',
                'reason':None}
    result['endpoint_norm_scope'] = 'whole finite Fourier period at the captured endpoint'
    result['complete_for_difference_matter_error'] = all(
        result['matter'][name]['value'] is not None for name in MATTER_INPUT_NAMES)
    return result


def contract_signed_pair(kwargs, negative_covariance=None, *, bits=120):
    """Contract signed sectors separately; norm isometry does not equate C's."""
    from .nsc_ks_current_field_cone import whole_cone_difference_matter_error
    positive = whole_cone_difference_matter_error(**kwargs,bits=bits)
    result = {'positive':positive, 'negative':None, 'covered_energy_sectors':[1],
              'source_accuracy_included':False, 'physical_local_gate':'OPEN'}
    if negative_covariance is not None:
        negative_kwargs = {**kwargs,'source_norm':source_norm_upper(negative_covariance,bits=bits)}
        result['negative'] = whole_cone_difference_matter_error(**negative_kwargs,bits=bits)
        result['covered_energy_sectors'].append(-1)
    with ctx.workprec(bits):
        for key in ('N','beta'):
            bound = restored_upper(positive[key])
            if result['negative'] is not None:bound += restored_upper(result['negative'][key])
            result[key] = exact_upper(bound)
    result['scope'] = 'numerical insertion for selected source rows; not full-family coverage'
    return result
