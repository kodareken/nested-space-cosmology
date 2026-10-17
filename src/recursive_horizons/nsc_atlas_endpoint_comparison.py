"""Three-grid diagnostics and an explicitly scaled Jordan wall observable.

This module performs no evolution or fit. Resolution spreads describe the
supplied numerical results, not confidence intervals or continuum bounds.
"""
import math


FIDUCIAL_LENGTH_METRES = 10000.0
ENDPOINT_REFERENCE_RADIUS = 0.55
ENDPOINT_SOURCE_X = 0.4815042082576364


def ssprk3_diagnostic_increment(width, stage_values):
    """Integrate a diagnostic with the actual RK weights.

    Stage order is left, right predictor, midpoint, accepted endpoint.
    The accepted endpoint has zero weight in this RK update.
    """
    if len(stage_values) != 4 or not math.isfinite(width) or width <= 0:
        raise ValueError('positive step and four finite stage values required')
    if not all(math.isfinite(v) for v in stage_values):
        raise ValueError('positive step and four finite stage values required')
    return width*(stage_values[0]+stage_values[1]+4*stage_values[2])/6


def jordan_radius(radius_e, scalar_c, alpha):
    """Invariant conversion using c, not a canonical-field approximation."""
    values = (radius_e, scalar_c, alpha)
    if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
        raise ValueError('finite real radius, invariant scalar and alpha required')
    if radius_e <= 0 or alpha < 0:
        raise ValueError('positive Einstein radius and nonnegative alpha required')
    # hypot avoids overflow in the sum of squares for representable outputs.
    result = math.hypot(radius_e, math.sqrt(alpha/6)*scalar_c)
    if not math.isfinite(result):
        raise ValueError('Jordan radius is not representable')
    return result


def refinement(gaps):
    """Report measured contraction without assigning an asymptotic order."""
    if len(gaps) != 3 or any(type(v) not in (float, int) or not math.isfinite(v)
                            or v < 0 for v in gaps):
        raise ValueError('three finite nonnegative gaps required')
    return dict(values=list(gaps), strictly_decreasing=gaps[0] > gaps[1] > gaps[2],
                fine_over_coarse=[b/a if a else None for a,b in zip(gaps, gaps[1:])],
                continuum_error_bound=False)


def scaled_radius(radii, *, length_metres=FIDUCIAL_LENGTH_METRES):
    """Freeze a fiducial scale and report the observed three-grid envelope.

    The fine-grid value is the numerical prediction. The symmetric range
    covers all three supplied values. It is deliberately named a resolution
    envelope: unresolved common input, quadrature and model errors are not
    numerically estimated by varying the grid.
    """
    if (len(radii) != 3 or any(type(r) not in (float,int) or not math.isfinite(r)
                             or r <= 0 for r in radii)):
        raise ValueError('three positive finite Jordan radii required')
    if (type(length_metres) not in (float,int) or not math.isfinite(length_metres)
            or length_metres <= 0):
        raise ValueError('positive finite physical length required')
    values = [length_metres*r for r in radii]
    fine = values[-1]
    half_width = max(abs(v-fine) for v in values)
    reference = length_metres*ENDPOINT_REFERENCE_RADIUS
    if not all(math.isfinite(v) for v in (*values, half_width, reference)):
        raise ValueError('physical scale conversion is not representable')
    return dict(
        observable='Jordan material-wall areal radius at the fixed incoming-label endpoint',
        definition='R_J = sqrt(Area_J/(4*pi)) = circumference_J/(2*pi)',
        event=dict(affine_source_x=ENDPOINT_SOURCE_X, normalized_V=1.0,
                   reference_radius_over_ell=ENDPOINT_REFERENCE_RADIUS),
        scale=dict(ell_metres=length_metres, role='declared_fiducial_input',
                   fitted_to_output=False, measured_astrophysical_scale=False),
        value_metres=fine, grid_values_metres=values,
        quiet_reference_metres=reference,
        difference_from_quiet_reference_metres=fine-reference,
        resolution_half_width_metres=half_width,
        resolution_range_metres=[fine-half_width, fine+half_width],
        circumference_metres=2*math.pi*fine,
        uncertainty=dict(kind='observed_three_grid_resolution_envelope',
            coverage='all three numerical endpoint values with the same input data',
            confidence_level=None, rigorous_error_enclosure=False,
            unresolved=['common source interpolation error',
                        'common radial quadrature and floating-point error',
                        'error beyond the tested grid family',
                        'physical preparation and model identification'],
            physical_scale_rule='an independently measured relative ell error rescales R_J by the same relative amount'),
        measurement_protocol='A spherical array comoving with the material wall measures its intrinsic area or equatorial proper circumference on the cross-section selected by the normalized incoming null label. This is a local wall observable; the calculation does not supply its reception by a distant or central detector.',
        observational_validation=False)
