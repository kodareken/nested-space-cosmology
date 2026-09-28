"""Allocate the stored compact spectral source to its canonical origins.

No compact eigenproblem or source integral runs here. The free KK mass map
and flat/warped moments are retained inputs. The conversion is Euclidean;
it does not invent a real-time state for a spectral subtraction.
"""
from math import log, pi
from types import SimpleNamespace

from .nsc_compact_matching import zero_limit


def _coefficients(q):
    """Apply the existing one-Dirac heat coefficient map linearly."""
    q0, q1, q2 = (q[k] for k in ('Q0', 'Q1', 'Q2'))
    norm = (4*pi)**2
    return {'V_Dirac': 2*q2/norm, 'A_Dirac': q1/(6*norm),
            'C_gauge_Dirac': q0/(3*norm), 'C_Weyl_Dirac': -q0/(40*norm),
            'C_Euler_Dirac': 11*q0/(720*norm),
            'C_boxR_Dirac': -q0/(60*norm), 'C_R2_Dirac': 0.}


def allocate(record):
    """Three origins of H_nu = h_warped - E1(y/nu²), from a pinned record.

    The unwarped weight is one light Dirac plus two Dirac fields at each
    n>=1, m_n=n*pi/ell. 'canonical' describes these fields; their raw
    finite-proper-time modulus is still distinct from a canonical state IF.
    """
    parameters = record['conventions']['parameters']
    cutoff, interval = parameters['Lambda'], parameters['ell']
    flat = record['flat_moment_check']
    warped = record['refined_moments']
    # The existing analytic zero-limit formula needs only these fields;
    # constructing a Galerkin compact operator is neither needed nor done.
    unwarped = SimpleNamespace(cutoff=cutoff, interval=interval,
                              path=0., amplitude=parameters['amplitude'])
    zero = zero_limit(unwarped, cutoff)
    tower = {'Q0': zero['matched_weight_at_zero'],
             'Q1': flat['analytic_Q1'] - cutoff**2,
             'Q2': flat['analytic_Q2'] - cutoff**4/2}
    anchor = record['matched_coefficients'][0]
    measure = {
        'Q0': anchor['matched_weight_at_zero']
              - 2*log(cutoff/anchor['matching_cutoff']) - tower['Q0'],
        'Q1': warped['Q1']['weight'] - flat['analytic_Q1'],
        'Q2': warped['Q2']['weight'] - flat['analytic_Q2'],
    }
    rows = []
    for old in record['matched_coefficients']:
        nu = old['matching_cutoff']
        conversion = {'Q0': 2*log(cutoff/nu),
                      'Q1': cutoff**2-nu**2,
                      'Q2': (cutoff**4-nu**4)/2}
        columns = {
            name: {'moments': dict(q), 'coefficients': _coefficients(q)}
            for name, q in [('massive_canonical_tower', tower),
                            ('light_endpoint_conversion', conversion),
                            ('warp_covariant_conversion', measure)]
        }
        summed = {key: sum(part['moments'][key] for part in columns.values())
                  for key in tower}
        summed_coefficients = {
            key: sum(part['coefficients'][key] for part in columns.values())
            for key in _coefficients(tower)}
        old_moments = {'Q0': old['matched_weight_at_zero'],
                       'Q1': old['Q1_complement'], 'Q2': old['Q2_complement']}
        rows.append({
            'matching_cutoff': nu, 'origins': columns,
            'sum_moments': summed, 'sum_coefficients': summed_coefficients,
            'stored_moments': old_moments,
            'stored_coefficients': {key: old[key] for key in summed_coefficients},
            'moment_reconstruction_residuals':
                {key: summed[key]-old_moments[key] for key in summed},
            'coefficient_reconstruction_residuals':
                {key: summed_coefficients[key]-old[key] for key in summed_coefficients},
        })
    return {
        'parameters': parameters,
        'canonical_fields': {'mass_formula': 'm_n=n*pi/ell',
                             'massless_Dirac_copies': 1, 'positive_level_Dirac_copies': 2,
                             'interval': 'ell=2 L_star, fixed chiral endpoint domain',
                             'constant_size_free_light_heavy_link': '0'},
        'canonical_Q0_tail_bound': zero['zero_limit_tail_bound'],
        'canonical_Q1_Q2_origin': 'stored analytic flat-tower moments; no reintegration',
        'rows': rows,
        'warp_path_moment_derivatives': {
            'Q1': warped['Q1']['d_path'], 'Q2': warped['Q2']['d_path'],
            'origin': 'stored warped derivatives; free canonical masses and fixed canonical state have no warp-path dependence'},
    }
