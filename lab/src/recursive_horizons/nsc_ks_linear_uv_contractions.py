"""Constraint contractions of the owned linear fixed-transfer projector pair.

This is an endpoint coefficient identity, not a finite-history tail bound.
The E-linear vertex requires the next diagonal projector coefficient.
"""
import sympy as sp


def contraction_identity():
    a, r = sp.symbols('a r', positive=True)
    m, ell = sp.symbols('m ell', real=True)
    what, x, y = sp.symbols('what q3_01 q3_10')  # complex Fourier amplitudes
    S1 = sp.Matrix([[0, 1], [1, 0]])
    S2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    S3 = sp.diag(1, -1)
    P0 = sp.diag(1, 0)
    # Reused from endpoint_pair_identity in nsc_incoming_fixed_transfer.
    P1 = a * (m * S1 - ell / r * S2) / 2
    Q2 = -ell * a**2 * what / (4 * r**2) * S1
    diagonal = m * ell * a**3 * what / (4 * r**2) * S3
    Q3 = diagonal + sp.Matrix([[0, x], [y, 0]])
    pair = P0 * Q3 + Q3 * P0 - Q3 + P1 * Q2 + Q2 * P1
    H0, H1 = -m * S1 + ell / r * S2, -S3 / a
    expressions = {
        'N_E_minus1': sp.trace(H1 * Q2),
        'N_E_minus2': sp.trace(H0 * Q2 + H1 * Q3),
        'beta_E_minus1': sp.trace(Q2),
        'beta_E_minus2': sp.trace(Q3),
    }
    residuals = {k: str(sp.simplify(v)) for k, v in expressions.items()}
    residuals['projector_pair_E_minus3'] = [str(sp.simplify(v)) for v in pair]
    if any(v != '0' for k, v in residuals.items() if k != 'projector_pair_E_minus3') or set(residuals['projector_pair_E_minus3']) != {'0'}:
        raise ArithmeticError('linear endpoint contraction identity failed')
    missing_diagonal = sp.simplify(sp.trace(H0 * Q2))
    if missing_diagonal == 0:
        raise ArithmeticError('control must detect omission of the forced diagonal')
    return {
        'fixed_transfer': 'Eo=E+omega/2; Ei=E-omega/2; omega fixed as E tends to positive infinity',
        'Q2': '-ell*a^2*what*sigma1/(4*r^2)',
        'forced_diagonal_Q3': 'm*ell*a^3*what*sigma3/(4*r^2)',
        'unfixed_Q3_entries': ['q3_01', 'q3_10'],
        'residuals': residuals,
        'omitted_diagonal_false_N_coefficient': str(missing_diagonal),
        'angular_pairing_needed_for_these_zeros': False,
        'complex_Fourier_amplitudes_preserved': True,
        'midpoint_offset_order': 'omega shifts P1/E first at E^-2, hence enters P_o Q+Q P_i first at E^-4',
        'order_conclusion': 'linear vacuum N,beta insertions are O(E^-3) at fixed omega after the E^-2 coefficient cancels',
        'remainder_input': 'owned Q expansion through E^-3 with O(E^-4) pair remainder, multiplied by an at-most-E vertex',
        'full_source_addition': 'the unchanged coherent thermal source remainder must also be bounded; it is not set to zero',
        'scope': {'linear_response_at_reference': True, 'finite_history_bound': False,
                  'uniform_in_axial_transfer': False, 'numeric_tail_constant': None,
                  'pointwise_on_I_tail_bound': None, 'physical_local_gate': 'OPEN'},
    }
