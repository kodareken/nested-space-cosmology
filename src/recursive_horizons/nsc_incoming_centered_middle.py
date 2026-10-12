"""Centered order16 enclosure for the retained finite middle-band source."""
import mpmath as mp

from .nsc_incoming_centered_order24 import trimmed_defect_jets
from .nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _geometry_jets, _up_float
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_projector_energy_bound import RefinedDefectGeometry


def centered_middle_bound(config, channel, *, intervals=128, precision=40, progress=None):
    """Retain the original16-term projector; refine only its error enclosure.

    The group14 order24 representation has its own separate certificate.
    This owner exposes the other retained groups without selecting geometry,
    redefining source data, or rerunning any mode equation.
    """
    group = channel.get('index')
    if group not in range(1, 33) or group == 14 or intervals != 128:
        raise ValueError('retained non-LLL group other than14 and fixed128-cell enclosure required')
    geometry = RefinedDefectGeometry(config, intervals=intervals, precision=precision)
    signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    with _precision(precision):
        _, _, mass, angular, _ = _parameters(channel)
        totals = {s: [mp.iv.mpf(0) for _ in range(17)] for s in signs}
        natural = {s: [mp.iv.mpf(0) for _ in range(17)] for s in signs}
        for index, (rho, _, _) in enumerate(geometry.cells):
            center = mp.iv.mpf((_lo(rho)+_hi(rho))/2)
            if not _lo(rho) <= _lo(center) <= _hi(center) <= _hi(rho):
                raise ArithmeticError('Taylor center outside its interval')
            whole_geometry = _geometry_jets(rho, 20)
            point_geometry = _geometry_jets(center, 20)
            for sign in signs:
                whole = trimmed_defect_jets(*whole_geometry, mass, sign*angular, order=16)
                point = trimmed_defect_jets(*point_geometry, mass, sign*angular, order=16)
                for j, (w, p) in enumerate(zip(whole, point)):
                    enclosed = centered_coefficient_enclosure(w, p, rho-center, 4)
                    totals[sign][j] += geometry.weight*mp.iv.mpf(_hi(abs(enclosed)))/intervals
                    natural[sign][j] += geometry.weight*mp.iv.mpf(_hi(abs(w[0])))/intervals
            if progress is not None and (index+1)%16 == 0: progress(index+1)
        return {'group': group, 'lower': 320. if group in (10, 11, 12, 31, 32) else 160.,
                'intervals': intervals, 'physical_Riccati_order': 16, 'centered_remainder_depth': 4,
                'directed_precision': precision, 'first_inverse_power': 16, 'last_inverse_power': 32,
                'per_sign': [{'sign': s, 'I16_to_I32_upper': [_up_float(v) for v in totals[s]],
                             'natural_I16_to_I32_upper': [_up_float(v) for v in natural[s]]} for s in signs],
                'endpoint_weight_upper': _up_float(geometry.weight),
                'horizon_enclosure': [mp.nstr(_lo(geometry.horizon),55), mp.nstr(_hi(geometry.horizon),55)],
                'initial_projector_difference': 0.,
                'initial_scope': 'declared affine-horizon limit; finite-delta input artifacts unchanged',
                'finite_offset_initial_projector_term_bound': None, 'finite_offset_archive_accuracy_certified': False,
                'physical_state_or_approximation_changed': False,
                'bound_method': 'trimmed order16 jets and centered degree4 interval remainder, fixed128 cells'}
