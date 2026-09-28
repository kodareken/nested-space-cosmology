"""Explicit retained-group order24 enclosure; never an order16 certificate."""
import mpmath as mp

from .nsc_incoming_centered_order24 import CenteredOrder24Geometry, trimmed_defect_jets
from .nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _geometry_jets, _up_float


def retained_order24_bound(config, channel, *, intervals=128, depth=4, precision=40, progress=None):
    """Same collar/operator/state with an explicitly different numerical order.

    This wrapper is reusable for an individually authorized retained group.
    Its record producer calculates only group12; source correction is separate.
    """
    group = channel.get('index')
    if isinstance(group, bool) or group not in range(1, 33) or intervals != 128 or depth != 4:
        raise ValueError('retained non-LLL group, fixed128 cells and depth4 required')
    geometry = CenteredOrder24Geometry(config, precision=precision)
    signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    with _precision(precision):
        _, _, mass, angular, _ = _parameters(channel)
        totals = {s: [mp.iv.mpf(0) for _ in range(25)] for s in signs}
        natural = {s: [mp.iv.mpf(0) for _ in range(25)] for s in signs}
        for index, (rho, center) in enumerate(geometry.cells):
            whole_geometry, point_geometry = _geometry_jets(rho, 28), _geometry_jets(center, 28)
            for sign in signs:
                whole = trimmed_defect_jets(*whole_geometry, mass, sign*angular, order=24)
                point = trimmed_defect_jets(*point_geometry, mass, sign*angular, order=24)
                for j, (w, p) in enumerate(zip(whole, point)):
                    enclosed = centered_coefficient_enclosure(w, p, rho-center, 4)
                    totals[sign][j] += geometry.weight*mp.iv.mpf(_hi(abs(enclosed)))/128
                    natural[sign][j] += geometry.weight*mp.iv.mpf(_hi(abs(w[0])))/128
            if progress is not None and (index+1)%16 == 0: progress(index+1)
        return {'group': group, 'lower': 320. if group in (10, 11, 12, 31, 32) else 160.,
                'intervals': 128, 'physical_Riccati_order': 24, 'centered_remainder_depth': 4,
                'first_inverse_power': 24, 'last_inverse_power': 48, 'directed_precision': precision,
                'per_sign': [{'sign': s, 'coefficient_integral_upper': [_up_float(v) for v in totals[s]],
                             'natural_coefficient_integral_upper': [_up_float(v) for v in natural[s]]}
                            for s in signs],
                'endpoint_weight_upper': _up_float(geometry.weight),
                'horizon_enclosure': [mp.nstr(_lo(geometry.horizon),55), mp.nstr(_hi(geometry.horizon),55)],
                'initial_projector_difference': 0.,
                'initial_scope': 'same declared affine-horizon vacuum; explicit numerical order24 only',
                'finite_offset_initial_projector_term_bound': None,
                'finite_offset_archive_accuracy_certified': False,
                'source_correction_required': True, 'source_correction_applied': False,
                'physical_state_or_action_changed': False,
                'bound_method': 'frozen trimmed24 jets, centered degree4 interval remainder, fixed128 cells'}
