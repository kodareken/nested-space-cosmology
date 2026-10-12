"""Explicit order24 projector remainder; original order16 records are intact."""
import mpmath as mp

from .nsc_incoming_vacuum_tail_bound import (
    _precision, _lo, _hi, _range, _up_float, _geometry_jets, _defect_coefficients,
)
from .nsc_incoming_projector_energy_bound import RefinedDefectGeometry
from .nsc_incoming_middle_bound import _parameters


class HigherOrderDefectGeometry(RefinedDefectGeometry):
    """A new numerical vacuum representation on the unchanged geometry."""
    def __init__(self, config, *, intervals=128, order=24, precision=40):
        if order != 24:
            raise ValueError('this version owns numerical order24 only')
        super().__init__(config, intervals=intervals, precision=precision)
        self.order = order
        with _precision(precision):
            self.cells = [(rho, *_geometry_jets(rho, order)) for rho, _, _ in self.cells]

    def bound_channel(self, channel, lower):
        group = channel.get('index')
        if group != 14 or lower != 160.:
            raise ValueError('this new order24 pilot owns group14 and its existing endpoint160 only')
        with _precision(self.precision):
            _, _, mass, angular, _ = _parameters(channel)
            rows = []
            for sign in (1, -1):
                integrals = [mp.iv.mpf(0) for _ in range(self.order+1)]
                for _, A, sphere in self.cells:
                    values = _defect_coefficients(A, sphere, mass, sign*angular, order=self.order)
                    for i, value in enumerate(values):
                        integrals[i] += self.weight*mp.iv.mpf(_hi(abs(value)))/self.intervals
                rows.append({'sign': sign, 'coefficient_integral_upper': [_up_float(v) for v in integrals]})
            return {'group': group, 'lower': float(lower), 'intervals': self.intervals,
                    'physical_Riccati_order': self.order, 'first_inverse_power': self.order,
                    'last_inverse_power': 2*self.order, 'directed_precision': self.precision,
                    'per_sign': rows, 'initial_projector_difference': 0.,
                    'initial_scope': 'same declared affine-horizon vacuum; new numerical order24 only',
                    'finite_offset_archive_accuracy_certified': False,
                    'horizon_enclosure': [mp.nstr(_lo(self.horizon), 55), mp.nstr(_hi(self.horizon), 55)],
                    'old_order16_source_replaced': False,
                    'bound_method': 'same recurrence with24 numerical terms and positive directed radial sums'}
