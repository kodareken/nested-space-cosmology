"""Node-wise local-action owner for a supplied homogeneous KS history.

This module applies the already locked four-dimensional NSC local ledger.  It
does not select a history, add an interface action, or turn endpoint forces
into constants.  Einstein includes its existing GHY completion, Maxwell and
Weyl are bulk terms, and Euler/box-R are retained as endpoint functionals.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi

import numpy as np


CHANNELS = (
    "einstein_ghy",
    "maxwell_bulk",
    "weyl_bulk",
    "euler_endpoint",
    "boxR_endpoint",
)
FIELDS = ("lapse", "beta", "q_ADM", "r")


@dataclass(frozen=True)
class GeneralKSLocalInducedHistory:
    """Locked local action evaluated on ``N,beta,q_ADM,r`` history nodes."""

    einstein_coefficient: float
    gauge_coefficient: float
    weyl_coefficient: float
    euler_coefficient: float
    boxR_coefficient: float
    magnetic_flux: int

    def __post_init__(self):
        values = (
            self.einstein_coefficient,
            self.gauge_coefficient,
            self.weyl_coefficient,
            self.euler_coefficient,
            self.boxR_coefficient,
        )
        if not all(isfinite(value) for value in values):
            raise ValueError("finite locked local coefficients required")
        if self.einstein_coefficient <= 0 or self.gauge_coefficient <= 0:
            raise ValueError("positive Einstein and gauge coefficients required")
        if isinstance(self.magnetic_flux, bool) or not isinstance(self.magnetic_flux, int):
            raise ValueError("integer magnetic flux required")

    @staticmethod
    def _arrays(history, *, allow_complex=False):
        names = ("time", "lapse", "shift", "a_parallel", "radius")
        dtype = complex if allow_complex else float
        time, lapse, beta, radial, radius = (
            np.asarray(getattr(history, name), dtype=dtype) for name in names
        )
        count = len(time)
        if count < 7 or any(value.shape != (count,) for value in (lapse, beta, radial, radius)):
            raise ValueError("at least seven equal-length KS history nodes required")
        if not np.isfinite(np.stack((time, lapse, beta, radial, radius))).all():
            raise ValueError("finite KS history required")
        if np.any(np.diff(time.real) <= 0):
            raise ValueError("strictly increasing history time required")
        if min(lapse.real) <= 0 or min(radial.real) <= 0 or min(radius.real) <= 0:
            raise ValueError("positive lapse and scale factors required")
        return time.real, lapse, beta, radial, radius

    @staticmethod
    def _derivative(values, time):
        return np.gradient(values, time, edge_order=2)

    def _channel_actions(self, time, lapse, beta, radial, radius):
        # Homogeneous scalar invariants do not depend on the shift.  It remains
        # an explicit fourth field so its zero variation is returned node-wise.
        del beta
        derivative = lambda value: self._derivative(value, time)
        lapse_dot = derivative(lapse)
        radial_dot = derivative(radial)
        radius_dot = derivative(radius)
        radial_ddot = derivative(radial_dot)
        radius_ddot = derivative(radius_dot)

        section_a = (
            radial_ddot-radial_dot*lapse_dot/lapse
        )/(lapse*lapse*radial)
        section_b = radial_dot*radius_dot/(lapse*lapse*radial*radius)
        section_c = (
            radius_ddot-radius_dot*lapse_dot/lapse
        )/(lapse*lapse*radius)
        section_d = (1+(radius_dot/lapse)**2)/(radius*radius)
        scalar = -2*(section_a+2*section_b+2*section_c+section_d)
        weyl_squared = 4*(section_a-section_b-section_c+section_d)**2/3

        # Existing spherical reduction after the existing four-dimensional GHY
        # completion.  Per unit homogeneous axial coordinate.
        einstein_density = 8*pi*self.einstein_coefficient*(
            lapse*radial
            -2*radius*radius_dot*radial_dot/lapse
            -radial*radius_dot*radius_dot/lapse
        )
        maxwell_density = (
            -2*pi*self.gauge_coefficient*self.magnetic_flux**2
            * lapse*radial/(radius*radius)
        )
        weyl_density = (
            -4*pi*self.weyl_coefficient*lapse*radial*radius*radius*weyl_squared
        )

        integrate = lambda value: np.trapezoid(value, time)
        euler_primitive = 32*pi*(radial_dot/lapse)*(1+(radius_dot/lapse)**2)
        scalar_dot = derivative(scalar)
        box_primitive = 4*pi*radial*radius*radius*scalar_dot/lapse
        return np.asarray((
            integrate(einstein_density),
            integrate(maxwell_density),
            integrate(weyl_density),
            -self.euler_coefficient*(euler_primitive[-1]-euler_primitive[0]),
            -self.boxR_coefficient*(box_primitive[-1]-box_primitive[0]),
        ))

    def channel_actions(self, history):
        arrays = self._arrays(history)
        values = self._channel_actions(*arrays)
        return {name: float(value.real) for name, value in zip(CHANNELS, values)}

    def action_gradients(self, history, *, complex_step=1e-28):
        """Return ``dS_channel/d(field_node)`` for all four ADM histories."""
        if not isfinite(complex_step) or complex_step <= 0:
            raise ValueError("positive finite complex step required")
        time, lapse, beta, radial, radius = self._arrays(history)
        fields = [lapse.astype(complex), beta.astype(complex),
                  radial.astype(complex), radius.astype(complex)]
        output = {channel: {field: np.zeros(len(time)) for field in FIELDS}
                  for channel in CHANNELS}
        for field_index, field_name in enumerate(FIELDS):
            if field_name == "beta":
                continue
            for node in range(len(time)):
                varied = [value.copy() for value in fields]
                varied[field_index][node] += 1j*complex_step
                values = self._channel_actions(time, *varied)
                for channel, value in zip(CHANNELS, values):
                    output[channel][field_name][node] = value.imag/complex_step
        return output

    def evaluate(self, history, *, stationarity_tolerance=3e-11):
        """Return composable actions, node forces and the open endpoint gate."""
        if not isfinite(stationarity_tolerance) or stationarity_tolerance <= 0:
            raise ValueError("positive finite stationarity tolerance required")
        actions = self.channel_actions(history)
        gradients = self.action_gradients(history)
        forces = {
            channel: {field: (-values).tolist() for field, values in rows.items()}
            for channel, rows in gradients.items()
        }
        weyl_endpoint = {
            field: [float(gradients["weyl_bulk"][field][0]),
                    float(gradients["weyl_bulk"][field][-1])]
            for field in FIELDS
        }
        endpoint_norm = max(abs(value) for pair in weyl_endpoint.values() for value in pair)
        return {
            "actions": actions,
            "total_action": float(sum(actions.values())),
            "action_gradients": {
                channel: {field: values.tolist() for field, values in rows.items()}
                for channel, rows in gradients.items()
            },
            "action_forces": forces,
            "weyl_endpoint_completion_residual": {
                "components": weyl_endpoint,
                "maximum_absolute": float(endpoint_norm),
                "tolerance": stationarity_tolerance,
                "below_tolerance": bool(endpoint_norm <= stationarity_tolerance),
                "completion_present": False,
            },
            "beta_variation_maximum_absolute": float(max(
                abs(value) for channel in CHANNELS
                for value in gradients[channel]["beta"]
            )),
        }
