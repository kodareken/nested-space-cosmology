"""Rectangular physical-metric trust region for the corrected retarded gate."""
from dataclasses import dataclass

import numpy as np

from .nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily


def flatten_constraints(gradient, jacobian):
    gradient = np.asarray(gradient, float)
    jacobian = np.asarray(jacobian, float)
    if gradient.ndim != 2 or gradient.shape[1] != 2:
        raise ValueError("gradient must have shape (nodes,2)")
    if jacobian.ndim != 3 or jacobian.shape[1:] != gradient.shape:
        raise ValueError("jacobian must have shape (directions,nodes,2)")
    residual = np.concatenate((gradient[:, 0], gradient[:, 1]))
    matrix = np.vstack((jacobian[:, :, 0].T, jacobian[:, :, 1].T))
    return residual, matrix


def joint_merit(gradient):
    gradient = np.asarray(gradient, float)
    if gradient.ndim != 2 or gradient.shape[1] != 2 or not np.isfinite(gradient).all():
        raise ValueError("finite (nodes,2) gradient required")
    return float(np.max(np.abs(gradient)))


def physical_metric(family, principal_matrix, *, samples=2001):
    """H3(w) plus H1(U), scaled by the certified geometric principal columns."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    principal = np.asarray(principal_matrix, float)
    if principal.shape != (2, 2) or not np.isfinite(principal).all():
        raise ValueError("finite 2x2 principal matrix required")
    n = len(family.coefficients[0])
    if isinstance(samples, bool) or int(samples) != samples or samples < 2:
        raise ValueError("at least two metric samples required")
    z = np.linspace(family.center - family.axial_outer,
                    family.center + family.axial_outer, int(samples))
    dz = float(z[1] - z[0])
    eye = np.eye(n)
    values = np.empty((n, 4, len(z)))
    for index in range(n):
        function = LocalAxialFunction(tuple(eye[index]), family.center,
                                      family.axial_inner, family.axial_outer)
        for order in range(4):
            values[index, order] = [function(point, order) for point in z]
    gram = lambda orders: sum((values[:, order] @ values[:, order].T
                               for order in orders), np.zeros((n, n))) * dz
    metric = np.zeros((2 * n, 2 * n))
    metric[:n, :n] = float(principal[:, 1] @ principal[:, 1]) * gram(range(4))
    metric[n:, n:] = float(principal[:, 0] @ principal[:, 0]) * gram(range(2))
    eigenvalues = np.linalg.eigvalsh(metric)
    if eigenvalues[0] <= 0 or not np.isfinite(eigenvalues).all():
        raise ArithmeticError("physical history metric must be positive definite")
    return metric


def metric_norm(vector, metric):
    vector, metric = np.asarray(vector, float), np.asarray(metric, float)
    if metric.shape != (len(vector), len(vector)):
        raise ValueError("metric shape must match vector")
    value = float(vector @ metric @ vector)
    if value < -1e-14 or not np.isfinite(value):
        raise ArithmeticError("finite nonnegative metric norm required")
    return float(np.sqrt(max(0., value)))


def levenberg_step(residual, jacobian, metric, radius):
    """Minimize the linear least-squares model inside a physical metric ball."""
    residual = np.asarray(residual, float).ravel()
    jacobian = np.asarray(jacobian, float)
    metric = np.asarray(metric, float)
    radius = float(radius)
    if (jacobian.ndim != 2 or jacobian.shape[0] != len(residual)
            or metric.shape != (jacobian.shape[1], jacobian.shape[1])):
        raise ValueError("residual, rectangular Jacobian and metric shapes differ")
    if radius <= 0 or not np.isfinite(radius):
        raise ValueError("positive finite trust radius required")
    normal = jacobian.T @ jacobian
    slope = jacobian.T @ residual

    def solve(mu):
        step = np.linalg.solve(normal + mu * metric, -slope)
        return step, metric_norm(step, metric)

    try:
        step, size = solve(0.)
    except np.linalg.LinAlgError:
        step, size = solve(np.finfo(float).eps)
    multiplier = 0.
    if size > radius:
        lower, upper = 0., 1.
        while solve(upper)[1] > radius:
            upper *= 4.
            if not np.isfinite(upper):
                raise ArithmeticError("failed to bracket trust multiplier")
        for _ in range(64):
            multiplier = (lower + upper) / 2
            step, size = solve(multiplier)
            if size > radius:
                lower = multiplier
            else:
                upper = multiplier
        multiplier = upper
        step, size = solve(multiplier)
    predicted = residual + jacobian @ step
    return {
        "step": step,
        "predicted_residual": predicted,
        "metric_norm": size,
        "trust_radius": radius,
        "levenberg_multiplier": multiplier,
        "predicted_merit": float(np.max(np.abs(predicted))),
        "parent_merit": float(np.max(np.abs(residual))),
    }


def broyden_update(jacobian, step, measured_change):
    jacobian = np.asarray(jacobian, float)
    step = np.asarray(step, float).ravel()
    measured_change = np.asarray(measured_change, float).ravel()
    if jacobian.shape != (len(measured_change), len(step)):
        raise ValueError("Broyden shapes differ")
    denominator = float(step @ step)
    if denominator <= np.finfo(float).tiny:
        raise ValueError("nonzero accepted step required for Broyden update")
    return jacobian + np.outer(measured_change - jacobian @ step, step) / denominator


def acceptance(parent_gradient, candidate_gradient, predicted_residual, *,
               radius_positive, preparation_supported):
    parent = np.asarray(parent_gradient, float)
    candidate = np.asarray(candidate_gradient, float)
    predicted = np.asarray(predicted_residual, float).ravel()
    parent_flat = np.concatenate((parent[:, 0], parent[:, 1]))
    candidate_flat = np.concatenate((candidate[:, 0], candidate[:, 1]))
    if predicted.shape != parent_flat.shape or candidate.shape != parent.shape:
        raise ValueError("acceptance residual shapes differ")
    parent_merit = float(np.max(np.abs(parent_flat)))
    candidate_merit = float(np.max(np.abs(candidate_flat)))
    predicted_merit = float(np.max(np.abs(predicted)))
    predicted_reduction = parent_merit - predicted_merit
    measured_reduction = parent_merit - candidate_merit
    ratio = (measured_reduction / predicted_reduction
             if predicted_reduction > 0 else float("-inf"))
    accepted = bool(radius_positive and preparation_supported
                    and measured_reduction > 0 and predicted_reduction > 0)
    return {
        "accepted": accepted,
        "parent_merit": parent_merit,
        "predicted_merit": predicted_merit,
        "candidate_merit": candidate_merit,
        "predicted_reduction": predicted_reduction,
        "measured_reduction": measured_reduction,
        "prediction_ratio": ratio,
        "radius_positive": bool(radius_positive),
        "preparation_supported": bool(preparation_supported),
    }


def update_radius(radius, ratio, accepted, *, minimum=1e-12, maximum=None):
    radius = float(radius)
    maximum = radius * 8 if maximum is None else float(maximum)
    if not accepted or ratio < .25:
        return max(minimum, radius / 4)
    if ratio > .75:
        return min(maximum, radius * 2)
    return radius


@dataclass(frozen=True)
class RefreshPolicy:
    accepted_since_refresh: int = 0
    poor_ratios_since_refresh: int = 0

    def update(self, *, accepted, prediction_ratio):
        accepted_count = self.accepted_since_refresh + int(bool(accepted))
        poor = self.poor_ratios_since_refresh + int(
            (not np.isfinite(prediction_ratio)) or prediction_ratio < .25)
        refresh = accepted_count >= 5 or poor >= 2
        return RefreshPolicy(0, 0) if refresh else RefreshPolicy(accepted_count, poor), refresh
