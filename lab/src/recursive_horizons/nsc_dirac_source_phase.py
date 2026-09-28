"""Bounded numerical owner of the proven formal coefficient f_s.

This evaluates the boxed geometry integral from
nsc-dirac-source-phase-transport.md and its full history tangent. It does
not evolve a source, insert a stress or contact term into N or beta, or
prove the uniform remainder / ordered coincidence limit.
"""
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from types import MappingProxyType
import json

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq

from .nsc_compatible_history_geometry import (
    CompatibleIncomingMetric,
    _plateau,
    _profile as _direction_profile,
)
from .nsc_ks_profile_identity import profile_description, profile_identity
from .nsc_ks_source_envelope import RHO_UP_MIN, _as_metric, _radius_lower_bound, continuum_speed_distance
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_pg_ks_metric_pullback import reference_chart


SOURCE_SIGNS = (1, -1)
RHO_SIGMA = 1.0
DEFAULT_GAUSS_NODES = 16
PHYSICAL_LOCAL_GATE = (
    'OPEN: formal geometric f_s only; uniform remainder and ordered '
    'source/coincidence limit are parent-owned; no assembly correction')


def _real_vector(value, name):
    raw = np.asarray(value)
    if np.iscomplexobj(raw):
        raise ValueError(name + ' must be real')
    try:
        result = np.array(raw, dtype=float, copy=True)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError('finite real ' + name + ' required') from error
    if result.ndim == 0:
        result = result.reshape(1)
    if result.ndim != 1 or result.size == 0 or not np.isfinite(result).all():
        raise ValueError('nonempty finite one-dimensional ' + name + ' required')
    return result


def _real_scalar(value, name):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError('explicit finite real ' + name + ' required')
    result = np.asarray(value)
    if result.shape != () or np.iscomplexobj(result):
        raise ValueError('scalar real ' + name + ' required')
    number = float(result)
    if not np.isfinite(number):
        raise ValueError('finite ' + name + ' required')
    return number


def _gauss_count(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError('positive integer Gauss-Legendre count required')
    count = int(value)
    if count < 1:
        raise ValueError('positive integer Gauss-Legendre count required')
    return count


def _normal_coordinate(rho):
    """s=T(rho)-T(1) on the original chart; independent of the radius family."""
    return chart_coordinates(float(rho))[0] - chart_coordinates(RHO_SIGMA)[0]


@lru_cache(maxsize=4096)
def _characteristic_distance(rho):
    """D(rho)=int_1^rho a^{-2} d rho of the original chart; no amplitude tangent."""
    rho = float(rho)
    reference_chart(rho)
    if rho == RHO_SIGMA:
        return 0.0
    if rho < RHO_SIGMA:
        raise ValueError('characteristic distance starts at the incoming rho=1')
    return float(continuum_speed_distance(rho, RHO_SIGMA))


def _hex(value):
    return float(value).hex()


def coefficient_binding(family, z, *, angular, rho_up, gauss_nodes=DEFAULT_GAUSS_NODES):
    """Identity of the analytic profiles, windows, chart data and quadrature options."""
    z = _real_vector(z, 'target z')
    angular = _real_scalar(angular, 'angular ell')
    rho_up = _real_scalar(rho_up, 'rho_up')
    gauss_nodes = _gauss_count(gauss_nodes)
    description = {
        'schema': 'NSC-FORMAL-SOURCE-PHASE-COEFFICIENT-v1',
        'profile': profile_description(family, include_normal_window=True),
        'rho_sigma': _hex(RHO_SIGMA),
        'rho_up': _hex(rho_up),
        'angular': _hex(angular),
        'z': [_hex(v) for v in z],
        'source_signs': [int(s) for s in SOURCE_SIGNS],
        'source_sign_axis': 0,
        'amplitude_axis': 1,
        'z_axis': -1,
        'amplitude_order': 'CompatibleIncomingMetric amplitudes; LocalIncomingFamily is w then U',
        'gauss_nodes': gauss_nodes,
        'integral': 'ell^2/2 int_1^rho_up (r_g^{-2}-r_ref^{-2})(rho, z-s D(rho)) d rho',
        'chart': 'original T(rho), a(rho); D=int a^{-2} d rho; r_ref=sqrt(1+rho^2)',
        'cache': False,
    }
    digest = sha256(json.dumps(description, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    return digest, description


def _transition_radii(metric, rho_up):
    s_up = _normal_coordinate(rho_up)
    widths = {float(d.inner_radius) for d in metric.directions}
    widths.update(float(d.outer_radius) for d in metric.directions)
    roots = []
    for width in sorted(widths):
        if not s_up < -width < 0:
            continue
        root = brentq(lambda rho, w=width: _normal_coordinate(rho) + w,
                      RHO_SIGMA, rho_up, xtol=2e-14, rtol=1e-12)
        roots.append(float(root))
    return tuple(roots)


def _panels(metric, rho_up):
    points = [RHO_SIGMA]
    for rho in _transition_radii(metric, rho_up):
        if rho > points[-1] + 1e-12 and rho < rho_up - 1e-12:
            points.append(rho)
    points.append(float(rho_up))
    return tuple((points[i], points[i + 1]) for i in range(len(points) - 1))


def _gauss_nodes(panels, count):
    nodes, weights = [], []
    x, w = leggauss(count)
    for left, right in panels:
        half = 0.5 * (right - left)
        if half <= 0:
            continue
        nodes.append(half * x + 0.5 * (left + right))
        weights.append(half * w)
    if not nodes:
        raise ValueError('quadrature panels must cover (1, rho_up]')
    return np.concatenate(nodes), np.concatenate(weights)


def _basis_jets(direction, normal, space):
    """Analytic (delta r, delta r_z) of one compatible radius direction."""
    space = np.asarray(space, dtype=float)
    window = _plateau(np.full(space.shape, float(normal)),
                      direction.inner_radius, direction.outer_radius)
    value = np.zeros(space.shape, dtype=float)
    deriv = np.zeros(space.shape, dtype=float)
    active = window != 0
    if not np.any(active):
        return value, deriv
    w0 = _direction_profile(direction.w, space[active], 0, 'w')
    w1 = _direction_profile(direction.w, space[active], 1, 'w')
    u0 = _direction_profile(direction.U, space[active], 0, 'U')
    u1 = _direction_profile(direction.U, space[active], 1, 'U')
    s = float(normal)
    cubic = s * s * s / 6.0
    scale = window[active]
    value[active] = scale * (s * w0 + cubic * u0)
    deriv[active] = scale * (s * w1 + cubic * u1)
    if not np.isfinite(value).all() or not np.isfinite(deriv).all():
        raise ValueError('radius direction is not finite')
    return value, deriv


def _geometry(metric, rho, z_char):
    """Actual r, r_z and every amplitude basis at the original-chart point."""
    r_ref = float(np.sqrt(1.0 + rho * rho))
    normal = _normal_coordinate(rho)
    n_amp = len(metric.amplitudes)
    basis = np.empty((n_amp,) + z_char.shape, dtype=float)
    basis_z = np.empty_like(basis)
    for index, direction in enumerate(metric.directions):
        basis[index], basis_z[index] = _basis_jets(direction, normal, z_char)
    amplitudes = np.asarray(metric.amplitudes, dtype=float)
    delta_r = np.tensordot(amplitudes, basis, axes=(0, 0))
    radius = r_ref + delta_r
    radius_z = np.tensordot(amplitudes, basis_z, axes=(0, 0))
    if not np.isfinite(radius).all() or np.any(radius <= 0):
        raise ValueError('actual radius must stay finite and positive')
    if not np.isfinite(radius_z).all():
        raise ArithmeticError('nonfinite analytic radius z-derivative')
    return radius, radius_z, r_ref, basis, basis_z, delta_r


def _assert_endpoint_radius(metric, rho, z, distance):
    for sign in SOURCE_SIGNS:
        z_char = z - sign * distance
        _geometry(metric, rho, z_char)


@dataclass(frozen=True)
class FormalSourcePhaseCoefficient:
    """f_s, f_s,z and amplitude tangents on both original source characteristics.

    Axis 0 is SOURCE_SIGNS=(+1,-1). Amplitude axis follows the bound metric,
    including zero-amplitude directions. Error bounds are missing on purpose.
    """
    signs: tuple
    z: object
    f: object
    f_z: object
    delta_f: object
    delta_f_z: object
    delta_n2: object
    angular: float
    rho_up: float
    characteristic_distance: float
    amplitudes: tuple
    profile_identity: str
    binding: str
    quadrature: object
    coefficient_error_bound: object
    remainder_error_bound: object
    physical_error_bound: object
    physical_local_gate: str
    assembly_correction: bool


def _freeze(array):
    result = np.array(array, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def formal_source_phase_coefficient(
        family, z, *, angular, rho_up, gauss_nodes=DEFAULT_GAUSS_NODES):
    """Integrate f_s and every retained history tangent in one rho quadrature.

    ell=0 returns exact zeros without quadrature. Zero amplitudes still keep
    their direction tangents. Characteristic D and a come from the original
    chart only; they are not varied with the radius history.
    """
    z = _real_vector(z, 'target z')
    angular = _real_scalar(angular, 'angular ell')
    rho_up = _real_scalar(rho_up, 'rho_up')
    gauss_nodes = _gauss_count(gauss_nodes)
    if rho_up < RHO_UP_MIN:
        raise ValueError('caller rho_up>=1.03 required')
    reference_chart(rho_up)
    local, metric = _as_metric(family)
    if not isinstance(metric, CompatibleIncomingMetric) or not metric.directions:
        raise TypeError('LocalIncomingFamily or CompatibleIncomingMetric pure-radius geometry required')
    binding, description = coefficient_binding(
        family, z, angular=angular, rho_up=rho_up, gauss_nodes=gauss_nodes)
    identity = profile_identity(family, include_normal_window=True)
    # profile_identity requires the owned analytic profiles, so the existing
    # coefficient bound applies globally, including between quadrature nodes.
    radius_lower, global_bound = _radius_lower_bound(local, metric, None, None)
    if not global_bound or radius_lower <= 0:
        raise ValueError('positive global radius bound required')
    if _normal_coordinate(rho_up) > -max(d.outer_radius for d in metric.directions):
        raise ValueError('history and all tangent directions must be flat at the fixed upstream slice')
    n_sign, n_amp, n_z = len(SOURCE_SIGNS), len(metric.amplitudes), z.size
    zeros = np.zeros((n_sign, n_z), dtype=float)
    zero_tangents = np.zeros((n_sign, n_amp, n_z), dtype=float)
    distance_up = _characteristic_distance(rho_up)
    quadrature = {
        'rule': 'composite Gauss-Legendre on normal-window split panels',
        'gauss_nodes': gauss_nodes,
        'panels': [],
        'split_rho': [],
        'n_nodes': 0,
        'skipped': False,
        'convergence_indicator': None,
        'certified_quadrature_error_bound': None,
        'radius_lower_bound': radius_lower,
    }

    def result(f, f_z, delta_f, delta_f_z, skipped):
        quadrature['skipped'] = bool(skipped)
        f, f_z = _freeze(f), _freeze(f_z)
        delta_f, delta_f_z = _freeze(delta_f), _freeze(delta_f_z)
        signs = np.asarray(SOURCE_SIGNS, dtype=float)
        return FormalSourcePhaseCoefficient(
            signs=SOURCE_SIGNS,
            z=_freeze(z),
            f=f,
            f_z=f_z,
            delta_f=delta_f,
            delta_f_z=delta_f_z,
            delta_n2=_freeze(signs[:, None] * f_z),
            angular=angular,
            rho_up=rho_up,
            characteristic_distance=distance_up,
            amplitudes=tuple(float(v) for v in metric.amplitudes),
            profile_identity=identity,
            binding=binding,
            quadrature=MappingProxyType({
                **quadrature,
                'binding_schema': description['schema'],
                'source_signs': list(SOURCE_SIGNS),
                'shapes': {
                    'f': list(f.shape),
                    'f_z': list(f_z.shape),
                    'delta_f': list(delta_f.shape),
                    'delta_f_z': list(delta_f_z.shape),
                    'delta_n2': [n_sign, n_z],
                },
            }),
            coefficient_error_bound=None,
            remainder_error_bound=None,
            physical_error_bound=None,
            physical_local_gate=PHYSICAL_LOCAL_GATE,
            assembly_correction=False,
        )

    if angular == 0.0:
        return result(zeros, zeros, zero_tangents, zero_tangents, True)

    panels = _panels(metric, rho_up)
    nodes, weights = _gauss_nodes(panels, gauss_nodes)
    quadrature['panels'] = [list(map(float, panel)) for panel in panels]
    quadrature['split_rho'] = [float(left) for left, _ in panels[1:]]
    quadrature['n_nodes'] = int(nodes.size)
    ell2 = angular * angular
    f, f_z = np.zeros_like(zeros), np.zeros_like(zeros)
    delta_f, delta_f_z = np.zeros_like(zero_tangents), np.zeros_like(zero_tangents)
    panel_mass = np.zeros(len(panels), dtype=float)
    panel_rights = np.fromiter((right for _, right in panels), dtype=float)
    endpoints = {RHO_SIGMA, float(rho_up)}
    endpoints.update(quadrature['split_rho'])
    for rho in sorted(endpoints):
        _assert_endpoint_radius(metric, float(rho), z, _characteristic_distance(float(rho)))

    for rho, weight in zip(nodes, weights):
        rho = float(rho)
        distance = _characteristic_distance(rho)
        panel_index = min(int(np.searchsorted(panel_rights, rho, side='right')), len(panels) - 1)
        for sign_index, sign in enumerate(SOURCE_SIGNS):
            z_char = z - sign * distance
            radius, radius_z, r_ref, basis, basis_z, delta_r = _geometry(metric, rho, z_char)
            # Stable difference of reciprocals: 1/r^2-1/r_ref^2 = -dr(2 r_ref+dr)/(r^2 r_ref^2).
            inv_diff = -delta_r * (2.0 * r_ref + delta_r) / (radius * radius * r_ref * r_ref)
            inv_r3 = 1.0 / (radius * radius * radius)
            inv_r4 = inv_r3 / radius
            f[sign_index] += weight * (ell2 / 2.0) * inv_diff
            f_z[sign_index] += weight * (-ell2) * radius_z * inv_r3
            delta_f[sign_index] += weight * (-ell2) * basis * inv_r3
            delta_f_z[sign_index] += weight * (-ell2) * (
                basis_z * inv_r3 - 3.0 * radius_z * basis * inv_r4)
            panel_mass[panel_index] += float(np.max(np.abs(weight * (ell2 / 2.0) * inv_diff)))

    if not all(np.isfinite(v).all() for v in (f, f_z, delta_f, delta_f_z)):
        raise ArithmeticError('nonfinite formal source-phase coefficient')
    if not np.any(metric.amplitudes):
        f[:] = 0.0
        f_z[:] = 0.0
    quadrature['panel_absolute_contribution_sums'] = [float(v) for v in panel_mass]
    return result(f, f_z, delta_f, delta_f_z, False)
