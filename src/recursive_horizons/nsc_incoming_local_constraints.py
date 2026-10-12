"""Local lapse/shift Euler forces of the existing NSC local actions.

The order-four Taylor data are a differentiation germ, not an evolved or
selected spacetime history. Complex Cauchy contours differentiate the owned
densities; they are numerical radii, not physical times or regulators.
Neither the physical covariance nor the reference band is evaluated here.
"""
from math import factorial

import numpy as np

from .nsc_incoming_cauchy_jets import IncomingCauchyJets, FIELDS
from .nsc_light_restoration_action import LightRestorationAction
from .nsc_spherical_local_history import LockedSphericalLocalAction
from .nsc_spatial_reference_symbol import SymbolJet, INDEX, INDICES


TWO_JETS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))


class _DensityBatch:
    """Identity integral: use existing action APIs to obtain their densities."""
    def __init__(self, count):
        self.time = np.empty(count)
        self.radius_coordinate = np.empty(1)

    @staticmethod
    def integral(value):
        return value


def _product(left, right):
    """Product of scalar value/first/second derivative arrays."""
    x, dx, ddx = left
    y, dy, ddy = right
    return (x*y, dx*y[..., None]+x[..., None]*dy,
            ddx*y[..., None, None]+x[..., None, None]*ddy
            +dx[..., :, None]*dy[..., None, :]
            +dy[..., :, None]*dx[..., None, :])


def raw_metric_two_jets(raw):
    """Map independent (N,beta,a,r) two-jets to the SAME canonical KS metric.

    Input shape (...,4,6); entries are actual derivatives, with Tz stored
    once. Varying that entry changes both symmetric Hessian entries once.
    """
    raw = np.asarray(raw)
    if raw.shape[-2:] != (4, 6) or not np.isfinite(raw).all():
        raise ValueError("four raw KS fields with six finite two-jets required")
    scalars = []
    for i in range(4):
        row = raw[..., i, :]
        dd = np.empty((*row.shape[:-1], 2, 2), dtype=raw.dtype)
        dd[..., 0, 0], dd[..., 0, 1] = row[..., 3], row[..., 4]
        dd[..., 1, 0], dd[..., 1, 1] = row[..., 4], row[..., 5]
        scalars.append((row[..., 0], row[..., 1:3], dd))
    N, beta, axial, radius = scalars
    nn, aa, bb = _product(N, N), _product(axial, axial), _product(beta, beta)
    aabb, aab = _product(aa, bb), _product(aa, beta)
    components = ((tuple(x-y for x, y in zip(nn, aabb)), tuple(-v for v in aab)),
                  (tuple(-v for v in aab), tuple(-v for v in aa)))
    shape = raw.shape[:-2]
    g = np.empty((*shape, 2, 2), dtype=raw.dtype)
    dg = np.empty((*shape, 2, 2, 2), dtype=raw.dtype)
    ddg = np.empty((*shape, 2, 2, 2, 2), dtype=raw.dtype)
    for i in range(2):
        for j in range(2):
            g[..., i, j], dg[..., :, i, j], ddg[..., :, :, i, j] = components[i][j]
    return g, radius[0], (dg, ddg, radius[1], radius[2])


def local_action_densities(raw, light, compact):
    """Unmodified light and locked-local actions, before spacetime integration.

    The result includes the compact Euler density as a cancellation control.
    The constant BoxR coefficient is a boundary functional, not a bulk term.
    """
    if not isinstance(light, LightRestorationAction) or not isinstance(compact, LockedSphericalLocalAction):
        raise ValueError("the existing light and locked compact action owners are required")
    if light.spectrum.charge != abs(compact.ledger['magnetic_flux']):
        raise ValueError("light and compact owners must use the same magnetic sector")
    raw = np.asarray(raw)
    shape = raw.shape[:-2]
    batch = raw.reshape(-1, 1, 4, 6)
    grid = _DensityBatch(len(batch))
    data = raw_metric_two_jets(batch)
    la, _ = light.actions(grid, *data, canonical_coordinates="KS")
    ca, _ = compact.actions(grid, *data)
    return {**{"light/"+k: v.reshape(shape) for k, v in la.items()},
            **{"compact/"+k: v.reshape(shape) for k, v in ca.items()}}


def scalar_taylor_coefficients(fields):
    """Read real scalar order-four SymbolJets; preserve every supplied jet."""
    if len(fields) != 4 or any(not isinstance(f, SymbolJet) for f in fields):
        raise ValueError("four raw KS SymbolJets required")
    count = fields[0].data.shape[1]
    coefficients = np.zeros((count, 4, 5, 5), dtype=np.longdouble)
    for f, field in enumerate(fields):
        if (not isinstance(field, SymbolJet) or field.physical != 4
                or field.data.shape[1:] != (count, 2, 2)
                or not np.isfinite(field.data).all()
                or np.max(abs(field.data.imag)) > 3e-12
                or np.max(abs(field.data-field.data[..., :1, :1]*np.eye(2))) > 3e-12):
            raise ValueError("real scalar metric jets through total coordinate order four required")
        if any(np.max(abs(field.data[i])) > 3e-12 for i, (_, _, k) in enumerate(INDICES) if k):
            raise ValueError("metric jets cannot depend on reference momentum")
        for t in range(5):
            for z in range(5-t):
                coefficients[:, f, t, z] = field.data[INDEX[t, z, 0], :, 0, 0].real
    if np.any(coefficients[:, (0, 2, 3), 0, 0] <= 0):
        raise ValueError("positive intrinsic lapse, axial scale and radius required")
    return coefficients


def taylor_two_jets(coefficients, time, space):
    """Evaluate only the two-jets of the supplied finite Taylor germ."""
    c = np.asarray(coefficients)
    if c.ndim != 4 or c.shape[1:] != (4, 5, 5):
        raise ValueError("batched raw Taylor coefficients of shape (n,4,5,5) required")
    T, Z = np.broadcast_arrays(np.asarray(time, np.clongdouble), np.asarray(space, np.clongdouble))
    output = np.zeros((len(c), *T.shape, 4, 6), dtype=np.clongdouble)
    for j, (dt, dz) in enumerate(TWO_JETS):
        for t in range(dt, 5):
            for z in range(dz, 5-t):
                factor = factorial(t)*factorial(z)/(factorial(t-dt)*factorial(z-dz))
                value = factor*T**(t-dt)*Z**(z-dz)
                output[..., j] += c[:, :, t, z].reshape(len(c), *([1]*T.ndim), 4)*value[..., None]
    return output


def _cauchy_euler(coefficients, light, compact, fields, samples, coordinate_radius, jet_radius):
    pi = np.arccos(np.longdouble(-1))
    angles = 2*pi*np.arange(samples, dtype=np.longdouble)/samples
    roots = np.exp(np.clongdouble(1j)*angles)
    gradients = None
    for B, name in enumerate(fields):
        for j, (dt, dz) in enumerate(TWO_JETS):
            rt = roots if dt else np.zeros(1, dtype=np.clongdouble)
            rz = roots if dz else np.zeros(1, dtype=np.clongdouble)
            T, Z, J = np.meshgrid(coordinate_radius*rt, coordinate_radius*rz,
                                  jet_radius*roots, indexing="ij")
            raw = taylor_two_jets(coefficients, T, Z)
            raw[..., FIELDS.index(name), j] += J
            density = local_action_densities(raw, light, compact)
            if gradients is None:
                gradients = {k: np.zeros((len(coefficients), len(fields)), np.clongdouble) for k in density}
            # Cauchy coefficient of jet_probe^1 T^dt z^dz; converting the
            # coordinate coefficient to a derivative supplies both factorials.
            wt = rt**(-dt) if dt else np.ones(1)
            wz = rz**(-dz) if dz else np.ones(1)
            weight = wt[:, None, None]*wz[None, :, None]/roots[None, None, :]
            factor = ((-1)**(dt+dz)*factorial(dt)*factorial(dz)
                      /(jet_radius*coordinate_radius**(dt+dz)))
            for k, value in density.items():
                gradients[k][:, B] += factor*np.mean(value*weight, axis=(1, 2, 3))
    return gradients


def local_euler_gradients(fields, light, compact, *, varied=("N", "beta"), samples=12,
                          coordinate_radius=.03125, jet_radius=.015625):
    """Local Euler coefficients of the OWNED densities for independent fields.

    E_B = sum_alpha (-partial)^alpha partial L/partial(B_alpha).
    Mixed Tz is one independent raw jet, so it has no extra factor of two.
    Output gradients are per dT dz, not node derivatives or integrated probes.
    Scalar gauges are imposed only on the base data, AFTER variation.
    """
    varied = tuple(varied)
    if not varied or len(set(varied)) != len(varied) or any(v not in FIELDS for v in varied):
        raise ValueError("distinct named raw KS fields required")
    if isinstance(samples, bool) or not isinstance(samples, int) or samples < 8:
        raise ValueError("at least eight complex Cauchy samples required")
    if not np.isfinite([coordinate_radius, jet_radius]).all() or min(coordinate_radius, jet_radius) <= 0:
        raise ValueError("positive finite numerical differentiation radii required")
    coefficients = scalar_taylor_coefficients(fields)
    return _cauchy_euler(coefficients, light, compact, varied, samples,
                         np.longdouble(coordinate_radius), np.longdouble(jet_radius))


def incoming_local_constraints(domain, light, compact, *, samples=12, refined_samples=16,
                               coordinate_radius=.03125, jet_radius=.015625,
                               numerical_tolerance=3e-11):
    """Actual local lapse/shift forces on supplied SAME-surface normal data.

    This returns the local action part of both constraints. It cannot claim
    full constraints: the physical/reference source and band remainder are
    deliberately absent. Local force = minus the action Euler coefficient.
    """
    if not isinstance(domain, IncomingCauchyJets):
        raise ValueError("the same intrinsic incoming Cauchy-jet domain is required")
    domain.validate()
    if refined_samples <= samples or not np.isfinite(numerical_tolerance) or numerical_tolerance <= 0:
        raise ValueError("a finer Cauchy resolution and positive tolerance required")
    options = dict(coordinate_radius=coordinate_radius, jet_radius=jet_radius)
    base = local_euler_gradients(domain.fields, light, compact, samples=samples, **options)
    fine = local_euler_gradients(domain.fields, light, compact, samples=refined_samples, **options)
    gradients = {k: v[0].real.astype(float) for k, v in fine.items()}
    errors = {k: np.maximum(abs(fine[k][0]-base[k][0]), abs(fine[k][0].imag)).astype(float) for k in fine}
    active = tuple(k for k in gradients if k != "compact/euler_bulk_diagnostic")
    light_total = sum(v for k, v in gradients.items() if k.startswith("light/"))
    compact_total = sum(gradients[k] for k in active if k.startswith("compact/"))
    euler = float(np.max(abs(gradients["compact/euler_bulk_diagnostic"])))
    maximum_error = max(float(np.max(v)) for v in errors.values())
    total = light_total+compact_total
    return {
        "constraint_order": ("N", "beta"), "channel_action_gradients": gradients,
        "channel_forces": {k: -v for k, v in gradients.items()},
        "light_action_gradient": light_total, "compact_action_gradient": compact_total,
        "local_action_gradient": total, "local_force": -total,
        "channel_numerical_indicators": errors, "maximum_numerical_indicator": maximum_error,
        "Euler_bulk_identity_residual": euler, "numerical_tolerance": numerical_tolerance,
        "numerically_resolved": max(maximum_error, euler) <= numerical_tolerance,
        "differentiation": {"samples": samples, "refined_samples": refined_samples,
                            "coordinate_radius": coordinate_radius, "jet_radius": jet_radius,
                            "dtype": np.dtype(np.clongdouble).name},
        "scope": {"local_lapse_shift_forces_evaluated": True,
                  "physical_state_reference_or_band_included": False,
                  "light_geometry_included_once": True, "compact_action_included_once": True,
                  "physical_history_or_normal_data_selected": False,
                  "constraint_roots_or_four_metric_equations_solved": False,
                  "endpoint_forces_evaluated": False,
                  "boundary": "Euler and constant BoxR have no bulk force; endpoint variations remain separate",
                  "normalization": "delta S = integral dT dz sum_B E_B delta B; local_force=-E_B"},
    }
