"""Identifiability of finite covariant terms on a smooth static axial cell.

The coefficients are unknown properties of ONE completed measure, not weights
added to an independently counted gravitational action. The static functional
F=4*pi*int N*q*r^2 I dx defines the coefficient/sign convention; S_fin=-int F dt.
Curvature uses signature (+---) and R^a_bcd=partial_c Gamma^a_db-... .
No quantum state, absolute Dirac stress or stationary geometry is supplied here.
"""
from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
from itertools import product

import numpy as np
import sympy as sp

from .nsc_shape_response import StaticAxialMetric, smooth_metric

BASIS = ("M4", "M2R", "C2", "R2", "E4", "boxR")
FIELDS = ("lapse", "radial_scale", "sphere_radius")


@lru_cache(maxsize=1)
def _symbolic_basis():
    jets = tuple(sp.symbols(name + "0:3", positive=True) for name in ("N", "q", "r"))
    (n, nx, nxx), (q, qx, _), (r, rx, rxx) = jets
    mass = sp.Symbol("M", positive=True)
    a = (nxx - nx * qx / q) / (n * q**2)
    b = nx * rx / (n * q**2 * r)
    c = (rxx - rx * qx / q) / (q**2 * r)
    d = (1 - rx**2 / q**2) / r**2
    scalar = 2 * (a + 2*b + 2*c - d)
    ricci = (a + 2*b)**2 + (a + 2*c)**2 + 2*(b + c - d)**2
    riemann = 4*(a*a + 2*b*b + 2*c*c + d*d)
    invariants = (mass**4, mass**2*scalar,
                  sp.Rational(4, 3)*(a-b-c-d)**2, scalar**2,
                  8*(-a*d + 2*b*c))
    densities = tuple(n*q*r*r*value for value in invariants)
    return jets, mass, (a, b, c, d), scalar, ricci, riemann, densities


def _total_derivative(expression, jets):
    return sum(sp.diff(expression, row[j])*row[j+1] for row in jets for j in range(2))


@lru_cache(maxsize=1)
def exact_identities():
    """Jet-algebra curvature, total derivative, local Weyl and unit identities."""
    jets, mass, sections, scalar, ricci, riemann, densities = _symbolic_basis()
    (n, nx, _), (q, _, _), (r, rx, _) = jets
    sigma = sp.symbols("sigma0:3")
    variations = tuple((sigma[0]*row[0], sigma[0]*row[1]+sigma[1]*row[0],
                        sigma[0]*row[2]+2*sigma[1]*row[1]+sigma[2]*row[0]) for row in jets)
    weyl = lambda expression: sp.expand(sum(sp.diff(expression, z)*v
        for row, delta in zip(jets, variations) for z, v in zip(row, delta)))
    primitive = -8*nx/q*(1-rx**2/q**2)
    a, b, c, d = sections
    rigid = []
    # Fixed coordinate interval: N,q,r all acquire s; M fixed. Density
    # weights are 4,2,0,0,0, including the time-coordinate normalization.
    for density, weight in zip(densities, (4, 2, 0, 0, 0)):
        rigid.append(str(sp.simplify(weyl(density).subs({sigma[1]: 0, sigma[2]: 0})
                                     -weight*sigma[0]*density)))
    return {
        "proper_radial_derivative": "D=(1/q)*partial_x",
        "sections": {"A": "D(D(N))/N", "B": "D(N)*D(r)/(N*r)",
                     "C": "D(D(r))/r", "D": "(1-D(r)^2)/r^2"},
        "R": "2*(A+2*B+2*C-D)",
        "Ricci2": "(A+2*B)^2+(A+2*C)^2+2*(B+C-D)^2",
        "Riemann2": "4*(A^2+2*B^2+2*C^2+D^2)",
        "C2": "4*(A-B-C-D)^2/3", "E4": "8*(-A*D+2*B*C)",
        "weyl_contraction_residual": str(sp.simplify(riemann-2*ricci+scalar**2/3
                                                      -sp.Rational(4, 3)*(a-b-c-d)**2)),
        "euler_contraction_residual": str(sp.simplify(riemann-4*ricci+scalar**2
                                                       -8*(-a*d+2*b*c))),
        "euler_density_primitive": "N*q*r^2*E4=partial_x[-8*N'/q*(1-r'^2/q^2)]",
        "euler_primitive_residual": str(sp.simplify(densities[4]-_total_derivative(primitive, jets))),
        "box_density_primitive": "N*q*r^2*box(R)=-partial_x[N*r^2*R'/q]",
        "pointwise_C2_local_weyl_residual": str(sp.simplify(weyl(densities[2]))),
        "rigid_weyl_density_residuals": rigid,
        "rigid_weyl_density_weights_at_fixed_M": [4, 2, 0, 0, 0, 0],
        "coefficient_matching_not_done": True,
    }


@lru_cache(maxsize=1)
def direct_curvature_identities():
    """Independent Christoffel/Riemann contractions of the full 4D metric.

    No sectional-curvature formula is used in constructing the connection or
    contractions. The symbolic comparison retains arbitrary N(x), q(x), r(x).
    """
    t, x, theta, phi = sp.symbols("t x theta phi", real=True)
    coordinates = (t, x, theta, phi)
    n, q, r = (sp.Function(name)(x) for name in ("N", "q", "r"))
    diagonal = (n*n, -q*q, -r*r, -r*r*sp.sin(theta)**2)
    connection = {}
    for i, j, k in product(range(4), repeat=3):
        value = ((sp.diff(diagonal[i], coordinates[j]) if i == k else 0)
                 +(sp.diff(diagonal[i], coordinates[k]) if i == j else 0)
                 -(sp.diff(diagonal[j], coordinates[i]) if j == k else 0))/(2*diagonal[i])
        if value != 0:
            connection[i, j, k] = sp.simplify(value)
    gamma = lambda i, j, k: connection.get((i, j, k), sp.S.Zero)
    tensor = {}
    for a, b, c, d in product(range(4), repeat=4):
        value = sp.diff(gamma(a, d, b), coordinates[c])-sp.diff(gamma(a, c, b), coordinates[d])
        value += sum(gamma(a, c, e)*gamma(e, d, b)-gamma(a, d, e)*gamma(e, c, b) for e in range(4))
        if value != 0:
            value = sp.simplify(value)
            if value != 0:
                tensor[a, b, c, d] = value
    ricci_tensor = {(b, d): sp.simplify(sum(tensor.get((a, b, a, d), 0) for a in range(4)))
                    for b, d in product(range(4), repeat=2)}
    scalar = sp.simplify(sum(ricci_tensor[i, i]/diagonal[i] for i in range(4)))
    ricci = sp.simplify(sum(value**2/(diagonal[i]*diagonal[j])
                            for (i, j), value in ricci_tensor.items()))
    riemann = sp.simplify(sum(diagonal[a]*value**2/(diagonal[b]*diagonal[c]*diagonal[d])
                              for (a, b, c, d), value in tensor.items()))
    jets, _, _, expected_scalar, expected_ricci, expected_riemann, _ = _symbolic_basis()
    substitutions = {jet: sp.diff(field, x, order)
                     for row, field in zip(jets, (n, q, r)) for order, jet in enumerate(row)}
    return {name: str(sp.simplify(observed-expected.subs(substitutions)))
            for name, observed, expected in (("R_residual", scalar, expected_scalar),
                                             ("Ricci2_residual", ricci, expected_ricci),
                                             ("Riemann2_residual", riemann, expected_riemann))}


@lru_cache(maxsize=1)
def exact_throat_sensitivities():
    """Exact local stress responses, obtained after general metric variation."""
    jets, mass, _, _, _, _, densities = _symbolic_basis()
    extended = tuple(row+sp.symbols(name+"3:5") for row, name in zip(jets, ("N", "q", "r")))
    derivative = lambda f: sum(sp.diff(f, row[j])*row[j+1] for row in extended for j in range(4))
    radius, wave = sp.symbols("a k", positive=True)
    values = ((1, 0, 0, 0, 0), (1, 0, 0, 0, 0),
              (radius, 0, 1/radius, 0, -3/radius**3-4*wave**2/radius))
    substitutions = {symbol: value for row, vals in zip(extended, values) for symbol, value in zip(row, vals)}
    output = {}
    for basis, density in zip(BASIS, densities):
        el = [sp.simplify((sp.diff(density, row[0])-derivative(sp.diff(density, row[1]))
                          +derivative(derivative(sp.diff(density, row[2])))).subs(substitutions)) for row in jets]
        rho, px, pt = el[0]/radius**2, -el[1]/radius**2, -el[2]/(2*radius)
        output[basis] = {name: str(sp.factor(value)) for name, value in
                         (("rho", rho), ("p_x", px), ("p_perp", pt), ("null", rho+px))}
    output["boxR"] = dict.fromkeys(("rho", "p_x", "p_perp", "null"), "0")
    return output


def periodic_derivative(values, length, order=1):
    """Real Fourier differentiation on an even periodic grid.

    Both Nyquist multipliers are zero, so the real first derivative is
    antisymmetric and its square is the symmetric second derivative.
    These adjoints, rather than a continuum product rule assumed on a lattice,
    define the exact gradient of the discretized integral.
    """
    values = np.asarray(values)
    if order not in (1, 2):
        raise ValueError("only first and second derivatives are required")
    frequency = 2*np.pi*np.fft.fftfreq(values.shape[-1], length/values.shape[-1])
    multiplier = (1j*frequency)**order
    multiplier[len(frequency)//2] = 0
    return np.fft.ifft(np.fft.fft(values, axis=-1)*multiplier, axis=-1).real


@lru_cache(maxsize=1)
def _density_functions():
    jets, mass, _, scalar, _, _, densities = _symbolic_basis()
    flat = tuple(z for row in jets for z in row)+(mass,)
    functions = []
    for density in densities:
        derivatives = [sp.diff(density, z) for row in jets for z in row]
        functions.append(sp.lambdify(flat, [density]+derivatives, "numpy", cse=True))
    return functions, sp.lambdify(flat, scalar, "numpy", cse=True)


def finite_response(metric, normalization_mass=1.):
    """All six finite terms' energies and N,q,r Euler-Lagrange densities.

    The final two terms have zero integrated bulk response on the closed cell.
    E4 is nevertheless differentiated before cancellation as a numerical
    control; boxR's gradient follows from its exact total-derivative identity.
    Gradients are densities per dx, not derivatives per lattice node.
    """
    if not np.isfinite(normalization_mass) or normalization_mass <= 0:
        raise ValueError("positive finite normalization mass required")
    values = []
    for field in FIELDS:
        f = getattr(metric, field)
        values += [f, periodic_derivative(f, metric.length), periodic_derivative(f, metric.length, 2)]
    functions, scalar_function = _density_functions()
    evaluated = [np.stack([np.broadcast_to(a, (metric.points,)) for a in function(*values, normalization_mass)])
                 for function in functions]
    densities = np.stack([item[0] for item in evaluated])
    gradients = np.empty((5, 3, metric.points))
    for i, item in enumerate(evaluated):
        for j in range(3):
            p0, p1, p2 = item[1+3*j:4+3*j]
            gradients[i, j] = 4*np.pi*(p0-periodic_derivative(p1, metric.length)
                                        +periodic_derivative(p2, metric.length, 2))
    n, q, r = (getattr(metric, field) for field in FIELDS)
    scalar = np.broadcast_to(scalar_function(*values, normalization_mass), (metric.points,))
    box_density = -periodic_derivative(n*r*r/q*periodic_derivative(scalar, metric.length), metric.length)
    densities = np.vstack((densities, box_density))
    gradients = np.concatenate((gradients, np.zeros((1, 3, metric.points))), axis=0)
    rho = gradients[:, 0]/(4*np.pi*q*r*r)
    px = -gradients[:, 1]/(4*np.pi*n*r*r)
    pt = -gradients[:, 2]/(8*np.pi*n*q*r)
    return {"energy": 4*np.pi*metric.spacing*np.sum(densities, axis=1),
            "densities": densities, "gradients": gradients, "R": scalar,
            "rho": rho, "p_x": px, "p_perp": pt, "null": rho+px,
            "weyl_density": n*gradients[:, 0]+q*gradients[:, 1]+r*gradients[:, 2],
            "box_density": box_density}


def independent_energy(metric, normalization_mass=1.):
    """Scalar energy path, using contracted curvature instead of jet gradients."""
    n, q, r = (getattr(metric, field) for field in FIELDS)
    dx = lambda f: periodic_derivative(f, metric.length)
    dxx = lambda f: periodic_derivative(f, metric.length, 2)
    a = (dxx(n)-dx(n)*dx(q)/q)/(n*q*q)
    b = dx(n)*dx(r)/(n*q*q*r)
    c = (dxx(r)-dx(r)*dx(q)/q)/(q*q*r)
    d = (1-dx(r)**2/q**2)/r**2
    scalar = 2*(a+2*b+2*c-d)
    ricci = (a+2*b)**2+(a+2*c)**2+2*(b+c-d)**2
    kretschmann = 4*(a*a+2*b*b+2*c*c+d*d)
    invariants = np.stack((np.full(metric.points, normalization_mass**4), normalization_mass**2*scalar,
                           kretschmann-2*ricci+scalar**2/3, scalar**2,
                           kretschmann-4*ricci+scalar**2, -dx(n*r*r/q*dx(scalar))/(n*q*r*r)))
    return 4*np.pi*metric.spacing*np.sum(invariants*n*q*r*r, axis=1)


def response_probes(metric):
    """Independent clock, areal-radius and spatial-volume-preserving shapes."""
    result = []
    for sector in ("lapse", "radius", "volume_preserving_shape"):
        for harmonic in range(3):
            f = np.cos(2*np.pi*harmonic*metric.x/metric.length)
            direction = np.zeros((3, metric.points))
            if sector == "lapse":
                direction[0] = metric.lapse*f
            elif sector == "radius":
                direction[2] = metric.sphere_radius*f
            else:
                direction[1] = metric.radial_scale*f
                direction[2] = -.5*metric.sphere_radius*f
            result.append((f"{sector}_cos{harmonic}", direction))
    return result


def response_matrix(metric, response=None):
    response = finite_response(metric) if response is None else response
    probes = response_probes(metric)
    matrix = np.array([metric.spacing*np.einsum("if,kif->k", delta, response["gradients"])
                       for _, delta in probes])
    norms = np.linalg.norm(matrix[:, :4], axis=0)
    scales = np.r_[norms, [1., 1.]]
    scaled = matrix/scales
    singular = np.linalg.svd(scaled, compute_uv=False)
    threshold = 1e-8
    return {"probe_names": [name for name, _ in probes], "basis": list(BASIS),
            "matrix": matrix.tolist(), "column_scales": scales.tolist(),
            "column_normalized_singular_values": singular.tolist(),
            "absolute_singular_value_rank_threshold": threshold,
            "rank": int(np.sum(singular > threshold)),
            "bulk_column_normalized_condition_number": float(singular[0]/singular[3]) if singular[3] > threshold else None}


def general_control_metric(points=128):
    metric = smooth_metric(2., points)
    angle = 2*np.pi*metric.x/metric.length
    return replace(metric, lapse=np.exp(.08*np.cos(angle)+.03*np.sin(2*angle)),
                   radial_scale=np.exp(.06*np.sin(angle)-.02*np.cos(2*angle)))


def ward_controls(metric, response=None):
    response = finite_response(metric) if response is None else response
    gradient = response["gradients"]
    n, q, r = (getattr(metric, field) for field in FIELDS)
    dx = lambda f: periodic_derivative(f, metric.length)
    diffeomorphism = dx(n)*gradient[:, 0]+dx(r)*gradient[:, 2]-q*dx(gradient[:, 1])
    expected_weyl = np.array([16*np.pi*response["densities"][0],
                              8*np.pi*response["densities"][1], np.zeros(metric.points),
                              -48*np.pi*response["box_density"], np.zeros(metric.points), np.zeros(metric.points)])
    xi = .2*np.sin(2*np.pi*metric.x/metric.length)+.05*np.cos(4*np.pi*metric.x/metric.length)
    gauge_direction = np.array([xi*dx(n), xi*dx(q)+q*dx(xi), xi*dx(r)])
    projection = metric.spacing*np.einsum("if,kif->k", gauge_direction, gradient)
    return {"maximum_local_radial_diffeomorphism_residual_by_basis": np.max(np.abs(diffeomorphism), axis=1).tolist(),
            "integrated_radial_gauge_response_by_basis": projection.tolist(),
            "maximum_local_weyl_residual_by_basis": np.max(np.abs(response["weyl_density"]-expected_weyl), axis=1).tolist(),
            "global_lapse_homogeneity_residual_by_basis": (metric.spacing*np.sum(n*gradient[:, 0], axis=1)-response["energy"]).tolist()}
