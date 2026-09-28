"""Formal action-first reference phase from the existing NSC projector symbol.

The frame diagonalizes a reference band; it is NOT physical Cauchy transport,
V_c, or an EndpointBranchJets completion. The finite heat-scheme conversion
and momentum-boundary allocation remain explicit subsequent operations.
"""
import numpy as np

from .nsc_spatial_reference_symbol import SymbolJet, star_order, I, ZERO


def zero(n): return SymbolJet.constant(0).trim(4-n, 4-n)
def dagger(A): return SymbolJet(A.data.swapaxes(-1, -2).conj(), A.physical, A.momentum)
def series_dagger(A): return [dagger(a) for a in A]


def product(A, B, order=4):
    """All coefficients of the graded Weyl product, preserving matrix order."""
    result = []
    for n in range(order+1):
        value = zero(n)
        for i in range(min(n+1, len(A))):
            for j in range(min(n-i+1, len(B))):
                value += star_order(A[i], B[j], n-i-j)
        result.append(value)
    return result


def add(A, B): return [a+b for a, b in zip(A, B)]
def subtract(A, B): return [a-b for a, b in zip(A, B)]
def values(A): return np.stack(np.broadcast_arrays(*[a.value for a in A]))
def band_trace(A, target): return np.einsum('fij,nfji->nf', target, values(A))


def exact_massless_reference(projector, masses, angular, momenta):
    """Reuse the owned chiral LLL identity, never assign it to massive modes.

    In a gapped k-sign chart H is exactly diagonal, P is constant and all
    higher adiabatic coefficients vanish. This avoids differentiating the
    roundoff of a generic sqrt(h^2)/h cancellation. The stored certificate
    and physical covariance remain unchanged; the generic-point delta is
    reported explicitly.
    """
    m, l, k = map(np.asarray, (masses, angular, momenta))
    mask = (m == 0) & (l == 0)
    if np.any(k[mask] == 0): raise ValueError('gap-zero LLL point is outside this chart')
    copies = [SymbolJet(p.data, p.physical, p.momentum) for p in projector['jets']]
    maximum = 0.
    for j in np.flatnonzero(mask):
        Pi = np.diag([1., 0.]) if k[j] < 0 else np.diag([0., 1.])
        h = projector['hamiltonian'].data[:, j]
        commutator = h@Pi-Pi@h
        if np.max(abs(commutator)) > 3e-11:
            raise ValueError('LLL reference requires the exact inherited chiral operator')
        for n, p in enumerate(copies):
            expected = Pi if n == 0 else np.zeros((2, 2))
            maximum = max(maximum, float(np.max(abs(p.value[j]-expected))))
            p.data[:, j] = 0.
            p.data[ZERO, j] = expected
    result = dict(projector)
    result['jets'], result['orders'] = copies, values(copies)
    result['LLL_generic_point_delta'] = maximum
    result['LLL_chiral_operator_residual'] = 0.
    return result


def band_frame(projector, target):
    """Polar star-normalization of S=Pi P+(1-Pi)(1-P).

    Target Pi is a fixed local band chart. Its overlap with P0 must be
    positive. For the massless zero-angular sector the two gapped k-sign
    charts are kept separate; the gap-zero point is not bridged here.
    """
    P, H = projector['jets'], projector['hamiltonian']
    target = np.asarray(target, complex)
    if target.shape != P[0].value.shape or np.max(abs(target@target-target)) > 3e-11:
        raise ValueError('fixed rank-one target projector for every reference band required')
    if np.max(abs(target-target.swapaxes(-1, -2).conj())) > 3e-11 or np.max(abs(np.trace(target, axis1=-2, axis2=-1)-1)) > 3e-11:
        raise ValueError('Hermitian rank-one band target required')
    Pi = SymbolJet.constant(target)
    complement = [1-P[0]]+[-p for p in P[1:]]
    S = [Pi*p+(1-Pi)*q for p, q in zip(P, complement)]
    A = product(series_dagger(S), S)
    # The exact leading identity is S0^dagger S0=Tr(Pi P0) I.
    overlap = np.einsum('fij,cfji->cf', target, P[0].data)
    scalar = SymbolJet(overlap[..., None, None]*I)
    if np.any(scalar.value[:, 0, 0].real <= 1e-12):
        raise ValueError('zero-overlap band chart; no arbitrary normalization floor')
    root = scalar.power(.5)
    Z = [scalar.power(-.5)]
    for n in range(1, 5):
        residual = product(product(Z, A, n), Z, n)[n]
        Z.append((-residual/(2*root)).trim(4-n, 4-n))
    U = product(S, Z); Ud = series_dagger(U)
    UHUd = product(product(U, [H]), Ud)
    derivative = [u.derivative(t=1) for u in U[:4]]
    connection = [zero(0)]+[1j*x for x in product(derivative, Ud, 3)]
    effective = add(UHUd, connection)
    identity = [SymbolJet.constant(np.broadcast_to(I, target.shape))]+[zero(n) for n in range(1, 5)]
    wanted = [Pi]+[zero(n) for n in range(1, 5)]
    checks = {
        'unitarity_left': subtract(product(Ud, U), identity),
        'unitarity_right': subtract(product(U, Ud), identity),
        'projector_intertwining': subtract(product(product(U, P), Ud), wanted),
        'band_decoupling': [h*Pi-Pi*h for h in effective],
        'effective_Hermiticity': subtract(effective, series_dagger(effective)),
        'normalization_commutator': subtract(product(Z, A), product(A, Z)),
    }
    raw, scaled = {}, {}
    scales = 1+np.max(abs(values(effective)), axis=(-2, -1))
    for name, error in checks.items():
        raw[name] = np.max(abs(values(error)), axis=(-2, -1))
        # Frame identities normalized by their constructed series coefficients;
        # energy identities by the corresponding energy coefficients.
        denominator = scales if name in ('band_decoupling', 'effective_Hermiticity') else 1+np.max(abs(values(U)), axis=(-2, -1))**2
        scaled[name] = raw[name]/denominator
    leading = np.max(abs(A[0].data-scalar.data))
    return {'U': U, 'effective': effective, 'connection': connection, 'target': target,
            'band_energy': band_trace(effective, target), 'raw_residuals': raw,
            'scaled_residuals': scaled, 'leading_overlap_identity': float(leading),
            'interpretation': 'formal band-frame action density; no physical state or regulator assigned'}


def direction_series(rows, weights):
    """Differentiate a real metric family; never conjugate a complex-step metric."""
    return [sum(weight*row[n] for weight, row in zip(weights, rows)) for n in range(5)]


def action_variation_identity(base_projector, frame, delta_H, delta_U, delta_effective):
    """Retain the local time derivative and star-trace cyclicity defect.

    delta h = U star deltaH star U^dagger + [D,h]_star + i epsilon dT D,
    D=deltaU star U^dagger. Thus delta tr(Pi h) equals the reference vertex
    plus the explicit time and space/momentum terms. No finite-k cyclic trace
    is presumed and no term is dropped as an alleged interface counterforce.
    """
    U, Ud, target = frame['U'], series_dagger(frame['U']), frame['target']
    D = product(delta_U, Ud)
    push = product(product(U, [delta_H]), Ud)
    commutator = subtract(product(D, frame['effective']), product(frame['effective'], D))
    time = [zero(0)]+[1j*d.derivative(t=1) for d in D[:4]]
    rhs = add(add(push, commutator), time)
    matrix_residual = values(subtract(delta_effective, rhs))
    derivative = band_trace(delta_effective, target)
    physical_vertex = np.trace(values(product(base_projector['jets'], [delta_H])), axis1=-2, axis2=-1)
    temporal = band_trace(time, target)
    cyclic_defect = band_trace(push, target)-physical_vertex+band_trace(commutator, target)
    residual = derivative-physical_vertex-temporal-cyclic_defect
    denominator = 1+abs(derivative)+abs(physical_vertex)+abs(temporal)+abs(cyclic_defect)
    return {
        'action_derivative_density': derivative, 'projector_vertex_density': physical_vertex,
        'temporal_derivative_density': temporal, 'star_trace_defect_density': cyclic_defect,
        'identity_absolute': abs(residual), 'identity_scaled': abs(residual)/denominator,
        'conjugation_absolute': np.max(abs(matrix_residual), axis=(-2, -1)),
        'conjugation_scaled': np.max(abs(matrix_residual), axis=(-2, -1))/(1+np.max(abs(values(delta_effective)), axis=(-2, -1))+np.max(abs(values(rhs)), axis=(-2, -1))),
    }
