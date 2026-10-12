"""Common near-diagonal KS pairing for the selected canonical subtraction.

Physical source frequencies and reference Weyl momenta label different
representations. They are integrated separately into the SAME equal-KS-time
kernel before pairing with the first-order metric vertex. A finite numerical
separation is a convergence device, not a new physical cutoff prescription.
This module neither evolves a sampled-frequency surrogate nor supplies the
missing physical spectral panels, local allocations or an absolute stress.
"""
from dataclasses import dataclass

import numpy as np

from .nsc_reference_band_action import action_variation_identity, product, values


_S1 = np.array([[0., 1.], [1., 0.]], complex)
_S2 = np.array([[0., -1j], [1j, 0.]], complex)
_S3 = np.diag([1., -1.]).astype(complex)
_I = np.eye(2, dtype=complex)


@dataclass(frozen=True)
class SpectrumDeclaration:
    """Explicit coverage of supplied numerical data, never inferred from nodes.

    ``complete_spectrum`` requires named convergence evidence. The declaration
    is a caller's evidence claim; a finite array alone cannot establish UV
    completeness. Control or partial data remain useful for kernel checks.
    """
    scope: str
    provenance: str
    convergence_evidence: str = ""

    def __post_init__(self):
        if self.scope not in ("control_only", "partial_spectrum", "complete_spectrum"):
            raise ValueError("explicit control_only, partial_spectrum or complete_spectrum scope required")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("resolved spectral input provenance required")
        if not isinstance(self.convergence_evidence, str):
            raise ValueError("spectral convergence evidence must be named text")
        if self.scope == "complete_spectrum" and not self.convergence_evidence.strip():
            raise ValueError("complete spectrum requires explicit spectral convergence evidence")


@dataclass(frozen=True)
class KSKernelSample:
    """Canonical 2x2 kernel and its separation derivative at one common split."""
    kernel: np.ndarray
    separation_derivative: np.ndarray
    separation: np.ndarray
    declaration: SpectrumDeclaration
    origin: str


def _quadrature(nodes, weights, name):
    if np.iscomplexobj(nodes) or np.iscomplexobj(weights):
        raise ValueError(f"real {name} nodes and quadrature weights required")
    x, w = np.asarray(nodes, float), np.asarray(weights, float)
    if (x.ndim != 1 or not len(x) or w.shape != x.shape
            or not np.isfinite(x).all() or not np.isfinite(w).all() or np.any(w <= 0)):
        raise ValueError(f"finite {name} nodes and positive quadrature weights required")
    return x, w


def _split(separation, shape):
    if np.iscomplexobj(separation):
        raise ValueError("real equal-time spatial separation required")
    try:
        eta = np.broadcast_to(np.asarray(separation, float), shape)
    except ValueError as error:
        raise ValueError("separation must broadcast to the common kernel sample shape") from error
    if not np.isfinite(eta).all():
        raise ValueError("finite equal-time spatial separation required")
    return eta


def physical_source_kernel(plus, minus, plus_z, minus_z, source_covariance,
                           frequencies, weights, *, separation, declaration):
    """Restrict resolved three-source fields to equal KS time at z+/-eta/2.

    Fields have shape (..., nE, 2, 3), covariance (nE, 3, 3), and source
    quadrature dE/(2*pi). The columns contain their existing source phases.
    The same canonical Pauli frame is used at both endpoints. Two-component
    seed covariances and integrated overlap matrices are not field columns.
    No relation between incoming E and instantaneous local k is assumed.
    """
    if not isinstance(declaration, SpectrumDeclaration):
        raise ValueError("physical spectrum completeness declaration required")
    E, w = _quadrature(frequencies, weights, "incoming source-frequency")
    fields = [np.asarray(p, complex) for p in (plus, minus, plus_z, minus_z)]
    expected = (len(E), 2, 3)
    if (fields[0].ndim < 3 or fields[0].shape[-3:] != expected
            or any(p.shape != fields[0].shape or not np.isfinite(p).all() for p in fields)):
        raise ValueError("resolved (..., nE, 2, 3) KS source columns and spatial derivatives required")
    C = np.asarray(source_covariance, complex)
    if (C.shape != (len(E), 3, 3) or not np.isfinite(C).all()
            or np.max(abs(C-C.swapaxes(-1, -2).conj())) > 3e-11):
        raise ValueError("the inherited Hermitian three-source covariance is required")
    eigenvalues = np.linalg.eigvalsh(C)
    if eigenvalues.min() < -3e-11 or eigenvalues.max() > 1+3e-11:
        raise ValueError("physical source covariance must satisfy CAR")
    p, m, pz, mz = fields
    def contract(left, right):
        return np.einsum('...eis,est,...ejt,e->...ij', left, C, right.conj(), w/(2*np.pi), optimize=True)
    matrix = contract(p, m)
    derivative = .5*(contract(pz, m)-contract(p, mz))
    eta = _split(separation, matrix.shape[:-2])
    return KSKernelSample(matrix, derivative, eta, declaration, "physical_source")


def reference_symbol_kernel(projector_orders, momenta, weights, *, separation, declaration):
    """Fourier/Weyl kernel of exactly the retained orders zero through four.

    ``projector_orders`` has shape (5, ..., nk, 2, 2). This formal symbol is
    not a physical Gaussian covariance; no pointwise CAR constraint is added.
    The k weights are numerical Fourier quadrature, not incoming E weights or
    a declaration of a new physical regulator. A finite reference sum does
    not by itself establish the distributional UV/coincidence limit.
    """
    if not isinstance(declaration, SpectrumDeclaration):
        raise ValueError("reference spectrum completeness declaration required")
    k, w = _quadrature(momenta, weights, "local Weyl-momentum")
    P = np.asarray(projector_orders, complex)
    if (P.ndim < 4 or P.shape[0] != 5 or P.shape[-3:] != (len(k), 2, 2)
            or not np.isfinite(P).all()):
        raise ValueError("all five finite order<=4 reference symbol arrays required")
    if np.max(abs(P-P.swapaxes(-1, -2).conj())) > 3e-10:
        raise ValueError("Hermitian formal reference symbol orders required")
    eta = _split(separation, P.shape[1:-3])
    phase = np.exp(1j*eta[..., None]*k)*w/(2*np.pi)
    symbol = np.sum(P, axis=0)
    matrix = np.einsum('...k,...kij->...ij', phase, symbol, optimize=True)
    derivative = np.einsum('...k,...kij->...ij', 1j*k*phase, symbol, optimize=True)
    return KSKernelSample(matrix, derivative, eta, declaration, "reference_symbol")


@dataclass(frozen=True)
class KSVertexCoefficients:
    """Multiplication M_B and momentum V_B in delta H=OpW(M_B+V_B*k)."""
    multiplication: np.ndarray
    momentum: np.ndarray


def raw_ks_vertex_coefficients(metric, mass, angular, *, envelopes):
    """Actual endpoint vertices, ordered (N, beta, a, r).

    ``metric`` ends in four real entries (N,beta,a,r). ``envelopes``
    broadcasts to that same shape and contains each actual raw-field
    variation at this endpoint. Call separately at plus and minus endpoints.
    The Weyl symbol already includes symmetric differential ordering.
    """
    if np.iscomplexobj(metric) or np.iscomplexobj(envelopes):
        raise ValueError("real physical metric and variation envelopes required")
    g = np.asarray(metric, float)
    if (g.ndim < 1 or g.shape[-1] != 4 or not np.isfinite(g).all()
            or np.any(g[..., [0, 2, 3]] <= 0)
            or not np.isfinite([mass, angular]).all()):
        raise ValueError("finite raw KS metric with positive N,a,r and finite channel labels required")
    try:
        f = np.broadcast_to(np.asarray(envelopes, float), g.shape)
    except ValueError as error:
        raise ValueError("four endpoint variation envelopes must broadcast to the metric") from error
    if not np.isfinite(f).all():
        raise ValueError("finite endpoint variation envelopes required")
    N, _, a, r = np.moveaxis(g, -1, 0)
    M = np.zeros((*g.shape, 2, 2), complex); V = np.zeros_like(M)
    M[..., 0, :, :] = -mass*_S1+angular/r[..., None, None]*_S2
    M[..., 3, :, :] = -(N*angular/r**2)[..., None, None]*_S2
    V[..., 0, :, :] = _S3/a[..., None, None]
    V[..., 1, :, :] = -_I
    V[..., 2, :, :] = -(N/a**2)[..., None, None]*_S3
    return KSVertexCoefficients(M*f[..., None, None], V*f[..., None, None])


@dataclass(frozen=True)
class SymmetricKSPairing:
    """Subtracted canonical action density; band/local terms remain separate."""
    vertex_density: np.ndarray
    action_gradient_density: np.ndarray
    separation: np.ndarray
    physical_declaration: SpectrumDeclaration
    reference_declaration: SpectrumDeclaration

    def integrate(self, spacetime_weights):
        """Integrate with supplied dT dz weights, dPGtime drho/a0 if in PG.

        No actual-a measure or 4*pi*r^2 projection is inserted into this
        canonical action trace. Angular/compact multiplicities belong to the
        common channel inventory and must be restored once by the caller.
        """
        if np.iscomplexobj(spacetime_weights):
            raise ValueError("real canonical spacetime quadrature weights required")
        w = np.asarray(spacetime_weights, float)
        if (w.shape != self.vertex_density.shape[:-1] or not np.isfinite(w).all()
                or np.any(w <= 0)):
            raise ValueError("positive spacetime weights on the exact common sample grid required")
        axes = tuple(range(w.ndim))
        return np.sum(w[..., None]*self.action_gradient_density, axis=axes)

    def require_complete_trace(self, *, coincidence_convergence_evidence, channel_inventory_evidence):
        """Reject a complete-source promotion of retained controls or panels.

        Passing validates explicit declarations for this canonical pairing;
        it does not assert that the separate band/local terms were assembled
        or that the resulting tensor solves the metric equations.
        """
        if any(d.scope != "complete_spectrum" for d in
               (self.physical_declaration, self.reference_declaration)):
            raise ValueError("control-only or partial spectra cannot support a complete source/trace claim")
        for evidence in (coincidence_convergence_evidence, channel_inventory_evidence):
            if not isinstance(evidence, str) or not evidence.strip():
                raise ValueError("named coincidence convergence and complete channel inventory evidence required")


def symmetric_pairing_density(physical, reference, plus_vertices, minus_vertices, *,
                              compact_spatial_boundary):
    """Pair {delta H,C-P}/2 at a common split, keeping the Gaussian minus sign.

    After spatial integration the exact bulk density is
    Re tr(Mbar D_eta - i Vbar partial_eta D_eta). The omitted boundary is
    Re[-i tr((Vplus-Vminus)D_eta)/4] at the spatial ends. This API accepts
    only compact boundary-identical vertices and one constant split across
    the integration grid. A spatially varying split would introduce extra
    derivatives. Reference time/star terms are NOT removed.
    """
    if compact_spatial_boundary is not True:
        raise ValueError("noncompact spatial vertices require explicit endpoint terms; omitted endpoints are not supported")
    if (not isinstance(physical, KSKernelSample) or not isinstance(reference, KSKernelSample)
            or physical.origin != "physical_source" or reference.origin != "reference_symbol"):
        raise ValueError("resolved physical and Weyl reference kernels required, not a seed or supplied force")
    if (physical.kernel.shape != reference.kernel.shape
            or not np.array_equal(physical.separation, reference.separation)):
        raise ValueError("state/reference kernels must share the same KS sample grid and separation")
    if physical.separation.size and not np.all(physical.separation == physical.separation.flat[0]):
        raise ValueError("one constant separation across the trace grid is required; variable-split derivative terms are not omitted")
    shape = (*physical.kernel.shape[:-2], 4, 2, 2)
    for vertices in (plus_vertices, minus_vertices):
        if (not isinstance(vertices, KSVertexCoefficients)
                or vertices.multiplication.shape != shape or vertices.momentum.shape != shape):
            raise ValueError("all four actual endpoint KS vertex coefficients on the common grid required")
    M = .5*(plus_vertices.multiplication+minus_vertices.multiplication)
    V = .5*(plus_vertices.momentum+minus_vertices.momentum)
    D = physical.kernel-reference.kernel
    Deta = physical.separation_derivative-reference.separation_derivative
    vertex = np.einsum('...bij,...ji->...b', M, D)-1j*np.einsum('...bij,...ji->...b', V, Deta)
    return SymmetricKSPairing(vertex.real, -vertex.real, physical.separation,
                              physical.declaration, reference.declaration)


@dataclass(frozen=True)
class ReferenceBandRemainder:
    """Actual band variation minus the symmetric projector/vertex insertion."""
    orders: np.ndarray
    symmetric_vertex_orders: np.ndarray
    action_derivative_orders: np.ndarray
    temporal_orders: np.ndarray
    star_exchange_orders: np.ndarray
    antisymmetric_exchange_orders: np.ndarray
    decomposition_residual: np.ndarray
    band_identity_residual: np.ndarray

    def fourier_pairing(self, momenta, weights, *, separation):
        """Same finite split and dk/(2*pi); retain orders <=4 before summing.

        No finite-k cyclicity or boundary cancellation is presumed. The
        caller integrates this positive subtraction-action remainder in
        dT dz and adds it to the negative common vertex action gradient.
        """
        k, w = _quadrature(momenta, weights, "reference remainder Weyl-momentum")
        if self.orders.ndim != 2 or self.orders.shape != (5, len(k)):
            raise ValueError("one local order<=4 band remainder on the supplied momentum nodes required")
        eta = np.asarray(separation)
        eta = _split(eta, eta.shape)
        phase = np.exp(1j*eta[..., None]*k)*w/(2*np.pi)
        return np.einsum('...k,k->...', phase, np.sum(self.orders, axis=0)).real


def reference_band_remainder(base_projector, frame, delta_H, delta_U, delta_effective):
    """Reuse the exact formal band-action identity, including its star terms.

    Inputs are the owned SymbolJet projector/frame and their actual variation,
    not an arbitrary supplied correction. If a_B=delta tr(Pi*h), this returns
    a_B - tr(P star deltaH + deltaH star P)/2 through formal order four.
    """
    if (len(base_projector['jets']) != 5 or len(delta_U) != 5
            or len(delta_effective) != 5):
        raise ValueError("all five owned formal projector/frame variation orders required")
    identity = action_variation_identity(base_projector, frame, delta_H, delta_U, delta_effective)
    left = identity['projector_vertex_density']
    right = np.trace(values(product([delta_H], base_projector['jets'])), axis1=-2, axis2=-1)
    symmetric = .5*(left+right)
    antisymmetric = .5*(left-right)
    remainder = identity['action_derivative_density']-symmetric
    decomposed = (identity['temporal_derivative_density']
                  +identity['star_trace_defect_density']+antisymmetric)
    return ReferenceBandRemainder(
        remainder, symmetric, identity['action_derivative_density'],
        identity['temporal_derivative_density'], identity['star_trace_defect_density'],
        antisymmetric, abs(remainder-decomposed), identity['identity_absolute'])
