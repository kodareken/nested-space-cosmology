"""Conditional normal geometric data on the existing rho=1 KS Cauchy slice.

The intrinsic metric, lapse, shift, normal and canonical spinor frame stay
fixed. Only normal derivatives of a and r, including mixed spatial jets,
are free. These local Taylor data select neither a history nor a duration.
The same abstract canonical C0 can therefore be held fixed on this surface;
its Hadamard regularity and origin in a new smooth parent history do not
follow from that identification.

The UV screen below changes the EXISTING fourth-order reference and pairs
its difference with all four instantaneous metric vertices. It does not
integrate a state, solve constraints, or add the separately owned LLL local
allocation. Finite momentum samples are not an all-momentum convergence proof.
"""
from dataclasses import dataclass
from math import factorial

import numpy as np

from .nsc_spatial_reference_symbol import (
    SymbolJet, background_scalar, reference_projector, INDEX, INDICES, I, SIGMA,
)
from .nsc_reference_band_action import exact_massless_reference, product, values


FIELDS = ("N", "beta", "a", "r")
NORMAL_ENTRIES = tuple((t, z) for t in range(1, 5) for z in range(5-t))


@dataclass(frozen=True)
class IncomingNormalJetChange:
    """Change in the derivative dT^normal_order dz^spatial_order of a or r.

    ``delta`` is a derivative value, not a Taylor coefficient. Division by
    both factorials occurs once when the SymbolJet is constructed.
    """
    field: str
    normal_order: int
    spatial_order: int
    delta: float

    def __post_init__(self):
        if self.field not in ("a", "r"):
            raise ValueError("only normal a/r data are free; lapse, shift and frame remain fixed")
        if (isinstance(self.normal_order, bool) or isinstance(self.spatial_order, bool)
                or not isinstance(self.normal_order, (int, np.integer))
                or not isinstance(self.spatial_order, (int, np.integer))
                or (self.normal_order, self.spatial_order) not in NORMAL_ENTRIES):
            raise ValueError("normal order >=1 and total T/z order <=4 are required")
        if np.ndim(self.delta) != 0 or np.iscomplexobj(self.delta) or not np.isfinite(self.delta):
            raise ValueError("finite real conditional normal-jet change required")


def reference_incoming_metric_jets():
    """Reuse the original chart's d rho/dT=-a0 at rho=1 through order four.

    T=0 is only the local Taylor origin on this SAME Cauchy surface. It is
    not the old neck slice and does not assign an elapsed physical time.
    """
    rho = SymbolJet.constant(1.)
    for n in range(4):
        _, axial, _ = background_scalar(rho)
        rho.data[INDEX[n+1, 0, 0]] = -axial.data[INDEX[n, 0, 0]]/(n+1)
    _, axial, radius = background_scalar(rho)
    return (SymbolJet.constant(1.), SymbolJet.constant(0.), axial, radius)


@dataclass(frozen=True)
class IncomingCauchyJets:
    baseline: tuple
    fields: tuple
    changes: tuple

    def validate(self):
        """Reject hidden intrinsic, lapse/shift, momentum or frame changes."""
        if len(self.baseline) != 4 or len(self.fields) != 4:
            raise ValueError("the four raw KS metric jets are required")
        owned = reference_incoming_metric_jets()
        for name, original, field, expected in zip(FIELDS, self.baseline, self.fields, owned):
            if (not isinstance(original, SymbolJet) or original.data.shape != expected.data.shape
                    or not np.isfinite(original.data).all()
                    or np.max(abs(original.data-expected.data)) > 3e-12):
                raise ValueError("the baseline must be the existing rho=1 reference-chart jets")
            if (not isinstance(field, SymbolJet) or field.data.shape[1:] != (1, 2, 2)
                    or field.physical != 4 or field.momentum != 4
                    or not np.isfinite(field.data).all()):
                raise ValueError("finite order-four local scalar metric jets required")
            scalar = field.data[..., :1, :1]*I
            if max(np.max(abs(field.data-scalar)), np.max(abs(field.data.imag))) > 3e-12:
                raise ValueError("physical incoming metric jets must be real scalars")
            difference = field.data-original.data
            for i, (t, z, k) in enumerate(INDICES):
                if k and np.max(abs(field.data[i])) > 3e-12:
                    raise ValueError("metric jets cannot depend on reference momentum")
                if (name in ("N", "beta") or t == 0) and np.max(abs(difference[i])) > 3e-12:
                    raise ValueError("incoming intrinsic metric, lapse/shift and normal/frame must remain unchanged")
        if any(self.fields[i].value[0, 0, 0].real <= 0 for i in (0, 2, 3)):
            raise ValueError("positive intrinsic N,a,r required")

    def changed_normal_entries(self):
        """Explicit before/after derivative values for the altered local data."""
        self.validate()
        rows = []
        for name in ("a", "r"):
            base, field = self.baseline[FIELDS.index(name)], self.fields[FIELDS.index(name)]
            for t, z in NORMAL_ENTRIES:
                factor = factorial(t)*factorial(z)
                old = float(base.data[INDEX[t, z, 0], 0, 0, 0].real*factor)
                new = float(field.data[INDEX[t, z, 0], 0, 0, 0].real*factor)
                if new != old:
                    rows.append({"field": name, "normal_order": t, "spatial_order": z,
                                 "baseline_derivative": old, "derivative": new, "delta": new-old})
        return rows

    def normal_geometry(self):
        """Local axial/sphere expansions and spatial derivatives, not solutions."""
        self.validate()
        _, _, a, r = self.fields
        scalar = lambda jet: float(jet.value[0, 0, 0].real)
        return {"k_parallel": scalar(a.derivative(t=1)/a),
                "k_perp": scalar(r.derivative(t=1)/r),
                "dz_k_parallel": scalar((a.derivative(t=1)/a).derivative(z=1)),
                "dz_k_perp": scalar((r.derivative(t=1)/r).derivative(z=1)),
                "intrinsic_a": scalar(a), "intrinsic_r": scalar(r),
                "canonical_density_factor": float(scalar(r)*np.sqrt(scalar(a))),
                "normal": "n=partial_T on the unchanged rho=1 surface"}


def incoming_cauchy_jets(changes=()):
    """Construct twenty free normal derivative slots; default is the baseline.

    No preferred nonzero values are supplied. Setting N=1,beta=0 fixes the
    background gauge only: all four lapse/shift/spatial metric equations
    must still be varied and checked independently by a constraint owner.
    """
    changes = tuple(changes)
    if any(not isinstance(change, IncomingNormalJetChange) for change in changes):
        raise ValueError("explicit IncomingNormalJetChange entries required")
    keys = [(c.field, c.normal_order, c.spatial_order) for c in changes]
    if len(keys) != len(set(keys)):
        raise ValueError("each normal Taylor derivative may be supplied only once")
    baseline = reference_incoming_metric_jets()
    fields = tuple(SymbolJet(f.data, f.physical, f.momentum) for f in baseline)
    for change in changes:
        fields[FIELDS.index(change.field)].data[INDEX[change.normal_order, change.spatial_order, 0]] += (
            change.delta/(factorial(change.normal_order)*factorial(change.spatial_order))*I)
    domain = IncomingCauchyJets(baseline, fields, changes)
    domain.validate()
    return domain


def _raw_vertex_jets(fields, momenta, masses, angular):
    """The existing raw KS Hamiltonian vertices, including their Weyl jets."""
    N, _, a, r = fields
    k = SymbolJet.variable(momenta, 2)
    m, ell = SymbolJet.constant(masses), SymbolJet.constant(angular)
    return (-m.right(SIGMA[0])+(ell/r).right(SIGMA[1])+(k/a).right(SIGMA[2]),
            -k, -(N*k/(a*a)).right(SIGMA[2]), -(N*ell/(r*r)).right(SIGMA[1]))


def _symmetric_reference_vertices(projector, vertices):
    traces = []
    for vertex in vertices:
        left = values(product(projector['jets'], [vertex]))
        right = values(product([vertex], projector['jets']))
        traces.append(.5*np.trace(left+right, axis1=-2, axis2=-1))
    return np.stack(traces, axis=-1)


def incoming_reference_uv_difference(domain, channels, momentum_magnitudes):
    """Bounded high-|k| screen of the fixed-C0 conditional initial source.

    ``channels`` is a sequence of dictionaries with name, mass and angular
    labels from the retained inventory; no angular/copy multiplicity is added.
    Positive increasing magnitudes are evaluated at BOTH momentum signs.
    Since H and its pointwise raw vertices on the surface are unchanged,
    these are the state-independent changes of the positive reference-action
    vertex. The canonical (C0-P) moment changes by their negative.

    The output preserves individual formal orders, signs, and four vertices.
    It does not include the band's time/boundary remainder or local-action
    change. A decreasing sampled tail is a plausibility screen, not a proof
    of Hadamard regularity, full stress convergence or parent preparation.
    """
    if not isinstance(domain, IncomingCauchyJets):
        raise ValueError("the same-surface conditional incoming metric-jet domain is required")
    domain.validate()
    kabs = np.asarray(momentum_magnitudes)
    if (np.iscomplexobj(kabs) or kabs.ndim != 1 or len(kabs) < 3
            or not np.isfinite(kabs).all() or np.any(kabs <= 0) or np.any(np.diff(kabs) <= 0)):
        raise ValueError("at least three positive increasing numerical momentum magnitudes required")
    kabs = kabs.astype(float)
    channels = tuple(channels)
    if not channels or len({c['name'] for c in channels}) != len(channels):
        raise ValueError("distinct named inherited channel controls required")
    masses = np.array([c['mass'] for c in channels], float)
    angular = np.array([c['angular'] for c in channels], float)
    if not np.isfinite(masses+angular).all() or np.any(masses < 0):
        raise ValueError("finite inherited masses and signed angular labels required")
    nk, nc = len(kabs), len(channels)
    k = np.tile(np.r_[-kabs, kabs], nc)
    m, ell = np.repeat(masses, 2*nk), np.repeat(angular, 2*nk)
    base = exact_massless_reference(reference_projector(*domain.baseline, k, m, ell), m, ell, k)
    new = exact_massless_reference(reference_projector(*domain.fields, k, m, ell), m, ell, k)
    old_vertices = _raw_vertex_jets(domain.baseline, k, m, ell)
    new_vertices = _raw_vertex_jets(domain.fields, k, m, ell)
    before = _symmetric_reference_vertices(base, old_vertices)
    after = _symmetric_reference_vertices(new, new_vertices)
    delta = (after-before).reshape(5, nc, 2, nk, 4)
    projector_difference = (new['orders']-base['orders']).reshape(5, nc, 2, nk, 2, 2)
    real_delta = delta.real.sum(axis=0)
    envelope = np.max(abs(real_delta), axis=1)
    paired = real_delta.sum(axis=1)
    # Zero/roundoff channels are marked unresolved rather than assigning a
    # fake decay exponent. The caller decides the relevant accuracy scale.
    power = np.zeros((nc, nk-1, 4))
    resolved = (envelope[:, :-1] > 1e-24) & (envelope[:, 1:] > 1e-24)
    ratio = np.divide(envelope[:, :-1], envelope[:, 1:], out=np.ones_like(power), where=resolved)
    np.divide(np.log(ratio), np.log(kabs[1:]/kabs[:-1])[None, :, None], out=power, where=resolved)
    formal = max(float(np.max(entry['scaled'])) for p in (base, new)
                 for order in p['residuals'] for entry in order.values())
    return {
        'channel_names': [c['name'] for c in channels], 'momenta': k.reshape(nc, 2, nk),
        'momentum_magnitudes': kabs, 'vertex_order': FIELDS,
        'projector_difference_orders': projector_difference,
        'reference_vertex_difference_orders': delta,
        'reference_vertex_difference': real_delta,
        'canonical_subtracted_moment_difference': -real_delta,
        'signed_pair_reference_difference': paired,
        'unsigned_vertex_envelope': envelope,
        'sampled_decay_power': power, 'decay_power_resolved': resolved,
        'initial_hamiltonian_difference': float(np.max(abs(new['hamiltonian'].value-base['hamiltonian'].value))),
        'initial_raw_vertex_difference': max(float(np.max(abs(a.value-b.value))) for a, b in zip(new_vertices, old_vertices)),
        'initial_projector_zero_order_difference': float(np.max(abs(projector_difference[0]))),
        'maximum_formal_reference_residual': formal,
        'maximum_symmetric_vertex_imaginary_part': float(np.max(abs(delta.imag))),
        'changed_normal_entries': domain.changed_normal_entries(),
        'scope': {
            'same_intrinsic_surface_normal_frame': True, 'fixed_abstract_canonical_C0': True,
            'physical_history_or_duration_selected': False, 'constraints_solved': False,
            'stress_integrated': False, 'band_boundary_remainder_evaluated': False,
            'full_UV_or_Hadamard_certificate': False, 'global_horizon_preparation_preserved': None,
            'LLL_projector': 'exact constant chiral reference; no duplicated geometric allocation',
            'LLL_local_allocation_owner': 'nsc_lll_geometric_history.py, counted separately once',
            'LLL_local_allocation_recomputed': False,
            'interpretation': 'sampled reference-vertex UV difference for a conditional Cauchy-data experiment',
        },
    }
