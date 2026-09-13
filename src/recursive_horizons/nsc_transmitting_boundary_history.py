"""Bind the existing transmitting-action declaration to the locked NSC records.

This owner executes identifiability and endpoint assembly. It does not supply
an absent scattering action, select a witness as physical, or solve a metric.
The finite covariance probe only asks whether the unselected interface family
can change the state in the fixed, recorded basis.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np


FIELDS = ("N", "beta", "q_ADM", "r")
LOCAL_NAMES = ("lapse", "beta", "q_ADM", "r")


@dataclass(frozen=True)
class EndpointVariation:
    """Action covector on explicitly identified nodal endpoint coordinates.

Both signs are inherited from dS/dg, including the owner's time orientation.
No second outward-normal sign or stress/volume conversion is applied here.
This object is not a continuum boundary momentum without a boundary pullback.
"""

    basis_id: str
    components: tuple[tuple[float, float], ...]

    def __post_init__(self):
        values = np.asarray(self.components)
        if not self.basis_id or values.shape != (4, 2) or not np.isfinite(values).all():
            raise ValueError("finite four-field/two-endpoint action covector required")

    def to_dict(self):
        return {name: list(pair) for name, pair in zip(FIELDS, self.components)}

    @property
    def maximum_absolute(self):
        return float(np.max(np.abs(self.components)))


def match_endpoint_variations(local, boundary, *, tolerance=3e-11):
    """Add derivatives only on an identical pulled-back endpoint basis.

None is an unevaluated derivative, never an implicit zero or fitted negative
of the local vector. A numerical match alone does not certify provenance.
"""
    if not isfinite(tolerance) or tolerance <= 0:
        raise ValueError("positive finite tolerance required")
    if boundary is None:
        return {
            "status": "OPEN", "evaluated": False,
            "boundary_components": None, "mismatch_components": None,
            "maximum_absolute": None, "tolerance": tolerance,
            "reason": "Gamma_rest has no evaluated metric derivative and no transmitting-to-KS endpoint pullback",
        }
    if local.basis_id != boundary.basis_id:
        raise ValueError("boundary derivative must use the same endpoint basis and history")
    total = np.asarray(local.components) + np.asarray(boundary.components)
    norm = float(np.max(np.abs(total)))
    return {
        "status": "PASS" if norm <= tolerance else "FAIL", "evaluated": True,
        "boundary_components": boundary.to_dict(),
        "mismatch_components": {field: row.tolist() for field, row in zip(FIELDS, total)},
        "maximum_absolute": norm, "tolerance": tolerance,
        "reason": "numerical covector sum only; action provenance is a separate gate",
    }


def covariance_freedom_probe(arrays):
    """One infinitesimal first/last-frequency mixing direction per channel.

At V=I, dV=iX, X=X^dagger, ||X||_F=1. The isometry tangent is
dV^dagger+dV=0, while the fixed-basis state tangent is i[X,C].
Only the four affected spin entries are formed. No temporal propagation,
stress integration, optimization or new occupation is introduced.
"""
    offsets = np.asarray(arrays["sample_offsets"], dtype=int)
    covariance = np.asarray(arrays["covariance_seed"], dtype=complex)
    frequency = np.asarray(arrays["frequency"], dtype=float)
    if offsets[0] != 0 or offsets[-1] != len(covariance) or np.any(np.diff(offsets) < 2):
        raise ValueError("at least two retained frequencies per authenticated channel required")
    x = np.zeros((4, 4), dtype=complex)
    x[:2, 2:] = np.eye(2)/2
    x[2:, :2] = np.eye(2)/2
    dv = 1j*x
    rows = []
    for channel, (start, stop) in enumerate(zip(offsets[:-1], offsets[1:])):
        c = np.zeros((4, 4), dtype=complex)
        c[:2, :2] = covariance[start]
        c[2:, 2:] = covariance[stop-1]
        dc = dv @ c + c @ dv.conj().T
        rows.append({
            "channel": channel, "frequency_indices": [int(start), int(stop-1)],
            "frequencies": [float(frequency[start]), float(frequency[stop-1])],
            "isometry_tangent_maximum_absolute": float(np.max(np.abs(dv.conj().T+dv))),
            "fixed_basis_covariance_tangent_frobenius": float(np.linalg.norm(dc)),
            "covariance_tangent_trace_absolute": float(abs(np.trace(dc))),
            "covariance_tangent_hermiticity_maximum_absolute": float(np.max(np.abs(dc-dc.conj().T))),
        })
    tol = 3e-11
    return {
        "definition": "dV=iX; first/last-frequency spin-preserving X per channel; ||X||_F=1; dC=i[X,C_seed]",
        "rows": rows,
        "maximum_isometry_tangent_residual": max(r["isometry_tangent_maximum_absolute"] for r in rows),
        "maximum_state_tangent": max(r["fixed_basis_covariance_tangent_frobenius"] for r in rows),
        "channels_with_state_tangent_above_tolerance": sum(r["fixed_basis_covariance_tangent_frobenius"] > tol for r in rows),
        "tolerance": tol,
        "interpretation": "some allowed maps change C in the fixed source basis; this is not a claim about a selected physical kernel or stress",
        "quotient_dimension_after_gauge_and_state_stabilizers": None,
        "physical_map_selected": False,
    }


class ModeResolvedTransmittingBoundaryHistoryAction:
    """Executable binding audit of the existing Gamma_rest declaration."""

    def selection(self, interface, compact_boundary, arrays):
        matching = compact_boundary["metric_boundary_matching"]
        if compact_boundary["scope"]["transmitting_throat_boundary_action_derived"]:
            raise ValueError("a new transmitting-action dependency requires a fresh binding audit")
        if matching["higher_order_or_transmission_metric_variation_solved"]:
            raise ValueError("new endpoint variation must be bound explicitly")
        if interface["physical_selector"]["selected"]:
            raise ValueError("selected kernel requires its authenticated variational residual")
        probe = covariance_freedom_probe(arrays)
        ranks = interface["rank_nonuniqueness_certificate"]
        dimension = sum(int(2*n)**2 for n in np.diff(arrays["sample_offsets"]))
        if dimension != ranks["admissible_family_real_dimension_and_linearized_nullity"]:
            raise ValueError("channel dimensions differ from the locked interface family")
        return {
            "status": "OPEN", "physical_Vc": None,
            "selection_residual": {"evaluated": False, "value": None, "tolerance": 3e-11},
            "admissible_family_real_dimension": dimension,
            "additional_selector_rank": None,
            "same_action_equation": matching["common_interface_equation"],
            "missing": [
                "mode-resolved transmitting link/embedding functional in the locked channel bases",
                "its variation with respect to channel scattering data and the history metric",
            ],
            "covariance_freedom_probe": probe,
            "reason": "the existing action inventory declares Gamma_rest but supplies no functional selecting V_c; norm conservation alone leaves the recorded family",
        }

    def endpoint(self, local, *, local_record_sha256):
        gradients = local["diagnostic_evaluation"]["action_gradients"]["weyl_bulk"]
        components = tuple((float(gradients[name][0]), float(gradients[name][-1])) for name in LOCAL_NAMES)
        basis_id = "KS-local-nodal-endpoints:"+local_record_sha256
        covector = EndpointVariation(basis_id, components)
        stored = local["residuals"]["weyl_endpoint_completion"]
        extraction = max(abs(a-b) for pair, name in zip(components, LOCAL_NAMES)
                         for a, b in zip(pair, stored["components"][name]))
        if extraction != 0 or covector.maximum_absolute != stored["maximum_absolute"]:
            raise ValueError("endpoint object does not reproduce the authenticated nodal gradient")
        count = len(gradients["r"])
        # Keep the entire already computed gradient; the end-node projection is
        # not a substitute for higher-derivative boundary jets or their pullback.
        return {
            "status": "OPEN",
            "basis": {
                "id": basis_id, "fields": list(FIELDS), "local_fields": list(LOCAL_NAMES),
                "node_indices": [0, count-1], "history_nodes": count,
                "orientation": "increasing stored time; signs already included in dS/dg",
                "units": "dimensionless recorded action per unit homogeneous axial coordinate; derivatives with respect to raw metric nodes",
                "history": local["owner"]["diagnostic_history"],
                "physical_history": False, "derivative_kind": "discrete action gradient, not a continuum boundary momentum",
            },
            "local_components": covector.to_dict(),
            "local_maximum_absolute": covector.maximum_absolute,
            "local_full_nodal_gradient": {field: gradients[name] for field, name in zip(FIELDS, LOCAL_NAMES)},
            "binding_residual": {"maximum_absolute": extraction, "tolerance": 0.0, "status": "PASS"},
            "match": match_endpoint_variations(covector, None, tolerance=stored["tolerance"]),
            "endpoint_pullback": {
                "transmitting_embedding_to_KS_nodes": None,
                "allowed_endpoint_and_higher_jet_variations": None,
                "reason": "fixed endpoint values do not specify the boundary jets of a higher-derivative action; the transmitting boundary is not the diagnostic temporal end node",
            },
            "independent_first_jet_derivative": {
                "available": False,
                "components": {field: [None, None] for field in FIELDS},
                "reason": "the existing nodal schema varies metric values, not independently specified normal derivatives; a boundary-jet decomposition and its admissible variations are not supplied",
            },
            "local_value_erased_or_cancelled": False,
        }


def compose_full_extended_gate(base_gate, selection, endpoint):
    """Assemble evaluated and missing residual channels without solving them."""
    if not base_gate["components_expose_residuals"]:
        raise ValueError("the three imported components must expose their residuals")
    blockers = []
    if selection["status"] != "PASS":
        blockers.append("delta_Gamma_rest/dV_c: physical transmitting kernel selector")
    if endpoint["match"]["status"] != "PASS":
        blockers.append("pullback(delta_Gamma_rest/dh) to the local Weyl endpoint variation basis")
    if not blockers:
        # A selected interface is necessary but not sufficient for stationarity.
        blockers.append("assembled history stationarity residuals for g_star,C_star")
    return {
        "decision": "OPEN", "blocking_residuals": blockers,
        "component_construction_pass": True,
        "residual_surface": {
            "fields": list(FIELDS),
            "bulk_reference_and_local": "authenticated component records on supplied histories",
            "kernel_selection": selection["selection_residual"],
            "endpoint_matching": endpoint["match"],
            "history_variation": {field: None for field in FIELDS},
        },
        "old_homogeneous_nonexistence": base_gate["killed_homogeneous_regression_pass"],
        "old_homogeneous_result_reused_not_rerun": True,
        "extended_existence_claimed": False, "extended_nonexistence_claimed": False,
        "optimizer_started": False, "finite_stress": None, "nulls": None,
        "updated_constraints": None, "metric_timestep_started": False,
        "coupled_evolution_reopened": False,
    }
