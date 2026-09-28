"""Nine-component Gate Budget v5: fill only from authenticated PASS enclosures.

Allocations are frozen. A component bound is written only from an
authenticated applicable PASS enclosure. Predecessor v4 bytes stay
immutable. The physical local gate remains OPEN until every component is
finite, the directed total is at most 2e-11 componentwise, and at least
1e-11 residual reserve remains.
"""
from __future__ import annotations

from fractions import Fraction
import math

import numpy as np

from .nsc_ks_assembly_arithmetic import (
    DECLARED_POSITIVE_FAMILIES,
    DECLARED_SIGNED_FAMILIES,
    signed_family_coverage,
)
from .nsc_ks_continuous_constraint_bound import VERIFICATION_NODE_COUNT
from .nsc_ks_continuous_constraint_enclosure import (
    PROFILE_IDENTITY,
    STATE_LAW,
    INTERVAL_LABEL,
    _sha256_text,
    authenticate_evidence_bytes,
    binary_rational,
)
from .nsc_local_gate_certificate_v2 import ERROR_COMPONENTS


SCHEMA = "NSC-KS-GATE-ERROR-BUDGET-v5"
PREDECESSOR_SCHEMA = "NSC-KS-GATE-ERROR-BUDGET-v4"
CONSTRAINTS = ("N", "beta")
ALLOCATION = {
    "field_space_time": 5e-12,
    "changed_history_UV_tail": 5e-12,
    "baseline_low_subgap": 3e-12,
    "upstream": 2e-12,
    "energy_interpolation": 1e-12,
    "covered_regions": 1e-12,
    "phase_value": 1e-12,
    "between_node": 1e-12,
    "arithmetic": 1e-12,
}
ALLOCATION_TOTAL = 2e-11
RESIDUAL_RESERVE = 1e-11
PHYSICAL_TOLERANCE = 3e-11
PASS_TOKEN = "PASS"


def _finite_pair(value, name):
    array = np.asarray(value, dtype=float)
    if array.shape != (2,) or not np.isfinite(array).all() or np.any(array < 0):
        raise ValueError(name + " must be a finite nonnegative N,beta pair")
    return (float(array[0]), float(array[1]))


def _status_token(status):
    if not isinstance(status, str) or not status:
        return ""
    return status.split(":", 1)[0]


def empty_components():
    return {
        name: {
            "allocation_per_component": allocation,
            "bound": None,
            "status": "OPEN: no authenticated PASS enclosure for the current history",
        }
        for name, allocation in ALLOCATION.items()
    }


def _require_node_count(name, enclosure, *, node_count=VERIFICATION_NODE_COUNT):
    declared = enclosure.get("node_count")
    if declared is None:
        raise ValueError(name + " node_count required")
    try:
        declared = int(declared)
    except (TypeError, ValueError) as error:
        raise ValueError(name + " node_count required") from error
    if declared == 129:
        raise ValueError(name + " 129-node phase cannot fill the 257-node budget")
    if declared != int(node_count) or declared != VERIFICATION_NODE_COUNT:
        raise ValueError(name + " requires 257 verification nodes")
    return declared


def fill_component(name, enclosure, *,
                   profile_identity=PROFILE_IDENTITY,
                   state_law=STATE_LAW,
                   root=None,
                   node_count=VERIFICATION_NODE_COUNT):
    """Admit one component bound only from an authenticated applicable PASS."""
    if name not in ALLOCATION:
        raise ValueError("unknown budget component: " + name)
    if enclosure is None:
        raise ValueError(name + " enclosure is missing")
    if not isinstance(enclosure, dict):
        raise TypeError(name + " enclosure must be a mapping")
    if enclosure.get("certificate_use") is False:
        raise ValueError(name + " diagnostic enclosure cannot fill the budget")
    if enclosure.get("certificate_use") is not True:
        raise ValueError(name + " certificate_use is required")
    if enclosure.get("numerical_indicators_used_as_bounds") is True:
        raise ValueError(name + " samples/indicators cannot fill the budget")
    if enclosure.get("samples_used_as_bound") is True:
        raise ValueError(name + " samples cannot fill the budget")
    if _status_token(enclosure.get("status")) != PASS_TOKEN:
        raise ValueError(name + " fills only from an authenticated PASS enclosure")
    bound = enclosure.get("bound")
    pair = _finite_pair(bound, name)
    declared_profile = enclosure.get("profile_identity")
    if declared_profile is None:
        raise ValueError(name + " profile identity required")
    if declared_profile != profile_identity:
        raise ValueError(name + " profile identity mismatch")
    declared_law = enclosure.get("state_law")
    if declared_law is None:
        raise ValueError(name + " state law required")
    if declared_law != state_law:
        raise ValueError(name + " state law mismatch")
    declared_nodes = _require_node_count(name, enclosure, node_count=node_count)
    evidence = enclosure.get("evidence")
    if evidence is None:
        raise ValueError(name + " evidence path and sha256 required")
    relative, digest, _counts = authenticate_evidence_bytes(
        evidence, root=root, name=name + " evidence")
    excess = tuple(max(pair[index] - ALLOCATION[name], 0.0) for index in range(2))
    row = {
        "allocation_per_component": ALLOCATION[name],
        "bound": [pair[0], pair[1]],
        "status": enclosure.get("status"),
        "allocation_excess": list(excess),
        "scope": enclosure.get("scope"),
        "profile_identity": declared_profile,
        "state_law": declared_law,
        "node_count": declared_nodes,
        "evidence": {
            **evidence,
            "path": relative,
            "sha256": digest,
        },
    }
    if enclosure.get("exact_upper_evidence") is not None:
        row["exact_upper_evidence"] = enclosure["exact_upper_evidence"]
    return row


def directed_component_sum(components):
    if set(components) != set(ERROR_COMPONENTS):
        raise ValueError("exactly the nine named components required")
    missing = [name for name in ERROR_COMPONENTS if components[name].get("bound") is None]
    if missing:
        return None, missing
    total = [Fraction(0), Fraction(0)]
    for name in ERROR_COMPONENTS:
        pair = _finite_pair(components[name]["bound"], name)
        for index in range(2):
            exact = binary_rational(pair[index])
            total[index] += exact
            widened = float(total[index])
            while binary_rational(widened) < total[index]:
                nxt = math.nextafter(widened, math.inf)
                if nxt == widened:
                    raise ValueError("directed component sum overflow")
                widened = nxt
            total[index] = binary_rational(widened)
    return (float(total[0]), float(total[1])), []


def admit_budget(record, *, family_inventory=None):
    """Strict v5 admission used by the successor search. OPEN records fail."""
    if not isinstance(record, dict):
        raise ValueError("Gate-2 budget must be an object")
    if record.get("schema") != SCHEMA:
        raise ValueError("Gate-2 requires NSC-KS-GATE-ERROR-BUDGET-v5")
    components = record.get("components")
    if not isinstance(components, dict) or set(components) != set(ERROR_COMPONENTS):
        raise ValueError("Gate-2 requires exactly the nine named components")
    if record.get("allocation_changed") is True:
        raise ValueError("allocations may not be moved")
    if record.get("allocation_total") != ALLOCATION_TOTAL:
        raise ValueError("allocation total is frozen at 2e-11")
    if float(record.get("residual_reserve", float("nan"))) < RESIDUAL_RESERVE:
        raise ValueError("residual reserve must be at least 1e-11")
    if record.get("physical_tolerance") != PHYSICAL_TOLERANCE:
        raise ValueError("physical tolerance is frozen at 3e-11")
    if record.get("constraint_order") != list(CONSTRAINTS):
        raise ValueError("constraint order is N then beta")
    if record.get("state_law") != STATE_LAW:
        raise ValueError("switched state law binding changed")
    if record.get("numerical_indicators_used_as_bounds") is not False:
        raise ValueError("indicators cannot be used as bounds")
    if record.get("profile_identity") != PROFILE_IDENTITY:
        raise ValueError("budget profile identity mismatch")
    for name, allocation in ALLOCATION.items():
        row = components[name]
        if float(row.get("allocation_per_component", float("nan"))) != allocation:
            raise ValueError("allocation moved for " + name)
    if record.get("physical_EXISTENCE_certificate"):
        raise ValueError("a budget is not an EXISTENCE certificate")
    if record.get("physical_NONEXISTENCE_certificate"):
        raise ValueError("a budget is not a NON_EXISTENCE certificate")
    if record.get("missing_components") not in ([], ()):
        raise ValueError("Gate-2 budget still has missing components")
    total, missing = directed_component_sum(components)
    if missing:
        raise ValueError("missing components: " + ", ".join(missing))
    declared = np.asarray(record.get("full_error_sum"), dtype=float)
    if (declared.shape != (2,) or not np.isfinite(declared).all()
            or np.any(declared < np.asarray(total))):
        raise ValueError("full_error_sum does not enclose the nine components")
    if np.any(declared > ALLOCATION_TOTAL):
        raise ValueError("componentwise total exceeds 2e-11")
    if PHYSICAL_TOLERANCE - float(np.max(declared)) < RESIDUAL_RESERVE:
        raise ValueError("residual reserve does not fit 3e-11")
    coverage = signed_family_coverage(
        family_inventory if family_inventory is not None
        else record.get("source_inventory"))
    if (coverage["authenticated"] is not True
            or coverage["covered_positive_families"] != DECLARED_POSITIVE_FAMILIES
            or coverage["covered_signed_families"] != DECLARED_SIGNED_FAMILIES
            or coverage["all_declared_source_paths_closed"] is not True):
        raise ValueError("Gate-2 lacks authenticated 60/120 source coverage")
    if record.get("covered_positive_families") != coverage["covered_positive_families"]:
        raise ValueError("covered_positive_families does not match authenticated coverage")
    if record.get("covered_signed_families") != coverage["covered_signed_families"]:
        raise ValueError("covered_signed_families does not match authenticated coverage")
    if _status_token(record.get("status")) not in ("PASS", "CLOSED", "ENCLOSED"):
        raise ValueError("Gate-2 budget status is not closed")
    return {
        "authorized": True,
        "component_sum": list(map(float, declared)),
        "allocation_total": ALLOCATION_TOTAL,
        "residual_reserve": RESIDUAL_RESERVE,
        "physical_tolerance": PHYSICAL_TOLERANCE,
    }


def compose_open_budget(*, predecessor=None, pass_enclosures=None,
                        profile_identity=PROFILE_IDENTITY,
                        source_identity=None,
                        predecessor_sha256=None,
                        root=None,
                        family_inventory=None):
    """Build the current v5 record. Unfilled components stay null."""
    if source_identity is not None:
        _sha256_text(source_identity, "source_identity")
    coverage = signed_family_coverage(family_inventory)
    components = empty_components()
    predecessor_filled = []
    if predecessor is not None:
        if predecessor.get("schema") != PREDECESSOR_SCHEMA:
            raise ValueError("predecessor must be NSC-KS-GATE-ERROR-BUDGET-v4")
        if predecessor.get("allocation") not in (None, ALLOCATION):
            if predecessor.get("allocation_changed") is True:
                raise ValueError("predecessor allocations moved")
        for name, allocation in ALLOCATION.items():
            row = predecessor["components"][name]
            if float(row["allocation_per_component"]) != allocation:
                raise ValueError("predecessor allocation moved for " + name)
        # Predecessor enclosed rows are recorded as evidence only. They are
        # not PASS enclosures of this successor, so they do not fill v5.
        predecessor_filled = [
            name for name in ERROR_COMPONENTS
            if predecessor["components"][name].get("bound") is not None
        ]
    filled = []
    enclosures = pass_enclosures or {}
    extra = sorted(set(enclosures) - set(ERROR_COMPONENTS))
    if extra:
        raise ValueError("unknown PASS enclosure: " + extra[0])
    for name, enclosure in enclosures.items():
        components[name] = fill_component(
            name, enclosure, profile_identity=profile_identity, state_law=STATE_LAW,
            root=root)
        filled.append(name)
    missing = [name for name in ERROR_COMPONENTS if components[name]["bound"] is None]
    total, _missing = directed_component_sum(components)
    partial = None
    if missing:
        known = []
        for name in ERROR_COMPONENTS:
            bound = components[name].get("bound")
            if bound is not None:
                known.append(_finite_pair(bound, name))
        if known:
            partial = [0.0, 0.0]
            for pair in known:
                for index in range(2):
                    widened = math.nextafter(partial[index] + pair[index], math.inf)
                    partial[index] = widened
    status = (
        "OPEN: nine-component v5 composer; components fill only from "
        "authenticated PASS enclosures; "
        + str(len(missing)) + " remain null"
    )
    return {
        "schema": SCHEMA,
        "status": status,
        "constraint_order": list(CONSTRAINTS),
        "state_law": STATE_LAW,
        "profile_identity": profile_identity,
        "source_identity": source_identity,
        "interval": INTERVAL_LABEL,
        "components": components,
        "missing_components": missing,
        "filled_from_pass": filled,
        "predecessor_schema": PREDECESSOR_SCHEMA,
        "predecessor_enclosed_components": predecessor_filled,
        "predecessor_sha256": predecessor_sha256,
        "predecessor_bounds_are_not_v5_pass_fills": True,
        "partial_known_error_sum": partial,
        "full_error_sum": None if missing else list(total),
        "allocation_total": ALLOCATION_TOTAL,
        "residual_reserve": RESIDUAL_RESERVE,
        "physical_tolerance": PHYSICAL_TOLERANCE,
        "allocation_changed": False,
        "declared_positive_families": coverage["declared_positive_families"],
        "declared_signed_families": coverage["declared_signed_families"],
        "covered_positive_families": coverage["covered_positive_families"],
        "covered_signed_families": coverage["covered_signed_families"],
        "source_inventory_authenticated": coverage["authenticated"] is True,
        "all_declared_source_paths_closed": coverage[
            "all_declared_source_paths_closed"] is True,
        "numerical_indicators_used_as_bounds": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "named_gap": "pass_enclosures_missing_for_null_v5_components",
        "arithmetic_scope": (
            "remaining geometry, field and assembly arithmetic; "
            "phase arithmetic is already in phase_value"),
        "phase_continuous_remainder_included": False,
        "v4_bytes_preserved": True,
    }


def validate_open_budget(record):
    if record.get("schema") != SCHEMA:
        raise ValueError("NSC-KS-GATE-ERROR-BUDGET-v5 required")
    if record.get("allocation_changed") is True:
        raise ValueError("allocations may not be moved")
    if record.get("v4_bytes_preserved") is not True:
        raise ValueError("existing v4 bytes must stay preserved")
    if set(record.get("components", {})) != set(ERROR_COMPONENTS):
        raise ValueError("exactly the nine named components required")
    for name, allocation in ALLOCATION.items():
        row = record["components"][name]
        if float(row["allocation_per_component"]) != allocation:
            raise ValueError("allocation moved for " + name)
    if record.get("physical_EXISTENCE_certificate"):
        raise ValueError("OPEN budget may not be reported as EXISTENCE")
    if record.get("physical_NONEXISTENCE_certificate"):
        raise ValueError("OPEN budget may not be reported as NON_EXISTENCE")
    if not str(record.get("status", "")).startswith("OPEN"):
        if record.get("missing_components"):
            raise ValueError("a budget with null components remains OPEN")
    missing = [name for name in ERROR_COMPONENTS
               if record["components"][name].get("bound") is None]
    if missing != list(record.get("missing_components", [])):
        raise ValueError("missing_components does not match null bounds")
    if missing and record.get("full_error_sum") is not None:
        raise ValueError("full_error_sum stays null while components are null")
    if missing and record.get("all_declared_source_paths_closed") is True:
        raise ValueError("source paths cannot be closed while components are null")
    if record.get("declared_positive_families") != DECLARED_POSITIVE_FAMILIES:
        raise ValueError("declared positive family count is 60")
    if record.get("declared_signed_families") != DECLARED_SIGNED_FAMILIES:
        raise ValueError("declared signed family count is 120")
    if record.get("source_inventory_authenticated") is not True:
        if record.get("covered_positive_families") is not None:
            raise ValueError(
                "unauthenticated coverage cannot report covered_positive_families")
        if record.get("covered_signed_families") is not None:
            raise ValueError(
                "unauthenticated coverage cannot report covered_signed_families")
    if record.get("predecessor_bounds_are_not_v5_pass_fills") is not True:
        raise ValueError("v4 enclosed rows are not silent v5 PASS fills")
    return record
