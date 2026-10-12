"""Strict 257-node / 256-cell continuous enclosure of both incoming constraints.

The certification grid on I is the owned Chebyshev-Lobatto set
``LocalIncomingFamily.collocation_nodes(257)``. Production input is an
outward nodal |E_N|, |E_beta| enclosure at every node and a directed
derivative upper bound on every cell and component. Samples are
diagnostics. Missing cells, a reordered grid, an inconsistent derivative
bound, a hash/profile/state-law mismatch, NaN/negative values or partial
coverage cannot fill the between-node budget slot.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
from types import MappingProxyType

import numpy as np

from .nsc_ks_continuous_constraint_bound import (
    CELL_COUNT,
    VERIFICATION_NODE_COUNT,
    verification_componentwise_continuous_max,
)
from .nsc_ks_evaluation_binding import (
    array_digest,
    file_digest,
    hex_float,
    repository_relative,
)
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_source_envelope import physical_incoming_interval
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_local_gate_certificate_v2 import ERROR_COMPONENTS

REPO_ROOT = Path(__file__).resolve().parents[2]


SCHEMA = "NSC-KS-CONTINUOUS-CONSTRAINT-ENCLOSURE-v1"
PROFILE_IDENTITY = (
    "0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0")
STATE_LAW = "C_Sigma[g]=U_g C_up U_g^dagger"
INTERVAL_LABEL = "I=S(1)+[.12,.18]"
CONSTRAINTS = ("N", "beta")
GRID_OWNER = "LocalIncomingFamily.collocation_nodes"
CANDIDATE_OWNER = "results/development/nsc-ks-gate-history-lm-broyden.json"
EVALUATOR_OWNERS = (
    "src/recursive_horizons/nsc_ks_history_evaluator.py",
    "src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py",
    "src/recursive_horizons/nsc_ks_value_evaluator.py",
    "src/recursive_horizons/nsc_ks_batched_constraints.py",
)
HEX64 = tuple("0123456789abcdef")
MISSING_PRODUCTION_INPUTS = (
    "outward nodal |E_N|, |E_beta| enclosures at all 257 verification nodes",
    "directed full-residual derivative uppers on all 256 cells for N and beta",
    "interval/directed abs majorants of remaining contraction arrays "
    "(delta, C_src, reference columns, vertices) at the assembly nodes",
    "family-summand interval enclosures for all 120 signed families",
    "authenticated original weights, signs, degeneracies and coherences "
    "for all 120 signed families from the source inventory",
    "interval enclosures of baseline, geometry and matter terms; edge phase "
    "arithmetic is already enclosed by phase_value and is not counted again",
    "authenticated PASS records for field_space_time, changed_history_UV_tail, "
    "baseline_low_subgap, upstream, between_node and arithmetic",
)


def _sha256_text(value, name):
    if (not isinstance(value, str) or len(value) != 64
            or any(char not in HEX64 for char in value)):
        raise ValueError("lowercase SHA-256 " + name + " required")
    return value


def _nonneg_fraction(value, name):
    if value is None or isinstance(value, bool):
        raise ValueError("explicit nonnegative " + name + " required")
    if isinstance(value, dict) and set(value) == {"numerator", "denominator"}:
        value = Fraction(int(value["numerator"]), int(value["denominator"]))
    elif isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("finite nonnegative " + name + " required")
        value = Fraction(str(value))
    else:
        try:
            value = Fraction(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError("finite nonnegative " + name + " required") from error
    if value < 0:
        raise ValueError(name + " must be nonnegative")
    return value


def serialize_fraction(value):
    if value is None:
        return None
    number = value if isinstance(value, Fraction) else _nonneg_fraction(value, "fraction")
    return {"numerator": str(number.numerator), "denominator": str(number.denominator)}


def binary_rational(value):
    """Exact rational of a finite binary64 value."""
    rounded = float(value)
    if not math.isfinite(rounded):
        raise ValueError("finite binary64 required")
    return Fraction(*rounded.as_integer_ratio())


def outward_float(value):
    """Least finite binary64 that is >= the exact nonnegative rational."""
    number = value if isinstance(value, Fraction) else _nonneg_fraction(value, "outward")
    rounded = float(number)
    if not math.isfinite(rounded):
        raise ValueError("finite outward float required")
    while binary_rational(rounded) < number:
        nxt = math.nextafter(rounded, math.inf)
        if nxt == rounded:
            raise ValueError("outward float overflow")
        rounded = nxt
    return rounded


def authenticate_file_digest(path, digest, *, root=None, name="evidence"):
    """Require a repository-relative path whose bytes match sha256."""
    root = REPO_ROOT if root is None else Path(root)
    relative = repository_relative(path, root)
    _sha256_text(digest, name)
    target = root / relative
    if not target.is_file():
        raise ValueError(name + " file missing: " + relative)
    actual = file_digest(target)
    if actual != digest:
        raise ValueError(name + " sha256 does not match file bytes: " + relative)
    return relative, actual


def evidence_grid_node_counts(payload):
    """Collect declared grid sizes from an evidence mapping, if present."""
    if not isinstance(payload, dict):
        return ()
    counts = []
    for value in (payload.get("node_count"), payload.get("verification_node_count")):
        if value is not None:
            counts.append(int(value))
    binding = payload.get("binding")
    if isinstance(binding, dict) and binding.get("node_count") is not None:
        counts.append(int(binding["node_count"]))
    return tuple(counts)


def authenticate_evidence_bytes(evidence, *, root=None, name="evidence"):
    """Match evidence sha256 to file bytes and reject a 129-node grid."""
    if not isinstance(evidence, dict):
        raise ValueError(name + " path and sha256 required")
    path = evidence.get("path")
    digest = evidence.get("sha256")
    if not path or not digest:
        raise ValueError(name + " path and sha256 required")
    relative, actual = authenticate_file_digest(
        path, digest, root=root, name=name)
    target = (REPO_ROOT if root is None else Path(root)) / relative
    try:
        payload = json.loads(target.read_bytes().decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        return relative, actual, ()
    counts = evidence_grid_node_counts(payload)
    for count in counts:
        if count == 129:
            raise ValueError(name + " 129-node phase cannot fill the 257-node budget")
        if count != VERIFICATION_NODE_COUNT:
            raise ValueError(name + " grid is not the 257-node verification set")
    return relative, actual, counts


def serialize_lipschitz(result):
    return {
        "nodal_maximum_upper": serialize_fraction(result["nodal_maximum_upper"]),
        "continuous_maximum_upper": serialize_fraction(
            result["continuous_maximum_upper"]),
        "between_node_remainder_upper": serialize_fraction(
            result["between_node_remainder_upper"]),
        "cell_maximum_uppers": [
            serialize_fraction(item) for item in result["cell_maximum_uppers"]],
        "outward_binary64": {
            "nodal_maximum_upper": outward_float(result["nodal_maximum_upper"]),
            "continuous_maximum_upper": outward_float(
                result["continuous_maximum_upper"]),
            "between_node_remainder_upper": outward_float(
                result["between_node_remainder_upper"]),
        },
        "derivative_bound_required": True,
        "dense_sampling_used_as_bound": False,
    }


@dataclass(frozen=True, slots=True)
class EnclosureBinding:
    profile_identity: str
    state_law: str
    source_identity: str
    interval: tuple[float, float]
    dependency_hashes: tuple[tuple[str, str], ...]

    def __post_init__(self):
        _sha256_text(self.profile_identity, "profile_identity")
        _sha256_text(self.source_identity, "source_identity")
        if self.state_law != STATE_LAW:
            raise ValueError("switched state law binding changed")
        left, right = (float(self.interval[0]), float(self.interval[1]))
        if (len(self.interval) != 2 or not math.isfinite(left)
                or not math.isfinite(right) or right <= left):
            raise ValueError("connected positive-length interval I required")
        if not self.dependency_hashes:
            raise ValueError("dependency hashes required")
        names = [name for name, _ in self.dependency_hashes]
        if len(set(names)) != len(names):
            raise ValueError("unique dependency hashes are required")
        for name, digest in self.dependency_hashes:
            if not isinstance(name, str) or not name or ".." in name:
                raise ValueError("canonical dependency path required")
            _sha256_text(digest, name)

    def as_record(self):
        return {
            "profile_identity": self.profile_identity,
            "state_law": self.state_law,
            "source_identity": self.source_identity,
            "interval": [hex_float(self.interval[0]), hex_float(self.interval[1])],
            "interval_label": INTERVAL_LABEL,
            "dependency_hashes": {
                name: digest for name, digest in self.dependency_hashes},
        }


def require_positive_interval(interval):
    bounds = np.asarray(interval, float)
    if bounds.shape != (2,) or not np.isfinite(bounds).all() or bounds[1] <= bounds[0]:
        raise ValueError("connected positive-length interval I required")
    return (float(bounds[0]), float(bounds[1]))


def require_interval(interval=None, *, grid_interval=None, physical=True):
    """Bind the declared physical I. Grid coverage uses the family interval."""
    expected = physical_incoming_interval()
    if interval is None:
        if not physical:
            raise ValueError("explicit interval required")
        bounds = np.asarray(expected, float)
    else:
        bounds = np.asarray(require_positive_interval(interval), float)
        if physical and not np.allclose(bounds, expected, atol=3e-12, rtol=0):
            raise ValueError("interval is the fixed I=S(1)+[.12,.18]")
    if grid_interval is not None:
        grid = np.asarray(require_positive_interval(grid_interval), float)
        if physical and not np.allclose(grid, expected, atol=3e-12, rtol=0):
            raise ValueError("family interval must be the declared I")
        return (float(grid[0]), float(grid[1]))
    return (float(bounds[0]), float(bounds[1]))


def require_family(family, *, profile=PROFILE_IDENTITY):
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required as the candidate/grid owner")
    identity = profile_identity(family, include_normal_window=True)
    if identity != profile:
        raise ValueError("profile identity mismatch")
    interval = require_interval(family.interval, grid_interval=family.interval)
    return family, identity, interval


def require_verification_nodes(nodes, interval, family=None, *, physical=None):
    """Exactly 257 strictly increasing nodes covering I; the family owns the grid."""
    array = np.asarray(nodes, float)
    if array.ndim != 1 or array.shape != (VERIFICATION_NODE_COUNT,):
        raise ValueError("exactly 257 verification nodes required")
    if not np.isfinite(array).all():
        raise ValueError("finite verification nodes required")
    if np.any(np.diff(array) <= 0):
        raise ValueError("strictly increasing verification nodes required")
    if physical is None:
        physical = family is not None
    if family is not None:
        owned = np.asarray(family.collocation_nodes(VERIFICATION_NODE_COUNT), float)
        if not np.array_equal(array, owned):
            raise ValueError(
                "verification nodes must be LocalIncomingFamily.collocation_nodes(257)")
        require_interval(interval, grid_interval=family.interval, physical=physical)
    else:
        left, right = require_interval(interval, physical=physical)
        if array[0] != left or array[-1] != right:
            raise ValueError("verification nodes must cover I with both endpoints")
    frozen = np.array(array, dtype=float, copy=True)
    frozen.setflags(write=False)
    return frozen


def _pair_array(values, shape, name):
    if values is None:
        raise ValueError("explicit " + name + " required; missing cells are not zero")
    if isinstance(values, dict):
        raise ValueError(name + " must be a complete array; partial maps are rejected")
    array = np.asarray(values, object)
    if array.shape != shape:
        raise ValueError(name + " must have shape " + str(shape))
    result = np.empty(shape, dtype=object)
    for index in np.ndindex(shape):
        item = array[index]
        if item is None:
            raise ValueError("missing " + name + " at " + str(index))
        result[index] = _nonneg_fraction(item, name)
    return result


def require_nodal_abs_enclosures(values):
    return _pair_array(values, (VERIFICATION_NODE_COUNT, 2), "nodal absolute-value enclosure")


def require_cell_derivative_uppers(values):
    return _pair_array(values, (CELL_COUNT, 2), "cell derivative upper")


def authenticate_binding(binding, *, family=None, source_hashes=None, root=None):
    if not isinstance(binding, EnclosureBinding):
        raise TypeError("EnclosureBinding required")
    if family is not None:
        require_family(family, profile=binding.profile_identity)
    actual = {}
    for name, digest in binding.dependency_hashes:
        _, hashed = authenticate_file_digest(
            name, digest, root=root, name=name)
        actual[name] = hashed
    if source_hashes is not None:
        expected = dict(binding.dependency_hashes)
        for name, digest in expected.items():
            if name not in source_hashes or source_hashes[name] != digest:
                raise ValueError("dependency hash mismatch: " + name)
            if actual[name] != digest:
                raise ValueError("dependency hash mismatch: " + name)
        extra = sorted(set(source_hashes) - set(expected))
        if extra:
            raise ValueError("undeclared dependency hash: " + extra[0])
    return binding


def enclose_continuous_constraints(nodes, nodal_abs, derivative_uppers, *,
                                   interval, family=None, binding=None,
                                   samples=None, samples_used_as_bound=False,
                                   source_hashes=None, root=None):
    """Directed Lipschitz enclosure of both raw constraints on the 257-node grid.

    ``samples`` may be stored as diagnostics. They never replace nodal
    absolute-value enclosures or cell derivative bounds.
    """
    if samples_used_as_bound:
        raise ValueError("samples are not a continuous bound")
    if binding is not None:
        authenticate_binding(
            binding, family=family, source_hashes=source_hashes, root=root)
        interval = binding.interval
    grid = require_verification_nodes(
        nodes, interval, family=family,
        physical=family is not None or binding is not None)
    values = require_nodal_abs_enclosures(nodal_abs)
    derivatives = require_cell_derivative_uppers(derivative_uppers)
    if samples is not None:
        sample = np.asarray(samples, float)
        if sample.shape != (VERIFICATION_NODE_COUNT, 2) or not np.isfinite(sample).all():
            raise ValueError("diagnostic samples must be finite with shape (257, 2)")
        if np.any(sample < 0):
            raise ValueError("diagnostic samples must be nonnegative")
    enclosed = verification_componentwise_continuous_max(grid, values, derivatives)
    return {
        "schema": SCHEMA,
        "node_count": VERIFICATION_NODE_COUNT,
        "cell_count": CELL_COUNT,
        "nodes_hex": [hex_float(value) for value in grid],
        "nodes_digest": array_digest(grid),
        "grid_owner": GRID_OWNER,
        "N": serialize_lipschitz(enclosed["N"]),
        "beta": serialize_lipschitz(enclosed["beta"]),
        "samples_are_not_proof": True,
        "samples_used_as_bound": False,
        "dense_sampling_used_as_bound": False,
        "certificate_use": True,
    }


def production_between_node_component(*, nodal_abs=None, derivative_uppers=None,
                                      nodes=None, interval=None, family=None,
                                      binding=None, samples=None,
                                      source_hashes=None, root=None):
    """Fill the between-node slot only from a complete directed enclosure."""
    missing = []
    if nodal_abs is None:
        missing.append(MISSING_PRODUCTION_INPUTS[0])
    if derivative_uppers is None:
        missing.append(MISSING_PRODUCTION_INPUTS[1])
    if missing:
        return {
            "bound": None,
            "status": (
                "OPEN: 257-node / 256-cell between-node enclosure lacks "
                "production nodal absolute-value and derivative inputs"),
            "missing_inputs": tuple(missing),
            "samples_are_not_proof": True,
            "certificate_use": False,
            "physical_EXISTENCE_certificate": False,
            "physical_NONEXISTENCE_certificate": False,
        }
    enclosed = enclose_continuous_constraints(
        nodes, nodal_abs, derivative_uppers, interval=interval,
        family=family, binding=binding, samples=samples,
        source_hashes=source_hashes, root=root)
    return {
        "bound": [
            enclosed["N"]["outward_binary64"]["between_node_remainder_upper"],
            enclosed["beta"]["outward_binary64"]["between_node_remainder_upper"],
        ],
        "exact_upper_evidence": {
            "N": enclosed["N"]["between_node_remainder_upper"],
            "beta": enclosed["beta"]["between_node_remainder_upper"],
        },
        "continuous_maximum_upper": {
            "N": enclosed["N"]["continuous_maximum_upper"],
            "beta": enclosed["beta"]["continuous_maximum_upper"],
        },
        "status": "PASS: directed 257-node / 256-cell Lipschitz enclosure",
        "missing_inputs": (),
        "enclosure": enclosed,
        "samples_are_not_proof": True,
        "certificate_use": True,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
    }


def manufactured_linear_enclosures(nodes):
    """Exact |f| and |f'| for the affine map on I. Diagnostic only."""
    grid = np.asarray(nodes, float)
    left, right = float(grid[0]), float(grid[-1])
    width = Fraction(str(right)) - Fraction(str(left))
    if width <= 0:
        raise ValueError("positive manufactured interval required")
    slope = 1 / width
    nodal = np.empty((VERIFICATION_NODE_COUNT, 2), dtype=object)
    derivatives = np.empty((CELL_COUNT, 2), dtype=object)
    for index, z in enumerate(grid):
        t = (Fraction(str(float(z))) - Fraction(str(left))) / width
        nodal[index, 0] = abs(t)
        nodal[index, 1] = abs(1 - t)
    for cell in range(CELL_COUNT):
        derivatives[cell, 0] = slope
        derivatives[cell, 1] = slope
    return MappingProxyType({"nodal_abs": nodal, "derivative_uppers": derivatives})


def diagnostic_linear_enclosure(family, binding, *, source_hashes=None, root=None):
    """Finite synthetic enclosure on the owned grid. certificate_use is false."""
    family, _identity, interval = require_family(
        family, profile=binding.profile_identity)
    nodes = require_verification_nodes(
        family.collocation_nodes(VERIFICATION_NODE_COUNT), interval, family=family)
    manufactured = manufactured_linear_enclosures(nodes)
    enclosed = enclose_continuous_constraints(
        nodes, manufactured["nodal_abs"], manufactured["derivative_uppers"],
        interval=interval, family=family, binding=binding,
        source_hashes=source_hashes, root=root)
    enclosed = dict(enclosed)
    enclosed["certificate_use"] = False
    enclosed["scope"] = (
        "synthetic Lipschitz control on LocalIncomingFamily.collocation_nodes(257); "
        "not a residual enclosure and not a budget component")
    enclosed["status"] = (
        "DIAGNOSTIC: manufactured affine |f| on the owned 257-node grid")
    return enclosed


def validate_continuous_enclosure_record(record):
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("NSC-KS-CONTINUOUS-CONSTRAINT-ENCLOSURE-v1 required")
    if record.get("profile_identity") != PROFILE_IDENTITY:
        raise ValueError("enclosure profile identity mismatch")
    if record.get("state_law") != STATE_LAW:
        raise ValueError("enclosure state law mismatch")
    if record.get("verification_node_count") != VERIFICATION_NODE_COUNT:
        raise ValueError("exactly 257 verification nodes required")
    if record.get("cell_count") != CELL_COUNT:
        raise ValueError("exactly 256 cells required")
    if record.get("grid_owner") != GRID_OWNER:
        raise ValueError("LocalIncomingFamily.collocation_nodes owns the grid")
    production = record.get("production")
    if not isinstance(production, dict):
        raise ValueError("production enclosure block required")
    if production.get("between_node") not in (None,):
        if production.get("certificate_use") is not True:
            raise ValueError("a production between-node bound requires certificate_use")
    else:
        if production.get("certificate_use") is True:
            raise ValueError("null production bound cannot claim certificate_use")
    if production.get("arithmetic") is not None:
        raise ValueError("arithmetic is owned by the assembly-arithmetic component")
    if record.get("physical_EXISTENCE_certificate"):
        raise ValueError("OPEN enclosure may not be reported as EXISTENCE")
    if record.get("physical_NONEXISTENCE_certificate"):
        raise ValueError("OPEN enclosure may not be reported as NON_EXISTENCE")
    if record.get("samples_are_not_proof") is not True:
        raise ValueError("samples are not a continuous bound")
    diagnostic = record.get("diagnostic")
    if diagnostic is not None:
        if diagnostic.get("certificate_use") is not False:
            raise ValueError("diagnostic enclosure cannot be used as a certificate")
        if diagnostic.get("status", "").split(":", 1)[0] in ("PASS", "EXISTENCE"):
            raise ValueError("diagnostic enclosure cannot claim PASS or EXISTENCE")
    status = record.get("status")
    if not isinstance(status, str) or not status.startswith("OPEN"):
        raise ValueError("production continuous enclosure remains OPEN")
    missing = record.get("missing_scientific_inputs")
    if not isinstance(missing, list) or not missing:
        raise ValueError("missing scientific inputs must remain explicit")
    if set(record.get("error_components_not_filled", ())) - set(ERROR_COMPONENTS):
        raise ValueError("unknown error-component name")
    return record
