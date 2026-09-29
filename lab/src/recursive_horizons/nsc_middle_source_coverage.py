"""Resumable source bounds for the original group14 angular+1 middle window.

Each completed positive row reuses the reviewed correction equation, its
nonzero horizon initial remainder, the finite thermal/coherent occupation
bound, and a separate archived comparison of the actual negative partner.
Missing rows stay out of the aggregate. Nothing here fills a full upstream
budget or replaces the saved source.
"""
from hashlib import sha256
from io import BytesIO
import json
import math
from pathlib import Path
import time

import numpy as np
from flint import arb, ctx

from .evidence_io import publish_exclusive_directory
from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_ks_evaluation_binding import implementation_hashes
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from .nsc_subgap_source_covariance import original_covariance_distance_bounds
from .nsc_vacuum_source_correction import (
    VacuumCorrection, capture_correction, validate_correction,
)
from .nsc_vacuum_source_remainder import VacuumSourceExpansion


SCHEMA = "NSC-MIDDLE-SOURCE-COVERAGE-v1"
ROW_SCHEMA = "NSC-MIDDLE-SOURCE-ROW-v1"
GROUP = 14
ANGULAR_SIGN = 1
PANEL = "group14/mid24_1"
ENERGY_LOWER = 16.0
ENERGY_UPPER = 32.0
EXPECTED_POSITIVE_ROWS = 48
START = -18.0
ORDER = 8
BITS = 192
METRIC_TERMS = 48
DEGREE = 12
RTOL = 1e-9
ATOL = 1e-21
MAX_STEP = 0.1
EQUATION = "e_y=A_E(y) e-r_M(y)"
ROW_STATUS = (
    "ENCLOSED: one original middle-energy preparation and its signed partner; "
    "full gate OPEN"
)
OPEN_STATUS = (
    "OPEN: group14 angular+1 middle-energy source coverage is incomplete; "
    "completed rows only; full upstream budget not filled"
)
WINDOW_STATUS = (
    "ENCLOSED: original group14 angular+1 middle energies 16<=E<32 and their "
    "signed partners; full upstream budget OPEN"
)
PAYLOAD_NAMES = (
    "negative_columns", "negative_covariance", "positive_columns",
    "positive_covariance", "trace",
)
OWNERS = (
    "docs/nsc-middle-source-coverage.md",
    "scripts/derive_nsc_middle_source_coverage.py",
    "scripts/derive_nsc_vacuum_source_correction.py",
    "src/recursive_horizons/nsc_middle_source_coverage.py",
    "src/recursive_horizons/nsc_vacuum_source_correction.py",
    "src/recursive_horizons/nsc_vacuum_source_remainder.py",
    "tests/test_nsc_middle_source_coverage.py",
)


def encode_record(record):
    return (json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def checkpoint_name(panel, row):
    if not isinstance(panel, str) or isinstance(row, bool) or type(row) is not int or row < 0:
        raise ValueError("panel and row identity required")
    slug = panel.replace("/", "__")
    if (not slug or not slug.isascii() or not all(character.isalnum() or character == "_"
                                                   for character in slug)):
        raise ValueError("unsafe panel checkpoint name")
    return f"row-{slug}-{row:05d}"


def checkpoint_relative(panel, row):
    return f"rows/{checkpoint_name(panel, row)}"


def in_middle_window(energy):
    value = float(energy)
    if not math.isfinite(value):
        raise ValueError("finite middle-window energy required")
    return ENERGY_LOWER <= value < ENERGY_UPPER


def complemented_endpoint(endpoint):
    """Signed partner map used by the reviewed correction, not a copied scalar bound."""
    if len(endpoint) != 3:
        raise ValueError("three-component Bloch endpoint required")
    return (endpoint[0], -endpoint[1], -endpoint[2])


def _interval(value):
    return {"lower": exact_upper(value.lower()), "upper": exact_upper(value.upper())}


def _frozen(value, dtype):
    array = np.array(value, dtype=dtype, copy=True, order="C")
    array.setflags(write=False)
    return array


def _budgets(row_budget, cpu_budget):
    if isinstance(row_budget, bool) or type(row_budget) is not int or row_budget < 0:
        raise ValueError("explicit nonnegative integer row budget required")
    if (isinstance(cpu_budget, bool) or not isinstance(cpu_budget, (int, float))
            or not math.isfinite(float(cpu_budget)) or float(cpu_budget) <= 0):
        raise ValueError("explicit positive finite CPU budget required")
    return row_budget, float(cpu_budget)


def _remember(seen, panel, row, energy_hex):
    key = (panel, row)
    energy_key = (panel, energy_hex)
    if key in seen or energy_key in seen:
        raise ValueError("duplicate source identity")
    seen.add(key)
    seen.add(energy_key)
    return key


class MiddleRow:
    """One original positive fiber plus the archived opposite-angular partner."""

    def __init__(self, **fields):
        self.__dict__.update(fields)
        for name in ("positive_columns", "negative_columns", "positive_covariance",
                     "negative_covariance"):
            array = getattr(self, name)
            if array.flags.writeable:
                raise ValueError("source arrays must be frozen before a row is stored")


def _row_from_batches(positive, partner, local, config):
    if positive.original_panel != PANEL or positive.group != GROUP or positive.angular_sign != 1:
        raise ValueError("middle window must stay on original group14/mid24_1")
    if partner.angular_sign != -1 or partner.energy_sign != -1:
        raise ValueError("actual opposite-angular negative partner required")
    if partner.original_panel != positive.original_panel or partner.rows != positive.rows:
        raise ValueError("negative partner must share the original panel rows")
    if positive.mass != np.pi/2 or positive.angular != np.sqrt(5):
        raise ValueError("original group14 mass and angular label required")
    if partner.mass != positive.mass or partner.angular != -positive.angular or partner.rho_up != positive.rho_up:
        raise ValueError("signed mass, angular label or rho changed")
    sl = slice(3*local, 3*local+3)
    fiber = np.asarray(positive.source.energies[sl], float)
    energy = float(fiber[0])
    if not np.all(fiber == energy) or not in_middle_window(energy):
        raise ValueError("middle-window three-column energy fiber required")
    negative_fiber = np.asarray(partner.source.energies[sl], float)
    if not np.array_equal(negative_fiber, -fiber):
        raise ValueError("negative partner energies must equal -E")
    positive_columns = _frozen(positive.initial_columns[:, sl], np.complex128)
    negative_columns = _frozen(partner.initial_columns[:, sl], np.complex128)
    positive_covariance = _frozen(positive.source.covariance[sl, sl], np.complex128)
    negative_covariance = _frozen(partner.source.covariance[sl, sl], np.complex128)
    if (np.array_equal(positive_columns, negative_columns)
            and np.array_equal(positive_covariance, negative_covariance)):
        raise ValueError("negative partner is a copy of the positive source")
    if (positive.source.digest == partner.source.digest
            or positive.preparation_digest == partner.preparation_digest):
        raise ValueError("duplicate source identity")
    return MiddleRow(
        panel=positive.original_panel, row=int(positive.rows[0]+local),
        batch_rows=(int(positive.rows[0]), int(positive.rows[1])), local_index=int(local),
        energy=energy, negative_energy=float(negative_fiber[0]), rho_up=float(positive.rho_up),
        mass=float(positive.mass), angular=float(positive.angular),
        positive_angular_sign=1, negative_angular_sign=-1,
        positive_source_digest=positive.source.digest,
        negative_source_digest=partner.source.digest,
        positive_preparation_digest=positive.preparation_digest,
        negative_preparation_digest=partner.preparation_digest,
        positive_columns=positive_columns, negative_columns=negative_columns,
        positive_covariance=positive_covariance, negative_covariance=negative_covariance,
        kappa=float(config["surface_gravity"]), omega=float(config["omega"]),
        horizon_rho=config["horizon_rho"])


def middle_catalogue(archive):
    """Original group14 angular+1 rows with 16<=E<32 and their negative partners."""
    if not isinstance(archive, RetainedUpstreamArchive):
        raise TypeError("retained upstream archive required")
    config = archive.meta["config"]
    for name in ("surface_gravity", "omega", "horizon_rho"):
        if name not in config:
            raise ValueError("original source scales required")
    entries = archive.family_entries((GROUP, ANGULAR_SIGN))
    positives = []
    for batch, _channel in entries:
        if batch.energy_sign <= 0:
            continue
        fibers = np.asarray(batch.source.energies, float).reshape(-1, 3)
        for local, fiber in enumerate(fibers):
            if in_middle_window(float(fiber[0])):
                positives.append((batch, local))
    if len(positives) != EXPECTED_POSITIVE_ROWS:
        raise ValueError("original group14 middle-energy row count changed")
    seen = set()
    rows = []
    for batch, local in positives:
        partners = [item for item, _channel in entries
                    if item.energy_sign < 0 and item.original_panel == batch.original_panel
                    and item.rows == batch.rows]
        if len(partners) != 1:
            raise ValueError("actual negative partner required")
        row = _row_from_batches(batch, partners[0], local, config)
        _remember(seen, row.panel, row.row, float(row.energy).hex())
        rows.append(row)
    rows.sort(key=lambda row: (row.energy, row.panel, row.row))
    if [row.row for row in rows] != list(range(EXPECTED_POSITIVE_ROWS)):
        raise ValueError("middle rows must be the original contiguous indices")
    if not in_middle_window(rows[0].energy) or not in_middle_window(rows[-1].energy):
        raise ValueError("middle energy window changed")
    return tuple(rows)


def dependency_identity(archive_root, archive):
    source_hashes = implementation_hashes(Path(archive_root), owners=OWNERS)
    inputs = dict(archive.input_hashes)
    if not source_hashes or not inputs:
        raise ValueError("source and dependency identity required")
    return {"source_hashes": source_hashes, "input_hashes": inputs}


def _flow(spec):
    # Same enclosure as the reviewed single-row driver: the float labels are
    # unioned with the exact mass and angular values, and only then expanded.
    expansion = VacuumSourceExpansion(
        spec.horizon_rho, arb(spec.mass).union(arb.pi()/2),
        arb(spec.angular).union(arb(5).sqrt()), order=ORDER, bits=BITS,
        metric_terms=METRIC_TERMS)
    flow = VacuumCorrection(expansion, spec.energy)
    target = float(expansion.target_distance(spec.rho_up).log().mid())
    return expansion, flow, target


def _physical(columns, covariance, validation, thermal):
    bounds = original_covariance_distance_bounds(columns, covariance, validation)
    lower = max(arb(0), (bounds["lower"]-thermal).lower())
    upper = (bounds["upper"]+thermal).upper()
    return bounds, lower, upper


def _comparisons(spec, proof, thermal):
    flipped = dict(proof)
    flipped["endpoint"] = complemented_endpoint(proof["endpoint"])
    plain, _lower, _upper = _physical(
        spec.negative_columns, spec.negative_covariance, proof, thermal)
    signed, _lower, _upper = _physical(
        spec.negative_columns, spec.negative_covariance, flipped, thermal)
    # The positive Bloch endpoint is the wrong comparison for the negative columns.
    if (exact_upper(plain["lower"]) == exact_upper(signed["lower"])
            and exact_upper(plain["upper"]) == exact_upper(signed["upper"])):
        raise ValueError("signed complement did not change the archived negative comparison")
    rows = []
    for sign, validation, columns, covariance, source_digest, preparation_digest, angular in (
            (1, proof, spec.positive_columns, spec.positive_covariance,
             spec.positive_source_digest, spec.positive_preparation_digest, 1),
            (-1, flipped, spec.negative_columns, spec.negative_covariance,
             spec.negative_source_digest, spec.negative_preparation_digest, -1)):
        _distance, lower, upper = _physical(columns, covariance, validation, thermal)
        rows.append({
            "angular_sign": angular, "energy_sign": sign,
            "physical_source_operator_error_lower": exact_upper(lower),
            "physical_source_operator_error_upper": exact_upper(upper),
            "preparation_digest": preparation_digest, "source_digest": source_digest,
        })
    if rows[0]["energy_sign"] == rows[1]["energy_sign"] or rows[0]["source_digest"] == rows[1]["source_digest"]:
        raise ValueError("duplicate source identity")
    return rows


def _prove(spec, expansion, flow, target, trace, nfev, dependencies):
    horizon = flow.initial_error(START)
    if not horizon.is_finite() or not horizon > 0:
        raise ArithmeticError("nonzero finite horizon initial remainder required")
    proof = validate_correction(flow, trace, START, spec.rho_up, degree=DEGREE)
    if not proof["bloch_error"].is_finite() or not proof["bloch_error"] >= horizon:
        raise ArithmeticError("vacuum error dropped the horizon remainder")
    # Constant zero keeps the curve's vacuum error from being added twice.
    # The occupation law is still the finite thermal/coherent remainder.
    thermal = expansion.covariance_error(
        0, spec.energy, kappa=spec.kappa, omega=spec.omega)["thermal_coherent_error"]
    if not thermal.is_finite() or not thermal > 0:
        raise ArithmeticError("finite positive thermal/coherent remainder required")
    comparisons = _comparisons(spec, proof, thermal)
    payload = deterministic_npz_bytes({
        "trace": np.ascontiguousarray(trace, dtype=np.float64),
        "positive_columns": spec.positive_columns,
        "negative_columns": spec.negative_columns,
        "positive_covariance": spec.positive_covariance,
        "negative_covariance": spec.negative_covariance,
    })
    relative = checkpoint_relative(spec.panel, spec.row)
    record = {
        "schema": ROW_SCHEMA, "coverage_schema": SCHEMA, "status": ROW_STATUS,
        "group": GROUP, "panel": spec.panel, "row": spec.row,
        "batch_rows": [spec.batch_rows[0], spec.batch_rows[1]],
        "local_index": spec.local_index,
        "energy_hex": float(spec.energy).hex(),
        "negative_energy_hex": float(spec.negative_energy).hex(),
        "angular_sign": 1, "negative_angular_sign": -1,
        "rho_up_hex": float(spec.rho_up).hex(),
        "horizon_rho_hex": float(spec.horizon_rho).hex(),
        "kappa_hex": float(spec.kappa).hex(), "omega_hex": float(spec.omega).hex(),
        "start_log_delta_hex": float(START).hex(), "end_log_delta_hex": float(target).hex(),
        "expansion_order": ORDER, "defect_degree": DEGREE, "bits": BITS,
        "metric_terms": METRIC_TERMS, "rtol_hex": float(RTOL).hex(),
        "atol_hex": float(ATOL).hex(), "max_step_hex": float(MAX_STEP).hex(),
        "equation": EQUATION, "cells": int(proof["cells"]), "nfev": int(nfev),
        "horizon_initial_remainder_upper": exact_upper(horizon),
        "initial_bloch_error_upper": exact_upper(proof["initial_error"]),
        "normalized_defect_integral_upper": exact_upper(proof["normalized_defect_integral"]),
        "endpoint_bridge_upper": exact_upper(proof["endpoint_bridge"]),
        "vacuum_bloch_error_upper": exact_upper(proof["bloch_error"]),
        "thermal_coherent_operator_error_upper": exact_upper(thermal),
        "validated_vacuum_bloch_center": [_interval(value) for value in proof["endpoint"]],
        "source_comparisons": comparisons,
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "covariance_errors_unweighted": True,
        "initial_true_correction_assumed_zero": False,
        "horizon_initial_remainder_nonzero": True,
        "negative_comparison_uses_complemented_endpoint": True,
        "positive_scalar_bound_copied_to_negative": False,
        "replay_uses_saved_witness": True, "archived_source_replaced": False,
        "source_occupations_changed": False, "original_preparation_ode_rerun": False,
        "proof_correction_ode_evolved": True, "source_quadrature_error_included": False,
        "all_source_families": False, "other_energies_filled_from_this_row": False,
        "physical_upstream_budget_component": None, "changed_history_C_M": None,
        "physical_local_gate": "OPEN",
        "payload": {"path": f"{relative}/witness.npz", "sha256": sha256(payload).hexdigest(),
                    "bytes": len(payload)},
        "source_hashes": dependencies["source_hashes"],
        "input_hashes": dependencies["input_hashes"],
    }
    return record, payload


def _require(record, field, expected, message):
    if record.get(field) != expected:
        raise ValueError(message)


def _require_settings(record, spec):
    message = "settings changed"
    _require(record, "schema", ROW_SCHEMA, "middle-source checkpoint schema required")
    _require(record, "coverage_schema", SCHEMA, "middle-source checkpoint schema required")
    _require(record, "group", GROUP, message)
    _require(record, "angular_sign", 1, message)
    _require(record, "negative_angular_sign", -1, message)
    _require(record, "start_log_delta_hex", float(START).hex(), message)
    _require(record, "expansion_order", ORDER, message)
    _require(record, "defect_degree", DEGREE, message)
    _require(record, "bits", BITS, message)
    _require(record, "metric_terms", METRIC_TERMS, message)
    _require(record, "rtol_hex", float(RTOL).hex(), message)
    _require(record, "atol_hex", float(ATOL).hex(), message)
    _require(record, "max_step_hex", float(MAX_STEP).hex(), message)
    _require(record, "equation", EQUATION, message)
    _require(record, "rho_up_hex", float(spec.rho_up).hex(), message)
    _require(record, "horizon_rho_hex", float(spec.horizon_rho).hex(), message)
    _require(record, "kappa_hex", float(spec.kappa).hex(), message)
    _require(record, "omega_hex", float(spec.omega).hex(), message)
    _require(record, "initial_true_correction_assumed_zero", False, message)
    _require(record, "positive_scalar_bound_copied_to_negative", False, message)


def _require_identity(record, spec):
    message = "checkpoint identity does not match the original middle row"
    _require(record, "panel", spec.panel, message)
    _require(record, "row", spec.row, message)
    _require(record, "local_index", spec.local_index, message)
    _require(record, "batch_rows", [spec.batch_rows[0], spec.batch_rows[1]], message)
    _require(record, "energy_hex", float(spec.energy).hex(), message)
    _require(record, "negative_energy_hex", float(spec.negative_energy).hex(), message)
    _require(record, "positive_source_digest", spec.positive_source_digest, "original source changed")
    _require(record, "negative_source_digest", spec.negative_source_digest, "original source changed")
    _require(record, "positive_preparation_digest", spec.positive_preparation_digest,
             "original source changed")
    _require(record, "negative_preparation_digest", spec.negative_preparation_digest,
             "original source changed")


def _require_dependencies(record, dependencies):
    if (record.get("source_hashes") != dependencies["source_hashes"]
            or record.get("input_hashes") != dependencies["input_hashes"]):
        raise ValueError("dependency identity changed")


def _load_payload(record, payload, spec):
    described = record.get("payload")
    if not isinstance(described, dict):
        raise ValueError("missing payload")
    digest = sha256(payload).hexdigest()
    if described.get("sha256") != digest or described.get("bytes") != len(payload):
        raise ValueError("payload hash changed")
    expected_path = f"{checkpoint_relative(spec.panel, spec.row)}/witness.npz"
    if described.get("path") != expected_path:
        raise ValueError("payload hash changed")
    try:
        with np.load(BytesIO(payload), allow_pickle=False) as data:
            if set(data.files) != set(PAYLOAD_NAMES):
                raise ValueError("missing payload")
            loaded = {name: np.array(data[name], copy=True) for name in PAYLOAD_NAMES}
    except ValueError:
        raise
    except Exception as error:
        raise ValueError("missing payload") from error
    originals = {
        "positive_columns": spec.positive_columns,
        "negative_columns": spec.negative_columns,
        "positive_covariance": spec.positive_covariance,
        "negative_covariance": spec.negative_covariance,
    }
    for name, value in originals.items():
        if not np.array_equal(loaded[name], value):
            raise ValueError("original source changed")
    trace = loaded["trace"]
    if trace.ndim != 2 or trace.shape[1] != 26 or not np.isfinite(trace).all():
        raise ValueError("missing payload")
    return trace


def authenticate_checkpoint(spec, dirname, record, payload, dependencies):
    """Replay one saved witness. This never calls the correction solver."""
    if dirname != checkpoint_name(spec.panel, spec.row):
        raise ValueError("checkpoint name does not match row identity")
    if not isinstance(record, dict):
        raise ValueError("middle-source checkpoint schema required")
    _require_settings(record, spec)
    _require_identity(record, spec)
    _require_dependencies(record, dependencies)
    trace = _load_payload(record, payload, spec)
    nfev = record.get("nfev")
    if isinstance(nfev, bool) or type(nfev) is not int or nfev < 0:
        raise ValueError("settings changed")
    with ctx.workprec(BITS):
        expansion, flow, target = _flow(spec)
        rebuilt, rebuilt_payload = _prove(
            spec, expansion, flow, target, trace, nfev, dependencies)
    if rebuilt_payload != payload or encode_record(rebuilt) != encode_record(record):
        raise ValueError("saved witness changed")
    return rebuilt


def _canonical_directory(path, *, create):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError("absolute coverage path required")
    if create:
        path.mkdir(parents=True, exist_ok=True)
    elif not path.exists():
        return path.resolve(strict=False)
    return path.resolve(strict=True)


def publish_checkpoint(output_root, record, payload):
    """Atomically publish one new row directory. An existing row is left untouched."""
    name = checkpoint_name(record["panel"], record["row"])
    described = record.get("payload") if isinstance(record, dict) else None
    if not isinstance(described, dict):
        raise ValueError("missing payload")
    digest = sha256(payload).hexdigest()
    if described.get("sha256") != digest or described.get("bytes") != len(payload):
        raise ValueError("payload hash changed")
    if described.get("path") != f"rows/{name}/witness.npz":
        raise ValueError("payload hash changed")
    root = _canonical_directory(output_root, create=True)
    (root/"rows").mkdir(exist_ok=True)
    publish_exclusive_directory(root, f"rows/{name}", {
        "record.json": encode_record(record), "witness.npz": payload,
    })


def _read_checkpoint(path):
    record_path = path/"record.json"
    payload_path = path/"witness.npz"
    if not record_path.is_file() or not payload_path.is_file():
        raise ValueError("missing payload")
    raw = record_path.read_bytes()
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError("missing payload") from error
    if encode_record(record) != raw:
        raise ValueError("settings changed")
    return record, payload_path.read_bytes()


def load_checkpoints(output_root):
    root = _canonical_directory(output_root, create=False)
    rows = root/"rows"
    if not rows.exists():
        return []
    if not rows.is_dir():
        raise ValueError("missing payload")
    loaded = []
    seen = set()
    for path in sorted(rows.iterdir(), key=lambda item: item.name):
        if not path.is_dir() or not path.name.startswith("row-"):
            raise ValueError("unexpected coverage entry")
        record, payload = _read_checkpoint(path)
        panel, row, energy_hex = _row_key(record)
        _remember(seen, panel, row, energy_hex)
        loaded.append((path.name, record, payload))
    return loaded


def _row_key(record):
    if not isinstance(record, dict):
        raise ValueError("checkpoint row identity is incomplete")
    panel, row, energy = record.get("panel"), record.get("row"), record.get("energy_hex")
    if (not isinstance(panel, str) or isinstance(row, bool) or type(row) is not int
            or not isinstance(energy, str)):
        raise ValueError("checkpoint row identity is incomplete")
    return panel, row, energy


def _summary(record):
    return {
        "panel": record["panel"], "row": record["row"],
        "energy_hex": record["energy_hex"],
        "negative_energy_hex": record["negative_energy_hex"],
        "checkpoint": record["payload"]["path"].rsplit("/", 1)[0],
        "cells": record["cells"], "nfev": record["nfev"],
        "equation": record["equation"],
        "horizon_initial_remainder_upper": record["horizon_initial_remainder_upper"],
        "vacuum_bloch_error_upper": record["vacuum_bloch_error_upper"],
        "thermal_coherent_operator_error_upper": record["thermal_coherent_operator_error_upper"],
        "source_comparisons": record["source_comparisons"],
        "physical_upstream_budget_component": None, "physical_local_gate": "OPEN",
    }


def aggregate(catalogue, completed, *, mode, row_budget, cpu_budget, cpu_seconds,
              capture_cpu_seconds, new_rows_captured):
    by_key = {(record["panel"], record["row"]): record for record in completed}
    if len(by_key) != len(completed):
        raise ValueError("duplicate source identity")
    ordered = [by_key[(spec.panel, spec.row)] for spec in catalogue if (spec.panel, spec.row) in by_key]
    missing = [{"panel": spec.panel, "row": spec.row, "energy_hex": float(spec.energy).hex()}
               for spec in catalogue if (spec.panel, spec.row) not in by_key]
    worst = None
    with ctx.workprec(BITS):
        for record in ordered:
            for comparison in record["source_comparisons"]:
                upper = restored_upper(comparison["physical_source_operator_error_upper"])
                worst = upper if worst is None else max(worst, upper)
    complete = len(ordered) == len(catalogue) and not missing
    return {
        "schema": SCHEMA,
        "status": WINDOW_STATUS if complete else OPEN_STATUS,
        "group": GROUP, "angular_sign": ANGULAR_SIGN, "panel": PANEL,
        "energy_window": {"lower": ENERGY_LOWER, "lower_included": True,
                          "upper": ENERGY_UPPER, "upper_included": False},
        "expected_positive_rows": EXPECTED_POSITIVE_ROWS,
        "expected_negative_partners": EXPECTED_POSITIVE_ROWS,
        "completed_positive_rows": len(ordered),
        "completed_negative_partners": len(ordered),
        "completed_signed_rows": 2*len(ordered),
        "missing_positive_rows": len(missing),
        "coverage_complete": complete,
        "interrupted": not complete,
        "mode": mode,
        "rows": [_summary(record) for record in ordered],
        "missing": missing,
        "maximum_completed_source_operator_error_upper":
            None if worst is None else exact_upper(worst),
        "covariance_errors_unweighted": True,
        "source_quadrature_error_included": False,
        "all_source_families": False,
        "archived_source_replaced": False,
        "source_occupations_changed": False,
        "other_energies_filled_from_completed_rows": False,
        "physical_upstream_budget_component": None,
        "changed_history_C_M": None,
        "physical_local_gate": "OPEN",
        "row_budget": row_budget, "cpu_budget_seconds": cpu_budget,
        "new_rows_captured": new_rows_captured,
        "cpu_seconds": cpu_seconds, "capture_cpu_seconds": capture_cpu_seconds,
        "replay_uses_saved_witness": True,
    }


def _authenticate_saved(output_root, catalogue, dependencies):
    by_key = {(spec.panel, spec.row): spec for spec in catalogue}
    authenticated = []
    for dirname, record, payload in load_checkpoints(output_root):
        panel, row, _energy = _row_key(record)
        spec = by_key.get((panel, row))
        if spec is None:
            raise ValueError("checkpoint identity is not an original middle row")
        authenticated.append(authenticate_checkpoint(
            spec, dirname, record, payload, dependencies))
    return authenticated


def cover(archive_root, output_root, *, mode, row_budget=None, cpu_budget=None):
    """Validate saved rows before any new capture. Check mode never solves."""
    if mode not in {"resume", "check"}:
        raise ValueError("coverage mode must be resume or check")
    if mode == "check":
        if row_budget is not None or cpu_budget is not None:
            raise ValueError("check replays saved witnesses and does not accept a solve budget")
    else:
        row_budget, cpu_budget = _budgets(row_budget, cpu_budget)
    started = time.process_time()
    archive_root = Path(archive_root).resolve(strict=True)
    output_root = Path(output_root)
    if not output_root.is_absolute():
        output_root = archive_root/output_root
    archive = RetainedUpstreamArchive(archive_root)
    catalogue = middle_catalogue(archive)
    dependencies = dependency_identity(archive_root, archive)
    if mode == "check" and not Path(output_root).exists():
        return aggregate(catalogue, [], mode=mode, row_budget=None, cpu_budget=None,
                         cpu_seconds=time.process_time()-started, capture_cpu_seconds=0.0,
                         new_rows_captured=0)
    output_root = _canonical_directory(output_root, create=mode == "resume")
    # Saved witnesses are proved before a budget is allowed to start a solve.
    authenticated = _authenticate_saved(output_root, catalogue, dependencies)
    done = {(record["panel"], record["row"]) for record in authenticated}
    captured = 0
    capture_started = time.process_time()
    if mode == "resume":
        for spec in catalogue:
            if (spec.panel, spec.row) in done:
                continue
            if captured >= row_budget:
                break
            remaining = cpu_budget-(time.process_time()-capture_started)
            if remaining <= 0:
                break
            try:
                record, payload = _capture(spec, dependencies, cpu_limit=remaining)
            except TimeoutError:
                break
            publish_checkpoint(output_root, record, payload)
            authenticated.append(record)
            done.add((spec.panel, spec.row))
            captured += 1
    return aggregate(catalogue, authenticated, mode=mode, row_budget=row_budget,
                     cpu_budget=cpu_budget, cpu_seconds=time.process_time()-started,
                     capture_cpu_seconds=time.process_time()-capture_started if mode == "resume" else 0.0,
                     new_rows_captured=captured)


def _capture(spec, dependencies, *, cpu_limit):
    if not math.isfinite(cpu_limit) or cpu_limit <= 0:
        raise TimeoutError("vacuum-correction proof-witness CPU budget exhausted")
    with ctx.workprec(BITS):
        expansion, flow, target = _flow(spec)
        trace, nfev = capture_correction(
            flow, START, target, rtol=RTOL, atol=ATOL, max_step=MAX_STEP,
            cpu_limit=float(cpu_limit))
        return _prove(spec, expansion, flow, target, trace, int(nfev), dependencies)
