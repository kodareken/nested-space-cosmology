"""Resumable direct-vacuum bounds for one original low-energy window.

The transport is the reviewed homogeneous Bloch owner. Correction forcing is
not used. Each completed positive row keeps the archived distance and the
analytic occupation separate, then records each signed total as their sum
once. Missing rows stay out of the aggregate. This finite window does not
close the upstream budget or the local gate.
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
from .nsc_direct_vacuum_source import (
    DirectVacuumBloch, archived_vacuum_distance, capture_direct_vacuum,
    compare_signed_archives, validate_direct_vacuum,
)
from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_ks_evaluation_binding import implementation_hashes
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


SCHEMA = "NSC-DIRECT-SOURCE-WINDOW-v1"
ROW_SCHEMA = "NSC-DIRECT-SOURCE-ROW-v1"
GROUP = 14
ANGULAR_SIGN = 1
LOW16_PANEL = "group14/low16_1"
LOW32_PANEL = "group14/low32_1"
ENERGY_LOWER = 2.0
ENERGY_UPPER = 16.0
MIN_ENERGY = 2.071557525252559
MAX_ENERGY = 15.9788018699833
EXPECTED_POSITIVE_ROWS = 87
LOW16_ROWS = tuple(range(16, 48))
LOW32_ROWS = tuple(range(41, 96))
START = -18.0
BITS = 192
FRAME_ORDER = 16
METRIC_TERMS = 48
DEGREE = 12
MAX_STEP = 0.0125
# Fixed inside capture_direct_vacuum. This owner does not pass or override them.
RTOL = 2e-13
ATOL = 2e-15
ANALYTIC_RADIUS = "0.1"
METHOD = "DOP853"
EQUATION = "n_y=2*cross(h(y),n)"
MASS_ENCLOSURE = "archived-float union pi/2"
ANGULAR_ENCLOSURE = "archived-float union sqrt(5)"
ROW_STATUS = (
    "OPEN: one original direct-vacuum row and its negative partner; "
    "finite window only; local gate open"
)
OPEN_STATUS = (
    "OPEN: group14 angular+1 direct-vacuum window 2<=E<16 is incomplete; "
    "completed rows only; no missing row is entered as zero; "
    "full upstream budget and local gate remain open"
)
WINDOW_STATUS = (
    "OPEN: finite group14 angular+1 direct-vacuum window 2<=E<16 and the "
    "signed partners are recorded; full upstream budget and local gate remain open"
)
PAYLOAD_NAMES = (
    "negative_columns", "negative_covariance", "positive_columns",
    "positive_covariance", "trace",
)
OWNERS = (
    "docs/nsc-direct-source-window.md",
    "scripts/derive_nsc_direct_source_window.py",
    "src/recursive_horizons/nsc_direct_source_window.py",
    "src/recursive_horizons/nsc_direct_vacuum_source.py",
    "tests/test_nsc_direct_source_window.py",
)
_CAPTURE_CONSTANTS = capture_direct_vacuum.__code__.co_consts


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


def in_direct_window(energy):
    value = float(energy)
    if not math.isfinite(value):
        raise ValueError("finite direct-window energy required")
    return ENERGY_LOWER <= value < ENERGY_UPPER


def _hex(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("finite source value required")
    return number.hex()


def _fiber_hex(values):
    fiber = np.asarray(values, float)
    if fiber.shape != (3,) or not np.isfinite(fiber).all():
        raise ValueError("three-entry energy fiber required")
    return [_hex(item) for item in fiber]


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


def _require_owner_tolerances():
    if (METHOD not in _CAPTURE_CONSTANTS or RTOL not in _CAPTURE_CONSTANTS
            or ATOL not in _CAPTURE_CONSTANTS):
        raise ValueError("direct-owner capture tolerances changed")


class DirectRow:
    """One original positive fiber plus the archived opposite-angular partner."""

    def __init__(self, **fields):
        self.__dict__.update(fields)
        for name in ("positive_columns", "negative_columns", "positive_covariance",
                     "negative_covariance"):
            array = getattr(self, name)
            if array.flags.writeable:
                raise ValueError("source arrays must be frozen before a row is stored")


def _row_from_batches(positive, partner, local, config, archive):
    if positive.group != GROUP or positive.angular_sign != 1:
        raise ValueError("direct window must stay on original group14 angular+1")
    if partner.angular_sign != -1 or partner.energy_sign != -1:
        raise ValueError("actual opposite-angular negative partner required")
    if partner.original_panel != positive.original_panel or partner.rows != positive.rows:
        raise ValueError("negative partner must share the original panel rows")
    if positive.mass != np.pi/2 or positive.angular != np.sqrt(5):
        raise ValueError("original group14 mass and angular label required")
    if (partner.mass != positive.mass or partner.angular != -positive.angular
            or partner.rho_up != positive.rho_up):
        raise ValueError("signed mass, angular label or rho changed")
    row_index = int(positive.rows[0] + local)
    sl = slice(3*local, 3*local+3)
    fiber = np.asarray(positive.source.energies[sl], float)
    negative_fiber = np.asarray(partner.source.energies[sl], float)
    if fiber.shape != (3,) or not np.all(fiber == fiber[0]) or not in_direct_window(float(fiber[0])):
        raise ValueError("direct-window three-column energy fiber required")
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
    column_weights = np.asarray(positive.source.column_weights[sl], float)
    negative_weights = np.asarray(partner.source.column_weights[sl], float)
    if column_weights.shape != (3,) or not np.array_equal(column_weights, negative_weights):
        raise ValueError("signed column weights changed")
    if not np.all(column_weights == column_weights[0]) or not column_weights[0] > 0:
        raise ValueError("repeated positive column weight required")
    raw_weight = float(archive._arrays[positive.original_panel + "/weights"][row_index])
    if not math.isfinite(raw_weight) or raw_weight <= 0:
        raise ValueError("positive archived quadrature weight required")
    if float(column_weights[0]) != float(np.sqrt(raw_weight / (2 * np.pi))):
        raise ValueError("stored column weight does not match the archived quadrature weight")
    return DirectRow(
        panel=positive.original_panel, row=row_index,
        batch_rows=(int(positive.rows[0]), int(positive.rows[1])), local_index=int(local),
        energy=float(fiber[0]), negative_energy=float(negative_fiber[0]),
        energy_fiber=tuple(float(item) for item in fiber),
        negative_energy_fiber=tuple(float(item) for item in negative_fiber),
        rho_up=float(positive.rho_up), mass=float(positive.mass),
        angular=float(positive.angular), negative_angular=float(partner.angular),
        positive_angular_sign=1, negative_angular_sign=-1,
        positive_source_digest=positive.source.digest,
        negative_source_digest=partner.source.digest,
        positive_preparation_digest=positive.preparation_digest,
        negative_preparation_digest=partner.preparation_digest,
        positive_columns=positive_columns, negative_columns=negative_columns,
        positive_covariance=positive_covariance, negative_covariance=negative_covariance,
        weight=raw_weight, column_weight=float(column_weights[0]),
        kappa=float(config["surface_gravity"]), omega=float(config["omega"]),
        horizon_rho=float(config["horizon_rho"]))


def _require_census(rows):
    if len(rows) != EXPECTED_POSITIVE_ROWS:
        raise ValueError("original group14 direct-window row count changed")
    low16 = tuple(sorted(row.row for row in rows if row.panel == LOW16_PANEL))
    low32 = tuple(sorted(row.row for row in rows if row.panel == LOW32_PANEL))
    if low16 != LOW16_ROWS:
        raise ValueError("group14/low16_1 direct-window rows changed")
    if low32 != LOW32_ROWS:
        raise ValueError("group14/low32_1 direct-window rows changed")
    if {row.panel for row in rows} != {LOW16_PANEL, LOW32_PANEL}:
        raise ValueError("direct-window panels changed")
    energies = [row.energy for row in rows]
    if energies != sorted(energies) or len(set(energies)) != len(energies):
        raise ValueError("duplicate source identity")
    if rows[0].energy != MIN_ENERGY or rows[-1].energy != MAX_ENERGY:
        raise ValueError("direct-window energy endpoints changed")
    if any(not in_direct_window(row.energy) for row in rows):
        raise ValueError("direct-window energy endpoints changed")


def direct_catalogue(archive):
    """Original group14 angular+1 rows with 2<=E<16 and their negative partners."""
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
            if in_direct_window(float(fiber[0])):
                positives.append((batch, local))
    seen = set()
    rows = []
    for batch, local in positives:
        partners = [item for item, _channel in entries
                    if item.energy_sign < 0 and item.original_panel == batch.original_panel
                    and item.rows == batch.rows]
        if len(partners) != 1:
            raise ValueError("actual negative partner required")
        row = _row_from_batches(batch, partners[0], local, config, archive)
        _remember(seen, row.panel, row.row, _hex(row.energy))
        rows.append(row)
    rows.sort(key=lambda row: (row.energy, row.panel, row.row))
    _require_census(tuple(rows))
    return tuple(rows)


def dependency_identity(archive_root, archive):
    source_hashes = implementation_hashes(Path(archive_root), owners=OWNERS)
    inputs = dict(archive.input_hashes)
    if not source_hashes or not inputs:
        raise ValueError("source and dependency identity required")
    return {"source_hashes": source_hashes, "archive_digests": inputs}


def _model(spec):
    # Same enclosure as the reviewed direct-vacuum row: float labels union the
    # exact group14 mass and angular values. The archived arrays stay unchanged.
    mass = arb(spec.mass).union(arb.pi()/2)
    angular = arb(spec.angular).union(arb(5).sqrt())
    return DirectVacuumBloch(
        spec.energy, mass, angular, spec.horizon_rho, bits=BITS,
        frame_order=FRAME_ORDER, metric_terms=METRIC_TERMS,
        analytic_radius=ANALYTIC_RADIUS)


def _pack_occupation(occupation):
    required = {
        "coherence_retained": True,
        "signed_source_norm_equal": True,
        "requires_normalized_sewing": True,
        "source_occupation_law_changed": False,
        "field_or_quadrature_error_included": False,
        "physical_local_gate": "OPEN",
    }
    for key, expected in required.items():
        if occupation.get(key) != expected:
            raise ValueError("source occupation law changed")
    allowance = occupation["operator_distance_upper"]
    if not allowance.is_finite() or not allowance > 0:
        raise ArithmeticError("positive analytic occupation allowance required")
    return {
        "operator_distance_upper": exact_upper(allowance),
        "horizon_occupation_upper": exact_upper(occupation["horizon_occupation_upper"]),
        "incoming_occupation_upper": exact_upper(occupation["incoming_occupation_upper"]),
        "incoming_gap_used": bool(occupation["incoming_gap_used"]),
        "coherence_retained": True,
        "signed_source_norm_equal": True,
        "requires_normalized_sewing": True,
        "source_occupation_law_changed": False,
        "field_or_quadrature_error_included": False,
        "physical_local_gate": "OPEN",
    }


def _comparison(sign, angular, fiber, source_digest, preparation_digest, archived,
               occupation_upper):
    if not archived["lower"].is_finite() or not archived["upper"].is_finite():
        raise ArithmeticError("finite archived distance required")
    if not archived["upper"] >= archived["lower"]:
        raise ArithmeticError("archived distance bounds are reversed")
    # Sum the serialized dyadics once: [L, U] + [0, Occ] = [L, U+Occ].
    # Occ is not written back into the archived interval.
    archived_lower = exact_upper(archived["lower"])
    archived_upper = exact_upper(archived["upper"])
    total_upper = (restored_upper(archived_upper)+restored_upper(occupation_upper)).upper()
    if not total_upper > restored_upper(archived_upper):
        raise ArithmeticError("signed total must add the occupation allowance once")
    return {
        "energy_sign": int(sign),
        "angular_sign": int(angular),
        "energy_fiber_hex": _fiber_hex(fiber),
        "source_digest": source_digest,
        "preparation_digest": preparation_digest,
        "archived_distance_lower": archived_lower,
        "archived_distance_upper": archived_upper,
        "signed_total_lower": archived_lower,
        "signed_total_upper": exact_upper(total_upper),
    }


def _prove(spec, model, trace, nfev, dependencies):
    with ctx.workprec(BITS):
        return _prove_at_precision(spec, model, trace, nfev, dependencies)


def _prove_at_precision(spec, model, trace, nfev, dependencies):
    trace = np.ascontiguousarray(trace, dtype=np.float64)
    proof = validate_direct_vacuum(model, trace, START, spec.rho_up, degree=DEGREE)
    initial, bloch = proof["initial_error"], proof["bloch_error"]
    if (not initial.is_finite() or not bloch.is_finite() or not initial > 0
            or not bloch >= initial):
        raise ArithmeticError("nonzero finite direct-vacuum Bloch error required")
    if int(proof["cells"]) != len(trace):
        raise ArithmeticError("Bloch cell count does not match the saved trace")
    compared = compare_signed_archives(
        model, proof, spec.positive_columns, spec.positive_covariance,
        spec.negative_columns, spec.negative_covariance,
        kappa=spec.kappa, omega=spec.omega)
    plain = archived_vacuum_distance(
        spec.negative_columns, spec.negative_covariance, proof, complement=False)
    negative = compared["archived_negative"]
    if (exact_upper(plain["lower"]) == exact_upper(negative["lower"])
            and exact_upper(plain["upper"]) == exact_upper(negative["upper"])):
        raise ValueError("signed complement did not change the archived negative comparison")
    packed_occupation = _pack_occupation(compared["occupation"])
    occupation_upper = packed_occupation["operator_distance_upper"]
    comparisons = [
        _comparison(1, 1, spec.energy_fiber, spec.positive_source_digest,
                    spec.positive_preparation_digest, compared["archived_positive"],
                    occupation_upper),
        _comparison(-1, -1, spec.negative_energy_fiber, spec.negative_source_digest,
                    spec.negative_preparation_digest, negative, occupation_upper),
    ]
    if comparisons[0]["source_digest"] == comparisons[1]["source_digest"]:
        raise ValueError("duplicate source identity")
    end_hex = _hex(model.target_log_distance(spec.rho_up).mid())
    payload = deterministic_npz_bytes({
        "trace": trace,
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
        "energy_sign": 1, "negative_energy_sign": -1,
        "angular_sign": 1, "negative_angular_sign": -1,
        "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
        "negative_energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
        "mass_hex": _hex(spec.mass), "angular_hex": _hex(spec.angular),
        "negative_angular_hex": _hex(spec.negative_angular),
        "mass_enclosure": MASS_ENCLOSURE, "angular_enclosure": ANGULAR_ENCLOSURE,
        "rho_up_hex": _hex(spec.rho_up), "horizon_rho_hex": _hex(spec.horizon_rho),
        "kappa_hex": _hex(spec.kappa), "omega_hex": _hex(spec.omega),
        "weight_hex": _hex(spec.weight), "column_weight_hex": _hex(spec.column_weight),
        "weight_applied_to_error": False,
        "start_log_delta_hex": _hex(START), "end_log_delta_hex": end_hex,
        "bits": BITS, "frame_order": FRAME_ORDER, "metric_terms": METRIC_TERMS,
        "analytic_radius": ANALYTIC_RADIUS, "defect_degree": DEGREE,
        "rtol_hex": _hex(RTOL), "atol_hex": _hex(ATOL), "max_step_hex": _hex(MAX_STEP),
        "method": METHOD, "equation": EQUATION,
        "cells": int(proof["cells"]), "nfev": int(nfev),
        "initial_bloch_error_upper": exact_upper(initial),
        "defect_integral_upper": exact_upper(proof["defect_integral"]),
        "endpoint_bridge_upper": exact_upper(proof["endpoint_bridge"]),
        "vacuum_bloch_error_upper": exact_upper(bloch),
        "validated_vacuum_bloch_endpoint": [
            {"lower": exact_upper(value.lower()), "upper": exact_upper(value.upper())}
            for value in proof["endpoint"]],
        "occupation": packed_occupation,
        "signed_comparisons": comparisons,
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "covariance_errors_unweighted": True,
        "occupation_added_once": True,
        "occupation_included_in_archived_distance": False,
        "negative_comparison_uses_complemented_endpoint": True,
        "positive_scalar_bound_copied_to_negative": False,
        "correction_forcing_used": False,
        "replay_uses_saved_trace": True,
        "archived_source_replaced": False,
        "source_occupations_changed": False,
        "initial_archived_arrays_changed": False,
        "original_preparation_ode_rerun": False,
        "direct_vacuum_trace_saved": True,
        "source_quadrature_error_included": False,
        "all_source_families": False,
        "other_positive_angular_family_included": False,
        "other_energies_filled_from_this_row": False,
        "physical_upstream_budget_component": None,
        "changed_history_C_M": None,
        "physical_local_gate": "OPEN",
        "payload": {"path": f"{relative}/witness.npz", "sha256": sha256(payload).hexdigest(),
                    "bytes": len(payload)},
        "source_hashes": dependencies["source_hashes"],
        "archive_digests": dependencies["archive_digests"],
    }
    return record, payload


def _require(record, field, expected, message):
    if field not in record:
        raise ValueError("missing bound")
    if record[field] != expected:
        raise ValueError(message)


def _require_bindings(record, spec, dependencies):
    if not isinstance(record, dict):
        raise ValueError("corrupt row directory")
    settings = {
        "schema": ROW_SCHEMA, "coverage_schema": SCHEMA, "group": GROUP,
        "angular_sign": 1, "negative_angular_sign": -1,
        "energy_sign": 1, "negative_energy_sign": -1,
        "start_log_delta_hex": _hex(START), "bits": BITS, "frame_order": FRAME_ORDER,
        "metric_terms": METRIC_TERMS, "analytic_radius": ANALYTIC_RADIUS,
        "defect_degree": DEGREE, "rtol_hex": _hex(RTOL), "atol_hex": _hex(ATOL),
        "max_step_hex": _hex(MAX_STEP), "method": METHOD, "equation": EQUATION,
        "mass_enclosure": MASS_ENCLOSURE, "angular_enclosure": ANGULAR_ENCLOSURE,
        "rho_up_hex": _hex(spec.rho_up), "horizon_rho_hex": _hex(spec.horizon_rho),
        "kappa_hex": _hex(spec.kappa), "omega_hex": _hex(spec.omega),
        "correction_forcing_used": False,
        "positive_scalar_bound_copied_to_negative": False,
        "occupation_added_once": True,
        "occupation_included_in_archived_distance": False,
        "weight_applied_to_error": False,
        "covariance_errors_unweighted": True,
    }
    for field, expected in settings.items():
        message = ("direct-source checkpoint schema required"
                   if field in {"schema", "coverage_schema"} else "settings changed")
        _require(record, field, expected, message)
    identity = {
        "panel": spec.panel, "row": spec.row, "local_index": spec.local_index,
        "batch_rows": [spec.batch_rows[0], spec.batch_rows[1]],
        "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
        "negative_energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
        "mass_hex": _hex(spec.mass), "angular_hex": _hex(spec.angular),
        "negative_angular_hex": _hex(spec.negative_angular),
        "weight_hex": _hex(spec.weight), "column_weight_hex": _hex(spec.column_weight),
    }
    for field, expected in identity.items():
        _require(record, field, expected, "changed bound")
    sources = {
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
    }
    for field, expected in sources.items():
        _require(record, field, expected, "original source changed")
    if "source_hashes" not in record or "archive_digests" not in record:
        raise ValueError("missing bound")
    if (record["source_hashes"] != dependencies["source_hashes"]
            or record["archive_digests"] != dependencies["archive_digests"]):
        raise ValueError("changed bound")


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
            files = list(data.files)
            loaded = {name: np.array(data[name], copy=True) for name in files}
    except Exception as error:
        raise ValueError("corrupt row directory") from error
    if set(files) != set(PAYLOAD_NAMES):
        raise ValueError("missing payload")
    originals = {
        "positive_columns": spec.positive_columns,
        "negative_columns": spec.negative_columns,
        "positive_covariance": spec.positive_covariance,
        "negative_covariance": spec.negative_covariance,
    }
    for name, value in originals.items():
        if not np.array_equal(loaded[name], value):
            raise ValueError("original source changed")
    trace = np.ascontiguousarray(loaded["trace"], dtype=np.float64)
    if (trace.ndim != 2 or trace.shape[1] != 26 or len(trace) == 0
            or not np.isfinite(trace).all()):
        raise ValueError("missing payload")
    return trace


def authenticate_checkpoint(spec, dirname, record, payload, dependencies):
    """Replay one saved trace. This never calls the direct-vacuum solver."""
    if dirname != checkpoint_name(spec.panel, spec.row):
        raise ValueError("checkpoint name does not match row identity")
    _require_bindings(record, spec, dependencies)
    trace = _load_payload(record, payload, spec)
    nfev = record.get("nfev")
    if "nfev" not in record:
        raise ValueError("missing bound")
    if isinstance(nfev, bool) or type(nfev) is not int or nfev < 0:
        raise ValueError("settings changed")
    with ctx.workprec(BITS):
        model = _model(spec)
        end_hex = _hex(model.target_log_distance(spec.rho_up).mid())
        _require(record, "end_log_delta_hex", end_hex, "settings changed")
        rebuilt, rebuilt_payload = _prove(spec, model, trace, nfev, dependencies)
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


def _refuse_foreign_campaign(output_root):
    parts = Path(output_root).resolve(strict=False).parts
    marker = ("results", "development", "nsc-middle-source-coverage-v1")
    if any(parts[index:index+3] == marker for index in range(len(parts)-2)):
        raise ValueError("refusing to write the frozen middle-source campaign")


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
    _refuse_foreign_campaign(output_root)
    root = _canonical_directory(output_root, create=True)
    (root/"rows").mkdir(exist_ok=True)
    publish_exclusive_directory(root, f"rows/{name}", {
        "record.json": encode_record(record), "witness.npz": payload,
    })


def _read_checkpoint(path):
    names = {item.name for item in path.iterdir()}
    if names != {"record.json", "witness.npz"}:
        if not {"record.json", "witness.npz"} <= names:
            raise ValueError("missing payload")
        raise ValueError("corrupt row directory")
    raw = (path/"record.json").read_bytes()
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError("corrupt row directory") from error
    if not isinstance(record, dict) or encode_record(record) != raw:
        raise ValueError("corrupt row directory")
    return record, (path/"witness.npz").read_bytes()


def _row_key(record):
    if not isinstance(record, dict):
        raise ValueError("corrupt row directory")
    panel, row, fiber = record.get("panel"), record.get("row"), record.get("energy_fiber_hex")
    if "panel" not in record or "row" not in record or "energy_fiber_hex" not in record:
        raise ValueError("missing bound")
    if (not isinstance(panel, str) or isinstance(row, bool) or type(row) is not int
            or not isinstance(fiber, list) or len(fiber) != 3
            or not all(isinstance(item, str) for item in fiber)):
        raise ValueError("changed bound")
    return panel, row, fiber[0]


def load_checkpoints(output_root):
    root = _canonical_directory(output_root, create=False)
    rows = root/"rows"
    if not rows.exists():
        return []
    if not rows.is_dir():
        raise ValueError("corrupt row directory")
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


def _summary(record):
    return {
        "panel": record["panel"], "row": record["row"],
        "energy_fiber_hex": record["energy_fiber_hex"],
        "negative_energy_fiber_hex": record["negative_energy_fiber_hex"],
        "checkpoint": record["payload"]["path"].rsplit("/", 1)[0],
        "cells": record["cells"], "nfev": record["nfev"],
        "equation": record["equation"],
        "vacuum_bloch_error_upper": record["vacuum_bloch_error_upper"],
        "occupation": record["occupation"],
        "signed_comparisons": record["signed_comparisons"],
        "weight_hex": record["weight_hex"],
        "column_weight_hex": record["column_weight_hex"],
        "weight_applied_to_error": False,
        "correction_forcing_used": False,
        "physical_upstream_budget_component": None,
        "physical_local_gate": "OPEN",
    }


def _comparison_upper(record):
    comparisons = record.get("signed_comparisons")
    if not isinstance(comparisons, list) or len(comparisons) != 2:
        raise ValueError("missing bound")
    if [item.get("energy_sign") for item in comparisons] != [1, -1]:
        raise ValueError("missing bound")
    uppers = []
    for item in comparisons:
        if not isinstance(item, dict) or "signed_total_upper" not in item:
            raise ValueError("missing bound")
        uppers.append(restored_upper(item["signed_total_upper"]))
    return uppers


def aggregate(catalogue, completed, *, mode, row_budget, cpu_budget, cpu_seconds,
              capture_cpu_seconds, new_rows_captured):
    by_key = {(record["panel"], record["row"]): record for record in completed}
    if len(by_key) != len(completed):
        raise ValueError("duplicate source identity")
    catalogue_keys = {(spec.panel, spec.row) for spec in catalogue}
    if any(key not in catalogue_keys for key in by_key):
        raise ValueError("checkpoint identity is not an original direct-source row")
    ordered = [by_key[(spec.panel, spec.row)] for spec in catalogue
               if (spec.panel, spec.row) in by_key]
    missing = [{"panel": spec.panel, "row": spec.row, "energy_hex": _hex(spec.energy)}
               for spec in catalogue if (spec.panel, spec.row) not in by_key]
    worst = None
    with ctx.workprec(BITS):
        for record in ordered:
            for upper in _comparison_upper(record):
                worst = upper if worst is None else max(worst, upper)
    complete = len(ordered) == len(catalogue) and not missing
    return {
        "schema": SCHEMA,
        "status": WINDOW_STATUS if complete else OPEN_STATUS,
        "group": GROUP, "angular_sign": ANGULAR_SIGN,
        "energy_window": {"lower": ENERGY_LOWER, "lower_included": True,
                          "upper": ENERGY_UPPER, "upper_included": False},
        "panels": [LOW16_PANEL, LOW32_PANEL],
        "expected_positive_rows": EXPECTED_POSITIVE_ROWS,
        "expected_negative_partners": EXPECTED_POSITIVE_ROWS,
        "completed_positive_rows": len(ordered),
        "completed_negative_partners": len(ordered),
        "completed_signed_rows": 2*len(ordered),
        "missing_positive_rows": len(missing),
        "coverage_complete": complete,
        "interrupted": not complete,
        "finite_window_closes_gate": False,
        "mode": mode,
        "rows": [_summary(record) for record in ordered],
        "missing": missing,
        "maximum_completed_signed_total_upper": None if worst is None else exact_upper(worst),
        "covariance_errors_unweighted": True,
        "weight_applied_to_error": False,
        "occupation_added_once": True,
        "correction_forcing_used": False,
        "source_quadrature_error_included": False,
        "all_source_families": False,
        "other_positive_angular_family_included": False,
        "archived_source_replaced": False,
        "source_occupations_changed": False,
        "initial_archived_arrays_changed": False,
        "other_energies_filled_from_completed_rows": False,
        "physical_upstream_budget_component": None,
        "changed_history_C_M": None,
        "physical_local_gate": "OPEN",
        "row_budget": row_budget, "cpu_budget_seconds": cpu_budget,
        "new_rows_captured": new_rows_captured,
        "cpu_seconds": cpu_seconds, "capture_cpu_seconds": capture_cpu_seconds,
        "replay_uses_saved_trace": True,
    }


def _authenticate_saved(output_root, catalogue, dependencies):
    by_key = {(spec.panel, spec.row): spec for spec in catalogue}
    authenticated = []
    for dirname, record, payload in load_checkpoints(output_root):
        panel, row, _energy = _row_key(record)
        spec = by_key.get((panel, row))
        if spec is None:
            raise ValueError("checkpoint identity is not an original direct-source row")
        authenticated.append(authenticate_checkpoint(
            spec, dirname, record, payload, dependencies))
    return authenticated


def _resolve_output(archive_root, output_root):
    output_root = Path(output_root)
    if not output_root.is_absolute():
        output_root = Path(archive_root)/output_root
    return output_root.resolve(strict=False)


def cover(archive_root, output_root, *, mode, row_budget=None, cpu_budget=None):
    """Validate saved rows before any new capture. Check mode never solves."""
    if mode not in {"resume", "check"}:
        raise ValueError("coverage mode must be resume or check")
    _require_owner_tolerances()
    if mode == "check":
        if row_budget is not None or cpu_budget is not None:
            raise ValueError("check replays saved witnesses and does not accept a solve budget")
    else:
        row_budget, cpu_budget = _budgets(row_budget, cpu_budget)
    started = time.process_time()
    archive_root = Path(archive_root).resolve(strict=True)
    output_root = _resolve_output(archive_root, output_root)
    _refuse_foreign_campaign(output_root)
    archive = RetainedUpstreamArchive(archive_root)
    catalogue = direct_catalogue(archive)
    dependencies = dependency_identity(archive_root, archive)
    if mode == "check" and not output_root.exists():
        return aggregate(catalogue, [], mode=mode, row_budget=None, cpu_budget=None,
                         cpu_seconds=time.process_time()-started, capture_cpu_seconds=0.0,
                         new_rows_captured=0)
    output_root = _canonical_directory(output_root, create=mode == "resume")
    # Every saved row is proved before a budget is allowed to start a solve.
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
        raise TimeoutError("direct vacuum Bloch CPU budget exhausted")
    with ctx.workprec(BITS):
        model = _model(spec)
        trace, nfev = capture_direct_vacuum(
            model, START, spec.rho_up, max_step=MAX_STEP, cpu_limit=float(cpu_limit))
        return _prove(spec, model, trace, int(nfev), dependencies)
