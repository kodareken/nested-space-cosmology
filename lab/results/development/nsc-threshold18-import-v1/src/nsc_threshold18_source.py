"""Actual18 direct-vacuum recorder for group14/low32_1 rows 32 through 40.

The homogeneous Bloch transport, horizon frame, and signed distance are the
reviewed direct-vacuum owner. This module does not solve a new source law
and it does not extend the frozen 2<=E<16 window. Occupation is added once.
Missing rows stay identities with a null aggregate, never a fabricated zero.
"""
from hashlib import sha256
from io import BytesIO
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import zipfile

import numpy as np
from flint import arb, ctx

_SUCCESSOR_ROOT = Path(__file__).resolve().parents[1]
_PIN_ROOT = _SUCCESSOR_ROOT / "pin"
_LAB = _PIN_ROOT / "lab"
_PIN_SRC = _LAB / "src"
_PARENT_LAB = _SUCCESSOR_ROOT.parent / "lab"
if str(_PIN_SRC) not in sys.path:
    sys.path.insert(0, str(_PIN_SRC))

from recursive_horizons.evidence_io import publish_exclusive_directory
from recursive_horizons.nsc_direct_vacuum_source import (
    DirectVacuumBloch, archived_vacuum_distance, capture_direct_vacuum,
    compare_signed_archives, validate_direct_vacuum,
)
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import file_digest, implementation_hashes
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_signed_state import (
    SOURCE_COMPLEMENT_TOLERANCE, s3_conjugate, source_complement_residual,
)
import recursive_horizons.nsc_direct_vacuum_source as _vacuum_owner
import recursive_horizons.nsc_source_occupation_enclosure as _occupation_owner


SCHEMA = "NSC-THRESHOLD18-ACTUAL-SOURCE-v1"
ROW_SCHEMA = "NSC-THRESHOLD18-ACTUAL-ROW-v1"
BASE_COMMIT = "eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9"
BASE_REPOSITORY = "https://github.com/kodareken/nested-space-cosmology"
GROUP = 14
ANGULAR_SIGN = 1
PANEL = "group14/low32_1"
ROWS = tuple(range(32, 41))
MIN_ENERGY_HEX = "0x1.92f9814c33825p+0"
MAX_ENERGY_HEX = "0x1.fae56f00f6a85p+0"
EXPECTED_POSITIVE_ROWS = 9
PUBLISHED_COVERAGE = "results/development/nsc-ks-source-operator-majorant-v3.json"
PUBLISHED_SIGNED_ROWS = 868
START = -18.0
BITS = 192
FRAME_ORDER = 16
METRIC_TERMS = 48
DEGREE = 12
MAX_STEP = 0.0125
RTOL = 2e-13
ATOL = 2e-15
ANALYTIC_RADIUS = "0.1"
METHOD = "DOP853"
EQUATION = "n_y=2*cross(h(y),n)"
MASS_ENCLOSURE = "archived-float union pi/2"
ANGULAR_ENCLOSURE = "archived-float union sqrt(5)"
ADJUDICATION = "c82b600b"
ROW_STATUS = (
    "OPEN: one original direct-vacuum row and its negative partner "
    "in pi/2<=E<2; finite actual18 window only; local gate open"
)
OPEN_STATUS = (
    "OPEN: group14 angular+1 actual18 window pi/2<=E<2 is incomplete; "
    "completed rows only; no missing row is entered as zero; "
    "full upstream budget and local gate remain open"
)
WINDOW_STATUS = (
    "OPEN: finite group14 angular+1 actual18 window pi/2<=E<2 and the "
    "signed partners are recorded; full upstream budget and local gate remain open"
)
PAYLOAD_NAMES = (
    "negative_columns", "negative_covariance", "positive_columns",
    "positive_covariance", "trace",
)
METHOD_OWNERS = (
    "src/recursive_horizons/evidence_io.py",
    "src/recursive_horizons/nsc_direct_vacuum_source.py",
    "src/recursive_horizons/nsc_ks_ball_trajectory.py",
    "src/recursive_horizons/nsc_ks_evaluation_binding.py",
    "src/recursive_horizons/nsc_ks_retained_upstream_archive.py",
    "src/recursive_horizons/nsc_ks_signed_state.py",
    "src/recursive_horizons/nsc_metric_horizon_frame.py",
    "src/recursive_horizons/nsc_source_occupation_enclosure.py",
    "src/recursive_horizons/nsc_subgap_source_covariance.py",
    "src/recursive_horizons/nsc_vacuum_source_remainder.py",
)
SUCCESSOR_FILES = (
    "docs/nsc-threshold18-actual-source.md",
    "scripts/derive_nsc_threshold18_source.py",
    "src/nsc_threshold18_source.py",
)
REVIEWED_OWNERS = (
    "src/recursive_horizons/nsc_direct_vacuum_source.py",
    "src/recursive_horizons/nsc_source_occupation_enclosure.py",
)
FROZEN_CAMPAIGNS = (
    ("results", "development", "nsc-middle-source-coverage-v1"),
    ("results", "development", "nsc-direct-source-window-v1"),
)
_CAPTURE_CONSTANTS = capture_direct_vacuum.__code__.co_consts


def encode_record(record):
    return (json.dumps(record, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def deterministic_npz_bytes(arrays):
    """Byte-stable uncompressed npz, same layout as the pinned Cauchy writer."""
    output = BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            value = np.ascontiguousarray(arrays[name])
            if value.dtype.hasobject:
                raise TypeError("object arrays are forbidden")
            stream = BytesIO()
            np.lib.format.write_array(stream, value, allow_pickle=False)
            info = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info, stream.getvalue())
    return output.getvalue()


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


def in_actual18(energy):
    value = float(energy)
    if not math.isfinite(value):
        raise ValueError("finite actual18 energy required")
    return math.pi / 2 <= value < 2.0


def in_fixed_direct_window(energy):
    """The frozen direct-window predicate, 2<=E<16. It does not own these rows."""
    value = float(energy)
    if not math.isfinite(value):
        raise ValueError("finite direct-window energy required")
    return 2.0 <= value < 16.0


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
    owner = Path(_vacuum_owner.__file__).resolve()
    if not owner.is_relative_to(_PIN_SRC.resolve()):
        raise ValueError("direct vacuum owner is not the pinned source")
    occupation = Path(_occupation_owner.__file__).resolve()
    if not occupation.is_relative_to(_PIN_SRC.resolve()):
        raise ValueError("occupation owner is not the pinned source")


def _require_clean_base(lab):
    pin = Path(lab).resolve().parent
    if pin.resolve() != _PIN_ROOT.resolve():
        raise ValueError("authoritative base must be the pinned checkout")
    head = subprocess.check_output(
        ["git", "-C", str(pin), "rev-parse", "HEAD"], text=True).strip()
    if head != BASE_COMMIT:
        raise ValueError("authoritative base changed")
    dirty = subprocess.check_output(
        ["git", "-C", str(pin), "status", "--porcelain"], text=True)
    if dirty.strip():
        raise ValueError("authoritative base is dirty")
    return head


def _reviewed_snapshot_match():
    matched = {}
    for relative in REVIEWED_OWNERS:
        pinned = _LAB / relative
        parent = _PARENT_LAB / relative
        pinned_bytes = pinned.read_bytes()
        parent_bytes = parent.read_bytes()
        if pinned_bytes != parent_bytes:
            raise ValueError("reviewed snapshot differs from the pinned direct-vacuum owner")
        matched[relative] = sha256(pinned_bytes).hexdigest()
    return matched


class ActualRow:
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
        raise ValueError("actual18 must stay on original group14 angular+1")
    if partner.angular_sign != -1 or partner.energy_sign != -1:
        raise ValueError("actual opposite-angular negative partner required")
    if partner.original_panel != positive.original_panel or partner.rows != positive.rows:
        raise ValueError("negative partner must share the original panel rows")
    if positive.mass != np.pi / 2 or positive.angular != np.sqrt(5):
        raise ValueError("original group14 mass and angular label required")
    if (partner.mass != positive.mass or partner.angular != -positive.angular
            or partner.rho_up != positive.rho_up):
        raise ValueError("signed mass, angular label or rho changed")
    row_index = int(positive.rows[0] + local)
    sl = slice(3 * local, 3 * local + 3)
    fiber = np.asarray(positive.source.energies[sl], float)
    negative_fiber = np.asarray(partner.source.energies[sl], float)
    if fiber.shape != (3,) or not np.all(fiber == fiber[0]) or not in_actual18(float(fiber[0])):
        raise ValueError("actual18 three-column energy fiber required")
    if in_fixed_direct_window(float(fiber[0])):
        raise ValueError("fixed direct window cannot claim an actual18 row")
    if not np.array_equal(negative_fiber, -fiber):
        raise ValueError("negative partner energies must equal -E")
    positive_columns = _frozen(positive.initial_columns[:, sl], np.complex128)
    negative_columns = _frozen(partner.initial_columns[:, sl], np.complex128)
    positive_covariance = _frozen(positive.source.covariance[sl, sl], np.complex128)
    negative_covariance = _frozen(partner.source.covariance[sl, sl], np.complex128)
    if not np.array_equal(negative_columns, s3_conjugate(positive_columns)):
        raise ValueError("negative columns are not the S3 conjugate of the positive columns")
    residual = source_complement_residual(positive_covariance, negative_covariance)
    if not math.isfinite(residual) or residual > SOURCE_COMPLEMENT_TOLERANCE:
        raise ValueError("negative source differs from the fixed horizon source law")
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
    return ActualRow(
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
        horizon_rho=float(config["horizon_rho"]),
        source_complement_residual=float(residual))


def _published_coverage(lab):
    path = Path(lab) / PUBLISHED_COVERAGE
    payload = json.loads(path.read_text())
    rows = payload.get("source_coverage", {}).get("rows")
    covered_count = payload.get("source_coverage", {}).get("covered_signed_rows")
    if not isinstance(rows, list) or covered_count != PUBLISHED_SIGNED_ROWS or len(rows) != PUBLISHED_SIGNED_ROWS:
        raise ValueError("published source coverage census changed")
    covered = set()
    for item in rows:
        if not isinstance(item, dict):
            raise ValueError("published source coverage census changed")
        key = (item.get("panel"), item.get("row"), item.get("energy_sign"))
        if key in covered or not isinstance(key[0], str) or type(key[1]) is not int or key[2] not in (1, -1):
            raise ValueError("published source coverage census changed")
        covered.add(key)
    return covered, file_digest(path), payload.get("input_hashes")


def _require_census(rows, covered):
    if len(rows) != EXPECTED_POSITIVE_ROWS:
        raise ValueError("original group14 actual18 row count changed")
    if tuple(row.row for row in rows) != ROWS:
        raise ValueError("group14/low32_1 actual18 rows changed")
    if {row.panel for row in rows} != {PANEL}:
        raise ValueError("actual18 panels changed")
    energies = [row.energy for row in rows]
    if energies != sorted(energies) or len(set(energies)) != len(energies):
        raise ValueError("duplicate source identity")
    if _hex(rows[0].energy) != MIN_ENERGY_HEX or _hex(rows[-1].energy) != MAX_ENERGY_HEX:
        raise ValueError("actual18 energy endpoints changed")
    if any(not in_actual18(row.energy) or in_fixed_direct_window(row.energy) for row in rows):
        raise ValueError("actual18 energy endpoints changed")
    signed = {(row.panel, row.row, sign) for row in rows for sign in (1, -1)}
    if len(signed) != 2 * EXPECTED_POSITIVE_ROWS or signed & covered:
        raise ValueError("actual18 row is already in the published source coverage")


def actual_catalogue(archive, lab):
    """Original group14 angular+1 rows with pi/2<=E<2 and their negative partners."""
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
            if in_actual18(float(fiber[0])):
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
    covered, _digest, _inputs = _published_coverage(lab)
    _require_census(tuple(rows), covered)
    return tuple(rows)


def _successor_hashes():
    hashes = {}
    for relative in SUCCESSOR_FILES:
        path = _SUCCESSOR_ROOT / relative
        if not path.is_file():
            raise ValueError("missing successor source " + relative)
        hashes[relative] = sha256(path.read_bytes()).hexdigest()
    return hashes


def _method_hashes(lab):
    hashes = implementation_hashes(Path(lab), owners=METHOD_OWNERS)
    if not hashes or any(name not in hashes for name in METHOD_OWNERS):
        raise ValueError("source and dependency identity required")
    return hashes


def dependency_identity(lab, archive):
    _covered, published_digest, published_inputs = _published_coverage(lab)
    if not isinstance(published_inputs, dict):
        raise ValueError("published source coverage census changed")
    shared = set(archive.input_hashes) & set(published_inputs)
    if not shared:
        raise ValueError("baseline input hash changed")
    for name in shared:
        if archive.input_hashes[name] != published_inputs[name]:
            raise ValueError("baseline input hash changed")
    family = archive._families.get((GROUP, ANGULAR_SIGN))
    if not isinstance(family, dict) or not isinstance(family.get("source_hashes"), dict):
        raise ValueError("baseline family source hashes required")
    if not family["source_hashes"]:
        raise ValueError("baseline family source hashes required")
    source_hashes = _method_hashes(lab)
    successor = _successor_hashes()
    inputs = dict(archive.input_hashes)
    if not source_hashes or not inputs or not successor:
        raise ValueError("source and dependency identity required")
    return {
        "source_hashes": source_hashes,
        "successor_hashes": successor,
        "archive_digests": inputs,
        "baseline_family_source_hashes": dict(family["source_hashes"]),
        "published_coverage_digest": published_digest,
        "published_signed_rows": PUBLISHED_SIGNED_ROWS,
        "reviewed_snapshot_sha256": _reviewed_snapshot_match(),
        "base_commit": _require_clean_base(lab),
        "base_repository": BASE_REPOSITORY,
        "covered_keys": _covered,
    }


def _model(spec):
    mass = arb(spec.mass).union(arb.pi() / 2)
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
    if allowance is None or not allowance.is_finite() or not allowance > 0:
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
    if archived.get("lower") is None or archived.get("upper") is None:
        raise ValueError("missing bound")
    if not archived["lower"].is_finite() or not archived["upper"].is_finite():
        raise ArithmeticError("finite archived distance required")
    if not archived["upper"] >= archived["lower"]:
        raise ArithmeticError("archived distance bounds are reversed")
    archived_lower = exact_upper(archived["lower"])
    archived_upper = exact_upper(archived["upper"])
    total_upper = (restored_upper(archived_upper) + restored_upper(occupation_upper)).upper()
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
        "weight_applied_to_error": False,
    }


def _prove(spec, model, trace, nfev, dependencies):
    with ctx.workprec(BITS):
        return _prove_at_precision(spec, model, trace, nfev, dependencies)


def _prove_at_precision(spec, model, trace, nfev, dependencies):
    trace = np.ascontiguousarray(trace, dtype=np.float64)
    proof = validate_direct_vacuum(model, trace, START, spec.rho_up, degree=DEGREE)
    initial, bloch = proof["initial_error"], proof["bloch_error"]
    if (initial is None or bloch is None or not initial.is_finite() or not bloch.is_finite()
            or not initial > 0 or not bloch >= initial):
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
    endpoint = proof["endpoint"]
    manual = dict(proof)
    manual["endpoint"] = (endpoint[0], -endpoint[1], -endpoint[2])
    manual_distance = archived_vacuum_distance(
        spec.negative_columns, spec.negative_covariance, manual, complement=False)
    if (exact_upper(manual_distance["lower"]) != exact_upper(negative["lower"])
            or exact_upper(manual_distance["upper"]) != exact_upper(negative["upper"])):
        raise ValueError("negative endpoint is not (nx,-ny,-nz)")
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
    if comparisons[0]["archived_distance_upper"] == comparisons[1]["archived_distance_upper"]:
        raise ValueError("positive scalar bound was copied onto the negative partner")
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
        "source_complement_residual_hex": _hex(spec.source_complement_residual),
        "negative_covariance_law": "Cminus=I-conj(Cplus)",
        "negative_column_law": "S3 column conjugation",
        "negative_density_endpoint": "(nx,-ny,-nz)",
        "negative_endpoint_adjudication": ADJUDICATION,
        "negative_endpoint_adjudication_flipped": False,
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
        "fixed_direct_window_claims_this_row": False,
        "other_energies_filled_from_this_row": False,
        "subgap108_batch_included": False,
        "physical_upstream_budget_component": None,
        "changed_history_C_M": None,
        "physical_local_gate": "OPEN",
        "payload": {"path": f"{relative}/witness.npz", "sha256": sha256(payload).hexdigest(),
                    "bytes": len(payload)},
        "source_hashes": dependencies["source_hashes"],
        "successor_hashes": dependencies["successor_hashes"],
        "archive_digests": dependencies["archive_digests"],
        "baseline_family_source_hashes": dependencies["baseline_family_source_hashes"],
        "published_coverage_digest": dependencies["published_coverage_digest"],
        "reviewed_snapshot_sha256": dependencies["reviewed_snapshot_sha256"],
        "base_commit": dependencies["base_commit"],
        "base_repository": dependencies["base_repository"],
    }
    return record, payload


def _require(record, field, expected, message):
    if field not in record or record[field] is None and expected is not None:
        raise ValueError("missing bound")
    if record[field] != expected:
        raise ValueError(message)


def _binding_expectations(spec, dependencies):
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
        "negative_covariance_law": "Cminus=I-conj(Cplus)",
        "negative_column_law": "S3 column conjugation",
        "negative_density_endpoint": "(nx,-ny,-nz)",
        "negative_endpoint_adjudication": ADJUDICATION,
        "negative_endpoint_adjudication_flipped": False,
        "fixed_direct_window_claims_this_row": False,
        "subgap108_batch_included": False,
        "base_commit": BASE_COMMIT,
        "base_repository": BASE_REPOSITORY,
        "panel": spec.panel, "row": spec.row, "local_index": spec.local_index,
        "batch_rows": [spec.batch_rows[0], spec.batch_rows[1]],
        "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
        "negative_energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
        "mass_hex": _hex(spec.mass), "angular_hex": _hex(spec.angular),
        "negative_angular_hex": _hex(spec.negative_angular),
        "weight_hex": _hex(spec.weight), "column_weight_hex": _hex(spec.column_weight),
        "source_complement_residual_hex": _hex(spec.source_complement_residual),
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "source_hashes": dependencies["source_hashes"],
        "successor_hashes": dependencies["successor_hashes"],
        "archive_digests": dependencies["archive_digests"],
        "baseline_family_source_hashes": dependencies["baseline_family_source_hashes"],
        "published_coverage_digest": dependencies["published_coverage_digest"],
        "reviewed_snapshot_sha256": dependencies["reviewed_snapshot_sha256"],
    }
    return settings


def binding_record(spec, dependencies):
    """Identity fields for refusal tests. This is not a solved witness."""
    return _binding_expectations(spec, dependencies)


def _require_bindings(record, spec, dependencies):
    if not isinstance(record, dict):
        raise ValueError("corrupt row directory")
    expected = _binding_expectations(spec, dependencies)
    for field, value in expected.items():
        if field in {"schema", "coverage_schema"}:
            message = "actual18 checkpoint schema required"
        elif field in {"source_hashes", "successor_hashes", "archive_digests",
                       "baseline_family_source_hashes", "published_coverage_digest",
                       "reviewed_snapshot_sha256", "base_commit", "base_repository"}:
            message = "changed bound"
        elif field.endswith("digest") or field.endswith("sha256") or "source" in field:
            message = "original source changed" if "digest" in field and "published" not in field else "changed bound"
        elif field in {"panel", "row", "local_index", "batch_rows", "energy_fiber_hex",
                       "negative_energy_fiber_hex", "mass_hex", "angular_hex",
                       "negative_angular_hex", "weight_hex", "column_weight_hex",
                       "source_complement_residual_hex"}:
            message = "changed bound"
        else:
            message = "settings changed"
        if field in {"positive_source_digest", "negative_source_digest",
                     "positive_preparation_digest", "negative_preparation_digest"}:
            message = "original source changed"
        _require(record, field, value, message)


def _load_payload(record, payload, spec):
    described = record.get("payload")
    if not isinstance(described, dict):
        raise ValueError("missing payload")
    if described.get("sha256") is None or described.get("bytes") is None or described.get("path") is None:
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
        if name not in loaded or loaded[name] is None:
            raise ValueError("missing payload")
        if not np.array_equal(loaded[name], value):
            raise ValueError("original source changed")
    if "trace" not in loaded or loaded["trace"] is None:
        raise ValueError("missing payload")
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
    if "nfev" not in record or record.get("nfev") is None:
        raise ValueError("missing bound")
    nfev = record.get("nfev")
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


def _refuse_foreign_output(output_root):
    path = Path(output_root).resolve(strict=False)
    blocked = (
        _PIN_ROOT.resolve(),
        _PARENT_LAB.resolve(),
        (_SUCCESSOR_ROOT / ".venv").resolve(),
    )
    for root in blocked:
        if path == root or root in path.parents:
            raise ValueError("refusing to write pinned source, the parent snapshot, or the validation environment")
    parts = path.parts
    for marker in FROZEN_CAMPAIGNS:
        if any(parts[index:index + len(marker)] == marker for index in range(len(parts) - len(marker) + 1)):
            raise ValueError("refusing to write a frozen source campaign")


def publish_checkpoint(output_root, record, payload):
    """Atomically publish one new row directory. An existing row is left untouched."""
    name = checkpoint_name(record["panel"], record["row"])
    described = record.get("payload") if isinstance(record, dict) else None
    if not isinstance(described, dict):
        raise ValueError("missing payload")
    if described.get("sha256") is None or described.get("bytes") is None:
        raise ValueError("missing payload")
    digest = sha256(payload).hexdigest()
    if described.get("sha256") != digest or described.get("bytes") != len(payload):
        raise ValueError("payload hash changed")
    if described.get("path") != f"rows/{name}/witness.npz":
        raise ValueError("payload hash changed")
    _refuse_foreign_output(output_root)
    root = _canonical_directory(output_root, create=True)
    (root / "rows").mkdir(exist_ok=True)
    publish_exclusive_directory(root, f"rows/{name}", {
        "record.json": encode_record(record), "witness.npz": payload,
    })


def _read_checkpoint(path):
    names = {item.name for item in path.iterdir()}
    if names != {"record.json", "witness.npz"}:
        if not {"record.json", "witness.npz"} <= names:
            raise ValueError("missing payload")
        raise ValueError("corrupt row directory")
    raw = (path / "record.json").read_bytes()
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError("corrupt row directory") from error
    if not isinstance(record, dict) or encode_record(record) != raw:
        raise ValueError("corrupt row directory")
    return record, (path / "witness.npz").read_bytes()


def _row_key(record):
    if not isinstance(record, dict):
        raise ValueError("corrupt row directory")
    panel, row, fiber = record.get("panel"), record.get("row"), record.get("energy_fiber_hex")
    if panel is None or row is None or fiber is None:
        raise ValueError("missing bound")
    if (not isinstance(panel, str) or isinstance(row, bool) or type(row) is not int
            or not isinstance(fiber, list) or len(fiber) != 3
            or not all(isinstance(item, str) for item in fiber)):
        raise ValueError("changed bound")
    return panel, row, fiber[0]


def load_checkpoints(output_root):
    root = _canonical_directory(output_root, create=False)
    rows = root / "rows"
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
        "negative_density_endpoint": record["negative_density_endpoint"],
        "negative_endpoint_adjudication_flipped": False,
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
        if not isinstance(item, dict) or item.get("signed_total_upper") is None:
            raise ValueError("missing bound")
        if item.get("weight_applied_to_error") is not False:
            raise ValueError("weight was applied to an unweighted signed error")
        uppers.append(restored_upper(item["signed_total_upper"]))
    return uppers


def _closure(dependencies):
    return {
        "base_commit": dependencies["base_commit"],
        "base_repository": dependencies["base_repository"],
        "reviewed_snapshot_sha256": dependencies["reviewed_snapshot_sha256"],
        "published_coverage_digest": dependencies["published_coverage_digest"],
        "published_signed_rows": PUBLISHED_SIGNED_ROWS,
        "source_hashes": dependencies["source_hashes"],
        "successor_hashes": dependencies["successor_hashes"],
        "archive_digests": dependencies["archive_digests"],
        "baseline_family_source_hashes": dependencies["baseline_family_source_hashes"],
        "myrsa_used": False,
        "stale_nuc_checkout_used": False,
        "authoritative_transport": "github-pin",
    }


def aggregate(catalogue, completed, dependencies, *, mode, row_budget, cpu_budget, cpu_seconds,
              capture_cpu_seconds, new_rows_captured):
    by_key = {(record["panel"], record["row"]): record for record in completed}
    if len(by_key) != len(completed):
        raise ValueError("duplicate source identity")
    catalogue_keys = {(spec.panel, spec.row) for spec in catalogue}
    if any(key not in catalogue_keys for key in by_key):
        raise ValueError("checkpoint identity is not an original actual18 row")
    ordered = [by_key[(spec.panel, spec.row)] for spec in catalogue
               if (spec.panel, spec.row) in by_key]
    missing = [{"panel": spec.panel, "row": spec.row, "energy_hex": _hex(spec.energy),
                "signed_total_upper": None, "archived_distance_upper": None}
               for spec in catalogue if (spec.panel, spec.row) not in by_key]
    if any(item["signed_total_upper"] is not None or item["archived_distance_upper"] is not None
           for item in missing):
        raise ValueError("missing actual18 row was fabricated")
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
        "panel": PANEL,
        "rows_included": list(ROWS),
        "energy_window": {
            "lower": "pi/2", "lower_included": True,
            "lower_pi_over_2_hex": _hex(math.pi / 2),
            "upper": 2.0, "upper_included": False,
            "minimum_row_energy_hex": MIN_ENERGY_HEX,
            "maximum_row_energy_hex": MAX_ENERGY_HEX,
        },
        "fixed_direct_window": {"lower": 2.0, "lower_included": True,
                                "upper": 16.0, "upper_included": False,
                                "claims_these_rows": False},
        "expected_positive_rows": EXPECTED_POSITIVE_ROWS,
        "expected_negative_partners": EXPECTED_POSITIVE_ROWS,
        "expected_signed_rows": 2 * EXPECTED_POSITIVE_ROWS,
        "completed_positive_rows": len(ordered),
        "completed_negative_partners": len(ordered),
        "completed_signed_rows": 2 * len(ordered),
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
        "negative_density_endpoint": "(nx,-ny,-nz)",
        "negative_covariance_law": "Cminus=I-conj(Cplus)",
        "negative_column_law": "S3 column conjugation",
        "negative_endpoint_adjudication": ADJUDICATION,
        "negative_endpoint_adjudication_flipped": False,
        "correction_forcing_used": False,
        "source_quadrature_error_included": False,
        "all_source_families": False,
        "fixed_direct_window_claims_these_rows": False,
        "subgap108_batch_included": False,
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
        "closure": _closure(dependencies),
        "mac_integration": {
            "lab_root": "/Users/admin/Documents/BlackHoles-Infinity/lab",
            "recorder": "src/recursive_horizons/nsc_threshold18_source.py",
            "script": "scripts/derive_nsc_threshold18_source.py",
            "tests": "tests/test_nsc_threshold18_source.py",
            "doc": "docs/nsc-threshold18-actual-source.md",
            "checkpoints": "results/development/nsc-threshold18-actual18-v1",
            "primary_edited": False,
        },
    }


def _authenticate_saved(output_root, catalogue, dependencies):
    by_key = {(spec.panel, spec.row): spec for spec in catalogue}
    authenticated = []
    for dirname, record, payload in load_checkpoints(output_root):
        panel, row, _energy = _row_key(record)
        spec = by_key.get((panel, row))
        if spec is None:
            raise ValueError("checkpoint identity is not an original actual18 row")
        authenticated.append(authenticate_checkpoint(
            spec, dirname, record, payload, dependencies))
    return authenticated


def _resolve_output(output_root):
    output_root = Path(output_root)
    if not output_root.is_absolute():
        output_root = _SUCCESSOR_ROOT / output_root
    return output_root.resolve(strict=False)


def load_context(lab_root=None):
    _require_owner_tolerances()
    lab = Path(lab_root) if lab_root is not None else _LAB
    lab = lab.resolve(strict=True)
    _require_clean_base(lab)
    archive = RetainedUpstreamArchive(lab)
    catalogue = actual_catalogue(archive, lab)
    dependencies = dependency_identity(lab, archive)
    return lab, catalogue, dependencies


def cover(output_root, *, mode, row_budget=None, cpu_budget=None, lab_root=None):
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
    output_root = _resolve_output(output_root)
    _refuse_foreign_output(output_root)
    _lab, catalogue, dependencies = load_context(lab_root)
    if mode == "check" and not output_root.exists():
        return aggregate(catalogue, [], dependencies, mode=mode, row_budget=None, cpu_budget=None,
                         cpu_seconds=time.process_time() - started, capture_cpu_seconds=0.0,
                         new_rows_captured=0)
    output_root = _canonical_directory(output_root, create=mode == "resume")
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
            remaining = cpu_budget - (time.process_time() - capture_started)
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
    return aggregate(
        catalogue, authenticated, dependencies, mode=mode, row_budget=row_budget,
        cpu_budget=cpu_budget, cpu_seconds=time.process_time() - started,
        capture_cpu_seconds=time.process_time() - capture_started if mode == "resume" else 0.0,
        new_rows_captured=captured)


def _capture(spec, dependencies, *, cpu_limit):
    if not math.isfinite(cpu_limit) or cpu_limit <= 0:
        raise TimeoutError("direct vacuum Bloch CPU budget exhausted")
    with ctx.workprec(BITS):
        model = _model(spec)
        trace, nfev = capture_direct_vacuum(
            model, START, spec.rho_up, max_step=MAX_STEP, cpu_limit=float(cpu_limit))
        return _prove(spec, model, trace, int(nfev), dependencies)


