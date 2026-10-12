"""Resumable witnesses for the uncovered subgap rows below mass.

The phase transport, horizon frame, Bloch witness and signed covariance
comparison are the existing owners. This module does not define a solver,
a tube, or a separation scale. Standard rows use the row-15 control. The
original subgap8/14_1 row 7 uses the documented wider tube and 16 defect
subdivisions. A row the existing control rejects is stored as a failure,
not as a zero bound and not under a retuned control.
"""
from contextlib import contextmanager
from hashlib import sha256
from io import BytesIO
import json
import math
from pathlib import Path
import resource
import signal
import sys
import time

import numpy as np
from flint import arb, ctx

from .evidence_io import publish_exclusive_directory, publish_exclusive_file
from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_ks_evaluation_binding import file_digest, hex_float, implementation_hashes
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_massive_jost_mixed_transport import subgap_mixed_phase_transport_bound
from .nsc_massive_jost_modes import solve_jost
from .nsc_massive_jost_transport_bound import (
    original_dense_solution, phase_transport_segments)
from .nsc_metric_horizon_frame import metric_horizon_frame, reflection_from_phase
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from .nsc_paired_horizon_preparation import PairedHorizonSeedMap
from .nsc_subgap_row_witness import (
    BITS, BLOCH_DEGREE, BLOCH_MAX_STEP, BLOCH_Y_START, DEFECT_SUBDIVISIONS,
    DEGREE, FRAME_ORDER, INNER_REQUEST, METRIC_TERMS, PAYLOAD_NAMES,
    PHASE_COLUMNS, SOLVER, TUBE,
    PhasePreparation, bloch_record, finalize_bloch, load_witness_arrays,
    phase_prefix_from_dense, phase_record, prepare_phase, replay_mode,
    requested_inner_radius, require_exact_source, require_outer_node,
    require_same_phase_segments, require_segments_match_prefix,
    segments_of_prefix, validate_phase_prefix)
from .nsc_subgap_source_covariance import BlochSource, capture_bloch, initial_bloch


SCHEMA = "NSC-SUBGAP-SOURCE-WINDOW-v1"
ROW_SCHEMA = "NSC-SUBGAP-SOURCE-ROW-v1"
FAILURE_SCHEMA = "NSC-SUBGAP-SOURCE-WINDOW-FAILURE-v1"
FREEZE_SCHEMA = "NSC-SUBGAP-SOURCE-WINDOW-PROOF-FREEZE-v1"
GROUP = 14
ANGULAR_SIGN = 1
FAMILY_SIGNED_ROWS = 1168
EXPECTED_POSITIVE_ROWS = 54
EXPECTED_SIGNED_ROWS = 108
LOW16_PANEL = "group14/low16_1"
LOW32_PANEL = "group14/low32_1"
SUBGAP_PANEL = "subgap8/14_1"
# Contiguous archived row intervals. Row 0 and row 15 of low16 stay outside.
PANEL_ROWS = (
    (LOW16_PANEL, tuple(range(1, 15))),
    (LOW32_PANEL, tuple(range(0, 32))),
    (SUBGAP_PANEL, tuple(range(0, 8))),
)
ENERGY_ENDPOINTS = {
    (LOW16_PANEL, 1): "0x1.c60a99e906500p-8",
    (LOW16_PANEL, 14): "0x1.f1cfab30b7cd8p-3",
    (LOW32_PANEL, 0): "0x1.010cf92bcf8b8p-2",
    (LOW32_PANEL, 31): "0x1.ff79836a183a4p-1",
    (SUBGAP_PANEL, 0): "0x1.02e6bb940d290p+0",
    (SUBGAP_PANEL, 7): "0x1.8f38f9b035a88p+0",
}
EXCLUDED_LOW16 = (0, 15)
EXCLUDED_LOW32 = tuple(range(32, 41))
EXCLUDED_LOW32_ENERGY = {
    32: "0x1.92f9814c33825p+0",
    40: "0x1.fae56f00f6a85p+0",
}
PI_OVER_TWO_HEX = "0x1.921fb54442d18p+0"
# Direct-vacuum window already authenticated elsewhere. Not selected here.
DIRECT_LOW16 = tuple(range(16, 48))
DIRECT_LOW32 = tuple(range(41, 96))
NEAR_THRESHOLD = (SUBGAP_PANEL, 7)
WIDE_TUBE = "0.001"
WIDE_SUBDIVISIONS = 16
STANDARD_PROFILE = "row15-mixed-transport"
WIDE_PROFILE = "near-threshold-wider-tube"
PHASE_SEPARATION = 1e-4
LOOSE_COVARIANCE = 1e-3
PROOF_RELATIVE = "src/recursive_horizons/nsc_subgap_source_window.py"
MIXED_PILOT = "scripts/pilot_nsc_massive_jost_mixed_transport.py"
NEAR_PILOT = "scripts/pilot_nsc_near_threshold_source.py"
ROW_STATUS = (
    "OPEN: one original subgap row and its actual negative partner; "
    "unweighted covariance errors; local gate open"
)
OPEN_STATUS = (
    "OPEN: uncovered subgap window is incomplete; completed witnesses only; "
    "failed rows are not entered as zero; full upstream budget and local gate remain open"
)
WINDOW_STATUS = (
    "OPEN: finite uncovered subgap window and the actual negative partners "
    "are recorded; full upstream budget and local gate remain open"
)
FAILURE_STATUS = (
    "FAILED: the established subgap control rejected this row; "
    "no covariance bound is entered and the control was not retuned"
)
PROTECTED_PREFIXES = (
    "results/development/nsc-subgap-row15-upstream-v1.json",
    "results/development/artifacts/nsc-subgap-row15-upstream-v1.npz",
    "results/development/artifacts",
    "results/development/nsc-direct-source-window-v1",
    "results/development/nsc-middle-source-coverage-v1",
    "results/development/nsc-ks-source-operator-majorant-v1.json",
    "results/development/nsc-ks-source-operator-majorant-v2.json",
    "results/development/nsc-ks-source-operator-majorant-v3.json",
    "results/development/nsc-ks-source-operator-majorant-v4.json",
)
ALLOWED_ROOT_NAMES = {"proof-freeze.json", "rows", "failures", "capture-log.txt"}


class EstablishedControlRejected(ArithmeticError):
    """The existing solver or enclosure rejected the row. Not a zero bound."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = dict(details or {})


class SubgapRow:
    """One archived positive row and its opposite-angular negative partner."""

    __slots__ = (
        "panel", "row", "batch_rows", "local_index", "energy", "negative_energy",
        "energy_fiber", "negative_energy_fiber", "mass", "angular",
        "negative_angular", "rho_up", "horizon_rho", "kappa", "omega",
        "positive_columns", "negative_columns", "positive_covariance",
        "negative_covariance", "positive_source_digest", "negative_source_digest",
        "positive_preparation_digest", "negative_preparation_digest", "weight",
        "column_weight", "background")

    def __init__(self, **kwargs):
        for name in self.__slots__:
            setattr(self, name, kwargs[name])


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


def failure_relative(panel, row):
    return f"failures/{checkpoint_name(panel, row)}"


def settings_for(panel, row):
    """Established controls only. Callers cannot substitute a tube."""
    if (panel, row) == NEAR_THRESHOLD:
        return {
            "transport_profile": WIDE_PROFILE,
            "tube": WIDE_TUBE,
            "defect_subdivisions": WIDE_SUBDIVISIONS,
        }
    return {
        "transport_profile": STANDARD_PROFILE,
        "tube": TUBE,
        "defect_subdivisions": DEFECT_SUBDIVISIONS,
    }


def _require_control_constants():
    expected_solver = {
        "order": 8, "radial_collar": 1e-10, "rtol": 2e-13, "atol": 2e-15}
    if (TUBE != "0.00001" or DEFECT_SUBDIVISIONS != 4 or SOLVER != expected_solver
            or DEGREE != 16 or METRIC_TERMS != 48 or BITS != 192
            or FRAME_ORDER != 16 or BLOCH_DEGREE != 12 or BLOCH_Y_START != -18.0
            or BLOCH_MAX_STEP != 0.025 or WIDE_TUBE != "0.001"
            or WIDE_SUBDIVISIONS != 16 or PHASE_SEPARATION != 1e-4
            or LOOSE_COVARIANCE != 1e-3 or INNER_REQUEST != 1.01e-4):
        raise ValueError("established subgap control constants changed")
    if tuple(PHASE_COLUMNS)[:5] != ("y0", "y1", "theta0", "theta1", "F0"):
        raise ValueError("phase prefix layout changed")


def _require_documented_controls(archive_root):
    """Lock the pilots that defined the two successful controls. Do not import
    the near-threshold pilot: importing it parses a required command flag.
    """
    root = Path(archive_root)
    near = (root / NEAR_PILOT).read_text()
    mixed = (root / MIXED_PILOT).read_text()
    if 'tube="0.001"' not in near or "defect_subdivisions=16" not in near:
        raise ValueError("near-threshold wider tube control changed")
    if ("LOOSE_COVARIANCE = 1e-3" not in mixed or "defect_subdivisions=4" not in mixed
            or "phase_float < 1e-4" not in mixed):
        raise ValueError("row15 mixed-transport control changed")


def _budgets(row_budget, cpu_budget):
    if isinstance(row_budget, bool) or type(row_budget) is not int or row_budget < 0:
        raise ValueError("explicit nonnegative integer row budget required")
    if (isinstance(cpu_budget, bool) or not isinstance(cpu_budget, (int, float))
            or not math.isfinite(float(cpu_budget)) or float(cpu_budget) <= 0):
        raise ValueError("explicit positive finite CPU budget required")
    return row_budget, float(cpu_budget)


def _hex(value):
    return hex_float(value)


def _fiber_hex(values):
    fiber = np.asarray(values, float)
    if fiber.shape != (3,) or not np.isfinite(fiber).all():
        raise ValueError("three-entry energy fiber required")
    return [_hex(item) for item in fiber]


def _frozen(value, dtype):
    array = np.array(value, dtype=dtype, copy=True, order="C")
    array.setflags(write=False)
    return array


def _remember(seen, panel, row, energy_hex):
    key = (panel, row)
    energy_key = (panel, energy_hex)
    if key in seen or energy_key in seen:
        raise ValueError("duplicate source identity")
    seen.add(key)
    seen.add(energy_key)
    return key


def _pi_over_two():
    value = float(np.pi / 2)
    if value.hex() != PI_OVER_TWO_HEX:
        raise ValueError("pi/2 hex changed")
    return value


def _row_from_batches(positive, partner, local, config, archive, background, horizon_rho):
    if positive.group != GROUP or positive.angular_sign != ANGULAR_SIGN:
        raise ValueError("subgap window must stay on original group14 angular+1")
    if partner.angular_sign != -1 or partner.energy_sign != -1:
        raise ValueError("actual opposite-angular negative partner required")
    if partner.original_panel != positive.original_panel or partner.rows != positive.rows:
        raise ValueError("negative partner must share the original panel rows")
    if positive.mass != np.pi / 2 or positive.angular != np.sqrt(5.0):
        raise ValueError("original group14 mass and angular label required")
    if (partner.mass != positive.mass or partner.angular != -positive.angular
            or partner.rho_up != positive.rho_up):
        raise ValueError("signed mass, angular label or rho changed")
    row_index = int(positive.rows[0] + local)
    sl = slice(3 * local, 3 * local + 3)
    fiber = np.asarray(positive.source.energies[sl], float)
    negative_fiber = np.asarray(partner.source.energies[sl], float)
    if fiber.shape != (3,) or not np.all(fiber == fiber[0]):
        raise ValueError("three-column energy fiber required")
    if not np.array_equal(negative_fiber, -fiber):
        raise ValueError("negative partner energies must equal -E")
    energy = float(fiber[0])
    if not math.isfinite(energy) or energy <= 0 or not positive.mass > energy:
        raise ValueError("subgap positive mass required")
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
    return SubgapRow(
        panel=positive.original_panel, row=row_index,
        batch_rows=(int(positive.rows[0]), int(positive.rows[1])),
        local_index=int(local), energy=energy, negative_energy=float(negative_fiber[0]),
        energy_fiber=tuple(float(item) for item in fiber),
        negative_energy_fiber=tuple(float(item) for item in negative_fiber),
        rho_up=float(positive.rho_up), mass=float(positive.mass),
        angular=float(positive.angular), negative_angular=float(partner.angular),
        positive_source_digest=positive.source.digest,
        negative_source_digest=partner.source.digest,
        positive_preparation_digest=positive.preparation_digest,
        negative_preparation_digest=partner.preparation_digest,
        positive_columns=positive_columns, negative_columns=negative_columns,
        positive_covariance=positive_covariance, negative_covariance=negative_covariance,
        weight=raw_weight, column_weight=float(column_weights[0]),
        kappa=float(config["surface_gravity"]), omega=float(config["omega"]),
        horizon_rho=float(horizon_rho), background=background)


def _indexed_positive_rows(entries):
    rows = {}
    for batch, _channel in entries:
        if batch.energy_sign <= 0:
            continue
        fibers = np.asarray(batch.source.energies, float).reshape(-1, 3)
        for local, fiber in enumerate(fibers):
            key = (batch.original_panel, int(batch.rows[0] + local))
            if key in rows:
                raise ValueError("duplicate original source row")
            rows[key] = (batch, local, float(fiber[0]))
    return rows


def _require_census(rows, entries):
    signed = sum(len(batch.source.energies) // 3 for batch, _channel in entries)
    if signed != FAMILY_SIGNED_ROWS:
        raise ValueError("family14_1 census changed")
    if len(rows) != EXPECTED_POSITIVE_ROWS:
        raise ValueError("uncovered subgap row count changed")
    by_panel = {panel: [] for panel, _wanted in PANEL_ROWS}
    for row in rows:
        if row.panel not in by_panel:
            raise ValueError("subgap window panels changed")
        by_panel[row.panel].append(row.row)
    for panel, wanted in PANEL_ROWS:
        if tuple(by_panel[panel]) != wanted:
            raise ValueError(panel + " uncovered rows changed")
    if [row.panel for row in rows] != [panel for panel, wanted in PANEL_ROWS for _row in wanted]:
        raise ValueError("subgap window order changed")
    pi_over_two = _pi_over_two()
    indexed = _indexed_positive_rows(entries)
    for row in rows:
        if row.energy >= pi_over_two or not row.mass > row.energy:
            raise ValueError("subgap positive mass required")
        if (row.panel, row.row) in {(LOW16_PANEL, item) for item in EXCLUDED_LOW16}:
            raise ValueError("row 0 and row 15 stay outside this window")
        if row.panel == LOW32_PANEL and row.row in EXCLUDED_LOW32:
            raise ValueError("low32 rows at or above pi/2 stay outside this window")
        if ((row.panel == LOW16_PANEL and row.row in DIRECT_LOW16)
                or (row.panel == LOW32_PANEL and row.row in DIRECT_LOW32)):
            raise ValueError("direct-window rows are not part of this subgap census")
    for key, expected in ENERGY_ENDPOINTS.items():
        matched = [row for row in rows if (row.panel, row.row) == key]
        if len(matched) != 1 or matched[0].energy.hex() != expected:
            raise ValueError("archived energy endpoint changed: " + key[0])
    for index in EXCLUDED_LOW16:
        key = (LOW16_PANEL, index)
        if key not in indexed or key in {(row.panel, row.row) for row in rows}:
            raise ValueError("preserved low16 row is missing from the archive")
    excluded_energies = []
    for index in EXCLUDED_LOW32:
        key = (LOW32_PANEL, index)
        if key not in indexed:
            raise ValueError("excluded low32 row is missing from the archive")
        _batch, _local, energy = indexed[key]
        if energy < pi_over_two or not energy < 2.0:
            raise ValueError("excluded low32 energy is not between pi/2 and 2")
        excluded_energies.append(energy)
        if key in {(row.panel, row.row) for row in rows}:
            raise ValueError("low32 rows at or above pi/2 stay outside this window")
    if (excluded_energies[0].hex() != EXCLUDED_LOW32_ENERGY[32]
            or excluded_energies[-1].hex() != EXCLUDED_LOW32_ENERGY[40]):
        raise ValueError("excluded low32 energy endpoints changed")
    rhos = {row.rho_up for row in rows}
    if len(rhos) != 1 or not next(iter(rhos)) >= 1.03:
        raise ValueError("upstream archive slice changed")


def subgap_catalogue(archive):
    """The 54 uncovered positive rows and their archived negative partners.

    Census failures raise before a caller is allowed to solve.
    """
    if not isinstance(archive, RetainedUpstreamArchive):
        raise TypeError("retained upstream archive required")
    _require_control_constants()
    config = archive.meta["config"]
    for name in ("surface_gravity", "omega", "horizon_rho"):
        if name not in config:
            raise ValueError("original source scales required")
    preparation = PairedHorizonSeedMap(
        config["horizon_rho"], config["surface_gravity"], config["omega"],
        config["horizon_offset"], config["scattering_tolerance"], config["outer_floor"])
    horizon_rho = float(preparation.horizon_rho)
    if horizon_rho != float(config["horizon_rho"]):
        raise ValueError("horizon label changed")
    background = preparation.background
    entries = archive.family_entries((GROUP, ANGULAR_SIGN))
    positives = {}
    partners = {}
    for batch, _channel in entries:
        key = (batch.original_panel, tuple(batch.rows))
        if batch.energy_sign > 0:
            if key in positives:
                raise ValueError("duplicate original source row")
            positives[key] = batch
        elif batch.energy_sign < 0:
            if key in partners:
                raise ValueError("duplicate original source row")
            partners[key] = batch
    selected = []
    seen = set()
    for panel, wanted in PANEL_ROWS:
        for row_index in wanted:
            owners = [key for key in positives if key[0] == panel and key[1][0] <= row_index < key[1][1]]
            if len(owners) != 1 or owners[0] not in partners:
                raise ValueError("actual negative partner required")
            batch = positives[owners[0]]
            local = row_index - int(batch.rows[0])
            row = _row_from_batches(
                batch, partners[owners[0]], local, config, archive, background, horizon_rho)
            if row.row != row_index or row.panel != panel:
                raise ValueError("subgap row identity changed")
            _remember(seen, row.panel, row.row, _hex(row.energy))
            selected.append(row)
    _require_census(tuple(selected), entries)
    return tuple(selected)


def dependency_identity(archive_root, archive):
    _require_documented_controls(archive_root)
    root = Path(archive_root)
    proof = PROOF_RELATIVE
    source_hashes = implementation_hashes(root, owners=(proof,))
    inputs = dict(archive.input_hashes)
    controls = {
        MIXED_PILOT: file_digest(root / MIXED_PILOT),
        NEAR_PILOT: file_digest(root / NEAR_PILOT),
    }
    if not source_hashes or not inputs or proof not in source_hashes:
        raise ValueError("source and dependency identity required")
    return {
        "proof": proof,
        "proof_sha256": source_hashes[proof],
        "source_hashes": source_hashes,
        "archive_digests": inputs,
        "control_hashes": controls,
    }


def _prepare_wider_tube(prefix, *, energy, mass, angular, background, outer_radius,
                        horizon_rho, kappa, rho_up):
    """Same frame and Bloch setup as ``prepare_phase``, with the documented tube.

    The row-15 owner hardcodes its tube. Only subgap8/14_1 row 7 reaches this
    path, and only with tube 0.001 and 16 subdivisions.
    """
    prefix = validate_phase_prefix(prefix)
    require_outer_node(prefix, outer_radius, horizon_rho)
    mode = replay_mode(
        prefix, energy=energy, mass=mass, angular=angular, background=background,
        outer_radius=outer_radius)
    inner = requested_inner_radius(horizon_rho)
    segments = phase_transport_segments(
        original_dense_solution(mode.run), horizon_rho, inner)
    require_segments_match_prefix(segments, prefix)
    with ctx.workprec(BITS):
        mass_hull = arb(mass).union(arb.pi() / 2)
        angular_hull = arb(angular).union(arb(5).sqrt())
        bound = subgap_mixed_phase_transport_bound(
            energy, mass_hull, angular_hull, background, mode=mode,
            inner_radius=inner, degree=DEGREE, metric_terms=METRIC_TERMS,
            defect_subdivisions=WIDE_SUBDIVISIONS, tube=WIDE_TUBE, bits=BITS)
    if float(segments[-1].y_end).hex() != bound["inner_y_node_hex"]:
        raise ArithmeticError("transported node does not match the captured endpoint")
    if int(bound["cells"]) != len(segments):
        raise ArithmeticError("transported cell count does not match the saved prefix")
    if int(bound["defect_subdivisions"]) != WIDE_SUBDIVISIONS:
        raise ValueError("transport subdivisions do not match the established control")
    phase_error = restored_upper(bound["phase_error_inner_upper"])
    with ctx.workprec(BITS):
        rho = arb(horizon_rho) + arb(float(segments[-1].y_end)).exp()
        theta = arb(float(segments[-1].theta_end)) + arb(0, phase_error)
        frame = metric_horizon_frame(
            energy, mass_hull, angular_hull, horizon_rho, order=FRAME_ORDER, bits=BITS)
        reflection, distance, tail = reflection_from_phase(frame, rho, theta)
        initial = initial_bloch(frame, reflection, kappa, BLOCH_Y_START)
        model = BlochSource(frame.q, energy, mass_hull, angular_hull, bits=BITS)
        target = (frame.q - arb.pi() / 2 - arb(rho_up).atan()).log()
    return PhasePreparation(
        bound=bound, segments=segments, mass_hull=mass_hull, angular_hull=angular_hull,
        frame=frame, reflection=reflection, distance=distance, tail=tail,
        initial=initial, model=model, target=target)


def prepare_recorded_phase(prefix, *, energy, mass, angular, background, outer_radius,
                           horizon_rho, kappa, rho_up, defect_subdivisions, tube):
    """Dispatch to the row-15 owner or the one documented wider tube."""
    if isinstance(defect_subdivisions, bool) or type(defect_subdivisions) is not int:
        raise ValueError("established defect subdivisions required")
    if isinstance(tube, bool) or not isinstance(tube, str):
        raise ValueError("established phase tube required")
    _require_control_constants()
    if defect_subdivisions == DEFECT_SUBDIVISIONS and tube == TUBE:
        prepared = prepare_phase(
            prefix, energy=energy, mass=mass, angular=angular, background=background,
            outer_radius=outer_radius, horizon_rho=horizon_rho, kappa=kappa,
            rho_up=rho_up)
        if int(prepared.bound["defect_subdivisions"]) != DEFECT_SUBDIVISIONS:
            raise ValueError("transport subdivisions do not match the established control")
        return prepared
    if defect_subdivisions == WIDE_SUBDIVISIONS and tube == WIDE_TUBE:
        return _prepare_wider_tube(
            prefix, energy=energy, mass=mass, angular=angular, background=background,
            outer_radius=outer_radius, horizon_rho=horizon_rho, kappa=kappa,
            rho_up=rho_up)
    raise ValueError("retuning an established subgap control is refused")


def _separated(record):
    phase = float(restored_upper(record["phase_error_inner_upper"]))
    positive = float(restored_upper(record["positive_covariance_error_upper"]))
    negative = float(restored_upper(record["negative_covariance_error_upper"]))
    return (phase < PHASE_SEPARATION and positive < LOOSE_COVARIANCE
            and negative < LOOSE_COVARIANCE)


def _comparisons(spec, positive_upper, negative_upper):
    return [
        {
            "energy_sign": 1, "angular_sign": 1,
            "energy_hex": _hex(spec.energy),
            "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
            "source_digest": spec.positive_source_digest,
            "preparation_digest": spec.positive_preparation_digest,
            "covariance_error_upper": positive_upper,
        },
        {
            "energy_sign": -1, "angular_sign": -1,
            "energy_hex": _hex(spec.negative_energy),
            "energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
            "source_digest": spec.negative_source_digest,
            "preparation_digest": spec.negative_preparation_digest,
            "covariance_error_upper": negative_upper,
        },
    ]


def assemble_witness(spec, prefix, trace, prepared, validated, positive_error,
                     negative_error, *, outer_radius, phase_nfev, bloch_nfev,
                     dependencies):
    """Serialize one finished witness. The phase proof has already run once."""
    settings = settings_for(spec.panel, spec.row)
    arrays = {
        "phase_prefix": np.ascontiguousarray(prefix, dtype=np.float64),
        "bloch_trace": np.ascontiguousarray(trace, dtype=np.float64),
        "positive_columns": np.ascontiguousarray(spec.positive_columns),
        "negative_columns": np.ascontiguousarray(spec.negative_columns),
        "positive_covariance": np.ascontiguousarray(spec.positive_covariance),
        "negative_covariance": np.ascontiguousarray(spec.negative_covariance),
        "outer_radius": np.array([float(outer_radius)], dtype=np.float64),
        "phase_nfev": np.array([int(phase_nfev)], dtype=np.int64),
        "bloch_nfev": np.array([int(bloch_nfev)], dtype=np.int64),
    }
    if set(arrays) != set(PAYLOAD_NAMES):
        raise ValueError("complete row witness payload required")
    raw = deterministic_npz_bytes(arrays)
    relative = checkpoint_relative(spec.panel, spec.row)
    positive_upper = exact_upper(positive_error)
    negative_upper = exact_upper(negative_error)
    record = {
        "schema": ROW_SCHEMA,
        "coverage_schema": SCHEMA,
        "status": ROW_STATUS,
        "group": GROUP,
        "panel": spec.panel,
        "row": spec.row,
        "batch_rows": [spec.batch_rows[0], spec.batch_rows[1]],
        "local_index": spec.local_index,
        "source_column_slice": [3 * spec.local_index, 3 * spec.local_index + 3],
        "energy_sign": 1,
        "negative_energy_sign": -1,
        "angular_sign": 1,
        "negative_angular_sign": -1,
        "energy_hex": _hex(spec.energy),
        "negative_energy_hex": _hex(spec.negative_energy),
        "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
        "negative_energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
        "mass_hex": _hex(spec.mass),
        "positive_angular_hex": _hex(spec.angular),
        "negative_angular_hex": _hex(spec.negative_angular),
        "rho_up_hex": _hex(spec.rho_up),
        "horizon_rho_hex": _hex(spec.horizon_rho),
        "kappa_hex": _hex(spec.kappa),
        "omega_hex": _hex(spec.omega),
        "outer_radius_hex": _hex(outer_radius),
        "inner_request_offset_hex": _hex(INNER_REQUEST),
        "weight_hex": _hex(spec.weight),
        "column_weight_hex": _hex(spec.column_weight),
        "transport_profile": settings["transport_profile"],
        "tube": settings["tube"],
        "defect_subdivisions": settings["defect_subdivisions"],
        "solver": {
            "method": "DOP853",
            "order": SOLVER["order"],
            "radial_collar_hex": _hex(SOLVER["radial_collar"]),
            "rtol_hex": _hex(SOLVER["rtol"]),
            "atol_hex": _hex(SOLVER["atol"]),
        },
        "transport": {
            "degree": DEGREE,
            "metric_terms": METRIC_TERMS,
            "defect_subdivisions": settings["defect_subdivisions"],
            "tube": settings["tube"],
            "bits": BITS,
            "frame_order": FRAME_ORDER,
        },
        "bloch_settings": {
            "y_start_hex": _hex(BLOCH_Y_START),
            "max_step_hex": _hex(BLOCH_MAX_STEP),
            "degree": BLOCH_DEGREE,
            "bits": BITS,
        },
        "separation_scales": {
            "phase_upper_hex": _hex(PHASE_SEPARATION),
            "covariance_upper_hex": _hex(LOOSE_COVARIANCE),
        },
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "phase_columns": list(PHASE_COLUMNS),
        "phase_prefix_rows": int(len(prefix)),
        "bloch_trace_rows": int(len(trace)),
        "phase_nfev": int(phase_nfev),
        "bloch_nfev": int(bloch_nfev),
        "covariance_error_unweighted": True,
        "quadrature_weight_applied": False,
        "negative_error_copied_from_positive": False,
        "negative_complement_applied": True,
        "phase_component_only": True,
        "amplitude_component_retained": False,
        "dop853_dense_output_reconstructed": True,
        "replay_uses_saved_phase_and_bloch_witness": True,
        "preparation_jost_rerun_on_replay": False,
        "bloch_ode_rerun_on_replay": False,
        "source_columns_replaced": False,
        "horizon_sewing_error": None,
        "frame_error_included_in_initial_bloch": True,
        "separate_frame_error_addition_required": False,
        "physical_rho1_source_error": None,
        "n_beta_aggregate": None,
        "physical_upstream_budget_component": None,
        "all_source_families": False,
        "method_too_loose": False,
        "physical_local_gate": "OPEN",
        "proof_sha256": dependencies["proof_sha256"],
        "source_hashes": dependencies["source_hashes"],
        "archive_digests": dependencies["archive_digests"],
        "control_hashes": dependencies["control_hashes"],
    }
    record.update(phase_record(prepared))
    record.update(bloch_record(validated, positive_error, negative_error))
    record["signed_comparisons"] = _comparisons(
        spec, record["positive_covariance_error_upper"],
        record["negative_covariance_error_upper"])
    if record["phase_cells"] != record["phase_prefix_rows"]:
        raise ArithmeticError("transported cell count does not match the saved prefix")
    if record["bloch_cells"] != record["bloch_trace_rows"]:
        raise ArithmeticError("validated Bloch cells do not match the saved trace")
    if int(record["transport"]["defect_subdivisions"]) != int(prepared.bound["defect_subdivisions"]):
        raise ValueError("transport subdivisions do not match the established control")
    if record["positive_covariance_error_upper"] != positive_upper:
        raise ArithmeticError("positive covariance bound changed")
    if record["negative_covariance_error_upper"] != negative_upper:
        raise ArithmeticError("negative covariance bound changed")
    if not _separated(record):
        raise EstablishedControlRejected(
            "method too loose under the established separation scales",
            details={
                "phase_error_inner_upper": record["phase_error_inner_upper"],
                "positive_covariance_error_upper": record["positive_covariance_error_upper"],
                "negative_covariance_error_upper": record["negative_covariance_error_upper"],
            })
    record["payload"] = {
        "path": f"{relative}/witness.npz",
        "sha256": sha256(raw).hexdigest(),
        "bytes": len(raw),
    }
    return record, raw


def prove_from_witness(spec, prefix, trace, *, outer_radius, phase_nfev, bloch_nfev,
                       dependencies, reject):
    """Validate one saved phase prefix and Bloch trace. This does not solve."""
    settings = settings_for(spec.panel, spec.row)
    try:
        prepared = prepare_recorded_phase(
            prefix, energy=spec.energy, mass=spec.mass, angular=spec.angular,
            background=spec.background, outer_radius=outer_radius,
            horizon_rho=spec.horizon_rho, kappa=spec.kappa, rho_up=spec.rho_up,
            defect_subdivisions=settings["defect_subdivisions"], tube=settings["tube"])
        validated, positive_error, negative_error = finalize_bloch(
            prepared, trace, spec.positive_columns, spec.positive_covariance,
            spec.negative_columns, spec.negative_covariance)
    except EstablishedControlRejected:
        raise
    except ArithmeticError as exc:
        if reject:
            raise EstablishedControlRejected(str(exc)) from exc
        raise
    return assemble_witness(
        spec, prefix, trace, prepared, validated, positive_error, negative_error,
        outer_radius=outer_radius, phase_nfev=phase_nfev, bloch_nfev=bloch_nfev,
        dependencies=dependencies)


def capture_row(spec, dependencies, *, cpu_limit):
    """One positive Jost solve, one phase proof and one Bloch capture.

    Replay does not call this. The phase owner runs once, on the saved prefix.
    """
    if not math.isfinite(cpu_limit) or cpu_limit <= 0:
        raise TimeoutError("subgap source-window CPU budget exhausted")
    settings = settings_for(spec.panel, spec.row)
    try:
        mode = solve_jost(
            spec.background, spec.energy, spec.mass, spec.angular, **SOLVER)
    except ArithmeticError as exc:
        raise EstablishedControlRejected(str(exc)) from exc
    dense = original_dense_solution(mode.run)
    inner = requested_inner_radius(spec.horizon_rho)
    prefix = phase_prefix_from_dense(dense, spec.horizon_rho, inner)
    require_same_phase_segments(
        phase_transport_segments(dense, spec.horizon_rho, inner),
        segments_of_prefix(prefix, spec.horizon_rho, inner))
    outer_radius = float(mode.outer_radius)
    try:
        prepared = prepare_recorded_phase(
            prefix, energy=spec.energy, mass=spec.mass, angular=spec.angular,
            background=spec.background, outer_radius=outer_radius,
            horizon_rho=spec.horizon_rho, kappa=spec.kappa, rho_up=spec.rho_up,
            defect_subdivisions=settings["defect_subdivisions"], tube=settings["tube"])
        trace, bloch_nfev = capture_bloch(
            prepared.model, prepared.initial, BLOCH_Y_START,
            float(prepared.target.mid()), max_step=BLOCH_MAX_STEP)
        validated, positive_error, negative_error = finalize_bloch(
            prepared, trace, spec.positive_columns, spec.positive_covariance,
            spec.negative_columns, spec.negative_covariance)
    except EstablishedControlRejected:
        raise
    except ArithmeticError as exc:
        raise EstablishedControlRejected(str(exc)) from exc
    return assemble_witness(
        spec, prefix, trace, prepared, validated, positive_error, negative_error,
        outer_radius=outer_radius, phase_nfev=int(mode.run.nfev),
        bloch_nfev=int(bloch_nfev), dependencies=dependencies)


def _failure_record(spec, dependencies, exc):
    details = getattr(exc, "details", {}) or {}
    message = str(exc)
    if not message:
        raise ValueError("failed row requires the rejecting control's message")
    for name in (
            "phase_error_inner_upper", "positive_covariance_error_upper",
            "negative_covariance_error_upper"):
        if name in details and details[name] is None:
            raise ValueError("a failed row must not replace an error with zero")
    settings = settings_for(spec.panel, spec.row)
    return {
        "schema": FAILURE_SCHEMA,
        "coverage_schema": SCHEMA,
        "status": FAILURE_STATUS,
        "group": GROUP,
        "panel": spec.panel,
        "row": spec.row,
        "energy_hex": _hex(spec.energy),
        "negative_energy_hex": _hex(spec.negative_energy),
        "transport_profile": settings["transport_profile"],
        "tube": settings["tube"],
        "defect_subdivisions": settings["defect_subdivisions"],
        "exception_type": type(exc).__name__,
        "exception_message": message,
        "phase_error_inner_upper": details.get("phase_error_inner_upper"),
        "positive_covariance_error_upper": details.get("positive_covariance_error_upper"),
        "negative_covariance_error_upper": details.get("negative_covariance_error_upper"),
        "entered_as_zero": False,
        "covariance_bound_used_in_aggregate": False,
        "retuned": False,
        "method_too_loose": message.startswith("method too loose"),
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "weight_hex": _hex(spec.weight),
        "quadrature_weight_applied": False,
        "proof_sha256": dependencies["proof_sha256"],
        "physical_upstream_budget_component": None,
        "physical_local_gate": "OPEN",
    }


def _require_bindings(record, spec, dependencies):
    settings = settings_for(spec.panel, spec.row)
    expected = {
        "schema": ROW_SCHEMA,
        "coverage_schema": SCHEMA,
        "panel": spec.panel,
        "row": spec.row,
        "group": GROUP,
        "energy_hex": _hex(spec.energy),
        "negative_energy_hex": _hex(spec.negative_energy),
        "energy_fiber_hex": _fiber_hex(spec.energy_fiber),
        "negative_energy_fiber_hex": _fiber_hex(spec.negative_energy_fiber),
        "mass_hex": _hex(spec.mass),
        "positive_angular_hex": _hex(spec.angular),
        "negative_angular_hex": _hex(spec.negative_angular),
        "rho_up_hex": _hex(spec.rho_up),
        "transport_profile": settings["transport_profile"],
        "tube": settings["tube"],
        "defect_subdivisions": settings["defect_subdivisions"],
        "positive_source_digest": spec.positive_source_digest,
        "negative_source_digest": spec.negative_source_digest,
        "positive_preparation_digest": spec.positive_preparation_digest,
        "negative_preparation_digest": spec.negative_preparation_digest,
        "weight_hex": _hex(spec.weight),
        "column_weight_hex": _hex(spec.column_weight),
        "quadrature_weight_applied": False,
        "negative_error_copied_from_positive": False,
        "covariance_error_unweighted": True,
        "method_too_loose": False,
        "physical_local_gate": "OPEN",
        "proof_sha256": dependencies["proof_sha256"],
    }
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError("saved subgap witness does not match the archived row")
    if record.get("source_hashes") != dependencies["source_hashes"]:
        raise ValueError("frozen proof code changed")
    if record.get("archive_digests") != dependencies["archive_digests"]:
        raise ValueError("archive binding changed")
    if record.get("control_hashes") != dependencies["control_hashes"]:
        raise ValueError("established control changed")
    if record.get("preparation_jost_rerun_on_replay") is not False:
        raise ValueError("replay must not solve")
    if record.get("bloch_ode_rerun_on_replay") is not False:
        raise ValueError("replay must not solve")


def _load_payload(record, payload, spec):
    described = record.get("payload")
    if (not isinstance(described, dict) or described.get("path") != f"{checkpoint_relative(spec.panel, spec.row)}/witness.npz"
            or described.get("bytes") != len(payload)
            or described.get("sha256") != sha256(payload).hexdigest()):
        raise ValueError("row witness payload hash changed")
    with np.load(BytesIO(payload), allow_pickle=False) as saved:
        arrays = {name: saved[name].copy() for name in saved.files}
    loaded = load_witness_arrays(arrays)
    require_exact_source(
        loaded, spec.positive_columns, spec.negative_columns,
        spec.positive_covariance, spec.negative_covariance)
    return loaded


def authenticate_checkpoint(spec, dirname, record, payload, dependencies):
    """Replay one saved witness. This never calls the Jost solver."""
    if dirname != checkpoint_name(spec.panel, spec.row):
        raise ValueError("checkpoint name does not match row identity")
    _require_bindings(record, spec, dependencies)
    loaded = _load_payload(record, payload, spec)
    rebuilt, rebuilt_payload = prove_from_witness(
        spec, loaded["phase_prefix"], loaded["bloch_trace"],
        outer_radius=loaded["outer_radius"], phase_nfev=loaded["phase_nfev"],
        bloch_nfev=loaded["bloch_nfev"], dependencies=dependencies, reject=False)
    if rebuilt_payload != payload or encode_record(rebuilt) != encode_record(record):
        raise ValueError("saved witness changed")
    return rebuilt


def _authenticate_failure(spec, dirname, record, dependencies):
    if dirname != checkpoint_name(spec.panel, spec.row):
        raise ValueError("checkpoint name does not match row identity")
    settings = settings_for(spec.panel, spec.row)
    if (record.get("schema") != FAILURE_SCHEMA or record.get("panel") != spec.panel
            or record.get("row") != spec.row or record.get("energy_hex") != _hex(spec.energy)
            or record.get("negative_energy_hex") != _hex(spec.negative_energy)
            or record.get("transport_profile") != settings["transport_profile"]
            or record.get("tube") != settings["tube"]
            or record.get("defect_subdivisions") != settings["defect_subdivisions"]
            or record.get("entered_as_zero") is not False
            or record.get("covariance_bound_used_in_aggregate") is not False
            or record.get("retuned") is not False
            or record.get("quadrature_weight_applied") is not False
            or record.get("physical_local_gate") != "OPEN"
            or record.get("proof_sha256") != dependencies["proof_sha256"]
            or record.get("positive_source_digest") != spec.positive_source_digest
            or record.get("negative_source_digest") != spec.negative_source_digest
            or not isinstance(record.get("exception_message"), str)
            or not record["exception_message"]
            or record.get("exception_type") not in {"EstablishedControlRejected", "ArithmeticError"}):
        raise ValueError("failed subgap row does not match the archived identity")
    for name in (
            "phase_error_inner_upper", "positive_covariance_error_upper",
            "negative_covariance_error_upper"):
        value = record.get(name)
        if value is None:
            continue
        upper = restored_upper(value)
        if not upper.is_finite() or not upper >= 0:
            raise ValueError("failed row bound is not a real nonnegative diagnostic")
    return record


def _canonical_directory(path, *, create):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError("absolute coverage path required")
    if create:
        path.mkdir(parents=True, exist_ok=True)
    elif not path.exists():
        return path.resolve(strict=False)
    return path.resolve(strict=True)


def _refuse_foreign_campaign(archive_root, output_root):
    root = Path(archive_root).resolve(strict=True)
    resolved = Path(output_root).resolve(strict=False)
    if not resolved.is_relative_to(root):
        raise ValueError("coverage output must stay inside the laboratory")
    relative = resolved.relative_to(root).as_posix()
    for banned in PROTECTED_PREFIXES:
        if relative == banned or relative.startswith(banned + "/"):
            raise ValueError("refusing to write a historical record")


def _require_root_names(output_root):
    if not output_root.exists():
        return
    names = {item.name for item in output_root.iterdir()}
    if not names <= ALLOWED_ROOT_NAMES:
        raise ValueError("unexpected coverage entry")


def freeze_payload(archive_root):
    _require_control_constants()
    _require_documented_controls(archive_root)
    digest = file_digest(Path(archive_root) / PROOF_RELATIVE)
    return {
        "schema": FREEZE_SCHEMA,
        "proof": PROOF_RELATIVE,
        "sha256": digest,
        "positive_rows": EXPECTED_POSITIVE_ROWS,
        "signed_rows": EXPECTED_SIGNED_ROWS,
        "standard_profile": settings_for(LOW16_PANEL, 1),
        "near_threshold_profile": settings_for(*NEAR_THRESHOLD),
        "panels": [
            {"panel": panel, "first_row": rows[0], "last_row": rows[-1], "count": len(rows)}
            for panel, rows in PANEL_ROWS
        ],
        "energy_endpoints": {f"{panel}:{row}": value for (panel, row), value in ENERGY_ENDPOINTS.items()},
        "preserved_low16_rows": list(EXCLUDED_LOW16),
        "excluded_low32_rows": list(EXCLUDED_LOW32),
        "phase_separation_hex": _hex(PHASE_SEPARATION),
        "loose_covariance_hex": _hex(LOOSE_COVARIANCE),
        "solver": {
            "order": SOLVER["order"],
            "radial_collar_hex": _hex(SOLVER["radial_collar"]),
            "rtol_hex": _hex(SOLVER["rtol"]),
            "atol_hex": _hex(SOLVER["atol"]),
        },
    }


def ensure_proof_freeze(archive_root, output_root):
    """Write the proof hash once, before any row record. A later edit is refused."""
    payload = encode_record(freeze_payload(archive_root))
    path = Path(output_root) / "proof-freeze.json"
    if not path.exists():
        publish_exclusive_file(output_root, "proof-freeze.json", payload)
    elif path.read_bytes() != payload:
        raise ValueError("frozen proof code changed after the census lock")
    return json.loads(payload)


def publish_checkpoint(output_root, record, payload):
    name = checkpoint_name(record["panel"], record["row"])
    described = record.get("payload")
    if not isinstance(described, dict):
        raise ValueError("missing payload")
    digest = sha256(payload).hexdigest()
    if described.get("sha256") != digest or described.get("bytes") != len(payload):
        raise ValueError("payload hash changed")
    if described.get("path") != f"rows/{name}/witness.npz":
        raise ValueError("payload hash changed")
    root = _canonical_directory(output_root, create=True)
    (root / "rows").mkdir(exist_ok=True)
    publish_exclusive_directory(root, f"rows/{name}", {
        "record.json": encode_record(record),
        "witness.npz": payload,
    })


def publish_failure(output_root, record):
    name = checkpoint_name(record["panel"], record["row"])
    if record.get("entered_as_zero") is not False:
        raise ValueError("a failed row must not be entered as zero")
    root = _canonical_directory(output_root, create=True)
    (root / "failures").mkdir(exist_ok=True)
    publish_exclusive_directory(root, f"failures/{name}", {
        "failure.json": encode_record(record),
    })


def _read_json_file(path, name):
    raw = (path / name).read_bytes()
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError("corrupt row directory") from exc
    if not isinstance(record, dict) or encode_record(record) != raw:
        raise ValueError("corrupt row directory")
    return record, raw


def _row_key(record):
    panel, row = record.get("panel"), record.get("row")
    if (not isinstance(panel, str) or isinstance(row, bool) or type(row) is not int):
        raise ValueError("changed bound")
    return panel, row


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
        names = {item.name for item in path.iterdir()}
        if names != {"record.json", "witness.npz"}:
            raise ValueError("corrupt row directory")
        record, _raw = _read_json_file(path, "record.json")
        panel, row = _row_key(record)
        if path.name != checkpoint_name(panel, row):
            raise ValueError("checkpoint name does not match row identity")
        fiber = record.get("energy_fiber_hex")
        if not isinstance(fiber, list) or len(fiber) != 3:
            raise ValueError("changed bound")
        _remember(seen, panel, row, fiber[0])
        loaded.append((path.name, record, (path / "witness.npz").read_bytes()))
    return loaded


def load_failures(output_root):
    root = _canonical_directory(output_root, create=False)
    failures = root / "failures"
    if not failures.exists():
        return []
    if not failures.is_dir():
        raise ValueError("corrupt row directory")
    loaded = []
    seen = set()
    for path in sorted(failures.iterdir(), key=lambda item: item.name):
        if not path.is_dir() or not path.name.startswith("row-"):
            raise ValueError("unexpected coverage entry")
        names = {item.name for item in path.iterdir()}
        if names != {"failure.json"}:
            raise ValueError("corrupt row directory")
        record, _raw = _read_json_file(path, "failure.json")
        panel, row = _row_key(record)
        if path.name != checkpoint_name(panel, row):
            raise ValueError("checkpoint name does not match row identity")
        _remember(seen, panel, row, record.get("energy_hex"))
        loaded.append((path.name, record))
    return loaded


def _summary(record):
    return {
        "panel": record["panel"],
        "row": record["row"],
        "checkpoint": record["payload"]["path"].rsplit("/", 1)[0],
        "energy_hex": record["energy_hex"],
        "negative_energy_hex": record["negative_energy_hex"],
        "transport_profile": record["transport_profile"],
        "tube": record["tube"],
        "defect_subdivisions": record["defect_subdivisions"],
        "phase_cells": record["phase_cells"],
        "bloch_cells": record["bloch_cells"],
        "phase_error_inner_upper": record["phase_error_inner_upper"],
        "initial_bloch_error_upper": record["initial_bloch_error_upper"],
        "endpoint_bridge_upper": record["endpoint_bridge_upper"],
        "positive_covariance_error_upper": record["positive_covariance_error_upper"],
        "negative_covariance_error_upper": record["negative_covariance_error_upper"],
        "positive_source_digest": record["positive_source_digest"],
        "negative_source_digest": record["negative_source_digest"],
        "positive_preparation_digest": record["positive_preparation_digest"],
        "negative_preparation_digest": record["negative_preparation_digest"],
        "weight_hex": record["weight_hex"],
        "column_weight_hex": record["column_weight_hex"],
        "signed_comparisons": record["signed_comparisons"],
        "quadrature_weight_applied": False,
        "negative_error_copied_from_positive": False,
        "physical_upstream_budget_component": None,
        "physical_local_gate": "OPEN",
    }


def _failure_summary(record):
    return {
        "panel": record["panel"],
        "row": record["row"],
        "checkpoint": failure_relative(record["panel"], record["row"]),
        "energy_hex": record["energy_hex"],
        "negative_energy_hex": record["negative_energy_hex"],
        "transport_profile": record["transport_profile"],
        "tube": record["tube"],
        "defect_subdivisions": record["defect_subdivisions"],
        "exception_type": record["exception_type"],
        "exception_message": record["exception_message"],
        "phase_error_inner_upper": record.get("phase_error_inner_upper"),
        "positive_covariance_error_upper": record.get("positive_covariance_error_upper"),
        "negative_covariance_error_upper": record.get("negative_covariance_error_upper"),
        "entered_as_zero": False,
        "covariance_bound_used_in_aggregate": False,
        "retuned": False,
        "method_too_loose": record.get("method_too_loose"),
        "physical_local_gate": "OPEN",
    }


def aggregate(catalogue, completed, failures, *, mode, row_budget, cpu_budget,
              cpu_seconds, capture_cpu_seconds, new_rows_captured, new_rows_failed,
              budget_exhausted):
    by_key = {(record["panel"], record["row"]): record for record in completed}
    failed_key = {(record["panel"], record["row"]): record for record in failures}
    if len(by_key) != len(completed) or len(failed_key) != len(failures):
        raise ValueError("duplicate source identity")
    if set(by_key) & set(failed_key):
        raise ValueError("a row cannot be both witnessed and failed")
    catalogue_keys = {(spec.panel, spec.row) for spec in catalogue}
    if any(key not in catalogue_keys for key in set(by_key) | set(failed_key)):
        raise ValueError("checkpoint identity is not an uncovered subgap row")
    ordered = [by_key[(spec.panel, spec.row)] for spec in catalogue
               if (spec.panel, spec.row) in by_key]
    failed = [failed_key[(spec.panel, spec.row)] for spec in catalogue
              if (spec.panel, spec.row) in failed_key]
    missing = [{"panel": spec.panel, "row": spec.row, "energy_hex": _hex(spec.energy)}
               for spec in catalogue
               if (spec.panel, spec.row) not in by_key and (spec.panel, spec.row) not in failed_key]
    if len(ordered) + len(failed) + len(missing) != len(catalogue):
        raise ValueError("subgap census does not add up")
    if any(item.get("entered_as_zero") for item in failed):
        raise ValueError("a failed row must not be entered as zero")
    complete = len(ordered) == len(catalogue) and not missing and not failed
    return {
        "schema": SCHEMA,
        "status": WINDOW_STATUS if complete else OPEN_STATUS,
        "group": GROUP,
        "angular_sign": ANGULAR_SIGN,
        "panels": [panel for panel, _rows in PANEL_ROWS],
        "expected_positive_rows": EXPECTED_POSITIVE_ROWS,
        "expected_negative_partners": EXPECTED_POSITIVE_ROWS,
        "expected_signed_rows": EXPECTED_SIGNED_ROWS,
        "completed_positive_rows": len(ordered),
        "completed_negative_partners": len(ordered),
        "completed_signed_rows": 2 * len(ordered),
        "failed_positive_rows": len(failed),
        "failed_signed_rows_not_bounded": 2 * len(failed),
        "missing_positive_rows": len(missing),
        "coverage_complete": complete,
        "interrupted": not complete,
        "finite_window_closes_gate": False,
        "failed_rows_entered_as_zero": False,
        "missing_rows_entered_as_zero": False,
        "mode": mode,
        "rows": [_summary(record) for record in ordered],
        "failed": [_failure_summary(record) for record in failed],
        "missing": missing,
        "covariance_errors_unweighted": True,
        "quadrature_weight_applied": False,
        "negative_bounds_supplied_separately": True,
        "source_columns_replaced": False,
        "all_source_families": False,
        "other_positive_angular_family_included": False,
        "physical_upstream_budget_component": None,
        "physical_local_gate": "OPEN",
        "row_budget": row_budget,
        "cpu_budget_seconds": cpu_budget,
        "new_rows_captured": new_rows_captured,
        "new_rows_failed": new_rows_failed,
        "budget_exhausted": budget_exhausted,
        "cpu_seconds": cpu_seconds,
        "capture_cpu_seconds": capture_cpu_seconds,
        "replay_uses_saved_phase_and_bloch_witness": True,
        "preparation_jost_rerun_on_replay": False,
    }


def _authenticate_saved(output_root, catalogue, dependencies):
    by_key = {(spec.panel, spec.row): spec for spec in catalogue}
    authenticated = []
    for dirname, record, payload in load_checkpoints(output_root):
        panel, row = _row_key(record)
        spec = by_key.get((panel, row))
        if spec is None:
            raise ValueError("checkpoint identity is not an uncovered subgap row")
        authenticated.append(authenticate_checkpoint(
            spec, dirname, record, payload, dependencies))
    failed = []
    for dirname, record in load_failures(output_root):
        panel, row = _row_key(record)
        spec = by_key.get((panel, row))
        if spec is None:
            raise ValueError("checkpoint identity is not an uncovered subgap row")
        failed.append(_authenticate_failure(spec, dirname, record, dependencies))
    return authenticated, failed


def _resolve_output(archive_root, output_root):
    output_root = Path(output_root)
    if not output_root.is_absolute():
        output_root = Path(archive_root) / output_root
    return output_root.resolve(strict=False)


@contextmanager
def _cpu_limit(seconds):
    """Interrupt the existing solver when the remaining capture budget is spent.

    ``RLIMIT_CPU`` is whole seconds. Only the soft limit is lowered. This
    platform refuses to raise a hard limit again, so lowering it would make
    the next row unable to restore the budget.
    """
    if not math.isfinite(seconds) or seconds <= 0:
        raise TimeoutError("subgap source-window CPU budget exhausted")
    usage = resource.getrusage(resource.RUSAGE_SELF)
    consumed = usage.ru_utime + usage.ru_stime
    soft = max(1, int(math.ceil(consumed + float(seconds))))
    old_limit = resource.getrlimit(resource.RLIMIT_CPU)
    hard = old_limit[1]
    if hard != resource.RLIM_INFINITY and soft >= hard:
        raise TimeoutError("subgap source-window CPU budget exhausted")
    old_handler = signal.getsignal(signal.SIGXCPU)

    def handler(signum, frame):
        raise TimeoutError("subgap source-window CPU budget exhausted")

    signal.signal(signal.SIGXCPU, handler)
    try:
        resource.setrlimit(resource.RLIMIT_CPU, (soft, hard))
        try:
            yield
        finally:
            resource.setrlimit(resource.RLIMIT_CPU, old_limit)
    finally:
        signal.signal(signal.SIGXCPU, old_handler)


def _emit(event):
    print(json.dumps(event, sort_keys=True), file=sys.stderr, flush=True)


def cover(archive_root, output_root, *, mode, row_budget=None, cpu_budget=None,
          emit=None):
    """Validate saved witnesses before any new capture. Check mode never solves."""
    if mode not in {"resume", "check"}:
        raise ValueError("coverage mode must be resume or check")
    if mode == "check":
        if row_budget is not None or cpu_budget is not None:
            raise ValueError("check replays saved witnesses and does not accept a solve budget")
    else:
        row_budget, cpu_budget = _budgets(row_budget, cpu_budget)
    started = time.process_time()
    archive_root = Path(archive_root).resolve(strict=True)
    output_root = _resolve_output(archive_root, output_root)
    _refuse_foreign_campaign(archive_root, output_root)
    _require_root_names(output_root)
    archive = RetainedUpstreamArchive(archive_root)
    # Census raises before a freeze file or a solve.
    catalogue = subgap_catalogue(archive)
    dependencies = dependency_identity(archive_root, archive)
    if mode == "check" and not output_root.exists():
        return aggregate(
            catalogue, [], [], mode=mode, row_budget=None, cpu_budget=None,
            cpu_seconds=time.process_time() - started, capture_cpu_seconds=0.0,
            new_rows_captured=0, new_rows_failed=0, budget_exhausted=False)
    if mode == "resume":
        output_root = _canonical_directory(output_root, create=True)
        _require_root_names(output_root)
        ensure_proof_freeze(archive_root, output_root)
    elif (output_root / "proof-freeze.json").exists():
        ensure_proof_freeze(archive_root, output_root)
    output_root = _canonical_directory(output_root, create=False)
    authenticated, failed = _authenticate_saved(output_root, catalogue, dependencies)
    done = {(record["panel"], record["row"]) for record in authenticated}
    failed_keys = {(record["panel"], record["row"]) for record in failed}
    captured = 0
    failed_now = 0
    budget_exhausted = False
    capture_started = time.process_time()
    reporter = emit or _emit
    if mode == "resume":
        for spec in catalogue:
            key = (spec.panel, spec.row)
            if key in done or key in failed_keys:
                continue
            if captured + failed_now >= row_budget:
                break
            remaining = cpu_budget - (time.process_time() - capture_started)
            if remaining <= 0:
                budget_exhausted = True
                break
            try:
                with _cpu_limit(remaining):
                    record, payload = capture_row(spec, dependencies, cpu_limit=remaining)
            except TimeoutError:
                budget_exhausted = True
                reporter({
                    "event": "budget_exhausted", "panel": spec.panel, "row": spec.row,
                    "published": False,
                    "capture_cpu_seconds": time.process_time() - capture_started,
                })
                break
            except EstablishedControlRejected as exc:
                failure = _failure_record(spec, dependencies, exc)
                publish_failure(output_root, failure)
                failed.append(failure)
                failed_keys.add(key)
                failed_now += 1
                reporter({
                    "event": "failed", "panel": spec.panel, "row": spec.row,
                    "exception_type": type(exc).__name__,
                    "exception_message": str(exc),
                    "entered_as_zero": False,
                    "capture_cpu_seconds": time.process_time() - capture_started,
                })
                continue
            publish_checkpoint(output_root, record, payload)
            authenticated.append(record)
            done.add(key)
            captured += 1
            reporter({
                "event": "captured", "panel": spec.panel, "row": spec.row,
                "profile": record["transport_profile"],
                "phase_cells": record["phase_cells"],
                "capture_cpu_seconds": time.process_time() - capture_started,
            })
    return aggregate(
        catalogue, authenticated, failed, mode=mode, row_budget=row_budget,
        cpu_budget=cpu_budget, cpu_seconds=time.process_time() - started,
        capture_cpu_seconds=(time.process_time() - capture_started) if mode == "resume" else 0.0,
        new_rows_captured=captured, new_rows_failed=failed_now,
        budget_exhausted=budget_exhausted)
