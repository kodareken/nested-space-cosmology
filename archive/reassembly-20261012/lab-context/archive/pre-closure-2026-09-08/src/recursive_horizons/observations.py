"""Compact validation of published GW170817/GRB 170817A timing inputs.

This dependency-free module validates four small public metadata products and
reconstructs the collaboration's *published input arithmetic*.  It is not a
strain or gamma-ray light-curve reanalysis, and it cannot separate an intrinsic
source emission lag from a propagation delay.
"""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
from math import isfinite, sqrt
from numbers import Real
from pathlib import Path
import re
from typing import Any, Final


SPEED_OF_LIGHT_M_S: Final[float] = 299_792_458.0
MPC_M: Final[float] = 3.085_677_581_491_367_3e22

_LVC_UTC_RE = re.compile(r"TRIGGER_TIME:\s*[^\{]+\{(?P<utc>\d{2}:\d{2}:\d{2}\.\d{6})\}\s*UT")
_GCN_UTC_RE = re.compile(r"At (?P<utc>\d{2}:\d{2}:\d{2}\.\d{2}) UT")
_GCN_TRIGGER_RE = re.compile(r"trigger (?P<met>\d+) / (?P<name>\d+)")


class ObservationalInputError(ValueError):
    """Raised when a compact local observational input is malformed or inconsistent."""


def _repository_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _finite(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ObservationalInputError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ObservationalInputError(f"{name} must be finite")
    return result


def _positive(name: str, value: Real) -> float:
    result = _finite(name, value)
    if result <= 0.0:
        raise ObservationalInputError(f"{name} must be positive")
    return result


def load_published_constraints(path: Path | None = None) -> dict[str, Any]:
    """Load the immutable, compact published-constraint record."""

    source = path or _repository_root() / "collected-data" / "published-constraints.json"
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ObservationalInputError(f"cannot read published constraints: {source}") from error
    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ObservationalInputError("published constraints must be schema version 1")
    if not isinstance(payload.get("GW170817_GRB170817A"), dict):
        raise ObservationalInputError("published constraints lack GW170817_GRB170817A")
    return payload


def file_sha256(path: Path) -> str:
    """Return a file digest without treating a matching hash as scientific validation."""

    try:
        return sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ObservationalInputError(f"cannot read input file: {path}") from error


def parse_fits_primary_header(path: Path) -> dict[str, Any]:
    """Parse standard 80-byte FITS primary-header cards without third-party code."""

    try:
        file_bytes = path.read_bytes()
    except OSError as error:
        raise ObservationalInputError(f"cannot read FITS input: {path}") from error
    if len(file_bytes) < 80:
        raise ObservationalInputError("FITS primary header is incomplete")

    cards: dict[str, Any] = {}
    saw_end = False
    # FITS headers occupy one or more 2880-byte blocks.  The local GBM file's
    # primary header spans two blocks, so stopping after the first would omit END.
    for start in range(0, len(file_bytes) - 79, 80):
        card = file_bytes[start : start + 80].decode("ascii", errors="strict")
        keyword = card[:8].strip()
        if keyword == "END":
            saw_end = True
            break
        if not keyword or card[8:10] != "= ":
            continue
        raw_value = card[10:80].strip()
        cards[keyword] = _parse_fits_value(raw_value)
    if not saw_end:
        raise ObservationalInputError("FITS primary header has no END card")
    return cards


def _parse_fits_value(raw_value: str) -> Any:
    if raw_value.startswith("'"):
        closing = raw_value.find("'", 1)
        if closing < 0:
            raise ObservationalInputError("unterminated FITS string card")
        return raw_value[1:closing].rstrip()
    raw_value = raw_value.split("/", 1)[0].strip()
    if raw_value == "T":
        return True
    if raw_value == "F":
        return False
    normalized = raw_value.replace("D", "E")
    try:
        return float(normalized) if any(mark in normalized for mark in ".Ee") else int(normalized)
    except ValueError as error:
        raise ObservationalInputError(f"unsupported FITS card value: {raw_value!r}") from error


def parse_lvc_notice_trigger_utc(text: str) -> str:
    """Return the repeated UTC trigger time from an LVC notice."""

    matches = _LVC_UTC_RE.findall(text)
    if not matches or len(set(matches)) != 1:
        raise ObservationalInputError("LVC notice must contain one consistent trigger UTC")
    return matches[0]


def parse_gcn_21520_trigger(text: str) -> dict[str, str]:
    """Extract the rounded UTC and GBM trigger identifiers from GCN 21520."""

    utc = _GCN_UTC_RE.search(text)
    trigger = _GCN_TRIGGER_RE.search(text)
    if utc is None or trigger is None:
        raise ObservationalInputError("GCN 21520 lacks the Fermi UTC or trigger identifiers")
    return {"utc": utc.group("utc"), "met_integer": trigger.group("met"), "trigger_name": trigger.group("name")}


def _utc_datetime(utc_timestamp: str) -> datetime:
    try:
        stamp = datetime.fromisoformat(utc_timestamp.replace("Z", "+00:00"))
    except ValueError as error:
        raise ObservationalInputError("UTC timestamp must be ISO-8601") from error
    if stamp.tzinfo != timezone.utc:
        raise ObservationalInputError("UTC timestamp must explicitly use UTC")
    return stamp


def utc_seconds(utc_timestamp: str) -> float:
    """Return seconds after midnight from an ISO UTC timestamp ending in ``Z``."""

    stamp = _utc_datetime(utc_timestamp)
    return stamp.hour * 3600.0 + stamp.minute * 60.0 + stamp.second + stamp.microsecond / 1_000_000.0


def reconstructed_geocentric_delay_seconds(
    fermi_trigger_utc: str,
    gw_geocentric_utc: str,
    onset_before_trigger_seconds: Real,
    fermi_arrival_before_geocenter_seconds: Real,
) -> float:
    """Reconstruct ``t_gamma,onset,geo - t_GW,merger,geo`` from published inputs."""

    onset = _positive("onset_before_trigger_seconds", onset_before_trigger_seconds)
    fermi_to_geocenter = _positive(
        "fermi_arrival_before_geocenter_seconds", fermi_arrival_before_geocenter_seconds
    )
    trigger = _utc_datetime(fermi_trigger_utc)
    merger = _utc_datetime(gw_geocentric_utc)
    result = (trigger - merger).total_seconds() - onset + fermi_to_geocenter
    if not isfinite(result):
        raise ObservationalInputError("reconstructed geocentric delay is not finite")
    return result


def combined_timing_uncertainty_seconds(
    gw_uncertainty_seconds: Real, gamma_onset_uncertainty_seconds: Real
) -> float:
    """Combine the declared independent published timing-input uncertainties."""

    gw_uncertainty = _positive("gw_uncertainty_seconds", gw_uncertainty_seconds)
    onset_uncertainty = _positive("gamma_onset_uncertainty_seconds", gamma_onset_uncertainty_seconds)
    return sqrt(gw_uncertainty**2 + onset_uncertainty**2)


def propagation_delta_from_gamma_minus_gw_delay(
    gamma_minus_gw_delay_seconds: Real, distance_mpc: Real
) -> float:
    """Solve ``delta=(v_GW-v_EM)/v_EM`` for a stated propagation-only delay.

    The input delay is ``t_gamma-t_GW`` after an explicitly chosen source-lag
    model.  It is not an inference from the observed lag alone.
    """

    delay = _finite("gamma_minus_gw_delay_seconds", gamma_minus_gw_delay_seconds)
    distance = _positive("distance_mpc", distance_mpc)
    travel_time = distance * MPC_M / SPEED_OF_LIGHT_M_S
    if not isfinite(travel_time) or travel_time <= delay:
        raise ObservationalInputError("distance/delay combination has no positive photon travel time")
    result = delay / (travel_time - delay)
    if not isfinite(result) or result <= -1.0:
        raise ObservationalInputError("derived relative speed shift is outside the physical range")
    return result


def published_speed_bounds(constraint: dict[str, Any]) -> dict[str, float]:
    """Reconstruct published conservative speed-bound arithmetic from its assumptions."""

    timing = _require_mapping(constraint, "published_timing")
    speed = _require_mapping(constraint, "published_speed_constraint")
    observed_delay = _finite("published_geocentric_delay_seconds", timing["published_geocentric_delay_seconds"])
    distance = _positive("conservative_distance_mpc", speed["conservative_distance_mpc"])
    source_lag = _positive("source_lag_for_lower_bound_seconds", speed["source_lag_for_lower_bound_seconds"])
    return {
        "lower": propagation_delta_from_gamma_minus_gw_delay(observed_delay - source_lag, distance),
        "upper": propagation_delta_from_gamma_minus_gw_delay(observed_delay, distance),
    }


def _require_mapping(payload: dict[str, Any], key: str) -> dict[str, Any]:
    value = payload.get(key)
    if not isinstance(value, dict):
        raise ObservationalInputError(f"missing mapping: {key}")
    return value


def _rounded_one_significant(value: float) -> float:
    """Round a nonzero finite value to one significant decimal digit."""

    if not isfinite(value) or value == 0.0:
        raise ObservationalInputError("cannot round a zero or non-finite speed bound")
    exponent = int(f"{abs(value):e}".split("e")[1])
    return round(value, -exponent)


def validate_gw170817_inputs(
    data_directory: Path | None = None, constraints_path: Path | None = None
) -> dict[str, Any]:
    """Validate compact source provenance and reconstruct published timing arithmetic.

    A successful result confirms only that the frozen small source files agree
    with the declared published inputs and arithmetic.  It does not reanalyse
    raw strain/light curves or remove the intrinsic source-lag assumption.
    """

    constraints = load_published_constraints(constraints_path)
    record = constraints["GW170817_GRB170817A"]
    if not isinstance(record, dict):  # guarded by load_published_constraints
        raise ObservationalInputError("GW170817 record must be an object")
    source_dir = data_directory or _repository_root() / "collected-data" / "gw170817"
    sources = _require_mapping(record, "local_sources")
    for filename, metadata in sources.items():
        if not isinstance(metadata, dict) or not isinstance(metadata.get("sha256"), str):
            raise ObservationalInputError(f"source metadata malformed for {filename}")
        if file_sha256(source_dir / filename) != metadata["sha256"]:
            raise ObservationalInputError(f"source hash mismatch: {filename}")

    gwosc = json.loads((source_dir / "gwosc-event-v1.json").read_text(encoding="utf-8"))
    if gwosc.get("name") != "GW170817" or gwosc.get("grace_id") != "G298048":
        raise ObservationalInputError("GWOSC metadata does not identify GW170817/G298048")
    if gwosc.get("detectors") != ["H1", "L1", "V1"]:
        raise ObservationalInputError("GWOSC detector list differs from frozen input")
    if abs(_finite("GWOSC gps", gwosc.get("gps")) - 1187008882.4) > 1.0e-6:
        raise ObservationalInputError("GWOSC GPS value differs from frozen event metadata")

    fits = parse_fits_primary_header(source_dir / "glg_tcat_all_bn170817529_v03.fit")
    if fits.get("OBJECT") != "GRB170817529" or fits.get("DATE-OBS") != "2017-08-17T12:38:54":
        raise ObservationalInputError("GBM FITS header does not identify the expected event/observation")
    timing = _require_mapping(record, "published_timing")
    if abs(_finite("TRIGTIME", fits.get("TRIGTIME")) - _finite("fermi_trigger_met_seconds", timing["fermi_trigger_met_seconds"])) > 1.0e-6:
        raise ObservationalInputError("GBM FITS TRIGTIME differs from published constraint record")

    lvc_utc = parse_lvc_notice_trigger_utc((source_dir / "G298048.lvc").read_text(encoding="utf-8"))
    if lvc_utc != "12:41:04.445710":
        raise ObservationalInputError("LVC notice trigger UTC differs from frozen notice")
    gcn = parse_gcn_21520_trigger((source_dir / "GCN-21520.txt").read_text(encoding="utf-8"))
    if gcn != {"utc": "12:41:06.47", "met_integer": "524666471", "trigger_name": "170817529"}:
        raise ObservationalInputError("GCN 21520 trigger fields differ from frozen circular")
    if not str(timing["fermi_trigger_utc"]).startswith("2017-08-17T" + gcn["utc"]):
        raise ObservationalInputError("precise Fermi trigger UTC does not agree with GCN rounding")

    reconstructed_delay = reconstructed_geocentric_delay_seconds(
        str(timing["fermi_trigger_utc"]),
        str(timing["gw_geocentric_utc"]),
        timing["gamma_onset_before_trigger_seconds"],
        timing["fermi_arrival_before_geocenter_seconds"],
    )
    uncertainty = combined_timing_uncertainty_seconds(
        timing["gw_geocentric_uncertainty_seconds"], timing["gamma_onset_uncertainty_seconds"]
    )
    bounds = published_speed_bounds(record)
    declared_interval = _require_mapping(record, "published_speed_constraint")["delta_v_over_v_em_interval"]
    if not isinstance(declared_interval, list) or len(declared_interval) != 2:
        raise ObservationalInputError("declared speed interval must have two bounds")
    if [_rounded_one_significant(bounds["lower"]), _rounded_one_significant(bounds["upper"])] != declared_interval:
        raise ObservationalInputError("reconstructed speed bounds do not round to the published interval")
    return {
        "classification": "published_timing_input_reconstruction_not_raw_reanalysis",
        "gwosc_gps_seconds": gwosc["gps"],
        "lvc_notice_utc": lvc_utc,
        "gbm_object": fits["OBJECT"],
        "gbm_trigger_met_seconds": fits["TRIGTIME"],
        "gcn_trigger_utc_rounded": gcn["utc"],
        "reconstructed_geocentric_delay_seconds": reconstructed_delay,
        "combined_input_uncertainty_seconds": uncertainty,
        "published_geocentric_delay_seconds": timing["published_geocentric_delay_seconds"],
        "published_geocentric_delay_uncertainty_seconds": timing["published_geocentric_delay_uncertainty_seconds"],
        "reconstructed_delta_v_over_v_em": bounds,
        "published_rounded_delta_v_over_v_em_interval": declared_interval,
        "limitation": "published timing-input reconstruction only; intrinsic source emission lag remains an assumption",
    }


def reproduce_gw170817_timing(
    data_directory: Path | None = None, constraints_path: Path | None = None
) -> dict[str, Any]:
    """Return the versioned, JSON-safe compact timing-reconstruction artifact."""

    source = constraints_path or _repository_root() / "collected-data" / "published-constraints.json"
    constraints = load_published_constraints(source)
    event = _require_mapping(constraints, "GW170817_GRB170817A")
    primary = _require_mapping(event, "primary_paper")
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "artifact": "GW170817_GRB170817A_published_timing_input_reconstruction",
        "classification": "published_timing_input_reconstruction_not_raw_reanalysis",
        "primary_paper": primary,
        "published_constraints": {
            "filename": source.name,
            "sha256": file_sha256(source),
        },
        "validation": validate_gw170817_inputs(data_directory, source),
        "nonclaims": [
            "Not a raw gravitational-wave strain or gamma-ray light-curve reanalysis.",
            "Does not remove the intrinsic source-emission-lag assumption.",
            "Does not measure a raw dimensionful local c or an inaccessible parent-domain speed.",
            "Does not establish an external origin for dark energy or a black-hole-to-child transition.",
        ],
    }
