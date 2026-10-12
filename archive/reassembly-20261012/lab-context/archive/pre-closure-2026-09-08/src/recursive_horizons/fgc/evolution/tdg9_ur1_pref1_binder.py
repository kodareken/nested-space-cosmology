"""Independent post-result binder for the completed TDG9 UR1 terminal.

The one-time live path authenticates the immutable UR1 authority commit and
every committed delta blob, reads the canonical two-leaf raw namespace with
no-following race-checked opens, independently decodes the published ``u:R``
envelope fields, and reclassifies D01/D12 through the design-only TDG6
classifier.  It deliberately does not import UR1 or AC1 runners or
authorities, TI2 runner or binder code, the runtime envelope implementation,
or LOC1 row reconstruction.  Ordinary verification consumes only the tracked
compact certificate and is raw/store/shadow blind.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import stat
import subprocess
import tomllib
from typing import Any, Mapping, NoReturn

from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    classify_tdg6_channel,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG9-UR1-PREF1"
CLASSIFICATION = (
    "independently_bound_retry3_u_R_envelope_lower_owner_terminal"
)
CONFIG_PATH = "configs/fgc/fgc-1-tdg9-ur1-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg9-ur1-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg9-ur1-pref1.md"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"

AUTHORITY_COMMIT = "17d8c4cf728ffd18ce6e137874501e48b35fbb32"
AUTHORITY_PARENT = "22f7dc5e181af45eb43f6f19dee23fbc29f89c0c"
RAW_NAMESPACE = "runs/fgc-2-sf1/tdg9-ur1/retry3-ur-envelope-owner"
RAW_SCHEMA = "UR1-raw-v1"
RAW_CLASSIFICATION = (
    "completed_u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity"
)
RUNNER_ID = "FGC-1-TDG9-UR1-RUN1"
UR1_ARTIFACT_ID = "FGC-1-TDG9-UR1-FRZ1"
RAW_MANIFEST_SHA256 = (
    "2b82a93a171eea852fb48770924a7a28b812c120cc5af5316d2a29a03dbd7db9"
)
RAW_TERMINAL_SHA256 = (
    "545a163a650810ecaceb11bbd26fd9e8adb3559635e26e8dad35d88460a3f06f"
)
COMPACT_RESULT_SHA256 = (
    "27c2aef9b033092285b868f583441842de99f61941ff57aa7e62749cc8d91300"
)

STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = (
    "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
)
AC1_PREF1_RESULT_SHA256 = (
    "173183871b2d750f98e9758cdc5f702be40b0bbb73d10e707672e0e68752935d"
)
EXPERIMENT_LABEL = (
    "retry3_SSPRK3_on_inherited_SBP4_u_R_row_hash_and_"
    "zero_lower_envelope_owner_diagnostic"
)
TABLEAU_RUNTIME_SELECTOR = "second_order_diagonal_norm_SBP_plus_SSPRK3"
MEMBER_KEY = "RK4-2049"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
RETRY = 3
WIDTH_HEX = "0x1.aaa9612df8000p-11"
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
TRANSACTION_SHA256 = (
    "7a90163d2eb4252fa7a1bbdf55d12f92a9127ed7cba5a37c28456e72d732b22b"
)
HISTORICAL_JOURNAL_SHA256 = (
    "0d21f650ccc46db778fd81fbebf5acab6394e1c2f7972b7940e24c76af9e1f59"
)
PUBLISHED_CHANNEL = "u:R"
ROW_HASH_DOMAIN = "TDG9-AR1-BINARY64-HERMITE-ROWS-v1"
TI2_RETRY3_U_R_ROW_SHA256 = (
    "b4a48c7bf72e015f3f3d9a340ea550d608514fb8cd886562265158c33c0b3969"
)
TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS = "sufficient_contraction_pass"
AC1_U_R_D01_LOWER_HEX = "0x0.0p+0"
AC1_U_R_D01_UPPER_HEX = "0x1.2d198e246e459p-38"
AC1_U_R_D12_LOWER_HEX = "0x0.0p+0"
AC1_U_R_D12_UPPER_HEX = "0x1.2d1c31bb91376p-38"
AC1_U_R_CLASS = "order_inconclusive"
OWNER_NAME = (
    "coefficient_plus_arithmetic_debit_clips_positive_raw_maximum"
)
ZERO_LOWER_OWNERS = (
    "raw_candidate_maximum_is_zero",
    OWNER_NAME,
    "downward_binary64_rounding_of_positive_exact_lower",
)
EXACT_LOWER_FORMULA = (
    "max(0,Fraction(raw_candidate_maximum)-"
    "Fraction(coefficient_construction_debit)-"
    "Fraction(outward_arithmetic_debit))"
)
D01_CONSTRUCTION_OVER_RAW = Fraction(654019454794891, 25865933815808)
D12_CONSTRUCTION_OVER_RAW = Fraction(163504863789045, 5012309316608)
WORK_BUDGET = {
    "SSPRK3_records_per_proposal": 4,
    "maximum_stage_and_endpoint_RHS_records": 28,
    "published_channel": PUBLISHED_CHANNEL,
    "published_channel_count": 1,
    "published_row_count": OWNED_ROW_COUNT,
    "resource_escalation_authorized": False,
    "retry": 3,
    "retry_4_authorized": False,
    "retry_5_authorized": False,
    "retry_count": 1,
    "fourth_width_authorized": False,
    "shadow_paths": 7,
    "shadow_proposals": 7,
}
EXPECTED_REPLAY_RECEIPT = {
    "accepted_time_hex": ACCEPTED_TIME_HEX,
    "attempted_width_hex": WIDTH_HEX,
    "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
    "historical_journal_sha256": HISTORICAL_JOURNAL_SHA256,
    "member_key": MEMBER_KEY,
    "retry": RETRY,
    "state_sha256": PHYSICAL_STATE_SHA256,
    "transaction_sha256": TRANSACTION_SHA256,
}
_EXACT_KEYS = frozenset({"binary64_hex", "numerator", "denominator"})
_RATIONAL_KEYS = frozenset({"numerator", "denominator"})
_ENVELOPE_KEYS = frozenset(
    {
        "ambiguous_discriminant_count",
        "bernstein_certification_slack",
        "certified_continuous_upper_bound",
        "coefficient_construction_debit",
        "endpoint_candidate_count",
        "outward_arithmetic_debit",
        "polynomial_count",
        "raw_candidate_maximum",
        "real_interior_root_count",
    }
)
_OWNER_KEYS = frozenset(
    {
        "bernstein_certification_slack",
        "bernstein_certification_slack_not_part_of_lower_clip",
        "coefficient_construction_debit",
        "exact_clipped",
        "exact_unclipped",
        "outward_arithmetic_debit",
        "owner",
        "raw_candidate_maximum",
        "reproduced_stored_lower",
        "source",
    }
)
_INTERVAL_KEYS = frozenset(
    {
        "envelope",
        "lower_bound",
        "owned_row_count",
        "subinterval_count",
        "upper_bound",
        "zero_lower_owner",
    }
)
_U_R_KEYS = frozenset(
    {
        "D01",
        "D01_D12_lower_owners_equal",
        "D12",
        "admission_passed",
        "channel",
        "classification",
        "matches_sealed_AC1_production_intervals",
        "order_threshold_passed",
        "order_threshold_resolved",
        "temporal_retry_permitted",
    }
)
_ROW_BASE_KEYS = frozenset(
    {
        "algorithm",
        "domain",
        "expected_sha256",
        "matched",
        "observed_sha256",
        "row_count",
    }
)
_ROW_MATCH_KEYS = _ROW_BASE_KEYS | {
    "related_TI2_retry3_u_R_radius_free_complete_C_class",
    "related_not_replacement_production_admission",
}
_DEBIT_SOURCE_KEYS = frozenset(
    {
        "coefficient_construction_debit",
        "outward_arithmetic_debit",
        "raw_candidate_maximum",
        "term",
    }
)
_COMPLETED_EXTRA_KEYS = frozenset(
    {
        "SSPRK3_stage_and_endpoint_record_count",
        "replay_receipt",
        "row_stream",
        "shadow_path_count",
        "shadow_proposal_count",
        "u_R",
    }
)
_MAX_RAW_LEAF_BYTES = 8 * 1024 * 1024
_MAX_STORE_LEAF_BYTES = 128 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")
_D01_HEX = {
    "arithmetic": "0x1.6b2daa9000016p-90",
    "bernstein": "0x0.0p+0",
    "certified": AC1_U_R_D01_UPPER_HEX,
    "construction": "0x1.2969e2a3be458p-38",
    "lower": AC1_U_R_D01_LOWER_HEX,
    "raw": "0x1.78661cc000000p-43",
    "upper": AC1_U_R_D01_UPPER_HEX,
}
_D12_HEX = {
    "arithmetic": "0x1.b9799e940001bp-91",
    "bernstein": "0x0.0p+0",
    "certified": AC1_U_R_D12_UPPER_HEX,
    "construction": "0x1.2969e2a67fea0p-38",
    "lower": AC1_U_R_D12_LOWER_HEX,
    "raw": "0x1.23c13aa500000p-43",
    "upper": AC1_U_R_D12_UPPER_HEX,
}
_D01_COUNTS = {
    "ambiguous_discriminant_count": 0,
    "endpoint_candidate_count": 8176,
    "owned_row_count": OWNED_ROW_COUNT,
    "polynomial_count": 4088,
    "real_interior_root_count": 2023,
    "subinterval_count": 2,
}
_D12_COUNTS = {
    "ambiguous_discriminant_count": 0,
    "endpoint_candidate_count": 16352,
    "owned_row_count": OWNED_ROW_COUNT,
    "polynomial_count": 8176,
    "real_interior_root_count": 4355,
    "subinterval_count": 4,
}
_AUTHORITY_BLOBS = (
    (
        "Makefile",
        "207174051d6a3a3ce70f19efc7c122b452ba71761906924d283525badcea4e64",
    ),
    (
        "README.md",
        "73ff82f9c064a0036e912e880c3825783ad77113021e0d749db04cb03477b50f",
    ),
    (
        "configs/fgc/fgc-1-tdg9-ur1-frz1.toml",
        "6aaa5ae0b91753ac86e36235629ab4ed00ddfe6e4c83b69ad2594d755dcbc303",
    ),
    (
        "docs/claim-ledger.md",
        "40861eaef270d47aa93f6ba9b0b610deabcaa25de8134dfde71941b5defcb085",
    ),
    (
        "docs/fgc-runtime-matrix.md",
        "124711908ef0a2b1b3ca54aeec72960d423c5002419161a8523ee04157d58f61",
    ),
    (
        "docs/fgc-tdg9-ur1-frz1.md",
        "aba34cd72182f795c969397c18acbd5b680b035711e58ea9ff1ff1651a354d2e",
    ),
    (
        "docs/research-roadmap.md",
        "daa2b02d2401057f521d3cf78a9f0e2d2ed2cccb94c180556f38e5dad1dff325",
    ),
    (
        "results/README.md",
        "ce620016d40dd1f5908e9b5f8ceb6d6e2630e20efae36f3a5367a59a733a1bbb",
    ),
    (
        "results/fgc-1-tdg9-ur1-frz1.json",
        "6dad298d01e98080bf3dcd7bde3e4aa969efdda790d97f8d68b40231cee9283f",
    ),
    (
        "scripts/check_repo.py",
        "cec3bf1d79a065a7e6e3c4254f3c4847eb39b573663ebf24099de479cf6df191",
    ),
    (
        "scripts/reproduce_fgc_tdg9_ur1_frz1.py",
        "f67c5eb612f92640e8555117f7026a48d3c7c2193e2cc8feed54a56ed155040f",
    ),
    (
        "scripts/run_fgc_tdg9_ur1.py",
        "f660ba40df882d6ded87d8b0c020f8328fa50fa26a0a61f7f5686febe32c5680",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_ur1_authority.py",
        "7ab52c8e910ae0b7bf9ceea415e9e12eedcf09afc67c6986cfc7588d882e015f",
    ),
    (
        "tests/test_check_repo_tdg9_ur1_frz1.py",
        "afae3ce3d06bf3296a5d86101860924083bd81612cec7aa429b6ae10ca9363ca",
    ),
    (
        "tests/test_fgc_tdg9_ur1_authority.py",
        "7e74735f07f27da97a42ff76e9d65eb672d14c411cfce160096d2392ab6e9336",
    ),
    (
        "tests/test_fgc_tdg9_ur1_runner.py",
        "931b677b48c8c8e91ae0d025dd0a8dba37d09e88507cdd48a8a34c6154badd22",
    ),
)


class TDG9UR1PREF1Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG9UR1PREF1Error(stop_id, detail)


def canonical_result(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                sort_keys=True,
                indent=2,
                ensure_ascii=True,
                allow_nan=False,
            )
            + "\n"
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG9UR1PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError(f"duplicate key {key}")
        answer[key] = value
    return answer


def _json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(
            raw.decode("ascii"),
            object_pairs_hook=_pairs,
            parse_constant=lambda item: (_ for _ in ()).throw(ValueError(item)),
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG9UR1PREF1Error("PREF1_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF1_JSON_DRIFT", f"{label} is not canonical pretty JSON")
    return value


def _git(root: Path, *arguments: str) -> bytes:
    environment = dict(
        os.environ,
        LC_ALL="C",
        LANG="C",
        GIT_NO_REPLACE_OBJECTS="1",
        GIT_CONFIG_NOSYSTEM="1",
    )
    try:
        return subprocess.run(
            (
                "git",
                "--no-replace-objects",
                "--no-optional-locks",
                "-c",
                "core.fsmonitor=false",
                "-c",
                "core.untrackedCache=false",
                *arguments,
            ),
            cwd=root,
            env=environment,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except subprocess.CalledProcessError as exc:
        raise TDG9UR1PREF1Error("PREF1_GIT_DRIFT", arguments) from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
    )


def _open_directory(root: Path, relative: str) -> int:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _stop("PREF1_PATH_UNSAFE", relative)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in candidate.parts:
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            if not stat.S_ISDIR(os.fstat(child).st_mode):
                os.close(child)
                _stop("PREF1_PATH_UNSAFE", relative)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _read_leaf_at(
    directory_fd: int, name: str, *, maximum: int, label: str
) -> bytes:
    if not name or "/" in name or name in {".", ".."}:
        _stop("PREF1_PATH_UNSAFE", label)
    before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > maximum
    ):
        _stop("PREF1_LEAF_UNSAFE", label)
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        active = os.fstat(descriptor)
        if _identity(active) != _identity(before):
            _stop("PREF1_LEAF_RACED", label)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _stop("PREF1_LEAF_SHORT_READ", label)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _stop("PREF1_LEAF_GREW", label)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    final = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    if _identity(before) != _identity(after) or _identity(before) != _identity(final):
        _stop("PREF1_LEAF_CHANGED", label)
    return b"".join(chunks)


def _raw_pair(root: Path) -> tuple[bytes, bytes]:
    descriptor = _open_directory(root, RAW_NAMESPACE)
    try:
        before = os.fstat(descriptor)
        with os.scandir(descriptor) as entries:
            names = tuple(sorted(item.name for item in entries))
        if names != ("manifest.json", "terminal.json"):
            _stop("PREF1_RAW_TREE_DRIFT", names)
        manifest_raw = _read_leaf_at(
            descriptor,
            "manifest.json",
            maximum=1 << 20,
            label=f"{RAW_NAMESPACE}/manifest.json",
        )
        terminal_raw = _read_leaf_at(
            descriptor,
            "terminal.json",
            maximum=_MAX_RAW_LEAF_BYTES,
            label=f"{RAW_NAMESPACE}/terminal.json",
        )
        with os.scandir(descriptor) as entries:
            final_names = tuple(sorted(item.name for item in entries))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if names != final_names or _identity(before) != _identity(after):
        _stop("PREF1_RAW_TREE_CHANGED", RAW_NAMESPACE)
    if (
        sha256(manifest_raw).hexdigest() != RAW_MANIFEST_SHA256
        or sha256(terminal_raw).hexdigest() != RAW_TERMINAL_SHA256
    ):
        _stop("PREF1_RAW_HASH_DRIFT", RAW_NAMESPACE)
    return manifest_raw, terminal_raw


def _raw_snapshot(
    root: Path,
) -> tuple[bytes, bytes, dict[str, Any], dict[str, Any]]:
    manifest_raw, terminal_raw = _raw_pair(root)
    return (
        manifest_raw,
        terminal_raw,
        _json(manifest_raw, "raw manifest"),
        _json(terminal_raw, "raw terminal"),
    )


def _snapshot_store(root: Path) -> tuple[int, str]:
    """Reproduce the sealed stack-order store digest with no-following reads."""

    base_fd = _open_directory(root, STORE_PATH)
    digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
    count = 0
    stack: list[tuple[int, str]] = [(base_fd, "")]
    try:
        while stack:
            directory_fd, prefix = stack.pop()
            try:
                before = os.fstat(directory_fd)
                with os.scandir(directory_fd) as entries:
                    items = sorted(entries, key=lambda item: item.name)
                initial_names = tuple(item.name for item in items)
                children: list[tuple[int, str]] = []
                for item in items:
                    metadata = item.stat(follow_symlinks=False)
                    relative = f"{prefix}/{item.name}" if prefix else item.name
                    if item.is_symlink():
                        _stop("PREF1_STORE_SYMLINK", relative)
                    if stat.S_ISDIR(metadata.st_mode):
                        child_fd = os.open(
                            item.name,
                            os.O_RDONLY
                            | getattr(os, "O_DIRECTORY", 0)
                            | getattr(os, "O_NOFOLLOW", 0),
                            dir_fd=directory_fd,
                        )
                        active = os.fstat(child_fd)
                        if _identity(metadata) != _identity(active):
                            os.close(child_fd)
                            _stop("PREF1_STORE_DIRECTORY_RACED", relative)
                        children.append((child_fd, relative))
                    elif stat.S_ISREG(metadata.st_mode):
                        raw = _read_leaf_at(
                            directory_fd,
                            item.name,
                            maximum=_MAX_STORE_LEAF_BYTES,
                            label=relative,
                        )
                        digest.update(
                            (
                                relative
                                + "\0"
                                + sha256(raw).hexdigest()
                                + "\n"
                            ).encode("ascii")
                        )
                        count += 1
                    else:
                        _stop("PREF1_STORE_FOREIGN_TYPE", relative)
                with os.scandir(directory_fd) as entries:
                    final_names = tuple(sorted(item.name for item in entries))
                after = os.fstat(directory_fd)
                if initial_names != final_names or _identity(before) != _identity(
                    after
                ):
                    for child_fd, _ in children:
                        os.close(child_fd)
                    _stop("PREF1_STORE_DIRECTORY_CHANGED", prefix or ".")
                stack.extend(children)
            finally:
                os.close(directory_fd)
    except Exception:
        for directory_fd, _ in stack:
            try:
                os.close(directory_fd)
            except OSError:
                pass
        raise
    return count, digest.hexdigest()


def _mapping(value: object, keys: frozenset[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        _stop("PREF1_SCHEMA_DRIFT", (label, sorted(keys)))
    return value


def _encode_exact(hex_value: str) -> dict[str, str]:
    rebuilt = float.fromhex(hex_value)
    fraction = Fraction.from_float(rebuilt)
    encoded = {
        "binary64_hex": rebuilt.hex(),
        "denominator": str(fraction.denominator),
        "numerator": str(fraction.numerator),
    }
    if (
        not math.isfinite(rebuilt)
        or rebuilt.hex() != hex_value
        or Fraction.from_float(float.fromhex(encoded["binary64_hex"])) != fraction
    ):
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", hex_value)
    return encoded


def _encode_rational(value: Fraction) -> dict[str, str]:
    if not isinstance(value, Fraction) or isinstance(value, bool):
        _stop("PREF1_EXACT_TYPE_DRIFT", type(value).__name__)
    encoded = {
        "denominator": str(value.denominator),
        "numerator": str(value.numerator),
    }
    rebuilt = Fraction(int(encoded["numerator"]), int(encoded["denominator"]))
    if rebuilt != value or str(rebuilt.numerator) != encoded["numerator"]:
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", encoded)
    return encoded


def _fraction_from_exact(value: object, *, label: str) -> Fraction:
    encoded = _mapping(value, _EXACT_KEYS, label)
    hex_value = encoded["binary64_hex"]
    numerator = encoded["numerator"]
    denominator = encoded["denominator"]
    if (
        not isinstance(hex_value, str)
        or not isinstance(numerator, str)
        or not isinstance(denominator, str)
    ):
        _stop("PREF1_EXACT_TYPE_DRIFT", label)
    try:
        rebuilt = float.fromhex(hex_value)
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise TDG9UR1PREF1Error("PREF1_EXACT_PARSE_DRIFT", label) from exc
    if (
        type(rebuilt) is not float
        or not math.isfinite(rebuilt)
        or rebuilt.hex() != hex_value
        or str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
        or Fraction.from_float(rebuilt) != fraction
    ):
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", label)
    return fraction


def _fraction_from_rational(value: object, *, label: str) -> Fraction:
    encoded = _mapping(value, _RATIONAL_KEYS, label)
    numerator = encoded["numerator"]
    denominator = encoded["denominator"]
    if not isinstance(numerator, str) or not isinstance(denominator, str):
        _stop("PREF1_EXACT_TYPE_DRIFT", label)
    try:
        fraction = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise TDG9UR1PREF1Error("PREF1_EXACT_PARSE_DRIFT", label) from exc
    if (
        str(fraction.numerator) != numerator
        or str(fraction.denominator) != denominator
    ):
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", label)
    return fraction


def downward_binary64_lower(exact: Fraction) -> float:
    """Reproduce stored-lower downward rounding without runtime envelope code."""

    if not isinstance(exact, Fraction) or isinstance(exact, bool) or exact < 0:
        _stop("PREF1_EXACT_LOWER_DRIFT", exact)
    lower = float(exact)
    if not math.isfinite(lower):
        _stop("PREF1_EXACT_LOWER_DRIFT", lower)
    if Fraction.from_float(lower) > exact:
        lower = math.nextafter(lower, 0.0)
    if not math.isfinite(lower) or lower < 0.0:
        _stop("PREF1_EXACT_LOWER_DRIFT", lower)
    return lower


def name_zero_lower_owner(
    *,
    raw: Fraction,
    unclipped: Fraction,
    reproduced: float,
) -> str:
    """Name the exact term that clips a stored lower bound to zero."""

    if reproduced != 0.0:
        _stop("PREF1_ZERO_LOWER_DRIFT", reproduced)
    if raw == 0:
        return "raw_candidate_maximum_is_zero"
    if unclipped <= 0:
        return OWNER_NAME
    if unclipped > 0 and reproduced == 0.0:
        return "downward_binary64_rounding_of_positive_exact_lower"
    _stop("PREF1_ZERO_LOWER_UNCLASSIFIED", (str(raw), str(unclipped)))


def reproduce_envelope_lower(
    *,
    raw: Fraction,
    construction: Fraction,
    arithmetic: Fraction,
) -> dict[str, Any]:
    """Clip by construction plus arithmetic debit; Bernstein is excluded."""

    for name, value in (
        ("raw", raw),
        ("construction", construction),
        ("arithmetic", arithmetic),
    ):
        if not isinstance(value, Fraction) or isinstance(value, bool) or value < 0:
            _stop("PREF1_ENVELOPE_FIELD_DRIFT", name)
    unclipped = raw - construction - arithmetic
    clipped = max(Fraction(0), unclipped)
    stored = downward_binary64_lower(clipped)
    construction_margin = raw - construction
    arithmetic_margin = raw - arithmetic
    if raw == 0:
        _stop("PREF1_RAW_MAXIMUM_ZERO", raw)
    return {
        "owner": name_zero_lower_owner(
            raw=raw, unclipped=unclipped, reproduced=stored
        ),
        "exact_unclipped": unclipped,
        "exact_clipped": clipped,
        "reproduced_stored_lower": stored,
        "construction_over_raw": construction / raw,
        "construction_margin": construction_margin,
        "arithmetic_margin": arithmetic_margin,
        "construction_alone_sufficient": construction_margin <= 0,
        "arithmetic_alone_sufficient": arithmetic_margin <= 0,
    }


def reduce_ur1_terminal(
    *,
    rows_matched: bool,
    intervals_match_ac1: bool,
    owners_equal: bool,
) -> str:
    """Closed four-arm reduction of a completed UR1 diagnostic terminal."""

    if not rows_matched:
        return "row_stream_identity_mismatch"
    if not intervals_match_ac1:
        return "replayed_production_intervals_differ_from_sealed_AC1"
    if owners_equal:
        return RAW_CLASSIFICATION
    return "completed_u_R_rows_match_TI2_but_D01_D12_lower_owners_differ"


def _count(value: object, *, label: str) -> int:
    if type(value) is not int or isinstance(value, bool) or value < 0:
        _stop("PREF1_COUNT_DRIFT", label)
    return value


def _interval_replay_payload(
    hexes: Mapping[str, str],
    counts: Mapping[str, int],
    expected_ratio: Fraction,
) -> dict[str, Any]:
    if counts["endpoint_candidate_count"] != 2 * counts["polynomial_count"]:
        _stop("PREF1_ENDPOINT_COUNT_DRIFT", counts)
    if counts["polynomial_count"] != (
        counts["subinterval_count"] * counts["owned_row_count"]
    ):
        _stop("PREF1_POLYNOMIAL_COUNT_DRIFT", counts)
    raw = Fraction.from_float(float.fromhex(hexes["raw"]))
    construction = Fraction.from_float(float.fromhex(hexes["construction"]))
    arithmetic = Fraction.from_float(float.fromhex(hexes["arithmetic"]))
    certified = Fraction.from_float(float.fromhex(hexes["certified"]))
    upper = Fraction.from_float(float.fromhex(hexes["upper"]))
    reproduced = reproduce_envelope_lower(
        raw=raw, construction=construction, arithmetic=arithmetic
    )
    if (
        reproduced["construction_over_raw"] != expected_ratio
        or reproduced["owner"] != OWNER_NAME
        or reproduced["construction_alone_sufficient"] is not True
        or reproduced["arithmetic_alone_sufficient"] is not False
        or certified < raw
        or upper != certified
        or reproduced["reproduced_stored_lower"] != 0.0
    ):
        _stop("PREF1_INTERVAL_REPLAY_DRIFT", expected_ratio)
    return {
        "ambiguous_discriminant_count": counts["ambiguous_discriminant_count"],
        "arithmetic_alone_sufficient": False,
        "arithmetic_margin": _encode_rational(reproduced["arithmetic_margin"]),
        "bernstein_certification_slack": _encode_exact(hexes["bernstein"]),
        "bernstein_certification_slack_not_part_of_lower_clip": True,
        "certified_at_least_raw": True,
        "certified_continuous_upper_bound": _encode_exact(hexes["certified"]),
        "coefficient_construction_debit": _encode_exact(hexes["construction"]),
        "construction_alone_sufficient": True,
        "construction_margin": _encode_rational(reproduced["construction_margin"]),
        "construction_over_raw": _encode_rational(expected_ratio),
        "endpoint_candidate_count": counts["endpoint_candidate_count"],
        "endpoint_candidate_identity_verified": True,
        "exact_clipped": _encode_rational(reproduced["exact_clipped"]),
        "exact_unclipped": _encode_rational(reproduced["exact_unclipped"]),
        "lower_bound": _encode_exact(hexes["lower"]),
        "outward_arithmetic_debit": _encode_exact(hexes["arithmetic"]),
        "owned_row_count": counts["owned_row_count"],
        "owner": OWNER_NAME,
        "polynomial_count": counts["polynomial_count"],
        "polynomial_count_identity_verified": True,
        "raw_candidate_maximum": _encode_exact(hexes["raw"]),
        "real_interior_root_count": counts["real_interior_root_count"],
        "reproduced_stored_lower": _encode_exact(hexes["lower"]),
        "subinterval_count": counts["subinterval_count"],
        "upper_bound": _encode_exact(hexes["upper"]),
        "upper_equals_certified_upper": True,
    }


def _expected_independent_replay() -> dict[str, Any]:
    d01 = _interval_replay_payload(
        _D01_HEX, _D01_COUNTS, D01_CONSTRUCTION_OVER_RAW
    )
    d12 = _interval_replay_payload(
        _D12_HEX, _D12_COUNTS, D12_CONSTRUCTION_OVER_RAW
    )
    return {
        "D01": d01,
        "D12": d12,
        "D01_D12_lower_owners_equal": True,
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "arithmetic_alone_sufficient": False,
        "bernstein_certification_slack_part_of_lower_clip": False,
        "bernstein_certification_slack_recorded": True,
        "both_owners": OWNER_NAME,
        "classifier": {
            "admission_passed": False,
            "channel": PUBLISHED_CHANNEL,
            "classification": AC1_U_R_CLASS,
            "order_threshold_passed": None,
            "order_threshold_resolved": False,
            "temporal_retry_permitted": True,
        },
        "construction_alone_sufficient": True,
        "production_intervals": {
            "D01": [AC1_U_R_D01_LOWER_HEX, AC1_U_R_D01_UPPER_HEX],
            "D12": [AC1_U_R_D12_LOWER_HEX, AC1_U_R_D12_UPPER_HEX],
            "matches_sealed_AC1": True,
        },
        "raw_terminal_class_matches_independent_reduction": True,
        "replay_receipt": dict(EXPECTED_REPLAY_RECEIPT),
        "row_identity": {
            "algorithm": "sha256",
            "domain": ROW_HASH_DOMAIN,
            "matched": True,
            "related_TI2_retry3_u_R_radius_free_complete_C_class": (
                TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
            ),
            "related_not_replacement_production_admission": True,
            "row_count": OWNED_ROW_COUNT,
            "sha256": TI2_RETRY3_U_R_ROW_SHA256,
        },
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
    }


def _expected_conclusion() -> dict[str, Any]:
    return {
        "PDE_or_candidate_state_opened": False,
        "arithmetic_alone_sufficient_to_clip": False,
        "both_zero_lowers_owned_by_coefficient_plus_arithmetic_debit": True,
        "classification_remains_order_inconclusive": True,
        "construction_alone_sufficient_to_clip": True,
        "independent_method_agreement": False,
        "licenses_only_later_prospectively_frozen_envelope_or_admission_design": True,
        "physics_inference_permitted": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "u_R_rows_match_TI2_and_zero_lowers_owned_by_named_envelope_quantity": True,
    }


def _expected_config() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "raw": {
            "namespace": RAW_NAMESPACE,
            "schema": RAW_SCHEMA,
            "classification": RAW_CLASSIFICATION,
            "runner_id": RUNNER_ID,
            "authority_commit": AUTHORITY_COMMIT,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "terminal_sha256": RAW_TERMINAL_SHA256,
            "leaf_count": 2,
        },
        "authority": {
            "parent_commit": AUTHORITY_PARENT,
            "blob_count": len(_AUTHORITY_BLOBS),
            "blobs": [
                {"path": path, "sha256": digest}
                for path, digest in _AUTHORITY_BLOBS
            ],
        },
        "store": {
            "path": STORE_PATH,
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "snapshot_sha256": SEALED_STORE_SHA256,
            "snapshot_algorithm": "TDG9-LOC1-STORE-SNAPSHOT-v1",
        },
        "selection": {
            "experiment_label": EXPERIMENT_LABEL,
            "tableau_selector": "SSPRK3",
            "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
            "actual_spatial_operator": "inherited_RK4_2049_SBP4",
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "member_key": MEMBER_KEY,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "retry": RETRY,
            "attempted_width_hex": WIDTH_HEX,
            "published_channel": PUBLISHED_CHANNEL,
            "row_hash_domain": ROW_HASH_DOMAIN,
            "TI2_retry3_u_R_row_sha256": TI2_RETRY3_U_R_ROW_SHA256,
            "TI2_retry3_u_R_related_complete_C_class": (
                TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
            ),
            "envelope_owner": "Binary64CubicEnvelope",
            "classifier_owner": "classify_tdg6_channel",
        },
        "work_budget": dict(WORK_BUDGET),
        "decision": {
            "exact_lower_formula": EXACT_LOWER_FORMULA,
            "bernstein_certification_slack_recorded": True,
            "bernstein_certification_slack_part_of_lower_clip": False,
            "stored_lower_reproduced_by_downward_binary64_rounding": True,
            "zero_lower_owner": OWNER_NAME,
            "D01_construction_over_raw": _encode_rational(
                D01_CONSTRUCTION_OVER_RAW
            ),
            "D12_construction_over_raw": _encode_rational(
                D12_CONSTRUCTION_OVER_RAW
            ),
            "construction_alone_sufficient": True,
            "arithmetic_alone_sufficient": False,
            "required_design_classifier": AC1_U_R_CLASS,
            "related_TI2_complete_C_class_is_not_replacement_admission": True,
            "successor_remedy_selected": False,
        },
        "scope": {
            "one_time_live_raw_and_store_authentication": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_store_blind": True,
            "ordinary_verifier_shadow_blind": True,
            "UR1_runner_decision_code_imported": False,
            "AC1_runner_or_authority_imported": False,
            "TI2_runner_or_binder_imported": False,
            "runtime_envelope_implementation_imported": False,
            "LOC1_row_reconstruction_imported": False,
            "seven_proposals_reconstructed": False,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "campaign_store_mutation_authorized": False,
            "raw_namespace_mutation_authorized": False,
            "temporal_retry_admission_authorized": False,
            "PDE_state_commit_authorized": False,
            "successor_remedy_selected": False,
            "common_event_authorized": False,
            "GR0_calibration_authorized": False,
            "candidate_branches_authorized": False,
            "mechanism_result_authorized": False,
            "physical_result_authorized": False,
            "envelope_or_admission_runtime_change_authorized": False,
        },
        "claims": {
            "UR1_raw_result_independently_bound": True,
            "u_R_rows_match_TI2": True,
            "D01_D12_lower_owners_equal": True,
            "zero_lowers_owned_by_coefficient_plus_arithmetic_debit": True,
            "classification_remains_order_inconclusive": True,
            "construction_alone_sufficient": True,
            "arithmetic_alone_sufficient": False,
            "production_method_earned": False,
            "production_SSPRK3_comparator": False,
            "independent_method_agreement": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
            "common_event_completed": False,
            "GR0_calibration_completed": False,
            "candidate_execution_authorized": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
            "retained_EFT_evolution_authorized": False,
            "physical_transition_claim_authorized": False,
        },
    }


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG9UR1PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
    if value != _expected_config():
        _stop("PREF1_CONFIG_DRIFT", "typed PREF1 config differs")
    return value


def expected_evidence(config_raw: bytes) -> dict[str, Any]:
    config = _config(config_raw)
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "raw": dict(config["raw"]),
        "authority": dict(config["authority"]),
        "store": dict(config["store"]),
        "selection": dict(config["selection"]),
        "work_budget": dict(config["work_budget"]),
        "decision": dict(config["decision"]),
        "scope": dict(config["scope"]),
        "claims": dict(config["claims"]),
    }


def expected_bound_payload(config_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    sealed = {
        "leaf_count": SEALED_STORE_LEAF_COUNT,
        "sha256": SEALED_STORE_SHA256,
    }
    return {
        **expected,
        "raw_binding": {
            "canonical_duplicate_free_schema_verified": True,
            "leaf_count": 2,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "raw_namespace_unchanged": True,
            "terminal_sha256": RAW_TERMINAL_SHA256,
        },
        "authority_binding": {
            "all_committed_blob_hashes_verified": True,
            "blob_count": len(_AUTHORITY_BLOBS),
            "commit": AUTHORITY_COMMIT,
            "parent_commit": AUTHORITY_PARENT,
        },
        "store_binding": {
            "binder_snapshot_after": dict(sealed),
            "binder_snapshot_before": dict(sealed),
            "raw_snapshot_after": dict(sealed),
            "raw_snapshot_before": dict(sealed),
            "store_unchanged": True,
        },
        "independent_replay": _expected_independent_replay(),
        "conclusion": _expected_conclusion(),
    }


def expected_compact_result(config_raw: bytes) -> dict[str, Any]:
    return {
        "artifact_id": ARTIFACT_ID,
        "artifact_payload": expected_bound_payload(config_raw),
    }


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT:
        _stop("PREF1_AUTHORITY_DRIFT", manifest.get("authority_commit"))
    if _git(root, "cat-file", "-t", AUTHORITY_COMMIT).strip() != b"commit":
        _stop("PREF1_AUTHORITY_DRIFT", "authority object is not a commit")
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [AUTHORITY_PARENT.encode("ascii")]:
        _stop("PREF1_AUTHORITY_LINEAGE_DRIFT", parents)
    if len(_AUTHORITY_BLOBS) != 16:
        _stop("PREF1_AUTHORITY_BLOB_COUNT_DRIFT", len(_AUTHORITY_BLOBS))
    for relative, expected in _AUTHORITY_BLOBS:
        observed = sha256(
            _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        ).hexdigest()
        if observed != expected:
            _stop("PREF1_AUTHORITY_BLOB_DRIFT", relative)


def _expected_manifest() -> dict[str, Any]:
    return {
        "AC1_PREF1_result_sha256": AC1_PREF1_RESULT_SHA256,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": UR1_ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "experiment_label": EXPERIMENT_LABEL,
        "output_leaves": ["manifest.json", "terminal.json"],
        "production_SSPRK3_comparator": False,
        "published_channel": PUBLISHED_CHANNEL,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "tableau_selector": "SSPRK3",
    }


def _expected_terminal_base() -> dict[str, Any]:
    sealed = {
        "leaf_count": SEALED_STORE_LEAF_COUNT,
        "sha256": SEALED_STORE_SHA256,
    }
    return {
        "GR0_calibration_completed": False,
        "PDE_state_committed": False,
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": UR1_ARTIFACT_ID,
        "attempted_width_hex": WIDTH_HEX,
        "authority_commit": AUTHORITY_COMMIT,
        "candidate_branch_opened": False,
        "classification": RAW_CLASSIFICATION,
        "common_event_completed": False,
        "completed_result_licenses_only_later_prospectively_frozen_envelope_or_admission_design": True,
        "experiment_label": EXPERIMENT_LABEL,
        "fine_path_committed": False,
        "independent_method_agreement": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "published_channel": PUBLISHED_CHANNEL,
        "related_TI2_complete_C_class_is_not_replacement_production_admission": True,
        "retry": RETRY,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "store_snapshot_after": sealed,
        "store_snapshot_before": sealed,
        "store_unchanged": True,
        "successor_remedy_selected": False,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "tableau_selector": "SSPRK3",
        "temporal_retry_admission_called": False,
    }


def _validate_replay_receipt(value: object) -> dict[str, Any]:
    receipt = _mapping(value, frozenset(EXPECTED_REPLAY_RECEIPT), "replay_receipt")
    observed = {key: receipt[key] for key in EXPECTED_REPLAY_RECEIPT}
    if observed != EXPECTED_REPLAY_RECEIPT:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", observed)
    return dict(EXPECTED_REPLAY_RECEIPT)


def _bind_interval(
    value: object,
    expected: Mapping[str, Any],
    *,
    label: str,
) -> dict[str, Any]:
    payload = _mapping(value, _INTERVAL_KEYS, label)
    envelope = _mapping(payload.get("envelope"), _ENVELOPE_KEYS, f"{label}.envelope")
    owner = _mapping(
        payload.get("zero_lower_owner"), _OWNER_KEYS, f"{label}.owner"
    )
    polynomial_count = _count(
        envelope.get("polynomial_count"), label=f"{label}.polynomial_count"
    )
    endpoint_count = _count(
        envelope.get("endpoint_candidate_count"),
        label=f"{label}.endpoint_candidate_count",
    )
    interior = _count(
        envelope.get("real_interior_root_count"),
        label=f"{label}.real_interior_root_count",
    )
    ambiguous = _count(
        envelope.get("ambiguous_discriminant_count"),
        label=f"{label}.ambiguous_discriminant_count",
    )
    subinterval_count = _count(
        payload.get("subinterval_count"), label=f"{label}.subinterval_count"
    )
    owned_row_count = _count(
        payload.get("owned_row_count"), label=f"{label}.owned_row_count"
    )
    if endpoint_count != 2 * polynomial_count:
        _stop("PREF1_ENDPOINT_COUNT_DRIFT", label)
    if polynomial_count != subinterval_count * owned_row_count:
        _stop("PREF1_POLYNOMIAL_COUNT_DRIFT", label)
    raw = _fraction_from_exact(
        envelope.get("raw_candidate_maximum"), label=f"{label}.raw"
    )
    construction = _fraction_from_exact(
        envelope.get("coefficient_construction_debit"),
        label=f"{label}.construction",
    )
    arithmetic = _fraction_from_exact(
        envelope.get("outward_arithmetic_debit"), label=f"{label}.arithmetic"
    )
    bernstein = _fraction_from_exact(
        envelope.get("bernstein_certification_slack"),
        label=f"{label}.bernstein",
    )
    certified = _fraction_from_exact(
        envelope.get("certified_continuous_upper_bound"),
        label=f"{label}.certified",
    )
    lower = _fraction_from_exact(
        payload.get("lower_bound"), label=f"{label}.lower"
    )
    upper = _fraction_from_exact(
        payload.get("upper_bound"), label=f"{label}.upper"
    )
    if any(
        item < 0
        for item in (raw, construction, arithmetic, bernstein, certified, lower, upper)
    ):
        _stop("PREF1_NEGATIVE_ENVELOPE_DRIFT", label)
    if certified < raw or upper != certified or lower != 0:
        _stop("PREF1_CERTIFIED_UPPER_DRIFT", label)
    owner_raw = _fraction_from_exact(
        owner.get("raw_candidate_maximum"), label=f"{label}.owner.raw"
    )
    owner_construction = _fraction_from_exact(
        owner.get("coefficient_construction_debit"),
        label=f"{label}.owner.construction",
    )
    owner_arithmetic = _fraction_from_exact(
        owner.get("outward_arithmetic_debit"),
        label=f"{label}.owner.arithmetic",
    )
    owner_bernstein = _fraction_from_exact(
        owner.get("bernstein_certification_slack"),
        label=f"{label}.owner.bernstein",
    )
    if (
        owner_raw != raw
        or owner_construction != construction
        or owner_arithmetic != arithmetic
        or owner_bernstein != bernstein
    ):
        _stop("PREF1_OWNER_ENVELOPE_DRIFT", label)
    if owner.get("owner") not in ZERO_LOWER_OWNERS:
        _stop("PREF1_OWNER_ENUM_DRIFT", owner.get("owner"))
    if owner.get("bernstein_certification_slack_not_part_of_lower_clip") is not True:
        _stop("PREF1_BERNSTEIN_CLIP_FLAG_DRIFT", label)
    reproduced = reproduce_envelope_lower(
        raw=raw, construction=construction, arithmetic=arithmetic
    )
    unclipped = _fraction_from_rational(
        owner.get("exact_unclipped"), label=f"{label}.unclipped"
    )
    clipped = _fraction_from_rational(
        owner.get("exact_clipped"), label=f"{label}.clipped"
    )
    stored = _fraction_from_exact(
        owner.get("reproduced_stored_lower"), label=f"{label}.reproduced"
    )
    if (
        unclipped != reproduced["exact_unclipped"]
        or clipped != reproduced["exact_clipped"]
        or stored != Fraction.from_float(reproduced["reproduced_stored_lower"])
        or stored != lower
        or owner.get("owner") != reproduced["owner"]
        or reproduced["owner"] != OWNER_NAME
        or reproduced["construction_alone_sufficient"] is not True
        or reproduced["arithmetic_alone_sufficient"] is not False
        or reproduced["construction_over_raw"]
        != _fraction_from_rational(
            expected["construction_over_raw"], label=f"{label}.ratio"
        )
    ):
        _stop("PREF1_OWNER_RECOMPUTE_DRIFT", label)
    source = owner.get("source")
    if not isinstance(source, Mapping) or set(source) != _DEBIT_SOURCE_KEYS:
        _stop("PREF1_OWNER_SOURCE_DRIFT", label)
    if source.get("term") != (
        "coefficient_construction_debit_plus_outward_arithmetic_debit"
    ):
        _stop("PREF1_OWNER_SOURCE_DRIFT", source.get("term"))
    if (
        _fraction_from_exact(
            source.get("raw_candidate_maximum"), label=f"{label}.source.raw"
        )
        != raw
        or _fraction_from_exact(
            source.get("coefficient_construction_debit"),
            label=f"{label}.source.construction",
        )
        != construction
        or _fraction_from_exact(
            source.get("outward_arithmetic_debit"),
            label=f"{label}.source.arithmetic",
        )
        != arithmetic
    ):
        _stop("PREF1_OWNER_SOURCE_IDENTITY_DRIFT", label)
    summary = {
        "ambiguous_discriminant_count": ambiguous,
        "arithmetic_alone_sufficient": False,
        "arithmetic_margin": _encode_rational(reproduced["arithmetic_margin"]),
        "bernstein_certification_slack": _encode_exact(
            envelope["bernstein_certification_slack"]["binary64_hex"]
        ),
        "bernstein_certification_slack_not_part_of_lower_clip": True,
        "certified_at_least_raw": True,
        "certified_continuous_upper_bound": _encode_exact(
            envelope["certified_continuous_upper_bound"]["binary64_hex"]
        ),
        "coefficient_construction_debit": _encode_exact(
            envelope["coefficient_construction_debit"]["binary64_hex"]
        ),
        "construction_alone_sufficient": True,
        "construction_margin": _encode_rational(reproduced["construction_margin"]),
        "construction_over_raw": dict(expected["construction_over_raw"]),
        "endpoint_candidate_count": endpoint_count,
        "endpoint_candidate_identity_verified": True,
        "exact_clipped": _encode_rational(reproduced["exact_clipped"]),
        "exact_unclipped": _encode_rational(reproduced["exact_unclipped"]),
        "lower_bound": _encode_exact(payload["lower_bound"]["binary64_hex"]),
        "outward_arithmetic_debit": _encode_exact(
            envelope["outward_arithmetic_debit"]["binary64_hex"]
        ),
        "owned_row_count": owned_row_count,
        "owner": OWNER_NAME,
        "polynomial_count": polynomial_count,
        "polynomial_count_identity_verified": True,
        "raw_candidate_maximum": _encode_exact(
            envelope["raw_candidate_maximum"]["binary64_hex"]
        ),
        "real_interior_root_count": interior,
        "reproduced_stored_lower": _encode_exact(
            owner["reproduced_stored_lower"]["binary64_hex"]
        ),
        "subinterval_count": subinterval_count,
        "upper_bound": _encode_exact(payload["upper_bound"]["binary64_hex"]),
        "upper_equals_certified_upper": True,
    }
    if summary != expected:
        _stop("PREF1_INTERVAL_SUMMARY_DRIFT", label)
    return dict(expected)


def _bind_row_stream(value: object) -> dict[str, Any]:
    payload = _mapping(value, _ROW_MATCH_KEYS, "row_stream")
    observed = payload.get("observed_sha256")
    if (
        payload.get("algorithm") != "sha256"
        or payload.get("domain") != ROW_HASH_DOMAIN
        or payload.get("row_count") != OWNED_ROW_COUNT
        or payload.get("expected_sha256") != TI2_RETRY3_U_R_ROW_SHA256
        or not isinstance(observed, str)
        or len(observed) != 64
        or any(character not in _HEX for character in observed)
        or payload.get("matched") is not True
        or observed != TI2_RETRY3_U_R_ROW_SHA256
        or payload.get("related_TI2_retry3_u_R_radius_free_complete_C_class")
        != TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
        or payload.get("related_not_replacement_production_admission") is not True
    ):
        _stop("PREF1_ROW_STREAM_DRIFT", observed)
    return {
        "algorithm": "sha256",
        "domain": ROW_HASH_DOMAIN,
        "matched": True,
        "related_TI2_retry3_u_R_radius_free_complete_C_class": (
            TI2_RETRY3_U_R_RELATED_COMPLETE_C_CLASS
        ),
        "related_not_replacement_production_admission": True,
        "row_count": OWNED_ROW_COUNT,
        "sha256": TI2_RETRY3_U_R_ROW_SHA256,
    }


def _bind_u_R(
    value: object, expected_replay: Mapping[str, Any]
) -> tuple[dict[str, Any], bool, bool]:
    payload = _mapping(value, _U_R_KEYS, "u_R")
    if payload.get("channel") != PUBLISHED_CHANNEL:
        _stop("PREF1_CHANNEL_DRIFT", payload.get("channel"))
    d01 = _bind_interval(payload.get("D01"), expected_replay["D01"], label="D01")
    d12 = _bind_interval(payload.get("D12"), expected_replay["D12"], label="D12")
    outer = CertifiedMagnitudeInterval(
        _fraction_from_exact(payload["D01"]["lower_bound"], label="D01.lower"),
        _fraction_from_exact(payload["D01"]["upper_bound"], label="D01.upper"),
    )
    finest = CertifiedMagnitudeInterval(
        _fraction_from_exact(payload["D12"]["lower_bound"], label="D12.lower"),
        _fraction_from_exact(payload["D12"]["upper_bound"], label="D12.upper"),
    )
    try:
        decision = classify_tdg6_channel(outer, finest)
    except (TypeError, ValueError) as exc:
        raise TDG9UR1PREF1Error("PREF1_RECLASSIFY_DRIFT", PUBLISHED_CHANNEL) from exc
    if (
        decision.classification != payload.get("classification")
        or decision.admission_passed is not payload.get("admission_passed")
        or decision.temporal_retry_permitted
        is not payload.get("temporal_retry_permitted")
        or decision.order_threshold_resolved
        is not payload.get("order_threshold_resolved")
        or decision.order_threshold_passed is not payload.get("order_threshold_passed")
        or decision.classification != AC1_U_R_CLASS
        or payload.get("admission_passed") is not False
        or decision.finest_pair_debit != finest.upper
    ):
        _stop("PREF1_RECLASSIFY_DRIFT", decision.classification)
    matches_ac1 = (
        payload["D01"]["lower_bound"]["binary64_hex"] == AC1_U_R_D01_LOWER_HEX
        and payload["D01"]["upper_bound"]["binary64_hex"] == AC1_U_R_D01_UPPER_HEX
        and payload["D12"]["lower_bound"]["binary64_hex"] == AC1_U_R_D12_LOWER_HEX
        and payload["D12"]["upper_bound"]["binary64_hex"] == AC1_U_R_D12_UPPER_HEX
        and payload.get("classification") == AC1_U_R_CLASS
    )
    if payload.get("matches_sealed_AC1_production_intervals") is not matches_ac1:
        _stop("PREF1_AC1_INTERVAL_DRIFT", matches_ac1)
    owners_equal = d01["owner"] == d12["owner"]
    if payload.get("D01_D12_lower_owners_equal") is not owners_equal:
        _stop("PREF1_OWNER_EQUALITY_DRIFT", owners_equal)
    if (
        d01["owner"] != OWNER_NAME
        or d12["owner"] != OWNER_NAME
        or owners_equal is not True
        or matches_ac1 is not True
    ):
        _stop("PREF1_OWNER_SET_DRIFT", (d01["owner"], d12["owner"]))
    return dict(expected_replay["classifier"]), matches_ac1, owners_equal


def bind_raw_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    expected = expected_bound_payload(config_raw)
    root = Path(os.path.abspath(os.fspath(repository)))
    metadata = root.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _stop("PREF1_REPOSITORY_UNSAFE", root)
    manifest_raw, terminal_raw, manifest, terminal = _raw_snapshot(root)
    if manifest != _expected_manifest():
        _stop("PREF1_RAW_MANIFEST_DRIFT", "manifest schema differs")
    _authenticate_authority(root, manifest)
    before_store = _snapshot_store(root)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SHA256)
    if before_store != sealed_store:
        _stop("PREF1_STORE_IDENTITY_DRIFT", before_store)
    base = _expected_terminal_base()
    if set(terminal) != set(base) | _COMPLETED_EXTRA_KEYS:
        _stop("PREF1_TERMINAL_SCHEMA_DRIFT", sorted(terminal))
    if any(terminal.get(name) != value for name, value in base.items()):
        _stop("PREF1_TERMINAL_NONCLAIM_DRIFT", terminal.get("classification"))
    receipt = _validate_replay_receipt(terminal.get("replay_receipt"))
    if (
        terminal.get("shadow_path_count") != 7
        or terminal.get("shadow_proposal_count") != 7
        or terminal.get("SSPRK3_stage_and_endpoint_record_count") != 28
    ):
        _stop("PREF1_SHADOW_COUNT_DRIFT", "completed terminal counts")
    row_identity = _bind_row_stream(terminal.get("row_stream"))
    replay = expected["independent_replay"]
    _classifier, matches_ac1, owners_equal = _bind_u_R(
        terminal.get("u_R"), replay
    )
    independent_class = reduce_ur1_terminal(
        rows_matched=True,
        intervals_match_ac1=matches_ac1,
        owners_equal=owners_equal,
    )
    if (
        independent_class != RAW_CLASSIFICATION
        or terminal.get("classification") != independent_class
        or receipt != EXPECTED_REPLAY_RECEIPT
        or row_identity != replay["row_identity"]
    ):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", independent_class)
    after_store = _snapshot_store(root)
    manifest_after, terminal_after = _raw_pair(root)
    if (
        after_store != before_store
        or manifest_after != manifest_raw
        or terminal_after != terminal_raw
    ):
        _stop("PREF1_INPUT_MUTATED_DURING_BIND", "raw/store input changed")
    return expected


def build_pref1_result(
    config_raw: bytes, repository: Path, *, live: bool = True
) -> dict[str, Any]:
    payload = (
        bind_raw_result(config_raw, repository)
        if live
        else expected_evidence(config_raw)
    )
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = expected_compact_result(config_raw)
    if sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256:
        _stop("PREF1_COMPACT_DRIFT", "compact result hash differs")
    result = _json(result_raw, "compact result")
    if result != expected:
        _stop("PREF1_COMPACT_DRIFT", "compact payload differs")
    return result


__all__ = (
    "ARTIFACT_ID",
    "COMPACT_RESULT_SHA256",
    "CONFIG_PATH",
    "D01_CONSTRUCTION_OVER_RAW",
    "D12_CONSTRUCTION_OVER_RAW",
    "OWNER_DOCUMENT",
    "OWNER_NAME",
    "RESULT_PATH",
    "TDG9UR1PREF1Error",
    "bind_raw_result",
    "build_pref1_result",
    "canonical_result",
    "downward_binary64_lower",
    "expected_bound_payload",
    "expected_compact_result",
    "expected_evidence",
    "name_zero_lower_owner",
    "reduce_ur1_terminal",
    "reproduce_envelope_lower",
    "validate_compact_result",
)
