"""Independent post-result binder for the completed TDG10 QA2 terminal.

The explicit live path authenticates the frozen authority commit, reads the
two-leaf raw namespace through no-follow/race-checked descriptors, restores
the retry-4 and retry-5 predecessors with lower-level campaign mechanics, and
independently rebuilds fourteen SSPRK3-on-inherited-SBP4 shadow proposals.  It
then reconstructs all thirty-six complete-state channel streams with exact
Fraction arithmetic and two independent local-extrema routes.  It imports no
QA2 decision code and never serializes or admits an evolved state.

Ordinary compact verification is intentionally raw-, store-, shadow-, and
Git-blind.  It checks the tracked certificate's canonical schema, exact
interval evidence, all-of reductions, and frozen nonclaims.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict, dataclass, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import platform
import stat
import subprocess
import sys
import tomllib
from typing import Any, Mapping, NoReturn, Sequence

import numpy as np

from . import hlt16_campaign_runtime as campaign_runtime
from . import proto15_runtime as p15
from . import tdg6_temporal_admission_runtime as tdg6
from .hlt16_campaign_store import HLT16CampaignStore
from .numerical_engine import COMPARATOR_METHOD, PRIMARY_METHOD, array_content_sha256
from .proto19_gr0_static_factory import build_static_gr0_shells
from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_ORDER_SQUARED_MULTIPLIER,
    classify_tdg6_channel,
    three_halves_order_passes_squared,
)
from .tdg9_exact_temporal_arithmetic import (
    exact_hermite_coefficients,
    restrict_exact_cubic_to_half,
    subtract_exact_cubics,
)
from .tdg9_local_extrema import (
    EVALUATOR_ID as PRIMARY_LOCALIZER_ID,
    LocalCubic,
    localize_absolute_maximum,
)
from .tdg9_local_extrema import _stationary_intervals as _primary_stationary_intervals
from .tdg9_local_extrema_independent_v2 import (
    EVALUATOR_ID as INDEPENDENT_LOCALIZER_ID,
    IndependentLocalCubicV2,
    RootIsolationInconclusive,
    localize_absolute_maximum_independently_v2,
)


SCHEMA_VERSION = 1
ARTIFACT_ID = "FGC-1-TDG10-QA2-PREF1"
CLASSIFICATION = (
    "independently_bound_retry4_5_ssprk3_sbp4_exact_complete_C_"
    "two_width_nonpass_terminal"
)
CONFIG_PATH = "configs/fgc/fgc-1-tdg10-qa2-pref1.toml"
RESULT_PATH = "results/fgc-1-tdg10-qa2-pref1.json"
OWNER_DOCUMENT = "docs/fgc-tdg10-qa2-pref1.md"
TARGET_PROTOCOL = "FGC-2-SF1-PROTO18"
# Filled by the integration owner after the one authorized live bind creates
# the no-clobber compact result.  Until then semantic/canonical validation is
# fully active, while no absent artifact hash is invented.
COMPACT_RESULT_SHA256: str | None = (
    "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6"
)
ENVIRONMENT = {
    "python_implementation": "CPython",
    "python_version": "3.14.3",
    "numpy_version": "2.5.1",
    "system": "Darwin",
    "machine": "arm64",
    "byteorder": "little",
}

AUTHORITY_COMMIT = "615144fe72b402835acd6dd89d28f89cbe690cce"
AUTHORITY_PARENT = "a7915925a23ceaff6a491b8d08c572a70d242985"
RAW_NAMESPACE = "runs/fgc-2-sf1/tdg10-qa2/retries4-5-ssprk3-sbp4-exact-complete-c"
RAW_SCHEMA = "FGC-1-TDG10-QA2-raw-v1"
RAW_CLASSIFICATION = "completed_one_or_more_width_all_channel_nonpass_no_state_advance"
BOTH_PASS_CLASS = "completed_both_width_all_channel_pass_no_state_advance"
INCONCLUSIVE_CLASS = "exact_route_or_resource_inconclusive_no_state_advance"
PREMISE_STOP_CLASS = "shadow_or_scientific_premise_stop_no_state_advance"
RUNNER_ID = "FGC-1-TDG10-QA2-RUN1"
QA2_ARTIFACT_ID = "FGC-1-TDG10-QA2-FRZ1"
RAW_MANIFEST_SHA256 = "8638ac8caf2fb4adc0ff419e2147080e78e13e4cb3d5640f38fd841a261fd1eb"
RAW_TERMINAL_SHA256 = "252810f4a388ab87e382c33c7cdd7fa576be67a317752828a89c247def49867c"

STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
SEALED_STORE_LEAF_COUNT = 115
SEALED_STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
TABLEAU_RUNTIME_SELECTOR = COMPARATOR_METHOD
INTERVAL_OWNER = "exact_radius_free_complete_C_dual_rational_localizer"
MEMBER_KEY = "RK4-2049"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
MEMBER_DESCRIPTOR_SHA256 = (
    "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
)
COORDINATES_SHA256 = "2a347b314de5e1e94ac0c999db283e03498070ece6762c53252fa41c1e04a646"
GRID_SPACING_HEX = "0x1.0000000000000p-4"
OUTER_RADIUS_HEX = "0x1.0000000000000p+7"
PRIMARY_REFINEMENT_DEPTH = 160
PER_CHANNEL_MAXIMUM_CANDIDATES_D01 = 16352
PER_CHANNEL_MAXIMUM_CANDIDATES_D12 = 32704
D01_POLYNOMIAL_COUNT = 2 * OWNED_ROW_COUNT
D12_POLYNOMIAL_COUNT = 4 * OWNED_ROW_COUNT
CHANNEL_ORDER = TDG6_COMPLETE_STATE_CHANNELS
SUFFICIENT_CONDITION = "8*U12^2 <= L01^2"
PRIMARY_EVALUATOR_ID = "tdg9_loc1_derivative_monotone_bisection_v1"
INDEPENDENT_EVALUATOR_ID = "tdg9_loc2_endpoint_deflated_adaptive_discriminant_v2"
TDG10_EVALUATOR_ID = "tdg10_exact_complete_c_dual_route_v1"
_ROW_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-ROW-STREAM-v1\n"
_COEFFICIENT_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-COMBINED-v1\n"
_SURVIVOR_KEY_HASH_DOMAIN = b"TDG10-EXACT-COMPLETE-C-SURVIVOR-KEY-STREAM-v1\n"
_STATIONARY_COUNT_DOMAIN = b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n"
_D12_STREAM_DOMAIN = b"TDG10-QA2-PREF1-RAW-D12-UPPER-STREAM-v1\n"

WIDTHS: tuple[dict[str, Any], ...] = (
    {
        "retry": 4,
        "prior_retry_count": 3,
        "generation": 10,
        "forbidden_predecessor_generation": 11,
        "checkpoint_sha256": (
            "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187"
        ),
        "checkpoint_raw_sha256": (
            "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b"
        ),
        "checkpoint_journal_sequence": 12,
        "checkpoint_journal_sha256": (
            "d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4"
        ),
        "checkpoint_journal_raw_sha256": (
            "dd02233dd5e4cf28becbc4fe29259d0879005d806b70493aa2281255081eb219"
        ),
        "historical_rejection_sequence": 13,
        "historical_rejection_sha256": (
            "1bd21553f3b6e0af613997da094a46d4dcd966a0e4346692bae62b05046862ca"
        ),
        "historical_rejection_raw_sha256": (
            "26542403a17b05664d5a447b88860186d5f554b1f4d555624493ed40b941694c"
        ),
        "attempted_width_hex": "0x1.aaa9612df8000p-12",
        "transaction_sha256": (
            "abab435c8222de0a86368e9562a2949578de38fe822c6032319410b49e4c6f6c"
        ),
        "historical_journal_sha256": (
            "fc68ce01980cb7a66f82de59bbd59d5e53ace08537b6b0dbfe0be1e070e6912f"
        ),
        "failed_channels": ("u:alpha", "u:lambda", "u:R"),
        "width_class": "one_or_more_channel_nonpass",
        "d12_upper_stream_sha256": (
            "e59c564eafc023495932a5f2578255c90332c01b323786cd52499546d80c54f4"
        ),
    },
    {
        "retry": 5,
        "prior_retry_count": 4,
        "generation": 11,
        "forbidden_predecessor_generation": 12,
        "checkpoint_sha256": (
            "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8"
        ),
        "checkpoint_raw_sha256": (
            "786bdcb835e043c65ddd33fad3a2af0e3d33e768dc26a9f9f00d2757b0da3ea5"
        ),
        "checkpoint_journal_sequence": 14,
        "checkpoint_journal_sha256": (
            "0810018134df958b8dab32bd9c0ecd7ebfc4376fbc87869a8bee7dd02eb1e199"
        ),
        "checkpoint_journal_raw_sha256": (
            "a1882d904091591d7fa3ccf3eab04aa748161115f97c15abd57e535e6fa01f3d"
        ),
        "historical_rejection_sequence": 15,
        "historical_rejection_sha256": (
            "475014d918952adef032d71b81073773364a37fa795491e739aae609f77e8537"
        ),
        "historical_rejection_raw_sha256": (
            "1997342415f74f330ce6b45ca6fd550ad0e5e674b85148fea5047ff82ac3f152"
        ),
        "attempted_width_hex": "0x1.aaa9612df0000p-13",
        "transaction_sha256": (
            "69fd069034fdd3e2a80a5f3fa989b639b6659105a8ff3e8742e32e78739fd461"
        ),
        "historical_journal_sha256": (
            "b602d12ab00757090734ac63042f172a372154d51b41c7d286c7c0a758b62221"
        ),
        "failed_channels": ("u:alpha", "u:lambda", "u:R"),
        "width_class": "one_or_more_channel_nonpass",
        "d12_upper_stream_sha256": (
            "94add3b1a6dc90041827d915a881d124b035925698a7c7ad49eb57ed9553ed99"
        ),
    },
)

_AUTHORITY_BLOBS: tuple[tuple[str, str], ...] = (
    ("Makefile", "5e84d1ef3a6cf951ffad13f0b6d39b974accaaab54124e959b10d3231389ea64"),
    ("README.md", "e5e6304f90730735ef162988721d3f26b5b74d7c316ca3caae577d65fecda80d"),
    (
        "configs/fgc/fgc-1-tdg10-qa2-frz1.toml",
        "5d32295113d8c675774319b3313426140b0ca4902f956541182bbc40c2949323",
    ),
    (
        "docs/claim-ledger.md",
        "dcae54e4ef16bb9c808dfe659f28099daaeced28090624e8db2faee4cc631ae4",
    ),
    (
        "docs/fgc-runtime-matrix.md",
        "c890f748e0ba82033f12956e8ee1971f467cbfc74e94a70634a4cdb917600cb2",
    ),
    (
        "docs/fgc-tdg10-qa2-frz1.md",
        "2f703fc2e5e1fb64f775c7eb02c136663d5d25af6d2a54e9df8a015884b92daa",
    ),
    (
        "docs/research-roadmap.md",
        "4d76b03b2a4ff901283494d73cc4394ddc5003bc00328e085d893ac08ee0152a",
    ),
    (
        "results/README.md",
        "fbd974b2bb87eac6ac0330661b713898f445cde5fdc5ffd9c63fbd93494fdf36",
    ),
    (
        "results/fgc-1-tdg10-qa2-frz1.json",
        "a2e1e7345b4fa9b2a4810dd0385ee592e0c1dccdb11e2a626ecbb9de9215868e",
    ),
    (
        "scripts/check_repo.py",
        "41d461e9191e8bacb2fba66c8e39b1cdab76a7ff0cda01b206dc86df4e016136",
    ),
    (
        "scripts/reproduce_fgc_tdg10_qa2_frz1.py",
        "c3f768c354206f2e54a76b2fcbf5be59b2bb2780816a7b886750ae1d61af29fa",
    ),
    (
        "scripts/run_fgc_tdg10_qa2.py",
        "9893dbfeddeadd93f4b72a93fc21563649edda6fcf46238d6d616622168a024e",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg10_qa2_authority.py",
        "5c1344691976315f06f306f412fbf018124888745ebb1aee1a67d7b4c1f0f690",
    ),
    (
        "tests/test_check_repo_tdg10_qa2_frz1.py",
        "68ef3d991002028eba4cdabad94d7eb46f84eb79dcccdcdca9451ceaf426b42a",
    ),
    (
        "tests/test_fgc_tdg10_qa2_authority.py",
        "404ed37e0dafb8db18b944f07b5e8f24ac5a8b01f04975c64b63913f49f15d8c",
    ),
    (
        "tests/test_fgc_tdg10_qa2_runner.py",
        "e747434429b2d306783a0f1e2c9b0429bdd5633b158d7dfd47c4d643b1a807ab",
    ),
)

_LIVE_DEPENDENCY_BLOBS: tuple[tuple[str, str], ...] = (
    (
        "src/recursive_horizons/fgc/evolution/hlt16_campaign_runtime.py",
        "ff60118a431aa99ea3b69bc8dfc3a38f1cde2d7bd2b7a9e9f3e7269a4d498d72",
    ),
    (
        "src/recursive_horizons/fgc/evolution/proto15_runtime.py",
        "73f07ae26e90642bc6e5158d9a85bfa499d8082c13e082fc8f19e968d8e624bc",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
    ),
    (
        "src/recursive_horizons/fgc/evolution/hlt16_campaign_store.py",
        "4625441983d0e014bf7d8b1bfdee512e2781ccf88213e2ccae3474317b22ee7e",
    ),
    (
        "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
    ),
    (
        "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py",
        "735545419a39e3fa5417c020d688cfa4b0ce4770446ae1083d0f4b41e039a5ad",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py",
        "87d9a90e272a89479690c0d432f2d5cd979d88841615aaaeb9b49a2c0cca651e",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_exact_temporal_arithmetic.py",
        "2c9a0c1a9c6560ed787a605ef508ede7471953ceb222527889c5622ea6d8e82d",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_local_extrema.py",
        "b97bca8c49313ad28322a1236143e03c5a7f0d3b6bcdb2940965e531a122b28a",
    ),
    (
        "src/recursive_horizons/fgc/evolution/tdg9_local_extrema_independent_v2.py",
        "750c38bb6bad4b37d946f66708055f3dc0835c807244e9c25c6e87e4e87e0511",
    ),
)

_PREF1_PRE_RESULT_PATHS = tuple(
    sorted(
        (
            "Makefile",
            "README.md",
            CONFIG_PATH,
            "docs/claim-ledger.md",
            "docs/fgc-runtime-matrix.md",
            OWNER_DOCUMENT,
            "docs/research-roadmap.md",
            "results/README.md",
            "scripts/check_repo.py",
            "scripts/reproduce_fgc_tdg10_qa2_pref1.py",
            "src/recursive_horizons/fgc/evolution/tdg10_qa2_pref1_binder.py",
            "tests/test_check_repo_tdg10_qa2_pref1.py",
            "tests/test_fgc_tdg10_qa2_pref1_binder.py",
        )
    )
)
_ALLOWED_ADDITIONAL_UNTRACKED = (".qdrant-initialized",)

_MAX_RAW_LEAF_BYTES = 512 * 1024 * 1024
_MAX_STORE_LEAF_BYTES = 512 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")
_RATIONAL_KEYS = frozenset({"numerator", "denominator"})
_LOCALIZATION_CLASSES = frozenset(
    {"unique_maximum", "nonunique_or_interval_inconclusive"}
)
_RAW_TERMINAL_CLASSES = frozenset(
    {RAW_CLASSIFICATION, BOTH_PASS_CLASS, INCONCLUSIVE_CLASS, PREMISE_STOP_CLASS}
)


class TDG10QA2PREF1Error(ValueError):
    """The compact contract, raw terminal, store, or independent replay differs."""

    def __init__(self, stop_id: str, detail: object) -> None:
        self.stop_id = str(stop_id)
        self.detail = " ".join(str(detail).split())[:640]
        super().__init__(f"{self.stop_id}: {self.detail}")


def _stop(stop_id: str, detail: object) -> NoReturn:
    raise TDG10QA2PREF1Error(stop_id, detail)


def _canonical(value: object) -> bytes:
    try:
        return json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("ascii")
    except (TypeError, ValueError) as exc:
        raise TDG10QA2PREF1Error("PREF1_CANONICAL_DRIFT", exc) from exc


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
        raise TDG10QA2PREF1Error("PREF1_COMPACT_DRIFT", exc) from exc


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
        raise TDG10QA2PREF1Error("PREF1_JSON_DRIFT", label) from exc
    if not isinstance(value, dict) or canonical_result(value) != raw:
        _stop("PREF1_JSON_DRIFT", f"{label} is not canonical pretty JSON")
    return value


def _git(root: Path, *arguments: str) -> bytes:
    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update({"GIT_OPTIONAL_LOCKS": "0", "LC_ALL": "C", "LANG": "C"})
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
    except (OSError, subprocess.CalledProcessError) as exc:
        raise TDG10QA2PREF1Error("PREF1_GIT_DRIFT", arguments) from exc


def _git_path_stream(root: Path, *arguments: str) -> tuple[str, ...]:
    raw = _git(root, *arguments)
    if not raw:
        return ()
    chunks = raw.split(b"\0")
    if chunks[-1] != b"":
        _stop("PREF1_GIT_DRIFT", "noncanonical Git path stream")
    try:
        return tuple(item.decode("utf-8") for item in chunks[:-1])
    except UnicodeDecodeError as exc:
        raise TDG10QA2PREF1Error("PREF1_GIT_DRIFT", "non-UTF-8 Git path") from exc


def _identity(metadata: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        metadata.st_nlink,
    )


def _open_directory(root: Path, relative: str) -> int:
    candidate = Path(relative)
    if candidate.is_absolute() or any(
        part in {"", ".", ".."} for part in candidate.parts
    ):
        _stop("PREF1_PATH_UNSAFE", relative)
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        _stop("PREF1_PATH_UNSAFE", root)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        if _identity(os.fstat(descriptor)) != _identity(root_before):
            _stop("PREF1_PATH_RACED", root)
        for part in candidate.parts:
            before = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _stop("PREF1_PATH_UNSAFE", relative)
            child = os.open(
                part,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=descriptor,
            )
            active = os.fstat(child)
            after = os.stat(part, dir_fd=descriptor, follow_symlinks=False)
            if _identity(before) != _identity(active) or _identity(before) != _identity(
                after
            ):
                os.close(child)
                _stop("PREF1_PATH_RACED", relative)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _read_leaf_at(directory_fd: int, name: str, *, maximum: int, label: str) -> bytes:
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


def _read_repo_leaf(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    parent = candidate.parent.as_posix()
    descriptor = (
        os.open(
            root,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        if parent == "."
        else _open_directory(root, parent)
    )
    try:
        return _read_leaf_at(
            descriptor, candidate.name, maximum=maximum, label=relative
        )
    finally:
        os.close(descriptor)


def _raw_tree(root: Path) -> dict[str, bytes]:
    descriptor = _open_directory(root, RAW_NAMESPACE)
    try:
        before = os.fstat(descriptor)
        with os.scandir(descriptor) as entries:
            names = tuple(sorted(item.name for item in entries))
        if names != ("manifest.json", "terminal.json"):
            _stop("PREF1_RAW_TREE_DRIFT", names)
        blobs = {
            name: _read_leaf_at(
                descriptor,
                name,
                maximum=_MAX_RAW_LEAF_BYTES,
                label=f"{RAW_NAMESPACE}/{name}",
            )
            for name in names
        }
        with os.scandir(descriptor) as entries:
            final_names = tuple(sorted(item.name for item in entries))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if names != final_names or _identity(before) != _identity(after):
        _stop("PREF1_RAW_TREE_CHANGED", RAW_NAMESPACE)
    observed = {name: sha256(raw).hexdigest() for name, raw in blobs.items()}
    expected = {
        "manifest.json": RAW_MANIFEST_SHA256,
        "terminal.json": RAW_TERMINAL_SHA256,
    }
    if observed != expected:
        _stop("PREF1_RAW_HASH_DRIFT", observed)
    return blobs


def _raw_snapshot(
    root: Path,
) -> tuple[dict[str, bytes], dict[str, Any], dict[str, Any]]:
    blobs = _raw_tree(root)
    return (
        blobs,
        _json(blobs["manifest.json"], "raw manifest"),
        _json(blobs["terminal.json"], "raw terminal"),
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
                        if _identity(metadata) != _identity(os.fstat(child_fd)):
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
                            (relative + "\0" + sha256(raw).hexdigest() + "\n").encode(
                                "ascii"
                            )
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


def _environment() -> dict[str, str]:
    observed = {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
    }
    if observed != ENVIRONMENT:
        _stop("PREF1_ENVIRONMENT_DRIFT", observed)
    return observed


def _width_config(width: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "retry": width["retry"],
        "prior_retry_count": width["prior_retry_count"],
        "generation": width["generation"],
        "forbidden_predecessor_generation": width["forbidden_predecessor_generation"],
        "checkpoint_sha256": width["checkpoint_sha256"],
        "checkpoint_raw_sha256": width["checkpoint_raw_sha256"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "checkpoint_journal_sha256": width["checkpoint_journal_sha256"],
        "checkpoint_journal_raw_sha256": width["checkpoint_journal_raw_sha256"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "historical_rejection_sha256": width["historical_rejection_sha256"],
        "historical_rejection_raw_sha256": width["historical_rejection_raw_sha256"],
        "attempted_width_hex": width["attempted_width_hex"],
        "transaction_sha256": width["transaction_sha256"],
        "historical_journal_sha256": width["historical_journal_sha256"],
        "failed_channels": list(width["failed_channels"]),
        "width_class": width["width_class"],
        "d12_upper_stream_sha256": width["d12_upper_stream_sha256"],
    }


def _expected_config() -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "target_protocol": TARGET_PROTOCOL,
        "owner_document": OWNER_DOCUMENT,
        "authority": {
            "authority_commit": AUTHORITY_COMMIT,
            "parent_commit": AUTHORITY_PARENT,
            "blob_count": len(_AUTHORITY_BLOBS),
            "blobs": [
                {"path": path, "sha256": digest} for path, digest in _AUTHORITY_BLOBS
            ],
        },
        "live_dependency_authority": {
            "blob_count": len(_LIVE_DEPENDENCY_BLOBS),
            "blobs": [
                {"path": path, "sha256": digest}
                for path, digest in _LIVE_DEPENDENCY_BLOBS
            ],
        },
        "live_working_boundary": {
            "authority_commit": AUTHORITY_COMMIT,
            "staged_path_count": 0,
            "pre_result_path_count": len(_PREF1_PRE_RESULT_PATHS),
            "pre_result_paths": list(_PREF1_PRE_RESULT_PATHS),
            "allowed_additional_untracked": list(_ALLOWED_ADDITIONAL_UNTRACKED),
        },
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
        "store": {
            "path": STORE_PATH,
            "leaf_count": SEALED_STORE_LEAF_COUNT,
            "snapshot_sha256": SEALED_STORE_SHA256,
            "snapshot_algorithm": "TDG9-LOC1-STORE-SNAPSHOT-v1",
        },
        "environment": dict(ENVIRONMENT),
        "selection": {
            "tableau_selector": "SSPRK3",
            "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
            "actual_spatial_operator": "inherited_RK4_2049_SBP4",
            "interval_owner": INTERVAL_OWNER,
            "classifier_owner": "classify_tdg6_channel",
            "member_key": MEMBER_KEY,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "selected_retries": [4, 5],
            "per_channel_maximum_candidates_D01": (PER_CHANNEL_MAXIMUM_CANDIDATES_D01),
            "per_channel_maximum_candidates_D12": (PER_CHANNEL_MAXIMUM_CANDIDATES_D12),
            "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
            "channel_order": list(CHANNEL_ORDER),
            "widths": [_width_config(width) for width in WIDTHS],
        },
        "work_budget": {
            "SSPRK3_records_per_proposal": 4,
            "channel_count_per_width": 18,
            "fourth_width_authorized": False,
            "maximum_stage_and_endpoint_RHS_records": 56,
            "owned_row_count": OWNED_ROW_COUNT,
            "per_channel_maximum_candidates_D01": (PER_CHANNEL_MAXIMUM_CANDIDATES_D01),
            "per_channel_maximum_candidates_D12": (PER_CHANNEL_MAXIMUM_CANDIDATES_D12),
            "primary_refinement_depth": PRIMARY_REFINEMENT_DEPTH,
            "resource_escalation_authorized": False,
            "retry_count": 2,
            "shadow_paths": 14,
            "shadow_paths_per_width": 7,
            "shadow_proposals": 14,
            "shadow_proposals_per_width": 7,
        },
        "decision": {
            "any_width_nonpass_rejects_exact_remedy_on_tested_neighborhood": True,
            "both_widths_all_18_pass_required_for_robustness": True,
            "classifier_owner": "classify_tdg6_channel",
            "complete_admission_is_all_of": True,
            "diagnostic_fine_endpoint_serialized": False,
            "equality_passes": True,
            "interval_owner": INTERVAL_OWNER,
            "production_method_earned": False,
            "state_advance_authorized": False,
            "sufficient_condition": SUFFICIENT_CONDITION,
            "successor_remedy_selected": False,
            "width_robustness_passed": False,
        },
        "scope": {
            "QA1_decision_owner_imported": False,
            "QA2_authority_imported": False,
            "QA2_runner_decision_code_imported": False,
            "TI2_authority_imported": False,
            "TI2_runner_imported": False,
            "campaign_store_mutation_authorized": False,
            "diagnostic_endpoint_serialization_authorized": False,
            "exact_complete_C_admission_imported": False,
            "exact_complete_C_runtime_imported": False,
            "fourteen_proposals_reconstructed": True,
            "lower_level_exact_arithmetic_and_dual_localizers_reused": True,
            "one_time_live_raw_and_store_authentication": True,
            "ordinary_verifier_git_blind": True,
            "ordinary_verifier_raw_blind": True,
            "ordinary_verifier_shadow_blind": True,
            "ordinary_verifier_store_blind": True,
            "production_admission_called": False,
            "production_commit_called": False,
            "raw_namespace_mutation_authorized": False,
            "state_serialization_called": False,
            "temporal_retry_admission_authorized": False,
            "two_width_predecessors_reconstructed": True,
        },
        "claims": {
            "FGCQR_authorized": False,
            "GR0_calibration_completed": False,
            "SGBL_authorized": False,
            "accepted_state": False,
            "both_widths_all_eighteen_channels_admitted": False,
            "candidate_execution_authorized": False,
            "common_event_completed": False,
            "exact_remedy_rejected_on_tested_retry_neighborhood": True,
            "external_publication_authorized": False,
            "independent_method_agreement": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
            "physical_transition_claim_authorized": False,
            "production_SSPRK3_comparator": False,
            "production_method_earned": False,
            "retained_EFT_evolution_authorized": False,
            "state_advance_authorized": False,
            "successor_remedy_selected": False,
            "width_robustness_passed": False,
        },
    }


def expected_config() -> dict[str, Any]:
    return json.loads(_canonical(_expected_config()).decode("ascii"))


def _config(raw: bytes) -> dict[str, Any]:
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise TDG10QA2PREF1Error("PREF1_CONFIG_DRIFT", exc) from exc
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
        "live_dependency_authority": dict(config["live_dependency_authority"]),
        "live_working_boundary": dict(config["live_working_boundary"]),
        "store": dict(config["store"]),
        "environment": dict(config["environment"]),
        "selection": dict(config["selection"]),
        "work_budget": dict(config["work_budget"]),
        "decision": dict(config["decision"]),
        "scope": dict(config["scope"]),
        "claims": dict(config["claims"]),
    }


def _sealed_store() -> dict[str, Any]:
    return {"leaf_count": SEALED_STORE_LEAF_COUNT, "sha256": SEALED_STORE_SHA256}


def _authenticate_authority(root: Path, manifest: Mapping[str, Any]) -> None:
    if manifest.get("authority_commit") != AUTHORITY_COMMIT:
        _stop("PREF1_AUTHORITY_DRIFT", manifest.get("authority_commit"))
    try:
        top = Path(_git(root, "rev-parse", "--show-toplevel").decode().strip())
        if not os.path.samefile(root, top):
            _stop("PREF1_AUTHORITY_DRIFT", "repository is not the worktree root")
    except (OSError, UnicodeDecodeError) as exc:
        raise TDG10QA2PREF1Error(
            "PREF1_AUTHORITY_DRIFT", "Git worktree root differs"
        ) from exc
    if _git(root, "cat-file", "-t", AUTHORITY_COMMIT).strip() != b"commit":
        _stop("PREF1_AUTHORITY_DRIFT", "authority object is not a commit")
    parents = _git(root, "show", "-s", "--format=%P", AUTHORITY_COMMIT).split()
    if parents != [AUTHORITY_PARENT.encode("ascii")]:
        _stop("PREF1_AUTHORITY_LINEAGE_DRIFT", parents)
    if _git(root, "for-each-ref", "--format=%(refname)", "refs/replace").strip():
        _stop("PREF1_AUTHORITY_DRIFT", "Git replace refs are present")
    expected_paths = tuple(path for path, _ in _AUTHORITY_BLOBS)
    if len(_AUTHORITY_BLOBS) != 16 or len(set(expected_paths)) != 16:
        _stop("PREF1_AUTHORITY_BLOB_COUNT_DRIFT", len(_AUTHORITY_BLOBS))
    observed_paths = tuple(
        sorted(
            _git(
                root,
                "diff-tree",
                "--root",
                "--no-commit-id",
                "--name-only",
                "-r",
                AUTHORITY_COMMIT,
            )
            .decode("utf-8")
            .splitlines()
        )
    )
    if observed_paths != tuple(sorted(expected_paths)):
        _stop("PREF1_AUTHORITY_DELTA_DRIFT", observed_paths)
    for relative, expected in _AUTHORITY_BLOBS:
        observed = sha256(
            _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        ).hexdigest()
        if observed != expected:
            _stop("PREF1_AUTHORITY_BLOB_DRIFT", relative)


def _authenticate_live_dependencies(root: Path) -> None:
    if len(_LIVE_DEPENDENCY_BLOBS) != 10 or len(
        {path for path, _ in _LIVE_DEPENDENCY_BLOBS}
    ) != len(_LIVE_DEPENDENCY_BLOBS):
        _stop("PREF1_DEPENDENCY_INVENTORY_DRIFT", len(_LIVE_DEPENDENCY_BLOBS))
    head = _git(root, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != AUTHORITY_COMMIT:
        _stop("PREF1_DEPENDENCY_GIT_DRIFT", head)
    paths = tuple(path for path, _ in _LIVE_DEPENDENCY_BLOBS)
    for arguments in (
        ("diff", "--name-only", "--", *paths),
        ("diff", "--cached", "--name-only", "--", *paths),
        ("ls-files", "--others", "--exclude-standard", "--", *paths),
    ):
        if _git(root, *arguments).strip():
            _stop("PREF1_DEPENDENCY_GIT_DRIFT", arguments[0])
    for relative, expected in _LIVE_DEPENDENCY_BLOBS:
        committed = _git(root, "show", f"{AUTHORITY_COMMIT}:{relative}")
        live = _read_repo_leaf(root, relative, 64 << 20)
        if (
            sha256(committed).hexdigest() != expected
            or sha256(live).hexdigest() != expected
            or live != committed
        ):
            _stop("PREF1_DEPENDENCY_BLOB_DRIFT", relative)


def _authenticate_live_working_boundary(root: Path) -> None:
    if (
        len(_PREF1_PRE_RESULT_PATHS) != 13
        or len(set(_PREF1_PRE_RESULT_PATHS)) != 13
        or RESULT_PATH in _PREF1_PRE_RESULT_PATHS
    ):
        _stop("PREF1_WORKING_BOUNDARY_DRIFT", "pre-result inventory")
    head = _git(root, "rev-parse", "--verify", "HEAD^{commit}").decode().strip()
    if head != AUTHORITY_COMMIT:
        _stop("PREF1_WORKING_BOUNDARY_DRIFT", head)
    unstaged = _git_path_stream(root, "diff", "--name-only", "-z", "--")
    staged = _git_path_stream(root, "diff", "--cached", "--name-only", "-z", "--")
    untracked = _git_path_stream(
        root, "ls-files", "--others", "--exclude-standard", "-z", "--"
    )
    task_untracked = tuple(
        item for item in untracked if item not in _ALLOWED_ADDITIONAL_UNTRACKED
    )
    combined = (*unstaged, *task_untracked)
    observed = tuple(sorted(combined))
    if (
        staged
        or len(combined) != len(set(combined))
        or observed != (_PREF1_PRE_RESULT_PATHS)
    ):
        _stop(
            "PREF1_WORKING_BOUNDARY_DRIFT",
            {"staged": staged, "observed": observed, "untracked": untracked},
        )
    for relative in _PREF1_PRE_RESULT_PATHS:
        _read_repo_leaf(root, relative, 64 << 20)


def _expected_manifest() -> dict[str, Any]:
    return {
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": QA2_ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "diagnostic_fine_endpoint_serialized": False,
        "diagnostic_qualification_only": True,
        "output_leaves": ["manifest.json", "terminal.json"],
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "selected_retries": [4, 5],
        "state_advance_authorized": False,
        "tableau_selector": "SSPRK3",
        "target_protocol": TARGET_PROTOCOL,
        "width_robustness_passed_requires_both_widths_all_18_pass": True,
    }


def _terminal_required_nonclaims() -> dict[str, Any]:
    return {
        "FGCQR_authorized": False,
        "GR0_calibration_authorized": False,
        "GR0_calibration_completed": False,
        "PDE_state_commit_authorized": False,
        "PDE_state_committed": False,
        "QA1_licenses_QA2_method_design_never_old_member_adoption": True,
        "RA1_authorized": False,
        "SGBL_authorized": False,
        "campaign_state_write_authorized": False,
        "candidate_branch_opened": False,
        "candidate_execution_authorized": False,
        "common_event_authorized": False,
        "common_event_completed": False,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "diagnostic_fine_endpoint_serialized": False,
        "diagnostic_qualification_only": True,
        "external_publication_authorized": False,
        "fine_path_commit_authorized": False,
        "fine_path_committed": False,
        "fourth_width_authorized": False,
        "independent_method_agreement_earned": False,
        "mechanism_result_authorized": False,
        "mechanism_result_earned": False,
        "new_protocol_id_authorized": False,
        "old_member_adoption_authorized": False,
        "physical_result_authorized": False,
        "physical_result_earned": False,
        "physical_transition_claim_authorized": False,
        "production_SSPRK3_comparator_earned": False,
        "production_method_earned": False,
        "push_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "retry_3_authorized": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "target_protocol": TARGET_PROTOCOL,
        "temporal_retry_admission_called": False,
    }


def _expected_terminal_base(classification: str, *, robustness: bool) -> dict[str, Any]:
    if classification not in _RAW_TERMINAL_CLASSES:
        _stop("PREF1_TERMINAL_CLASS_DRIFT", classification)
    if robustness is not (classification == BOTH_PASS_CLASS):
        _stop("PREF1_ROBUSTNESS_DRIFT", (classification, robustness))
    sealed = _sealed_store()
    return {
        **_terminal_required_nonclaims(),
        "actual_spatial_operator": "inherited_RK4_2049_SBP4",
        "artifact_id": QA2_ARTIFACT_ID,
        "authority_commit": AUTHORITY_COMMIT,
        "channel_order": list(CHANNEL_ORDER),
        "classification": classification,
        "runner_id": RUNNER_ID,
        "schema": RAW_SCHEMA,
        "selected_retries": [4, 5],
        "store_snapshot_after": sealed,
        "store_snapshot_before": sealed,
        "store_unchanged": True,
        "tableau_runtime_selector": TABLEAU_RUNTIME_SELECTOR,
        "tableau_selector": "SSPRK3",
        "width_robustness_passed": robustness,
    }


def _mapping(value: object, keys: frozenset[str], label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        _stop("PREF1_SCHEMA_DRIFT", (label, sorted(keys)))
    return value


def _digest(value: object, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in _HEX for character in value)
    ):
        _stop("PREF1_DIGEST_DRIFT", label)
    return value


def _count(value: object, *, label: str, allow_zero: bool = False) -> int:
    if type(value) is not int or value < (0 if allow_zero else 1):
        _stop("PREF1_COUNT_DRIFT", label)
    return value


def _encode_rational(value: Fraction) -> dict[str, str]:
    if not isinstance(value, Fraction) or isinstance(value, bool):
        _stop("PREF1_EXACT_TYPE_DRIFT", type(value).__name__)
    encoded = {
        "denominator": str(value.denominator),
        "numerator": str(value.numerator),
    }
    rebuilt = Fraction(int(encoded["numerator"]), int(encoded["denominator"]))
    if rebuilt != value:
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", encoded)
    return encoded


def _fraction_from_rational(value: object, *, label: str) -> Fraction:
    encoded = _mapping(value, _RATIONAL_KEYS, label)
    numerator = encoded["numerator"]
    denominator = encoded["denominator"]
    if not isinstance(numerator, str) or not isinstance(denominator, str):
        _stop("PREF1_EXACT_TYPE_DRIFT", label)
    try:
        answer = Fraction(int(numerator), int(denominator))
    except (ValueError, ZeroDivisionError) as exc:
        raise TDG10QA2PREF1Error("PREF1_EXACT_PARSE_DRIFT", label) from exc
    if str(answer.numerator) != numerator or str(answer.denominator) != denominator:
        _stop("PREF1_EXACT_ROUNDTRIP_DRIFT", label)
    return answer


def _class_name(value: object) -> str:
    kind = type(value)
    return f"{kind.__module__}.{kind.__qualname__}"


def _callable_name(value: object) -> str:
    return (
        f"{getattr(value, '__module__', type(value).__module__)}."
        f"{getattr(value, '__qualname__', type(value).__qualname__)}"
    )


def _fingerprint_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if type(value) is float:
        if not math.isfinite(value):
            _stop("PREF1_FINGERPRINT_DRIFT", "nonfinite binary64")
        return {"binary64_hex": value.hex()}
    if is_dataclass(value) and not isinstance(value, type):
        return _fingerprint_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _fingerprint_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_fingerprint_exact(item) for item in value]
    if value is None or type(value) in {str, int, bool}:
        return value
    _stop("PREF1_FINGERPRINT_DRIFT", type(value).__name__)


def _tracer_history_hash(member: object) -> str:
    arrays = [
        np.asarray(member.tracers.labels, dtype=np.float64),
        np.asarray(member.tracers.positions, dtype=np.float64),
        np.asarray(member.tracers.proper_times, dtype=np.float64),
    ]
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_proper_times
    )
    arrays.extend(
        np.asarray(row, dtype=np.float64) for row in member.tracers.event_fields
    )
    return array_content_sha256(*arrays)


def _member_fingerprint(member: object, descriptor_sha256: str) -> dict[str, object]:
    ledger = member.temporal_ledger
    if ledger is None:
        _stop("PREF1_REPLAY_DRIFT", "temporal ledger absent")
    return {
        "member_key": member.key,
        "method_label": member.method_label,
        "source_integrator": member.integrator_id,
        "point_count": member.point_count,
        "accepted_time_hex": float(member.time).hex(),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
        "coordinates_sha256": array_content_sha256(member.initial.grid.coordinates),
        "grid_spacing_hex": float(member.initial.grid.spacing).hex(),
        "outer_radius_hex": float(member.initial.grid.maximum).hex(),
        "descriptor_sha256": descriptor_sha256,
        "operator_class": _class_name(member.operator),
        "projector_callable": _callable_name(member.projector),
        "transaction_class": _class_name(member.transaction),
        "tracer_class": _class_name(member.tracers),
        "ledger_class": _class_name(ledger),
        "transaction_sha256": sha256(
            _canonical(
                _fingerprint_exact(
                    {
                        "monitor": asdict(member.transaction.state),
                        "causal": asdict(member.transaction.causal_state),
                        "ledger": asdict(ledger),
                        "step_index": member.step_index,
                        "transaction_serial": member.transaction_serial,
                    }
                )
            )
        ).hexdigest(),
        "tracer_history_sha256": _tracer_history_hash(member),
    }


def _expected_replay_receipt(width: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "member_key": MEMBER_KEY,
        "retry": width["retry"],
        "predecessor_generation": width["generation"],
        "checkpoint_journal_sequence": width["checkpoint_journal_sequence"],
        "historical_rejection_sequence": width["historical_rejection_sequence"],
        "attempted_width_hex": width["attempted_width_hex"],
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": width["transaction_sha256"],
        "historical_journal_sha256": width["historical_journal_sha256"],
        "ledger_owner": "TDG6TemporalLedger",
        "tracer_owner": "NormalFlowTracers",
        "transaction_owner": "GR0RuntimeStageTransaction",
    }


def width_spec(retry: int) -> dict[str, Any]:
    for width in WIDTHS:
        if width["retry"] == retry:
            return dict(width)
    _stop("PREF1_RETRY_DRIFT", retry)


def _journal_evidence(root: Path, width: Mapping[str, Any]) -> Mapping[str, Any]:
    relative = (
        f"{STORE_PATH}/journal/"
        f"{int(width['historical_rejection_sequence']):020d}-"
        f"{width['historical_rejection_sha256']}.journal"
    )
    raw = _read_repo_leaf(root, relative, 8 << 20)
    if sha256(raw).hexdigest() != width["historical_rejection_raw_sha256"]:
        _stop("PREF1_JOURNAL_HASH_DRIFT", relative)
    try:
        value = json.loads(raw.decode("ascii"), object_pairs_hook=_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise TDG10QA2PREF1Error("PREF1_JOURNAL_EVIDENCE_DRIFT", relative) from exc
    if (
        not isinstance(value, Mapping)
        or value.get("sequence") != width["historical_rejection_sequence"]
        or value.get("record_sha256") != width["historical_rejection_sha256"]
        or value.get("kind") != "tdg6_rejection"
    ):
        _stop("PREF1_JOURNAL_EVIDENCE_DRIFT", width["retry"])
    evidence = value.get("payload", {}).get("evidence")
    if not isinstance(evidence, Mapping):
        _stop("PREF1_JOURNAL_EVIDENCE_DRIFT", "evidence absent")
    if sha256(_canonical(evidence)).hexdigest() != width["historical_journal_sha256"]:
        _stop("PREF1_HISTORICAL_JOURNAL_DRIFT", width["retry"])
    return evidence


def _verify_checkpoint(root: Path, width: Mapping[str, Any]) -> None:
    generation = int(width["generation"])
    relative = (
        f"{STORE_PATH}/checkpoints/{generation:020d}-{width['checkpoint_sha256']}.json"
    )
    raw = _read_repo_leaf(root, relative, 8 << 20)
    if sha256(raw).hexdigest() != width["checkpoint_raw_sha256"]:
        _stop("PREF1_CHECKPOINT_RAW_DRIFT", relative)
    value = json.loads(raw.decode("ascii"), object_pairs_hook=_pairs)
    if (
        not isinstance(value, Mapping)
        or value.get("generation") != generation
        or generation == width["forbidden_predecessor_generation"]
        or value.get("checkpoint_sha256") != width["checkpoint_sha256"]
        or value.get("journal_sequence") != width["checkpoint_journal_sequence"]
        or value.get("journal_tip_sha256") != width["checkpoint_journal_sha256"]
    ):
        _stop("PREF1_CHECKPOINT_DRIFT", width["retry"])
    journal_relative = (
        f"{STORE_PATH}/journal/"
        f"{int(width['checkpoint_journal_sequence']):020d}-"
        f"{width['checkpoint_journal_sha256']}.journal"
    )
    journal_raw = _read_repo_leaf(root, journal_relative, 8 << 20)
    if sha256(journal_raw).hexdigest() != width["checkpoint_journal_raw_sha256"]:
        _stop("PREF1_CHECKPOINT_JOURNAL_DRIFT", journal_relative)
    journal = json.loads(journal_raw.decode("ascii"), object_pairs_hook=_pairs)
    if (
        not isinstance(journal, Mapping)
        or journal.get("sequence") != width["checkpoint_journal_sequence"]
        or journal.get("record_sha256") != width["checkpoint_journal_sha256"]
        or journal.get("kind") == "tdg6_rejection"
    ):
        _stop("PREF1_CHECKPOINT_JOURNAL_DRIFT", width["retry"])


def _restore_shadow(
    root: Path, width: Mapping[str, Any]
) -> tuple[tdg6.TDG6PreparedGR0Compositor, dict[str, Any]]:
    generation = int(width["generation"])
    if generation == 9 or generation == width["forbidden_predecessor_generation"]:
        _stop("PREF1_GENERATION_DRIFT", (width["retry"], generation))
    _verify_checkpoint(root, width)
    _journal_evidence(root, width)
    store = HLT16CampaignStore(root / STORE_PATH)
    shells = build_static_gr0_shells(root)
    checkpoint = store.authenticated_checkpoint_at_generation(generation)
    if (
        checkpoint.sha256 != width["checkpoint_sha256"]
        or checkpoint.generation != generation
        or checkpoint.journal_sequence != width["checkpoint_journal_sequence"]
    ):
        _stop("PREF1_CHECKPOINT_DRIFT", width["retry"])
    state = checkpoint.members.get(MEMBER_KEY)
    if state is None:
        _stop("PREF1_MEMBER_ABSENT", MEMBER_KEY)
    member = deepcopy(shells[MEMBER_KEY])
    try:
        campaign_runtime.restore_member_with_overlay(
            store, checkpoint, member, key=MEMBER_KEY
        )
        cursor = p15.Proto15Cursor(dict(state.cursor))
        cursor.validate()
    except Exception as exc:
        raise TDG10QA2PREF1Error("PREF1_RESTORE_DRIFT", width["retry"]) from exc
    ledger = member.temporal_ledger
    if (
        ledger is None
        or member.key != MEMBER_KEY
        or member.method_label != "RK4"
        or member.integrator_id != PRIMARY_METHOD
        or member.point_count != POINT_COUNT
        or float(member.time).hex() != ACCEPTED_TIME_HEX
        or ledger.current_macro_step_temporal_retry_count != width["prior_retry_count"]
        or cursor.mode != "RETRY_PENDING"
        or state.pending_owner != "temporal"
    ):
        _stop("PREF1_PREDECESSOR_DRIFT", width["retry"])
    fingerprint = _member_fingerprint(member, state.descriptor_sha256)
    fixed = {
        "member_key": MEMBER_KEY,
        "method_label": "RK4",
        "source_integrator": PRIMARY_METHOD,
        "point_count": POINT_COUNT,
        "accepted_time_hex": ACCEPTED_TIME_HEX,
        "state_sha256": PHYSICAL_STATE_SHA256,
        "coordinates_sha256": COORDINATES_SHA256,
        "grid_spacing_hex": GRID_SPACING_HEX,
        "outer_radius_hex": OUTER_RADIUS_HEX,
        "descriptor_sha256": MEMBER_DESCRIPTOR_SHA256,
        "transaction_sha256": width["transaction_sha256"],
    }
    if any(fingerprint.get(name) != expected for name, expected in fixed.items()):
        _stop("PREF1_FINGERPRINT_DRIFT", width["retry"])
    before = dict(fingerprint)
    try:
        prepared = tdg6.prepare_tdg6_gr0_compositor(
            method=COMPARATOR_METHOD,
            time=member.time,
            step_size=float.fromhex(str(width["attempted_width_hex"])),
            state=member.state,
            rhs=member.operator,
            projector=member.projector,
            transaction=member.transaction,
            tracers=member.tracers,
            coordinates=member.initial.grid.coordinates,
            temporal_ledger=ledger,
            previous_step_index=member.step_index,
            previous_transaction_serial=member.transaction_serial,
        )
    except Exception as exc:
        raise TDG10QA2PREF1Error("PREF1_SHADOW_PREMISE_STOP", width["retry"]) from exc
    after = _member_fingerprint(member, state.descriptor_sha256)
    paths = (prepared.outer, prepared.medium, prepared.fine)
    proposals = tuple(attempt.proposal for path in paths for attempt in path.attempts)
    if (
        after != before
        or prepared.method != COMPARATOR_METHOD
        or prepared.initial_state_sha256 != PHYSICAL_STATE_SHA256
        or tuple(len(path.attempts) for path in paths) != (1, 2, 4)
        or len(proposals) != 7
        or any(proposal.method != COMPARATOR_METHOD for proposal in proposals)
        or any(len(proposal.stages) != 4 for proposal in proposals)
    ):
        _stop("PREF1_ONE_VARIABLE_IDENTITY_DRIFT", width["retry"])
    receipt = {
        **_expected_replay_receipt(width),
        "accepted_time_hex": fingerprint["accepted_time_hex"],
        "state_sha256": fingerprint["state_sha256"],
        "descriptor_sha256": fingerprint["descriptor_sha256"],
        "transaction_sha256": fingerprint["transaction_sha256"],
    }
    expected_receipt = _expected_replay_receipt(width)
    if receipt != expected_receipt:
        _stop("PREF1_REPLAY_RECEIPT_DRIFT", receipt)
    return prepared, expected_receipt


@dataclass(frozen=True, slots=True)
class _IndependentDifference:
    lower: Fraction
    upper: Fraction
    polynomial_count: int
    candidate_count: int
    co_maximizer_count: int
    localization_classification: str
    coefficient_stream_sha256: str
    survivor_key_stream_sha256: str
    primary_stationary_count_stream_sha256: str
    independent_stationary_count_stream_sha256: str
    primary_evaluator_id: str
    independent_evaluator_id: str
    maximum_candidates: int
    routes_agree: bool = True


@dataclass(frozen=True, slots=True)
class _IndependentAdmission:
    d01: _IndependentDifference
    d12: _IndependentDifference
    decision: object
    row_count: int
    row_stream_sha256: str
    combined_coefficient_stream_sha256: str
    sufficient_pass_left: Fraction
    sufficient_pass_right: Fraction
    sufficient_contraction_pass: bool
    sufficient_contraction_failure: bool
    threshold_inconclusive: bool
    maximum_candidates_D01: int
    maximum_candidates_D12: int
    refinement_depth: int
    evaluator_id: str = TDG10_EVALUATOR_ID


@dataclass(frozen=True, slots=True)
class _IndependentChannel:
    channel: str
    evidence: _IndependentAdmission


@dataclass(frozen=True, slots=True)
class _IndependentAssessment:
    channel_evidence: tuple[_IndependentChannel, ...]
    complete_admission_passed: bool
    failed_channels: tuple[str, ...]


def _binary64(value: object, *, label: str) -> float:
    if type(value) is float:
        answer = value
    elif type(value) is np.float64:
        answer = float(value)
    else:
        _stop("PREF1_BINARY64_DRIFT", label)
    if not math.isfinite(answer):
        _stop("PREF1_BINARY64_DRIFT", f"{label} is nonfinite")
    return answer


def _same_binary64(left: object, right: object) -> bool:
    first = _binary64(left, label="left binary64")
    second = _binary64(right, label="right binary64")
    return np.float64(first).tobytes() == np.float64(second).tobytes()


def _state_sha256(state: object) -> str:
    try:
        return array_content_sha256(state.u, state.p, state.q)
    except Exception as exc:
        raise TDG10QA2PREF1Error("PREF1_STATE_DRIFT", type(state).__name__) from exc


def _rhs_sha256(rhs: object) -> str:
    try:
        return array_content_sha256(rhs.du, rhs.dp, rhs.dq)
    except Exception as exc:
        raise TDG10QA2PREF1Error("PREF1_RHS_DRIFT", type(rhs).__name__) from exc


def _owned_state(state: object) -> np.ndarray:
    arrays = tuple(np.asarray(getattr(state, name)) for name in ("u", "p", "q"))
    if any(
        value.dtype != np.dtype("<f8")
        or value.shape != (POINT_COUNT, 6)
        or not value.flags.c_contiguous
        or not np.isfinite(value).all()
        for value in arrays
    ):
        _stop("PREF1_OWNED_STATE_DRIFT", [value.shape for value in arrays])
    answer = np.stack(tuple(value[1:-4] for value in arrays), axis=0)
    if answer.shape != (3, OWNED_ROW_COUNT, 6):
        _stop("PREF1_OWNED_STATE_DRIFT", answer.shape)
    return answer


def _owned_rhs(rhs: object) -> np.ndarray:
    arrays = tuple(np.asarray(getattr(rhs, name)) for name in ("du", "dp", "dq"))
    if any(
        value.dtype != np.dtype("<f8")
        or value.shape != (POINT_COUNT, 6)
        or not value.flags.c_contiguous
        or not np.isfinite(value).all()
        for value in arrays
    ):
        _stop("PREF1_OWNED_RHS_DRIFT", [value.shape for value in arrays])
    answer = np.stack(tuple(value[1:-4] for value in arrays), axis=0)
    if answer.shape != (3, OWNED_ROW_COUNT, 6):
        _stop("PREF1_OWNED_RHS_DRIFT", answer.shape)
    return answer


def _endpoint_rhs(proposal: object) -> tuple[object, object]:
    stages = tuple(getattr(proposal, "stages", ()))
    if len(stages) != 4 or tuple(
        getattr(item, "stage_name", None) for item in stages
    ) != ("ssprk3_s0", "ssprk3_s1", "ssprk3_s2", "candidate_endpoint"):
        _stop("PREF1_STAGE_RECORD_DRIFT", len(stages))
    first, last = stages[0], stages[-1]
    if (
        not _same_binary64(first.time, proposal.initial_time)
        or not _same_binary64(last.time, proposal.final_time)
        or _state_sha256(first.state) != _state_sha256(proposal.initial_state)
        or _state_sha256(last.state) != _state_sha256(proposal.candidate_state)
    ):
        _stop("PREF1_ENDPOINT_RHS_DRIFT", proposal.method)
    return first.rhs, last.rhs


def _validate_prepared_shadow(
    prepared: object, width_specification: Mapping[str, Any]
) -> float:
    if getattr(prepared, "method", None) != COMPARATOR_METHOD:
        _stop("PREF1_TABLEAU_DRIFT", getattr(prepared, "method", None))
    start = _binary64(prepared.initial_time, label="prepared.initial_time")
    final = _binary64(prepared.final_time, label="prepared.final_time")
    width = final - start
    if (
        not math.isfinite(width)
        or width <= 0.0
        or not _same_binary64(
            width, float.fromhex(width_specification["attempted_width_hex"])
        )
    ):
        _stop("PREF1_SHADOW_WIDTH_DRIFT", width)
    if prepared.initial_state_sha256 != PHYSICAL_STATE_SHA256:
        _stop("PREF1_INITIAL_STATE_DRIFT", prepared.initial_state_sha256)
    initial_rhs: list[str] = []
    for path, level, count in (
        (prepared.outer, "outer", 1),
        (prepared.medium, "medium", 2),
        (prepared.fine, "fine", 4),
    ):
        attempts = tuple(path.attempts)
        if (
            path.level != level
            or path.method != COMPARATOR_METHOD
            or len(attempts) != count
            or path.initial_state_sha256 != PHYSICAL_STATE_SHA256
            or not _same_binary64(path.initial_time, start)
            or not _same_binary64(path.final_time, final)
        ):
            _stop("PREF1_SHADOW_PATH_DRIFT", level)
        step = width / count
        boundaries = tuple(start + index * step for index in range(count + 1))
        if not _same_binary64(boundaries[-1], final):
            _stop("PREF1_SHADOW_BOUNDARY_DRIFT", level)
        previous_state: str | None = None
        previous_rhs: str | None = None
        for index, (attempt, left, right) in enumerate(
            zip(attempts, boundaries[:-1], boundaries[1:], strict=True)
        ):
            if attempt.accepted is None or attempt.retry is not None:
                _stop("PREF1_SHADOW_ATTEMPT_DRIFT", (level, index))
            proposal = attempt.proposal
            left_rhs, right_rhs = _endpoint_rhs(proposal)
            if (
                proposal.method != COMPARATOR_METHOD
                or not _same_binary64(proposal.initial_time, left)
                or not _same_binary64(proposal.final_time, right)
                or not _same_binary64(right - left, step)
            ):
                _stop("PREF1_SHADOW_BOUNDARY_DRIFT", (level, index))
            state_digest = _state_sha256(proposal.initial_state)
            if index == 0:
                if state_digest != PHYSICAL_STATE_SHA256:
                    _stop("PREF1_SHADOW_INITIAL_DRIFT", level)
                initial_rhs.append(_rhs_sha256(left_rhs))
            elif (
                state_digest != previous_state or _rhs_sha256(left_rhs) != previous_rhs
            ):
                _stop("PREF1_SHADOW_CONTIGUITY_DRIFT", (level, index))
            previous_state = _state_sha256(proposal.candidate_state)
            previous_rhs = _rhs_sha256(right_rhs)
            if _state_sha256(attempt.accepted.state) != previous_state:
                _stop("PREF1_SHADOW_ACCEPTED_DRIFT", (level, index))
    if len(set(initial_rhs)) != 1:
        _stop("PREF1_SHARED_INITIAL_RHS_DRIFT", initial_rhs)
    _owned_state(prepared.outer.attempts[0].proposal.initial_state)
    return width


def _segment(
    proposal: object,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    left_rhs, right_rhs = _endpoint_rhs(proposal)
    width = _binary64(proposal.final_time, label="proposal.final_time") - _binary64(
        proposal.initial_time, label="proposal.initial_time"
    )
    if not math.isfinite(width) or width <= 0.0:
        _stop("PREF1_SEGMENT_WIDTH_DRIFT", width)
    return (
        _owned_state(proposal.initial_state),
        _owned_rhs(left_rhs),
        _owned_state(proposal.candidate_state),
        _owned_rhs(right_rhs),
        width,
    )


def _path_segments(
    path: object,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float], ...]:
    return tuple(_segment(attempt.proposal) for attempt in path.attempts)


def _channel_rows(
    prepared: object, channel: str, *, width: float
) -> tuple[tuple[object, object, object], ...]:
    outer = _path_segments(prepared.outer)
    medium = _path_segments(prepared.medium)
    fine = _path_segments(prepared.fine)
    if len(outer) != 1 or len(medium) != 2 or len(fine) != 4:
        _stop("PREF1_SHADOW_COUNT_DRIFT", channel)
    if (
        not _same_binary64(outer[0][4], width)
        or any(not _same_binary64(item[4], width / 2.0) for item in medium)
        or any(not _same_binary64(item[4], width / 4.0) for item in fine)
    ):
        _stop("PREF1_SEGMENT_WIDTH_DRIFT", channel)
    block_name, field_name = channel.split(":", 1)
    block = ("u", "p", "q").index(block_name)
    field = ("alpha", "v", "lambda", "R", "phi", "chi").index(field_name)

    def scalar(
        item: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float], row: int
    ) -> tuple[float, ...]:
        return (
            _binary64(item[0][block, row, field], label=f"{channel}.left"),
            _binary64(item[1][block, row, field], label=f"{channel}.left_rhs"),
            _binary64(item[2][block, row, field], label=f"{channel}.right"),
            _binary64(item[3][block, row, field], label=f"{channel}.right_rhs"),
            _binary64(item[4], label=f"{channel}.width"),
        )

    return tuple(
        (
            scalar(outer[0], row),
            tuple(scalar(item, row) for item in medium),
            tuple(scalar(item, row) for item in fine),
        )
        for row in range(OWNED_ROW_COUNT)
    )


def _row_line(ordinal: int, row: Sequence[object]) -> bytes:
    outer, medium, fine = row
    values: list[float] = []
    for segment in (outer, *medium, *fine):
        values.extend(_binary64(item, label=f"row[{ordinal}]") for item in segment)
    return (f"{ordinal}|" + "|".join(item.hex() for item in values) + "\n").encode(
        "ascii"
    )


def _feed_coefficients(hasher: object, ordinal: int, cubic: Sequence[Fraction]) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(f"{item.numerator}/{item.denominator}" for item in cubic)
        + "\n"
    )
    hasher.update(line.encode("ascii"))


def _metadata(level: str, row: int, subinterval: int) -> dict[str, object]:
    return {"level": level, "row_index": row, "subinterval": subinterval}


def _survivor_key(item: object) -> tuple[object, ...]:
    return (
        item.polynomial_ordinal,
        item.location,
        item.location_ordinal,
        tuple(item.metadata),
    )


def _survivor_digest(items: Sequence[object]) -> str:
    answer = sha256(_SURVIVOR_KEY_HASH_DOMAIN)
    for ordinal, location, location_ordinal, metadata in sorted(
        _survivor_key(item) for item in items
    ):
        encoded = ",".join(f"{key}={value}" for key, value in metadata)
        answer.update(
            f"{ordinal}|{location}|{location_ordinal}|{encoded}\n".encode("ascii")
        )
    return answer.hexdigest()


def _stationary_digest(cubics: Sequence[object]) -> str:
    answer = sha256(_STATIONARY_COUNT_DOMAIN)
    for ordinal, cubic in enumerate(cubics):
        roots = _primary_stationary_intervals(
            cubic.coefficients, depth=PRIMARY_REFINEMENT_DEPTH
        )
        answer.update(f"{ordinal}|{len(roots)}\n".encode("ascii"))
    return answer.hexdigest()


def _overlap(a: Fraction, b: Fraction, c: Fraction, d: Fraction) -> bool:
    return max(a, c) <= min(b, d)


def _candidate_overlap(first: object, second: object) -> bool:
    return _overlap(
        first.parameter_lower,
        first.parameter_upper,
        second.parameter_lower,
        second.parameter_upper,
    ) and _overlap(
        first.absolute_lower,
        first.absolute_upper,
        second.absolute_lower,
        second.absolute_upper,
    )


def _localize_independently(
    coefficients: Sequence[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ],
    *,
    level: str,
    coefficient_sha256: str,
    ceiling: int,
) -> _IndependentDifference:
    primary_cubics = tuple(
        LocalCubic(cubic, metadata) for cubic, metadata in coefficients
    )
    second_cubics = tuple(
        IndependentLocalCubicV2(cubic, metadata) for cubic, metadata in coefficients
    )
    try:
        primary = localize_absolute_maximum(
            primary_cubics,
            maximum_candidates=ceiling,
            refinement_depth=PRIMARY_REFINEMENT_DEPTH,
        )
        second = localize_absolute_maximum_independently_v2(
            second_cubics, maximum_candidates=ceiling
        )
    except RootIsolationInconclusive as exc:
        raise TDG10QA2PREF1Error(
            "PREF1_EXACT_RESOURCE_INCONCLUSIVE", f"{level}:{exc.reason}"
        ) from exc
    except RuntimeError as exc:
        raise TDG10QA2PREF1Error(
            "PREF1_EXACT_RESOURCE_INCONCLUSIVE", f"{level}:{exc}"
        ) from exc
    primary_stationary = _stationary_digest(primary_cubics)
    independent_stationary = second.stationary_count_stream_sha256
    first_items = {_survivor_key(item): item for item in primary.candidates}
    second_items = {_survivor_key(item): item for item in second.candidates}
    keys_equal = (
        len(first_items) == len(primary.candidates)
        and len(second_items) == len(second.candidates)
        and tuple(sorted(first_items)) == tuple(sorted(second_items))
    )
    route_facts = {
        "polynomial_count": primary.polynomial_count == second.polynomial_count,
        "candidate_count": primary.candidate_count == second.candidate_count,
        "classification": primary.classification == second.classification,
        "stationary_count": primary_stationary == independent_stationary,
        "survivor_keys": keys_equal,
        "survivor_intervals": keys_equal
        and all(
            _candidate_overlap(first_items[key], second_items[key])
            for key in first_items
        ),
        "global_interval": _overlap(
            primary.global_absolute_lower,
            primary.global_absolute_upper,
            second.global_absolute_lower,
            second.global_absolute_upper,
        ),
        "primary_tolerance_absent": primary.tolerance_used is False,
        "independent_tolerance_absent": second.tolerance_used is False,
        "independent_boundary_clipping_absent": (
            second.boundary_clipping_used is False
        ),
    }
    if not all(route_facts.values()):
        _stop("PREF1_EXACT_ROUTE_DISAGREEMENT", (level, route_facts))
    lower = max(primary.global_absolute_lower, second.global_absolute_lower)
    upper = min(primary.global_absolute_upper, second.global_absolute_upper)
    if lower < 0 or upper < lower:
        _stop("PREF1_EXACT_INTERVAL_DRIFT", level)
    return _IndependentDifference(
        lower=lower,
        upper=upper,
        polynomial_count=primary.polynomial_count,
        candidate_count=primary.candidate_count,
        co_maximizer_count=len(primary.candidates),
        localization_classification=primary.classification,
        coefficient_stream_sha256=coefficient_sha256,
        survivor_key_stream_sha256=_survivor_digest(primary.candidates),
        primary_stationary_count_stream_sha256=primary_stationary,
        independent_stationary_count_stream_sha256=independent_stationary,
        primary_evaluator_id=PRIMARY_LOCALIZER_ID,
        independent_evaluator_id=INDEPENDENT_LOCALIZER_ID,
        maximum_candidates=ceiling,
    )


def _assess_rows_independently(rows: Sequence[object]) -> _IndependentAdmission:
    row_hash = sha256(_ROW_HASH_DOMAIN)
    d01_hash = sha256(_COEFFICIENT_HASH_DOMAIN)
    d12_hash = sha256(_COEFFICIENT_HASH_DOMAIN)
    d01_coefficients: list[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ] = []
    d12_coefficients: list[
        tuple[tuple[Fraction, Fraction, Fraction, Fraction], dict[str, object]]
    ] = []
    for row_index, row in enumerate(rows):
        outer, medium, fine = row
        if len(outer) != 5 or len(medium) != 2 or len(fine) != 4:
            _stop("PREF1_ROW_SCHEMA_DRIFT", row_index)
        row_hash.update(_row_line(row_index, row))
        outer_cubic = exact_hermite_coefficients(outer)
        medium_cubics = tuple(exact_hermite_coefficients(item) for item in medium)
        fine_cubics = tuple(exact_hermite_coefficients(item) for item in fine)
        outer_width = Fraction(
            *_binary64(outer[4], label="outer width").as_integer_ratio()
        )
        if any(
            Fraction(*_binary64(item[4], label="medium width").as_integer_ratio())
            != outer_width / 2
            for item in medium
        ) or any(
            Fraction(*_binary64(item[4], label="fine width").as_integer_ratio())
            != outer_width / 4
            for item in fine
        ):
            _stop("PREF1_ROW_WIDTH_DRIFT", row_index)
        for half, child in enumerate(medium_cubics):
            cubic = subtract_exact_cubics(
                restrict_exact_cubic_to_half(outer_cubic, half=half), child
            )
            _feed_coefficients(d01_hash, len(d01_coefficients), cubic)
            d01_coefficients.append((cubic, _metadata("D01", row_index, half)))
        for parent_index, parent in enumerate(medium_cubics):
            for half in (0, 1):
                child_index = 2 * parent_index + half
                cubic = subtract_exact_cubics(
                    restrict_exact_cubic_to_half(parent, half=half),
                    fine_cubics[child_index],
                )
                _feed_coefficients(d12_hash, len(d12_coefficients), cubic)
                d12_coefficients.append(
                    (cubic, _metadata("D12", row_index, child_index))
                )
    if len(rows) != OWNED_ROW_COUNT:
        _stop("PREF1_ROW_COUNT_DRIFT", len(rows))
    d01 = _localize_independently(
        d01_coefficients,
        level="D01",
        coefficient_sha256=d01_hash.hexdigest(),
        ceiling=PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
    )
    d12 = _localize_independently(
        d12_coefficients,
        level="D12",
        coefficient_sha256=d12_hash.hexdigest(),
        ceiling=PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
    )
    decision = classify_tdg6_channel(
        CertifiedMagnitudeInterval(d01.lower, d01.upper),
        CertifiedMagnitudeInterval(d12.lower, d12.upper),
    )
    left = TDG6_ORDER_SQUARED_MULTIPLIER * d12.upper**2
    right = d01.lower**2
    exact_zero = d01.upper == 0 and d12.upper == 0
    sufficient_pass = (not exact_zero) and left <= right
    sufficient_failure = (
        (not exact_zero)
        and not sufficient_pass
        and TDG6_ORDER_SQUARED_MULTIPLIER * d12.lower**2 > d01.upper**2
    )
    combined = sha256(
        _COMBINED_HASH_DOMAIN
        + f"{d01.coefficient_stream_sha256}\n{d12.coefficient_stream_sha256}\n".encode(
            "ascii"
        )
    ).hexdigest()
    return _IndependentAdmission(
        d01=d01,
        d12=d12,
        decision=decision,
        row_count=len(rows),
        row_stream_sha256=row_hash.hexdigest(),
        combined_coefficient_stream_sha256=combined,
        sufficient_pass_left=left,
        sufficient_pass_right=right,
        sufficient_contraction_pass=sufficient_pass,
        sufficient_contraction_failure=sufficient_failure,
        threshold_inconclusive=(
            not exact_zero and not sufficient_pass and not sufficient_failure
        ),
        maximum_candidates_D01=PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        maximum_candidates_D12=PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        refinement_depth=PRIMARY_REFINEMENT_DEPTH,
    )


def _assess_prepared_independently(
    prepared: object, width_specification: Mapping[str, Any]
) -> _IndependentAssessment:
    width = _validate_prepared_shadow(prepared, width_specification)
    channels = tuple(
        _IndependentChannel(
            channel=channel,
            evidence=_assess_rows_independently(
                _channel_rows(prepared, channel, width=width)
            ),
        )
        for channel in CHANNEL_ORDER
    )
    failed = tuple(
        item.channel for item in channels if not item.evidence.decision.admission_passed
    )
    return _IndependentAssessment(
        channel_evidence=channels,
        complete_admission_passed=not failed,
        failed_channels=failed,
    )


def _difference_payload(evidence: object, *, level: str) -> dict[str, Any]:
    maximum = (
        PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        if level == "D01"
        else PER_CHANNEL_MAXIMUM_CANDIDATES_D12
    )
    polynomial_count = D01_POLYNOMIAL_COUNT if level == "D01" else D12_POLYNOMIAL_COUNT
    payload = {
        "candidate_count": _count(
            evidence.candidate_count, label=f"{level}.candidate_count"
        ),
        "co_maximizer_count": _count(
            evidence.co_maximizer_count, label=f"{level}.co_maximizer_count"
        ),
        "coefficient_stream_sha256": _digest(
            evidence.coefficient_stream_sha256,
            label=f"{level}.coefficient_stream_sha256",
        ),
        "independent_evaluator_id": evidence.independent_evaluator_id,
        "independent_stationary_count_stream_sha256": _digest(
            evidence.independent_stationary_count_stream_sha256,
            label=f"{level}.independent_stationary",
        ),
        "localization_classification": evidence.localization_classification,
        "lower": _encode_rational(evidence.lower),
        "maximum_candidates": evidence.maximum_candidates,
        "polynomial_count": evidence.polynomial_count,
        "primary_evaluator_id": evidence.primary_evaluator_id,
        "primary_stationary_count_stream_sha256": _digest(
            evidence.primary_stationary_count_stream_sha256,
            label=f"{level}.primary_stationary",
        ),
        "routes_agree": evidence.routes_agree,
        "survivor_key_stream_sha256": _digest(
            evidence.survivor_key_stream_sha256,
            label=f"{level}.survivor_key_stream_sha256",
        ),
        "upper": _encode_rational(evidence.upper),
    }
    if (
        payload["maximum_candidates"] != maximum
        or payload["polynomial_count"] != polynomial_count
        or payload["routes_agree"] is not True
        or payload["primary_evaluator_id"] != PRIMARY_EVALUATOR_ID
        or payload["independent_evaluator_id"] != INDEPENDENT_EVALUATOR_ID
        or payload["localization_classification"] not in _LOCALIZATION_CLASSES
        or payload["candidate_count"] > maximum
        or payload["co_maximizer_count"] > payload["candidate_count"]
    ):
        _stop("PREF1_DIFFERENCE_CONTRACT_DRIFT", level)
    return payload


def _reclassify_channel(
    *,
    channel: str,
    d01_lower: Fraction,
    d01_upper: Fraction,
    d12_lower: Fraction,
    d12_upper: Fraction,
) -> dict[str, Any]:
    try:
        outer = CertifiedMagnitudeInterval(d01_lower, d01_upper)
        finest = CertifiedMagnitudeInterval(d12_lower, d12_upper)
        decision = classify_tdg6_channel(outer, finest)
    except (TypeError, ValueError) as exc:
        raise TDG10QA2PREF1Error("PREF1_RECLASSIFY_DRIFT", channel) from exc
    left = TDG6_ORDER_SQUARED_MULTIPLIER * d12_upper**2
    right = d01_lower**2
    exact_zero = d01_upper == 0 and d12_upper == 0
    sufficient_holds = left <= right
    sufficient_pass = (not exact_zero) and sufficient_holds
    sufficient_failure = (
        (not exact_zero)
        and not sufficient_pass
        and TDG6_ORDER_SQUARED_MULTIPLIER * d12_lower**2 > d01_upper**2
    )
    threshold_inconclusive = (
        not exact_zero and not sufficient_pass and not sufficient_failure
    )
    if d01_lower > 0:
        order_pass = three_halves_order_passes_squared(
            outer_lower_squared=d01_lower**2,
            finest_upper_squared=d12_upper**2,
        )
        if order_pass is not sufficient_holds:
            _stop("PREF1_SUFFICIENT_IDENTITY_DRIFT", channel)
    return {
        "decision": decision,
        "left": left,
        "right": right,
        "sufficient_condition_holds": sufficient_holds,
        "sufficient_contraction_pass": sufficient_pass,
        "sufficient_contraction_failure": sufficient_failure,
        "threshold_inconclusive": threshold_inconclusive,
    }


def _channel_replay(item: object) -> dict[str, Any]:
    channel = item.channel
    evidence = item.evidence
    decision = evidence.decision
    classified = _reclassify_channel(
        channel=channel,
        d01_lower=evidence.d01.lower,
        d01_upper=evidence.d01.upper,
        d12_lower=evidence.d12.lower,
        d12_upper=evidence.d12.upper,
    )
    if (
        classified["decision"] != decision
        or classified["left"] != evidence.sufficient_pass_left
        or classified["right"] != evidence.sufficient_pass_right
        or classified["sufficient_contraction_pass"]
        is not evidence.sufficient_contraction_pass
        or classified["sufficient_contraction_failure"]
        is not evidence.sufficient_contraction_failure
        or classified["threshold_inconclusive"] is not evidence.threshold_inconclusive
        or evidence.row_count != OWNED_ROW_COUNT
        or evidence.maximum_candidates_D01 != PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        or evidence.maximum_candidates_D12 != PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        or evidence.refinement_depth != PRIMARY_REFINEMENT_DEPTH
        or evidence.evaluator_id != TDG10_EVALUATOR_ID
    ):
        _stop("PREF1_RECLASSIFY_DRIFT", channel)
    return {
        "D01": _difference_payload(evidence.d01, level="D01"),
        "D12": _difference_payload(evidence.d12, level="D12"),
        "admission_passed": decision.admission_passed,
        "candidate_evidence_available": True,
        "channel": channel,
        "classification": decision.classification,
        "combined_coefficient_stream_sha256": _digest(
            evidence.combined_coefficient_stream_sha256,
            label=f"{channel}.combined",
        ),
        "evaluator_id": TDG10_EVALUATOR_ID,
        "hash_evidence_available": True,
        "maximum_candidates_D01": PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        "maximum_candidates_D12": PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        "order_threshold_passed": decision.order_threshold_passed,
        "order_threshold_resolved": decision.order_threshold_resolved,
        "refinement_depth": PRIMARY_REFINEMENT_DEPTH,
        "route_evidence_available": True,
        "row_count": OWNED_ROW_COUNT,
        "row_stream_sha256": _digest(
            evidence.row_stream_sha256, label=f"{channel}.row_stream"
        ),
        "sufficient_condition": SUFFICIENT_CONDITION,
        "sufficient_condition_holds": classified["sufficient_condition_holds"],
        "sufficient_contraction_failure": evidence.sufficient_contraction_failure,
        "sufficient_contraction_pass": evidence.sufficient_contraction_pass,
        "sufficient_pass_left": _encode_rational(classified["left"]),
        "sufficient_pass_right": _encode_rational(classified["right"]),
        "temporal_retry_permitted": decision.temporal_retry_permitted,
        "threshold_inconclusive": evidence.threshold_inconclusive,
    }


def _validate_raw_exact(value: object, *, retry: int) -> Mapping[str, Any]:
    keys = frozenset(
        {
            "admission_is_all_of",
            "channel_count",
            "channel_order",
            "channels",
            "complete_admission_passed",
            "failed_channels",
            "independent_route_required",
            "interval_owner",
        }
    )
    exact = _mapping(value, keys, f"retry-{retry}.exact_complete_C")
    channels = exact.get("channels")
    if not isinstance(channels, list) or len(channels) != len(CHANNEL_ORDER):
        _stop("PREF1_CHANNEL_COUNT_DRIFT", retry)
    failed: list[str] = []
    for channel, name in zip(channels, CHANNEL_ORDER, strict=True):
        item = _mapping(
            channel,
            frozenset({"admission_passed", "channel", "classification", "d12_upper"}),
            f"retry-{retry}.{name}",
        )
        if item.get("channel") != name:
            _stop("PREF1_CHANNEL_ORDER_DRIFT", (retry, name, item.get("channel")))
        passed = item.get("admission_passed") is True
        expected_class = "resolved_order_pass" if passed else "resolved_order_failure"
        if item.get("classification") != expected_class:
            _stop("PREF1_CHANNEL_CLASS_DRIFT", (retry, name))
        upper = _fraction_from_rational(
            item.get("d12_upper"), label=f"retry-{retry}.{name}.d12_upper"
        )
        if upper < 0:
            _stop("PREF1_INTERVAL_DRIFT", (retry, name))
        if not passed:
            failed.append(name)
    complete = not failed
    if (
        exact.get("admission_is_all_of") is not True
        or exact.get("channel_count") != 18
        or exact.get("channel_order") != list(CHANNEL_ORDER)
        or exact.get("independent_route_required") is not True
        or exact.get("interval_owner") != INTERVAL_OWNER
        or exact.get("complete_admission_passed") is not complete
        or exact.get("failed_channels") != failed
    ):
        _stop("PREF1_ALL_OF_DRIFT", (retry, failed))
    return exact


def _validate_typed_stop(value: object, *, retry: int) -> Mapping[str, Any]:
    stop = _mapping(value, frozenset({"type", "detail"}), f"retry-{retry}.stop")
    if stop["type"] != "TDG6RefinementPathStop":
        _stop("PREF1_STOP_DRIFT", retry)
    if not isinstance(stop["detail"], str) or not stop["detail"]:
        _stop("PREF1_STOP_DRIFT", retry)
    return stop


def _validate_typed_inconclusive(value: object, *, retry: int) -> Mapping[str, Any]:
    item = _mapping(
        value,
        frozenset({"type", "reason", "detail"}),
        f"retry-{retry}.inconclusive",
    )
    if (
        item.get("type")
        not in {
            "ExactCompleteCResourceExhausted",
            "ExactCompleteCRouteDisagreement",
            "TDG10ExactCompleteCRuntimeClosed",
        }
        or not isinstance(item.get("reason"), str)
        or not item.get("reason")
        or not isinstance(item.get("detail"), str)
        or not item.get("detail")
    ):
        _stop("PREF1_INCONCLUSIVE_DRIFT", retry)
    return item


def _validate_raw_width_record(
    value: object, width: Mapping[str, Any]
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _stop("PREF1_WIDTH_SCHEMA_DRIFT", width["retry"])
    base = {
        "retry",
        "width_class",
        "predecessor_generation",
        "attempted_width_hex",
        "replay_receipt",
    }
    width_class = value.get("width_class")
    if width_class in {"all_18_channel_pass", "one_or_more_channel_nonpass"}:
        expected_keys = base | {
            "shadow_path_count",
            "shadow_proposal_count",
            "SSPRK3_stage_and_endpoint_record_count",
            "exact_complete_C",
            "complete_admission_passed",
        }
    elif width_class == "inconclusive":
        expected_keys = base | {
            "shadow_path_count",
            "shadow_proposal_count",
            "SSPRK3_stage_and_endpoint_record_count",
            "typed_inconclusive",
        }
    elif width_class == "premise_stop":
        expected_keys = base | {"typed_stop"}
    else:
        _stop("PREF1_WIDTH_CLASS_DRIFT", (width["retry"], width_class))
    if set(value) != expected_keys:
        _stop("PREF1_WIDTH_SCHEMA_DRIFT", (width["retry"], sorted(value)))
    if (
        value.get("retry") != width["retry"]
        or value.get("predecessor_generation") != width["generation"]
        or value.get("predecessor_generation")
        == width["forbidden_predecessor_generation"]
        or value.get("attempted_width_hex") != width["attempted_width_hex"]
        or value.get("replay_receipt") != _expected_replay_receipt(width)
    ):
        _stop("PREF1_WIDTH_IDENTITY_DRIFT", width["retry"])
    if width_class in {
        "all_18_channel_pass",
        "one_or_more_channel_nonpass",
        "inconclusive",
    } and (
        value.get("shadow_path_count") != 7
        or value.get("shadow_proposal_count") != 7
        or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
    ):
        _stop("PREF1_SHADOW_COUNT_DRIFT", width["retry"])
    if width_class in {"all_18_channel_pass", "one_or_more_channel_nonpass"}:
        exact = _validate_raw_exact(value.get("exact_complete_C"), retry=width["retry"])
        complete = exact.get("complete_admission_passed") is True
        expected_class = (
            "all_18_channel_pass" if complete else "one_or_more_channel_nonpass"
        )
        if (
            width_class != expected_class
            or value.get("complete_admission_passed") is not complete
        ):
            _stop("PREF1_WIDTH_CLASS_DRIFT", width["retry"])
    elif width_class == "inconclusive":
        _validate_typed_inconclusive(
            value.get("typed_inconclusive"), retry=width["retry"]
        )
    else:
        _validate_typed_stop(value.get("typed_stop"), retry=width["retry"])
    return value


def reduce_qa2_terminal(width_records: Sequence[Mapping[str, Any]]) -> str:
    if len(width_records) != 2 or tuple(
        item.get("retry") for item in width_records
    ) != (
        4,
        5,
    ):
        _stop("PREF1_WIDTH_ORDER_DRIFT", [item.get("retry") for item in width_records])
    classes = tuple(item.get("width_class") for item in width_records)
    allowed = {
        "all_18_channel_pass",
        "one_or_more_channel_nonpass",
        "inconclusive",
        "premise_stop",
    }
    if any(item not in allowed for item in classes):
        _stop("PREF1_WIDTH_CLASS_DRIFT", classes)
    if "premise_stop" in classes:
        return PREMISE_STOP_CLASS
    if "inconclusive" in classes:
        return INCONCLUSIVE_CLASS
    if "one_or_more_channel_nonpass" in classes:
        return RAW_CLASSIFICATION
    return BOTH_PASS_CLASS


def _validate_raw_terminal_shape(terminal: object) -> list[Mapping[str, Any]]:
    if not isinstance(terminal, Mapping):
        _stop("PREF1_TERMINAL_SCHEMA_DRIFT", type(terminal).__name__)
    classification = terminal.get("classification")
    if classification not in _RAW_TERMINAL_CLASSES:
        _stop("PREF1_TERMINAL_CLASS_DRIFT", classification)
    robustness = terminal.get("width_robustness_passed") is True
    base = _expected_terminal_base(str(classification), robustness=robustness)
    extra = {"widths"}
    if classification == PREMISE_STOP_CLASS:
        extra.add("typed_stop")
    elif classification == INCONCLUSIVE_CLASS:
        extra.add("typed_inconclusive")
    if set(terminal) != set(base) | extra:
        _stop("PREF1_TERMINAL_SCHEMA_DRIFT", sorted(terminal))
    if any(terminal.get(name) != expected for name, expected in base.items()):
        _stop("PREF1_TERMINAL_NONCLAIM_DRIFT", classification)
    raw_widths = terminal.get("widths")
    if not isinstance(raw_widths, list) or len(raw_widths) != 2:
        _stop("PREF1_WIDTH_SCHEMA_DRIFT", type(raw_widths).__name__)
    widths = [
        _validate_raw_width_record(value, width)
        for value, width in zip(raw_widths, WIDTHS, strict=True)
    ]
    reduced = reduce_qa2_terminal(widths)
    if classification != reduced:
        _stop("PREF1_TERMINAL_CLASS_DRIFT", (classification, reduced))
    if classification == PREMISE_STOP_CLASS:
        top = _mapping(
            terminal.get("typed_stop"),
            frozenset({"type", "detail", "retry"}),
            "terminal.typed_stop",
        )
        if top.get("retry") not in {4, 5}:
            _stop("PREF1_STOP_DRIFT", top)
        _validate_typed_stop(
            {"type": top.get("type"), "detail": top.get("detail")},
            retry=int(top["retry"]),
        )
        owner = next(
            item for item in widths if item.get("width_class") == "premise_stop"
        )
        expected_top = {**dict(owner["typed_stop"]), "retry": owner["retry"]}
        if dict(top) != expected_top:
            _stop("PREF1_STOP_DRIFT", "top/width mismatch")
    if classification == INCONCLUSIVE_CLASS:
        top = _mapping(
            terminal.get("typed_inconclusive"),
            frozenset({"type", "reason", "detail", "retry"}),
            "terminal.typed_inconclusive",
        )
        if top.get("retry") not in {4, 5}:
            _stop("PREF1_INCONCLUSIVE_DRIFT", top)
        _validate_typed_inconclusive(
            {
                "type": top.get("type"),
                "reason": top.get("reason"),
                "detail": top.get("detail"),
            },
            retry=int(top["retry"]),
        )
        owner = next(
            item for item in widths if item.get("width_class") == "inconclusive"
        )
        expected_top = {
            **dict(owner["typed_inconclusive"]),
            "retry": owner["retry"],
        }
        if dict(top) != expected_top:
            _stop("PREF1_INCONCLUSIVE_DRIFT", "top/width mismatch")
    return widths


def _compare_terminal_channel(live: Mapping[str, Any], raw: Mapping[str, Any]) -> None:
    channel = live["channel"]
    raw_upper = _fraction_from_rational(
        raw.get("d12_upper"), label=f"{channel}.raw.d12_upper"
    )
    live_upper = _fraction_from_rational(
        live["D12"]["upper"], label=f"{channel}.live.d12_upper"
    )
    if (
        raw.get("channel") != channel
        or raw.get("admission_passed") is not live["admission_passed"]
        or raw.get("classification") != live["classification"]
        or raw_upper != live_upper
    ):
        _stop("PREF1_TERMINAL_CHANNEL_DRIFT", channel)


def _d12_stream_digest(retry: int, channels: Sequence[Mapping[str, Any]]) -> str:
    answer = sha256(_D12_STREAM_DOMAIN)
    for channel in channels:
        encoded = _mapping(
            channel.get("D12"), frozenset(_DIFFERENCE_KEYS), f"retry-{retry}.D12"
        ).get("upper")
        rational = _mapping(
            encoded, _RATIONAL_KEYS, f"retry-{retry}.{channel.get('channel')}.upper"
        )
        answer.update(
            (
                f"{retry}|{channel.get('channel')}|"
                f"{rational['numerator']}/{rational['denominator']}\n"
            ).encode("ascii")
        )
    return answer.hexdigest()


def _recompute_width(
    prepared: tdg6.TDG6PreparedGR0Compositor,
    raw_width: Mapping[str, Any],
    width: Mapping[str, Any],
    live_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    try:
        assessment = _assess_prepared_independently(prepared, width)
    except Exception as exc:
        if isinstance(exc, TDG10QA2PREF1Error):
            raise
        raise TDG10QA2PREF1Error("PREF1_EXACT_RECOMPUTE_DRIFT", exc) from exc
    raw_exact = _validate_raw_exact(
        raw_width.get("exact_complete_C"), retry=width["retry"]
    )
    raw_channels = raw_exact["channels"]
    if tuple(item.channel for item in assessment.channel_evidence) != CHANNEL_ORDER:
        _stop("PREF1_CHANNEL_ORDER_DRIFT", width["retry"])
    channels: list[dict[str, Any]] = []
    failed: list[str] = []
    for item, raw_channel in zip(
        assessment.channel_evidence, raw_channels, strict=True
    ):
        live = _channel_replay(item)
        _compare_terminal_channel(live, raw_channel)
        if live["admission_passed"] is not True:
            failed.append(item.channel)
        channels.append(live)
    complete = not failed
    expected_width_class = (
        "all_18_channel_pass" if complete else "one_or_more_channel_nonpass"
    )
    stream = _d12_stream_digest(int(width["retry"]), channels)
    if (
        assessment.complete_admission_passed is not complete
        or tuple(assessment.failed_channels) != tuple(failed)
        or raw_exact.get("complete_admission_passed") is not complete
        or raw_exact.get("failed_channels") != failed
        or raw_width.get("complete_admission_passed") is not complete
        or raw_width.get("width_class") != expected_width_class
        or tuple(failed) != tuple(width["failed_channels"])
        or expected_width_class != width["width_class"]
        or stream != width["d12_upper_stream_sha256"]
        or live_receipt != _expected_replay_receipt(width)
    ):
        _stop("PREF1_ALL_OF_DRIFT", (width["retry"], failed, stream))
    return {
        "SSPRK3_stage_and_endpoint_record_count": 28,
        "attempted_width_hex": width["attempted_width_hex"],
        "channel_count": 18,
        "channels": channels,
        "complete_admission_passed": complete,
        "d12_upper_stream_sha256": stream,
        "failed_channels": failed,
        "independent_route_required": True,
        "interval_owner": INTERVAL_OWNER,
        "predecessor_generation": width["generation"],
        "raw_width_class_matches_independent_reduction": True,
        "replay_receipt": dict(live_receipt),
        "retry": width["retry"],
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "sufficient_condition": SUFFICIENT_CONDITION,
        "width_class": expected_width_class,
    }


def _conclusion(classification: str) -> dict[str, Any]:
    if classification not in _RAW_TERMINAL_CLASSES:
        _stop("PREF1_TERMINAL_CLASS_DRIFT", classification)
    both_pass = classification == BOTH_PASS_CLASS
    nonpass = classification == RAW_CLASSIFICATION
    return {
        "PDE_or_candidate_state_opened": False,
        "both_widths_all_eighteen_exact_complete_C_channels_pass": both_pass,
        "exact_remedy_rejected_on_tested_retry_neighborhood": nonpass,
        "independent_method_agreement": False,
        "physics_inference_permitted": False,
        "production_SSPRK3_comparator": False,
        "production_method_earned": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "state_advance_authorized": False,
        "successor_remedy_selected": False,
        "accepted_state": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "external_publication_authorized": False,
        "typed_inconclusive_only": classification == INCONCLUSIVE_CLASS,
        "typed_premise_stop_only": classification == PREMISE_STOP_CLASS,
        "width_robustness_passed": both_pass,
    }


def _static_bound_payload(config_raw: bytes) -> dict[str, Any]:
    expected = expected_evidence(config_raw)
    sealed = _sealed_store()
    return {
        **expected,
        "raw_binding": {
            "canonical_duplicate_free_schema_verified": True,
            "leaf_count": 2,
            "manifest_sha256": RAW_MANIFEST_SHA256,
            "raw_namespace_unchanged": True,
            "terminal_sha256": RAW_TERMINAL_SHA256,
            "two_leaf_tree_verified": True,
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
        "conclusion": _conclusion(RAW_CLASSIFICATION),
    }


def bind_raw_result(config_raw: bytes, repository: Path) -> dict[str, Any]:
    expected = _static_bound_payload(config_raw)
    root = Path(os.path.abspath(os.fspath(repository)))
    metadata = root.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _stop("PREF1_REPOSITORY_UNSAFE", root)
    _environment()
    blobs, manifest, terminal = _raw_snapshot(root)
    if manifest != _expected_manifest():
        _stop("PREF1_RAW_MANIFEST_DRIFT", "manifest schema differs")
    _authenticate_authority(root, manifest)
    _authenticate_live_working_boundary(root)
    _authenticate_live_dependencies(root)
    before_store = _snapshot_store(root)
    sealed_store = (SEALED_STORE_LEAF_COUNT, SEALED_STORE_SHA256)
    if before_store != sealed_store:
        _stop("PREF1_STORE_IDENTITY_DRIFT", before_store)
    raw_widths = _validate_raw_terminal_shape(terminal)
    if terminal.get("classification") != RAW_CLASSIFICATION:
        _stop("PREF1_TERMINAL_CLASS_DRIFT", terminal.get("classification"))
    replay_widths: list[dict[str, Any]] = []
    for width, raw_width in zip(WIDTHS, raw_widths, strict=True):
        if raw_width.get("width_class") not in {
            "all_18_channel_pass",
            "one_or_more_channel_nonpass",
        }:
            _stop("PREF1_TERMINAL_CLASS_DRIFT", raw_width.get("width_class"))
        prepared, receipt = _restore_shadow(root, width)
        replay_widths.append(_recompute_width(prepared, raw_width, width, receipt))
    independent_class = reduce_qa2_terminal(replay_widths)
    if independent_class != terminal.get("classification"):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", independent_class)
    after_store = _snapshot_store(root)
    blobs_after = _raw_tree(root)
    if after_store != before_store or blobs_after != blobs:
        _stop("PREF1_INPUT_MUTATED_DURING_BIND", "raw/store input changed")
    return {
        **expected,
        "independent_replay": {
            "SSPRK3_stage_and_endpoint_record_count": 56,
            "all_of_identity_verified": True,
            "channel_assessment_count": 36,
            "raw_terminal_class_matches_independent_reduction": True,
            "selected_retries": [4, 5],
            "shadow_path_count": 14,
            "shadow_proposal_count": 14,
            "terminal_classification": independent_class,
            "width_robustness_passed": independent_class == BOTH_PASS_CLASS,
            "widths": replay_widths,
        },
    }


def build_pref1_result(
    config_raw: bytes, repository: Path, *, live: bool = True
) -> dict[str, Any]:
    payload = (
        bind_raw_result(config_raw, repository)
        if live
        else _static_bound_payload(config_raw)
    )
    return {"artifact_id": ARTIFACT_ID, "artifact_payload": payload}


_DIFFERENCE_KEYS = {
    "candidate_count",
    "co_maximizer_count",
    "coefficient_stream_sha256",
    "independent_evaluator_id",
    "independent_stationary_count_stream_sha256",
    "localization_classification",
    "lower",
    "maximum_candidates",
    "polynomial_count",
    "primary_evaluator_id",
    "primary_stationary_count_stream_sha256",
    "routes_agree",
    "survivor_key_stream_sha256",
    "upper",
}
_CHANNEL_REPLAY_KEYS = {
    "D01",
    "D12",
    "admission_passed",
    "candidate_evidence_available",
    "channel",
    "classification",
    "combined_coefficient_stream_sha256",
    "evaluator_id",
    "hash_evidence_available",
    "maximum_candidates_D01",
    "maximum_candidates_D12",
    "order_threshold_passed",
    "order_threshold_resolved",
    "refinement_depth",
    "route_evidence_available",
    "row_count",
    "row_stream_sha256",
    "sufficient_condition",
    "sufficient_condition_holds",
    "sufficient_contraction_failure",
    "sufficient_contraction_pass",
    "sufficient_pass_left",
    "sufficient_pass_right",
    "temporal_retry_permitted",
    "threshold_inconclusive",
}


def _validate_difference(
    value: object, *, level: str
) -> tuple[Fraction, Fraction, str]:
    payload = _mapping(value, frozenset(_DIFFERENCE_KEYS), level)
    lower = _fraction_from_rational(payload.get("lower"), label=f"{level}.lower")
    upper = _fraction_from_rational(payload.get("upper"), label=f"{level}.upper")
    if lower < 0 or upper < lower:
        _stop("PREF1_INTERVAL_DRIFT", level)
    d01 = level.endswith(".D01") or level == "D01"
    maximum = (
        PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        if d01
        else PER_CHANNEL_MAXIMUM_CANDIDATES_D12
    )
    polynomial_count = D01_POLYNOMIAL_COUNT if d01 else D12_POLYNOMIAL_COUNT
    candidate_count = _count(
        payload.get("candidate_count"), label=f"{level}.candidate_count"
    )
    co_maximizers = _count(
        payload.get("co_maximizer_count"), label=f"{level}.co_maximizer_count"
    )
    classification = payload.get("localization_classification")
    primary_stationary = _digest(
        payload.get("primary_stationary_count_stream_sha256"),
        label=f"{level}.primary_stationary_count_stream_sha256",
    )
    independent_stationary = _digest(
        payload.get("independent_stationary_count_stream_sha256"),
        label=f"{level}.independent_stationary_count_stream_sha256",
    )
    coefficient_stream = _digest(
        payload.get("coefficient_stream_sha256"),
        label=f"{level}.coefficient_stream_sha256",
    )
    if (
        payload.get("maximum_candidates") != maximum
        or payload.get("polynomial_count") != polynomial_count
        or payload.get("routes_agree") is not True
        or payload.get("primary_evaluator_id") != PRIMARY_EVALUATOR_ID
        or payload.get("independent_evaluator_id") != INDEPENDENT_EVALUATOR_ID
        or classification not in _LOCALIZATION_CLASSES
        or candidate_count > maximum
        or co_maximizers > candidate_count
        or primary_stationary != independent_stationary
        or (classification == "unique_maximum" and co_maximizers != 1)
        or (
            classification == "nonunique_or_interval_inconclusive" and co_maximizers < 2
        )
    ):
        _stop("PREF1_DIFFERENCE_CONTRACT_DRIFT", level)
    for name in ("survivor_key_stream_sha256",):
        _digest(payload.get(name), label=f"{level}.{name}")
    return lower, upper, coefficient_stream


def _validate_rich_channel(value: object, *, expected_channel: str) -> dict[str, Any]:
    payload = _mapping(value, frozenset(_CHANNEL_REPLAY_KEYS), expected_channel)
    if payload.get("channel") != expected_channel:
        _stop(
            "PREF1_CHANNEL_ORDER_DRIFT",
            (expected_channel, payload.get("channel")),
        )
    d01_lower, d01_upper, d01_coefficient_stream = _validate_difference(
        payload.get("D01"), level=f"{expected_channel}.D01"
    )
    d12_lower, d12_upper, d12_coefficient_stream = _validate_difference(
        payload.get("D12"), level=f"{expected_channel}.D12"
    )
    classified = _reclassify_channel(
        channel=expected_channel,
        d01_lower=d01_lower,
        d01_upper=d01_upper,
        d12_lower=d12_lower,
        d12_upper=d12_upper,
    )
    decision = classified["decision"]
    left = _fraction_from_rational(
        payload.get("sufficient_pass_left"), label=f"{expected_channel}.left"
    )
    right = _fraction_from_rational(
        payload.get("sufficient_pass_right"), label=f"{expected_channel}.right"
    )
    combined = sha256(
        _COMBINED_HASH_DOMAIN
        + f"{d01_coefficient_stream}\n{d12_coefficient_stream}\n".encode("ascii")
    ).hexdigest()
    if (
        decision.classification != payload.get("classification")
        or decision.admission_passed is not payload.get("admission_passed")
        or decision.temporal_retry_permitted
        is not payload.get("temporal_retry_permitted")
        or decision.order_threshold_resolved
        is not payload.get("order_threshold_resolved")
        or decision.order_threshold_passed is not payload.get("order_threshold_passed")
        or classified["left"] != left
        or classified["right"] != right
        or classified["sufficient_condition_holds"]
        is not payload.get("sufficient_condition_holds")
        or classified["sufficient_contraction_pass"]
        is not payload.get("sufficient_contraction_pass")
        or classified["sufficient_contraction_failure"]
        is not payload.get("sufficient_contraction_failure")
        or classified["threshold_inconclusive"]
        is not payload.get("threshold_inconclusive")
        or payload.get("combined_coefficient_stream_sha256") != combined
        or payload.get("sufficient_condition") != SUFFICIENT_CONDITION
        or payload.get("row_count") != OWNED_ROW_COUNT
        or payload.get("maximum_candidates_D01") != PER_CHANNEL_MAXIMUM_CANDIDATES_D01
        or payload.get("maximum_candidates_D12") != PER_CHANNEL_MAXIMUM_CANDIDATES_D12
        or payload.get("refinement_depth") != PRIMARY_REFINEMENT_DEPTH
        or payload.get("evaluator_id") != TDG10_EVALUATOR_ID
        or payload.get("route_evidence_available") is not True
        or payload.get("candidate_evidence_available") is not True
        or payload.get("hash_evidence_available") is not True
    ):
        _stop("PREF1_RECLASSIFY_DRIFT", expected_channel)
    _digest(payload.get("row_stream_sha256"), label=f"{expected_channel}.row")
    return dict(payload)


def _validate_compact_width(
    value: object, width: Mapping[str, Any]
) -> Mapping[str, Any]:
    keys = frozenset(
        {
            "SSPRK3_stage_and_endpoint_record_count",
            "attempted_width_hex",
            "channel_count",
            "channels",
            "complete_admission_passed",
            "d12_upper_stream_sha256",
            "failed_channels",
            "independent_route_required",
            "interval_owner",
            "predecessor_generation",
            "raw_width_class_matches_independent_reduction",
            "replay_receipt",
            "retry",
            "shadow_path_count",
            "shadow_proposal_count",
            "sufficient_condition",
            "width_class",
        }
    )
    record = _mapping(value, keys, f"retry-{width['retry']}")
    if (
        record.get("retry") != width["retry"]
        or record.get("predecessor_generation") != width["generation"]
        or record.get("attempted_width_hex") != width["attempted_width_hex"]
        or record.get("replay_receipt") != _expected_replay_receipt(width)
        or record.get("shadow_path_count") != 7
        or record.get("shadow_proposal_count") != 7
        or record.get("SSPRK3_stage_and_endpoint_record_count") != 28
        or record.get("channel_count") != 18
        or record.get("independent_route_required") is not True
        or record.get("interval_owner") != INTERVAL_OWNER
        or record.get("sufficient_condition") != SUFFICIENT_CONDITION
        or record.get("raw_width_class_matches_independent_reduction") is not True
    ):
        _stop("PREF1_WIDTH_IDENTITY_DRIFT", width["retry"])
    channels = record.get("channels")
    if not isinstance(channels, list) or len(channels) != 18:
        _stop("PREF1_CHANNEL_COUNT_DRIFT", width["retry"])
    validated = [
        _validate_rich_channel(item, expected_channel=name)
        for item, name in zip(channels, CHANNEL_ORDER, strict=True)
    ]
    failed = [
        name
        for name, channel in zip(CHANNEL_ORDER, validated, strict=True)
        if channel["admission_passed"] is not True
    ]
    complete = not failed
    expected_class = (
        "all_18_channel_pass" if complete else "one_or_more_channel_nonpass"
    )
    stream = _d12_stream_digest(int(width["retry"]), validated)
    if (
        record.get("complete_admission_passed") is not complete
        or record.get("failed_channels") != failed
        or tuple(failed) != tuple(width["failed_channels"])
        or record.get("width_class") != expected_class
        or expected_class != width["width_class"]
        or record.get("d12_upper_stream_sha256") != stream
        or stream != width["d12_upper_stream_sha256"]
    ):
        _stop("PREF1_ALL_OF_DRIFT", (width["retry"], failed, stream))
    return record


def _validate_independent_replay(value: object) -> None:
    keys = frozenset(
        {
            "SSPRK3_stage_and_endpoint_record_count",
            "all_of_identity_verified",
            "channel_assessment_count",
            "raw_terminal_class_matches_independent_reduction",
            "selected_retries",
            "shadow_path_count",
            "shadow_proposal_count",
            "terminal_classification",
            "width_robustness_passed",
            "widths",
        }
    )
    replay = _mapping(value, keys, "independent_replay")
    widths = replay.get("widths")
    if not isinstance(widths, list) or len(widths) != 2:
        _stop("PREF1_WIDTH_SCHEMA_DRIFT", "compact widths")
    records = [
        _validate_compact_width(item, width)
        for item, width in zip(widths, WIDTHS, strict=True)
    ]
    classification = reduce_qa2_terminal(records)
    if (
        replay.get("SSPRK3_stage_and_endpoint_record_count") != 56
        or replay.get("all_of_identity_verified") is not True
        or replay.get("channel_assessment_count") != 36
        or replay.get("raw_terminal_class_matches_independent_reduction") is not True
        or replay.get("selected_retries") != [4, 5]
        or replay.get("shadow_path_count") != 14
        or replay.get("shadow_proposal_count") != 14
        or replay.get("terminal_classification") != classification
        or replay.get("width_robustness_passed")
        is not (classification == BOTH_PASS_CLASS)
        or classification != RAW_CLASSIFICATION
    ):
        _stop("PREF1_TERMINAL_CLASS_DRIFT", classification)


def validate_compact_result(config_raw: bytes, result_raw: bytes) -> dict[str, Any]:
    expected = _static_bound_payload(config_raw)
    if (
        COMPACT_RESULT_SHA256 is not None
        and sha256(result_raw).hexdigest() != COMPACT_RESULT_SHA256
    ):
        _stop("PREF1_COMPACT_HASH_DRIFT", COMPACT_RESULT_SHA256)
    result = _json(result_raw, "compact result")
    if (
        set(result) != {"artifact_id", "artifact_payload"}
        or result.get("artifact_id") != ARTIFACT_ID
        or not isinstance(result.get("artifact_payload"), Mapping)
    ):
        _stop("PREF1_COMPACT_DRIFT", "compact identity differs")
    payload = result["artifact_payload"]
    if set(payload) != set(expected) | {"independent_replay"}:
        _stop("PREF1_COMPACT_DRIFT", "compact payload schema")
    for key, expected_value in expected.items():
        if payload.get(key) != expected_value:
            _stop("PREF1_COMPACT_DRIFT", key)
    _validate_independent_replay(payload.get("independent_replay"))
    return result


__all__ = (
    "ARTIFACT_ID",
    "AUTHORITY_COMMIT",
    "AUTHORITY_PARENT",
    "BOTH_PASS_CLASS",
    "CLASSIFICATION",
    "COMPACT_RESULT_SHA256",
    "CONFIG_PATH",
    "ENVIRONMENT",
    "INCONCLUSIVE_CLASS",
    "OWNER_DOCUMENT",
    "PREMISE_STOP_CLASS",
    "RAW_CLASSIFICATION",
    "RESULT_PATH",
    "TDG10QA2PREF1Error",
    "WIDTHS",
    "bind_raw_result",
    "build_pref1_result",
    "canonical_result",
    "expected_config",
    "expected_evidence",
    "reduce_qa2_terminal",
    "validate_compact_result",
    "width_spec",
)
