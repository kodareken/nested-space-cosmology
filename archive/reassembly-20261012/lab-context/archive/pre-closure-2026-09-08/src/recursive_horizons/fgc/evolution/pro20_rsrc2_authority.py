"""Exact committed-image authority for one high-capacity RSRC2 event run."""

from __future__ import annotations

from pathlib import Path

from .pro20_rsrc2_profile import (
    ARTIFACT_ID as IMPLEMENTATION_ARTIFACT_ID,
    CAMPAIGN_ID,
    FREEZE_ARTIFACT_ID,
    PRODUCTION_CONTAINER,
    PRODUCTION_NAMESPACE,
    activate_rsrc2,
)


activate_rsrc2()

from . import pro20_rsrc1_authority as _base  # noqa: E402


CONFIG_PATH = "configs/fgc/fgc-1-pro20-ev1-rsrc2-frz1.toml"
OWNER_PATH = "docs/fgc-pro20-ev1-rsrc2-frz1.md"
RUNNER_PATH = "scripts/run_fgc_pro20_rsrc2.py"
AUTHORITY_PATH = "src/recursive_horizons/fgc/evolution/pro20_rsrc2_authority.py"
HARNESS_PATH = "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py"
ALLOWED_FREEZE_DELTA = (CONFIG_PATH, OWNER_PATH)
HISTORICAL_RUNS_FILE_COUNT = 1568
HISTORICAL_RUNS_BYTE_COUNT = 323080865
HISTORICAL_RUNS_DIGEST = (
    "9f1e0e85771fe12dcf2d6cc95af8af2434c95b16e8660bad3952b7ee0211c528"
)


def _configure_base() -> None:
    _base.ARTIFACT_ID = FREEZE_ARTIFACT_ID
    _base.CONFIG_PATH = CONFIG_PATH
    _base.OWNER_PATH = OWNER_PATH
    _base.RUNNER_PATH = RUNNER_PATH
    _base.AUTHORITY_PATH = AUTHORITY_PATH
    _base.HARNESS_PATH = HARNESS_PATH
    _base.OUTPUT_NAMESPACE = PRODUCTION_NAMESPACE
    _base.PRODUCTION_CONTAINER = PRODUCTION_CONTAINER
    _base.ALLOWED_FREEZE_DELTA = ALLOWED_FREEZE_DELTA
    _base.HISTORICAL_RUNS_FILE_COUNT = HISTORICAL_RUNS_FILE_COUNT
    _base.HISTORICAL_RUNS_BYTE_COUNT = HISTORICAL_RUNS_BYTE_COUNT
    _base.HISTORICAL_RUNS_DIGEST = HISTORICAL_RUNS_DIGEST
    _base._SOURCE_DOMAIN = b"FGC-1-PRO20-EV1-RSRC2-SOURCE-CLOSURE-v1\n"
    _base._STAGING_PREFIX = ".pro20-rsrc2-event1.evidence-io-stage-"


_configure_base()


def emit_config_bytes(
    root: Path, *, implementation_commit: str | None = None
) -> bytes:
    _configure_base()
    return _base.emit_config_bytes(root, implementation_commit=implementation_commit)


def check_tracked_config(root: Path) -> dict[str, object]:
    _configure_base()
    return _base.check_tracked_config(root)


def validate_authority_delta(
    *,
    repository_root: Path,
    authority_commit: str,
    implementation_sha256: str,
) -> dict[str, object]:
    _configure_base()
    receipt = dict(
        _base.validate_authority_delta(
            repository_root=repository_root,
            authority_commit=authority_commit,
            implementation_sha256=implementation_sha256,
        )
    )
    receipt["artifact_id"] = FREEZE_ARTIFACT_ID
    receipt["campaign_id"] = CAMPAIGN_ID
    receipt["output_namespace"] = PRODUCTION_NAMESPACE
    receipt["implementation_artifact_id"] = IMPLEMENTATION_ARTIFACT_ID
    receipt["rsrc2_namespace_absent"] = True
    return receipt


recommended_resource_ceilings = _base.recommended_resource_ceilings
environment_identity = _base.environment_identity
implementation_identity = _base.implementation_identity


__all__ = [
    "ALLOWED_FREEZE_DELTA",
    "ARTIFACT_ID",
    "AUTHORITY_PATH",
    "CONFIG_PATH",
    "HARNESS_PATH",
    "OWNER_PATH",
    "RUNNER_PATH",
    "check_tracked_config",
    "emit_config_bytes",
    "environment_identity",
    "implementation_identity",
    "recommended_resource_ceilings",
    "validate_authority_delta",
]


ARTIFACT_ID = FREEZE_ARTIFACT_ID
