"""Process-local RSRC2 identity for the constant-memory first-event rerun.

RSRC2 deliberately reuses the already-qualified RSRC1 scientific member,
attempt, protocol, and publisher implementation.  A dedicated parent and each
fresh child activate this profile before constructing any campaign object.  The
profile changes only provenance/schema names, the fresh namespace, the child
ceiling, and the executable paths; it does not change equations, grids,
methods, thresholds, retry caps, or accepted-state semantics.
"""

from __future__ import annotations

from types import MappingProxyType


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC2"
FREEZE_ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC2-FRZ1"
CAMPAIGN_ID = "FGC-2-SF1-PRO20-RSRC2-EVENT1"
PRODUCTION_CONTAINER = "runs/fgc-2-sf1/pro20-rsrc2-event1"
PRODUCTION_NAMESPACE = f"{PRODUCTION_CONTAINER}/event"
PARENT_MAX_CURRENT_RSS_BYTES = 1 * 1024**3
CHILD_MAX_PEAK_RSS_BYTES = 8 * 1024**3
HOST_RESERVE_BYTES = 2 * 1024**3
MIN_AVAILABLE_MEMORY_BYTES = (
    PARENT_MAX_CURRENT_RSS_BYTES
    + CHILD_MAX_PEAK_RSS_BYTES
    + HOST_RESERVE_BYTES
)


def _replace(module, name: str, old, new) -> None:
    current = getattr(module, name)
    if current == new:
        return
    if current != old:
        raise RuntimeError(f"RSRC2 profile found unexpected {module.__name__}.{name}")
    setattr(module, name, new)


def activate_rsrc2() -> None:
    """Activate the explicit RSRC2 profile in one isolated process."""

    from . import pro20_rsrc1_attempt as attempt
    from . import pro20_rsrc1_isolation as isolation
    from . import pro20_rsrc1_member as member
    from . import pro20_rsrc1_protocol as protocol
    from . import pro20_rsrc1_runtime as runtime
    from . import pro20_rsrc1_seed as seed
    from . import pro20_rsrc1_store as store

    _replace(attempt, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)

    _replace(member, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)
    _replace(
        member,
        "CAMPAIGN_ID",
        "FGC-2-SF1-PRO20-RSRC1-EVENT1",
        CAMPAIGN_ID,
    )

    _replace(seed, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)
    _replace(
        seed,
        "CAMPAIGN_ID",
        "FGC-2-SF1-PRO20-RSRC1-EVENT1",
        CAMPAIGN_ID,
    )
    _replace(
        seed,
        "SEED_REQUEST_SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-seed-request-v1",
        "FGC-1-PRO20-EV1-RSRC2-seed-request-v1",
    )
    _replace(
        seed,
        "SEED_RESPONSE_SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-seed-response-v1",
        "FGC-1-PRO20-EV1-RSRC2-seed-response-v1",
    )

    _replace(isolation, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)
    _replace(
        isolation,
        "REQUEST_SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-child-request-v1",
        "FGC-1-PRO20-EV1-RSRC2-child-request-v1",
    )
    _replace(
        isolation,
        "RESPONSE_SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-child-response-v1",
        "FGC-1-PRO20-EV1-RSRC2-child-response-v1",
    )
    _replace(
        isolation,
        "CHILD_MAX_PEAK_RSS_BYTES",
        4 * 1024**3,
        CHILD_MAX_PEAK_RSS_BYTES,
    )
    _replace(
        isolation,
        "PARENT_MAX_CURRENT_RSS_BYTES",
        1 * 1024**3,
        PARENT_MAX_CURRENT_RSS_BYTES,
    )
    _replace(
        isolation,
        "HOST_RESERVE_BYTES",
        4 * 1024**3,
        HOST_RESERVE_BYTES,
    )
    _replace(
        isolation,
        "MIN_AVAILABLE_MEMORY_BYTES",
        9 * 1024**3,
        MIN_AVAILABLE_MEMORY_BYTES,
    )
    isolation.RESOURCE_POLICY = MappingProxyType(
        {
            **dict(isolation.RESOURCE_POLICY),
            "child_max_peak_rss_bytes": CHILD_MAX_PEAK_RSS_BYTES,
            "parent_max_current_rss_bytes": PARENT_MAX_CURRENT_RSS_BYTES,
            "host_reserve_bytes": HOST_RESERVE_BYTES,
            "min_available_memory_bytes": MIN_AVAILABLE_MEMORY_BYTES,
        }
    )

    protocol_values = {
        "FREEZE_ARTIFACT_ID": FREEZE_ARTIFACT_ID,
        "ARTIFACT_ID": ARTIFACT_ID,
        "STORE_KIND_PRODUCTION": "production_rsrc2_event1",
        "PRODUCTION_NAMESPACE": PRODUCTION_NAMESPACE,
        "RECEIPT_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-authority-receipt-v1",
        "ROOT_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-campaign-root-v1",
        "JOURNAL_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-campaign-journal-v1",
        "CHECKPOINT_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-campaign-checkpoint-v1",
        "GENERATION_MANIFEST_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-generation-manifest-v1",
        "SESSION_OPEN_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-writer-session-open-v1",
        "SESSION_CLOSE_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-writer-session-close-v1",
        "TERMINAL_LOCK_SCHEMA": "FGC-1-PRO20-EV1-RSRC2-terminal-lock-v1",
    }
    for name, new in protocol_values.items():
        old = getattr(protocol, name)
        if isinstance(old, str):
            old = old.replace("RSRC2", "RSRC1").replace(
                "production_rsrc2_event1", "production_rsrc1_event1"
            ).replace("pro20-rsrc2-event1", "pro20-rsrc1-event1")
        _replace(protocol, name, old, new)

    _replace(
        store,
        "FREEZE_ARTIFACT_ID",
        "FGC-1-PRO20-EV1-RSRC1-FRZ1",
        FREEZE_ARTIFACT_ID,
    )
    _replace(
        store,
        "PRODUCTION_NAMESPACE",
        "runs/fgc-2-sf1/pro20-rsrc1-event1/event",
        PRODUCTION_NAMESPACE,
    )
    _replace(
        store,
        "WRITER_EXCLUSION_SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-writer-exclusion-v1",
        "FGC-1-PRO20-EV1-RSRC2-writer-exclusion-v1",
    )
    _replace(store, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)
    _replace(
        store,
        "PRODUCTION_CONTAINER",
        "runs/fgc-2-sf1/pro20-rsrc1-event1",
        PRODUCTION_CONTAINER,
    )

    _replace(runtime, "ARTIFACT_ID", "FGC-1-PRO20-EV1-RSRC1", ARTIFACT_ID)
    _replace(
        runtime,
        "FUTURE_FREEZE_ARTIFACT_ID",
        "FGC-1-PRO20-EV1-RSRC1-FRZ1",
        FREEZE_ARTIFACT_ID,
    )
    _replace(
        runtime,
        "SCHEMA",
        "FGC-1-PRO20-EV1-RSRC1-parent-runtime-v1",
        "FGC-1-PRO20-EV1-RSRC2-parent-runtime-v1",
    )
    _replace(
        runtime,
        "AUTHORITY_MODULE",
        "recursive_horizons.fgc.evolution.pro20_rsrc1_authority",
        "recursive_horizons.fgc.evolution.pro20_rsrc2_authority",
    )
    _replace(
        runtime,
        "CHILD_SCRIPT",
        "scripts/run_fgc_pro20_rsrc1_child.py",
        "scripts/run_fgc_pro20_rsrc2_child.py",
    )
    _replace(
        runtime,
        "OWNER_TOKEN",
        "pro20-rsrc1-first-event-one-shot",
        "pro20-rsrc2-first-event-one-shot",
    )
    _replace(
        runtime,
        "CAMPAIGN_ID",
        "FGC-2-SF1-PRO20-RSRC1-EVENT1",
        CAMPAIGN_ID,
    )
    _replace(
        runtime,
        "PRODUCTION_NAMESPACE",
        "runs/fgc-2-sf1/pro20-rsrc1-event1/event",
        PRODUCTION_NAMESPACE,
    )
    _replace(
        runtime,
        "PRODUCTION_CONTAINER",
        "runs/fgc-2-sf1/pro20-rsrc1-event1",
        PRODUCTION_CONTAINER,
    )
    runtime.IMPLEMENTATION_PATHS = (
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_protocol.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_store.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc1_authority.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc2_profile.py",
        "src/recursive_horizons/fgc/evolution/pro20_rsrc2_authority.py",
        "scripts/run_fgc_pro20_rsrc2_child.py",
        "scripts/run_fgc_pro20_rsrc2.py",
    )


__all__ = [
    "ARTIFACT_ID",
    "CAMPAIGN_ID",
    "CHILD_MAX_PEAK_RSS_BYTES",
    "FREEZE_ARTIFACT_ID",
    "HOST_RESERVE_BYTES",
    "MIN_AVAILABLE_MEMORY_BYTES",
    "PARENT_MAX_CURRENT_RSS_BYTES",
    "PRODUCTION_CONTAINER",
    "PRODUCTION_NAMESPACE",
    "activate_rsrc2",
]
