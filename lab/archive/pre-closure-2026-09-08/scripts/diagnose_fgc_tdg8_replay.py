#!/usr/bin/env python3
"""Run the frozen read-only TDG8 generation-eight replay discriminator."""
from __future__ import annotations

import argparse
import copy
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import hlt16_campaign_schema as schema  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_lifecycle as lifecycle  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_progression_attempt as attempt  # noqa: E402
from recursive_horizons.fgc.evolution import proto15_runtime as p15  # noqa: E402
from recursive_horizons.fgc.evolution import proto19_pref28_binder as pref28  # noqa: E402
from recursive_horizons.fgc.evolution import proto19_resume_authority as sid3  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg7_stage_safe_runtime as tdg7  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import build_static_gr0_shells  # noqa: E402


ARTIFACT = "FGC-1-TDG8-DIAG1"
FREEZE_COMMIT = "1496542cb1489e158d9edb2925fd54cad34d1f8d"
FREEZE_RESULT = Path("results/fgc-1-tdg8-frz1.json")
PREF28_CONFIG = Path("configs/fgc/fgc-1-pro19-pref28.toml")
RESULT = Path("results/fgc-1-tdg8-diag1.json")
CHECKPOINT_RELATIVE = Path(pref28.STORE_PATH) / (
    "checkpoints/00000000000000000008-"
    f"{sid3.GENERATION8_CHECKPOINT_SHA256}.json"
)


class TDG8DiagnosisError(RuntimeError):
    """The frozen discriminator could not classify the replay evidence."""


def _canonical(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def _pretty(value: object) -> bytes:
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


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _mapping_differences(
    live: object,
    durable: object,
    path: str = "$",
) -> tuple[list[dict[str, str]], list[str]]:
    containers: list[dict[str, str]] = []
    semantics: list[str] = []
    if isinstance(live, (tuple, list)) and isinstance(durable, (tuple, list)):
        if type(live) is not type(durable):
            containers.append(
                {
                    "path": path,
                    "live_type": type(live).__name__,
                    "durable_type": type(durable).__name__,
                }
            )
        if len(live) != len(durable):
            semantics.append(f"{path}.length")
            return containers, semantics
        for index, (left, right) in enumerate(zip(live, durable)):
            child_containers, child_semantics = _mapping_differences(
                left, right, f"{path}[{index}]"
            )
            containers.extend(child_containers)
            semantics.extend(child_semantics)
        return containers, semantics
    if isinstance(live, Mapping) and isinstance(durable, Mapping):
        if set(live) != set(durable):
            semantics.append(f"{path}.keys")
            return containers, semantics
        for key in sorted(live):
            child_containers, child_semantics = _mapping_differences(
                live[key], durable[key], f"{path}.{key}"
            )
            containers.extend(child_containers)
            semantics.extend(child_semantics)
        return containers, semantics
    if type(live) is not type(durable) or live != durable:
        semantics.append(path)
    return containers, semantics


def _load_generation_eight() -> schema.CampaignCheckpoint:
    raw = pref28._read_repository_leaf(ROOT, CHECKPOINT_RELATIVE.as_posix(), "TDG8 generation eight")
    value = pref28._json(raw, "TDG8 generation eight")
    checkpoint = schema.CampaignCheckpoint.from_mapping(value)
    if (
        checkpoint.generation != 8
        or checkpoint.sha256 != sid3.GENERATION8_CHECKPOINT_SHA256
        or checkpoint.journal_sequence != sid3.GENERATION8_JOURNAL_SEQUENCE
        or checkpoint.journal_tip_sha256 != sid3.GENERATION8_JOURNAL_SHA256
        or checkpoint.disposition != "nonterminal"
    ):
        raise TDG8DiagnosisError("generation-eight checkpoint differs")
    return checkpoint


def diagnose() -> dict[str, object]:
    freeze_raw = (ROOT / FREEZE_RESULT).read_bytes()
    freeze = json.loads(freeze_raw.decode("ascii"))
    if (
        freeze.get("artifact_id") != "FGC-1-TDG8-FRZ1"
        or freeze.get("gate_status") != "pass"
    ):
        raise TDG8DiagnosisError("TDG8 freeze is absent")
    pref28_config_raw = (ROOT / PREF28_CONFIG).read_bytes()
    before_store = pref28.bind_terminal_store(pref28_config_raw, ROOT)
    checkpoint = _load_generation_eight()
    key = sid3.GENERATION8_MEMBER_KEY
    state = checkpoint.members[key]
    cursor = schema._cursor(state.cursor)
    ledger = schema._ledger(state.ledger)
    templates = build_static_gr0_shells(ROOT)
    member = runtime.restore_member_with_overlay(
        HLT16CampaignStore(ROOT / pref28.STORE_PATH),
        checkpoint,
        copy.deepcopy(templates[key]),
        key=key,
    )
    boundary = member.snapshot()
    physical_sha256 = array_content_sha256(
        member.state.u, member.state.p, member.state.q
    )
    lifecycle.validate_retry_payload(
        cursor, ledger, evolution_state_sha256=physical_sha256
    )
    pending = cursor.payload["retry_successor_payload_or_none"]
    if not isinstance(pending, Mapping):
        raise TDG8DiagnosisError("generation-eight cursor is not retry-pending")
    durable = pending["TDG6_rejection_evidence"]
    if not isinstance(durable, Mapping):
        raise TDG8DiagnosisError("durable evidence is absent")
    prior = p15._prior_retry_ledger(ledger, durable)
    predecessor = p15._plan_mapping(pending["predecessor_plan"])
    prepared = attempt._prepare_initial(
        member,
        3 / 2,
        prior,
        float.fromhex(predecessor["requested_cap_hex"]),
    )
    seen: list[Mapping[str, object]] = []
    replay_outcome: tdg6.TDG6TemporalRetryRequired | None = None
    try:
        tdg7.require_tdg7_stage_safe_admission(
            prepared,
            transaction=member.transaction,
            tracers=member.tracers,
            temporal_ledger=prior,
            current_time=member.time,
            current_state=member.state,
            current_step_index=member.step_index,
            current_transaction_serial=member.transaction_serial,
            durable_rejection_sink=lambda item: seen.append(dict(item)),
        )
    except tdg6.TDG6TemporalRetryRequired as replayed:
        replay_outcome = replayed
    except tdg6.TDG6TemporalRetryExhausted as error:
        raise TDG8DiagnosisError("replay changed to temporal exhaustion") from error
    else:
        raise TDG8DiagnosisError("replay changed to admission")
    finally:
        member.restore(boundary)
    if len(seen) != 1:
        raise TDG8DiagnosisError("replay emitted a non-singleton evidence set")
    live = seen[0]
    containers, semantics = _mapping_differences(live, durable)
    normalized_live = p15._json_safe(live)
    canonical_equal = _canonical(normalized_live) == _canonical(dict(durable))
    raw_equal = live == dict(durable)
    expected_containers = [
        {
            "path": "$.channel_admissions",
            "live_type": "tuple",
            "durable_type": "list",
        },
        {
            "path": "$.failed_channels",
            "live_type": "tuple",
            "durable_type": "list",
        },
    ]
    if replay_outcome is None:
        raise TDG8DiagnosisError("replay outcome was not retained")
    successor = tdg7.replan_tdg7_after_tdg6_retry(prepared, replay_outcome)
    successor_equal = (
        p15._plan_mapping(successor.plan)
        == p15._plan_mapping(pending["successor_plan"])
    )
    ledger_equal = replay_outcome.updated_ledger == ledger
    after_store = pref28.bind_terminal_store(pref28_config_raw, ROOT)
    classification = (
        "representation_only_tuple_list_drift"
        if (
            raw_equal is False
            and canonical_equal
            and containers == expected_containers
            and not semantics
            and successor_equal
            and ledger_equal
            and before_store == after_store
        )
        else "semantic_or_replay_drift"
    )
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT,
        "project_version": "0.11.0",
        "target_protocol": "FGC-2-SF1-PROTO18",
        "classification": classification,
        "gate_status": "pass" if classification == "representation_only_tuple_list_drift" else "stop",
        "artifact_payload": {
            "freeze_commit": FREEZE_COMMIT,
            "freeze_result_sha256": _sha(freeze_raw),
            "source_checkpoint": {
                "generation": checkpoint.generation,
                "checkpoint_sha256": checkpoint.sha256,
                "journal_sequence": checkpoint.journal_sequence,
                "journal_tip_sha256": checkpoint.journal_tip_sha256,
                "member_key": key,
                "cursor_sha256": cursor.sha256,
                "physical_state_sha256": physical_sha256,
            },
            "comparison": {
                "raw_mapping_equal": raw_equal,
                "canonical_JSON_equal": canonical_equal,
                "container_differences": containers,
                "semantic_difference_paths": semantics,
                "replayed_ledger_equal": ledger_equal,
                "replayed_successor_plan_equal": successor_equal,
                "replay_classification": type(replay_outcome).__name__,
            },
            "store": {
                "locked_terminal_boundary_equal_before_after": before_store == after_store,
                "locked_campaign_mutated": False,
            },
            "decision": {
                "conditional_representation_repair_permitted": classification
                == "representation_only_tuple_list_drift",
                "permitted_expression": "p15._json_safe(seen[0]) != dict(evidence)",
                "TDG6_or_TDG7_rule_change_permitted": False,
                "successor_campaign_authorized": False,
            },
            "claims": {
                "root_cause_mechanism_reproduced": classification
                == "representation_only_tuple_list_drift",
                "PREF28_terminal_retroactively_reclassified": False,
                "repair_applied": False,
                "common_event_completed": False,
                "GR0_calibration_completed": False,
                "candidate_execution_authorized": False,
                "physical_result_earned": False,
            },
        },
    }
    return result


def validate_compact(raw: bytes) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise TDG8DiagnosisError("TDG8 diagnosis result is malformed") from error
    if not isinstance(value, dict) or _pretty(value) != raw:
        raise TDG8DiagnosisError("TDG8 diagnosis result is noncanonical")
    payload = value.get("artifact_payload")
    comparison = payload.get("comparison") if isinstance(payload, dict) else None
    decision = payload.get("decision") if isinstance(payload, dict) else None
    claims = payload.get("claims") if isinstance(payload, dict) else None
    if (
        value.get("artifact_id") != ARTIFACT
        or value.get("classification") != "representation_only_tuple_list_drift"
        or value.get("gate_status") != "pass"
        or not isinstance(comparison, dict)
        or comparison.get("canonical_JSON_equal") is not True
        or comparison.get("semantic_difference_paths") != []
        or not isinstance(decision, dict)
        or decision.get("conditional_representation_repair_permitted") is not True
        or decision.get("successor_campaign_authorized") is not False
        or not isinstance(claims, dict)
        or claims.get("repair_applied") is not False
        or claims.get("physical_result_earned") is not False
    ):
        raise TDG8DiagnosisError("TDG8 diagnosis claim boundary differs")
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    arguments = parser.parse_args()
    if arguments.run:
        result = diagnose()
        destination = ROOT / RESULT
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
        temporary.write_bytes(_pretty(result))
        os.replace(temporary, destination)
    else:
        result = validate_compact((ROOT / RESULT).read_bytes())
    print(
        json.dumps(
            {
                "artifact_id": result["artifact_id"],
                "classification": result["classification"],
                "gate_status": result["gate_status"],
            }
        )
    )
    return 0 if result["gate_status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
