#!/usr/bin/env python3
"""Build or verify the prospective TDG8 replay-evidence freeze."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import tomllib


ROOT = Path(__file__).resolve().parents[1]
CONFIG = Path("configs/fgc/fgc-1-tdg8-frz1.toml")
RESULT = Path("results/fgc-1-tdg8-frz1.json")
ARTIFACT = "FGC-1-TDG8-FRZ1"
PREDECESSOR = "0ba8b0fa5c9d683905ff92f901ab31718f4b1cb7"


class TDG8FreezeError(ValueError):
    """The prospective diagnosis freeze differs from its immutable inputs."""


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _canonical(value: object) -> bytes:
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


def _git_show(repository: Path, commit: str, path: str) -> bytes:
    environment = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("GIT_")
    }
    environment.update({"LC_ALL": "C", "LANG": "C"})
    completed = subprocess.run(
        [
            "git",
            "--no-replace-objects",
            "--no-optional-locks",
            "-c",
            "core.fsmonitor=false",
            "-c",
            "core.untrackedCache=false",
            "-C",
            str(repository),
            "show",
            f"{commit}:{path}",
        ],
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=environment,
    )
    if completed.returncode != 0:
        raise TDG8FreezeError(f"immutable binding is unavailable: {path}")
    return completed.stdout


def build_result(config_raw: bytes, repository: Path = ROOT) -> dict[str, object]:
    try:
        config = tomllib.loads(config_raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise TDG8FreezeError("TDG8 config is malformed") from error
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT
        or config.get("target_protocol") != "FGC-2-SF1-PROTO18"
        or config.get("classification")
        != "prospective_retry_evidence_representation_diagnosis_freeze"
    ):
        raise TDG8FreezeError("TDG8 identity differs")
    predecessor = config.get("predecessor")
    bindings = config.get("bindings")
    question = config.get("question")
    repair = config.get("conditional_repair")
    successor = config.get("successor")
    claims = config.get("claims")
    if not all(
        isinstance(value, dict)
        for value in (predecessor, bindings, question, repair, successor, claims)
    ):
        raise TDG8FreezeError("TDG8 section is absent")
    if (
        predecessor.get("commit") != PREDECESSOR
        or predecessor.get("artifact_id") != "FGC-1-PRO19-PREF28"
        or predecessor.get("terminal_generation") != 9
        or predecessor.get("diagnostic_generation") != 8
        or predecessor.get("affected_member_key") != "RK4-2049"
        or predecessor.get("event") != 23
        or predecessor.get("target_rational") != "3/2"
    ):
        raise TDG8FreezeError("TDG8 predecessor scope differs")
    paths = {
        str(predecessor["result_path"]): str(predecessor["result_sha256"]),
        str(bindings["attempt_path"]): str(bindings["attempt_sha256"]),
        str(bindings["lifecycle_path"]): str(bindings["lifecycle_sha256"]),
        str(bindings["persistence_path"]): str(bindings["persistence_sha256"]),
    }
    for path, expected in paths.items():
        if len(expected) != 64 or _sha(_git_show(repository, PREDECESSOR, path)) != expected:
            raise TDG8FreezeError(f"immutable binding differs: {path}")
    attempt_source = _git_show(repository, PREDECESSOR, str(bindings["attempt_path"]))
    persistence_source = _git_show(
        repository, PREDECESSOR, str(bindings["persistence_path"])
    )
    if b"seen[0] != dict(evidence)" not in attempt_source:
        raise TDG8FreezeError("frozen raw replay comparison is absent")
    if b"normalizing tuples to arrays" not in persistence_source:
        raise TDG8FreezeError("frozen canonical persistence owner is absent")
    if question.get("hypotheses") != [
        "representation_only_tuple_list_drift",
        "semantic_value_or_order_drift",
        "different_replay_classification_or_runtime_stop",
    ]:
        raise TDG8FreezeError("TDG8 outcome partition differs")
    if (
        repair.get("replacement_comparison")
        != "p15._json_safe(seen[0]) != dict(evidence)"
        or repair.get("preserve_scalar_types") is not True
        or repair.get("preserve_array_order") is not True
        or repair.get("preserve_all_TDG6_thresholds") is not True
        or repair.get("forbid_persisted_validator_weakening") is not True
        or repair.get("forbid_locked_campaign_mutation") is not True
    ):
        raise TDG8FreezeError("TDG8 conditional repair differs")
    if (
        successor.get("fresh_campaign_required") is not True
        or successor.get("fresh_namespace")
        != "runs/fgc-2-sf1/proto19/calibration"
        or successor.get("restart_source") != "authenticated_generation_zero"
        or successor.get("branch") != "GR-0"
        or successor.get("amplitude") != "3"
        or successor.get("start_rational") != "23/16"
        or successor.get("target_rational") != "3/2"
        or successor.get("event_count") != 1
        or successor.get("candidate_branches_forbidden") is not True
        or successor.get("execution_requires_separate_committed_image_authority")
        is not True
    ):
        raise TDG8FreezeError("TDG8 successor boundary differs")
    if claims != {
        "root_cause_localized": False,
        "repair_applied": False,
        "successor_campaign_authorized": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "candidate_execution_authorized": False,
        "physical_result_earned": False,
    }:
        raise TDG8FreezeError("TDG8 claim boundary differs")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT,
        "project_version": config["project_version"],
        "target_protocol": config["target_protocol"],
        "classification": config["classification"],
        "gate_status": "pass",
        "source_config_sha256": _sha(config_raw),
        "artifact_payload": {
            "predecessor": dict(predecessor),
            "immutable_bindings": [
                {"path": path, "sha256": digest}
                for path, digest in sorted(paths.items())
            ],
            "question": dict(question),
            "conditional_repair": dict(repair),
            "successor": dict(successor),
            "claims": dict(claims),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    config_raw = (ROOT / CONFIG).read_bytes()
    expected = _canonical(build_result(config_raw))
    if arguments.write:
        destination = ROOT / RESULT
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
        temporary.write_bytes(expected)
        os.replace(temporary, destination)
    elif (ROOT / RESULT).read_bytes() != expected:
        raise SystemExit("TDG8 freeze result differs")
    print(json.dumps({"artifact_id": ARTIFACT, "gate_status": "pass"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
