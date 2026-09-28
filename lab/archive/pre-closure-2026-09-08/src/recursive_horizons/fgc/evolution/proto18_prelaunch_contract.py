"""Static, fail-closed contract helpers for the PRO18 production prelaunch.

This module deliberately does not read campaign archives, event logs, or
candidate state.  It validates the compact freeze and checks only whether the
two future PROTO17 output roots already exist.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import subprocess
from typing import Any, Mapping


class Proto18PrelaunchStop(ValueError):
    """A static PRO18 premise is not true; no launch may proceed."""


PRO18_TRUE_CLAIMS = frozenset({
    "PROTO17_target_identity_and_hash_bindings_frozen",
    "HLT14_synthetic_result_preserved",
    "PRO18_static_production_prelaunch_overlay_frozen",
    "PRO18_two_source_container_and_fresh_output_boundary_frozen",
    "PREF26_implementation_or_design_authorized",
})

PRO18_FALSE_CLAIMS = frozenset({
    "PREF26_completed", "AUTH1_authorized", "HLT15_GEN1_authorized",
    "PREF27_authorized", "raw_archive_opened", "historical_event_log_opened",
    "future_namespace_created", "pretrajectory_operation_authorized",
    "trajectory_read", "candidate_execution_authorized",
    "physical_claim_authorized", "fresh_GR0_calibration_authorized",
    "fresh_GR0_dynamic_calibration_completed", "GR0_case_eligible",
    "classical_spherical_diagnostic_authorized", "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized", "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized", "physical_transition_claim_authorized",
    "global_continuation_authorized", "FGCQR_mechanism_rejected",
    "general_gradient_route_rejected", "singularity_resolution_derived",
    "child_domain_or_topology_derived", "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
})


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def canonical_sha256(value: object) -> str:
    return sha256(canonical_bytes(value)).hexdigest()


def duplicate_safe_json(payload: bytes, label: str) -> dict[str, Any]:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in pairs:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer
    try:
        value = json.loads(payload.decode("utf-8"), object_pairs_hook=reject_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto18PrelaunchStop(f"{label} is malformed or has duplicate keys") from error
    if not isinstance(value, dict) or canonical_bytes(value) != payload:
        raise Proto18PrelaunchStop(f"{label} is not canonical JSON")
    return value


def _git_show(root: Path, commit: str, relative: str) -> bytes:
    try:
        return subprocess.check_output(
            ["git", "show", f"{commit}:{relative}"], cwd=root, stderr=subprocess.DEVNULL
        )
    except subprocess.CalledProcessError as error:
        raise Proto18PrelaunchStop(f"immutable blob is unavailable: {relative}") from error


def safe_compact_file(root: Path, relative: str) -> Path:
    """Resolve one declared tracked file without following any symlink."""
    candidate_relative = Path(relative)
    if candidate_relative.is_absolute() or ".." in candidate_relative.parts:
        raise Proto18PrelaunchStop(f"unsafe compact-file path: {relative}")
    root_mode = root.lstat().st_mode
    if stat.S_ISLNK(root_mode) or not stat.S_ISDIR(root_mode):
        raise Proto18PrelaunchStop("repository root is a symlink or non-directory")
    candidate = root
    for part in candidate_relative.parts[:-1]:
        candidate /= part
        if not os.path.lexists(candidate):
            raise Proto18PrelaunchStop(f"live compact-file ancestor is absent: {relative}")
        mode = candidate.lstat().st_mode
        if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
            raise Proto18PrelaunchStop(
                f"live compact-file ancestor is a symlink or non-directory: {relative}"
            )
    candidate /= candidate_relative.parts[-1]
    if not os.path.lexists(candidate):
        raise Proto18PrelaunchStop(f"live compact file is absent: {relative}")
    mode = candidate.lstat().st_mode
    if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
        raise Proto18PrelaunchStop(
            f"live compact file is a symlink or non-regular file: {relative}"
        )
    return candidate


def verify_immutable_compact_files(root: Path, commit: str, files: Mapping[str, str]) -> dict[str, str]:
    """Bind every compact source to both immutable Git bytes and live bytes."""
    if len(commit) != 40 or commit.lower() != commit or any(c not in "0123456789abcdef" for c in commit):
        raise Proto18PrelaunchStop("immutable authority commit is not full lowercase SHA-1")
    try:
        subprocess.run(["git", "merge-base", "--is-ancestor", commit, "HEAD"], cwd=root, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as error:
        raise Proto18PrelaunchStop("immutable authority commit is not an ancestor of HEAD") from error
    observed: dict[str, str] = {}
    for relative, expected in files.items():
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise Proto18PrelaunchStop("compact-file binding is malformed")
        historic = _git_show(root, commit, relative)
        live_path = safe_compact_file(root, relative)
        live = live_path.read_bytes()
        digest = sha256(live).hexdigest()
        if historic != live or sha256(historic).hexdigest() != expected or digest != expected:
            raise Proto18PrelaunchStop(f"immutable/live compact binding differs: {relative}")
        observed[relative] = digest
    return observed


def inspect_future_output_roots_only(root: Path, roots: Mapping[str, str]) -> dict[str, bool]:
    """Observe just output-root existence; never traverse or open source inputs."""
    answer: dict[str, bool] = {}
    for name, relative in roots.items():
        if name in {"temporal_absence_observation_at_freeze", "only_output_root_existence_may_be_inspected"}:
            continue
        if not isinstance(relative, str) or not relative.startswith("runs/fgc-2-sf1/proto17/"):
            raise Proto18PrelaunchStop(f"unsafe future output root: {name}")
        path = root / relative
        # The output path itself may be absent, but every existing ancestor
        # between the repository and that path must be a real directory.  A
        # symlinked or non-directory ancestor could redirect the later atomic
        # genesis outside the frozen namespace without making the leaf exist.
        ancestor = root
        for part in Path(relative).parts[:-1]:
            ancestor /= part
            if not os.path.lexists(ancestor):
                continue
            try:
                ancestor_mode = ancestor.lstat().st_mode
            except OSError as error:
                raise Proto18PrelaunchStop("future output ancestor lstat failed") from error
            if stat.S_ISLNK(ancestor_mode) or not stat.S_ISDIR(ancestor_mode):
                raise Proto18PrelaunchStop("future output ancestor is a symlink or non-directory")
        # lexists detects a dangling symlink, while lstat lets us reject every
        # symlink or non-directory rather than following it.
        if os.path.lexists(path):
            try:
                mode = path.lstat().st_mode
            except OSError as error:
                raise Proto18PrelaunchStop("future output root lstat failed") from error
            if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
                raise Proto18PrelaunchStop("future output root is a symlink or non-directory")
            answer[name] = True
        else:
            answer[name] = False
    return answer


def require_exact_claims(claims: Mapping[str, Any]) -> None:
    required = PRO18_TRUE_CLAIMS | PRO18_FALSE_CLAIMS
    if set(claims) != required:
        raise Proto18PrelaunchStop("PRO18 claim schema differs")
    if any(claims[name] is not True for name in PRO18_TRUE_CLAIMS):
        raise Proto18PrelaunchStop("PRO18 earned static claim is absent")
    if any(claims[name] is not False for name in PRO18_FALSE_CLAIMS):
        raise Proto18PrelaunchStop(
            "PRO18 must preserve every production, trajectory, and physical nonclaim"
        )
