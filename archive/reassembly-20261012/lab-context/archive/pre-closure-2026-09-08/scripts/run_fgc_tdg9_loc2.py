#!/usr/bin/env python3
"""Run the bounded read-only TDG9 LOC2 localization diagnostic."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import errno
import json
import os
from pathlib import Path
import secrets
import sys
from typing import Any, Mapping, NoReturn, Sequence


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_loc1 as loc1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_loc2_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema_independent_v2 as independent  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import build_static_gr0_shells  # noqa: E402


RUNNER_ID = "FGC-1-TDG9-LOC2-RUN1"
RAW_SCHEMA = "FGC-1-TDG9-LOC2-raw-v1"


class LOC2RunnerError(RuntimeError):
    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise LOC2RunnerError(owner, code, detail)


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _json_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {"numerator": str(value.numerator), "denominator": str(value.denominator)}
    if is_dataclass(value) and not isinstance(value, type):
        return _json_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    _fail("serialization", "nonexact_value", type(value).__name__)


def _authority(root: Path, commit: str) -> authority.LOC2Authority:
    try:
        return authority.authorize(
            root,
            loc1._read(root, authority.CONFIG_PATH),
            loc1._read(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise LOC2RunnerError("authority", "rejected", exc) from exc


def _prepare_replay(
    root: Path,
    store: HLT16CampaignStore,
    shells: Mapping[str, object],
    replay: Mapping[str, object],
) -> tuple[Any, Mapping[str, object]]:
    try:
        return loc1._prepare(root, store, shells, replay)
    except loc1.LOC1RunnerError as exc:
        raise LOC2RunnerError(
            "replay",
            "sealed_LOC1_reconstruction_failed",
            (exc.owner, exc.code, exc.detail),
        ) from exc


def _cubics(
    rows: Sequence[object],
    level: str,
    component: str,
) -> tuple[list[primary.LocalCubic], list[independent.IndependentLocalCubicV2]]:
    first, _unused_v1 = loc1._cubics(rows, level, component)
    second = [
        independent.IndependentLocalCubicV2(item.coefficients, item.metadata)
        for item in first
    ]
    return first, second


def _candidate_key(item: object) -> tuple[object, ...]:
    return (
        item.polynomial_ordinal,  # type: ignore[attr-defined]
        item.location,  # type: ignore[attr-defined]
        item.location_ordinal,  # type: ignore[attr-defined]
        tuple(item.metadata),  # type: ignore[attr-defined]
    )


def _intervals_overlap(left: object, right: object) -> bool:
    return (
        max(left.parameter_lower, right.parameter_lower)  # type: ignore[attr-defined]
        <= min(left.parameter_upper, right.parameter_upper)  # type: ignore[attr-defined]
        and max(left.absolute_lower, right.absolute_lower)  # type: ignore[attr-defined]
        <= min(left.absolute_upper, right.absolute_upper)  # type: ignore[attr-defined]
    )


def _primary_stationary_count_digest(cubics: Sequence[primary.LocalCubic]) -> str:
    digest = sha256(b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n")
    for ordinal, cubic in enumerate(cubics):
        roots = primary._stationary_intervals(  # type: ignore[attr-defined]
            cubic.coefficients,
            depth=authority.PRIMARY_REFINEMENT_DEPTH,
        )
        digest.update(f"{ordinal}|{len(roots)}\n".encode())
    return digest.hexdigest()


def _require_route_agreement(
    first: object,
    second: object,
    primary_stationary_digest: str,
    detail: object,
) -> None:
    first_items = {_candidate_key(item): item for item in first.candidates}  # type: ignore[attr-defined]
    second_items = {_candidate_key(item): item for item in second.candidates}  # type: ignore[attr-defined]
    if (
        first.polynomial_count != second.polynomial_count  # type: ignore[attr-defined]
        or first.candidate_count != second.candidate_count  # type: ignore[attr-defined]
        or first.classification != second.classification  # type: ignore[attr-defined]
        or primary_stationary_digest != second.stationary_count_stream_sha256  # type: ignore[attr-defined]
        or len(first_items) != len(first.candidates)  # type: ignore[attr-defined]
        or len(second_items) != len(second.candidates)  # type: ignore[attr-defined]
        or tuple(sorted(first_items)) != tuple(sorted(second_items))
        or any(
            not _intervals_overlap(first_items[key], second_items[key])
            for key in first_items
        )
        or max(first.global_absolute_lower, second.global_absolute_lower)  # type: ignore[attr-defined]
        > min(first.global_absolute_upper, second.global_absolute_upper)  # type: ignore[attr-defined]
    ):
        _fail("localization", "independent_disagreement", detail)


def _localize_occurrence(
    rows: Sequence[object],
    retry: int,
    channel: str,
    pref1: Mapping[str, Any],
) -> dict[str, object]:
    row_hash = loc1._row_hash(rows)
    if row_hash != pref1["row_stream_sha256"]:
        _fail("replay", "PREF1_row_hash", (retry, channel))
    sealed = loc1.pref1_exact.assess_exact_temporal_refinement(
        rows,
        max_depth=ar1.MAX_DEPTH,
        max_nodes=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
    )
    sealed_hash = sha256(_canonical(_json_exact(sealed))).hexdigest()
    sealed_independent = loc1.pref1_exact_independent.assess_exact_temporal_refinement_independently(
        rows,
        max_depth=ar1.MAX_DEPTH,
        max_nodes=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR,
    )
    sealed_independent_hash = sha256(_canonical(_json_exact(sealed_independent))).hexdigest()
    if (
        sealed_hash != pref1["primary_evidence_sha256"]
        or sealed_independent_hash != pref1["independent_evidence_sha256"]
        or sealed.classification != "sufficient_contraction_failure"
        or sealed_independent.classification != sealed.classification
    ):
        _fail("replay", "PREF1_evidence_hash", (retry, channel))
    component_evidence: dict[str, dict[str, object]] = {}
    serialized_levels: dict[str, dict[str, object]] = {"D01": {}, "D12": {}}
    for component in ("value_V", "slope_S", "complete_C"):
        component_evidence[component] = {}
        for level, expected in (
            ("D01", 2 * authority.OWNED_ROW_COUNT),
            ("D12", 4 * authority.OWNED_ROW_COUNT),
        ):
            first_cubics, second_cubics = _cubics(rows, level, component)
            if len(first_cubics) != expected or len(second_cubics) != expected:
                _fail("localization", "polynomial_count", (retry, channel, component, level))
            first = primary.localize_absolute_maximum(
                first_cubics,
                maximum_candidates=4 * expected,
                refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
            )
            try:
                second = independent.localize_absolute_maximum_independently_v2(
                    second_cubics,
                    maximum_candidates=4 * expected,
                )
            except independent.RootIsolationInconclusive as exc:
                raise LOC2RunnerError("localization", exc.reason, (retry, channel, component, level, exc.detail)) from exc
            stationary_digest = _primary_stationary_count_digest(first_cubics)
            _require_route_agreement(
                first,
                second,
                stationary_digest,
                (retry, channel, component, level),
            )
            component_evidence[component][level] = first
            serialized_levels[level][component] = {
                "polynomial_count": expected,
                "primary_stationary_count_stream_sha256": stationary_digest,
                "primary": _json_exact(first),
                "independent_v2": _json_exact(second),
                "co_maximizer_count": len(first.candidates),
                "maximizer_input_decomposition": [
                    loc1._hermite_inputs(rows, level, item) for item in first.candidates
                ],
            }
        contraction = loc1._contraction_assessment(
            component_evidence[component]["D01"],
            component_evidence[component]["D12"],
        )
        component_evidence[component]["contraction"] = contraction
        for level in ("D01", "D12"):
            serialized_levels[level][component]["contraction"] = contraction
    if component_evidence["complete_C"]["contraction"]["classification"] != sealed.classification:  # type: ignore[index]
        _fail("localization", "complete_classification_drift", (retry, channel))
    return {
        "retry": retry,
        "channel": channel,
        "PREF1_row_stream_sha256": row_hash,
        "PREF1_primary_evidence_sha256": sealed_hash,
        "PREF1_independent_evidence_sha256": sealed_independent_hash,
        "levels": [
            {
                "level": level,
                "base_cubic_count": 2 * authority.OWNED_ROW_COUNT if level == "D01" else 4 * authority.OWNED_ROW_COUNT,
                "components": serialized_levels[level],
            }
            for level in ("D01", "D12")
        ],
        "component_ownership": loc1._component_ownership(component_evidence),
    }


def _open_output_parent(root: Path) -> tuple[int, str]:
    relative = Path(authority.OUTPUT_NAMESPACE)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in relative.parent.parts:
            try:
                child = os.open(
                    part,
                    os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                os.mkdir(part, mode=0o755, dir_fd=descriptor)
                child = os.open(
                    part,
                    os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            os.close(descriptor)
            descriptor = child
        return descriptor, relative.name
    except Exception:
        os.close(descriptor)
        raise


def _publish(root: Path, manifest: Mapping[str, object], terminal: Mapping[str, object]) -> None:
    parent_fd, target_name = _open_output_parent(root)
    stage_name = f"{authority.STAGING_PREFIX}{os.getpid()}-{secrets.token_hex(8)}"
    stage_fd: int | None = None
    try:
        try:
            os.stat(target_name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            _fail("output", "namespace_exists", authority.OUTPUT_NAMESPACE)
        os.mkdir(stage_name, mode=0o700, dir_fd=parent_fd)
        stage_fd = os.open(
            stage_name,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        for name, value in (("manifest.json", manifest), ("terminal.json", terminal)):
            leaf = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o644,
                dir_fd=stage_fd,
            )
            try:
                pending = memoryview(_pretty(value))
                while pending:
                    written = os.write(leaf, pending)
                    if written <= 0:
                        _fail("output", "short_write", name)
                    pending = pending[written:]
                os.fsync(leaf)
            finally:
                os.close(leaf)
        os.fsync(stage_fd)
        if os.uname().sysname != "Darwin":
            _fail("output", "exclusive_adoption_unavailable", os.uname().sysname)
        rename = getattr(ctypes.CDLL(None, use_errno=True), "renameatx_np", None)
        if rename is None:
            _fail("output", "exclusive_adoption_unavailable", "renameatx_np")
        rename.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
        rename.restype = ctypes.c_int
        if rename(parent_fd, os.fsencode(stage_name), parent_fd, os.fsencode(target_name), 0x00000004):
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        if stage_fd is not None:
            for name in ("manifest.json", "terminal.json"):
                try:
                    os.unlink(name, dir_fd=stage_fd)
                except FileNotFoundError:
                    pass
            os.close(stage_fd)
        try:
            os.rmdir(stage_name, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        raise
    finally:
        os.close(parent_fd)


def prospective_coefficient_depth_audit(root: Path) -> dict[str, object]:
    repository = root.resolve()
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    prepared: dict[int, Any] = {}
    for replay in ar1.REPLAYS:
        value, _historical = _prepare_replay(repository, store, shells, replay)
        prepared[int(replay["retry"])] = value
    maxima = {
        "coefficient_numerator_bits": 0,
        "coefficient_denominator_bits": 0,
        "derivative_numerator_bits": 0,
        "derivative_denominator_bits": 0,
        "discriminant_numerator_bits": 0,
        "discriminant_denominator_bits": 0,
        "primitive_height_bits": 0,
        "proof_bits": 0,
    }
    counts = {
        "base_cubics": 0,
        "component_cubics": 0,
        "endpoint_deflated": 0,
        "positive_rational_square": 0,
        "positive_nonsquare": 0,
        "negative_discriminant": 0,
        "zero_discriminant": 0,
    }
    proof_owner: tuple[object, ...] | None = None
    for retry, channel in authority.FAILED_OCCURRENCES:
        rows = tuple(loc1._rows(loc1._surface(prepared[retry]), channel))
        for level in ("D01", "D12"):
            base, _ = loc1._cubics(rows, level, "complete_C")
            counts["base_cubics"] += len(base)
            for component in ("value_V", "slope_S", "complete_C"):
                cubics, _ = loc1._cubics(rows, level, component)
                counts["component_cubics"] += len(cubics)
                for ordinal, cubic in enumerate(cubics):
                    for coefficient in cubic.coefficients:
                        maxima["coefficient_numerator_bits"] = max(maxima["coefficient_numerator_bits"], abs(coefficient.numerator).bit_length())
                        maxima["coefficient_denominator_bits"] = max(maxima["coefficient_denominator_bits"], coefficient.denominator.bit_length())
                    _, a1, a2, a3 = cubic.coefficients
                    derivative = a1, 2 * a2, 3 * a3
                    for coefficient in derivative:
                        maxima["derivative_numerator_bits"] = max(maxima["derivative_numerator_bits"], abs(coefficient.numerator).bit_length())
                        maxima["derivative_denominator_bits"] = max(maxima["derivative_denominator_bits"], coefficient.denominator.bit_length())
                    b0, b1, b2 = derivative
                    A, B, C = independent._primitive_integer_quadratic(b0, b1, b2)
                    height = max(abs(A), abs(B), abs(C), 1).bit_length()
                    proof = 2 * height + 8
                    if proof > maxima["proof_bits"]:
                        maxima["primitive_height_bits"] = height
                        maxima["proof_bits"] = proof
                        proof_owner = (retry, channel, level, component, ordinal)
                    if C == 0 or A + B + C == 0:
                        counts["endpoint_deflated"] += 1
                    discriminant = b1 * b1 - 4 * b2 * b0
                    maxima["discriminant_numerator_bits"] = max(maxima["discriminant_numerator_bits"], abs(discriminant.numerator).bit_length())
                    maxima["discriminant_denominator_bits"] = max(maxima["discriminant_denominator_bits"], discriminant.denominator.bit_length())
                    if discriminant < 0:
                        counts["negative_discriminant"] += 1
                    elif discriminant == 0:
                        counts["zero_discriminant"] += 1
                    elif independent._perfect_rational_square(discriminant) is None:
                        counts["positive_nonsquare"] += 1
                    else:
                        counts["positive_rational_square"] += 1
    return {
        "counts": counts,
        "maxima": maxima,
        "maximum_proof_owner": list(proof_owner or ()),
        "initial_refinement_bits": independent.INITIAL_REFINEMENT_BITS,
        "refinement_schedule": independent.REFINEMENT_SCHEDULE,
        "global_proof_bit_ceiling": independent.GLOBAL_PROOF_BIT_CEILING,
        "all_proof_caps_within_global_ceiling": maxima["proof_bits"] <= independent.GLOBAL_PROOF_BIT_CEILING,
        "localization_executed": False,
    }


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _authority(repository, authority_commit)
    failures = loc1._pref1_failure_map(repository)
    before = loc1._snapshot_store(repository)
    sealed_store = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    if before != sealed_store:
        _fail("provenance", "sealed_campaign_store_identity", (before, sealed_store))
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    prepared: dict[int, Any] = {}
    historical_hashes: dict[int, str] = {}
    for replay in ar1.REPLAYS:
        value, historical = _prepare_replay(repository, store, shells, replay)
        prepared[int(replay["retry"])] = value
        historical_hashes[int(replay["retry"])] = sha256(_canonical(historical)).hexdigest()
    occurrences: list[dict[str, object]] = []
    base_cubic_count = 0
    candidate_counts = {"value_V": 0, "slope_S": 0, "complete_C": 0}
    for retry, channel in authority.FAILED_OCCURRENCES:
        rows = tuple(loc1._rows(loc1._surface(prepared[retry]), channel))
        result = _localize_occurrence(rows, retry, channel, failures[(retry, channel)])
        occurrences.append(result)
        for level in result["levels"]:  # type: ignore[assignment]
            base_cubic_count += int(level["base_cubic_count"])  # type: ignore[index]
            for component in candidate_counts:
                candidate_counts[component] += int(level["components"][component]["primary"]["candidate_count"])  # type: ignore[index]
    total_candidates = sum(candidate_counts.values())
    if (
        base_cubic_count != receipt.base_cubic_count
        or any(value > authority.FAMILY_CANDIDATE_CEILING for value in candidate_counts.values())
        or total_candidates > authority.AGGREGATE_CANDIDATE_CEILING
    ):
        _fail("localization", "global_budget", (base_cubic_count, candidate_counts))
    after = loc1._snapshot_store(repository)
    if before != after:
        _fail("provenance", "campaign_store_mutated", (before, after))
    manifest = {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "ERR1_result_sha256": authority.ERR1_RESULT_SHA256,
        "base_cubic_count": authority.BASE_CUBIC_COUNT,
        "component_cubic_count": authority.COMPONENT_CUBIC_COUNT,
        "output_leaves": ["manifest.json", "terminal.json"],
    }
    terminal = {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": "bounded_exact_localization_completed",
        "authority_commit": authority_commit,
        "historical_TDG6_evidence_sha256": historical_hashes,
        "base_cubic_count": base_cubic_count,
        "component_candidate_counts": candidate_counts,
        "total_VSC_candidate_count": total_candidates,
        "occurrences": occurrences,
        "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]},
        "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]},
        "store_unchanged": True,
        "evaluator_inconclusive": False,
        "ULP_used_as_tolerance": False,
        "fourth_width_executed": False,
        "stage2_source_decomposition_executed": False,
        "SSPRK3_comparator_executed": False,
        "PDE_state_committed": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }
    _publish(repository, manifest, terminal)
    return terminal


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    arguments = parser.parse_args(argv)
    result = run(ROOT, authority_commit=arguments.authority_commit)
    print(json.dumps({"artifact_id": authority.ARTIFACT_ID, "classification": result["classification"], "base_cubic_count": result["base_cubic_count"], "candidate_branch_opened": False}, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
