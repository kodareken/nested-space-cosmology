#!/usr/bin/env python3
"""Construct once or compactly verify the premise-only LOC1 ERR1 diagnosis."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
ARTIFACT_ID = "FGC-1-TDG9-LOC1-ERR1"
AUTHORITY_COMMIT = "14e26c43edbc88003d68cbcf24681bb090d98534"
RESULT_PATH = "results/fgc-1-tdg9-loc1-err1.json"
STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
STORE_LEAF_COUNT = 115
STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
PREF2_RESULT_PATH = "results/fgc-1-tdg8-rcv3-pref2.json"
PREF2_RESULT_SHA256 = "7e676cb40fd83deeff0912603641fd93ec139f274287b1d092aa8eccd2c1231f"
PREF2_STORE_MANIFEST_SHA256 = "46b57bc38f3bbfbec70c4b6b11a089cb799dc73135af0d9dac0e3dc289b952d4"
LOC1_NAMESPACE = "runs/fgc-2-sf1/tdg9-loc1/rcv3-retries-3-5"
LOC1_STAGING_PREFIX = ".rcv3-retries-3-5.stage-"
SEALED = {
    "configs/fgc/fgc-1-tdg9-loc1-frz1.toml": "0d1ddd4825fba479f8b0b245790ee4cc712d5ba0d0878fde0e9cc7c914e3a099",
    "results/fgc-1-tdg9-loc1-frz1.json": "1e43d81eb077563e824b661b4a7049cdf13f96d14a710695f7ae74576d9ebd1b",
    "scripts/run_fgc_tdg9_loc1.py": "19d12813b7ee3fde6b044b224124d62909170fdfaffb728414a378984bac0996",
    "src/recursive_horizons/fgc/evolution/tdg9_local_extrema_independent.py": "d02fc7e719f962eb038e2704b004cff2b5ea11fa5f9b1584e82b3869a06f90c7",
    "docs/fgc-tdg9-loc1-frz1.md": "a43e262f639899c1cd4788a191ee84ba18a5bf3e4004e1439b2361569235f2e9",
}


class ERR1Error(RuntimeError):
    pass


def _canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            raise ERR1Error("duplicate compact JSON key")
        answer[key] = value
    return answer


def _expected() -> dict[str, object]:
    maximum = {
        "numerator": "7421011090394622581",
        "denominator": "10384593717069655257060992658440192",
    }
    factor = {
        "numerator": "21210020543942968745638655",
        "denominator": "338460656020607282663380637712778772392143197677711984273740183180495765112991409062496875745134225841966700556811959451779072",
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": "independent_root_enclosure_false_candidate_diagnosis",
        "gate_status": "pass",
        "authority_commit": AUTHORITY_COMMIT,
        "sealed_LOC1_sha256": dict(SEALED),
        "exception": {
            "owner": "localization",
            "code": "independent_disagreement",
            "detail": [3, "u:alpha", "value_V", "D01"],
            "observation_kind": "operator_observed_unserialized",
            "raw_exception_terminal_present": False,
        },
        "diagnosis": {
            "primary_candidate_count": 8176,
            "independent_candidate_count": 8544,
            "affected_cubic_count": 184,
            "false_extra_candidate_count": 368,
            "primary_classification": "nonunique_or_interval_inconclusive",
            "independent_classification": "nonunique_or_interval_inconclusive",
            "same_co_maximizer_count": 2,
            "same_exact_maximum": maximum,
            "co_maximizers": [
                {
                    "polynomial_ordinal": 284,
                    "location": "right_endpoint",
                    "location_ordinal": 0,
                    "parameter": "1",
                    "owned_row": 142,
                    "global_index": 143,
                    "radius_hex": "0x1.1e00000000000p+3",
                    "row_region": "smooth_interior",
                },
                {
                    "polynomial_ordinal": 285,
                    "location": "left_endpoint",
                    "location_ordinal": 0,
                    "parameter": "0",
                    "owned_row": 142,
                    "global_index": 143,
                    "radius_hex": "0x1.1e00000000000p+3",
                    "row_region": "smooth_interior",
                },
            ],
            "first_witness": {
                "polynomial_ordinal": 0,
                "component": "value_V",
                "level": "D01",
                "owned_row": 0,
                "global_index": 1,
                "subinterval": 0,
                "delta": factor,
                "polynomial_factorization": "V(t)=delta*t^2*(2*t-3)",
                "derivative_factorization": "V'(t)=6*delta*t*(t-1)",
                "exact_derivative_roots": ["0", "1"],
                "primary_interior_roots": [],
                "false_independent_root_intervals": [
                    {"root_ordinal": 0, "lower": "0", "upper": "1/2"},
                    {"root_ordinal": 1, "lower": "1/2", "upper": "1"},
                ],
            },
            "owner": "tdg9_local_extrema_independent_v1_fixed_absolute_sqrt_enclosure",
            "physical_or_model_breakdown": False,
        },
        "construction_evidence": {
            "LOC1_output_namespace_absent": True,
            "LOC1_staging_paths_absent": True,
            "LOC1_manifest_absent": True,
            "LOC1_terminal_absent": True,
            "campaign_store_leaf_count": STORE_LEAF_COUNT,
            "campaign_store_current_snapshot_sha256": STORE_SHA256,
            "PREF2_owner_result_path": PREF2_RESULT_PATH,
            "PREF2_owner_result_sha256": PREF2_RESULT_SHA256,
            "PREF2_store_manifest_sha256": PREF2_STORE_MANIFEST_SHA256,
            "current_snapshot_matches_sealed_PREF_lineage": True,
        },
        "scope": {
            "diagnosis_only": True,
            "ordinary_verification_store_blind": True,
            "LOC1_executed_successfully": False,
            "LOC2_executed": False,
            "PDE_state_committed": False,
            "fourth_width_executed": False,
            "stage2_source_decomposition_executed": False,
            "SSPRK3_comparator_executed": False,
            "candidate_branch_opened": False,
            "GR0_calibration_completed": False,
            "mechanism_result_earned": False,
            "physical_result_earned": False,
        },
    }


def _git_show(path: str) -> bytes:
    try:
        return subprocess.run(
            ("git", "show", f"{AUTHORITY_COMMIT}:{path}"),
            cwd=ROOT,
            check=True,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except subprocess.CalledProcessError as exc:
        raise ERR1Error("sealed LOC1 image is unavailable") from exc


def _snapshot_store() -> tuple[int, str]:
    base = ROOT / STORE_PATH
    metadata = base.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise ERR1Error("campaign store root is unsafe")
    digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
    count = 0
    stack = [base]
    while stack:
        directory = stack.pop()
        with os.scandir(directory) as entries:
            ordered = sorted(entries, key=lambda item: item.name)
        for item in ordered:
            observed = item.stat(follow_symlinks=False)
            if item.is_symlink():
                raise ERR1Error("campaign store contains a symlink")
            if stat.S_ISDIR(observed.st_mode):
                stack.append(Path(item.path))
                continue
            if not stat.S_ISREG(observed.st_mode):
                raise ERR1Error("campaign store contains a foreign leaf")
            path = Path(item.path)
            descriptor = os.open(
                path,
                os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            )
            try:
                active = os.fstat(descriptor)
                raw = b""
                while block := os.read(descriptor, 1 << 20):
                    raw += block
                after = os.fstat(descriptor)
            finally:
                os.close(descriptor)
            identity = (
                observed.st_dev,
                observed.st_ino,
                observed.st_size,
                observed.st_mtime_ns,
            )
            if (
                identity
                != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns)
                or identity
                != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
                or len(raw) != observed.st_size
            ):
                raise ERR1Error("campaign store leaf changed during snapshot")
            relative = path.relative_to(base).as_posix()
            digest.update((relative + "\0" + sha256(raw).hexdigest() + "\n").encode())
            count += 1
    return count, digest.hexdigest()


def _require_absence() -> None:
    target = ROOT / LOC1_NAMESPACE
    current = ROOT
    for part in target.relative_to(ROOT).parent.parts:
        current /= part
        try:
            metadata = current.lstat()
        except FileNotFoundError:
            return
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
            raise ERR1Error("LOC1 output ancestry is unsafe")
    try:
        target.lstat()
    except FileNotFoundError:
        pass
    else:
        raise ERR1Error("LOC1 output namespace exists")
    if any(item.name.startswith(LOC1_STAGING_PREFIX) for item in os.scandir(target.parent)):
        raise ERR1Error("LOC1 staging path exists")


def _reconstruct_diagnosis() -> dict[str, object]:
    """Rebuild the observed LOC1 tuple without importing the LOC2 evaluator."""

    from scripts import run_fgc_tdg9_loc1 as loc1
    from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1
    from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore
    from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (
        build_static_gr0_shells,
    )

    replay = next(item for item in ar1.REPLAYS if item["retry"] == 3)
    store = HLT16CampaignStore(ROOT / ar1.PREF2_STORE_PATH)
    prepared, _historical = loc1._prepare(
        ROOT,
        store,
        build_static_gr0_shells(ROOT),
        replay,
    )
    rows = tuple(loc1._rows(loc1._surface(prepared), "u:alpha"))
    primary_cubics, v1_cubics = loc1._cubics(rows, "D01", "value_V")
    primary = loc1.primary.localize_absolute_maximum(
        primary_cubics,
        maximum_candidates=4 * len(primary_cubics),
        refinement_depth=loc1.authority.PRIMARY_REFINEMENT_DEPTH,
    )
    v1 = loc1.independent.localize_absolute_maximum_independently(
        v1_cubics,
        maximum_candidates=4 * len(v1_cubics),
        refinement_bits=loc1.authority.INDEPENDENT_REFINEMENT_BITS,
    )
    affected = 0
    extras = 0
    first_witness: dict[str, object] | None = None
    for ordinal, (first, second) in enumerate(
        zip(primary_cubics, v1_cubics, strict=True)
    ):
        primary_roots = loc1.primary._stationary_intervals(
            first.coefficients,
            depth=loc1.authority.PRIMARY_REFINEMENT_DEPTH,
        )
        v1_roots = loc1.independent._roots(
            second.coefficients,
            loc1.authority.INDEPENDENT_REFINEMENT_BITS,
        )
        false_extras = len(v1_roots) - len(primary_roots)
        if false_extras <= 0:
            continue
        affected += 1
        extras += false_extras
        if first_witness is None:
            a0, a1, a2, a3 = first.coefficients
            delta = a3 / 2
            if a0 or a1 or a2 != -3 * delta:
                raise ERR1Error("first witness factorization differs")
            metadata = dict(first.metadata)
            first_witness = {
                "polynomial_ordinal": ordinal,
                "component": metadata["component"],
                "level": metadata["level"],
                "owned_row": metadata["owned_row"],
                "global_index": metadata["global_index"],
                "subinterval": metadata["subinterval"],
                "delta": {
                    "numerator": str(delta.numerator),
                    "denominator": str(delta.denominator),
                },
                "polynomial_factorization": "V(t)=delta*t^2*(2*t-3)",
                "derivative_factorization": "V'(t)=6*delta*t*(t-1)",
                "exact_derivative_roots": ["0", "1"],
                "primary_interior_roots": [],
                "false_independent_root_intervals": [
                    {
                        "root_ordinal": index,
                        "lower": str(lower),
                        "upper": str(upper),
                    }
                    for index, (lower, upper) in enumerate(v1_roots)
                ],
            }
    primary_maximizers = [
        {
            "polynomial_ordinal": item.polynomial_ordinal,
            "location": item.location,
            "location_ordinal": item.location_ordinal,
            "parameter": str(item.parameter_lower),
            "owned_row": dict(item.metadata)["owned_row"],
            "global_index": dict(item.metadata)["global_index"],
            "radius_hex": dict(item.metadata)["radius_hex"],
            "row_region": dict(item.metadata)["row_region"],
        }
        for item in primary.candidates
    ]
    v1_keys = [
        (item.polynomial_ordinal, item.location, item.location_ordinal)
        for item in v1.candidates
    ]
    primary_keys = [
        (item.polynomial_ordinal, item.location, item.location_ordinal)
        for item in primary.candidates
    ]
    if v1_keys != primary_keys or first_witness is None:
        raise ERR1Error("shared co-maximizer reconstruction differs")
    if not (
        primary.global_absolute_lower
        == primary.global_absolute_upper
        == v1.global_absolute_lower
        == v1.global_absolute_upper
    ):
        raise ERR1Error("shared exact maximum reconstruction differs")
    return {
        "primary_candidate_count": primary.candidate_count,
        "independent_candidate_count": v1.candidate_count,
        "affected_cubic_count": affected,
        "false_extra_candidate_count": extras,
        "primary_classification": primary.classification,
        "independent_classification": v1.classification,
        "same_co_maximizer_count": len(primary.candidates),
        "same_exact_maximum": {
            "numerator": str(primary.global_absolute_lower.numerator),
            "denominator": str(primary.global_absolute_lower.denominator),
        },
        "co_maximizers": primary_maximizers,
        "first_witness": first_witness,
    }


def _authenticate_pref_lineage() -> None:
    from recursive_horizons.fgc.evolution import tdg8_rcv3_pref2_binder as pref2

    config = (ROOT / pref2.CONFIG_PATH).read_bytes()
    result = (ROOT / PREF2_RESULT_PATH).read_bytes()
    if sha256(result).hexdigest() != PREF2_RESULT_SHA256:
        raise ERR1Error("sealed PREF2 owner result differs")
    regenerated = pref2.build_pref2_result(config, ROOT)
    if pref2.canonical_result(regenerated) != result:
        raise ERR1Error("current store does not match sealed PREF2 lineage")
    store = regenerated["artifact_payload"]["terminal_evidence"]["store"]
    if (
        store["leaf_count"] != STORE_LEAF_COUNT
        or store["manifest_sha256"] != PREF2_STORE_MANIFEST_SHA256
    ):
        raise ERR1Error("sealed PREF2 store identity differs")


def construct() -> dict[str, object]:
    head = subprocess.check_output(("git", "rev-parse", "HEAD"), cwd=ROOT, text=True).strip()
    if head != AUTHORITY_COMMIT:
        raise ERR1Error("ERR1 construction requires the sealed LOC1 authority commit")
    for path, digest in SEALED.items():
        if sha256(_git_show(path)).hexdigest() != digest:
            raise ERR1Error("sealed LOC1 hash differs")
    _require_absence()
    _authenticate_pref_lineage()
    if _snapshot_store() != (STORE_LEAF_COUNT, STORE_SHA256):
        raise ERR1Error("sealed campaign store identity differs")
    expected = _expected()
    reproduced = _reconstruct_diagnosis()
    if reproduced != {
        key: expected["diagnosis"][key]
        for key in (
            "primary_candidate_count",
            "independent_candidate_count",
            "affected_cubic_count",
            "false_extra_candidate_count",
            "primary_classification",
            "independent_classification",
            "same_co_maximizer_count",
            "same_exact_maximum",
            "co_maximizers",
            "first_witness",
        )
    }:
        raise ERR1Error("independent ERR1 reconstruction differs")
    return expected


def validate(raw: bytes) -> dict[str, object]:
    try:
        value = json.loads(raw, object_pairs_hook=_unique)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ERR1Error("invalid ERR1 JSON") from exc
    if value != _expected() or raw != _canonical(value):
        raise ERR1Error("ERR1 compact result differs")
    return value


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--verify-compact", action="store_true")
    arguments = parser.parse_args(argv)
    path = ROOT / RESULT_PATH
    if arguments.write:
        raw = _canonical(construct())
        path.write_bytes(raw)
        mode = "write"
    else:
        raw = path.read_bytes()
        validate(raw)
        mode = "verify_compact"
    print(json.dumps({"artifact_id": ARTIFACT_ID, "mode": mode, "result_sha256": sha256(raw).hexdigest(), "execution": False}, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
