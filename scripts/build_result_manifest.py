#!/usr/bin/env python3
"""Build the curated public result graph from the frozen frontier closure."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results" / "manifest.json"
SOURCE_COMMIT = "ff2cf2722b966589b98a61accdbb6cee819a58c7"
FOLLOW_UP_SOURCE_COMMIT = "5f38712ca01ddd71e715fd265088925a73369aba"
NSC3_SOURCE_COMMIT = "445d5b069adea8b5641384b0ee02f07fe9cfcc0a"
BOUNDARY_SOURCE_COMMIT = "76975387b9b03e323c1832237e0270dc9d49a382"
NSC3_INTRODUCED_COMMIT = "3490f19164eb9303915db41b8c90eca0f40a836e"
FRONTIER_PATH = "results/nsc-2-zeta1-recursion-map.json"
FOLLOW_UP_PATH = "results/nsc-2-zeta1-unit-closure-check.json"
HISTORICAL_COUNT = 58
RELEASE_SPEC_PATH = ROOT / "results/release-spec.json"
RELEASE_SPEC = json.loads(RELEASE_SPEC_PATH.read_text(encoding="utf-8"))
SCOPED_SPECS = {row["output"]: row for row in RELEASE_SPEC["scoped_follow_ups"]}
SCOPED_FOLLOW_UP_PATHS = tuple(SCOPED_SPECS)


IMPORTED = {
    "NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING",
    "NSC-1-GRAVITATING-BPS-OBSERVATION-LINK",
    "NSC-1-S-ONE-FERMIONIC-GEOMETRY-BIND",
    "NSC-1-S-ONE-OUTSIDE-BLACK-BIND",
    "NSC-1-S-ONE-PARTICLE-WAVE-IDENTITY",
    "NSC-1-S-ONE-SHADOW-MATTER-BIND",
    "NSC-1-S-ONE-SPECTRAL-RESOLUTION-WALL-BIND",
    "NSC-1-S-ONE-SPECTRAL-ROOM-BIND",
}

DIAGNOSTIC = {
    "NSC-1-S-ONE-DOUBLED-FLRW-CLOSURE",
    "NSC-1-S-ONE-EXACT-GEOMETRY-OUTSIDE-COMPATIBILITY",
    "NSC-1-S-ONE-INTERFACE-CHARACTERISTIC-SCAN",
    "NSC-1-S-ONE-REFLECTION-POSITIVITY-OBSTRUCTION",
    "NSC-1-S-ONE-RELATIVE-REFLECTION-TEST",
    "NSC-1-S-ONE-SHADOW-LENSING-OBSERVATION",
    "NSC-2-ZETA1-ANGULAR-TOWER",
    "NSC-2-ZETA1-ANOMALY-DECOMPOSITION",
    "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
}

SUPERSEDED = {
    "NSC-1-S-ONE-BOUNDARY-RETARDED-POLE",
    "NSC-1-S-ONE-CHILD-ORIENTATION",
    "NSC-1-S-ONE-GAP-TRANSITION-CLOSURE",
    "NSC-2-ZETA1-LOWEST-MODE",
    "NSC-2-ZETA1-ORBIFOLD-PARITY",
    "NSC-2-ZETA1-REGULATED-DETERMINANT",
    "NSC-2-ZETA1-WARPED-Y",
    "NSC-2-ZETA1-Y-BOUNDARY-SENSITIVITY",
}

HEADLINES = {
    "NSC-1-GRAVITATING-BPS-OBSERVATION-LINK": [
        ("/self_gravitating_predictions/full_field_maximum_mass_solar", 0.0, 0.0),
    ],
    "NSC-1-S-ONE-GAUSS-BONNET-KINETIC-RANK": [
        ("/result/critical_local_rank_all_samples/0", 0.0, 0.0),
    ],
    "NSC-1-S-ONE-INDUCED-INTERFACE-RANK": [
        ("/result/combined_interface_ranks/0", 0.0, 0.0),
    ],
    "NSC-1-S-ONE-INTERFACE-CHARACTERISTIC-SCAN": [
        ("/method/normalization_grid/count", 0.0, 0.0),
    ],
    "NSC-1-S-ONE-RELATIVE-REFLECTION-TEST": [
        ("/reflection_spectra/3/minimum_eigenvalue", 1e-8, 1e-4),
    ],
    "NSC-2-ZETA1-LOWEST-MODE": [
        ("/first_scale_estimate/zeta", 1e-8, 2e-5),
        ("/first_scale_estimate/mu", 1e-8, 2e-5),
    ],
    "NSC-2-ZETA1-REGULATED-DETERMINANT": [
        ("/current_scale_candidate/zeta", 1e-8, 2e-5),
        ("/current_scale_candidate/mu", 1e-8, 2e-5),
    ],
    "NSC-2-ZETA1-ORBIFOLD-PARITY": [
        ("/current_unwarped_orbifold_candidate/zeta", 1e-8, 2e-6),
        ("/current_unwarped_orbifold_candidate/mu", 1e-8, 2e-5),
    ],
    "NSC-2-ZETA1-WARPED-Y": [
        ("/current_warped_y_candidate/zeta", 1e-8, 5e-6),
        ("/current_warped_y_candidate/mu", 1e-8, 5e-5),
    ],
    "NSC-2-ZETA1-ANOMALY-DECOMPOSITION": [
        ("/prior_determinant_candidate/zeta", 1e-8, 5e-6),
        (
            "/prior_determinant_candidate/anomaly_compensated_child_link_derivative",
            1e-8,
            1e-4,
        ),
        ("/scan/minimum_physical_derivative", 1e-8, 1e-4),
        ("/scan/maximum_physical_derivative", 1e-8, 1e-5),
    ],
    "NSC-2-ZETA1-UNIT-CLOSURE-CHECK": [
        (
            "/comparison_at_frozen_old_candidate/corrected_subtracted_derivative",
            1e-8,
            1e-8,
        ),
        ("/corrected_determinant_diagnostic/zeta", 1e-8, 1e-8),
        ("/corrected_scan/maximum_subtracted_derivative", 1e-8, 1e-8),
    ],
    "NSC-3-GEOMETRIC-CHAIN": [
        ("/gap_families/0/continuum/band_edge", 1e-8, 1e-10),
        ("/gap_families/1/continuum/band_edge", 1e-8, 1e-10),
        ("/gap_families/2/continuum/band_edge", 1e-8, 1e-10),
    ],
}

PAPER_CLAIMS = {
    "NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING": ["imported-black-universe"],
    "NSC-1-S-ONE-LOCAL-TWO-SHEET-ANOMALY": ["shared-phi-identity"],
    "NSC-2-ZETA1-FOLIATION": ["global-spectral-foliation"],
    "NSC-2-ZETA1-SELF-ADJOINT-DOMAIN": ["self-adjoint-throat"],
    "NSC-2-ZETA1-ANOMALY-DECOMPOSITION": ["anomaly-removes-root"],
    "NSC-2-ZETA1-ANOMALY-OWNER-CORRECTION": ["recursive-tail-owner"],
    "NSC-2-ZETA1-RECURSION-MAP": ["unitary-recursion-map"],
    "NSC-2-ZETA1-UNIT-CLOSURE-CHECK": [
        "unit-consistent-mass-term",
        "finite-family-sign-theorem",
        "raw-versus-regulated-determinant",
    ],
    "NSC-3-REGULATED-RECURSION": [
        "finite-regulated-variations",
        "ultrastatic-frequency-factor",
        "nonelliptic-coordinate-time",
    ],
    "NSC-3-RADIAL-SPECTRUM": ["isolated-radial-gapless"],
    "NSC-3-GEOMETRIC-CHAIN": ["periodic-throat-gap", "stencil-derived-link"],
    "NSC-12-COMPACT-INTERACTION": ["compact-mode-overlaps"],
    "NSC-13-TORSION-UV-MAP": ["torsion-uv-finite-matching"],
    "NSC-14-FLOW-COMPATIBILITY": ["published-flow-compatibility"],
    "NSC-15-CHARGED-SECTOR": ["charged-sector-parity"],
}

OUTPUT_RE = re.compile(r"OUTPUT\s*=\s*ROOT\s*/\s*[\"']([^\"']+)[\"']")
RESULT_RE = re.compile(r"ROOT\s*/\s*[\"'](results/[^\"']+\.json)[\"']")
SCRIPT_RE = re.compile(r"ROOT\s*/\s*[\"'](scripts/[^\"']+\.py)[\"']")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def discover() -> tuple[dict[str, Path], dict[str, list[str]], dict[str, list[str]]]:
    generators: dict[str, Path] = {}
    dependencies: dict[str, list[str]] = {}
    source_dependencies: dict[str, list[str]] = {}
    for path in sorted((ROOT / "scripts").glob("*.py")):
        source = path.read_text(encoding="utf-8")
        if path.name == "check_nsc_boundary_response.py":
            output = "results/nsc-3-boundary-response.json"
            generators[output] = path
            dependencies[output] = []
            source_dependencies[output] = []
            continue
        match = OUTPUT_RE.search(source)
        if match is None:
            continue
        output = match.group(1)
        generators[output] = path
        dependencies[output] = sorted(
            {item for item in RESULT_RE.findall(source) if item != output}
        )
        source_dependencies[output] = sorted(set(SCRIPT_RE.findall(source)))
    for output, row in SCOPED_SPECS.items():
        generators[output] = ROOT / row["generator"]
        dependencies[output] = row["dependencies"]
        source_dependencies[output] = row["source_dependencies"]
    return generators, dependencies, source_dependencies


def closure(dependencies: dict[str, list[str]]) -> set[str]:
    seen: set[str] = set()

    def visit(path: str) -> None:
        if path in seen:
            return
        for dependency in dependencies.get(path, []):
            visit(dependency)
        seen.add(path)

    visit(FRONTIER_PATH)
    return seen


def order_outputs(paths: set[str], dependencies: dict[str, list[str]]) -> list[str]:
    remaining = set(paths)
    ordered: list[str] = []
    complete: set[str] = set()
    while remaining:
        ready = sorted(
            path
            for path in remaining
            if set(dependencies.get(path, [])) <= complete
        )
        if not ready:
            raise RuntimeError("result graph contains a cycle or undeclared dependency")
        ordered.extend(ready)
        complete.update(ready)
        remaining.difference_update(ready)
    return ordered


def category(artifact_id: str, generator_source: str) -> str:
    if artifact_id == "NSC-3-BOUNDARY-RESPONSE":
        return "repository_derived_numerical_result"
    if artifact_id == "NSC-2-ZETA1-RECURSION-MAP":
        return "current_frontier"
    if artifact_id in IMPORTED:
        return "imported_benchmark_reproduction"
    if artifact_id in DIAGNOSTIC:
        return "diagnostic_nonpass"
    if artifact_id in SUPERSEDED:
        return "superseded_candidate"
    if any(name in generator_source for name in ("numpy", "scipy", "mpmath")):
        return "repository_derived_numerical_result"
    return "repository_derived_exact_identity"



def build() -> dict[str, object]:
    generators, dependencies, source_dependencies = discover()
    missing_follow_ups = [path for path in SCOPED_FOLLOW_UP_PATHS if path not in generators]
    if missing_follow_ups:
        raise RuntimeError(f"follow-up generator is absent: {missing_follow_ups}")
    historical_dependencies = {
        path: deps
        for path, deps in dependencies.items()
        if path not in SCOPED_FOLLOW_UP_PATHS
    }
    paths = closure(historical_dependencies)
    paths.difference_update(SCOPED_FOLLOW_UP_PATHS)
    ordered = order_outputs(paths, historical_dependencies)
    if len(ordered) != HISTORICAL_COUNT:
        raise RuntimeError(
            f"expected {HISTORICAL_COUNT} historical public results, observed {len(ordered)}"
        )
    pending = list(SCOPED_FOLLOW_UP_PATHS)
    while pending:
        ready = [path for path in pending if set(dependencies[path]) <= set(ordered)]
        if not ready:
            raise RuntimeError("scoped result graph contains a cycle or absent dependency")
        # The explicit release order breaks ties, preserving published checkpoints.
        path = ready[0]
        ordered.append(path)
        pending.remove(path)

    steps: list[dict[str, object]] = []
    portable_outputs: set[str] = set()
    for index, output in enumerate(ordered, start=1):
        generator = generators.get(output)
        if generator is None:
            raise RuntimeError(f"no generator for {output}")
        result_path = ROOT / output
        if not result_path.is_file():
            raise RuntimeError(f"result is absent: {output}")
        value = json.loads(result_path.read_text(encoding="utf-8"))
        identity = SCOPED_SPECS.get(output, {}).get("identity_policy", {"kind": "terminal"})
        if identity["kind"] == "schema_gate":
            gate = value.get("gate")
            if (value.get("schema") != identity["schema"]
                    or value.get("classification") != identity["classification"]
                    or not isinstance(gate, dict) or set(gate) != set(identity["gate"])
                    or any(gate[key] is not expected for key, expected in identity["gate"].items())
                    or "artifact_id" in value or "terminal" in value):
                raise RuntimeError(f"invalid schema-gated result: {output}")
        elif identity["kind"] == "schema_status":
            if (value.get("schema") != identity["schema"]
                    or value.get("status") != identity["status"]
                    or "artifact_id" in value or "terminal" in value):
                raise RuntimeError(f"invalid schema-status result: {output}")
        elif identity["kind"] == "plateau_scope":
            if (value.get("schema") != identity["schema"]
                    or value.get("scope") != identity["scope"]
                    or value.get("observational_audit", {}).get("prediction_claim") is not False):
                raise RuntimeError(f"invalid scoped plateau record: {output}")
        elif value.get("terminal") is not True:
            raise RuntimeError(f"result is not terminal: {output}")
        artifact_id = SCOPED_SPECS.get(output, {}).get("artifact_id", value.get("artifact_id"))
        if not isinstance(artifact_id, str) or not artifact_id:
            raise RuntimeError(f"result has no declared artifact identity: {output}")
        generator_relative = str(generator.relative_to(ROOT))
        source = generator.read_text(encoding="utf-8")
        is_follow_up = output in SCOPED_FOLLOW_UP_PATHS
        declared_dependencies = dependencies.get(output, [])
        declared_source_dependencies = sorted(set(source_dependencies.get(output, [])))
        direct_numeric = any(name in source for name in ("numpy", "scipy", "mpmath"))
        dynamic_numeric = any(
            any(
                name in (ROOT / relative).read_text(encoding="utf-8")
                for name in ("numpy", "scipy", "mpmath")
            )
            for relative in declared_source_dependencies
        )
        inherited_numeric = any(
            dependency in portable_outputs for dependency in declared_dependencies
        )
        is_portable_numeric = (
            False
            if is_follow_up
            else direct_numeric or dynamic_numeric or inherited_numeric
        )
        if is_follow_up:
            policy = {
                "kind": "all_fields",
                "relative_tolerance": 1e-8,
                "absolute_tolerance": 1e-8,
            }
        elif is_portable_numeric:
            policy = {
                "kind": "portable_numeric",
                "relative_tolerance": 1e-8,
                "absolute_tolerance": 1e-10,
            }
            portable_outputs.add(output)
        else:
            policy = {"kind": "exact"}
        step = {
            "order": index,
            "artifact_id": artifact_id,
            "category": category(artifact_id, source),
            "generator": generator_relative,
            "generator_sha256": sha256(generator),
            "output": output,
            "output_sha256": sha256(result_path),
            "dependencies": declared_dependencies,
            "source_dependencies": declared_source_dependencies,
            "comparison_policy": policy,
            "headline_observables": [
                {
                    "pointer": pointer,
                    "relative_tolerance": relative_tolerance,
                    "absolute_tolerance": absolute_tolerance,
                }
                for pointer, relative_tolerance, absolute_tolerance in HEADLINES.get(
                    artifact_id, []
                )
            ],
            "paper_claim_ids": PAPER_CLAIMS.get(artifact_id, []),
        }
        if is_follow_up:
            row = SCOPED_SPECS[output]
            for key in ("comparison_policy", "generator_args", "json_format",
                        "follow_up_source_commit", "auxiliary_inputs", "identity_policy", "category"):
                if key in row:
                    step[key] = row[key]
            step["source_dependency_hashes"] = {
                relative: sha256(ROOT / relative)
                for relative in declared_source_dependencies
            }
            for auxiliary in row.get("auxiliary_inputs", []):
                if sha256(ROOT / auxiliary["path"]) != auxiliary["sha256"]:
                    raise RuntimeError(f"auxiliary provenance mismatch: {auxiliary['path']}")
        steps.append(step)

    return {
        "schema": "NSC-PUBLIC-RESULT-MANIFEST-v1",
        "release_version": RELEASE_SPEC["release_version"],
        "release_spec": "results/release-spec.json",
        "release_spec_sha256": sha256(RELEASE_SPEC_PATH),
        "release_source_commit": RELEASE_SPEC["source_commit"],
        "source_commit": SOURCE_COMMIT,
        "follow_up_source_commit": FOLLOW_UP_SOURCE_COMMIT,
        "nsc3_source_commit": NSC3_SOURCE_COMMIT,
        "boundary_source_commit": BOUNDARY_SOURCE_COMMIT,
        "nsc3_introduced_commit": NSC3_INTRODUCED_COMMIT,
        "frontier_artifact_id": "NSC-2-ZETA1-RECURSION-MAP",
        "frontier_output": FRONTIER_PATH,
        "follow_up_artifact_id": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
        "follow_up_output": FOLLOW_UP_PATH,
        "scoped_follow_up_outputs": ordered[HISTORICAL_COUNT:],
        "historical_result_count": HISTORICAL_COUNT,
        "measured_reference_runtime_seconds": 246.73,
        "result_count": len(steps),
        "steps": steps,
    }


def main() -> int:
    OUTPUT.write_text(
        json.dumps(build(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
