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
FRONTIER_PATH = "results/nsc-2-zeta1-recursion-map.json"

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
}

PAPER_CLAIMS = {
    "NSC-1-EXACT-BLACK-UNIVERSE-DEFOCUSING": ["imported-black-universe"],
    "NSC-1-S-ONE-LOCAL-TWO-SHEET-ANOMALY": ["shared-phi-identity"],
    "NSC-2-ZETA1-FOLIATION": ["global-spectral-foliation"],
    "NSC-2-ZETA1-SELF-ADJOINT-DOMAIN": ["self-adjoint-throat"],
    "NSC-2-ZETA1-ANOMALY-DECOMPOSITION": ["anomaly-removes-root"],
    "NSC-2-ZETA1-ANOMALY-OWNER-CORRECTION": ["recursive-tail-owner"],
    "NSC-2-ZETA1-RECURSION-MAP": ["unitary-recursion-map"],
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
        match = OUTPUT_RE.search(source)
        if match is None:
            continue
        output = match.group(1)
        generators[output] = path
        dependencies[output] = sorted(
            {item for item in RESULT_RE.findall(source) if item != output}
        )
        source_dependencies[output] = sorted(set(SCRIPT_RE.findall(source)))
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
    paths = closure(dependencies)
    ordered = order_outputs(paths, dependencies)
    if len(ordered) != 58:
        raise RuntimeError(f"expected 58 public results, observed {len(ordered)}")

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
        if value.get("terminal") is not True:
            raise RuntimeError(f"result is not terminal: {output}")
        artifact_id = str(value["artifact_id"])
        generator_relative = str(generator.relative_to(ROOT))
        source = generator.read_text(encoding="utf-8")
        direct_numeric = any(name in source for name in ("numpy", "scipy", "mpmath"))
        dynamic_numeric = any(
            any(
                name in (ROOT / relative).read_text(encoding="utf-8")
                for name in ("numpy", "scipy", "mpmath")
            )
            for relative in source_dependencies.get(output, [])
        )
        inherited_numeric = any(
            dependency in portable_outputs
            for dependency in dependencies.get(output, [])
        )
        is_portable_numeric = direct_numeric or dynamic_numeric or inherited_numeric
        policy = (
            {"kind": "portable_numeric", "relative_tolerance": 1e-8, "absolute_tolerance": 1e-10}
            if is_portable_numeric
            else {"kind": "exact"}
        )
        if is_portable_numeric:
            portable_outputs.add(output)
        steps.append(
            {
                "order": index,
                "artifact_id": artifact_id,
                "category": category(artifact_id, source),
                "generator": generator_relative,
                "generator_sha256": sha256(generator),
                "output": output,
                "output_sha256": sha256(result_path),
                "dependencies": dependencies.get(output, []),
                "source_dependencies": source_dependencies.get(output, []),
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
        )

    return {
        "schema": "NSC-PUBLIC-RESULT-MANIFEST-v1",
        "source_commit": SOURCE_COMMIT,
        "frontier_artifact_id": "NSC-2-ZETA1-RECURSION-MAP",
        "frontier_output": FRONTIER_PATH,
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
