#!/usr/bin/env python3
"""Regenerate the FGC-1-VAR1 covariant metric-variation certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import sys
import tomllib
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.metric_variation import (  # noqa: E402
    EXPANDED_GB_COEFFICIENTS,
    generate_algebraic_curvature_fixtures,
    metric_variation_certificate,
    required_var1_nonclaims,
)
from scripts.reproduce_fgc_action import load_config as load_action_config  # noqa: E402

DEFAULT_CONFIG = REPOSITORY / "configs" / "fgc" / "fgc-1-metric-variation.toml"
DEFAULT_OUTPUT = REPOSITORY / "results" / "fgc-1-metric-variation.json"
ACTION_SOURCE = REPOSITORY / "src" / "recursive_horizons" / "fgc" / "action.py"
VARIATION_SOURCE = (
    REPOSITORY / "src" / "recursive_horizons" / "fgc" / "metric_variation.py"
)
VARIATION_DOCUMENT = REPOSITORY / "docs" / "fgc-metric-variation.md"


@dataclass(frozen=True, slots=True)
class VariationConfig:
    schema_version: int
    artifact_id: str
    project_version: str
    metric_signature: str
    riemann_convention: str
    action_config: Path
    action_config_relative: str
    dimension: int
    algorithm: str
    seed: int
    fixture_count: int
    integer_bound: int
    kulkarni_nomizu_pair_count: int
    source_sha256: str


def _mapping(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a table")
    return value


def _exact_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    actual = set(value)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"{name} keys differ: missing={missing}, extra={extra}")


def _positive_integer(name: str, value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _repository_path(name: str, value: Any) -> tuple[Path, str]:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{name} must be a non-empty repository-relative path")
    relative = Path(value)
    if relative.is_absolute():
        raise ValueError(f"{name} must be repository relative")
    resolved = (REPOSITORY / relative).resolve()
    try:
        canonical = resolved.relative_to(REPOSITORY)
    except ValueError as exc:
        raise ValueError(f"{name} must stay inside the repository") from exc
    if canonical != relative:
        raise ValueError(f"{name} must be canonical and traversal free")
    if not resolved.is_file():
        raise ValueError(f"{name} does not exist")
    return resolved, canonical.as_posix()


def load_config(path: Path = DEFAULT_CONFIG) -> VariationConfig:
    """Load and strictly validate the frozen VAR1 convention fixture."""

    path = path.resolve()
    source = path.read_bytes()
    raw = _mapping("root", tomllib.loads(source.decode("utf-8")))
    _exact_keys(
        "root",
        raw,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "action_config",
            "fixture_generator",
        },
    )
    if raw["schema_version"] != 1:
        raise ValueError("schema_version must equal 1")
    if raw["artifact_id"] != "FGC-1-VAR1":
        raise ValueError("artifact_id must equal FGC-1-VAR1")
    if raw["project_version"] != "0.11.0":
        raise ValueError("project_version must equal 0.11.0")
    if raw["metric_signature"] != "-+++":
        raise ValueError("metric_signature must equal -+++")
    if raw["riemann_convention"] != "plus_partial_mu_gamma_nu":
        raise ValueError(
            "riemann_convention must equal plus_partial_mu_gamma_nu"
        )
    action_path, action_relative = _repository_path(
        "action_config", raw["action_config"]
    )
    if action_relative != "configs/fgc/fgc-1-action-gate.toml":
        raise ValueError("action_config must name the frozen ACT1 configuration")

    generator = _mapping("fixture_generator", raw["fixture_generator"])
    _exact_keys(
        "fixture_generator",
        generator,
        {
            "dimension",
            "algorithm",
            "seed",
            "fixture_count",
            "integer_bound",
            "kulkarni_nomizu_pair_count",
        },
    )
    if generator["dimension"] != 4:
        raise ValueError("fixture_generator.dimension must equal 4")
    if generator["algorithm"] != "lcg64_v1":
        raise ValueError("fixture_generator.algorithm must equal lcg64_v1")
    seed = generator["seed"]
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("fixture_generator.seed must be a nonnegative integer")
    fixture_count = _positive_integer(
        "fixture_generator.fixture_count", generator["fixture_count"]
    )
    integer_bound = _positive_integer(
        "fixture_generator.integer_bound", generator["integer_bound"]
    )
    pair_count = _positive_integer(
        "fixture_generator.kulkarni_nomizu_pair_count",
        generator["kulkarni_nomizu_pair_count"],
    )
    if pair_count != 2:
        raise ValueError(
            "fixture_generator.kulkarni_nomizu_pair_count must equal 2"
        )

    action = load_action_config(action_path)
    if action.project_version != raw["project_version"]:
        raise ValueError("ACT1 and VAR1 project versions differ")
    if action.metric_signature != raw["metric_signature"]:
        raise ValueError("ACT1 and VAR1 metric signatures differ")

    return VariationConfig(
        schema_version=1,
        artifact_id="FGC-1-VAR1",
        project_version="0.11.0",
        metric_signature="-+++",
        riemann_convention="plus_partial_mu_gamma_nu",
        action_config=action_path,
        action_config_relative=action_relative,
        dimension=4,
        algorithm="lcg64_v1",
        seed=seed,
        fixture_count=fixture_count,
        integer_bound=integer_bound,
        kulkarni_nomizu_pair_count=2,
        source_sha256=sha256(source).hexdigest(),
    )


def _source_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, object]:
    """Build the source-bound VAR1 result from exact rational fixtures."""

    config_path = config_path.resolve()
    config = load_config(config_path)
    action = load_action_config(config.action_config)
    fixtures = generate_algebraic_curvature_fixtures(
        config.seed, config.fixture_count, config.integer_bound
    )
    certificate = metric_variation_certificate(action.models, fixtures)
    checks = certificate["verified_exact_checks"]
    if not all(value is True for value in checks.values()):
        raise ValueError("VAR1 exact identity certificate did not pass")
    if certificate["aggregate"]["maximum_exact_residual"] != "0":
        raise ValueError("VAR1 aggregate residual is nonzero")
    nonclaims = required_var1_nonclaims()
    try:
        source_config = config_path.relative_to(REPOSITORY).as_posix()
    except ValueError:
        source_config = str(config_path)

    return {
        "schema_version": config.schema_version,
        "project_version": config.project_version,
        "artifact_id": config.artifact_id,
        "artifact": "fgc1_covariant_metric_variation_and_exact_identity_certificate",
        "classification": "covariant_metric_equation_derivation_and_crosschecks_not_principal_symbol_or_collapse",
        "development_status": "pre_publication_gate",
        "source_config": source_config,
        "source_configs_sha256": {
            source_config: config.source_sha256,
            config.action_config_relative: _source_digest(config.action_config),
        },
        "implementation_sources_sha256": {
            "src/recursive_horizons/fgc/action.py": _source_digest(ACTION_SOURCE),
            "src/recursive_horizons/fgc/metric_variation.py": _source_digest(
                VARIATION_SOURCE
            ),
            "scripts/reproduce_fgc_metric_variation.py": _source_digest(
                Path(__file__).resolve()
            ),
        },
        "derivation_document": "docs/fgc-metric-variation.md",
        "derivation_document_sha256": _source_digest(VARIATION_DOCUMENT),
        "fixture_contract": {
            "dimension": config.dimension,
            "metric_signature": config.metric_signature,
            "riemann_convention": config.riemann_convention,
            "algorithm": config.algorithm,
            "seed": config.seed,
            "fixture_count": config.fixture_count,
            "integer_bound": config.integer_bound,
            "kulkarni_nomizu_pair_count": config.kulkarni_nomizu_pair_count,
        },
        "variation_certificate": certificate,
        "primary_source_crosscheck": {
            "source": "Thaalba_et_al_arXiv_2306.01695_equations_1_and_3",
            "url": "https://arxiv.org/abs/2306.01695",
            "mapping": "f_T=alpha_T*phi^2/2_and_double_epsilon_term_equals_-8*P_acbd*nabla^c*nabla^d*f_T",
            "coefficient_and_sign_match": True,
            "unsafe_normalization_warning": "Lara_et_al_arXiv_2403.08705_equation_6_is_not_used_as_coefficient_authority",
        },
        "gate_status": {
            "act1_action_and_scalar_algebra_certificate_passed": True,
            "var1_covariant_metric_equation_derived_and_cross_checked": True,
            "full_fgc1_health_gate_passed": False,
            "publication_local_defocusing_gate_passed": False,
        },
        "nonclaims": nonclaims,
        "scope": "Derives the complete covariant bulk metric equation for the frozen action and cross-checks the Gauss-Bonnet sign, factor, expansion, trace, Noether curvature-contraction identity, nonminimal source, branch specializations, and exact automatic FLRW lapse variation. It does not numerically evaluate the differential Bianchi identity or full divergence, and does not derive boundary or junction terms, the spherical system, principal symbol, constraints, nonlinear health, collapse, affine defocusing, singularity resolution, child spacetime, dark sector, variable c, particle spectrum, or a theory of everything.",
        "expanded_gb_coefficient_order": list(EXPANDED_GB_COEFFICIENTS),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(
        record(arguments.config), indent=2, sort_keys=True, allow_nan=False
    ) + "\n"
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
