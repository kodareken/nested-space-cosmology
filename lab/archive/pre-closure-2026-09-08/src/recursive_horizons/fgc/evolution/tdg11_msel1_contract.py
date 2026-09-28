"""Prospective TDG11-MSEL1 data contract and no-trajectory algebra controls.

No Git, store, NumPy, proposer, or runner is imported here. A valid config is
not execution authority: the separately committed-image owner must authorize
the one-shot diagnostic. All thresholds and candidates precede measurement.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import json
import re
import tomllib

from recursive_horizons.evidence_io import canonical_json_bytes
from .tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    three_halves_order_passes_squared,
)


ARTIFACT_ID = "FGC-1-TDG11-MSEL1-FRZ1"
CLASSIFICATION = "prospective_three_route_temporal_selection_no_measurement"
BASE_COMMIT = "dbfe92df40bf18fd287663e0b9997a812fe4ccb2"
CONFIG_PATH = "configs/fgc/fgc-1-tdg11-msel1-frz1.toml"
RESULT_PATH = "results/fgc-1-tdg11-msel1-frz1.json"
OWNER_DOCUMENT = "docs/fgc-tdg11-msel1-frz1.md"
OUTPUT_NAMESPACE = "runs/fgc-2-sf1/tdg11-msel1"
STORE_PATH = "runs/fgc-2-sf1/tdg8-rcv3/calibration"
STORE_LEAF_COUNT = 115
STORE_SHA256 = "5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445"
MEMBER_KEY = "RK4-2049"
POINT_COUNT = 2049
OWNED_ROW_COUNT = 2044
ACCEPTED_TIME_HEX = "0x1.78554de5a30e0p+0"
PHYSICAL_STATE_SHA256 = (
    "3cd9f576d079595343e8cdaabd33dd0153763b88d1beecd7d85717bdba4b6c9a"
)
DESCRIPTOR_SHA256 = "77847126340e78c8ac300fac2795bceda724c15e0894e34444b0da42a84166b4"
COORDINATES_SHA256 = "2a347b314de5e1e94ac0c999db283e03498070ece6762c53252fa41c1e04a646"
ORIGINAL_ARITHMETIC_ID = "tdg11_recorded_legacy_binary64_rk_v1"
COMPENSATED_ARITHMETIC_ID = "tdg11_binary64_twofold_dot2_coherent_rk_v1"
CANDIDATES = (
    "exact_accumulation_reconstruction_with_debit",
    "compensated_coherent_RK_accumulation",
    "embedded_RK_quadrature_defect",
)

# These immutable references own the copied retry/origin/environment facts.
PINNED_REFERENCES = {
    "results/fgc-1-tdg10-qa2-pref1.json": "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6",
    "results/fgc-1-pro18-pref27.json": "102b385d058066d668e9e8541c974c736c583f6d5c4afaf9b99e95a16aabeb83",
    "configs/fgc/fgc-1-pro19-frz1.toml": "e4dac87e91d96802830418d624a4328b247bd6e4376a5b3de632c7f26ae66988",
    "configs/fgc/fgc-1-tdg10-qa1-frz1.toml": "0d171c3caea98c4efa9f44c90e666d8a665ecc99ab273f2470d57e5d297f64b3",
    "configs/fgc/fgc-1-tdg10-qa2-frz1.toml": "5d32295113d8c675774319b3313426140b0ca4902f956541182bbc40c2949323",
    "configs/fgc/fgc-runtime-matrix.toml": "66b760bc042e3bc24c67ba7b7e0769eed32448fb32f3caf7d570a83e6dfabe68",
}
_REPLAY_ROWS = (
    (
        3,
        9,
        10,
        11,
        "0x1.aaa9612df8000p-11",
        "eb6fddc480c94aa7ed15c80399fc2b26637693efef02399f267075fc6c258e56",
        "07572d4450829ff2e2bad65e217a917b0d31ba677a3b8fc20a3be43e849ba846",
        "5b533eb7009a9c7c3f353d26813f9cbbde8574afae77492929b435b4c701c5ff",
        "b50dc39ebd8bbdd3729d40d9e9ae4b223ff2ef19389d4ca74f9a316bac72202f",
        "341cd8cd328434d85774bcb452889ac1886c520a72e945a92337b72cf564161a",
        "e4f9ca42f2d014d4bdbd875a0f24a9bbb9c04159fb0d9e0e7de8d0d72db61fdf",
    ),
    (
        4,
        10,
        12,
        13,
        "0x1.aaa9612df8000p-12",
        "6dc263e7719c9a422e9fed1b81ac573f3127607a2285f3c55b095d9019f60187",
        "9abb59999809d22354d6b5f810bc592cd2673e5e4d354d0ec58e9ad53522dc4b",
        "d25d371b67cd227edc735479f719970b382443d9097ee718a4dd3888c0b6cca4",
        "dd02233dd5e4cf28becbc4fe29259d0879005d806b70493aa2281255081eb219",
        "1bd21553f3b6e0af613997da094a46d4dcd966a0e4346692bae62b05046862ca",
        "26542403a17b05664d5a447b88860186d5f554b1f4d555624493ed40b941694c",
    ),
    (
        5,
        11,
        14,
        15,
        "0x1.aaa9612df0000p-13",
        "11ec8a80d300de72075661a97b616aecead8ec78208bb541b93e62ca4c8124c8",
        "786bdcb835e043c65ddd33fad3a2af0e3d33e768dc26a9f9f00d2757b0da3ea5",
        "0810018134df958b8dab32bd9c0ecd7ebfc4376fbc87869a8bee7dd02eb1e199",
        "a1882d904091591d7fa3ccf3eab04aa748161115f97c15abd57e535e6fa01f3d",
        "475014d918952adef032d71b81073773364a37fa795491e739aae609f77e8537",
        "1997342415f74f330ce6b45ca6fd550ad0e5e674b85148fea5047ff82ac3f152",
    ),
)
_REPLAY_KEYS = (
    "retry",
    "generation",
    "journal_tip_sequence",
    "rejection_sequence",
    "width_hex",
    "checkpoint_sha256",
    "checkpoint_raw_sha256",
    "journal_tip_sha256",
    "journal_tip_raw_sha256",
    "rejection_sha256",
    "rejection_raw_sha256",
)

IMPLEMENTATION_PATHS = tuple(
    sorted(
        (
            "src/recursive_horizons/evidence_io.py",
            *(
                f"src/recursive_horizons/fgc/evolution/{name}.py"
                for name in (
                    "tdg11_msel1_contract",
                    "tdg11_compensated_rk",
                    "tdg11_rational_complete_c",
                    "tdg11_msel1_reconstruction",
                    "tdg11_msel1_runtime",
                    "tdg11_msel1_authority",
                )
            ),
            "scripts/run_fgc_tdg11_msel1.py",
            "scripts/reproduce_fgc_tdg11_msel1_frz1.py",
            *(
                f"tests/test_fgc_tdg11_{name}.py"
                for name in (
                    "compensated_rk",
                    "rational_complete_c",
                    "msel1_reconstruction",
                    "msel1_contract",
                    "msel1_runtime",
                    "msel1_authority",
                    "msel1_runner",
                )
            ),
        )
    )
)
MINIMUM_ENVIRONMENT = {
    "python_implementation": "CPython",
    "python_version": "3.14.3",
    "numpy_version": "2.5.1",
    "blas_name": "accelerate",
    "blas_version": "unknown",
    "system": "Darwin",
    "machine": "arm64",
    "byteorder": "little",
}
DELTA_PATHS = tuple(
    sorted(
        set(IMPLEMENTATION_PATHS)
        | {
            CONFIG_PATH,
            RESULT_PATH,
            OWNER_DOCUMENT,
            "README.md",
            "docs/active-code-map.md",
            "docs/claim-ledger.md",
            "docs/research-roadmap.md",
            "docs/fgc-runtime-matrix.md",
            "docs/evidence-io.md",
            "results/README.md",
            "configs/fgc/artifact-catalog.json",
            "scripts/build_artifact_catalog.py",
            "mk/current-foundation.mk",
            "scripts/check_repo.py",
            "scripts/repo_checks/core.py",
            "scripts/repo_checks/catalog.py",
            "scripts/repo_checks/temporal_selection.py",
            "tests/test_evidence_io.py",
            "tests/test_check_repo_tdg11_msel1_frz1.py",
            "tests/test_phase_minus1_artifact_catalog.py",
        }
    )
)
ENVIRONMENT_IMAGE_KEYS = (
    "kernel_release",
    "kernel_version",
    "macos_version",
    "python_executable",
    "python_executable_sha256",
    "numpy_extension",
    "numpy_extension_sha256",
)
RESOURCES = {
    "selected_widths": 3,
    "channels_per_width": 18,
    "candidate_width_channel_records": 162,
    "baseline_width_channel_records": 54,
    "path_families": 6,
    "level_paths": 18,
    "proposals": 42,
    "stage_and_endpoint_records": 210,
    "accepted_source_prechecks": 42,
    "static_shells": 6,
    "static_source_prechecks": 6,
    "per_channel_maximum_candidates_D01": 16352,
    "per_channel_maximum_candidates_D12": 32704,
    "primary_refinement_depth": 160,
    "raw_leaf_count": 2,
    "maximum_raw_bytes": 16 * 1024 * 1024,
    "maximum_wall_seconds": 14400,
    "maximum_rational_bits": 32768,
    "git_listing_entry_ceiling": 16384,
}
_NONCLAIMS = (
    "endpoint_serialization_authorized",
    "campaign_state_write_authorized",
    "accepted_state_advance_authorized",
    "production_method_selected",
    "production_admission_authorized",
    "temporal_retry_authorized",
    "threshold_fitting_authorized",
    "fourth_width_authorized",
    "diagnostic_endpoint_transplantation_authorized",
    "HLT15_origin_modified",
    "GR0_calibration_authorized",
    "SGBL_execution_authorized",
    "FGCQR_execution_authorized",
    "holdout_opened",
    "retained_EFT_evolution_authorized",
    "global_PDE_error_theorem_earned",
    "mechanism_result_earned",
    "physical_result_earned",
    "publication_authorized",
    "push_authorized",
)
_HEX = re.compile(r"[0-9a-f]{64}\Z")


class MSEL1ContractError(ValueError):
    """An input changes the prospective selection contract."""


def replay_specs() -> tuple[dict[str, object], ...]:
    return tuple(
        dict(zip(_REPLAY_KEYS, row, strict=True)) | {"prior_retry_count": row[0] - 1}
        for row in _REPLAY_ROWS
    )


def nonclaims() -> dict[str, bool]:
    return dict.fromkeys(_NONCLAIMS, False)


def _strict_equal(value: object, expected: object, path: str = "config") -> None:
    if type(value) is not type(expected):
        raise MSEL1ContractError(f"{path} has a different type")
    if type(expected) is dict:
        if value.keys() != expected.keys():
            raise MSEL1ContractError(f"{path} keys differ")
        for key in expected:
            _strict_equal(value[key], expected[key], f"{path}.{key}")
    elif type(expected) is list:
        if len(value) != len(expected):
            raise MSEL1ContractError(f"{path} length differs")
        for index, (actual, wanted) in enumerate(zip(value, expected, strict=True)):
            _strict_equal(actual, wanted, f"{path}[{index}]")
    elif value != expected:
        raise MSEL1ContractError(f"{path} differs")


def _environment_record(value: object) -> dict[str, str]:
    keys = set(MINIMUM_ENVIRONMENT) | set(ENVIRONMENT_IMAGE_KEYS)
    if type(value) is not dict or set(value) != keys:
        raise MSEL1ContractError("environment image keys differ")
    if any(
        type(item) is not str or not item or len(item) > 1024 or "\0" in item
        for item in value.values()
    ):
        raise MSEL1ContractError("environment image contains invalid text")
    for key, expected in MINIMUM_ENVIRONMENT.items():
        if value[key] != expected:
            raise MSEL1ContractError(f"environment {key} differs")
    for key in ("python_executable_sha256", "numpy_extension_sha256"):
        if not _HEX.fullmatch(value[key]):
            raise MSEL1ContractError(f"environment {key} is not a SHA-256")
    return dict(value)


def _implementation_records(value: object) -> list[dict[str, str]]:
    if type(value) is not list or len(value) != len(IMPLEMENTATION_PATHS):
        raise MSEL1ContractError("implementation binding count differs")
    answer = []
    for item, path in zip(value, IMPLEMENTATION_PATHS, strict=True):
        if (
            type(item) is not dict
            or set(item) != {"path", "sha256"}
            or type(item["path"]) is not str
            or item["path"] != path
            or type(item["sha256"]) is not str
            or not _HEX.fullmatch(item["sha256"])
        ):
            raise MSEL1ContractError("implementation path/hash inventory differs")
        answer.append(dict(item))
    return answer


def expected_config(
    *, environment: dict[str, str], implementation: list[dict[str, str]]
) -> dict[str, object]:
    """Build only the frozen metadata schema; this grants no run capability."""

    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": "0.11.0",
        "base_commit": BASE_COMMIT,
        "owner_document": OWNER_DOCUMENT,
        "output_namespace": OUTPUT_NAMESPACE,
        "predecessor_artifact_id": "FGC-1-TDG10-QA2-PREF1",
        "authority_delta_paths": list(DELTA_PATHS),
        "environment": _environment_record(environment),
        "implementation": _implementation_records(implementation),
        "references": [
            {"path": path, "sha256": digest}
            for path, digest in sorted(PINNED_REFERENCES.items())
        ],
        "source": {
            "store": STORE_PATH,
            "leaf_count": STORE_LEAF_COUNT,
            "snapshot_sha256": STORE_SHA256,
            "member": MEMBER_KEY,
            "point_count": POINT_COUNT,
            "owned_row_count": OWNED_ROW_COUNT,
            "accepted_time_hex": ACCEPTED_TIME_HEX,
            "physical_state_sha256": PHYSICAL_STATE_SHA256,
            "descriptor_sha256": DESCRIPTOR_SHA256,
            "coordinates_sha256": COORDINATES_SHA256,
            "cursor_mode": "RETRY_PENDING",
            "pending_owner": "temporal",
            "spatial_operator": "Proto12GR0EvolutionOperator_SBP4",
            "measurement_tableau": "RK4",
            "rejection_records_are_predecessors": False,
        },
        "replays": list(replay_specs()),
        "future_origin": {
            "artifact_id": "FGC-1-PRO18-PREF27",
            "time": "23/16",
            "namespace": "runs/fgc-2-sf1/proto17/calibration",
            "generation": 0,
            "checkpoint_sha256": "7c059dcc2197a9f57baf8012f4113c240917090171ee06cc639e9e8734603d2d",
            "receipt_sha256": "78ae89a9f74193658003cf51f8567c4088db80fd6fb86867c2d581a0c16d5baf",
            "store_tree_sha256": "95a1bb96570216b8007e32a0cc05f4c31f420052f1d90c1b7134fa27ce958585",
            "used_as_measurement_state": False,
            "diagnostic_endpoint_adoption_allowed": False,
        },
        "selection": {
            "candidate_precedence": list(CANDIDATES),
            "channels": list(TDG6_COMPLETE_STATE_CHANNELS),
            "all_widths_all_channels_required": True,
            "minimum_observed_order": "3/2",
            "squared_order_multiplier": 8,
            "equality_passes": True,
            "arithmetic_order": [ORIGINAL_ARITHMETIC_ID, COMPENSATED_ARITHMETIC_ID],
            "roundoff_debit": "max_rows_sum_steps_abs_exact_accumulation_defect",
            "raw_baseline_classifications_retained": True,
            "separately_bound_PREF1_required": True,
            "one_shot": True,
            "no_pass_with_unresolved_candidate_is_scientific_failure": False,
        },
        "resources": dict(RESOURCES),
        "nonclaims": nonclaims(),
    }


def validate_config(value: object) -> dict[str, object]:
    if type(value) is not dict:
        raise MSEL1ContractError("config must be a dictionary")
    expected = expected_config(
        environment=value.get("environment"), implementation=value.get("implementation")
    )
    _strict_equal(value, expected)
    return deepcopy(expected)


def parse_config(raw: bytes) -> dict[str, object]:
    if type(raw) is not bytes or len(raw) > 1024 * 1024:
        raise MSEL1ContractError("config is not bounded immutable bytes")
    try:
        value = tomllib.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise MSEL1ContractError("config is not valid TOML") from error
    return validate_config(value)


def _order_values(A: tuple[tuple[Q, ...], ...], b: tuple[Q, ...]) -> tuple[Q, ...]:
    c = tuple(sum(row, Q(0)) for row in A)
    Ac = tuple(sum((a * x for a, x in zip(row, c, strict=True)), Q(0)) for row in A)
    Ac2 = tuple(sum((a * x**2 for a, x in zip(row, c, strict=True)), Q(0)) for row in A)
    AAc = tuple(sum((a * x for a, x in zip(row, Ac, strict=True)), Q(0)) for row in A)

    def dot(xs):
        return sum((a * x for a, x in zip(b, xs, strict=True)), Q(0))

    return (
        sum(b, Q(0)),
        dot(c),
        dot(tuple(x * x for x in c)),
        dot(Ac),
        dot(tuple(x**3 for x in c)),
        dot(tuple(x * y for x, y in zip(c, Ac, strict=True))),
        dot(Ac2),
        dot(AAc),
    )


def mathematical_controls() -> dict[str, object]:
    """Exact finite-algebra controls, never a sampled trajectory theorem."""

    A4 = (
        (Q(0),) * 5,
        (Q(1, 2), Q(0), Q(0), Q(0), Q(0)),
        (Q(0), Q(1, 2), Q(0), Q(0), Q(0)),
        (Q(0), Q(0), Q(1), Q(0), Q(0)),
        (Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6), Q(0)),
    )
    A3 = ((Q(0),) * 3, (Q(1), Q(0), Q(0)), (Q(1, 4), Q(1, 4), Q(0)))
    weights = (
        ("RK4", A4, A4[-1], 4),
        ("RK4_embedded", A4, (Q(1, 6), Q(1, 3), Q(1, 3), Q(0), Q(1, 6)), 3),
        ("SSPRK3", A3, (Q(1, 6), Q(1, 6), Q(2, 3)), 3),
        ("SSPRK3_embedded", A3, (Q(1, 2), Q(1, 2), Q(0)), 2),
    )
    expected = (Q(1), Q(1, 2), Q(1, 3), Q(1, 6), Q(1, 4), Q(1, 8), Q(1, 12), Q(1, 24))
    counts = {2: 2, 3: 4, 4: 8}
    records = []
    for name, A, b, order in weights:
        actual = _order_values(A, b)
        satisfied = tuple(x == y for x, y in zip(actual, expected, strict=True))
        if not all(satisfied[: counts[order]]) or (
            order < 4 and all(satisfied[: counts[order + 1]])
        ):
            raise AssertionError("embedded/main tableau order differs")
        records.append(
            {
                "name": name,
                "order": order,
                "order_values": [str(x) for x in actual],
                "satisfied": list(satisfied),
            }
        )
    left = (Q(1), Q(0), Q(-3), Q(2))
    right = (Q(0), Q(0), Q(3), Q(-2))
    if tuple(a + b for a, b in zip(left, right, strict=True)) != (1, 0, 0, 0):
        raise AssertionError("Hermite partition of unity differs")

    def bernstein(coefficients):
        a0, a1, a2, a3 = coefficients
        return a0, a0 + a1 / 3, a0 + 2 * a1 / 3 + a2 / 3, a0 + a1 + a2 + a3

    if bernstein(left) != (1, 1, 0, 0) or bernstein(right) != (0, 0, 1, 1):
        raise AssertionError("zero-slope Hermite convex-hull proof differs")
    equality = three_halves_order_passes_squared(
        outer_lower_squared=Q(8), finest_upper_squared=Q(1)
    )
    below = three_halves_order_passes_squared(
        outer_lower_squared=Q(8) - Q(1, 2**52), finest_upper_squared=Q(1)
    )
    if not equality or below:
        raise AssertionError("three-halves threshold changed")
    return {
        "tableaux": records,
        "zero_slope_Hermite_Bernstein_controls": ["e0", "e0", "e1", "e1"],
        "endpoint_error_bound": "max(abs(e0),abs(e1))<=sum_steps_abs_delta",
        "telescoping_identity": "y_j-z_j=sum_i_before_j_delta_i",
        "exact_threshold_equality_passes": equality,
        "threshold_below_equality_fails": not below,
        "no_PDE_stability_theorem_inferred": True,
    }


def compact_freeze(config_raw: bytes) -> dict[str, object]:
    config = parse_config(config_raw)
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "classification": CLASSIFICATION,
        "project_version": "0.11.0",
        "owner_document": OWNER_DOCUMENT,
        "config_sha256": sha256(config_raw).hexdigest(),
        "base_commit": BASE_COMMIT,
        "config": config,
        "mathematical_controls": mathematical_controls(),
        "diagnostic_executed": False,
        "selection_result_earned": False,
        "committed_image_authority_required": True,
        "nonclaims": nonclaims(),
    }


def validate_compact(config_raw: bytes, result: object) -> dict[str, object]:
    expected = compact_freeze(config_raw)
    _strict_equal(result, expected, "compact")
    canonical_json_bytes(expected)
    return expected


def render_config(value: object) -> bytes:
    """Deterministic TOML for this closed metadata schema; no file writes."""

    checked = validate_config(value)
    lines: list[str] = []

    def scalar(item):
        if type(item) not in (str, int, bool, list):
            raise MSEL1ContractError("unsupported TOML metadata value")
        if type(item) is list and any(
            type(part) not in (str, int, bool) for part in item
        ):
            raise MSEL1ContractError("unsupported TOML scalar array")
        return json.dumps(item, ensure_ascii=True, allow_nan=False)

    def table(data, prefix):
        for key, item in data.items():
            if type(item) is dict or (
                type(item) is list and item and type(item[0]) is dict
            ):
                continue
            lines.append(f"{key} = {scalar(item)}")
        for key, item in data.items():
            name = f"{prefix}.{key}" if prefix else key
            if type(item) is dict:
                lines.extend(("", f"[{name}]"))
                table(item, name)
            elif type(item) is list and item and type(item[0]) is dict:
                for record in item:
                    lines.extend(("", f"[[{name}]]"))
                    table(record, name)

    table(checked, "")
    return ("\n".join(lines) + "\n").encode("ascii")
