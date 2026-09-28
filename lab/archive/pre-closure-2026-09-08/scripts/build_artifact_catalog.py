#!/usr/bin/env python3
"""Build or verify the deterministic artifact catalog and compact result index."""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import sys
import tomllib
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
CATALOG_PATH = ROOT / "configs/fgc/artifact-catalog.json"
RESULT_INDEX_PATH = ROOT / "results/README.md"
SCHEMA = "FGC-artifact-catalog-v2"
PREF1_ARTIFACT = "FGC-1-TDG11-MSEL1-PREF1"
IMP1_ARTIFACT = "FGC-1-TDG11-IMP1"
REC1_PREF1_ARTIFACT = "FGC-1-HLT17-SRCQ1-REC1-PREF1"
REC1_FRZ1_ARTIFACT = "FGC-1-HLT17-SRCQ1-REC1-FRZ1"
PRO20_FRZ1_ARTIFACT = "FGC-1-PRO20-EV1-FRZ1"
PRO20_PREF1_ARTIFACT = "FGC-1-PRO20-EV1-PREF1"
DEF1_FRZ1_ARTIFACT = "FGC-1-DEF1-STAB1-FRZ1"
DEF1_STAB1_ARTIFACT = "FGC-1-DEF1-STAB1"
DEF1_PREF1_ARTIFACT = "FGC-1-DEF1-STAB1-PREF1"
DEF1_FRZ1_SOURCE = "src/recursive_horizons/fgc/def1_stab1_frz1_certificate.py"
DEF1_FRZ1_REPRODUCER = "scripts/reproduce_fgc_def1_stab1_frz1.py"
DEF1_FRZ1_OWNER = "docs/fgc-def1-stab1-frz1.md"
DEF1_FRZ1_CONFIG = "configs/fgc/fgc-1-def1-stab1-frz1.toml"
DEF1_FRZ1_RESULT = "results/fgc-1-def1-stab1-frz1.json"
DEF1_FRZ1_CLASSIFICATION = (
    "candidate_blind_conversion_error_map_instrument_freeze_no_trajectory"
)
DEF1_PREF1_SOURCE = "src/recursive_horizons/fgc/def1_stab1_pref1_binder.py"
DEF1_PREF1_REPRODUCER = "scripts/reproduce_fgc_def1_stab1_pref1.py"
DEF1_PREF1_OWNER = "docs/fgc-def1-stab1-pref1.md"
DEF1_PREF1_CONFIG = "configs/fgc/fgc-1-def1-stab1-pref1.toml"
DEF1_PREF1_RESULT = "results/fgc-1-def1-stab1-pref1.json"
DEF1_PREF1_CLASSIFICATION = (
    "independently_bound_candidate_blind_error_map_readiness_only_no_trajectory"
)
DEF1_PREF1_CONCLUSION = (
    "PREF1 binds candidate-blind pre-holdout error-map readiness only."
)
SOL1_FRZ1_ARTIFACT = "FGC-1-SGB1-CTL1-SOL1-FRZ1"
SOL1_ARTIFACT = "FGC-1-SGB1-CTL1-SOL1"
SOL1_FRZ1_SOURCE = "src/recursive_horizons/fgc/sgb1_ctl1_sol1_frz1_certificate.py"
SOL1_FRZ1_REPRODUCER = "scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py"
SOL1_FRZ1_OWNER = "docs/fgc-sgb1-ctl1-sol1-frz1.md"
SOL1_FRZ1_CONFIG = "configs/fgc/fgc-1-sgb1-ctl1-sol1-frz1.toml"
SOL1_FRZ1_RESULT = "results/fgc-1-sgb1-ctl1-sol1-frz1.json"
SOL1_FRZ1_CLASSIFICATION = (
    "prospective_sol1_frozen_nominal_interval_inconclusive_no_health"
)
SOL1_PREF1_ARTIFACT = "FGC-1-SGB1-CTL1-SOL1-PREF1"
SOL1_PREF1_SOURCE = "src/recursive_horizons/fgc/sgb1_ctl1_sol1_pref1_binder.py"
SOL1_PREF1_REPRODUCER = "scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py"
SOL1_PREF1_OWNER = "docs/fgc-sgb1-ctl1-sol1-pref1.md"
SOL1_PREF1_CONFIG = "configs/fgc/fgc-1-sgb1-ctl1-sol1-pref1.toml"
SOL1_PREF1_RESULT = "results/fgc-1-sgb1-ctl1-sol1-pref1.json"
SOL1_PREF1_CLASSIFICATION = (
    "independently_bound_sol1_nominal_interval_inconclusive_no_health"
)
CURRENT_FRONTIER = PRO20_PREF1_ARTIFACT
NEXT_DESIGN = "FGC-1-PRO20-EV1-RSRC1"
RSRC1_CLASSIFICATION = "prospective_resource_isolation_core_no_authority"
RSRC1_OWNER = "docs/fgc-pro20-ev1-rsrc1.md"
RSRC1_SOURCES = [
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py",
]
RSRC1_REPRODUCERS = ["scripts/run_fgc_pro20_rsrc1_child.py"]
RSRC1_TESTS = [
    "tests/test_fgc_pro20_rsrc1_attempt.py",
    "tests/test_fgc_pro20_rsrc1_isolation.py",
    "tests/test_fgc_pro20_rsrc1_member.py",
    "tests/test_fgc_pro20_rsrc1_seed.py",
]
QA2_PREF1 = "FGC-1-TDG10-QA2-PREF1"
MSEL1_FRZ1 = "FGC-1-TDG11-MSEL1-FRZ1"
PREF1_SOURCE_PREFIX = "tdg11_msel1_pref1_"
PREF1_SOURCE_ROOT = "src/recursive_horizons/fgc/evolution/"
PREF1_REPRODUCER = "scripts/reproduce_fgc_tdg11_msel1_pref1.py"
PREF1_OWNER = "docs/fgc-tdg11-msel1-pref1.md"
IMP1_SOURCE_PREFIX = "tdg11_imp1_"
IMP1_SOURCE_ROOT = "src/recursive_horizons/fgc/evolution/"
IMP1_REPRODUCER = "scripts/reproduce_fgc_tdg11_imp1.py"
IMP1_OWNER = "docs/fgc-tdg11-imp1.md"
PROSPECTIVE_INSTRUMENTS = {
    "FGC-1-SGB1-CTL1": {
        "classification": "partial_outcome_blind_sgb_source_and_synthetic_runtime_instruments",
        "conclusion_or_nonclaim": (
            "Exact annular constraints, highest-time affinity and pointwise source/Jacobian "
            "controls plus matched initial-family/center instruments are implemented; "
            "a branch-owned covariant/exact point-principal comparison is implemented; "
            "declared-box source invertibility and a branch-owned synthetic source/time-method "
            "runtime with rollback are implemented. Prospective parametric source admission, "
            "validated initial-constraint ODE, exact branch controls, continuum-cone reference "
            "and constraint-continuity instruments are also implemented. The nominal initial "
            "compactness enclosure remains nonpassing through two finite resource ladders, and "
            "the matched-family principal feeder retains whole-cell acceleration/cone owners; evolving "
            "center and aggregate branch health stay closed. No execution authority or complete "
            "SGB1 certificate."
        ),
        "owner_documents": [
            "docs/fgc-sgb1-ctl1.md",
            "docs/fgc-sgb1-ctl1-admission.md",
            "docs/fgc-sgb1-ctl1-initial-health.md",
            "docs/fgc-sgb1-ctl1-cone.md",
            "docs/fgc-sgb1-ctl1-continuity.md",
            "docs/fgc-sgb1-ctl1-trap-refinement.md",
            "docs/fgc-sgb1-ctl1-family-principal.md",
            "docs/fgc-sgb1-ctl1-cell-admission.md",
            "docs/fgc-sgb1-ctl1-local-symmetrizer.md",
        ],
        "source_paths": [
            "src/recursive_horizons/fgc/sgb1_ctl1_constraints.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_source.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_family.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_center.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_principal.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_interval_health.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_runtime.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_admission.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_initial_health.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_controls.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_cone.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_continuity.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_family_principal.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_cell_admission.py",
            "src/recursive_horizons/fgc/sgb1_ctl1_local_symmetrizer.py",
        ],
        "test_paths": [
            "tests/test_fgc_sgb1_ctl1_constraints.py",
            "tests/test_fgc_sgb1_ctl1_source.py",
            "tests/test_fgc_sgb1_ctl1_time_affinity.py",
            "tests/test_fgc_sgb1_ctl1_family.py",
            "tests/test_fgc_sgb1_ctl1_center.py",
            "tests/test_fgc_sgb1_ctl1_principal.py",
            "tests/test_fgc_sgb1_ctl1_interval_health.py",
            "tests/test_fgc_sgb1_ctl1_runtime.py",
            "tests/test_fgc_sgb1_ctl1_admission.py",
            "tests/test_fgc_sgb1_ctl1_initial_health.py",
            "tests/test_fgc_sgb1_ctl1_controls.py",
            "tests/test_fgc_sgb1_ctl1_cone.py",
            "tests/test_fgc_sgb1_ctl1_continuity.py",
            "tests/test_fgc_sgb1_ctl1_trap_refinement.py",
            "tests/test_fgc_sgb1_ctl1_family_principal.py",
            "tests/test_fgc_sgb1_ctl1_cell_admission.py",
            "tests/test_fgc_sgb1_ctl1_local_symmetrizer.py",
            "tests/test_fgc_wave0_instrument_routing.py",
        ],
        "predecessor_ids": ["FGC-1-ACT1", "FGC-1-VAR1", "FGC-1-HYP1-RED1", "FGC-1-ID1-FAM1"],
    },
    "FGC-1-SGB1-CTL1-TRAP-REFINEMENT2": {
        "classification": "finite_product_box_no_trap_refinement_nonpass_no_health",
        "conclusion_or_nonclaim": (
            "The separately predeclared depth-14/16/18 product-box ladder remains "
            "interval-inconclusive on gap-free nominal A_chi=3 prefixes; no level "
            "tiles the support with C<1 and no branch-health or holdout gate opens."
        ),
        "owner_documents": ["docs/fgc-sgb1-ctl1-trap-refinement2.md"],
        "source_paths": [
            "src/recursive_horizons/fgc/sgb1_ctl1_trap_refinement2.py",
        ],
        "test_paths": ["tests/test_fgc_sgb1_ctl1_trap_refinement2.py"],
        "predecessor_ids": ["FGC-1-SGB1-CTL1"],
    },
    "FGC-1-SGB1-CTL1-TRAP-TAYLOR": {
        "classification": "order4_taylor_no_trap_enclosure_nonpass_no_health",
        "conclusion_or_nonclaim": (
            "The separately frozen order-4 interval-jet Taylor route remains "
            "interval-inconclusive on nominal A_chi=3 gap-free prefixes; it never "
            "tiles the support with C<1 and opens no branch-health or holdout gate."
        ),
        "owner_documents": ["docs/fgc-sgb1-ctl1-trap-taylor.md"],
        "source_paths": ["src/recursive_horizons/fgc/sgb1_ctl1_trap_taylor.py"],
        "test_paths": ["tests/test_fgc_sgb1_ctl1_trap_taylor.py"],
        "predecessor_ids": ["FGC-1-SGB1-CTL1-TRAP-REFINEMENT2"],
    },
    "FGC-1-SGB1-CTL1-TRAP-BARRIER": {
        "classification": "barrier_nonpass_not_trajectory_trapping",
        "conclusion_or_nonclaim": (
            "An exact D=1-C crossing-surface witness has outward D_r<0, so the "
            "whole-domain Nagumo barrier proof fails for nominal A_chi=3. This does "
            "not show the actual orbit traps and does not reject SGB-L."
        ),
        "owner_documents": ["docs/fgc-sgb1-ctl1-trap-barrier.md"],
        "source_paths": ["src/recursive_horizons/fgc/sgb1_ctl1_trap_barrier.py"],
        "test_paths": ["tests/test_fgc_sgb1_ctl1_trap_barrier.py"],
        "predecessor_ids": ["FGC-1-SGB1-CTL1-TRAP-TAYLOR"],
    },
    "FGC-1-SGB1-CTL1-SOL1": {
        "classification": "orbit_local_picard_lindelof_interval_inconclusive_no_health",
        "conclusion_or_nonclaim": (
            "The prospectively frozen A_chi=3 orbit-local continuation is interval-"
            "inconclusive because the authenticated prefix endpoint is wider than the "
            "declared rho=1/8 tube. The policy is not retuned; this is not on-orbit "
            "trapping, SGB-L rejection, branch health, or holdout evidence."
        ),
        "owner_documents": ["docs/fgc-sgb1-ctl1-sol1.md"],
        "source_paths": ["src/recursive_horizons/fgc/sgb1_ctl1_sol1.py"],
        "test_paths": ["tests/test_fgc_sgb1_ctl1_sol1.py"],
        "predecessor_ids": ["FGC-1-SGB1-CTL1-TRAP-BARRIER"],
    },
    "FGC-1-DEF1-STAB1": {
        "classification": "partial_outcome_blind_error_map_and_margin_algebra",
        "conclusion_or_nonclaim": (
            "Checked error/margin contracts, mass-flux algebra and whole-box local Q "
            "sensitivities, joint-box provider conversions, exact conversion-ownership "
            "maps and a non-promoting freeze-contract payload are implemented. "
            "The separate STAB1-FRZ1 compact freeze binds this candidate-blind contract. "
            "Independent STAB1-PREF1 now binds candidate-blind pre-holdout map readiness "
            "only (DEF1_error_map_passed with map_readiness_only). Real trajectory "
            "input-enclosure values and the nine DEF1 booleans remain later; no "
            "trajectory or physical result."
        ),
        "owner_documents": ["docs/fgc-def1-stab1.md", "docs/fgc-def1-stab1-geometry.md"],
        "source_paths": [
            "src/recursive_horizons/fgc/def1_stab1.py",
            "src/recursive_horizons/fgc/def1_geometry_error.py",
            "src/recursive_horizons/fgc/def1_stab1_providers.py",
            "src/recursive_horizons/fgc/def1_stab1_qualification.py",
            "src/recursive_horizons/fgc/def1_stab1_freeze_contract.py",
        ],
        "test_paths": [
            "tests/test_fgc_def1_stab1.py",
            "tests/test_fgc_def1_geometry_error.py",
            "tests/test_fgc_def1_stab1_providers.py",
            "tests/test_fgc_def1_stab1_qualification.py",
            "tests/test_fgc_def1_stab1_freeze_contract.py",
            "tests/test_fgc_wave0_instrument_routing.py",
        ],
        "predecessor_ids": ["FGC-1-DEF0-OBS1", "FGC-1-HYP1-RED1"],
    },
    "FGC-1-TDG11-C1R1": {
        "classification": "partial_exact_c1_representation_no_production_authority",
        "conclusion_or_nonclaim": (
            "A separately named integer/exponent implementation reproduces the sealed IMP1 "
            "C1 mathematical/reference wire on bounded controls. No full-size physical-source "
            "qualification, production selection, trajectory, or physical result."
        ),
        "owner_documents": ["docs/fgc-tdg11-c1r1.md"],
        "source_paths": [
            "src/recursive_horizons/fgc/evolution/tdg11_c1r1_ring.py",
            "src/recursive_horizons/fgc/evolution/tdg11_c1r1_enclosure.py",
            "src/recursive_horizons/fgc/evolution/tdg11_c1r1_runtime.py",
        ],
        "test_paths": [
            "tests/test_fgc_tdg11_c1r1_ring.py",
            "tests/test_fgc_tdg11_c1r1_enclosure.py",
            "tests/test_fgc_tdg11_c1r1_runtime.py",
        ],
        "predecessor_ids": ["FGC-1-TDG11-IMP1"],
    },
    "FGC-1-HLT17-MON17": {
        "classification": "partial_c1r1_finite_state_runtime_and_atomic_store",
        "conclusion_or_nonclaim": (
            "Fixed C1R1 prepare/replay/commit, immutable member/cursor checkpoints, and a "
            "fail-closed atomic synthetic store are implemented. Physical origin/source, "
            "environment/resource qualification, binder, and campaign authority remain open."
        ),
        "owner_documents": ["docs/fgc-hlt17-mon17.md"],
        "source_paths": [
            "src/recursive_horizons/fgc/evolution/hlt17_admission_runtime.py",
            "src/recursive_horizons/fgc/evolution/hlt17_runtime_member.py",
            "src/recursive_horizons/fgc/evolution/hlt17_imp1_cursor.py",
            "src/recursive_horizons/fgc/evolution/hlt17_imp1_bridge.py",
            "src/recursive_horizons/fgc/evolution/hlt17_member_codec.py",
            "src/recursive_horizons/fgc/evolution/hlt17_member_checkpoint.py",
            "src/recursive_horizons/fgc/evolution/protocol_v19.py",
            "src/recursive_horizons/fgc/evolution/hlt17_campaign_store.py",
        ],
        "test_paths": [
            "tests/test_fgc_hlt17_runtime_member.py",
            "tests/test_fgc_hlt17_imp1_cursor.py",
            "tests/test_fgc_hlt17_imp1_bridge.py",
            "tests/test_fgc_hlt17_member_codec.py",
            "tests/test_fgc_hlt17_member_checkpoint.py",
            "tests/test_fgc_hlt17_c1r1_integration.py",
            "tests/test_fgc_protocol_v19.py",
            "tests/test_fgc_hlt17_campaign_store.py",
        ],
        "predecessor_ids": [
            "FGC-1-TDG11-C1R1", "FGC-1-TDG11-IMP1",
            "FGC-1-HLT16-MON16", "FGC-1-PRO18-PREF27",
        ],
    },
}
TRANSIENT_TEST_ID = "FGC-1-PRO19-SID3-REAL-STORE-PREFLIGHT"
TRANSIENT_TEST_PATH = (
    "archive/historical-tests/test_fgc_pro19_sid3_real_store_preflight.py"
)
TRANSIENT_TEST_SHA256 = (
    "09d3515947530384063ad83c11faaea25ef7d43d527f23b15b50597067e63b3e"
)
CATALOG_RELATIVE = "configs/fgc/artifact-catalog.json"
ALLOWED_STATUSES = (
    "compact_active",
    "current_frontier",
    "historical_authority",
    "historical_runner_disabled",
    "prospective_authority",
    "superseded_but_hash_bound",
    "transient_test_archived",
)
RECORD_KEYS = (
    "artifact_id",
    "classification",
    "compact_verify_target_or_none",
    "conclusion_or_nonclaim",
    "config_paths",
    "explicit_live_target_or_none",
    "family",
    "may_execute_state_change",
    "owner_documents",
    "predecessor_ids",
    "raw_dependencies",
    "reproducer_paths",
    "result_paths",
    "source_paths",
    "status",
    "successor_ids",
    "test_paths",
)
CLASSIFICATION_KEYS = (
    "classification",
    "heading_artifact_id",
    "path",
    "sha256",
)
ARTIFACT_ID = re.compile(r"(?:FGC-[12]-[A-Z0-9]+(?:-[A-Z0-9]+)*|FGC-RUNTIME-MATRIX-\d+)\Z")
MENTION_ID = re.compile(
    r"(?<![A-Za-z0-9-])"
    r"(FGC-[12]-[A-Z0-9]+(?:-[A-Z0-9]+)*|FGC-RUNTIME-MATRIX-\d+)"
    r"(?![A-Za-z0-9-])"
)
HEADING_ID = re.compile(
    r"(?m)^#\s+(?P<id>(?:FGC-[12]-[A-Z0-9]+(?:-[A-Z0-9]+)*"
    r"|GMF-[A-Z0-9]+(?:-[A-Z0-9]+)*"
    r"|EC-[A-Z0-9]+(?:-[A-Z0-9]+)*))"
)
FGC_ID_IN_TEXT = re.compile(r"FGC-[12]-[A-Z0-9][A-Z0-9-]*")
RAW_PATH = re.compile(r"runs/fgc-2-sf1/[A-Za-z0-9_./-]+")
DECLARED_ID = re.compile(
    r"(?m)^(?:ARTIFACT_ID|ARTIFACT)\s*(?::[^=]+)?=\s*[\"']([^\"']+)[\"']"
)
MAKEFILE_TARGET = re.compile(r"(?m)^([A-Za-z0-9_.-]+):")
MAKE_LIST_BODY = re.compile(
    r"(?P<name>[A-Z0-9_]+) := \\\n(?P<body>(?:\t[^\n]+(?: \\)?\n)+)"
)
STATE_PREFIXES = ("run-", "resume-", "recover-")
AUTHORITY_TOKENS = frozenset(
    {
        "AUTH",
        "AUTH1",
        "FRZ",
        "FRZ1",
        "GEN",
        "GEN1",
        "IMP",
        "IMP1",
        "IMP2",
        "IMP3",
        "MON",
        "MON1",
        "MON2",
        "MON3",
        "MON4",
        "MON5",
        "MON6",
        "MON7",
        "MON8",
        "MON9",
        "MON10",
        "MON11",
        "MON12",
        "MON13",
        "MON14",
        "MON15",
        "MON16",
        "PLAN",
    }
)
PREDECESSOR_TOKENS = (
    "ancestor",
    "authority",
    "base",
    "consumed",
    "immediate",
    "input",
    "parent",
    "predecessor",
    "prior",
    "source",
)
PROTO_ID = re.compile(r"\AFGC-2-SF1-PROTO(\d+)\Z")
CLOSED_COMPACT_TARGETS = {
    "run-fgc-pro13-calibration": "verify-fgc-pro13-result",
    "run-fgc-pro14-calibration": "verify-fgc-pro14-result",
    "recover-fgc-pro19-sid1": "fgc-pro19-sid2-pref1",
    "resume-fgc-pro19-sid3-event1": "fgc-pro19-pref28",
    "run-fgc-pro19-event1": "fgc-pro19-pref28",
    "resume-fgc-pro19-event1": "fgc-pro19-pref28",
    "run-fgc-tdg8-successor-event1": "verify-fgc-tdg8-rcv3-pref2",
    "recover-fgc-tdg8-rcv1": "verify-fgc-tdg8-rcv3-pref2",
    "resume-fgc-tdg8-rcv1": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv2": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv3": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg8-rcv3-rec1": "verify-fgc-tdg8-rcv3-pref2",
    "run-fgc-tdg9-ar1": "verify-fgc-tdg9-ar1-pref1",
    "run-fgc-tdg9-loc1": "verify-fgc-tdg9-loc2-pref2",
    "run-fgc-tdg9-loc2": "verify-fgc-tdg9-loc2-pref2",
    "run-fgc-tdg9-ti1": "verify-fgc-tdg9-ti2-pref1",
    "run-fgc-tdg9-ti2": "verify-fgc-tdg9-ti2-pref1",
    "run-fgc-tdg9-ac1": "verify-fgc-tdg9-ac1-pref1",
    "run-fgc-tdg9-ur1": "verify-fgc-tdg9-ur1-pref1",
    "run-fgc-tdg10-qa1": "verify-fgc-tdg10-qa1-pref1",
    "run-fgc-tdg10-qa2": "verify-fgc-tdg10-qa2-pref1",
    "run-fgc-tdg11-msel1": "fgc-tdg11-msel1-pref1",
    "run-fgc-hlt17-srcq1-rec1": "fgc-hlt17-srcq1-rec1-pref1",
    "run-fgc-pro20-ev1": "fgc-pro20-ev1-pref1",
}


def _canonical(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("ascii")


def _tracked_paths() -> tuple[str, ...]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    chunks = completed.stdout.split(b"\0")
    if chunks[-1] != b"":
        raise ValueError("Git path stream is not NUL terminated")
    return tuple(sorted(item.decode("utf-8") for item in chunks[:-1]))


def _text(path: str) -> str | None:
    try:
        raw = (ROOT / path).read_bytes()
    except OSError:
        return None
    if b"\0" in raw:
        return None
    return raw.decode("utf-8", "replace")


def _structured(path: str) -> Mapping[str, Any] | None:
    try:
        raw = (ROOT / path).read_bytes()
        if path.endswith(".json"):
            value = json.loads(raw)
        elif path.endswith(".toml"):
            value = tomllib.loads(raw.decode("utf-8"))
        else:
            return None
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, tomllib.TOMLDecodeError):
        return None
    return value if isinstance(value, Mapping) else None


def _top_artifact_id(path: str) -> str | None:
    value = _structured(path)
    artifact = value.get("artifact_id") if value is not None else None
    if isinstance(artifact, str) and ARTIFACT_ID.fullmatch(artifact):
        return artifact
    return None


def _heading_artifact_id(text: str) -> str | None:
    matched = HEADING_ID.search(text)
    if matched is None:
        return None
    return matched.group("id")


def _walk_strings(
    value: object, path: tuple[str, ...] = ()
) -> Iterable[tuple[tuple[str, ...], str]]:
    if isinstance(value, Mapping):
        for key, item in value.items():
            yield from _walk_strings(item, (*path, str(key)))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from _walk_strings(item, (*path, str(index)))
    elif isinstance(value, str):
        yield path, value


def _slug(artifact_id: str) -> str:
    value = artifact_id.lower()
    for prefix in ("fgc-1-", "fgc-2-sf1-"):
        if value.startswith(prefix):
            return value[len(prefix) :]
    return value


def _family(artifact_id: str) -> str:
    if artifact_id in CLOSED_TARGETS:
        return "closed_historical_operation"
    if artifact_id == TRANSIENT_TEST_ID:
        return "archived_transient_test"
    if artifact_id.startswith("FGC-1-"):
        return "FGC-1-" + artifact_id.split("-")[2]
    if artifact_id.startswith("FGC-2-"):
        return "FGC-2-" + artifact_id.split("-")[2]
    if artifact_id.startswith("FGC-"):
        parts = artifact_id.split("-")
        return "-".join(parts[:-1]) if len(parts) > 2 else artifact_id
    return artifact_id


def _parse_make_list(source: str, name: str) -> tuple[str, ...]:
    for matched in MAKE_LIST_BODY.finditer(source):
        if matched.group("name") != name:
            continue
        return tuple(
            line.strip().removesuffix(" \\")
            for line in matched.group("body").splitlines()
            if line.strip()
        )
    raise ValueError(f"Make list {name} is absent")


def _make_sources() -> str:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    includes = re.findall(r"(?m)^include\s+(\S+)\s*$", makefile)
    parts = [makefile]
    for relative in includes:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(f"Makefile include is absent: {relative}")
        parts.append(path.read_text(encoding="utf-8"))
    return "\n".join(parts)


def _make_targets(source: str, closed: tuple[str, ...]) -> tuple[str, ...]:
    names = set(MAKEFILE_TARGET.findall(source))
    names.update(closed)
    return tuple(sorted(names))


def _owner_from_config_path(path: str, tracked: set[str]) -> str | None:
    stem = Path(path).stem
    candidates = [f"docs/{stem}.md"]
    if stem.startswith("fgc-1-"):
        candidates.append(f"docs/fgc-{stem[len('fgc-1-'):]}.md")
    if stem.startswith("fgc-2-"):
        candidates.append(f"docs/fgc-{stem[len('fgc-2-'):]}.md")
    for candidate in candidates:
        if candidate in tracked:
            return candidate
    return None


def _classification_of(value: Mapping[str, Any] | None) -> str | None:
    if value is None:
        return None
    current = value.get("classification")
    if isinstance(current, str) and current:
        return current
    payload = value.get("artifact_payload")
    if isinstance(payload, Mapping):
        nested = payload.get("classification")
        if isinstance(nested, str) and nested:
            return nested
    return None


def _first_sentence(text: str) -> str:
    collapsed = " ".join(text.split())
    return collapsed if collapsed else text


def _summary(record: Mapping[str, Any] | None, classification: str | None) -> str:
    if record is not None:
        payload = record.get("artifact_payload")
        for owner in (payload, record):
            if not isinstance(owner, Mapping):
                continue
            nonclaims = owner.get("nonclaims")
            if isinstance(nonclaims, list) and nonclaims and isinstance(nonclaims[0], str):
                return _first_sentence(nonclaims[0])
            conclusion = owner.get("conclusion")
            if isinstance(conclusion, str) and conclusion:
                return _first_sentence(conclusion)
        gate = record.get("gate_status")
        if isinstance(gate, str) and gate:
            return f"gate_status={gate}; see the owner artifact for scope and nonclaims."
    if classification:
        return (
            f"classification={classification}; see the owner artifact for scope and nonclaims."
        )
    return "Tracked artifact; see its owner/configuration for the exact claim boundary."


def _compact_verify_target(artifact_id: str, targets: set[str]) -> str | None:
    slug = _slug(artifact_id)
    for name in (f"verify-fgc-{slug}", f"fgc-{slug}"):
        if name in targets:
            return name
    matches = [name for name in targets if name.startswith(f"verify-fgc-{slug}-")]
    if not matches:
        matches = [
            name
            for name in targets
            if name.startswith(f"fgc-{slug}-")
            and not name.startswith(STATE_PREFIXES)
        ]
    if not matches:
        return None

    def rank(name: str) -> tuple[int, str]:
        if "pref" in name or name.endswith("-result"):
            return (0, name)
        if "compact" in name:
            return (1, name)
        if "runtime" in name or "post" in name:
            return (2, name)
        return (3, name)

    return sorted(matches, key=rank)[0]


def _authority_like(artifact_id: str) -> bool:
    tokens = artifact_id.upper().split("-")
    if tokens[-1].endswith("PLAN") or "PLAN" in tokens:
        return True
    return any(
        token in AUTHORITY_TOKENS or token.startswith(("AUTH", "FRZ", "MON", "IMP", "GEN"))
        for token in tokens
    )


def _proto_number(artifact_id: str) -> int | None:
    matched = PROTO_ID.fullmatch(artifact_id)
    return int(matched.group(1)) if matched else None


def _status(
    artifact_id: str,
    *,
    current_dev: set[str],
    latest_proto: int | None,
) -> str:
    if artifact_id == CURRENT_FRONTIER:
        return "current_frontier"
    if artifact_id == NEXT_DESIGN:
        return "prospective_authority"
    if artifact_id in CLOSED_TARGETS:
        return "historical_runner_disabled"
    if artifact_id == TRANSIENT_TEST_ID:
        return "transient_test_archived"
    if _authority_like(artifact_id):
        return "historical_authority"
    if artifact_id in current_dev:
        return "compact_active"
    proto = _proto_number(artifact_id)
    if proto is not None and latest_proto is not None:
        return "compact_active" if proto == latest_proto else "superseded_but_hash_bound"
    return "superseded_but_hash_bound"


def _blank_record(artifact_id: str) -> dict[str, Any]:
    return {
        "artifact_id": artifact_id,
        "classification": None,
        "compact_verify_target_or_none": None,
        "conclusion_or_nonclaim": (
            "Tracked artifact; see its owner/configuration for the exact claim boundary."
        ),
        "config_paths": [],
        "explicit_live_target_or_none": None,
        "family": _family(artifact_id),
        "may_execute_state_change": False,
        "owner_documents": [],
        "predecessor_ids": [],
        "raw_dependencies": [],
        "reproducer_paths": [],
        "result_paths": [],
        "source_paths": [],
        "status": "compact_active",
        "successor_ids": [],
        "test_paths": [],
    }


def _safe_raw_path(value: str) -> str | None:
    cleaned = value.rstrip(".,;:)")
    path = Path(cleaned)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        return None
    if not cleaned.startswith("runs/"):
        return None
    return cleaned


def _classify_unindexed(path: str, heading_id: str | None) -> str:
    if path.startswith("configs/fgc/"):
        return "auxiliary_bundle_without_artifact_id"
    if path.startswith("results/gmf-") or path.startswith("results/einstein-"):
        return "non_fgc_campaign_compact_result"
    if path.startswith("results/"):
        return "programme_compact_result"
    if heading_id and heading_id.startswith(("GMF-", "EC-")):
        return "non_fgc_campaign_owner_document"
    return "explicitly_classified_path"


def _explicit_classification(path: str, heading_id: str | None) -> dict[str, Any]:
    return {
        "classification": _classify_unindexed(path, heading_id),
        "heading_artifact_id": heading_id,
        "path": path,
        "sha256": sha256((ROOT / path).read_bytes()).hexdigest(),
    }


def _record_keys(record: dict[str, Any]) -> dict[str, Any]:
    return {key: record[key] for key in RECORD_KEYS}


CLOSED_TARGETS: tuple[str, ...] = ()


def build_catalog() -> dict[str, Any]:
    global CLOSED_TARGETS
    tracked = _tracked_paths()
    tracked_set = set(tracked)
    make_source = _make_sources()
    CLOSED_TARGETS = _parse_make_list(make_source, "CLOSED_HISTORICAL_STATE_TARGETS")
    current_dev_stems = _parse_make_list(make_source, "CURRENT_DEVELOPMENT_CERTIFICATES")
    targets = set(_make_targets(make_source, CLOSED_TARGETS))
    if set(CLOSED_TARGETS) != set(CLOSED_COMPACT_TARGETS):
        raise ValueError("closed-operation compact-verifier map differs")
    missing_compact_targets = sorted(set(CLOSED_COMPACT_TARGETS.values()) - targets)
    if missing_compact_targets:
        raise ValueError(
            f"closed-operation compact verifier is absent: {missing_compact_targets}"
        )

    relevant = tuple(
        path
        for path in tracked
        if path.startswith(
            (
                "configs/fgc/",
                "docs/",
                "results/",
                "scripts/",
                "src/recursive_horizons/fgc/evolution/",
                "tests/",
            )
        )
        and path != CATALOG_RELATIVE
        and path != "results/README.md"
    )
    texts = {path: text for path in relevant if (text := _text(path)) is not None}
    structured_cache = {
        path: _structured(path)
        for path in relevant
        if path.endswith((".json", ".toml"))
    }

    declarations: dict[str, list[str]] = defaultdict(list)
    heading_docs: dict[str, list[str]] = defaultdict(list)
    unindexed_headings: dict[str, str] = {}
    for path, text in texts.items():
        if path.startswith("docs/") and path.endswith(".md"):
            heading = _heading_artifact_id(text)
            if heading is None:
                continue
            if ARTIFACT_ID.fullmatch(heading):
                heading_docs[heading].append(path)
            else:
                unindexed_headings[path] = heading
        declared = _top_artifact_id(path)
        if declared is None and path.endswith(".py"):
            matched = DECLARED_ID.search(text)
            if matched and ARTIFACT_ID.fullmatch(matched.group(1)):
                declared = matched.group(1)
        if declared is not None:
            declarations[declared].append(path)

    artifact_ids = set(declarations)
    artifact_ids.update(heading_docs)
    artifact_ids.update(
        {CURRENT_FRONTIER, NEXT_DESIGN, TRANSIENT_TEST_ID, DEF1_PREF1_ARTIFACT}
    )
    artifact_ids.update(CLOSED_TARGETS)

    current_dev_ids = set()
    for stem in current_dev_stems:
        candidate = "FGC-1-" + stem.removeprefix("fgc-").upper()
        if candidate in artifact_ids:
            current_dev_ids.add(candidate)

    latest_proto = max(
        (number for artifact_id in artifact_ids if (number := _proto_number(artifact_id)) is not None),
        default=None,
    )

    records: dict[str, dict[str, Any]] = {
        artifact_id: _blank_record(artifact_id) for artifact_id in artifact_ids
    }

    for path in tracked:
        if path == CATALOG_RELATIVE:
            continue
        if path.startswith("configs/fgc/") and path.endswith((".toml", ".json")):
            artifact_id = _top_artifact_id(path)
            if artifact_id is not None:
                records[artifact_id]["config_paths"].append(path)
        elif path.startswith("results/") and path.endswith(".json"):
            artifact_id = _top_artifact_id(path)
            if artifact_id is not None:
                records[artifact_id]["result_paths"].append(
                    {
                        "path": path,
                        "sha256": sha256((ROOT / path).read_bytes()).hexdigest(),
                    }
                )

    for artifact_id, paths in heading_docs.items():
        records[artifact_id]["owner_documents"].extend(paths)

    for artifact_id, paths in declarations.items():
        record = records[artifact_id]
        for path in paths:
            if path.startswith("src/recursive_horizons/fgc/evolution/"):
                record["source_paths"].append(path)
            elif path.startswith("scripts/"):
                record["reproducer_paths"].append(path)
            elif path.startswith("tests/"):
                record["test_paths"].append(path)
            elif path.startswith("docs/") and path.endswith(".md"):
                record["owner_documents"].append(path)

    for artifact_id, record in records.items():
        for path in record["config_paths"]:
            owner = _owner_from_config_path(path, tracked_set)
            if owner is not None:
                record["owner_documents"].append(owner)
            value = structured_cache.get(path)
            if not isinstance(value, Mapping):
                continue
            documented = value.get("owner_document")
            if isinstance(documented, str) and documented in tracked_set:
                record["owner_documents"].append(documented)
            derived = value.get("derivation_document")
            if isinstance(derived, str) and derived in tracked_set:
                record["owner_documents"].append(derived)
            if record["classification"] is None:
                record["classification"] = _classification_of(value)
        for result in record["result_paths"]:
            value = structured_cache.get(result["path"])
            if not isinstance(value, Mapping):
                continue
            if record["classification"] is None:
                record["classification"] = _classification_of(value)
            derived = value.get("derivation_document")
            if isinstance(derived, str) and derived in tracked_set:
                record["owner_documents"].append(derived)
            payload = value.get("artifact_payload")
            if isinstance(payload, Mapping):
                derived = payload.get("derivation_document")
                if isinstance(derived, str) and derived in tracked_set:
                    record["owner_documents"].append(derived)

    for artifact_id, source_root, source_prefix, reproducer, owner, test_prefix, repo_test in (
        (
            PREF1_ARTIFACT,
            PREF1_SOURCE_ROOT,
            PREF1_SOURCE_PREFIX,
            PREF1_REPRODUCER,
            PREF1_OWNER,
            "test_fgc_tdg11_msel1_pref1_",
            "test_check_repo_tdg11_msel1_pref1.py",
        ),
        (
            IMP1_ARTIFACT,
            IMP1_SOURCE_ROOT,
            IMP1_SOURCE_PREFIX,
            IMP1_REPRODUCER,
            IMP1_OWNER,
            "test_fgc_tdg11_imp1_",
            "test_check_repo_tdg11_imp1.py",
        ),
    ):
        record = records[artifact_id]
        for path in tracked:
            name = Path(path).name
            if (
                path.startswith(source_root)
                and name.startswith(source_prefix)
                and path.endswith(".py")
            ):
                record["source_paths"].append(path)
            elif path == reproducer:
                record["reproducer_paths"].append(path)
            elif path.startswith("tests/") and (
                name.startswith(test_prefix) or name == repo_test
            ):
                record["test_paths"].append(path)
        if owner in tracked_set:
            record["owner_documents"].append(owner)

    mentions_by_id: dict[str, set[str]] = defaultdict(set)
    for path, text in texts.items():
        for match in MENTION_ID.findall(text):
            if match in records:
                mentions_by_id[match].add(path)

    for artifact_id, record in records.items():
        associated = set(
            record["config_paths"]
            + [item["path"] for item in record["result_paths"]]
            + record["owner_documents"]
            + record["source_paths"]
            + record["reproducer_paths"]
            + record["test_paths"]
        )
        mentioned = mentions_by_id.get(artifact_id, set()) - associated
        record["source_paths"].extend(
            path
            for path in mentioned
            if path.startswith("src/recursive_horizons/fgc/evolution/")
        )
        record["reproducer_paths"].extend(
            path
            for path in mentioned
            if path.startswith("scripts/reproduce_")
        )
        record["test_paths"].extend(
            path for path in mentioned if path.startswith("tests/")
        )
        associated = sorted(
            associated
            | set(record["source_paths"])
            | set(record["reproducer_paths"])
            | set(record["test_paths"])
        )
        raw_dependencies = {
            cleaned
            for path in associated
            for match in RAW_PATH.findall(texts.get(path, ""))
            if (cleaned := _safe_raw_path(match)) is not None
        }
        record["raw_dependencies"] = sorted(raw_dependencies)
        predecessors: set[str] = set()
        for path in record["config_paths"] + [item["path"] for item in record["result_paths"]]:
            value = structured_cache.get(path)
            if value is None:
                continue
            for key_path, string in _walk_strings(value):
                joined = ".".join(key_path).lower()
                if not any(token in joined for token in PREDECESSOR_TOKENS):
                    continue
                for candidate in FGC_ID_IN_TEXT.findall(string):
                    if candidate != artifact_id and candidate in records:
                        predecessors.add(candidate)
        record["predecessor_ids"] = sorted(predecessors)
        primary_result = None
        if record["result_paths"]:
            primary_result = structured_cache.get(record["result_paths"][0]["path"])
        record["conclusion_or_nonclaim"] = _summary(
            primary_result if isinstance(primary_result, Mapping) else None,
            record["classification"],
        )
        record["compact_verify_target_or_none"] = _compact_verify_target(
            artifact_id, targets
        )
        record["status"] = _status(
            artifact_id,
            current_dev=current_dev_ids,
            latest_proto=latest_proto,
        )

    records[NEXT_DESIGN]["classification"] = RSRC1_CLASSIFICATION
    records[NEXT_DESIGN]["conclusion_or_nonclaim"] = (
        "Prospective no-store per-member seed and attempt isolation is implemented. "
        "Production store, namespace, authority and live execution remain absent. Any new "
        "measurement needs a separately committed RSRC1-FRZ1."
    )
    records[NEXT_DESIGN]["predecessor_ids"] = [CURRENT_FRONTIER]
    records[NEXT_DESIGN]["status"] = "prospective_authority"
    records[NEXT_DESIGN]["may_execute_state_change"] = False
    records[NEXT_DESIGN]["explicit_live_target_or_none"] = None
    records[NEXT_DESIGN]["compact_verify_target_or_none"] = None
    records[NEXT_DESIGN]["config_paths"] = []
    records[NEXT_DESIGN]["result_paths"] = []
    records[NEXT_DESIGN]["reproducer_paths"] = list(RSRC1_REPRODUCERS)
    records[NEXT_DESIGN]["source_paths"] = list(RSRC1_SOURCES)
    records[NEXT_DESIGN]["test_paths"] = list(RSRC1_TESTS)
    records[NEXT_DESIGN]["owner_documents"] = [RSRC1_OWNER]
    records[NEXT_DESIGN]["raw_dependencies"] = []

    for artifact_id, specification in PROSPECTIVE_INSTRUMENTS.items():
        if artifact_id not in records:
            raise ValueError(f"prospective instrument owner is absent: {artifact_id}")
        record = records[artifact_id]
        if record["config_paths"] or record["result_paths"]:
            raise ValueError(f"{artifact_id} needs explicit catalog promotion before a certificate")
        record.update(specification)
        record["status"] = "prospective_authority"
        record["reproducer_paths"] = []
        record["raw_dependencies"] = []
        record["may_execute_state_change"] = False
        record["explicit_live_target_or_none"] = None
        record["compact_verify_target_or_none"] = None

    if SOL1_FRZ1_ARTIFACT not in records:
        raise ValueError("SOL1-FRZ1 compact freeze owner is absent")
    records[SOL1_FRZ1_ARTIFACT]["status"] = "compact_active"
    records[SOL1_FRZ1_ARTIFACT]["classification"] = SOL1_FRZ1_CLASSIFICATION
    records[SOL1_FRZ1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[SOL1_FRZ1_ARTIFACT]["predecessor_ids"]) | {SOL1_ARTIFACT}
    )
    records[SOL1_FRZ1_ARTIFACT]["source_paths"] = [SOL1_FRZ1_SOURCE]
    records[SOL1_FRZ1_ARTIFACT]["reproducer_paths"] = [SOL1_FRZ1_REPRODUCER]
    records[SOL1_FRZ1_ARTIFACT]["owner_documents"] = [SOL1_FRZ1_OWNER]
    records[SOL1_FRZ1_ARTIFACT]["test_paths"] = [
        "tests/test_fgc_sgb1_ctl1_sol1_frz1_certificate.py",
        "tests/test_check_repo_sgb1_ctl1_sol1_frz1.py",
    ]
    records[SOL1_FRZ1_ARTIFACT]["config_paths"] = [SOL1_FRZ1_CONFIG]
    records[SOL1_FRZ1_ARTIFACT]["raw_dependencies"] = []
    records[SOL1_FRZ1_ARTIFACT]["may_execute_state_change"] = False
    records[SOL1_FRZ1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[SOL1_FRZ1_ARTIFACT]["compact_verify_target_or_none"] = (
        "fgc-sgb1-ctl1-sol1-frz1"
    )
    records[SOL1_ARTIFACT]["successor_ids"] = sorted(
        set(records[SOL1_ARTIFACT]["successor_ids"]) | {SOL1_FRZ1_ARTIFACT}
    )
    if SOL1_PREF1_ARTIFACT not in records:
        raise ValueError("SOL1-PREF1 compact binder owner is absent")
    records[SOL1_FRZ1_ARTIFACT]["successor_ids"] = sorted(
        set(records[SOL1_FRZ1_ARTIFACT]["successor_ids"]) | {SOL1_PREF1_ARTIFACT}
    )
    records[SOL1_PREF1_ARTIFACT]["status"] = "compact_active"
    records[SOL1_PREF1_ARTIFACT]["classification"] = SOL1_PREF1_CLASSIFICATION
    records[SOL1_PREF1_ARTIFACT]["predecessor_ids"] = [SOL1_FRZ1_ARTIFACT]
    records[SOL1_PREF1_ARTIFACT]["source_paths"] = [SOL1_PREF1_SOURCE]
    records[SOL1_PREF1_ARTIFACT]["reproducer_paths"] = [SOL1_PREF1_REPRODUCER]
    records[SOL1_PREF1_ARTIFACT]["owner_documents"] = [SOL1_PREF1_OWNER]
    records[SOL1_PREF1_ARTIFACT]["test_paths"] = [
        "tests/test_fgc_sgb1_ctl1_sol1_pref1_binder.py",
        "tests/test_check_repo_sgb1_ctl1_sol1_pref1.py",
    ]
    records[SOL1_PREF1_ARTIFACT]["config_paths"] = [SOL1_PREF1_CONFIG]
    records[SOL1_PREF1_ARTIFACT]["raw_dependencies"] = []
    records[SOL1_PREF1_ARTIFACT]["may_execute_state_change"] = False
    records[SOL1_PREF1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[SOL1_PREF1_ARTIFACT]["compact_verify_target_or_none"] = (
        "fgc-sgb1-ctl1-sol1-pref1"
    )

    if DEF1_FRZ1_ARTIFACT not in records:
        raise ValueError("DEF1-STAB1-FRZ1 compact freeze owner is absent")
    records[DEF1_FRZ1_ARTIFACT]["status"] = "compact_active"
    records[DEF1_FRZ1_ARTIFACT]["classification"] = DEF1_FRZ1_CLASSIFICATION
    records[DEF1_FRZ1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["predecessor_ids"]) | {DEF1_STAB1_ARTIFACT}
    )
    records[DEF1_FRZ1_ARTIFACT]["source_paths"] = [DEF1_FRZ1_SOURCE]
    records[DEF1_FRZ1_ARTIFACT]["reproducer_paths"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["reproducer_paths"]) | {DEF1_FRZ1_REPRODUCER}
    )
    records[DEF1_FRZ1_ARTIFACT]["owner_documents"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["owner_documents"]) | {DEF1_FRZ1_OWNER}
    )
    records[DEF1_FRZ1_ARTIFACT]["test_paths"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["test_paths"])
        | {
            "tests/test_fgc_def1_stab1_frz1_certificate.py",
            "tests/test_check_repo_def1_stab1_frz1.py",
        }
    )
    records[DEF1_FRZ1_ARTIFACT]["config_paths"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["config_paths"]) | {DEF1_FRZ1_CONFIG}
    )
    records[DEF1_FRZ1_ARTIFACT]["raw_dependencies"] = []
    records[DEF1_FRZ1_ARTIFACT]["may_execute_state_change"] = False
    records[DEF1_FRZ1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[DEF1_FRZ1_ARTIFACT]["compact_verify_target_or_none"] = "fgc-def1-stab1-frz1"
    records[DEF1_FRZ1_ARTIFACT]["successor_ids"] = sorted(
        set(records[DEF1_FRZ1_ARTIFACT]["successor_ids"]) | {DEF1_PREF1_ARTIFACT}
    )

    if DEF1_PREF1_ARTIFACT not in records:
        records[DEF1_PREF1_ARTIFACT] = _blank_record(DEF1_PREF1_ARTIFACT)
    pref1_result_hash = sha256((ROOT / DEF1_PREF1_RESULT).read_bytes()).hexdigest()
    records[DEF1_PREF1_ARTIFACT]["status"] = "compact_active"
    records[DEF1_PREF1_ARTIFACT]["classification"] = DEF1_PREF1_CLASSIFICATION
    records[DEF1_PREF1_ARTIFACT]["conclusion_or_nonclaim"] = DEF1_PREF1_CONCLUSION
    records[DEF1_PREF1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[DEF1_PREF1_ARTIFACT]["predecessor_ids"]) | {DEF1_FRZ1_ARTIFACT}
    )
    records[DEF1_PREF1_ARTIFACT]["source_paths"] = [DEF1_PREF1_SOURCE]
    records[DEF1_PREF1_ARTIFACT]["reproducer_paths"] = sorted(
        set(records[DEF1_PREF1_ARTIFACT]["reproducer_paths"]) | {DEF1_PREF1_REPRODUCER}
    )
    records[DEF1_PREF1_ARTIFACT]["owner_documents"] = sorted(
        set(records[DEF1_PREF1_ARTIFACT]["owner_documents"]) | {DEF1_PREF1_OWNER}
    )
    records[DEF1_PREF1_ARTIFACT]["test_paths"] = sorted(
        set(records[DEF1_PREF1_ARTIFACT]["test_paths"])
        | {
            "tests/test_fgc_def1_stab1_pref1_binder.py",
            "tests/test_check_repo_def1_stab1_pref1.py",
        }
    )
    records[DEF1_PREF1_ARTIFACT]["config_paths"] = sorted(
        set(records[DEF1_PREF1_ARTIFACT]["config_paths"]) | {DEF1_PREF1_CONFIG}
    )
    records[DEF1_PREF1_ARTIFACT]["result_paths"] = [
        {"path": DEF1_PREF1_RESULT, "sha256": pref1_result_hash}
    ]
    records[DEF1_PREF1_ARTIFACT]["raw_dependencies"] = []
    records[DEF1_PREF1_ARTIFACT]["may_execute_state_change"] = False
    records[DEF1_PREF1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[DEF1_PREF1_ARTIFACT]["compact_verify_target_or_none"] = "fgc-def1-stab1-pref1"

    records[PREF1_ARTIFACT]["status"] = "compact_active"
    records[PREF1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[PREF1_ARTIFACT]["predecessor_ids"]) | {MSEL1_FRZ1}
    )
    records[PREF1_ARTIFACT]["may_execute_state_change"] = False
    records[PREF1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[PREF1_ARTIFACT]["compact_verify_target_or_none"] = "fgc-tdg11-msel1-pref1"
    records[PREF1_ARTIFACT]["successor_ids"] = sorted(
        set(records[PREF1_ARTIFACT]["successor_ids"]) | {IMP1_ARTIFACT}
    )
    records[IMP1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[IMP1_ARTIFACT]["predecessor_ids"]) | {PREF1_ARTIFACT}
    )

    records[REC1_PREF1_ARTIFACT]["status"] = "compact_active"
    records[REC1_PREF1_ARTIFACT]["conclusion_or_nonclaim"] = (
        "All six GR-0 origin members are independently bound through physical "
        "source/C1R1 qualification without state advance. This licenses only a "
        "separately frozen PRO20 first-event authority, not event execution, "
        "calibration, a candidate, mechanism, or physical claim."
    )
    records[REC1_PREF1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[REC1_PREF1_ARTIFACT]["predecessor_ids"])
        | {IMP1_ARTIFACT, REC1_FRZ1_ARTIFACT}
    )
    records[REC1_PREF1_ARTIFACT]["may_execute_state_change"] = False
    records[REC1_PREF1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[REC1_PREF1_ARTIFACT]["compact_verify_target_or_none"] = (
        "fgc-hlt17-srcq1-rec1-pref1"
    )
    records[REC1_PREF1_ARTIFACT]["test_paths"] = sorted(
        set(records[REC1_PREF1_ARTIFACT]["test_paths"])
        | {
            "tests/test_fgc_hlt17_srcq1_rec1_pref1_binder.py",
            "tests/test_check_repo_hlt17_srcq1_rec1_pref1.py",
        }
    )
    records[IMP1_ARTIFACT]["successor_ids"] = sorted(
        set(records[IMP1_ARTIFACT]["successor_ids"]) | {REC1_PREF1_ARTIFACT}
    )
    records[REC1_FRZ1_ARTIFACT]["successor_ids"] = sorted(
        set(records[REC1_FRZ1_ARTIFACT]["successor_ids"]) | {REC1_PREF1_ARTIFACT}
    )

    records[PRO20_FRZ1_ARTIFACT]["status"] = "historical_authority"
    records[PRO20_FRZ1_ARTIFACT]["classification"] = (
        "consumed_one_shot_first_event_authority"
    )
    records[PRO20_FRZ1_ARTIFACT]["conclusion_or_nonclaim"] = (
        "The corrected one-shot PRO20 first-event freeze was consumed. The independent "
        "PREF1 certificate binds a typed resource_exhausted terminal with partial "
        "accepted state and no common event. Resume, takeover, calibration and physics "
        "remain closed."
    )
    records[PRO20_FRZ1_ARTIFACT]["may_execute_state_change"] = False
    records[PRO20_FRZ1_ARTIFACT]["explicit_live_target_or_none"] = None
    records[PRO20_FRZ1_ARTIFACT]["compact_verify_target_or_none"] = "fgc-pro20-ev1-frz1"
    records[PRO20_FRZ1_ARTIFACT]["predecessor_ids"] = sorted(
        set(records[PRO20_FRZ1_ARTIFACT]["predecessor_ids"]) | {REC1_PREF1_ARTIFACT}
    )

    records[CURRENT_FRONTIER]["status"] = "current_frontier"
    records[CURRENT_FRONTIER]["classification"] = (
        "independently_bound_resource_exhausted_no_common_event"
    )
    records[CURRENT_FRONTIER]["conclusion_or_nonclaim"] = (
        "Independent certificate of the one authorized PRO20 first event. The event "
        "created the production campaign store and persisted 146 accepted fine states. "
        "The binder and certificate wrote nothing and serialized or adopted no "
        "additional endpoint. Typed resource_exhausted at SSPRK3-8193 because max RSS "
        "exceeded, with explicit partial accepted state and no common event. Not a "
        "scientific nonpass, calibration, PRO21/Wave2, SGB/DEF, holdout, candidate, "
        "mechanism or physics claim."
    )
    records[CURRENT_FRONTIER]["predecessor_ids"] = sorted(
        set(records[CURRENT_FRONTIER]["predecessor_ids"]) | {PRO20_FRZ1_ARTIFACT}
    )
    records[CURRENT_FRONTIER]["successor_ids"] = sorted(
        set(records[CURRENT_FRONTIER]["successor_ids"]) | {NEXT_DESIGN}
    )
    records[CURRENT_FRONTIER]["may_execute_state_change"] = False
    records[CURRENT_FRONTIER]["explicit_live_target_or_none"] = None
    records[CURRENT_FRONTIER]["compact_verify_target_or_none"] = "fgc-pro20-ev1-pref1"
    records[CURRENT_FRONTIER]["config_paths"] = sorted(
        set(records[CURRENT_FRONTIER]["config_paths"])
        | {"configs/fgc/fgc-1-pro20-ev1-pref1.toml"}
    )
    records[CURRENT_FRONTIER]["owner_documents"] = sorted(
        set(records[CURRENT_FRONTIER]["owner_documents"])
        | {"docs/fgc-pro20-ev1-pref1.md"}
    )
    records[CURRENT_FRONTIER]["reproducer_paths"] = sorted(
        set(records[CURRENT_FRONTIER]["reproducer_paths"])
        | {"scripts/reproduce_fgc_pro20_ev1_pref1.py"}
    )
    records[CURRENT_FRONTIER]["source_paths"] = sorted(
        set(records[CURRENT_FRONTIER]["source_paths"])
        | {
            "src/recursive_horizons/fgc/evolution/pro20_ev1_pref1_certificate.py",
            "src/recursive_horizons/fgc/evolution/pro20_ev1_pref1_binder.py",
        }
    )
    pref1_result = "results/fgc-1-pro20-ev1-pref1.json"
    records[CURRENT_FRONTIER]["result_paths"] = sorted(
        records[CURRENT_FRONTIER]["result_paths"]
        + (
            []
            if any(item["path"] == pref1_result for item in records[CURRENT_FRONTIER]["result_paths"])
            else [{
                "path": pref1_result,
                "sha256": "3ec879fdfd309864b913b036f6b8c08c4717cb40269a18749cc19b8c0ca7b09b",
            }]
        ),
        key=lambda item: item["path"],
    )
    records[CURRENT_FRONTIER]["test_paths"] = sorted(
        set(records[CURRENT_FRONTIER]["test_paths"])
        | {
            "tests/test_fgc_pro20_ev1_pref1_certificate.py",
            "tests/test_check_repo_pro20_ev1_pref1.py",
        }
    )
    records[PRO20_FRZ1_ARTIFACT]["successor_ids"] = sorted(
        set(records[PRO20_FRZ1_ARTIFACT]["successor_ids"]) | {CURRENT_FRONTIER}
    )

    if QA2_PREF1 in records:
        records[QA2_PREF1]["conclusion_or_nonclaim"] = (
            "Retry 4 and retry 5 each fail u:alpha, u:lambda, and u:R; this rejects only "
            "the tested SSPRK3-on-inherited-SBP4 plus exact-C remedy and authorizes no state advance."
        )
        records[QA2_PREF1]["may_execute_state_change"] = False
        records[QA2_PREF1]["explicit_live_target_or_none"] = None

    for target in CLOSED_TARGETS:
        record = records[target]
        record["classification"] = "closed_historical_state_changing_make_target"
        record["conclusion_or_nonclaim"] = (
            "Historical operation completed or terminal; state-changing rerun is not "
            "authorized; use the named compact verifier."
        )
        record["status"] = "historical_runner_disabled"
        record["may_execute_state_change"] = False
        record["explicit_live_target_or_none"] = None
        record["compact_verify_target_or_none"] = CLOSED_COMPACT_TARGETS[target]
        record["family"] = "closed_historical_operation"

    test_record = records[TRANSIENT_TEST_ID]
    test_record["classification"] = "archived_transient_live_fixture_preflight"
    test_record["conclusion_or_nonclaim"] = (
        "SID3 real-store preflight archived from tests/ after the live generation-eight "
        "fixture ceased to be the store tip; SHA-256 "
        f"{TRANSIENT_TEST_SHA256} is preserved; the file is not discovered by the "
        "canonical test runner."
    )
    test_record["status"] = "transient_test_archived"
    test_record["may_execute_state_change"] = False
    test_record["explicit_live_target_or_none"] = None
    test_record["compact_verify_target_or_none"] = None
    test_record["test_paths"] = [TRANSIENT_TEST_PATH]
    test_record["owner_documents"] = ["archive/historical-tests/README.md"]
    if "FGC-1-PRO19-SID3-AUTH1" in records:
        test_record["predecessor_ids"] = ["FGC-1-PRO19-SID3-AUTH1"]
    test_record["family"] = "archived_transient_test"

    for record in records.values():
        record["config_paths"] = sorted(set(record["config_paths"]))
        record["owner_documents"] = sorted(set(record["owner_documents"]))
        record["source_paths"] = sorted(set(record["source_paths"]))
        record["reproducer_paths"] = sorted(set(record["reproducer_paths"]))
        record["test_paths"] = sorted(set(record["test_paths"]))
        record["raw_dependencies"] = sorted(set(record["raw_dependencies"]))
        record["predecessor_ids"] = sorted(
            {
                item
                for item in record["predecessor_ids"]
                if item in records and item != record["artifact_id"]
            }
        )
        record["result_paths"] = sorted(record["result_paths"], key=lambda item: item["path"])
        record["successor_ids"] = []
        record["may_execute_state_change"] = False

    for artifact_id, record in records.items():
        for predecessor in record["predecessor_ids"]:
            records[predecessor]["successor_ids"] = sorted(
                set(records[predecessor]["successor_ids"]) | {artifact_id}
            )

    records[NEXT_DESIGN]["may_execute_state_change"] = False
    records[NEXT_DESIGN]["explicit_live_target_or_none"] = None
    records[NEXT_DESIGN]["compact_verify_target_or_none"] = None
    records[CURRENT_FRONTIER]["may_execute_state_change"] = False
    records[CURRENT_FRONTIER]["explicit_live_target_or_none"] = None
    records[CURRENT_FRONTIER]["successor_ids"] = sorted(
        set(records[CURRENT_FRONTIER]["successor_ids"]) | {NEXT_DESIGN}
    )
    records[PRO20_FRZ1_ARTIFACT]["may_execute_state_change"] = False
    records[PRO20_FRZ1_ARTIFACT]["explicit_live_target_or_none"] = None
    if MSEL1_FRZ1 in records:
        records[MSEL1_FRZ1]["may_execute_state_change"] = False
        records[MSEL1_FRZ1]["explicit_live_target_or_none"] = None
    if DEF1_FRZ1_ARTIFACT in records:
        records[DEF1_FRZ1_ARTIFACT]["status"] = "compact_active"
        records[DEF1_FRZ1_ARTIFACT]["classification"] = DEF1_FRZ1_CLASSIFICATION
        records[DEF1_FRZ1_ARTIFACT]["source_paths"] = [DEF1_FRZ1_SOURCE]
        records[DEF1_FRZ1_ARTIFACT]["raw_dependencies"] = []
        records[DEF1_FRZ1_ARTIFACT]["may_execute_state_change"] = False
        records[DEF1_FRZ1_ARTIFACT]["explicit_live_target_or_none"] = None
        records[DEF1_FRZ1_ARTIFACT]["compact_verify_target_or_none"] = (
            "fgc-def1-stab1-frz1"
        )
        records[DEF1_FRZ1_ARTIFACT]["successor_ids"] = sorted(
            set(records[DEF1_FRZ1_ARTIFACT]["successor_ids"]) | {DEF1_PREF1_ARTIFACT}
        )
    if DEF1_PREF1_ARTIFACT in records:
        records[DEF1_PREF1_ARTIFACT]["status"] = "compact_active"
        records[DEF1_PREF1_ARTIFACT]["classification"] = DEF1_PREF1_CLASSIFICATION
        records[DEF1_PREF1_ARTIFACT]["conclusion_or_nonclaim"] = DEF1_PREF1_CONCLUSION
        records[DEF1_PREF1_ARTIFACT]["source_paths"] = [DEF1_PREF1_SOURCE]
        records[DEF1_PREF1_ARTIFACT]["raw_dependencies"] = []
        records[DEF1_PREF1_ARTIFACT]["may_execute_state_change"] = False
        records[DEF1_PREF1_ARTIFACT]["explicit_live_target_or_none"] = None
        records[DEF1_PREF1_ARTIFACT]["compact_verify_target_or_none"] = (
            "fgc-def1-stab1-pref1"
        )

    explicit: list[dict[str, Any]] = []
    classified_paths: set[str] = set()
    for path in tracked:
        if path == CATALOG_RELATIVE or path == "results/README.md":
            continue
        if path.startswith("configs/fgc/") and path.endswith((".toml", ".json")):
            if _top_artifact_id(path) is None:
                explicit.append(_explicit_classification(path, None))
                classified_paths.add(path)
        elif path.startswith("results/") and path.endswith(".json"):
            if _top_artifact_id(path) is None:
                explicit.append(_explicit_classification(path, None))
                classified_paths.add(path)
    for path, heading in unindexed_headings.items():
        if path in classified_paths:
            continue
        explicit.append(_explicit_classification(path, heading))
    explicit.sort(key=lambda item: item["path"])

    return {
        "schema": SCHEMA,
        "allowed_statuses": list(ALLOWED_STATUSES),
        "frontier": {
            "current_artifact_id": CURRENT_FRONTIER,
            "next_design_artifact_id": NEXT_DESIGN,
            "state_advance_authorized": False,
        },
        "explicit_classifications": explicit,
        "artifacts": [_record_keys(records[key]) for key in sorted(records)],
    }


def build_result_index(catalog: Mapping[str, Any]) -> bytes:
    lines = [
        "# Result index",
        "",
        "Generated deterministically by `scripts/build_artifact_catalog.py` from tracked "
        "compact results and `configs/fgc/artifact-catalog.json`.",
        "",
        "Ordinary compact verification does not infer scientific promotion from this index.",
        "",
        "| Artifact ID | Class | SHA-256 | Owner | Conclusion / nonclaim |",
        "|---|---|---|---|---|",
    ]
    for artifact in catalog["artifacts"]:
        owners = artifact["owner_documents"]
        owner = "—"
        if owners:
            relative = "../" + owners[0]
            owner = f"[{Path(owners[0]).name}]({relative})"
        classification = artifact["classification"] or "—"
        conclusion = str(artifact["conclusion_or_nonclaim"]).replace("|", "\\|")
        for result in artifact["result_paths"]:
            result_link = Path(result["path"]).name
            lines.append(
                f"| `{artifact['artifact_id']}` | `{classification}` | "
                f"[`{result['sha256']}`]({result_link}) | {owner} | {conclusion} |"
            )
    for classified in catalog["explicit_classifications"]:
        path = str(classified["path"])
        if not path.startswith("results/") or not path.endswith(".json"):
            continue
        synthetic_id = "UNINDEXED-" + re.sub(
            r"[^A-Z0-9]+", "-", Path(path).stem.upper()
        ).strip("-")
        lines.append(
            f"| `{synthetic_id}` | `{classified['classification']}` | "
            f"[`{classified['sha256']}`]({Path(path).name}) | — | "
            "Explicitly classified legacy/programme result; see the result file "
            "and claim ledger for scope. |"
        )
    lines.extend(
        [
            "",
            f"Current frontier: `{catalog['frontier']['current_artifact_id']}`. "
            "Next design: "
            f"`{catalog['frontier']['next_design_artifact_id']}`. PREF1 independently "
            "binds a typed resource_exhausted stop: the event created the production "
            "campaign store and persisted 146 accepted fine states; the binder and "
            "certificate wrote nothing. Partial accepted state is explicit and no "
            "common event completed. RSRC1 now has prospective no-store per-member "
            "seed/attempt isolation, but no production store, namespace, authority or live run.",
            "",
        ]
    )
    return "\n".join(lines).encode("utf-8")


def _check_or_write(*, write: bool) -> int:
    catalog = build_catalog()
    catalog_raw = _canonical(catalog)
    index_raw = build_result_index(catalog)
    failures = []
    for path, expected in ((CATALOG_PATH, catalog_raw), (RESULT_INDEX_PATH, index_raw)):
        if write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(expected)
        elif not path.is_file():
            failures.append(f"missing generated artifact: {path.relative_to(ROOT)}")
        elif path.read_bytes() != expected:
            failures.append(f"generated artifact differs: {path.relative_to(ROOT)}")
    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    mode = "wrote" if write else "verified"
    print(
        f"{mode} artifact catalog: {len(catalog['artifacts'])} entries; "
        f"frontier={CURRENT_FRONTIER}; next={NEXT_DESIGN}"
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    return _check_or_write(write=args.write)


if __name__ == "__main__":
    raise SystemExit(main())
