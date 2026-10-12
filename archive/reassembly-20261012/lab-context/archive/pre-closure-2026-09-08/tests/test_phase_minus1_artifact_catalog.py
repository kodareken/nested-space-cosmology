"""Phase -1 artifact-catalog schema, coverage, and authority-routing invariants."""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import re
import unittest
from unittest.mock import patch

from scripts.build_artifact_catalog import (
    ALLOWED_STATUSES,
    CLOSED_COMPACT_TARGETS,
    CURRENT_FRONTIER,
    IMP1_ARTIFACT,
    PRO20_FRZ1_ARTIFACT,
    PRO20_PREF1_ARTIFACT,
    REC1_PREF1_ARTIFACT,
    IMP1_OWNER,
    IMP1_REPRODUCER,
    MSEL1_FRZ1,
    NEXT_DESIGN,
    RSRC1_CLASSIFICATION,
    RSRC1_REPRODUCERS,
    RSRC1_SOURCES,
    RSRC1_TESTS,
    QA2_PREF1,
    PREF1_ARTIFACT,
    PREF1_OWNER,
    PREF1_REPRODUCER,
    PROSPECTIVE_INSTRUMENTS,
    RECORD_KEYS,
    SCHEMA,
    TRANSIENT_TEST_ID,
    TRANSIENT_TEST_PATH,
    TRANSIENT_TEST_SHA256,
    build_catalog,
    build_result_index,
)
from scripts.repo_checks import catalog as catalog_checks
from scripts.repo_checks.catalog import check_artifact_catalog


ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "configs/fgc/artifact-catalog.json"
RESULT_INDEX = ROOT / "results/README.md"
CLOSED = ROOT / "mk/closed-live-targets.mk"
README = ROOT / "README.md"
ACTIVE_MAP = ROOT / "docs/active-code-map.md"
TDG8_OWNER = "docs/fgc-tdg8-frz1.md"
PRO20_OWNER = "docs/fgc-pro20-ev1-frz1.md"
PRO20_CONFIG = "configs/fgc/fgc-1-pro20-ev1-frz1.toml"
PRO20_SOURCES = {
    "src/recursive_horizons/fgc/evolution/pro20_ev1_authority.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_protocol.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_runtime.py",
    "src/recursive_horizons/fgc/evolution/pro20_ev1_store.py",
}
PRO20_TESTS = {
    "tests/test_fgc_pro20_ev1_authority.py",
    "tests/test_fgc_pro20_ev1_protocol.py",
    "tests/test_fgc_pro20_ev1_runtime.py",
    "tests/test_fgc_pro20_ev1_store.py",
}
PRO20_REPRODUCERS = {
    "scripts/reproduce_fgc_pro20_ev1_frz1.py",
    "scripts/run_fgc_pro20_ev1.py",
}

EXPECTED_CLOSED = (
    "run-fgc-pro13-calibration",
    "run-fgc-pro14-calibration",
    "recover-fgc-pro19-sid1",
    "resume-fgc-pro19-sid3-event1",
    "run-fgc-pro19-event1",
    "resume-fgc-pro19-event1",
    "run-fgc-tdg8-successor-event1",
    "recover-fgc-tdg8-rcv1",
    "resume-fgc-tdg8-rcv1",
    "run-fgc-tdg8-rcv2",
    "run-fgc-tdg8-rcv3",
    "run-fgc-tdg8-rcv3-rec1",
    "run-fgc-tdg9-ar1",
    "run-fgc-tdg9-loc1",
    "run-fgc-tdg9-loc2",
    "run-fgc-tdg9-ti1",
    "run-fgc-tdg9-ti2",
    "run-fgc-tdg9-ac1",
    "run-fgc-tdg9-ur1",
    "run-fgc-tdg10-qa1",
    "run-fgc-tdg10-qa2",
    "run-fgc-tdg11-msel1",
    "run-fgc-hlt17-srcq1-rec1",
    "run-fgc-pro20-ev1",
)
PREF1_SOURCES = (
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_binder.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_protocol.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_reconstruction.py",
    "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_localization.py",
)
PREF1_TESTS = (
    "tests/test_fgc_tdg11_msel1_pref1_binder.py",
    "tests/test_fgc_tdg11_msel1_pref1_protocol.py",
    "tests/test_fgc_tdg11_msel1_pref1_reconstruction.py",
    "tests/test_fgc_tdg11_msel1_pref1_localization.py",
    "tests/test_check_repo_tdg11_msel1_pref1.py",
)
IMP1_SOURCES = (
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_enclosure.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_ledger.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_runtime.py",
    "src/recursive_horizons/fgc/evolution/tdg11_imp1_qualification.py",
)
IMP1_TESTS = (
    "tests/test_fgc_tdg11_imp1_enclosure.py",
    "tests/test_fgc_tdg11_imp1_ledger.py",
    "tests/test_fgc_tdg11_imp1_runtime.py",
    "tests/test_fgc_tdg11_imp1_retry_limits.py",
    "tests/test_fgc_tdg11_imp1_reproduction.py",
    "tests/test_check_repo_tdg11_imp1.py",
)
IMP1_RESULT_SHA256 = "1c3b57d239fcbb62cce0d196b3c943dd7a94a34f613fdc8f05c7736c63e7edc1"
SCIENTIFIC_AUTHORITY_MARKERS = (
    "`AGENTS.md` — creative research orientation.",
    "[Claim ledger]",
    "[Research roadmap]",
    "runtime matrix",
    "Results index",
    "Versioned paper",
    "short project entrypoint and routing page",
)


def _closed_targets(source: str) -> tuple[str, ...]:
    match = re.search(
        r"CLOSED_HISTORICAL_STATE_TARGETS := \\\n(?P<body>(?:\t[^\n]+(?: \\)?\n)+)",
        source,
    )
    if match is None:
        raise AssertionError("closed target variable is absent")
    return tuple(
        line.strip().removesuffix(" \\")
        for line in match.group("body").splitlines()
        if line.strip()
    )


class PhaseMinusOneArtifactCatalogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.catalog = json.loads(CATALOG.read_text(encoding="ascii"))
        cls.by_id = {entry["artifact_id"]: entry for entry in cls.catalog["artifacts"]}
        cls.counts = Counter(entry["status"] for entry in cls.catalog["artifacts"])

    def test_repo_check_accepts_the_generated_catalog(self) -> None:
        self.assertEqual(check_artifact_catalog(), [])

    def test_schema_is_strict_and_sorted(self) -> None:
        self.assertEqual(self.catalog["schema"], SCHEMA)
        self.assertEqual(self.catalog["schema"], "FGC-artifact-catalog-v2")
        self.assertEqual(
            set(self.catalog["allowed_statuses"]),
            set(ALLOWED_STATUSES),
        )
        self.assertEqual(
            self.catalog["allowed_statuses"],
            sorted(ALLOWED_STATUSES),
        )
        self.assertEqual(
            set(ALLOWED_STATUSES),
            {
                "current_frontier",
                "compact_active",
                "prospective_authority",
                "historical_authority",
                "historical_runner_disabled",
                "superseded_but_hash_bound",
                "transient_test_archived",
            },
        )
        for entry in self.catalog["artifacts"]:
            self.assertEqual(set(entry), set(RECORD_KEYS))
            self.assertIn(entry["status"], ALLOWED_STATUSES)
            self.assertIsInstance(entry["family"], str)
            self.assertTrue(entry["family"])
            self.assertIsInstance(entry["conclusion_or_nonclaim"], str)
            self.assertTrue(entry["conclusion_or_nonclaim"])
            self.assertIsInstance(entry["may_execute_state_change"], bool)
            self.assertEqual(entry["owner_documents"], sorted(set(entry["owner_documents"])))
            self.assertEqual(entry["config_paths"], sorted(set(entry["config_paths"])))
            self.assertEqual(
                [item["path"] for item in entry["result_paths"]],
                sorted(item["path"] for item in entry["result_paths"]),
            )
            self.assertEqual(entry["predecessor_ids"], sorted(set(entry["predecessor_ids"])))
            self.assertEqual(entry["successor_ids"], sorted(set(entry["successor_ids"])))

    def test_every_status_is_machine_distinguished(self) -> None:
        for status in ALLOWED_STATUSES:
            self.assertGreater(self.counts[status], 0, status)
        self.assertEqual(self.counts["current_frontier"], 1)
        self.assertEqual(
            self.counts["prospective_authority"],
            1 + len(PROSPECTIVE_INSTRUMENTS),
        )
        self.assertEqual(self.counts["historical_runner_disabled"], 24)
        self.assertEqual(self.counts["transient_test_archived"], 1)

    def test_frontier_successor_and_only_prospective_first_event_routing(self) -> None:
        frontier = self.catalog["frontier"]
        current = self.by_id[CURRENT_FRONTIER]
        successor = self.by_id[NEXT_DESIGN]
        self.assertEqual(PREF1_ARTIFACT, "FGC-1-TDG11-MSEL1-PREF1")
        self.assertEqual(CURRENT_FRONTIER, "FGC-1-PRO20-EV1-PREF1")
        self.assertEqual(PRO20_PREF1_ARTIFACT, CURRENT_FRONTIER)
        self.assertEqual(NEXT_DESIGN, "FGC-1-PRO20-EV1-RSRC1")
        self.assertEqual(
            frontier["current_artifact_id"],
            "FGC-1-PRO20-EV1-PREF1",
        )
        self.assertEqual(frontier["next_design_artifact_id"], "FGC-1-PRO20-EV1-RSRC1")
        self.assertIs(frontier["state_advance_authorized"], False)
        self.assertEqual(current["status"], "current_frontier")
        self.assertEqual(successor["status"], "prospective_authority")
        self.assertEqual(
            current["classification"],
            "independently_bound_resource_exhausted_no_common_event",
        )
        self.assertEqual(
            successor["classification"],
            RSRC1_CLASSIFICATION,
        )
        self.assertIn(
            "no-store per-member seed and attempt",
            successor["conclusion_or_nonclaim"].lower(),
        )
        self.assertIn("RSRC1-FRZ1", successor["conclusion_or_nonclaim"])
        self.assertNotIn("numerical method", successor["conclusion_or_nonclaim"])
        self.assertIn(NEXT_DESIGN, current["successor_ids"])
        self.assertEqual(successor["predecessor_ids"], [CURRENT_FRONTIER])
        self.assertIs(current["may_execute_state_change"], False)
        self.assertIs(successor["may_execute_state_change"], False)
        self.assertIsNone(successor["explicit_live_target_or_none"])
        self.assertIsNone(successor["compact_verify_target_or_none"])
        self.assertEqual(successor["config_paths"], [])
        self.assertEqual(successor["result_paths"], [])
        self.assertEqual(set(successor["source_paths"]), set(RSRC1_SOURCES))
        self.assertEqual(set(successor["reproducer_paths"]), set(RSRC1_REPRODUCERS))
        self.assertEqual(set(successor["test_paths"]), set(RSRC1_TESTS))
        self.assertEqual(successor["raw_dependencies"], [])
        self.assertIn("docs/fgc-pro20-ev1-rsrc1.md", successor["owner_documents"])
        self.assertIsNone(current["explicit_live_target_or_none"])
        self.assertEqual(
            current["compact_verify_target_or_none"],
            "fgc-pro20-ev1-pref1",
        )
        self.assertEqual(
            [key for key, entry in self.by_id.items() if entry["may_execute_state_change"]],
            [],
        )
        self.assertEqual(
            [
                key
                for key, entry in self.by_id.items()
                if entry["explicit_live_target_or_none"] is not None
            ],
            [],
        )
        frz1 = self.by_id[PRO20_FRZ1_ARTIFACT]
        self.assertEqual(frz1["status"], "historical_authority")
        self.assertIs(frz1["may_execute_state_change"], False)
        self.assertIsNone(frz1["explicit_live_target_or_none"])
        self.assertEqual(frz1["compact_verify_target_or_none"], "fgc-pro20-ev1-frz1")
        self.assertEqual(frz1["config_paths"], [PRO20_CONFIG])
        self.assertTrue(PRO20_SOURCES.issubset(set(frz1["source_paths"])))
        self.assertTrue(PRO20_TESTS.issubset(set(frz1["test_paths"])))
        self.assertEqual(set(frz1["reproducer_paths"]), PRO20_REPRODUCERS)
        self.assertIn(PRO20_OWNER, frz1["owner_documents"])
        rec1 = self.by_id[REC1_PREF1_ARTIFACT]
        self.assertEqual(rec1["status"], "compact_active")
        self.assertEqual(
            rec1["classification"],
            "independently_bound_all_six_source_c1r1_qualified_no_state_advance",
        )
        qa2 = self.by_id[QA2_PREF1]
        self.assertNotEqual(qa2["status"], "current_frontier")
        self.assertIn("exact-C remedy", qa2["conclusion_or_nonclaim"])
        frozen = self.by_id[MSEL1_FRZ1]
        self.assertEqual(frozen["status"], "historical_authority")
        self.assertIs(frozen["may_execute_state_change"], False)
        self.assertIsNone(frozen["explicit_live_target_or_none"])
        self.assertEqual(frozen["compact_verify_target_or_none"], "fgc-tdg11-msel1-frz1")
        pref1 = self.by_id[PREF1_ARTIFACT]
        sol1_frz1 = self.by_id["FGC-1-SGB1-CTL1-SOL1-FRZ1"]
        self.assertEqual(sol1_frz1["status"], "compact_active")
        self.assertEqual(
            sol1_frz1["classification"],
            "prospective_sol1_frozen_nominal_interval_inconclusive_no_health",
        )
        self.assertEqual(
            sol1_frz1["compact_verify_target_or_none"],
            "fgc-sgb1-ctl1-sol1-frz1",
        )
        self.assertIs(sol1_frz1["may_execute_state_change"], False)
        self.assertIsNone(sol1_frz1["explicit_live_target_or_none"])
        self.assertEqual(sol1_frz1["raw_dependencies"], [])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1", sol1_frz1["predecessor_ids"])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1-PREF1", sol1_frz1["successor_ids"])
        sol1_pref1 = self.by_id["FGC-1-SGB1-CTL1-SOL1-PREF1"]
        self.assertEqual(sol1_pref1["status"], "compact_active")
        self.assertEqual(
            sol1_pref1["classification"],
            "independently_bound_sol1_nominal_interval_inconclusive_no_health",
        )
        self.assertEqual(
            sol1_pref1["compact_verify_target_or_none"],
            "fgc-sgb1-ctl1-sol1-pref1",
        )
        self.assertIs(sol1_pref1["may_execute_state_change"], False)
        self.assertIsNone(sol1_pref1["explicit_live_target_or_none"])
        self.assertEqual(sol1_pref1["raw_dependencies"], [])
        self.assertIn("FGC-1-SGB1-CTL1-SOL1-FRZ1", sol1_pref1["predecessor_ids"])
        def1_frz1 = self.by_id["FGC-1-DEF1-STAB1-FRZ1"]
        self.assertEqual(def1_frz1["status"], "compact_active")
        self.assertEqual(
            def1_frz1["classification"],
            "candidate_blind_conversion_error_map_instrument_freeze_no_trajectory",
        )
        self.assertEqual(
            def1_frz1["compact_verify_target_or_none"],
            "fgc-def1-stab1-frz1",
        )
        self.assertIs(def1_frz1["may_execute_state_change"], False)
        self.assertIsNone(def1_frz1["explicit_live_target_or_none"])
        self.assertEqual(def1_frz1["raw_dependencies"], [])
        self.assertIn("FGC-1-DEF1-STAB1", def1_frz1["predecessor_ids"])
        self.assertIn("FGC-1-DEF1-STAB1-PREF1", def1_frz1["successor_ids"])
        self.assertIn("FGC-1-DEF1-STAB1-FRZ1", self.by_id["FGC-1-DEF1-STAB1"]["successor_ids"])
        def1_pref1 = self.by_id["FGC-1-DEF1-STAB1-PREF1"]
        self.assertEqual(def1_pref1["status"], "compact_active")
        self.assertEqual(
            def1_pref1["classification"],
            "independently_bound_candidate_blind_error_map_readiness_only_no_trajectory",
        )
        self.assertEqual(
            def1_pref1["compact_verify_target_or_none"],
            "fgc-def1-stab1-pref1",
        )
        self.assertIs(def1_pref1["may_execute_state_change"], False)
        self.assertIsNone(def1_pref1["explicit_live_target_or_none"])
        self.assertEqual(def1_pref1["raw_dependencies"], [])
        self.assertIn("FGC-1-DEF1-STAB1-FRZ1", def1_pref1["predecessor_ids"])
        self.assertIn(PREF1_ARTIFACT, frozen["successor_ids"])
        self.assertIn(MSEL1_FRZ1, pref1["predecessor_ids"])
        self.assertEqual(pref1["status"], "compact_active")
        self.assertEqual(pref1["compact_verify_target_or_none"], "fgc-tdg11-msel1-pref1")
        self.assertIn(IMP1_ARTIFACT, pref1["successor_ids"])
        closed_run = self.by_id["run-fgc-tdg11-msel1"]
        self.assertEqual(closed_run["status"], "historical_runner_disabled")
        self.assertEqual(closed_run["compact_verify_target_or_none"], "fgc-tdg11-msel1-pref1")
        self.assertIs(closed_run["may_execute_state_change"], False)
        pro20_run = self.by_id["run-fgc-pro20-ev1"]
        self.assertEqual(pro20_run["status"], "historical_runner_disabled")
        self.assertEqual(pro20_run["compact_verify_target_or_none"], "fgc-pro20-ev1-pref1")
        self.assertIs(pro20_run["may_execute_state_change"], False)

    def test_existing_runners_do_not_make_history_current(self) -> None:
        for artifact_id, entry in self.by_id.items():
            if artifact_id == CURRENT_FRONTIER:
                continue
            self.assertNotEqual(entry["status"], "current_frontier")
            if any(
                Path(path).name.startswith("run_") for path in entry["reproducer_paths"]
            ) or any("run_fgc_" in path for path in entry["source_paths"]):
                self.assertNotEqual(entry["status"], "current_frontier")
                self.assertIs(entry["may_execute_state_change"], False)

    def test_pref1_and_imp1_keep_separate_owner_module_and_test_associations(self) -> None:
        for artifact_id, sources, tests, owner, reproducer in (
            (PREF1_ARTIFACT, PREF1_SOURCES, PREF1_TESTS, PREF1_OWNER, PREF1_REPRODUCER),
            (IMP1_ARTIFACT, IMP1_SOURCES, IMP1_TESTS, IMP1_OWNER, IMP1_REPRODUCER),
        ):
            with self.subTest(artifact_id=artifact_id):
                entry = self.by_id[artifact_id]
                self.assertEqual(entry["source_paths"], sorted(sources))
                for path in tests:
                    self.assertIn(path, entry["test_paths"])
                self.assertIn(owner, entry["owner_documents"])
                self.assertIn(reproducer, entry["reproducer_paths"])
        current = self.by_id[CURRENT_FRONTIER]
        self.assertEqual(
            current["config_paths"],
            ["configs/fgc/fgc-1-pro20-ev1-pref1.toml"],
        )
        self.assertEqual(
            current["result_paths"],
            [{
                "path": "results/fgc-1-pro20-ev1-pref1.json",
                "sha256": "3ec879fdfd309864b913b036f6b8c08c4717cb40269a18749cc19b8c0ca7b09b",
            }],
        )
        self.assertNotIn(PREF1_OWNER, current["owner_documents"])
        self.assertNotIn(PREF1_REPRODUCER, current["reproducer_paths"])

    def _check_mutated_catalog(self, catalog: dict[str, object]) -> list[str]:
        original_read = catalog_checks._read_unique_regular_bytes

        def read(relative: str, label: str) -> bytes:
            if relative == "configs/fgc/artifact-catalog.json":
                return json.dumps(catalog).encode("ascii")
            return original_read(relative, label)

        with patch.object(catalog_checks, "_read_unique_regular_bytes", side_effect=read):
            return check_artifact_catalog()

    def test_checker_rejects_missing_pref1_and_imp1_associations(self) -> None:
        for artifact_id, label, source, test in (
            (PREF1_ARTIFACT, "PREF1", PREF1_SOURCES[0], PREF1_TESTS[0]),
            (IMP1_ARTIFACT, "IMP1", IMP1_SOURCES[0], IMP1_TESTS[0]),
        ):
            for key, path in (("source_paths", source), ("test_paths", test)):
                with self.subTest(artifact_id=artifact_id, key=key):
                    catalog = deepcopy(self.catalog)
                    entry = next(
                        item for item in catalog["artifacts"]
                        if item["artifact_id"] == artifact_id
                    )
                    entry[key].remove(path)
                    self.assertIn(
                        f"artifact catalog {label} source/config/result/owner/test association differs",
                        self._check_mutated_catalog(catalog),
                    )

    def test_checker_rejects_absent_pref1_instead_of_skipping_it(self) -> None:
        catalog = deepcopy(self.catalog)
        catalog["artifacts"] = [
            entry for entry in catalog["artifacts"]
            if entry["artifact_id"] != PREF1_ARTIFACT
        ]
        failures = self._check_mutated_catalog(catalog)
        self.assertIn("artifact catalog compact-active PREF1 routing differs", failures)
        self.assertIn(
            "artifact catalog PREF1 source/config/result/owner/test association differs",
            failures,
        )

    def test_checker_rejects_focused_suite_as_imp1_compact_target(self) -> None:
        catalog = deepcopy(self.catalog)
        current = next(
            entry for entry in catalog["artifacts"]
            if entry["artifact_id"] == CURRENT_FRONTIER
        )
        current["compact_verify_target_or_none"] = (
            "verify-fgc-pro20-ev1-pref1"
        )
        self.assertIn(
            "artifact catalog current-to-next edge differs",
            self._check_mutated_catalog(catalog),
        )

    def test_exact_disabled_operations_match_closed_targets(self) -> None:
        closed = _closed_targets(CLOSED.read_text(encoding="utf-8"))
        self.assertEqual(closed, EXPECTED_CLOSED)
        disabled = tuple(
            artifact_id
            for artifact_id, entry in self.by_id.items()
            if entry["status"] == "historical_runner_disabled"
        )
        self.assertEqual(disabled, tuple(sorted(EXPECTED_CLOSED)))
        for target in EXPECTED_CLOSED:
            entry = self.by_id[target]
            self.assertEqual(entry["status"], "historical_runner_disabled")
            self.assertEqual(entry["family"], "closed_historical_operation")
            self.assertIsNone(entry["explicit_live_target_or_none"])
            self.assertEqual(
                entry["compact_verify_target_or_none"],
                CLOSED_COMPACT_TARGETS[target],
            )
            self.assertIs(entry["may_execute_state_change"], False)

    def test_transient_archive_record_is_hash_bound(self) -> None:
        entry = self.by_id[TRANSIENT_TEST_ID]
        self.assertEqual(entry["status"], "transient_test_archived")
        self.assertEqual(entry["family"], "archived_transient_test")
        self.assertIn(TRANSIENT_TEST_PATH, entry["test_paths"])
        digest = sha256((ROOT / TRANSIENT_TEST_PATH).read_bytes()).hexdigest()
        self.assertEqual(digest, TRANSIENT_TEST_SHA256)
        self.assertIn(TRANSIENT_TEST_SHA256, entry["conclusion_or_nonclaim"])
        self.assertFalse((ROOT / "tests/test_fgc_pro19_sid3_real_store_preflight.py").exists())

    def test_tdg8_owner_remains_indexed(self) -> None:
        entry = self.by_id["FGC-1-TDG8-FRZ1"]
        self.assertIn(TDG8_OWNER, entry["owner_documents"])
        self.assertTrue((ROOT / TDG8_OWNER).is_file())

    def test_result_index_stays_compact(self) -> None:
        text = RESULT_INDEX.read_text(encoding="utf-8")
        self.assertIn("| Artifact ID | Class | SHA-256 | Owner | Conclusion / nonclaim |", text)
        self.assertIn("`FGC-1-TDG10-QA2-PREF1`", text)
        self.assertIn("`FGC-1-TDG11-MSEL1-FRZ1`", text)
        self.assertIn("`FGC-1-TDG11-MSEL1-PREF1`", text)
        self.assertIn("`FGC-1-TDG11-IMP1`", text)
        self.assertIn("`FGC-1-PRO20-EV1-PREF1`", text)
        self.assertIn("`synthetic_runtime_qualification_passed_no_campaign`", text)
        self.assertIn("08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6", text)
        self.assertIn("542c195b77bf79b317dbebec88d08f2833d51cf171c99b87ea88bc3375ae9632", text)
        self.assertIn(IMP1_RESULT_SHA256, text)
        self.assertIn("3ec879fdfd309864b913b036f6b8c08c4717cb40269a18749cc19b8c0ca7b09b", text)
        self.assertIn("`UNINDEXED-CONTROLLED-MODEL`", text)
        self.assertIn("`UNINDEXED-NESTED-GRADIENT-SPECTRUM`", text)
        self.assertNotIn("This is the full artifact chronicle", text)
        self.assertIn(
            "Current frontier: `FGC-1-PRO20-EV1-PREF1`.", text
        )
        self.assertIn("Next design: `FGC-1-PRO20-EV1-RSRC1`.", text)
        self.assertIn("typed resource_exhausted stop", text)
        self.assertIn("146 accepted fine states", text)
        self.assertIn("prospective no-store per-member seed/attempt isolation", text)
        self.assertIn("no production store, namespace, authority or live run", text)
        generated = build_result_index(self.catalog).decode("utf-8")
        self.assertEqual(text, generated)
        self.assertIn("`FGC-1-TDG11-MSEL1-PREF1`", generated)
        self.assertIn("exact-C remedy", generated)

    def test_explicit_legacy_paths_are_hash_bound(self) -> None:
        for classified in self.catalog["explicit_classifications"]:
            path = ROOT / classified["path"]
            self.assertTrue(path.is_file())
            self.assertEqual(sha256(path.read_bytes()).hexdigest(), classified["sha256"])

    def test_generation_is_deterministic(self) -> None:
        first = build_catalog()
        second = build_catalog()
        self.assertEqual(first, second)
        encoded = (
            json.dumps(first, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
            + "\n"
        ).encode("ascii")
        self.assertEqual(CATALOG.read_bytes(), encoded)

    def test_readme_separates_agents_plan_and_scientific_authority(self) -> None:
        text = README.read_text(encoding="utf-8")
        self.assertGreaterEqual(len(text.splitlines()), 200)
        self.assertLessEqual(len(text.splitlines()), 300)
        self.assertIn("user-authored creative working-hypothesis orientation", text)
        self.assertIn("not evidence", text)
        self.assertIn("current implementation and finish-line execution plan", text)
        self.assertIn("validated routing and discovery", text)
        self.assertIn("does not outrank that scientific authority", text)
        self.assertNotIn("configs/fgc/artifact-catalog.json` — generated machine routing when the", text)
        for marker in SCIENTIFIC_AUTHORITY_MARKERS:
            self.assertIn(marker, text)
        agents = text.index("`AGENTS.md` — creative research orientation.")
        catalog = text.index("validated routing and discovery")
        self.assertLess(agents, catalog)

    def test_active_code_map_does_not_rank_the_catalog_above_science(self) -> None:
        text = ACTIVE_MAP.read_text(encoding="utf-8")
        self.assertIn("user-authored creative working-hypothesis orientation", text)
        self.assertIn("not evidence", text)
        self.assertIn("current implementation and finish-line execution plan", text)
        self.assertIn("validated routing and discovery", text)
        self.assertIn("not a higher scientific", text)
        self.assertIn("authority than the order above", text)
        self.assertIn("`AGENTS.md` — creative research orientation.", text)
        self.assertLess(
            text.index("`AGENTS.md` — creative research orientation."),
            text.index("validated routing and discovery"),
        )


if __name__ == "__main__":
    unittest.main()
