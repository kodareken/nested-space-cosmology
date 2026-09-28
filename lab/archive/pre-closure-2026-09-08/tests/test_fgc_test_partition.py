from __future__ import annotations

from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_test_partition as partition  # noqa: E402


class FGCTestPartitionTests(unittest.TestCase):
    def test_postlaunch_contract_lists_the_exact_three_prelaunch_test_ids(self) -> None:
        config = partition._load_config()
        sealed = config["partition"]["proto14_postlaunch"][
            "sealed_prelaunch_test_ids"
        ]
        self.assertEqual(
            sealed,
            [
                "tests.test_fgc_pro14_frz1_reproduction."
                "PRO14FreezeReproductionTests.test_canonical_record_reproduces",
                "tests.test_fgc_hlt12_mon12_reproduction."
                "FGCHLT12MON12ReproductionTests."
                "test_canonical_result_is_checked_when_sealed",
                "tests.test_fgc_gr0_campaign_runner_v14."
                "FGCGR0CampaignRunnerV14Tests."
                "test_frozen_run_plan_is_strict_and_candidate_closed",
            ],
        )

    def test_sealing_removes_exactly_the_named_cases(self) -> None:
        class FirstCase(unittest.TestCase):
            def test_first(self) -> None:
                pass

        class SecondCase(unittest.TestCase):
            def test_second(self) -> None:
                pass

        suite = unittest.TestSuite(
            (FirstCase("test_first"), SecondCase("test_second"))
        )
        sealed = (FirstCase("test_first").id(),)
        routed = partition._seal_prelaunch_tests(suite, sealed)
        self.assertEqual(
            [case.id() for case in partition._iter_cases(routed)],
            [SecondCase("test_second").id()],
        )

    def test_pref27_postlaunch_contract_lists_the_exact_temporal_test_ids(self) -> None:
        config = partition._load_config()
        sealed = config["partition"]["proto18_pref27_postlaunch"][
            "sealed_prelaunch_test_ids"
        ]
        prefix = "tests.test_fgc_proto18_pref27_binder.Pref27BinderTests."
        self.assertEqual(
            sealed,
            [
                prefix + "test_actual_canonical_store_is_exact_and_build_is_nonmutating",
                prefix + "test_arriving_foreign_root_wins_exclusive_rename_without_clobber",
                prefix + "test_authority_receipt_checkpoint_state_tree_and_claim_mutations_fail_closed",
                prefix + "test_byte_identical_copy_is_admitted_but_is_not_a_reusable_result",
                prefix + "test_config_scope_and_claim_promotions_fail_before_store_access",
                prefix + "test_content_address_filename_and_empty_directory_mutations_fail_closed",
                prefix + "test_partial_file_and_symlink_store_shapes_fail_closed_without_adoption",
                prefix + "test_partial_hlt15_store_availability_is_not_silently_accepted",
                prefix + "test_repository_root_nested_directory_and_leaf_symlinks_fail_closed",
                prefix + "test_semantically_promoted_receipt_and_authority_are_rejected_after_rehash",
                prefix + "test_special_store_entry_fails_closed_when_supported",
                prefix + "test_staging_residue_is_rejected_even_when_the_eight_store_leaves_are_valid",
                "tests.test_fgc_pro18_pref27_reproduction.Pref27ReproductionTests."
                "test_fixed_pref27_result_reproduces_without_advancing_canonical_store",
            ],
        )

    def test_unknown_sealed_test_id_fails_closed(self) -> None:
        class FirstCase(unittest.TestCase):
            def test_first(self) -> None:
                pass

        suite = unittest.TestSuite((FirstCase("test_first"),))
        with self.assertRaisesRegex(ValueError, "are absent"):
            partition._seal_prelaunch_tests(suite, ("tests.missing.test",))

    def test_complete_evolution_suite_seals_exactly_the_configured_ids(self) -> None:
        config = partition._load_config()
        modules = partition._select_modules(
            partition._module_names(),
            set(config["partition"]["runtime_bound_test_modules"]),
            "evolution-restart",
        )
        suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
        sealed = config["partition"]["proto14_postlaunch"][
            "sealed_prelaunch_test_ids"
        ]
        routed = partition._seal_prelaunch_tests(suite, sealed)
        all_ids = {case.id() for case in partition._iter_cases(suite)}
        routed_ids = {case.id() for case in partition._iter_cases(routed)}
        self.assertEqual(all_ids - routed_ids, set(sealed))
        self.assertEqual(suite.countTestCases() - routed.countTestCases(), 3)

    def test_complete_analytic_suite_seals_exactly_pref27_temporal_cases(self) -> None:
        config = partition._load_config()
        modules = partition._select_modules(
            partition._module_names(),
            set(config["partition"]["runtime_bound_test_modules"]),
            "analytic",
        )
        suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
        sealed = config["partition"]["proto18_pref27_postlaunch"][
            "sealed_prelaunch_test_ids"
        ]
        routed = partition._seal_prelaunch_tests(
            suite, sealed, seal_id="FGC-1-PRO18-PREF27"
        )
        all_ids = {case.id() for case in partition._iter_cases(suite)}
        routed_ids = {case.id() for case in partition._iter_cases(routed)}
        self.assertEqual(all_ids - routed_ids, set(sealed))
        self.assertEqual(suite.countTestCases() - routed.countTestCases(), 13)

    def test_cal11_postlaunch_verifier_is_check_only(self) -> None:
        config = partition._load_config()
        completed = SimpleNamespace(returncode=0, stdout="bound\n", stderr="")
        with mock.patch.object(
            partition.subprocess, "run", return_value=completed
        ) as run:
            partition._run_cal11_pref22_check(config)
        command = run.call_args.args[0]
        self.assertIn("--check", command)
        self.assertNotIn("--resume", command)
        self.assertNotIn("--run", command)
        self.assertEqual(run.call_args.kwargs["cwd"], REPOSITORY)

    def test_cal11_failure_prevents_postlaunch_routing(self) -> None:
        config = partition._load_config()
        failed = SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="terminal binding differs",
        )
        with mock.patch.object(partition.subprocess, "run", return_value=failed):
            with self.assertRaisesRegex(RuntimeError, "terminal binding differs"):
                partition._run_cal11_pref22_check(config)

    def test_ur1_pref1_postlaunch_contract_lists_the_exact_ids(self) -> None:
        config = partition._load_config()
        section = config["partition"]["tdg9_ur1_pref1_postlaunch"]
        self.assertEqual(
            set(section),
            {
                "ur1_pref1_reproducer",
                "sealed_prelaunch_test_ids",
                "evolution_restart_reconstruction_test_ids",
            },
        )
        self.assertEqual(
            section["ur1_pref1_reproducer"],
            "scripts/reproduce_fgc_tdg9_ur1_pref1.py",
        )
        self.assertEqual(
            section["sealed_prelaunch_test_ids"],
            [
                "tests.test_fgc_tdg8_rcv3_execution_authority."
                "TDG8RCV3ExecutionAuthorityTests."
                "test_live_projection_is_a_safe_generation_nine_retry_boundary",
                "tests.test_fgc_tdg8_rcv3_freeze.TDG8RCV3FreezeTests."
                "test_repair_sources_are_still_exact_commit_bytes",
                "tests.test_fgc_tdg8_rcv3_pref1_binder."
                "TDG8RCV3PREF1BinderTests."
                "test_live_projection_passes_independent_full_binding",
                "tests.test_fgc_tdg8_rcv3_rec1_authority."
                "TDG8RCV3REC1AuthorityTests."
                "test_exact_entry_and_two_phase_recovery_use_only_a_copy",
                "tests.test_fgc_tdg8_rcv3_rec1_authority."
                "TDG8RCV3REC1AuthorityTests."
                "test_owner_token_is_absent_from_tracked_contract",
                "tests.test_fgc_tdg8_rcv3_rec1_authority."
                "TDG8RCV3REC1AuthorityTests."
                "test_repeated_takeover_crash_cuts_remain_bounded",
                "tests.test_fgc_tdg8_rcv3_rec1_event_runner."
                "TDG8RCV3REC1RunnerTests."
                "test_real_copied_store_reconciles_to_clean_exact_gen10_without_PDE",
                "tests.test_fgc_tdg9_ac1_authority.AC1AuthorityTests."
                "test_output_absence_is_prelaunch_only",
                "tests.test_fgc_tdg9_ti2_authority.TI2AuthorityTests."
                "test_old_and_new_output_absence_is_prelaunch_only",
            ],
        )
        self.assertEqual(len(section["sealed_prelaunch_test_ids"]), 9)
        self.assertEqual(
            section["evolution_restart_reconstruction_test_ids"],
            [
                "tests.test_fgc_tdg9_loc1_err1.LOC1ERR1Tests."
                "test_sealed_tuple_is_reconstructed_without_loc2",
                "tests.test_fgc_tdg9_loc2_runner.LOC2RunnerTests."
                "test_retry3_value_D01_hash_bound_family_matches",
            ],
        )
        self.assertEqual(
            len(section["evolution_restart_reconstruction_test_ids"]), 2
        )
        self.assertFalse(
            set(section["sealed_prelaunch_test_ids"])
            & set(section["evolution_restart_reconstruction_test_ids"])
        )
        self.assertEqual(
            set(config["partition"]["proto18_pref27_postlaunch"]),
            {"sealed_prelaunch_test_ids"},
        )
        self.assertEqual(
            len(
                config["partition"]["proto18_pref27_postlaunch"][
                    "sealed_prelaunch_test_ids"
                ]
            ),
            13,
        )

    def test_ur1_pref1_owned_ids_are_analytic_not_runtime_bound(self) -> None:
        config = partition._load_config()
        bound = set(config["partition"]["runtime_bound_test_modules"])
        section = config["partition"]["tdg9_ur1_pref1_postlaunch"]
        owned = (
            section["sealed_prelaunch_test_ids"]
            + section["evolution_restart_reconstruction_test_ids"]
        )
        for test_id in owned:
            self.assertNotIn(test_id.split(".")[1], bound)

    def test_reconstruction_routing_removes_exactly_the_named_cases(self) -> None:
        class SiblingTests(unittest.TestCase):
            def test_compact_sibling(self) -> None:
                pass

            def test_reconstruction(self) -> None:
                pass

        suite = unittest.TestSuite(
            (
                SiblingTests("test_compact_sibling"),
                SiblingTests("test_reconstruction"),
            )
        )
        routed_id = SiblingTests("test_reconstruction").id()
        retained, extracted = partition._route_reconstruction_cases(
            suite, (routed_id,)
        )
        self.assertEqual(
            [case.id() for case in partition._iter_cases(retained)],
            [SiblingTests("test_compact_sibling").id()],
        )
        self.assertEqual([case.id() for case in extracted], [routed_id])
        self.assertEqual(suite.countTestCases() - retained.countTestCases(), 1)

    def test_unknown_reconstruction_test_id_fails_closed(self) -> None:
        class FirstCase(unittest.TestCase):
            def test_first(self) -> None:
                pass

        suite = unittest.TestSuite((FirstCase("test_first"),))
        with self.assertRaisesRegex(ValueError, "are absent"):
            partition._route_reconstruction_cases(
                suite, ("tests.missing.test_reconstruction",)
            )

    def test_append_reconstruction_cases_adds_exactly_the_extracted_tests(
        self,
    ) -> None:
        class AnalyticCase(unittest.TestCase):
            def test_reconstruction(self) -> None:
                pass

        class EvolutionCase(unittest.TestCase):
            def test_runtime(self) -> None:
                pass

        analytic = unittest.TestSuite((AnalyticCase("test_reconstruction"),))
        evolution = unittest.TestSuite((EvolutionCase("test_runtime"),))
        retained, extracted = partition._route_reconstruction_cases(
            analytic, (AnalyticCase("test_reconstruction").id(),)
        )
        self.assertEqual(retained.countTestCases(), 0)
        combined = partition._append_reconstruction_cases(evolution, extracted)
        self.assertEqual(
            [case.id() for case in partition._iter_cases(combined)],
            [
                EvolutionCase("test_runtime").id(),
                AnalyticCase("test_reconstruction").id(),
            ],
        )
        self.assertEqual(
            combined.countTestCases() - evolution.countTestCases(), 1
        )

    def test_append_does_not_duplicate_existing_evolution_cases(self) -> None:
        class FirstCase(unittest.TestCase):
            def test_first(self) -> None:
                pass

        case = FirstCase("test_first")
        suite = unittest.TestSuite((case,))
        with self.assertRaisesRegex(ValueError, "already present"):
            partition._append_reconstruction_cases(suite, (FirstCase("test_first"),))

    def test_complete_analytic_suite_routes_exactly_the_reconstruction_cases(
        self,
    ) -> None:
        config = partition._load_config()
        modules = partition._select_modules(
            partition._module_names(),
            set(config["partition"]["runtime_bound_test_modules"]),
            "analytic",
        )
        suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
        routed = config["partition"]["tdg9_ur1_pref1_postlaunch"][
            "evolution_restart_reconstruction_test_ids"
        ]
        retained, extracted = partition._route_reconstruction_cases(suite, routed)
        all_ids = {case.id() for case in partition._iter_cases(suite)}
        retained_ids = {case.id() for case in partition._iter_cases(retained)}
        extracted_ids = {case.id() for case in extracted}
        self.assertEqual(all_ids - retained_ids, set(routed))
        self.assertEqual(extracted_ids, set(routed))
        self.assertEqual(tuple(case.id() for case in extracted), tuple(routed))
        self.assertEqual(suite.countTestCases() - retained.countTestCases(), 2)

    def test_evolution_loads_exactly_two_analytic_reconstruction_cases(
        self,
    ) -> None:
        config = partition._load_config()
        extracted = partition._analytic_reconstruction_cases(
            config,
            set(config["partition"]["runtime_bound_test_modules"]),
        )
        routed = config["partition"]["tdg9_ur1_pref1_postlaunch"][
            "evolution_restart_reconstruction_test_ids"
        ]
        self.assertEqual(tuple(case.id() for case in extracted), tuple(routed))
        self.assertEqual(len(extracted), 2)

    def test_complete_analytic_suite_seals_exactly_ur1_pref1_prelaunch_cases(
        self,
    ) -> None:
        config = partition._load_config()
        modules = partition._select_modules(
            partition._module_names(),
            set(config["partition"]["runtime_bound_test_modules"]),
            "analytic",
        )
        suite = unittest.defaultTestLoader.loadTestsFromNames(modules)
        sealed = config["partition"]["tdg9_ur1_pref1_postlaunch"][
            "sealed_prelaunch_test_ids"
        ]
        routed = partition._seal_prelaunch_tests(
            suite, sealed, seal_id="FGC-1-TDG9-UR1-PREF1"
        )
        all_ids = {case.id() for case in partition._iter_cases(suite)}
        routed_ids = {case.id() for case in partition._iter_cases(routed)}
        self.assertEqual(all_ids - routed_ids, set(sealed))
        self.assertEqual(suite.countTestCases() - routed.countTestCases(), 9)

    def test_ur1_pref1_postlaunch_verifier_is_check_only_and_not_live(
        self,
    ) -> None:
        config = partition._load_config()
        completed = SimpleNamespace(returncode=0, stdout="bound\n", stderr="")
        with mock.patch.object(
            partition.subprocess, "run", return_value=completed
        ) as run:
            partition._run_ur1_pref1_check(config)
        command = run.call_args.args[0]
        self.assertEqual(len(command), 2)
        self.assertEqual(command[0], partition.sys.executable)
        self.assertEqual(
            Path(command[1]),
            REPOSITORY / "scripts/reproduce_fgc_tdg9_ur1_pref1.py",
        )
        self.assertNotIn("--live", command)
        self.assertNotIn("--run", command)
        self.assertNotIn("--status", command)
        self.assertFalse(any(item.lstrip("-") == "status" for item in command[1:]))
        self.assertEqual(run.call_args.kwargs["cwd"], REPOSITORY)

    def test_ur1_pref1_failure_prevents_postlaunch_routing(self) -> None:
        config = partition._load_config()
        class FirstCase(unittest.TestCase):
            def test_first(self) -> None:
                pass

        suite = unittest.TestSuite((FirstCase("test_first"),))
        failed = SimpleNamespace(
            returncode=1,
            stdout="",
            stderr="compact binding differs",
        )
        with mock.patch.object(partition.subprocess, "run", return_value=failed):
            with mock.patch.object(
                partition, "_seal_prelaunch_tests"
            ) as seal, mock.patch.object(
                partition, "_route_reconstruction_cases"
            ) as route:
                with self.assertRaisesRegex(
                    RuntimeError, "compact binding differs"
                ):
                    partition._apply_tdg9_ur1_pref1_analytic_routing(config, suite)
        seal.assert_not_called()
        route.assert_not_called()


if __name__ == "__main__":
    unittest.main()
