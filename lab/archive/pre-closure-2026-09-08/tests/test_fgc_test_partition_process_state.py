from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from scripts import run_fgc_test_partition as partition  # noqa: E402


MODULE = "tests.test_fgc_tdg11_imp1_enclosure"
REAL_CASE = (
    MODULE + ".LargeIntegerAndDerivedOverflowTests."
    "test_chunked_decimal_formatting_keeps_global_digit_limit"
)


class ProcessStateRoutingTests(unittest.TestCase):
    def test_split_preserves_every_identity_and_original_case(self) -> None:
        owned_type = type("Owned", (unittest.TestCase,), {
            "__module__": MODULE, "test_owned": lambda self: None,
        })
        ordinary_type = type("Ordinary", (unittest.TestCase,), {
            "__module__": __name__, "test_ordinary": lambda self: None,
        })
        owned, ordinary = owned_type("test_owned"), ordinary_type("test_ordinary")
        suite = unittest.TestSuite((unittest.TestSuite((owned,)), ordinary))
        retained, groups = partition._split_process_isolated_cases(suite)
        self.assertEqual(tuple(partition._iter_cases(retained)), (ordinary,))
        self.assertEqual(groups, {MODULE: (owned.id(),)})
        self.assertEqual(suite.countTestCases(), retained.countTestCases() + 1)
        self.assertEqual(tuple(partition._iter_cases(suite)), (owned, ordinary))

    def test_empty_groups_do_not_spawn_a_process(self) -> None:
        with mock.patch.object(partition.subprocess, "run") as run:
            self.assertEqual(partition._run_process_isolated_cases({}), 0)
        run.assert_not_called()

    def test_command_keeps_runtime_and_finite_limit_without_shell(self) -> None:
        completed = SimpleNamespace(returncode=0, stdout="", stderr="")
        with mock.patch.object(partition.subprocess, "run", return_value=completed) as run:
            with redirect_stdout(StringIO()):
                self.assertEqual(partition._run_process_isolated_cases({MODULE: (REAL_CASE,)}), 1)
        args, kwargs = run.call_args
        command = args[0]
        self.assertEqual(command[:5], [sys.executable, "-I", "-B", "-X", "int_max_str_digits=4300"])
        self.assertEqual(command[-2:], [str(ROOT), REAL_CASE])
        self.assertNotIn("shell", kwargs)
        self.assertFalse(kwargs["check"])
        self.assertEqual(kwargs["cwd"], ROOT)
        self.assertEqual(kwargs["timeout"], 600)

    def test_failure_output_is_retained_and_stops_the_gate(self) -> None:
        failed = SimpleNamespace(returncode=1, stdout="child stdout\n", stderr="child failure\n")
        output, errors = StringIO(), StringIO()
        with mock.patch.object(partition.subprocess, "run", return_value=failed):
            with redirect_stdout(output), redirect_stderr(errors):
                with self.assertRaisesRegex(RuntimeError, "isolated process tests failed"):
                    partition._run_process_isolated_cases({MODULE: (REAL_CASE,)})
        self.assertIn("child stdout", output.getvalue())
        self.assertIn("child failure", errors.getvalue())

    def test_timeout_is_not_treated_as_success(self) -> None:
        with mock.patch.object(partition.subprocess, "run", side_effect=subprocess.TimeoutExpired("test", 600)):
            with redirect_stdout(StringIO()):
                with self.assertRaises(subprocess.TimeoutExpired):
                    partition._run_process_isolated_cases({MODULE: (REAL_CASE,)})

    def test_unknown_foreign_duplicate_or_empty_group_is_refused_before_spawn(self) -> None:
        for groups in (
            {"tests.unknown": ("tests.unknown.Case.test",)},
            {MODULE: ("tests.foreign.Case.test",)},
            {MODULE: (REAL_CASE, REAL_CASE)},
            {MODULE: ()},
            {MODULE: [REAL_CASE]},
        ):
            with self.subTest(groups=groups), mock.patch.object(partition.subprocess, "run") as run:
                with self.assertRaises(ValueError):
                    partition._run_process_isolated_cases(groups)
                run.assert_not_called()

    def test_real_historical_import_does_not_disable_isolated_imp1_premise(self) -> None:
        before = sys.get_int_max_str_digits()
        try:
            # The historical module legitimately owns unbounded exact decimal
            # rendering. Re-importing may be cached, so reproduce its setting.
            from scripts import reproduce_fgc_hyp1_dom3_uhyp1  # noqa: F401
            sys.set_int_max_str_digits(0)
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                count = partition._run_process_isolated_cases({MODULE: (REAL_CASE,)})
            self.assertEqual(count, 1)
            self.assertEqual(sys.get_int_max_str_digits(), 0)
        finally:
            sys.set_int_max_str_digits(before)


if __name__ == "__main__":
    unittest.main()
