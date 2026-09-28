"""Small synthetic selection/publication tests, never historical campaign data."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from recursive_horizons import evidence_io as io
from recursive_horizons.fgc.evolution import tdg11_msel1_contract as contract
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.tdg11_rational_complete_c import (
    assess_rational_complete_c_rows,
)
from scripts import run_fgc_tdg11_msel1 as runner
from tests.test_fgc_tdg11_msel1_contract import fixture_config
from tests.test_fgc_tdg11_rational_complete_c import _smoothstep_row
from tests import test_fgc_tdg11_msel1_runtime as fixture


COMMIT = "a" * 40
CONFIG_SHA = "b" * 64
FREEZE_SHA = "c" * 64


def exact_evidence(a=3, b=1):
    return assess_rational_complete_c_rows(
        (_smoothstep_row(a, b),),
        expected_row_count=1,
        maximum_candidates_D01=contract.RESOURCES["per_channel_maximum_candidates_D01"],
        maximum_candidates_D12=contract.RESOURCES["per_channel_maximum_candidates_D12"],
        refinement_depth=contract.RESOURCES["primary_refinement_depth"],
    )


def width_records():
    result = [runner._new_width(spec) for spec in contract.replay_specs()]
    evidence = exact_evidence()
    embedded = SimpleNamespace(
        embedded_defect_bounds=(Q(9), Q(3), Q(1)),
        accumulation_bounds=(Q(1, 10), Q(1, 20), Q(1, 40)),
    )
    for width, spec in zip(result, contract.replay_specs(), strict=True):
        width["restored"] = {
            "checkpoint_sha256": spec["checkpoint_sha256"],
            "descriptor_sha256": contract.DESCRIPTOR_SHA256,
            "state_sha256": contract.PHYSICAL_STATE_SHA256,
            "fingerprint_sha256": "d" * 64,
        }
        for arithmetic in width["families"]:
            width["families"][arithmetic] = {
                "status": "complete",
                "family_sha256": "e" * 64,
                "accepted_source_prechecks": 7,
                "stage_and_endpoint_records": 35,
                "rhs_calls": 42,
            }
        for group in (width["baseline"], *width["candidates"].values()):
            group["status"] = "complete"
        for channel in runner.CHANNELS:
            width["baseline"]["channels"].append(
                runner._exact_channel(channel, evidence)
            )
            width["candidates"][contract.CANDIDATES[0]]["channels"].append(
                runner._exact_channel(channel, evidence, roundoff=(Q(0), Q(0), Q(1, 7)))
            )
            width["candidates"][contract.CANDIDATES[1]]["channels"].append(
                runner._exact_channel(channel, evidence)
            )
            width["candidates"][contract.CANDIDATES[2]]["channels"].append(
                runner._embedded_channel(channel, embedded)
            )
    return result


def terminal_record():
    widths = width_records()
    return {
        "schema": runner.RAW_SCHEMA,
        "artifact_id": contract.ARTIFACT_ID,
        "runner_id": runner.RUNNER_ID,
        "authority_commit": COMMIT,
        "config_sha256": CONFIG_SHA,
        "freeze_sha256": FREEZE_SHA,
        "environment": fixture_config()["environment"],
        "widths": widths,
        "store_snapshot_before": [contract.STORE_LEAF_COUNT, contract.STORE_SHA256],
        "store_snapshot_after": [contract.STORE_LEAF_COUNT, contract.STORE_SHA256],
        "accounting": {
            "accepted_source_prechecks": 42,
            "stage_and_endpoint_records": 210,
            "rhs_calls": 252,
            "static_shells": 6,
        },
        "global_stop": None,
        "elapsed_seconds_hex": (1.0).hex(),
        **runner.reduce_selection(widths, expected_rows=1),
        "nonclaims": contract.nonclaims(),
    }


def validate(terminal):
    return runner.validate_terminal(
        terminal,
        authority_commit=COMMIT,
        config_sha256=CONFIG_SHA,
        freeze_sha256=FREEZE_SHA,
        environment=fixture_config()["environment"],
        expected_rows=1,
    )


class SelectionContractTests(unittest.TestCase):
    def test_first_fully_passing_candidate_wins(self):
        terminal = terminal_record()
        self.assertEqual(validate(terminal), terminal)
        self.assertEqual(terminal["selected_candidate"], contract.CANDIDATES[0])
        self.assertEqual(terminal["classification"], runner.SELECTED_CLASS)
        self.assertFalse(terminal["licenses_only_separate_TDG11_IMP1"])
        self.assertTrue(terminal["independent_PREF1_required_before_IMP1"])
        self.assertTrue(all(value is False for value in terminal["nonclaims"].values()))

    def test_all_width_channel_gate_and_fixed_precedence(self):
        widths = width_records()
        failure = exact_evidence(2, 1)
        first = widths[2]["candidates"][contract.CANDIDATES[0]]["channels"]
        first[0] = runner._exact_channel(
            runner.CHANNELS[0], failure, roundoff=(Q(0),) * 3
        )
        reduced = runner.reduce_selection(widths, expected_rows=1)
        self.assertEqual(reduced["candidate_summaries"][0]["status"], "nonpass")
        self.assertEqual(reduced["selected_candidate"], contract.CANDIDATES[1])
        for width in widths:
            width["candidates"][contract.CANDIDATES[1]]["channels"][0] = (
                runner._exact_channel(runner.CHANNELS[0], failure)
            )
            bad = SimpleNamespace(
                embedded_defect_bounds=(Q(1), Q(1), Q(1)),
                accumulation_bounds=(Q(0),) * 3,
            )
            width["candidates"][contract.CANDIDATES[2]]["channels"][0] = (
                runner._embedded_channel(runner.CHANNELS[0], bad)
            )
        self.assertEqual(
            runner.reduce_selection(widths, expected_rows=1)["classification"],
            runner.NONPASS_CLASS,
        )

    def test_unmeasured_routes_are_inconclusive_not_failures(self):
        widths = [runner._new_width(spec) for spec in contract.replay_specs()]
        answer = runner.reduce_selection(widths, expected_rows=1)
        self.assertEqual(answer["classification"], runner.INCONCLUSIVE_CLASS)
        self.assertIsNone(answer["selected_candidate"])
        self.assertTrue(
            all(
                row["status"] == "inconclusive" for row in answer["candidate_summaries"]
            )
        )

    def test_saved_classifications_and_complete_debits_are_recomputed(self):
        record = deepcopy(
            width_records()[0]["candidates"][contract.CANDIDATES[0]]["channels"][0]
        )
        record["public_fine_debit"] = runner.exact_object(Q(1))
        with self.assertRaises(runner.MSEL1RunnerError):
            runner.channel_outcome(record, contract.CANDIDATES[0], expected_rows=1)
        record = deepcopy(
            width_records()[0]["candidates"][contract.CANDIDATES[1]]["channels"][0]
        )
        record["exact_complete_C"]["decision"]["admission_passed"] = False
        with self.assertRaises(ValueError):
            runner.channel_outcome(record, contract.CANDIDATES[1], expected_rows=1)
        record = deepcopy(
            width_records()[0]["candidates"][contract.CANDIDATES[2]]["channels"][0]
        )
        record["global_PDE_enclosure"] = True
        with self.assertRaises(runner.MSEL1RunnerError):
            runner.channel_outcome(record, contract.CANDIDATES[2], expected_rows=1)

    def test_hash_generation_width_promotion_and_partial_terminal_attacks(self):
        original = terminal_record()
        mutations = []
        item = deepcopy(original)
        item["authority_commit"] = "f" * 40
        mutations.append(item)
        item = deepcopy(original)
        item["config_sha256"] = "f" * 64
        mutations.append(item)
        item = deepcopy(original)
        item["widths"][0]["generation"] = 10
        mutations.append(item)
        item = deepcopy(original)
        item["widths"][2]["width_hex"] = item["widths"][1]["width_hex"]
        mutations.append(item)
        item = deepcopy(original)
        item["widths"][0]["restored"]["state_sha256"] = "f" * 64
        mutations.append(item)
        item = deepcopy(original)
        item["widths"][0]["families"][contract.ORIGINAL_ARITHMETIC_ID] = None
        mutations.append(item)
        item = deepcopy(original)
        item["widths"][0]["candidates"][contract.CANDIDATES[0]]["channels"].pop()
        mutations.append(item)
        item = deepcopy(original)
        item["nonclaims"]["accepted_state_advance_authorized"] = True
        mutations.append(item)
        item = deepcopy(original)
        item["accounting"]["rhs_calls"] = 251
        mutations.append(item)
        item = deepcopy(original)
        item["selected_candidate"] = contract.CANDIDATES[1]
        mutations.append(item)
        for index, item in enumerate(mutations):
            with (
                self.subTest(index=index),
                self.assertRaises(
                    ValueError
                    if index in (0, 1, 2, 3, 4, 7, 8, 9)
                    else runner.MSEL1RunnerError
                ),
            ):
                validate(item)

    def test_channel_reorder_and_extra_candidate_are_rejected(self):
        widths = width_records()
        widths[0]["candidates"][contract.CANDIDATES[1]]["channels"].reverse()
        with self.assertRaises(runner.MSEL1RunnerError):
            runner.reduce_selection(widths, expected_rows=1)
        widths = width_records()
        widths[0]["candidates"]["unfrozen"] = runner.empty_group()
        with self.assertRaises(runner.MSEL1RunnerError):
            runner.reduce_selection(widths, expected_rows=1)

    def test_exact_integer_codec_keeps_global_interpreter_limits_unchanged(self):
        import sys

        before = sys.get_int_max_str_digits()
        value = Q(2**20000 + 1, 2**18000 + 3)
        self.assertEqual(runner.exact_fraction(runner.exact_object(value)), value)
        self.assertEqual(sys.get_int_max_str_digits(), before)
        for value in (
            {"numerator": "02", "denominator": "1"},
            {"numerator": "2", "denominator": "4"},
            {"numerator": "-0", "denominator": "1"},
            {"numerator": "1", "denominator": "0"},
        ):
            with self.assertRaises(runner.MSEL1RunnerError):
                runner.exact_fraction(value)


class SyntheticInvocationTests(unittest.TestCase):
    def test_preflight_refusal_never_builds_a_shadow(self):
        fake = SimpleNamespace(
            authorize_execution=lambda *args, **kwargs: (_ for _ in ()).throw(
                runner.MSEL1RunnerError("closed")
            )
        )
        with patch.object(runner, "_authority_module", return_value=fake):
            self.assertFalse(
                runner.status(Path("/unavailable"), authority_commit=COMMIT)[
                    "safe_to_run"
                ]
            )
            with self.assertRaises(runner.MSEL1RunnerError):
                runner.run(Path("/unavailable"), authority_commit=COMMIT)

    def test_one_synthetic_run_publishes_only_terminal_and_manifest(self):
        template = SimpleNamespace(
            state=fixture._state(),
            operator=fixture._rhs(),
            projector=None,
            transaction=fixture._transaction(),
            tracers=fixture._Tracers(),
            initial=SimpleNamespace(
                grid=SimpleNamespace(coordinates=fixture._coordinates())
            ),
            time=0.0,
            step_index=0,
            transaction_serial=0,
        )
        state_hash = array_content_sha256(
            template.state.u, template.state.p, template.state.q
        )
        environment = fixture_config()["environment"]
        receipt = SimpleNamespace(
            authority_commit=COMMIT,
            config_sha256=CONFIG_SHA,
            result_sha256=FREEZE_SHA,
            environment=environment,
            store_snapshot=(contract.STORE_LEAF_COUNT, contract.STORE_SHA256),
        )
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            calls = []

            def authorize(repository, *, authority_commit):
                calls.append(authority_commit)
                if (repository / "diagnostic").exists():
                    raise runner.MSEL1RunnerError("namespace already used")
                return receipt

            def restore(repository, retry, *, templates):
                spec = next(
                    item for item in contract.replay_specs() if item["retry"] == retry
                )
                return SimpleNamespace(
                    member=deepcopy(template),
                    checkpoint_sha256=spec["checkpoint_sha256"],
                    descriptor_sha256=contract.DESCRIPTOR_SHA256,
                    fingerprint={"retry": retry, "synthetic": True},
                )

            fake = SimpleNamespace(
                authorize_execution=authorize, restore_predecessor=restore
            )
            with (
                patch.object(runner, "_authority_module", return_value=fake),
                patch.object(contract, "OWNED_ROW_COUNT", 4),
                patch.object(contract, "PHYSICAL_STATE_SHA256", state_hash),
                patch.object(contract, "OUTPUT_NAMESPACE", "diagnostic"),
                patch(
                    "recursive_horizons.fgc.evolution.proto19_gr0_static_factory.build_static_gr0_shells",
                    return_value={str(i): object() for i in range(6)},
                ),
            ):
                terminal = runner.run(root, authority_commit=COMMIT)
                self.assertEqual(
                    terminal["accounting"],
                    {
                        "static_shells": 6,
                        "accepted_source_prechecks": 42,
                        "stage_and_endpoint_records": 210,
                        "rhs_calls": 252,
                    },
                )
                self.assertEqual(
                    {item.name for item in (root / "diagnostic").iterdir()},
                    {"manifest.json", "terminal.json"},
                )
                manifest = io.load_canonical_json(
                    io.read_regular_file(root, "diagnostic/manifest.json")
                )
                from hashlib import sha256

                self.assertEqual(
                    manifest["terminal_sha256"],
                    sha256(
                        io.read_regular_file(root, "diagnostic/terminal.json")
                    ).hexdigest(),
                )
                self.assertIs(manifest["endpoint_serialized"], False)
                with self.assertRaises(runner.MSEL1RunnerError):
                    runner.run(root, authority_commit=COMMIT)
                self.assertEqual(calls, [COMMIT, COMMIT, COMMIT])


if __name__ == "__main__":
    unittest.main()
