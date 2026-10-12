"""Prospective metadata controls; real inputs are compact references only."""

from __future__ import annotations

import ast
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import tomllib
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import tdg11_msel1_contract as contract


ROOT = Path(__file__).resolve().parents[1]


def fixture_config():
    environment = dict(contract.MINIMUM_ENVIRONMENT)
    environment.update(
        dict.fromkeys(contract.ENVIRONMENT_IMAGE_KEYS, "synthetic_fixture")
    )
    for name in ("python_executable_sha256", "numpy_extension_sha256"):
        environment[name] = "0" * 64
    return contract.expected_config(
        environment=environment,
        implementation=[
            {"path": path, "sha256": "1" * 64} for path in contract.IMPLEMENTATION_PATHS
        ],
    )


class MSEL1ContractTests(unittest.TestCase):
    def test_config_roundtrip_is_deterministic_and_nonpromoting(self):
        config = fixture_config()
        raw = contract.render_config(config)
        self.assertEqual(contract.parse_config(raw), config)
        self.assertEqual(contract.render_config(contract.parse_config(raw)), raw)
        self.assertTrue(all(value is False for value in config["nonclaims"].values()))
        self.assertEqual(config["source"]["measurement_tableau"], "RK4")
        self.assertFalse(config["future_origin"]["used_as_measurement_state"])
        self.assertEqual(config["future_origin"]["time"], "23/16")

    def test_all_reference_hashes_match_unchanged_tracked_inputs(self):
        for relative, digest in contract.PINNED_REFERENCES.items():
            self.assertEqual(
                sha256((ROOT / relative).read_bytes()).hexdigest(), digest, relative
            )

    def test_retry_predecessors_are_not_their_rejections(self):
        qa1 = tomllib.loads(
            (ROOT / "configs/fgc/fgc-1-tdg10-qa1-frz1.toml").read_text()
        )["predecessor"]
        qa2 = tomllib.loads(
            (ROOT / "configs/fgc/fgc-1-tdg10-qa2-frz1.toml").read_text()
        )["predecessor"]
        for replay in contract.replay_specs():
            retry = replay["retry"]
            prior = qa1 if retry == 3 else qa2[f"retry_{retry}"]
            self.assertEqual(replay["generation"], prior["generation"])
            self.assertEqual(replay["checkpoint_sha256"], prior["checkpoint_sha256"])
            self.assertEqual(
                replay["checkpoint_raw_sha256"], prior["checkpoint_raw_sha256"]
            )
            self.assertEqual(
                replay["journal_tip_sequence"], prior["checkpoint_journal_sequence"]
            )
            self.assertEqual(
                replay["journal_tip_sha256"], prior["checkpoint_journal_sha256"]
            )
            self.assertEqual(
                replay["journal_tip_raw_sha256"], prior["checkpoint_journal_raw_sha256"]
            )
            self.assertEqual(
                replay["rejection_sequence"], replay["journal_tip_sequence"] + 1
            )
            self.assertEqual(replay["width_hex"], prior["attempted_width_binary64_hex"])
            rejection_prefix = (
                "historical_retry3_rejection" if retry == 3 else "historical_rejection"
            )
            self.assertEqual(
                replay["rejection_sha256"], prior[f"{rejection_prefix}_sha256"]
            )
            self.assertEqual(
                replay["rejection_raw_sha256"], prior[f"{rejection_prefix}_raw_sha256"]
            )
        self.assertEqual(
            tuple(item["generation"] for item in contract.replay_specs()), (9, 10, 11)
        )
        self.assertNotEqual(
            float.fromhex(contract.replay_specs()[2]["width_hex"]),
            float.fromhex(contract.replay_specs()[1]["width_hex"]) / 2,
        )

    def test_candidate_order_width_channel_and_budget_mutations_fail(self):
        for section, key, changed in (
            ("selection", "candidate_precedence", list(reversed(contract.CANDIDATES))),
            (
                "selection",
                "channels",
                list(reversed(fixture_config()["selection"]["channels"])),
            ),
            ("selection", "minimum_observed_order", "1"),
            ("selection", "squared_order_multiplier", 4),
            ("resources", "proposals", 41),
            ("resources", "stage_and_endpoint_records", 209),
            ("resources", "selected_widths", 4),
            ("source", "measurement_tableau", "SSPRK3"),
            ("future_origin", "used_as_measurement_state", True),
        ):
            value = fixture_config()
            value[section][key] = changed
            with (
                self.subTest(section=section, key=key),
                self.assertRaises(contract.MSEL1ContractError),
            ):
                contract.validate_config(value)
        for index in range(3):
            value = fixture_config()
            value["replays"][index]["generation"] += 1
            with self.assertRaises(contract.MSEL1ContractError):
                contract.validate_config(value)

    def test_bool_int_aliases_unknown_keys_and_promotion_are_rejected(self):
        for changed in (True, 1.0, "1"):
            value = fixture_config()
            value["schema_version"] = changed
            with self.assertRaises(contract.MSEL1ContractError):
                contract.validate_config(value)
        for key in contract.nonclaims():
            for changed in (True, 0, None):
                value = fixture_config()
                value["nonclaims"][key] = changed
                with self.assertRaises(contract.MSEL1ContractError):
                    contract.validate_config(value)
        value = fixture_config()
        value["unfrozen_override"] = True
        with self.assertRaises(contract.MSEL1ContractError):
            contract.validate_config(value)

    def test_environment_and_implementation_inventory_are_exact(self):
        for mutation in ("version", "hash", "path", "omitted", "reordered"):
            value = fixture_config()
            if mutation == "version":
                value["environment"]["python_version"] = "3.14.6"
            elif mutation == "hash":
                value["environment"]["numpy_extension_sha256"] = "not_a_hash"
            elif mutation == "path":
                value["implementation"][0]["path"] = "../foreign.py"
            elif mutation == "omitted":
                value["implementation"].pop()
            else:
                value["implementation"].reverse()
            with (
                self.subTest(mutation=mutation),
                self.assertRaises(contract.MSEL1ContractError),
            ):
                contract.validate_config(value)

    def test_compact_control_does_not_read_a_store_git_or_shadow(self):
        raw = contract.render_config(fixture_config())
        with (
            patch("builtins.open", side_effect=AssertionError("unexpected file read")),
            patch("subprocess.Popen", side_effect=AssertionError("unexpected process")),
        ):
            expected = contract.compact_freeze(raw)
            self.assertEqual(contract.validate_compact(raw, expected), expected)
        self.assertIs(expected["diagnostic_executed"], False)
        self.assertIs(expected["selection_result_earned"], False)
        altered = deepcopy(expected)
        altered["selection_result_earned"] = True
        with self.assertRaises(contract.MSEL1ContractError):
            contract.validate_compact(raw, altered)

    def test_exact_main_and_embedded_orders_and_hermite_controls(self):
        controls = contract.mathematical_controls()
        records = {item["name"]: item for item in controls["tableaux"]}
        self.assertEqual([records[name]["order"] for name in records], [4, 3, 3, 2])
        self.assertEqual(records["RK4_embedded"]["order_values"][-2:], ["7/72", "1/36"])
        self.assertEqual(records["SSPRK3_embedded"]["order_values"][2], "1/2")
        for name in (
            "exact_threshold_equality_passes",
            "threshold_below_equality_fails",
            "no_PDE_stability_theorem_inferred",
        ):
            self.assertIs(controls[name], True)
        for t in (Q(0), Q(1, 7), Q(1, 2), Q(2, 3), Q(1)):
            for e0, e1 in ((Q(-2), Q(3)), (Q(1), Q(-1)), (Q(0), Q(7))):
                value = e0 * (1 - 3 * t * t + 2 * t**3) + e1 * (3 * t * t - 2 * t**3)
                self.assertLessEqual(abs(value), max(abs(e0), abs(e1)))

    def test_strict_bytes_parser_and_data_copies(self):
        raw = contract.render_config(fixture_config())
        for invalid in (raw.decode(), b"[invalid", b"\xff", b"x" * (1024 * 1024 + 1)):
            with self.assertRaises(contract.MSEL1ContractError):
                contract.parse_config(invalid)
        copy = contract.validate_config(fixture_config())
        copy["nonclaims"]["physical_result_earned"] = True
        self.assertIs(contract.nonclaims()["physical_result_earned"], False)

    def test_contract_has_no_execution_or_raw_import_dependency(self):
        source = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_msel1_contract.py"
        imports = [
            node.module or ""
            for node in ast.walk(ast.parse(source.read_text()))
            if isinstance(node, ast.ImportFrom)
        ]
        for name in imports:
            self.assertFalse(
                any(
                    part in name
                    for part in (
                        "runner",
                        "authority",
                        "runtime",
                        "numpy",
                        "campaign",
                        "reconstruction",
                    )
                ),
                name,
            )


if __name__ == "__main__":
    unittest.main()
