"""Independent binding composition and adversarial compact provenance controls."""

from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons import evidence_io as io
from recursive_horizons.fgc.evolution import tdg11_msel1_contract as freeze
from recursive_horizons.fgc.evolution import tdg11_msel1_pref1_binder as b
from recursive_horizons.fgc.evolution import tdg11_msel1_pref1_protocol as p
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_localization import (
    assess_complete_c_independently,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_reconstruction import (
    reconstruct_channel_independently,
)
from tests import test_fgc_tdg11_msel1_pref1_protocol as synthetic
from tests import test_fgc_tdg11_msel1_pref1_reconstruction as reconstruction


ROOT = Path(__file__).resolve().parents[1]


def synthetic_terminal():
    """A wire-only, nonphysical fixture; no real state, RHS, or campaign input."""
    result = synthetic.terminal()
    for width in result["widths"]:
        for candidate, group in [
            ("baseline", width["baseline"]),
            *width["candidates"].items(),
        ]:
            if candidate == freeze.CANDIDATES[2]:
                continue
            for channel in group["channels"]:
                evidence = channel["exact_complete_C"]
                evidence["row_count"] = evidence["expected_row_count"] = (
                    freeze.OWNED_ROW_COUNT
                )
                for key, multiplier in (("d01", 2), ("d12", 4)):
                    count = multiplier * freeze.OWNED_ROW_COUNT
                    evidence[key].update(
                        polynomial_count=count,
                        candidate_count=2 * count,
                        survivor_count=2 * count,
                    )
    result.update(
        p.reduce_widths(result["widths"], expected_rows=freeze.OWNED_ROW_COUNT)
    )
    return result


def independent_widths(raw):
    widths = deepcopy(raw["widths"])
    for width in widths:
        for candidate, group in [
            ("baseline", width["baseline"]),
            *width["candidates"].items(),
        ]:
            if candidate != freeze.CANDIDATES[2]:
                for channel in group["channels"]:
                    channel["exact_complete_C"]["evaluator_id"] = (
                        b.INDEPENDENT_COMPLETE_C_ID
                    )
    return widths


def synthetic_bundle():
    terminal = synthetic_terminal()
    terminal_hash = sha256(io.canonical_json_bytes(terminal)).hexdigest()
    config = {
        "raw": {"terminal_sha256": terminal_hash},
        "environment": synthetic.ENVIRONMENT.copy(),
    }
    manifest = b._expected_manifest(config)
    config["raw"]["manifest_sha256"] = sha256(
        io.canonical_json_bytes(manifest)
    ).hexdigest()
    config_raw = b"synthetic fixture; never a production config\n"
    result = b._compact_payload(
        config_raw, config, terminal, independent_widths(terminal)
    )
    return config_raw, config, terminal, result


class PREF1BindingTests(unittest.TestCase):
    def test_large_base_composes_independent_reconstruction_and_localization(self):
        _paths, family = reconstruction._family(
            reconstruction.propose_step,
            reconstruction.ORIGINAL_ARITHMETIC_ID,
            reconstruction.PRIMARY_METHOD,
            reconstruction._state(u=reconstruction.LARGE_BASE),
            reconstruction._unit_u_rhs(),
        )
        data = reconstruct_channel_independently(family, "u:alpha")
        with patch.object(freeze, "OWNED_ROW_COUNT", reconstruction.OWNED_ROWS):
            baseline = b._channel_record(
                "baseline", "u:alpha", data, assess_complete_c_independently
            )
            corrected = b._channel_record(
                freeze.CANDIDATES[0], "u:alpha", data, assess_complete_c_independently
            )
            embedded = b._channel_record(
                freeze.CANDIDATES[2], "u:alpha", data, assess_complete_c_independently
            )
        self.assertEqual(
            p.fraction(baseline["exact_complete_C"]["d01"]["lower"]), Q(1, 36)
        )
        self.assertEqual(
            p.fraction(baseline["exact_complete_C"]["d12"]["upper"]), Q(1, 72)
        )
        self.assertEqual(p.fraction(corrected["exact_complete_C"]["d01"]["upper"]), 0)
        self.assertEqual(p.fraction(corrected["exact_complete_C"]["d12"]["upper"]), 0)
        self.assertEqual(p.fraction(corrected["public_fine_debit"]), Q(1, 4))
        self.assertEqual(p.fraction(embedded["public_fine_debit"]), Q(1, 4))
        self.assertEqual(
            corrected["exact_complete_C"]["evaluator_id"], b.INDEPENDENT_COMPLETE_C_ID
        )
        self.assertEqual(
            b._raw_channel(corrected, freeze.CANDIDATES[0])["exact_complete_C"][
                "evaluator_id"
            ],
            b.RAW_COMPLETE_C_ID,
        )

    def test_id_mapping_cannot_modify_math_or_masquerade_as_independence(self):
        raw = synthetic_terminal()
        widths = independent_widths(raw)
        self.assertEqual(b._raw_widths(widths), raw["widths"])
        self.assertEqual(
            widths[0]["baseline"]["channels"][0]["exact_complete_C"]["evaluator_id"],
            b.INDEPENDENT_COMPLETE_C_ID,
        )
        with self.assertRaises(b.TDG11PREF1Error):
            b._raw_widths(raw["widths"])
        changed = deepcopy(widths)
        changed[0]["baseline"]["channels"][0]["public_fine_debit"] = p.wire(Q(7))
        self.assertEqual(
            p.fraction(
                b._raw_widths(changed)[0]["baseline"]["channels"][0][
                    "public_fine_debit"
                ]
            ),
            7,
        )

    def test_compact_reconstructs_full_raw_commitment_without_live_operations(self):
        config_raw, config, _terminal, result = synthetic_bundle()
        raw = io.canonical_json_bytes(result) + b"\n"
        with (
            patch.object(b, "_tracked_inputs", return_value=config),
            patch.object(b, "_git", side_effect=AssertionError("Git")),
            patch.object(b, "_raw_snapshot", side_effect=AssertionError("raw")),
            patch.object(
                b, "_replay_independently", side_effect=AssertionError("shadow")
            ),
            patch("os.scandir", side_effect=AssertionError("store")),
        ):
            verified = b.validate_compact_result(config_raw, raw, ROOT)
        self.assertEqual(verified["classification"], p.SELECTED)
        self.assertTrue(verified["conclusion"]["licenses_only_separate_TDG11_IMP1"])
        self.assertFalse(
            verified["conclusion"]["production_method_selected_or_implemented"]
        )

    def test_self_consistent_changed_hash_receipt_still_fails_raw_commitment(self):
        config_raw, config, terminal, _result = synthetic_bundle()
        widths = independent_widths(terminal)
        widths[0]["baseline"]["channels"][0]["exact_complete_C"][
            "row_stream_sha256"
        ] = "a" * 64
        # All typed inequalities and the apparent selection remain consistent.
        forged = b._compact_payload(config_raw, config, terminal, widths)
        with (
            patch.object(b, "_tracked_inputs", return_value=config),
            self.assertRaisesRegex(b.TDG11PREF1Error, "full raw terminal commitment"),
        ):
            b.validate_compact_result(config_raw, io.canonical_json_bytes(forged), ROOT)

    def test_compact_metadata_and_state_promotion_mutations_fail(self):
        config_raw, config, _terminal, original = synthetic_bundle()
        for path, replacement in [
            (("authority", "commit"), "0" * 40),
            (("config_sha256",), "0" * 64),
            (("raw_commitment", "manifest_sha256"), "0" * 64),
            (("independent_replay", "source_store_unchanged"), 1),
            (("independent_replay", "complete_guarded_families_recomputed"), 5),
            (("conclusion", "diagnostic_endpoint_may_be_adopted"), True),
            (("conclusion", "physical_claim"), True),
        ]:
            value = deepcopy(original)
            cursor = value
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = replacement
            with (
                self.subTest(path=path),
                patch.object(b, "_tracked_inputs", return_value=config),
                self.assertRaises((b.TDG11PREF1Error, p.PREF1ProtocolError)),
            ):
                b.validate_compact_result(
                    config_raw, io.canonical_json_bytes(value), ROOT
                )

    def test_raw_namespace_is_exact_exclusive_and_hash_bound(self):
        _config_raw, config, terminal, _result = synthetic_bundle()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            parent = root / freeze.OUTPUT_NAMESPACE
            parent.mkdir(parents=True)
            (parent / "manifest.json").write_bytes(
                io.canonical_json_bytes(b._expected_manifest(config))
            )
            (parent / "terminal.json").write_bytes(io.canonical_json_bytes(terminal))
            b._raw_snapshot(root, config)
            (parent / "endpoint.json").write_text("{}")
            with self.assertRaises(b.TDG11PREF1Error):
                b._raw_snapshot(root, config)
            (parent / "endpoint.json").unlink()
            (parent / "terminal.json").write_text("{}")
            with self.assertRaises(b.TDG11PREF1Error):
                b._raw_snapshot(root, config)

    def test_config_pins_scope_and_canonical_encoding_fail_closed(self):
        raw = b.emit_config_bytes(ROOT)
        frozen_raw = (ROOT / freeze.CONFIG_PATH).read_bytes()
        config = b.parse_config(raw, frozen_raw)
        for path, replacement in [
            (("authority_commit",), "0" * 40),
            (("raw", "terminal_sha256"), "0" * 64),
            (("limits", "maximum_wall_seconds"), 14401),
            (("environment", "python_version"), "0.0.0"),
            (("candidate_precedence",), list(reversed(freeze.CANDIDATES))),
            (("schema_version",), True),
        ]:
            changed = deepcopy(config)
            cursor = changed
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = replacement
            with (
                self.subTest(path=path),
                self.assertRaises((b.TDG11PREF1Error, p.PREF1ProtocolError)),
            ):
                b.parse_config(b.render_config(changed), frozen_raw)
        with self.assertRaises(b.TDG11PREF1Error):
            b.parse_config(raw + b"\n", frozen_raw)

    def test_authority_authenticates_immutable_commit_without_head_equality(self):
        if not (ROOT / ".git").exists():
            self.skipTest("live immutable-Git proof requires repository metadata")
        original = b._git
        observed = []

        def query(root, operation):
            observed.append(operation)
            self.assertFalse(
                isinstance(operation, io.ResolveCommit) and operation.revision == "HEAD"
            )
            return original(root, operation)

        with patch.object(b, "_git", side_effect=query):
            receipt = b._authenticate_authority(ROOT)
        self.assertEqual(receipt, b._authority_receipt())
        self.assertTrue(any(isinstance(item, io.CommitParents) for item in observed))

    def test_authority_rejects_parent_delta_and_live_source_mutation(self):
        if not (ROOT / ".git").exists():
            self.skipTest("live immutable-Git proof requires repository metadata")
        original = b._git
        for operation_type in (io.CommitParents, io.InspectDelta, io.InspectTree):

            def query(root, operation, targeted=operation_type):
                result = original(root, operation)
                if isinstance(operation, targeted):
                    return ("0" * 40,) if targeted is io.CommitParents else result[:-1]
                return result

            with (
                self.subTest(operation=operation_type),
                patch.object(b, "_git", side_effect=query),
                self.assertRaises(b.TDG11PREF1Error),
            ):
                b._authenticate_authority(ROOT)
        read = io.read_regular_file

        def altered(root, relative, **kwargs):
            raw = read(root, relative, **kwargs)
            return (
                raw + b"\n# altered instrument\n"
                if relative.endswith("/numerical_engine.py")
                else raw
            )

        with (
            patch.object(io, "read_regular_file", side_effect=altered),
            self.assertRaises(b.TDG11PREF1Error),
        ):
            b._authenticate_authority(ROOT)

    def test_live_orchestration_rechecks_inputs_and_propagates_failed_replay(self):
        from recursive_horizons.fgc.evolution import tdg11_msel1_authority as authority

        config_raw, config, terminal, _result = synthetic_bundle()
        source = (freeze.STORE_LEAF_COUNT, freeze.STORE_SHA256)
        for failure in (None, b.TDG11PREF1ResourceStop("synthetic budget stop")):
            with self.subTest(failure=failure), ExitStack() as stack:
                stack.enter_context(
                    patch.object(b, "require_live_config", return_value=None)
                )
                stack.enter_context(
                    patch.object(b, "_tracked_inputs", return_value=config)
                )
                auth = stack.enter_context(
                    patch.object(
                        b,
                        "_authenticate_authority",
                        return_value=b._authority_receipt(),
                    )
                )
                snapshot = stack.enter_context(
                    patch.object(
                        b,
                        "_raw_snapshot",
                        return_value=({"synthetic": b"bytes"}, {}, terminal),
                    )
                )
                stack.enter_context(
                    patch.object(
                        authority,
                        "observe_environment",
                        return_value=config["environment"],
                    )
                )
                store = stack.enter_context(
                    patch.object(authority, "snapshot_store", return_value=source)
                )
                stack.enter_context(
                    patch.object(authority, "inspect_predecessors", return_value=())
                )
                replay = stack.enter_context(
                    patch.object(
                        b,
                        "_replay_independently",
                        side_effect=failure,
                        return_value=independent_widths(terminal),
                    )
                )
                if failure is None:
                    result = b.bind_raw_result(config_raw, ROOT)
                    self.assertEqual(result["classification"], p.SELECTED)
                    self.assertEqual(
                        (auth.call_count, snapshot.call_count, store.call_count),
                        (2, 2, 2),
                    )
                else:
                    with self.assertRaises(b.TDG11PREF1ResourceStop):
                        b.bind_raw_result(config_raw, ROOT)
                self.assertEqual(replay.call_count, 1)

    def test_live_config_drift_during_replay_prevents_a_result(self):
        from recursive_horizons.fgc.evolution import tdg11_msel1_authority as authority

        config_raw, config, terminal, _result = synthetic_bundle()
        for remove in (False, True):
            with (
                self.subTest(remove=remove),
                tempfile.TemporaryDirectory() as directory,
                ExitStack() as stack,
            ):
                root = Path(directory).resolve()
                path = root / b.CONFIG_PATH
                path.parent.mkdir(parents=True)
                path.write_bytes(config_raw)
                stack.enter_context(
                    patch.object(b, "_tracked_inputs", return_value=config)
                )
                stack.enter_context(
                    patch.object(
                        b,
                        "_authenticate_authority",
                        return_value=b._authority_receipt(),
                    )
                )
                stack.enter_context(
                    patch.object(
                        b,
                        "_raw_snapshot",
                        return_value=({"fixture": b"bytes"}, {}, terminal),
                    )
                )
                stack.enter_context(
                    patch.object(
                        authority,
                        "observe_environment",
                        return_value=config["environment"],
                    )
                )
                stack.enter_context(
                    patch.object(
                        authority,
                        "snapshot_store",
                        return_value=(freeze.STORE_LEAF_COUNT, freeze.STORE_SHA256),
                    )
                )
                stack.enter_context(
                    patch.object(authority, "inspect_predecessors", return_value=())
                )

                def change_config(*_args):
                    if remove:
                        path.unlink()
                    else:
                        path.write_bytes(config_raw + b"changed after capture")
                    return independent_widths(terminal)

                stack.enter_context(
                    patch.object(b, "_replay_independently", side_effect=change_config)
                )
                with self.assertRaisesRegex(b.TDG11PREF1Error, "live PREF1 config"):
                    b.bind_raw_result(config_raw, root)
                self.assertFalse((root / b.RESULT_PATH).exists())

    def test_cli_default_is_compact_and_existing_output_refuses_before_replay(self):
        spec = importlib.util.spec_from_file_location(
            "pref1_cli_control", ROOT / "scripts/reproduce_fgc_tdg11_msel1_pref1.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _config_raw, _config, _terminal, result = synthetic_bundle()
        with (
            patch.object(b, "verify_compact", return_value=result) as compact,
            patch.object(
                b,
                "bind_raw_result",
                side_effect=AssertionError("unexpected live replay"),
            ),
            patch.object(
                io,
                "publish_exclusive_file",
                side_effect=AssertionError("unexpected write"),
            ),
            redirect_stdout(StringIO()) as output,
        ):
            self.assertEqual(module.main([]), 0)
        compact.assert_called_once()
        self.assertIn('"mode":"compact_check"', output.getvalue())
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            target = root / b.RESULT_PATH
            target.parent.mkdir(parents=True)
            target.write_bytes(b"foreign output must remain")
            with (
                patch.object(module, "ROOT", root),
                patch.object(
                    b,
                    "bind_raw_result",
                    side_effect=AssertionError("unexpected live replay"),
                ),
                self.assertRaises(b.TDG11PREF1Error),
            ):
                module.main(["--bind", "--write-result"])
            self.assertEqual(target.read_bytes(), b"foreign output must remain")
