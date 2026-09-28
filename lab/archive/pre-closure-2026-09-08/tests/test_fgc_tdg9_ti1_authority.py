"""Focused prospective authority controls for TDG9 TI1."""

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]

from recursive_horizons.fgc.evolution import tdg9_ti1_authority as authority  # noqa: E402


class TI1AuthorityTests(unittest.TestCase):
    def test_config_freezes_one_variable_semantics_and_nonclaims(self) -> None:
        parsed = authority.parse_config(
            (ROOT / authority.CONFIG_PATH).read_bytes()
        )
        semantics = parsed["semantics"]
        self.assertEqual(semantics["tableau_selector"], "SSPRK3")
        self.assertEqual(
            semantics["actual_spatial_operator"], "inherited_RK4_2049_SBP4"
        )
        self.assertFalse(semantics["production_SSPRK3_comparator"])
        self.assertFalse(semantics["independent_method_agreement"])
        self.assertFalse(parsed["scope"]["PDE_state_commit_authorized"])
        self.assertEqual(parsed["work_budget"]["maximum_shadow_proposals"], 21)
        self.assertEqual(
            parsed["work_budget"]["maximum_stage_and_endpoint_RHS_records"], 84
        )

    def test_wrong_tableau_state_width_channel_or_threshold_is_rejected(self) -> None:
        raw = (ROOT / authority.CONFIG_PATH).read_bytes()
        mutations = (
            (b'tableau_selector = "SSPRK3"', b'tableau_selector = "RK4"'),
            (authority.PHYSICAL_STATE_SHA256.encode(), b"0" * 64),
            (b"0x1.aaa9612df8000p-11", b"0x1.aaa9612df0000p-11"),
            (b'channel = "q:R"', b'channel = "p:R"'),
            (b'minimum_observed_order = "3/2"', b'minimum_observed_order = "1"'),
        )
        for old, new in mutations:
            with self.subTest(old=old):
                changed = raw.replace(old, new, 1)
                self.assertNotEqual(changed, raw)
                with self.assertRaises(authority.TI1AuthorityError):
                    authority.parse_config(changed)

    def test_output_absence_rejects_existing_and_symlinked_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / authority.OUTPUT_NAMESPACE
            target.mkdir(parents=True)
            with self.assertRaisesRegex(
                authority.TI1AuthorityError, "output namespace already exists"
            ):
                authority.require_output_absent(root)
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaisesRegex(
                authority.TI1AuthorityError, "output parent is unsafe"
            ):
                authority.require_output_absent(root)

    def test_compact_rejects_claim_mutation_and_duplicate_key(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        raw = (ROOT / authority.RESULT_PATH).read_bytes()
        value = authority.validate_compact(config, raw)
        value["artifact_payload"]["claims"]["physical_result_earned"] = True
        with self.assertRaises(authority.TI1AuthorityError):
            authority.validate_compact(config, authority.canonical_pretty(value))
        duplicate = raw.replace(b'{\n  "artifact_id"', b'{\n  "artifact_id": "x",\n  "artifact_id"', 1)
        with self.assertRaises(authority.TI1AuthorityError):
            authority.validate_compact(config, duplicate)

    def test_delta_excludes_qdrant_and_raw_namespace(self) -> None:
        self.assertNotIn(".qdrant-initialized", authority.DELTA_PATHS)
        self.assertFalse(
            any(path.startswith("runs/") for path in authority.DELTA_PATHS)
        )
        self.assertNotIn(authority.OUTPUT_NAMESPACE, authority.DELTA_PATHS)


if __name__ == "__main__":
    unittest.main()
