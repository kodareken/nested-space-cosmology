"""Focused prospective-authority controls for TDG9 LOC1."""

from __future__ import annotations

from pathlib import Path
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_loc1_authority as authority  # noqa: E402


class LOC1AuthorityTests(unittest.TestCase):
    def test_frozen_counts_selection_and_nonclaims(self) -> None:
        self.assertEqual(authority.FAILED_OCCURRENCES, ((3, "u:alpha"), (3, "u:R"), (4, "u:alpha"), (4, "u:lambda"), (4, "u:R"), (5, "u:alpha"), (5, "u:v"), (5, "u:lambda"), (5, "u:R"), (5, "q:R")))
        self.assertEqual(authority.POLYNOMIALS_PER_OCCURRENCE, 12_264)
        self.assertEqual(authority.TOTAL_POLYNOMIALS, 122_640)
        self.assertEqual(authority.COMPONENT_CUBIC_COUNT, 122_640)
        self.assertEqual(authority.EXECUTED_COMPONENT_CUBIC_COUNT, 367_920)
        self.assertEqual(authority.COMPLETE_C_MAXIMUM_CANDIDATES, 490_560)
        self.assertEqual(authority.VALUE_V_MAXIMUM_CANDIDATES, 490_560)
        self.assertEqual(authority.SLOPE_S_MAXIMUM_CANDIDATES, 490_560)
        self.assertEqual(authority.TOTAL_VSC_MAXIMUM_CANDIDATES, 1_471_680)
        self.assertEqual(authority.GRID_SPACING_HEX, "0x1.0000000000000p-4")
        self.assertEqual(authority.OUTER_RADIUS_HEX, "0x1.0000000000000p+7")
        self.assertFalse(authority.SCOPE["ULP_used_as_tolerance"])
        self.assertFalse(authority.SCOPE["fourth_width_authorized"])
        self.assertFalse(authority.SCOPE["stage2_source_decomposition_authorized"])
        self.assertFalse(authority.SCOPE["SSPRK3_comparator_authorized"])

    def test_config_is_exact(self) -> None:
        parsed = authority.parse_config((ROOT / authority.CONFIG_PATH).read_bytes())
        self.assertEqual(parsed["work_budget"]["base_cubic_count"], 122_640)  # type: ignore[index]
        self.assertEqual(parsed["work_budget"]["component_polynomials_per_base_cubic"], 3)  # type: ignore[index]
        self.assertEqual(parsed["work_budget"]["executed_component_cubic_count"], 367_920)  # type: ignore[index]
        self.assertEqual(parsed["committed_image"]["base_commit"], authority.BASE_COMMIT)  # type: ignore[index]
        self.assertEqual(tuple(parsed["committed_image"]["delta_paths"]), authority.DELTA_PATHS)  # type: ignore[index]

    def test_prospective_image_requires_exact_direct_delta_and_absence(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        with (
            patch.object(authority, "_environment", return_value=authority.ENVIRONMENT),
            patch.object(authority, "_head", return_value=authority.BASE_COMMIT),
            patch.object(authority, "_working", return_value=(authority.DELTA_PATHS, (), (".qdrant-initialized",))),
            patch.object(authority, "_require_predecessor"),
            patch.object(authority, "_require_absent"),
        ):
            result = authority.build_prelaunch(config, ROOT)
        self.assertEqual(result["gate_status"], "pass")
        self.assertFalse(result["artifact_payload"]["diagnostic_executed"])  # type: ignore[index]

    def test_wrong_delta_and_existing_namespace_fail_closed(self) -> None:
        config = (ROOT / authority.CONFIG_PATH).read_bytes()
        with (
            patch.object(authority, "_environment", return_value=authority.ENVIRONMENT),
            patch.object(authority, "_head", return_value=authority.BASE_COMMIT),
            patch.object(authority, "_working", return_value=(authority.DELTA_PATHS[:-1], (), (".qdrant-initialized",))),
            self.assertRaisesRegex(authority.LOC1AuthorityError, "delta differs"),
        ):
            authority.build_prelaunch(config, ROOT)

    def test_absence_rejects_existing_and_broken_symlink_parents(self) -> None:
        for broken in (False, True):
            with self.subTest(broken=broken), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                (root / "runs").mkdir()
                target = root / "missing" if broken else root / "foreign"
                if not broken:
                    target.mkdir()
                os.symlink(target, root / "runs" / "fgc-2-sf1")
                with self.assertRaisesRegex(authority.LOC1AuthorityError, "output parent is unsafe"):
                    authority._require_absent(root)  # type: ignore[attr-defined]


if __name__ == "__main__":
    unittest.main()
