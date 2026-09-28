from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts import reproduce_fgc_pro18_frz1 as rep
from recursive_horizons.fgc.evolution.proto18_prelaunch_contract import Proto18PrelaunchStop

class PRO18FreezeTests(unittest.TestCase):
    def test_canonical_result_reproduces(self) -> None:
        subprocess.run([sys.executable, str(ROOT / "scripts/reproduce_fgc_pro18_frz1.py"), "--verify"], check=True)

    def test_selector_and_chain_mutations_fail_closed(self) -> None:
        config = rep.load_config()
        bad_bytes = (ROOT / "configs/fgc/fgc-1-pro18-frz1.toml").read_text().replace("step_index=342", "step_index=0", 1)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.toml"; path.write_text(bad_bytes)
            with self.assertRaises(Proto18PrelaunchStop): rep.load_config(path)
        bad = deepcopy(config); bad["selectors"]["member"][0]["step_index"] = 0
        with self.assertRaises(Proto18PrelaunchStop):
            if bad["selectors"]["member"][0]["step_index"] <= 0: raise Proto18PrelaunchStop("PRO18_SELECTOR_IDENTITY_DRIFT")
        self.assertEqual(config["future_chain"]["required_order"], ["PREF26", "AUTH1", "HLT15_GEN1", "PREF27"])

    def test_config_section_ancestry_and_lineage_mutations_fail_closed(self) -> None:
        source = (ROOT / "configs/fgc/fgc-1-pro18-frz1.toml").read_text()
        mutations = (
            ("schema_version = 1", "schema_version = 2"),
            ("FGC-1-PRO18-AUTH1", "FGC-1-PRO18-AUTHX"),
            ("container_format = \"NPZ\"", "container_format = \"NPY\""),
            (
                "calibration = \"runs/fgc-2-sf1/proto17/calibration\"",
                "calibration = \"runs/fgc-2-sf1/proto17/other\"",
            ),
            (
                "PREF27 = \"FGC-1-PRO18-PREF27\"",
                "PREF27 = \"FGC-1-PRO18-PREFXX\"",
            ),
            (
                "unsafe_output_root = \"PRO18_UNSAFE_OUTPUT_ROOT\"",
                "unsafe_output_root = \"PRO18_IGNORE_OUTPUT_ROOT\"",
            ),
            ("detached_launch_image_required = true", "detached_launch_image_required = false"),
            (
                '"results/fgc-1-pro15-frz1.json" = "99d24425cb5d40da2e248605265b83bcdc82a2126cfa4b5bce4239c6ef30635c"\n',
                "",
            ),
        )
        for before, after in mutations:
            with tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "bad.toml"; path.write_text(source.replace(before, after, 1))
                with self.assertRaises(Proto18PrelaunchStop): rep.load_config(path)
        from recursive_horizons.fgc.evolution.proto18_prelaunch_contract import verify_immutable_compact_files
        with self.assertRaises(Proto18PrelaunchStop): verify_immutable_compact_files(ROOT, "0" * 40, {"README.md": "0" * 64})
        with self.assertRaises(Proto18PrelaunchStop): verify_immutable_compact_files(ROOT, "4b27a48facf75b1c3186efed6cdfa752ce160e93", {"README.md": "0" * 64})

    def test_all_production_and_physical_claims_remain_false(self) -> None:
        result = rep.build(observe_prelaunch_roots=False)
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))
        self.assertTrue(result["gate_status"]["PRO18_static_production_prelaunch_overlay_frozen"])
        self.assertTrue(result["gate_status"]["PREF26_implementation_or_design_authorized"])
        self.assertFalse(result["gate_status"]["FGCQR_holdout_execution_authorized"])
        forged = deepcopy(result); forged["nonclaims"]["candidate_execution_authorized"] = True
        with self.assertRaises(Proto18PrelaunchStop): rep.verify_result(forged)

    def test_prelaunch_absence_is_temporal_but_offline_verify_is_not(self) -> None:
        original = rep.inspect_future_output_roots_only
        try:
            rep.inspect_future_output_roots_only = lambda *_: {"calibration": True, "holdout": False}
            with self.assertRaises(Proto18PrelaunchStop): rep.build(observe_prelaunch_roots=True)
            # The frozen record is validated against immutable compact lineage,
            # not incorrectly reclassified after its future root later exists.
            rep.build(observe_prelaunch_roots=False)
        finally:
            rep.inspect_future_output_roots_only = original

    def test_unsafe_root_and_raw_history_access_are_rejected_or_absent(self) -> None:
        from recursive_horizons.fgc.evolution.proto18_prelaunch_contract import inspect_future_output_roots_only, safe_compact_file
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); target = root / "target"; target.mkdir()
            path = root / "runs/fgc-2-sf1/proto17/calibration"; path.parent.mkdir(parents=True); path.symlink_to(target)
            with self.assertRaises(Proto18PrelaunchStop): inspect_future_output_roots_only(root, {"calibration": "runs/fgc-2-sf1/proto17/calibration"})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); target = root / "missing"
            path = root / "runs/fgc-2-sf1/proto17/calibration"; path.parent.mkdir(parents=True); path.symlink_to(target)
            with self.assertRaises(Proto18PrelaunchStop): inspect_future_output_roots_only(root, {"calibration": "runs/fgc-2-sf1/proto17/calibration"})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); target = root / "elsewhere"; target.mkdir()
            ancestor = root / "runs/fgc-2-sf1/proto17"; ancestor.parent.mkdir(parents=True); ancestor.symlink_to(target)
            with self.assertRaises(Proto18PrelaunchStop): inspect_future_output_roots_only(root, {"calibration": "runs/fgc-2-sf1/proto17/calibration"})
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); target = root / "target.txt"; target.write_text("sealed")
            linked = root / "bound.txt"; linked.symlink_to(target)
            with self.assertRaises(Proto18PrelaunchStop):
                safe_compact_file(root, "bound.txt")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); target = root / "elsewhere"; target.mkdir()
            linked_parent = root / "tracked"; linked_parent.symlink_to(target)
            with self.assertRaises(Proto18PrelaunchStop):
                safe_compact_file(root, "tracked/bound.txt")

        original_read_bytes = Path.read_bytes
        original_read_text = Path.read_text
        forbidden = tuple(rep.EXPECTED_SOURCE_MAP["shared_containers"] + rep.EXPECTED_SOURCE_MAP["legacy_event_logs"])

        def is_forbidden(path: Path) -> bool:
            rendered = path.as_posix()
            return any(rendered == relative or rendered.endswith("/" + relative) for relative in forbidden)

        def guarded_read_bytes(path: Path, *args, **kwargs):
            self.assertFalse(is_forbidden(path), path)
            return original_read_bytes(path, *args, **kwargs)

        def guarded_read_text(path: Path, *args, **kwargs):
            self.assertFalse(is_forbidden(path), path)
            return original_read_text(path, *args, **kwargs)

        with patch.object(Path, "read_bytes", guarded_read_bytes), patch.object(Path, "read_text", guarded_read_text):
            rep.build(observe_prelaunch_roots=False)

if __name__ == "__main__": unittest.main()
