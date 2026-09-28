from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
from recursive_horizons.fgc.evolution.proto17_hlt14_inputs import (
    HLT14AuthorityError, HLT14BundleError, HLT14NamespaceForbiddenError,
    LaunchAuthorityTuple, LaunchAuthorityVerifier, VerifiedBundleAdapter,
    canonical_bytes, canonical_sha256, shape_framed_array_content_sha256,
)
from recursive_horizons.fgc.evolution.numerical_engine import array_content_sha256
from recursive_horizons.fgc.evolution.proto17_pure_construction import build_genesis, digest, synthetic_genesis_fixture


NAMES = ("u", "p", "q", "tracer_positions", "tracer_proper_times", "event_proper_times", "event_fields")


def git(root: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=root, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True).stdout.strip()


class HLT14InputsTests(unittest.TestCase):
    def authority_repo(self) -> tuple[Path, LaunchAuthorityTuple]:
        temporary = TemporaryDirectory(); self.addCleanup(temporary.cleanup); root = Path(temporary.name)
        git(root, "init", "-q"); git(root, "config", "user.email", "test@example.invalid"); git(root, "config", "user.name", "test")
        module = root / "fixturepkg"; module.mkdir(); (module / "__init__.py").write_text("")
        sources = []
        for role in ("runtime_module", "adapter_module", "runner"):
            source = module / f"{role}.py"; source.write_text("VALUE = 1\n")
            sources.append((role, source))
        result_path = root / "authority.json"
        raw = synthetic_genesis_fixture()
        roles = {role: sha256(source.read_bytes()).hexdigest() for role, source in sources}
        raw["runtime_module_sha256"] = roles["runtime_module"]; raw["adapter_module_sha256"] = roles["adapter_module"]; raw["runner_sha256"] = roles["runner"]
        for field in ("member_descriptor_set_sha256", "genesis_checkpoint_sha256", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256"):
            raw[field] = "0" * 64
        genesis = build_genesis(raw, _validate_derived=False)["genesis_spec"]
        pinned = [{"path": f"fixturepkg/{role}.py", "sha256": roles[role], "module": f"fixturepkg.{role}", "role": role} for role, _source in sources]
        result = {"artifact_id": "FGC-1-HLT14-MON14", "genesis_spec": genesis, "pinned_sources": pinned}
        result_path.write_bytes(canonical_bytes(result)); git(root, "add", "."); git(root, "commit", "-qm", "authority")
        commit = git(root, "rev-parse", "HEAD"); sys.path.insert(0, str(root)); self.addCleanup(lambda: sys.path.remove(str(root)))
        return root, LaunchAuthorityTuple(commit, "authority.json", sha256(result_path.read_bytes()).hexdigest(), canonical_sha256(genesis))

    def test_authority_reads_exact_committed_blob_and_pinned_origin(self):
        root, launch = self.authority_repo(); verified = LaunchAuthorityVerifier(root).verify(launch)
        self.assertEqual(len(verified.genesis_spec["member_descriptors"]), 6); self.assertEqual(len(verified.pinned_sources), 3)

    def test_authority_mutations_fail_closed(self):
        root, launch = self.authority_repo(); verifier = LaunchAuthorityVerifier(root)
        for altered in (
            LaunchAuthorityTuple("f" * 40, launch.authorization_result_path, launch.authorization_result_sha256, launch.genesis_spec_sha256),
            LaunchAuthorityTuple(launch.authorization_commit, launch.authorization_result_path, "0" * 64, launch.genesis_spec_sha256),
            LaunchAuthorityTuple(launch.authorization_commit, launch.authorization_result_path, launch.authorization_result_sha256, "0" * 64),
        ):
            with self.assertRaises(HLT14AuthorityError): verifier.verify(altered)
        with self.assertRaises(HLT14AuthorityError):
            LaunchAuthorityTuple(launch.authorization_commit, "../authority.json", launch.authorization_result_sha256, launch.genesis_spec_sha256)
        (root / "fixturepkg/runtime_module.py").write_text("VALUE = 2\n")
        with self.assertRaises(HLT14AuthorityError): verifier.verify(launch)

    def test_self_referential_authority_payload_fails_closed(self):
        root, launch = self.authority_repo()
        raw = json.loads((root / "authority.json").read_text())
        raw["authorization_commit"] = launch.authorization_commit
        (root / "authority.json").write_bytes(canonical_bytes(raw))
        git(root, "add", "authority.json"); git(root, "commit", "-qm", "self reference")
        altered = LaunchAuthorityTuple(git(root, "rev-parse", "HEAD"), "authority.json", sha256((root / "authority.json").read_bytes()).hexdigest(), launch.genesis_spec_sha256)
        with self.assertRaises(HLT14AuthorityError): LaunchAuthorityVerifier(root).verify(altered)

    def bundle_descriptor(self, root: Path) -> tuple[Path, dict[str, object], dict[str, object]]:
        fixture = synthetic_genesis_fixture()
        raw_descriptor = fixture["member_descriptors"][0]
        point_count = raw_descriptor["point_count"]
        arrays = {"u": np.arange(point_count * 6, dtype="<f8").reshape(point_count, 6), "p": np.ones((point_count, 6), dtype="<f8"), "q": np.zeros((point_count, 6), dtype="<f8"), "tracer_positions": np.zeros(48, dtype="<f8"), "tracer_proper_times": np.ones(48, dtype="<f8"), "event_proper_times": np.zeros((24, 48), dtype="<f8"), "event_fields": np.zeros((24, 48, 6), dtype="<f8")}
        hashes = {name: sha256(np.ascontiguousarray(array).tobytes()).hexdigest() for name, array in arrays.items()}
        physical = array_content_sha256(arrays["u"], arrays["p"], arrays["q"]); restart = array_content_sha256(*(arrays[name] for name in NAMES))
        self.assertEqual(physical, shape_framed_array_content_sha256(arrays["u"], arrays["p"], arrays["q"]))
        self.assertEqual(restart, shape_framed_array_content_sha256(*(arrays[name] for name in NAMES)))
        state = raw_descriptor["state_object"]
        for item in state["physical_arrays"]:
            item["little_endian_c_bytes_sha256"] = hashes[item["logical_name"]]
        state["physical_state_sha256"] = physical; state["restart_payload_sha256"] = restart
        raw_descriptor["state_object_canonical_sha256"] = digest(state)
        meta = {"member_key": raw_descriptor["member_key"], "method": raw_descriptor["method"], "point_count": raw_descriptor["point_count"], "coordinate_time": state["coordinate_time"], "accepted_boundary_time": state["accepted_boundary_time"], "step_index": state["step_index"], "transaction_serial": state["transaction_serial"], "input_hash": state["input_hash"], "source_retry_count": state["source_retry_count"], "CFL_retry_count": state["CFL_retry_count"], "runtime_monitor_state": state["runtime_monitor_state"], "causal_state": state["causal_state"], "tracer_state": state["tracer_state"], "event_history_state": state["event_history_state"], "initial_TDG6_ledger": raw_descriptor["initial_TDG6_ledger"], "initial_TDG6_ledger_sha256": raw_descriptor["initial_TDG6_ledger_sha256"]}
        stored = {item["npz_storage_key"]: arrays[item["logical_name"]] for item in state["physical_arrays"]}
        path = root / "bundle.npz"; np.savez(path, **stored, metadata_utf8=np.frombuffer(canonical_bytes(meta), dtype=np.uint8))
        bundle = fixture["physical_input_manifest"]["bundle_members"][0]
        bundle["relative_input_path"] = path.name; bundle["raw_file_sha256"] = sha256(path.read_bytes()).hexdigest()
        fixture["physical_input_manifest_sha256"] = digest(fixture["physical_input_manifest"])
        for field in ("member_descriptor_set_sha256", "genesis_checkpoint_sha256", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256"):
            fixture[field] = "0" * 64
        genesis = build_genesis(fixture, _validate_derived=False)["genesis_spec"]
        descriptor = genesis["member_descriptors"][0]
        entry = genesis["physical_input_manifest"]["bundle_members"][0]
        return path, descriptor, entry

    def test_bundle_verifies_all_bytes_and_rich_identity(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary); path, descriptor, entry = self.bundle_descriptor(root)
            verified = VerifiedBundleAdapter(root).verify(descriptor, entry)
            self.assertEqual(verified.physical_state_sha256, descriptor["state_object"]["physical_state_sha256"])
            self.assertEqual(len(verified.array_hashes), 7)

    def test_manifest_root_is_bound_and_production_root_is_refused(self):
        with TemporaryDirectory() as temporary:
            authority_root = Path(temporary)
            fixture = synthetic_genesis_fixture()
            fixture["physical_input_manifest"]["relative_input_root"] = "inputs"
            adapter = VerifiedBundleAdapter.from_manifest_root(authority_root, fixture)
            self.assertEqual(adapter.input_root, (authority_root / "inputs").resolve())
            fixture["physical_input_manifest"]["relative_input_root"] = "runs/fgc-2-sf1/proto17/calibration"
            with self.assertRaises(HLT14NamespaceForbiddenError):
                VerifiedBundleAdapter.from_manifest_root(authority_root, fixture)

    def test_bundle_mutations_and_namespace_access_fail_closed(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary); path, descriptor, entry = self.bundle_descriptor(root); adapter = VerifiedBundleAdapter(root)
            for mutate in (
                lambda d, b: b.__setitem__("raw_file_sha256", "0" * 64),
                lambda d, b: d["state_object"]["physical_arrays"][0].__setitem__("shape", [3, 6]),
                lambda d, b: d["state_object"]["physical_arrays"][0].__setitem__("little_endian_c_bytes_sha256", "0" * 64),
                lambda d, b: d["state_object"].__setitem__("physical_state_sha256", "0" * 64),
                lambda d, b: b["required_npz_member_keys"].append("extra"),
            ):
                altered_descriptor, altered_entry = deepcopy(descriptor), deepcopy(entry); mutate(altered_descriptor, altered_entry)
                with self.assertRaises(HLT14BundleError): adapter.verify(altered_descriptor, altered_entry)
            with self.assertRaises(HLT14NamespaceForbiddenError): adapter.forbid_production_namespace(Path("runs/fgc-2-sf1/proto17/calibration"))
            # Duplicate and path-bearing members are rejected before NumPy load.
            bad = root / "bad.npz"
            with zipfile.ZipFile(bad, "w") as archive:
                for key in entry["required_npz_member_keys"]:
                    archive.writestr(f"dir/{key}.npy", b"x")
            altered_descriptor, altered_entry = deepcopy(descriptor), deepcopy(entry); altered_entry["relative_input_path"] = bad.name; altered_entry["raw_file_sha256"] = sha256(bad.read_bytes()).hexdigest()
            with self.assertRaises(HLT14BundleError): adapter.verify(altered_descriptor, altered_entry)


if __name__ == "__main__":
    unittest.main()
