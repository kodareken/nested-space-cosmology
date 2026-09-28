from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto18_pref27_binder import (  # noqa: E402
    Proto18Pref27Error,
    build_pref27,
    validate_pref27_config,
    verify_generation_zero_store,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import canonical, digest  # noqa: E402


CANONICAL_NAMESPACE = "runs/fgc-2-sf1/proto17/calibration"
CANONICAL_STORE = ROOT / CANONICAL_NAMESPACE


def immutable_manifest(path: Path) -> tuple[tuple[str, str, int], ...]:
    """A read-only type-and-byte inventory; never follows a symlink."""
    root_info = path.lstat()
    if os.path.islink(path):
        return ((".", "symlink", root_info.st_mode),)
    if path.is_file():
        return ((".", sha256(path.read_bytes()).hexdigest(), root_info.st_mode),)
    if not path.is_dir():
        return ((".", "other", root_info.st_mode),)
    rows: list[tuple[str, str, int]] = []
    pending = [path]
    while pending:
        parent = pending.pop()
        for entry in sorted(parent.iterdir(), reverse=True):
            info = entry.lstat()
            relative = entry.relative_to(path).as_posix()
            if os.path.islink(entry):
                rows.append((relative, "symlink", info.st_mode))
            elif entry.is_file():
                rows.append((relative, sha256(entry.read_bytes()).hexdigest(), info.st_mode))
            elif entry.is_dir():
                rows.append((relative, "directory", info.st_mode))
                pending.append(entry)
            else:
                rows.append((relative, "other", info.st_mode))
    return tuple(rows)


class Pref27BinderTests(unittest.TestCase):
    """All adversarial stores live in disposable local clones, never GEN1."""

    @classmethod
    def setUpClass(cls) -> None:
        if not CANONICAL_STORE.is_dir():  # pragma: no cover - concrete production prerequisite
            raise unittest.SkipTest("HLT15 canonical generation-zero store is absent")

    def clone_with_store(self, temporary: Path) -> Path:
        clone = temporary / "repository"
        subprocess.run(["git", "clone", "-q", "--no-hardlinks", str(ROOT), str(clone)], check=True)
        # The receipt's authority commit is the checked-out repository's exact
        # AUTH1 predecessor; preserve it in the temporary Git object database.
        subprocess.run(["git", "checkout", "-q", "--detach", "2e3373ae84e6f69ffc4339096ad163a014f3a985"], cwd=clone, check=True)
        destination = clone / CANONICAL_NAMESPACE
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(CANONICAL_STORE, destination, symlinks=True)
        return clone

    def assert_rejected_unchanged(self, root: Path, *, stop_id: str | None = None) -> None:
        store = root / CANONICAL_NAMESPACE
        before = immutable_manifest(store) if store.exists() or store.is_symlink() else ()
        with self.assertRaises(Proto18Pref27Error) as raised:
            verify_generation_zero_store(root)
        if stop_id is not None:
            self.assertEqual(raised.exception.stop_id, stop_id)
        after = immutable_manifest(store) if store.exists() or store.is_symlink() else ()
        self.assertEqual(after, before)

    def test_actual_canonical_store_is_exact_and_build_is_nonmutating(self) -> None:
        before = immutable_manifest(CANONICAL_STORE)
        evidence = verify_generation_zero_store(ROOT)
        self.assertEqual(len(evidence.leaves), 8)
        self.assertEqual(evidence.checkpoint["campaign_generation"], 0)
        self.assertFalse(evidence.checkpoint["terminal_lock"])
        self.assertFalse(evidence.receipt["claims"]["state_advanced"])
        self.assertFalse(evidence.receipt["claims"]["candidate_execution_authorized"])
        self.assertEqual(immutable_manifest(CANONICAL_STORE), before)

        artifact = build_pref27(ROOT)
        self.assertEqual(artifact["artifact_id"], "FGC-1-PRO18-PREF27")
        claims = artifact["artifact_payload"]["claims"]
        self.assertFalse(claims["state_advanced"])
        self.assertFalse(claims["candidate_execution_authorized"])
        self.assertFalse(claims["physical_transition_claim_authorized"])
        self.assertEqual(immutable_manifest(CANONICAL_STORE), before)

    def test_byte_identical_copy_is_admitted_but_is_not_a_reusable_result(self) -> None:
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            evidence = verify_generation_zero_store(clone)
            self.assertEqual(evidence.tree_sha256, verify_generation_zero_store(ROOT).tree_sha256)
            # PREF27's temporary controls must reject this exact copy as an
            # occupied namespace rather than treating it as resumable state.
            artifact = build_pref27(clone)
            payload = artifact["artifact_payload"]
            self.assertTrue(payload["claims"]["production_namespace_reuse_refusal_verified"])
            controls = payload["temporary_refusal_controls"]
            self.assertTrue(controls["byte_identical_actual_store_copy"])
            self.assertTrue(controls["byte_identical_actual_store_copy_no_staging_residue"])

    def test_arriving_foreign_root_wins_exclusive_rename_without_clobber(self) -> None:
        """The decisive race is a foreign store arriving after staging, not two twins."""
        before = immutable_manifest(CANONICAL_STORE)
        artifact = build_pref27(ROOT)
        controls = artifact["artifact_payload"]["temporary_refusal_controls"]
        self.assertTrue(controls["arriving_foreign_root_at_exclusive_rename"])
        self.assertTrue(controls["arriving_foreign_root_no_staging_residue"])
        self.assertTrue(controls["internally_self_consistent_foreign_store_no_staging_residue"])
        self.assertEqual(immutable_manifest(CANONICAL_STORE), before)

    def test_authority_receipt_checkpoint_state_tree_and_claim_mutations_fail_closed(self) -> None:
        mutators = {
            "authority": lambda clone: (clone / "results/fgc-1-pro18-auth1.json").write_bytes(b"{}"),
            "receipt": lambda clone: (clone / CANONICAL_NAMESPACE / "receipts/generation-zero.json").write_bytes(b"{}"),
            "checkpoint": lambda clone: next((clone / CANONICAL_NAMESPACE / "checkpoints").glob("*.json")).write_bytes(b"{}"),
            "state": lambda clone: next((clone / CANONICAL_NAMESPACE / "states").glob("*.json")).write_bytes(b"{}"),
            "extra_tree_leaf": lambda clone: (clone / CANONICAL_NAMESPACE / "unexpected").write_text("foreign", encoding="utf-8"),
        }
        for name, mutate in mutators.items():
            with self.subTest(name=name), TemporaryDirectory() as directory:
                clone = self.clone_with_store(Path(directory))
                mutate(clone)
                self.assert_rejected_unchanged(clone)

    def test_semantically_promoted_receipt_and_authority_are_rejected_after_rehash(self) -> None:
        for name, mutate, stop_id in (
            (
                "authority",
                lambda receipt: receipt["external_authority"].__setitem__(
                    "authorization_commit", "0" * 40
                ),
                "PREF27_AUTHORITY_GIT_INVALID",
            ),
            (
                "claim",
                lambda receipt: receipt["claims"].__setitem__("state_advanced", True),
                "PREF27_RECEIPT_DRIFT",
            ),
        ):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                clone = self.clone_with_store(Path(directory))
                receipt_path = clone / CANONICAL_NAMESPACE / "receipts/generation-zero.json"
                receipt = json.loads(receipt_path.read_text("utf-8"))
                mutate(receipt)
                bare = dict(receipt)
                bare.pop("receipt_sha256")
                receipt["receipt_sha256"] = digest(bare)
                receipt_path.write_bytes(canonical(receipt))
                self.assert_rejected_unchanged(clone, stop_id=stop_id)

    def test_content_address_filename_and_empty_directory_mutations_fail_closed(self) -> None:
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            state = next((clone / CANONICAL_NAMESPACE / "states").glob("*.json"))
            state.rename(state.with_name("0" * 64 + ".json"))
            self.assert_rejected_unchanged(clone, stop_id="PREF27_STORE_PATH_UNSAFE")
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            (clone / CANONICAL_NAMESPACE / "empty").mkdir()
            self.assert_rejected_unchanged(clone, stop_id="PREF27_STORE_TREE_DRIFT")

    def test_config_scope_and_claim_promotions_fail_before_store_access(self) -> None:
        config = tomllib.loads((ROOT / "configs/fgc/fgc-1-pro18-pref27.toml").read_text("utf-8"))
        for mutate in (
            lambda value: value["scope"].__setitem__("state_advance", True),
            lambda value: value["claims"].__setitem__("candidate_execution_authorized", True),
        ):
            with self.subTest(mutate=mutate):
                candidate = json.loads(json.dumps(config))
                mutate(candidate)
                with self.assertRaises(Proto18Pref27Error) as raised:
                    validate_pref27_config(candidate)
                self.assertEqual(raised.exception.stop_id, "PREF27_CONTRACT_DRIFT")

    def test_staging_residue_is_rejected_even_when_the_eight_store_leaves_are_valid(self) -> None:
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            (clone / CANONICAL_NAMESPACE).parent.joinpath(
                ".calibration.hlt15-stage-leaked"
            ).mkdir()
            with self.assertRaises(Proto18Pref27Error) as raised:
                build_pref27(clone)
            self.assertEqual(raised.exception.stop_id, "PREF27_STORE_TREE_DRIFT")

    def test_partial_file_and_symlink_store_shapes_fail_closed_without_adoption(self) -> None:
        for name, populate in (
            ("partial", lambda store: (store.mkdir(), (store / "partial").write_text("x", encoding="utf-8"))),
            ("file", lambda store: store.write_text("foreign", encoding="utf-8")),
            ("symlink", lambda store: os.symlink("/tmp", store)),
        ):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                clone = self.clone_with_store(Path(directory))
                shutil.rmtree(clone / CANONICAL_NAMESPACE)
                store = clone / CANONICAL_NAMESPACE
                populate(store)
                self.assert_rejected_unchanged(clone)

    def test_repository_root_nested_directory_and_leaf_symlinks_fail_closed(self) -> None:
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            linked_root = Path(directory) / "repository-link"
            os.symlink(clone, linked_root)
            with self.assertRaises(Proto18Pref27Error) as raised:
                verify_generation_zero_store(linked_root)
            self.assertEqual(raised.exception.stop_id, "PREF27_STORE_PATH_UNSAFE")

        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            states = clone / CANONICAL_NAMESPACE / "states"
            real_states = states.with_name("states-real")
            states.rename(real_states)
            os.symlink(real_states.name, states)
            self.assert_rejected_unchanged(clone, stop_id="PREF27_STORE_PATH_UNSAFE")

        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            state = next((clone / CANONICAL_NAMESPACE / "states").glob("*.json"))
            real_state = state.with_suffix(".real")
            state.rename(real_state)
            os.symlink(real_state.name, state)
            self.assert_rejected_unchanged(clone, stop_id="PREF27_STORE_PATH_UNSAFE")

    def test_special_store_entry_fails_closed_when_supported(self) -> None:
        if not hasattr(os, "mkfifo"):
            self.skipTest("platform has no portable FIFO constructor")
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            special = clone / CANONICAL_NAMESPACE / "special"
            try:
                os.mkfifo(special)
            except OSError as error:
                self.skipTest(f"FIFO unsupported in this temporary filesystem: {error}")
            self.assert_rejected_unchanged(clone, stop_id="PREF27_STORE_PATH_UNSAFE")

    def test_partial_hlt15_store_availability_is_not_silently_accepted(self) -> None:
        """PREF27 authenticates GEN1 leaves; it does not reopen old raw archives."""
        with TemporaryDirectory() as directory:
            clone = self.clone_with_store(Path(directory))
            missing = next((clone / CANONICAL_NAMESPACE / "states").glob("*.json"))
            missing.unlink()
            with self.assertRaises(Proto18Pref27Error):
                build_pref27(clone)


if __name__ == "__main__":
    unittest.main()
