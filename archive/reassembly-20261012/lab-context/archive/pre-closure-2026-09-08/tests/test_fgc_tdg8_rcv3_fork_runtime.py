"""Adversarial tests for the external RCV3 recovery-fork installer."""
from __future__ import annotations

import ast
from dataclasses import replace
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import stat
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore
from recursive_horizons.fgc.evolution import tdg8_rcv3_fork_runtime as rcv3


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / rcv3.SOURCE_RELATIVE
AUTHORITY = rcv3.RecoveryForkAuthority(
    authority_commit_sha="a" * 40,
    frz1_result_sha256="b" * 64,
    runtime_sha256="c" * 64,
    bootstrap_script_sha256="d" * 64,
)


def _tree(root: Path) -> tuple[tuple[object, ...], ...]:
    rows: list[tuple[object, ...]] = []
    for directory, names, files in os.walk(root, followlinks=False):
        base = Path(directory)
        for name in sorted((*names, *files)):
            path = base / name
            info = path.lstat()
            relative = path.relative_to(root).as_posix()
            digest = (
                sha256(path.read_bytes()).hexdigest()
                if stat.S_ISREG(info.st_mode)
                else None
            )
            rows.append(
                (
                    relative,
                    stat.S_IFMT(info.st_mode),
                    info.st_size,
                    digest,
                    os.readlink(path) if stat.S_ISLNK(info.st_mode) else None,
                )
            )
    return tuple(sorted(rows))


def _clone_source(repository: Path) -> Path:
    destination = repository / rcv3.SOURCE_RELATIVE
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SOURCE, destination)
    return destination


def _stages(repository: Path) -> tuple[Path, ...]:
    parent = repository / rcv3.WRAPPER_RELATIVE.parent
    return tuple(sorted(parent.glob(".tdg8-rcv3-*"))) if parent.exists() else ()


class TDG8RCV3ForkRuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.actual_before = _tree(SOURCE)
        cls.authenticated = rcv3.authenticate_recovery_fork_source(ROOT)
        cls.actual_after = _tree(SOURCE)

    def test_exact_terminal_source_and_fifty_leaf_lineage_are_read_only(self) -> None:
        source = self.authenticated
        self.assertEqual(self.actual_before, self.actual_after)
        self.assertEqual(len(source.source_leaves), 59)
        self.assertEqual(source.source_tree_sha256, rcv3.SOURCE_TREE_SHA256)
        self.assertEqual(len(source.projected_leaves), 50)
        self.assertEqual(source.projected_tree_sha256, rcv3.PROJECTED_TREE_SHA256)
        self.assertEqual(source.terminal_checkpoint.generation, 10)
        self.assertEqual(
            source.terminal_checkpoint.sha256,
            rcv3.TERMINAL_CHECKPOINT_SHA256,
        )
        self.assertEqual(source.terminal_checkpoint.disposition, "invalid_terminal")
        self.assertEqual(source.anchor_checkpoint.generation, 9)
        self.assertEqual(
            source.anchor_checkpoint.sha256,
            rcv3.ANCHOR_CHECKPOINT_SHA256,
        )
        self.assertEqual(source.anchor_checkpoint.disposition, "nonterminal")
        self.assertEqual(
            set(source.source_leaves) - set(source.projected_leaves),
            {leaf for leaf in source.source_leaves if leaf.path in rcv3.EXCLUDED_SOURCE_PATHS},
        )
        self.assertEqual(
            {leaf.path for leaf in source.source_leaves}
            - {leaf.path for leaf in source.projected_leaves},
            set(rcv3.EXCLUDED_SOURCE_PATHS),
        )

    def test_install_is_exact_restartable_and_idempotent_without_source_mutation(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            source = _clone_source(repository)
            before = _tree(source)
            installed = rcv3.install_recovery_fork(repository, authority=AUTHORITY)
            self.assertTrue(installed.installed_now)
            self.assertEqual(_tree(source), before)
            reopened = rcv3.install_recovery_fork(repository, authority=AUTHORITY)
            self.assertFalse(reopened.installed_now)
            self.assertEqual(reopened.receipt_sha256, installed.receipt_sha256)
            self.assertEqual(_tree(source), before)
            wrapper = repository / rcv3.WRAPPER_RELATIVE
            store = repository / rcv3.PROJECTED_STORE_RELATIVE
            self.assertEqual(
                len([path for path in store.rglob("*") if path.is_file()]),
                rcv3.PROJECTED_LEAF_COUNT,
            )
            self.assertEqual(list((store / "locks").iterdir()), [])
            checkpoint = HLT16CampaignStore(store).authenticated_checkpoint_at_generation(9)
            status = HLT16CampaignStore(store).inspect_recovery()
            self.assertEqual(checkpoint.sha256, rcv3.ANCHOR_CHECKPOINT_SHA256)
            self.assertTrue(status.safe_to_restart)
            self.assertFalse(status.terminal)
            receipt_path = repository / rcv3.RECEIPT_RELATIVE
            receipt = json.loads(receipt_path.read_text("ascii"))
            self.assertEqual(receipt["projection_id"], rcv3.PROJECTION_ID)
            self.assertEqual(
                receipt["kind"], "externally_identified_recovery_fork_projection"
            )
            self.assertEqual(
                receipt["external_install_authority"], AUTHORITY.as_mapping()
            )
            self.assertFalse(
                receipt["internal_hlt16_identity"]["new_internal_campaign_created"]
            )
            self.assertEqual(
                receipt["internal_hlt16_identity"]["campaign_id"],
                rcv3.INTERNAL_CAMPAIGN_ID,
            )
            self.assertEqual(
                sha256(rcv3.canonical({k: v for k, v in receipt.items() if k != "receipt_sha256"})).hexdigest(),
                receipt["receipt_sha256"],
            )
            self.assertEqual(_stages(repository), ())
            self.assertTrue(wrapper.is_dir())

    def test_production_identity_is_fixed_and_cannot_select_another_path(self) -> None:
        mutators = (
            replace(rcv3.PRODUCTION_SPEC, projection_id="foreign"),
            replace(rcv3.PRODUCTION_SPEC, source_relative="runs/foreign"),
            replace(rcv3.PRODUCTION_SPEC, wrapper_relative="runs/other"),
            replace(rcv3.PRODUCTION_SPEC, projected_store_relative="runs/other/calibration"),
            replace(rcv3.PRODUCTION_SPEC, receipt_relative="runs/other/receipt.json"),
        )
        for spec in mutators:
            with self.subTest(spec=spec), self.assertRaisesRegex(
                rcv3.TDG8RCV3ForkRuntimeError, "production projection spec differs"
            ):
                rcv3.authenticate_recovery_fork_source(ROOT, spec)

    def test_install_authority_is_required_and_each_field_is_validated(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            with self.assertRaises(TypeError):
                rcv3.install_recovery_fork(repository)  # type: ignore[call-arg]
            with self.assertRaisesRegex(
                rcv3.TDG8RCV3ForkRuntimeError, "authority type differs"
            ):
                rcv3.install_recovery_fork(
                    repository, authority=object()  # type: ignore[arg-type]
                )

            mutations = (
                replace(AUTHORITY, authority_commit_sha="a" * 39),
                replace(AUTHORITY, authority_commit_sha="A" * 40),
                replace(AUTHORITY, frz1_result_sha256="z" * 64),
                replace(AUTHORITY, runtime_sha256="c" * 63),
                replace(AUTHORITY, bootstrap_script_sha256="D" * 64),
            )
            for authority in mutations:
                with self.subTest(authority=authority), self.assertRaisesRegex(
                    rcv3.TDG8RCV3ForkRuntimeError, "is not hexadecimal"
                ):
                    rcv3.install_recovery_fork(repository, authority=authority)
            self.assertFalse((repository / rcv3.WRAPPER_RELATIVE).exists())

    def test_idempotent_reopen_requires_the_identical_authority_tuple(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            _clone_source(repository)
            installed = rcv3.install_recovery_fork(
                repository, authority=AUTHORITY
            )
            target = repository / rcv3.WRAPPER_RELATIVE
            before = _tree(target)
            different = replace(AUTHORITY, runtime_sha256="e" * 64)
            with self.assertRaisesRegex(
                rcv3.TDG8RCV3ForkRuntimeError, "installed bytes differ"
            ):
                rcv3.install_recovery_fork(repository, authority=different)
            self.assertEqual(_tree(target), before)
            reopened = rcv3.install_recovery_fork(
                repository, authority=AUTHORITY
            )
            self.assertFalse(reopened.installed_now)
            self.assertEqual(reopened.receipt_sha256, installed.receipt_sha256)

    def test_source_symlink_special_foreign_and_partial_trees_fail_closed(self) -> None:
        def symlink_leaf(source: Path) -> None:
            leaf = source / "locks/terminal.lock"
            raw = leaf.read_bytes()
            target = source.parent / "terminal-copy"
            target.write_bytes(raw)
            leaf.unlink()
            os.symlink(target, leaf)

        def symlink_directory(source: Path) -> None:
            real = source / "payloads-real"
            (source / "payloads").rename(real)
            os.symlink(real, source / "payloads")

        def special(source: Path) -> None:
            os.mkfifo(source / "locks/foreign.fifo")

        def foreign(source: Path) -> None:
            (source / "foreign").write_bytes(b"foreign")

        def partial(source: Path) -> None:
            (source / "journal" / (
                "00000000000000000011-"
                f"{rcv3.TERMINAL_JOURNAL_SHA256}.journal"
            )).unlink()

        for name, mutate in (
            ("leaf-symlink", symlink_leaf),
            ("directory-symlink", symlink_directory),
            ("special", special),
            ("foreign", foreign),
            ("partial", partial),
        ):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                repository = Path(directory)
                source = _clone_source(repository)
                mutate(source)
                before = _tree(source)
                with self.assertRaises(rcv3.TDG8RCV3ForkRuntimeError):
                    rcv3.authenticate_recovery_fork_source(repository)
                self.assertEqual(_tree(source), before)
                self.assertFalse((repository / rcv3.WRAPPER_RELATIVE).exists())

    def test_existing_foreign_partial_and_symlink_destinations_are_preserved(self) -> None:
        for name in ("foreign", "partial", "symlink"):
            with self.subTest(name=name), TemporaryDirectory() as directory, TemporaryDirectory() as outside:
                repository = Path(directory)
                _clone_source(repository)
                target = repository / rcv3.WRAPPER_RELATIVE
                target.parent.mkdir(parents=True, exist_ok=True)
                if name == "foreign":
                    target.mkdir()
                    (target / "sentinel").write_bytes(b"foreign")
                elif name == "partial":
                    (target / "calibration/checkpoints").mkdir(parents=True)
                    (target / "partial").write_bytes(b"partial")
                else:
                    os.symlink(outside, target)
                before = _tree(target.parent)
                with self.assertRaises(rcv3.TDG8RCV3ForkRuntimeError):
                    rcv3.install_recovery_fork(repository, authority=AUTHORITY)
                self.assertEqual(_tree(target.parent), before)
                self.assertEqual(_stages(repository), ())

    def test_foreign_destination_race_wins_without_clobber(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            _clone_source(repository)
            target = repository / rcv3.WRAPPER_RELATIVE
            original = rcv3.hlt15._rename_exclusive

            def race(stage: Path, destination: Path) -> None:
                destination.mkdir()
                (destination / "sentinel").write_bytes(b"arrived")
                original(stage, destination)

            with patch.object(rcv3.hlt15, "_rename_exclusive", side_effect=race):
                with self.assertRaisesRegex(
                    rcv3.TDG8RCV3ForkRuntimeError, "destination appeared"
                ):
                    rcv3.install_recovery_fork(repository, authority=AUTHORITY)
            self.assertEqual((target / "sentinel").read_bytes(), b"arrived")
            self.assertEqual(_stages(repository), ())

    def test_exact_concurrent_destination_is_accepted_idempotently(self) -> None:
        with TemporaryDirectory() as template_dir, TemporaryDirectory() as directory:
            template = Path(template_dir)
            _clone_source(template)
            rcv3.install_recovery_fork(template, authority=AUTHORITY)
            exact = template / rcv3.WRAPPER_RELATIVE

            repository = Path(directory)
            _clone_source(repository)
            target = repository / rcv3.WRAPPER_RELATIVE
            original = rcv3.hlt15._rename_exclusive

            def race(stage: Path, destination: Path) -> None:
                shutil.copytree(exact, destination)
                original(stage, destination)

            with patch.object(rcv3.hlt15, "_rename_exclusive", side_effect=race):
                result = rcv3.install_recovery_fork(repository, authority=AUTHORITY)
            self.assertFalse(result.installed_now)
            self.assertEqual(_tree(target), _tree(exact))
            self.assertEqual(_stages(repository), ())

    def test_faults_before_adoption_leave_no_destination_or_private_stage(self) -> None:
        for phase in (
            "after_stage_directories",
            "after_projected_leaves",
            "after_receipt",
            "after_stage_fsync",
            "after_staged_validation",
            "after_source_revalidation",
            "before_exclusive_rename",
        ):
            with self.subTest(phase=phase), TemporaryDirectory() as directory:
                repository = Path(directory)
                source = _clone_source(repository)
                before = _tree(source)

                def fail(observed: str) -> None:
                    if observed == phase:
                        raise RuntimeError(phase)

                with self.assertRaisesRegex(RuntimeError, phase):
                    rcv3.install_recovery_fork(
                        repository, authority=AUTHORITY, fault_hook=fail
                    )
                self.assertEqual(_tree(source), before)
                self.assertFalse((repository / rcv3.WRAPPER_RELATIVE).exists())
                self.assertEqual(_stages(repository), ())

    def test_post_adoption_fault_preserves_exact_idempotent_destination(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            source = _clone_source(repository)
            before = _tree(source)

            def fail(phase: str) -> None:
                if phase == "after_exclusive_rename":
                    raise RuntimeError(phase)

            with self.assertRaisesRegex(RuntimeError, "after_exclusive_rename"):
                rcv3.install_recovery_fork(
                    repository, authority=AUTHORITY, fault_hook=fail
                )
            self.assertTrue((repository / rcv3.WRAPPER_RELATIVE).is_dir())
            self.assertEqual(_tree(source), before)
            resumed = rcv3.install_recovery_fork(repository, authority=AUTHORITY)
            self.assertFalse(resumed.installed_now)
            self.assertEqual(_stages(repository), ())

    def test_stage_symlink_and_special_file_attacks_fail_before_adoption(self) -> None:
        def symlink_attack(stage: Path) -> None:
            payloads = stage / "calibration/payloads"
            real = stage / "calibration/payloads-real"
            payloads.rename(real)
            os.symlink(real, payloads)

        def special_attack(stage: Path) -> None:
            os.mkfifo(stage / "calibration/locks/foreign.fifo")

        for name, attack in (("symlink", symlink_attack), ("special", special_attack)):
            with self.subTest(name=name), TemporaryDirectory() as directory:
                repository = Path(directory)
                _clone_source(repository)

                def inject(phase: str) -> None:
                    if phase == "after_stage_directories":
                        stage = next((repository / rcv3.WRAPPER_RELATIVE.parent).glob(".tdg8-rcv3-*"))
                        attack(stage)

                with self.assertRaises(rcv3.TDG8RCV3ForkRuntimeError):
                    rcv3.install_recovery_fork(
                        repository, authority=AUTHORITY, fault_hook=inject
                    )
                self.assertFalse((repository / rcv3.WRAPPER_RELATIVE).exists())
                self.assertEqual(_stages(repository), ())

    def test_source_change_during_stage_is_detected_before_install(self) -> None:
        with TemporaryDirectory() as directory:
            repository = Path(directory)
            source = _clone_source(repository)

            def mutate(phase: str) -> None:
                if phase == "after_staged_validation":
                    terminal = source / "locks/terminal.lock"
                    terminal.write_bytes(terminal.read_bytes() + b"x")

            with self.assertRaisesRegex(
                rcv3.TDG8RCV3ForkRuntimeError, "source tree identity differs"
            ):
                rcv3.install_recovery_fork(
                    repository, authority=AUTHORITY, fault_hook=mutate
                )
            self.assertFalse((repository / rcv3.WRAPPER_RELATIVE).exists())
            self.assertEqual(_stages(repository), ())

    def test_runtime_imports_no_candidate_branch_or_evolution_operator(self) -> None:
        source = Path(rcv3.__file__).read_text("utf-8")
        tree = ast.parse(source)
        imports: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        joined = "\n".join(imports).lower()
        for forbidden in (
            "sgb", "fgc_qr", "fgcqr", "numerical_engine",
            "progression_attempt", "bounded_retry_runner",
        ):
            self.assertNotIn(forbidden, joined)


if __name__ == "__main__":
    unittest.main()
