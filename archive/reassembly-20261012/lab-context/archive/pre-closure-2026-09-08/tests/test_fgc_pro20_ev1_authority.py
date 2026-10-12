"""PRO20-EV1 exact-delta authority and nonmutating CLI controls."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from recursive_horizons.fgc.evolution import pro20_ev1_authority as authority
from recursive_horizons.fgc.evolution.pro20_ev1_authority import (
    ACTIVE_MAP_PATH,
    ALLOWED_FREEZE_DELTA,
    ARTIFACT_ID as AUTHORITY_ARTIFACT_ID,
    AUTHORITY_PATH,
    CATALOG_BUILDER_PATH,
    CATALOG_CHECKER_PATH,
    CATALOG_INDEX_PATH,
    CATALOG_OUTPUT_PATH,
    CATALOG_TEST_PATH,
    CHANGELOG_PATH,
    CLAIM_LEDGER_PATH,
    CONFIG_PATH,
    HARNESS_PATH,
    MAKE_PATH,
    MAKE_TEST_PATH,
    OWNER_PATH,
    README_PATH,
    ROADMAP_PATH,
    RUNNER_PATH,
    RUNTIME_MATRIX_PATH,
    PRO20EV1AuthorityError,
    authenticate_historical_source_tree,
    authenticate_rec1_pref1_compact,
    authenticate_static_inputs,
    environment_identity,
    expected_config,
    inventory_historical_runs,
    recommended_resource_ceilings,
    render_config,
    require_output_absent,
    validate_authority_delta,
)


ROOT = Path(__file__).resolve().parents[1]
PYTHON = "/opt/homebrew/Caskroom/miniconda/base/bin/python3"
ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"


def _run(root: Path, *arguments: str) -> str:
    completed = subprocess.run(
        arguments,
        cwd=root,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode:
        raise AssertionError(completed.stderr or completed.stdout)
    return completed.stdout.strip()


def _write(root: Path, relative: str, raw: bytes) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _production_tree_identity() -> tuple[tuple[object, ...], ...] | None:
    """Snapshot the canonical tree without assuming an unconsumed authority."""

    container = ROOT / authority.PRODUCTION_CONTAINER
    if not container.exists():
        return None
    paths = [container, *container.rglob("*")]
    return tuple(
        (
            "." if path == container else path.relative_to(container).as_posix(),
            metadata.st_dev,
            metadata.st_ino,
            metadata.st_mode,
            metadata.st_size,
            metadata.st_mtime_ns,
        )
        for path in sorted(paths, key=lambda item: item.as_posix())
        for metadata in (path.lstat(),)
    )


class _AuthorityRepository:
    def __init__(self, root: Path) -> None:
        self.root = root
        _run(root, "git", "init", "-q")
        _run(root, "git", "config", "user.email", "test@example.invalid")
        _run(root, "git", "config", "user.name", "PRO20 Test")
        self.harness = b"# implementation\n"
        self.runner = b"# runner\n"
        self.authority = b"# authority\n"
        _write(root, HARNESS_PATH, self.harness)
        _write(root, RUNNER_PATH, self.runner)
        _write(root, AUTHORITY_PATH, self.authority)
        _run(root, "git", "add", ".")
        _run(root, "git", "commit", "-q", "-m", "implementation")
        self.implementation_commit = _run(root, "git", "rev-parse", "HEAD")

    def freeze(self, *, environment: dict[str, str], source_sha256: str) -> str:
        hashes = {
            "implementation_sha256": _sha(self.harness),
            "runner_sha256": _sha(self.runner),
            "authority_sha256": _sha(self.authority),
        }
        config = expected_config(
            implementation_commit=self.implementation_commit,
            implementation_hashes=hashes,
            source_closure_sha256=source_sha256,
            environment=environment,
        )
        for relative in ALLOWED_FREEZE_DELTA:
            raw = (
                render_config(config)
                if relative == CONFIG_PATH
                else f"freeze:{relative}\n".encode("ascii")
            )
            _write(self.root, relative, raw)
        _run(self.root, "git", "add", ".")
        _run(self.root, "git", "commit", "-q", "-m", "freeze")
        return _run(self.root, "git", "rev-parse", "HEAD")


class PRO20EV1AuthorityTests(unittest.TestCase):
    def test_new_owners_do_not_import_historical_campaign_writers(self) -> None:
        forbidden = (
            "hlt16_campaign_runtime",
            "proto19_launch_authority",
            "proto19_cfl_continuation_authority",
            "run_fgc_gr0_calibration_v",
            "run_fgc_pro19_event1",
        )
        for relative in (HARNESS_PATH, AUTHORITY_PATH, RUNNER_PATH):
            source = (ROOT / relative).read_text(encoding="utf-8")
            with self.subTest(path=relative):
                for marker in forbidden:
                    self.assertNotIn(marker, source)

    def test_allowed_freeze_delta_is_exact_and_publicly_truthful(self) -> None:
        self.assertEqual(
            set(ALLOWED_FREEZE_DELTA),
            {
                CONFIG_PATH,
                OWNER_PATH,
                MAKE_PATH,
                CATALOG_BUILDER_PATH,
                CATALOG_OUTPUT_PATH,
                CATALOG_INDEX_PATH,
                CATALOG_CHECKER_PATH,
                CATALOG_TEST_PATH,
                MAKE_TEST_PATH,
                README_PATH,
                CHANGELOG_PATH,
                ACTIVE_MAP_PATH,
                CLAIM_LEDGER_PATH,
                RUNTIME_MATRIX_PATH,
                ROADMAP_PATH,
            },
        )
        resources = recommended_resource_ceilings()
        self.assertEqual(resources["max_total_wall_seconds"], 86400.0)
        self.assertEqual(resources["max_rss_bytes"], 4294967296)
        self.assertEqual(resources["max_namespace_bytes"], 17179869184)
        self.assertEqual(resources["min_free_disk_bytes"], 34359738368)
        self.assertEqual(resources["max_published_attempts"], 4096)
        self.assertFalse(resources["scientific_threshold_fitted"])

    def test_real_inputs_authenticate_and_output_state_is_fail_closed(self) -> None:
        rec1 = authenticate_rec1_pref1_compact(ROOT)
        self.assertEqual(rec1["artifact_id"], "FGC-1-HLT17-SRCQ1-REC1-PREF1")
        source = authenticate_historical_source_tree(ROOT)
        self.assertEqual(source["origin_capture_sha256"], authority.ORIGIN_CAPTURE_SHA256)
        self.assertEqual(
            source["generation1_bridge_content_id"],
            authority.GENERATION1_BRIDGE_CONTENT_ID,
        )
        self.assertEqual(authenticate_static_inputs(ROOT), dict(authority.STATIC_INPUT_SHA256))
        runs = inventory_historical_runs(ROOT)
        self.assertEqual(runs["digest"], authority.HISTORICAL_RUNS_DIGEST)
        if (ROOT / authority.PRODUCTION_CONTAINER).exists():
            with self.assertRaisesRegex(
                PRO20EV1AuthorityError, "output namespace already exists"
            ):
                require_output_absent(ROOT)
        else:
            require_output_absent(ROOT)

    def test_exact_successor_validation_and_full_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            repository = _AuthorityRepository(root)
            environment = environment_identity()
            source_sha = "6" * 64
            freeze = repository.freeze(
                environment=environment,
                source_sha256=source_sha,
            )
            implementation_sha = _sha(repository.harness)
            source_configurations = dict(authority.SOURCE_CONFIGURATION_SHA256_BY_MEMBER)
            with (
                patch.object(
                    authority,
                    "authenticate_rec1_pref1_compact",
                    return_value={"artifact_id": "REC1"},
                ),
                patch.object(
                    authority,
                    "authenticate_historical_source_tree",
                    return_value={
                        "generation1_bridge_content_id": authority.GENERATION1_BRIDGE_CONTENT_ID
                    },
                ),
                patch.object(
                    authority,
                    "authenticate_static_inputs",
                    return_value=dict(authority.STATIC_INPUT_SHA256),
                ),
                patch.object(authority, "authenticate_planck", return_value={}),
                patch.object(
                    authority,
                    "inventory_historical_runs",
                    return_value={"digest": authority.HISTORICAL_RUNS_DIGEST},
                ),
                patch.object(authority, "environment_identity", return_value=environment),
                patch.object(
                    authority,
                    "source_closure_identity",
                    return_value={"sha256": source_sha, "file_count": 3},
                ),
                patch.object(authority, "require_output_absent"),
                patch.object(authority, "_require_process_preflight"),
                patch.object(
                    authority,
                    "_require_disk_and_rss",
                    return_value={"rss_bytes": 1, "free_disk_bytes": 2**40},
                ),
            ):
                receipt = validate_authority_delta(
                    repository_root=root,
                    authority_commit=freeze,
                    implementation_sha256=implementation_sha,
                )
            self.assertEqual(receipt["implementation_commit"], repository.implementation_commit)
            self.assertEqual(receipt["exact_delta_paths"], list(ALLOWED_FREEZE_DELTA))
            self.assertEqual(receipt["source_configuration"], source_configurations)
            self.assertFalse(receipt["campaign_execution_authorized"])
            self.assertTrue(receipt["pro20_namespace_absent"])

    def test_extra_or_intervening_commit_refuses_exact_delta(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            repository = _AuthorityRepository(root)
            environment = environment_identity()
            repository.freeze(environment=environment, source_sha256="6" * 64)
            _write(root, "unexpected.txt", b"foreign\n")
            _run(root, "git", "add", ".")
            _run(root, "git", "commit", "-q", "-m", "intervening")
            head = _run(root, "git", "rev-parse", "HEAD")
            with self.assertRaisesRegex(PRO20EV1AuthorityError, "delta paths"):
                authority._require_committed_image(root, head)

    def test_output_and_staging_paths_refuse_without_cleanup(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw).resolve()
            parent = root / "runs" / "fgc-2-sf1"
            parent.mkdir(parents=True)
            require_output_absent(root)
            stage = parent / ".pro20-event1.evidence-io-stage-foreign"
            stage.mkdir()
            with self.assertRaisesRegex(PRO20EV1AuthorityError, "staging"):
                require_output_absent(root)
            self.assertTrue(stage.is_dir())

    def test_cli_status_and_missing_authority_are_structured_and_nonmutating(self) -> None:
        production_before = _production_tree_identity()
        commands = (
            ([PYTHON, "-I", "-B", "scripts/run_fgc_pro20_ev1.py"], 0),
            ([PYTHON, "-I", "-B", "scripts/run_fgc_pro20_ev1.py", "--run"], 2),
            (
                [
                    PYTHON,
                    "-I",
                    "-B",
                    "scripts/run_fgc_pro20_ev1.py",
                    "--status",
                    "--authority-commit",
                    "0" * 40,
                ],
                1,
            ),
        )
        for command, expected in commands:
            completed = subprocess.run(
                command,
                cwd=ROOT,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            with self.subTest(command=command):
                self.assertEqual(completed.returncode, expected)
                self.assertEqual(completed.stderr, "")
                payload = json.loads(completed.stdout)
                self.assertFalse(payload.get("physical_source_constructed", False))
                self.assertFalse(payload.get("campaign_execution_authorized", False))
                self.assertEqual(_production_tree_identity(), production_before)

    def test_config_renderer_is_canonical_and_does_not_write(self) -> None:
        tracked = ROOT / CONFIG_PATH
        before = tracked.read_bytes() if tracked.is_file() else None
        environment = environment_identity()
        config = expected_config(
            implementation_commit="a" * 40,
            implementation_hashes={
                "implementation_sha256": "1" * 64,
                "runner_sha256": "2" * 64,
                "authority_sha256": "3" * 64,
            },
            source_closure_sha256="4" * 64,
            environment=environment,
        )
        raw = render_config(config)
        self.assertEqual(authority._parse_config(raw), config)
        if before is None:
            self.assertFalse(tracked.exists())
        else:
            self.assertEqual(tracked.read_bytes(), before)
        self.assertEqual(ARTIFACT_ID, AUTHORITY_ARTIFACT_ID)
        self.assertEqual(ARTIFACT_ID, config["artifact_id"])


if __name__ == "__main__":
    unittest.main()
