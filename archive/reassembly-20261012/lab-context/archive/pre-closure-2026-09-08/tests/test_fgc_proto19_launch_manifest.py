"""Pure manifest-builder tests; no Git, store, solver, or arrays are opened."""
from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path
import tempfile
import unittest

from recursive_horizons.fgc.evolution import proto19_launch_manifest as manifest


ROLES = ("runtime", "config", "result", "documentation", "reproducer", "runner", "test")
MANIFEST_PATH = "configs/fgc/fgc-1-pro19-launch-authority.toml"


def write(root: Path, path: str, value: bytes) -> None:
    target = root / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(value)


def scaffold(*, pending: bool = True) -> tuple[tempfile.TemporaryDirectory[str], Path, bytes]:
    directory = tempfile.TemporaryDirectory()
    root = Path(directory.name)
    entries: list[tuple[str, str]] = []
    for index, role in enumerate(ROLES):
        path = f"tests/{role}-{index}.txt"
        write(root, path, f"{role} {index}\n".encode())
        entries.append((role, path))
    write(root, manifest.AUTH1_PATH, b"{}\n")
    write(root, "results/fgc-1-pro19-frz1.json", b"{}\n")
    mon16_predecessor = "tests/mon16-predecessor.txt"
    mon16_inventory = "tests/mon16-inventory.txt"
    write(root, mon16_predecessor, b"predecessor\n")
    write(root, mon16_inventory, b"inventory\n")
    write(
        root,
        manifest.MON16_CONFIG_PATH,
        (
            'artifact_id = "FGC-1-HLT16-MON16"\n\n'
            '[[predecessors]]\nartifact_id = "fixture"\n'
            f'path = "{mon16_predecessor}"\nsha256 = "{'a' * 64}"\n\n'
            '[[inventory]]\nrole = "source"\n'
            f'path = "{mon16_inventory}"\nsha256 = "{'b' * 64}"\n'
        ).encode(),
    )
    write(root, manifest.MON16_RESULT_PATH, b"{}\n")
    entries.extend((
        ("result", manifest.AUTH1_PATH),
        ("result", "results/fgc-1-pro19-frz1.json"),
        ("config", manifest.MON16_CONFIG_PATH),
        ("result", manifest.MON16_RESULT_PATH),
        ("config", mon16_predecessor),
        ("runtime", mon16_inventory),
    ))
    rows = [f'schema = "{manifest.MANIFEST_SCHEMA}"', "", "[environment]", 'source = "AUTH1"', f'auth1_result_path = "{manifest.AUTH1_PATH}"', ""]
    for role, path in entries:
        digest = manifest.PENDING_SHA256 if pending else sha256((root / path).read_bytes()).hexdigest()
        rows.extend(("[[authority_path]]", f'role = "{role}"', f'path = "{path}"', f'sha256 = "{digest}"', ""))
    return directory, root, "\n".join(rows).encode()


class Proto19LaunchManifestTests(unittest.TestCase):
    def test_refresh_is_deterministic_and_verify_requires_sealed_hashes(self) -> None:
        directory, root, raw = scaffold()
        self.addCleanup(directory.cleanup)
        with self.assertRaises(manifest.Proto19LaunchManifestError):
            manifest.verify_manifest(root, raw, manifest_path=MANIFEST_PATH)
        first = manifest.render_manifest(manifest.refresh_manifest(root, raw, manifest_path=MANIFEST_PATH))
        second = manifest.render_manifest(manifest.refresh_manifest(root, first, manifest_path=MANIFEST_PATH))
        self.assertEqual(first, second)
        verified = manifest.verify_manifest(root, first, manifest_path=MANIFEST_PATH)
        self.assertEqual({row["role"] for row in verified["authority_path"]}, set(ROLES))

    def test_missing_role_duplicate_path_and_self_binding_fail(self) -> None:
        directory, root, raw = scaffold(pending=False)
        self.addCleanup(directory.cleanup)
        text = raw.decode()
        self.assertRaises(manifest.Proto19LaunchManifestError, manifest.parse_manifest, text.replace('role = "runner"', 'role = "runtime"').encode(), manifest_path=MANIFEST_PATH, allow_pending=False)
        self.assertRaises(manifest.Proto19LaunchManifestError, manifest.parse_manifest, text.replace('path = "tests/runner-5.txt"', 'path = "tests/runtime-0.txt"').encode(), manifest_path=MANIFEST_PATH, allow_pending=False)
        self.assertRaises(manifest.Proto19LaunchManifestError, manifest.parse_manifest, text.replace('path = "tests/runner-5.txt"', f'path = "{MANIFEST_PATH}"').encode(), manifest_path=MANIFEST_PATH, allow_pending=False)

    def test_unsafe_path_and_symlink_are_rejected(self) -> None:
        directory, root, raw = scaffold(pending=False)
        self.addCleanup(directory.cleanup)
        text = raw.decode()
        self.assertRaises(manifest.Proto19LaunchManifestError, manifest.parse_manifest, text.replace('path = "tests/runtime-0.txt"', 'path = "../escape.txt"').encode(), manifest_path=MANIFEST_PATH, allow_pending=False)
        target = root / "tests/runtime-0.txt"
        replacement = root / "tests/replacement.txt"
        replacement.write_bytes(target.read_bytes())
        target.unlink()
        os.symlink(replacement, target)
        with self.assertRaises(manifest.Proto19LaunchManifestError):
            manifest.verify_manifest(root, raw, manifest_path=MANIFEST_PATH)

    def test_missing_binding_fails_refresh(self) -> None:
        directory, root, raw = scaffold()
        self.addCleanup(directory.cleanup)
        (root / "tests/runtime-0.txt").unlink()
        with self.assertRaises(manifest.Proto19LaunchManifestError):
            manifest.refresh_manifest(root, raw, manifest_path=MANIFEST_PATH)

    def test_every_mon16_config_path_must_be_bound(self) -> None:
        directory, root, raw = scaffold(pending=False)
        self.addCleanup(directory.cleanup)
        text = raw.decode()
        stanza = (
            '[[authority_path]]\nrole = "config"\n'
            'path = "tests/mon16-predecessor.txt"\n'
            f'sha256 = "{sha256((root / "tests/mon16-predecessor.txt").read_bytes()).hexdigest()}"\n\n'
        )
        self.assertIn(stanza, text)
        with self.assertRaisesRegex(manifest.Proto19LaunchManifestError, "MON16 configured"):
            manifest.verify_manifest(root, text.replace(stanza, "").encode(), manifest_path=MANIFEST_PATH)


if __name__ == "__main__":
    unittest.main()
