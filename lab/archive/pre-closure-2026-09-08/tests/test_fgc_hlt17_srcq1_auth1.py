from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import hlt17_srcq1_auth1 as authority  # noqa: E402


CAPS = (
    ("RK4-2049", "0x1.0000000000000p-4", "0x1.33345e70b5188p+0", "0x1.aaa90b0fb5c26p-8"),
    ("RK4-4097", "0x1.0000000000000p-5", "0x1.3333333aecf3ep+0", "0x1.aaaaaa9fefc9cp-9"),
    ("RK4-8193", "0x1.0000000000000p-6", "0x1.3333333333653p+0", "0x1.aaaaaaaaaa654p-10"),
    ("SSPRK3-4097", "0x1.0000000000000p-5", "0x1.3333333fe51c4p+0", "0x1.aaaaaa9908e70p-9"),
    ("SSPRK3-8193", "0x1.0000000000000p-6", "0x1.33333333349e5p+0", "0x1.aaaaaaaaa8b26p-10"),
    ("SSPRK3-16385", "0x1.0000000000000p-7", "0x1.3333333333346p+0", "0x1.aaaaaaaaaaa91p-11"),
)


def _cap_config():
    return {
        "member_cap": [
            {
                "member_key": key,
                "grid_spacing_hex": spacing,
                "previous_speed_upper_hex": speed,
                "cfl_maximum_hex": "0x1.0000000000000p-3",
                "requested_cap_hex": cap,
            }
            for key, spacing, speed, cap in CAPS
        ]
    }


class SRCQ1AuthorityImplementationTests(unittest.TestCase):
    def test_cap_mapping_is_exact_and_rejects_a_swapped_value(self) -> None:
        parsed = authority._cap_mapping(_cap_config())
        self.assertEqual(tuple(parsed), tuple(row[0] for row in CAPS))
        self.assertEqual(parsed["RK4-2049"], CAPS[0][3])
        altered = _cap_config()
        altered["member_cap"][0]["requested_cap_hex"] = CAPS[2][3]
        with self.assertRaisesRegex(authority.HLT17SRCQ1AuthorityError, "formula"):
            authority._cap_mapping(altered)

    def test_environment_binds_executable_numpy_extension_and_kernel(self) -> None:
        observed = authority.environment_identity()
        for name in (
            "executable_sha256", "numpy_extension_sha256", "kernel_release",
            "python_version", "numpy_version", "machine",
        ):
            self.assertTrue(observed[name])
        self.assertEqual(len(observed["executable_sha256"]), 64)
        self.assertEqual(len(observed["numpy_extension_sha256"]), 64)

    def test_source_closure_covers_committed_package_and_runner(self) -> None:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, check=True,
            capture_output=True, text=True,
        ).stdout.strip()
        record = authority.source_closure_identity(ROOT, head)
        self.assertGreater(record["file_count"], 100)
        paths = {item[0] for item in record["files"]}
        self.assertIn(authority.RUNNER_PATH, paths)
        self.assertIn("src/recursive_horizons/fgc/evolution/hlt17_srcq1_auth1.py", paths)
        self.assertEqual(len(record["sha256"]), 64)

    def test_live_validator_refuses_noncommit_and_missing_freeze(self) -> None:
        with self.assertRaises(authority.HLT17SRCQ1AuthorityError):
            authority.validate_authority_delta(
                repository_root=ROOT,
                authority_commit="HEAD",
                implementation_sha256="0" * 64,
            )
        with self.assertRaises(authority.HLT17SRCQ1AuthorityError):
            authority.validate_authority_delta(
                repository_root=ROOT,
                authority_commit="0" * 40,
                implementation_sha256="0" * 64,
            )


if __name__ == "__main__":
    unittest.main()
