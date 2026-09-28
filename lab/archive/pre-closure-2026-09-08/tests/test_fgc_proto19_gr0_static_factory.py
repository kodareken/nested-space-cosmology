"""Focused equivalence and fail-closed tests for the direct GR-0 factory."""

from __future__ import annotations

import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import (  # noqa: E402
    proto19_gr0_static_factory as factory,
)
from recursive_horizons.fgc.evolution.hlt16_member_codec import (  # noqa: E402
    _digest,
    _template_identity,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto14_runtime import (  # noqa: E402
    Proto14RunMember,
)
from recursive_horizons.fgc.evolution.protocol_v17 import MEMBER_KEYS  # noqa: E402


GENERATION8_CHECKPOINT = ROOT / (
    "runs/fgc-2-sf1/proto17/calibration/checkpoints/"
    "00000000000000000008-"
    "6d2b6e37e7c7cfd4ce64504353d886099003f959ea2aaee1a14fc778ce14fa46.json"
)


def _state_arrays(state):
    return state.u, state.p, state.q


def _tracer_arrays(tracers):
    return (
        tracers.labels,
        tracers.positions,
        tracers.proper_times,
        *tracers.event_proper_times,
        *tracers.event_fields,
    )


def _captured_static_inputs() -> dict[str, bytes]:
    return {
        relative: (ROOT / relative).read_bytes()
        for relative in factory.STATIC_INPUT_PATHS
    }


def _copy_static_inputs(root: Path) -> None:
    for relative in factory.STATIC_INPUT_PATHS:
        source = ROOT / relative
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)


class Proto19GR0StaticFactoryImportTests(unittest.TestCase):
    def test_source_has_no_executable_or_candidate_only_import(self) -> None:
        path = ROOT / (
            "src/recursive_horizons/fgc/evolution/proto19_gr0_static_factory.py"
        )
        tree = ast.parse(path.read_text("utf-8"))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imported.append(node.module or "")
        normalized = [
            name.lower().replace("_", "").replace("-", "") for name in imported
        ]
        self.assertFalse(
            any(name == "scripts" or name.startswith("scripts.") for name in imported)
        )
        self.assertFalse(
            any(marker in name for name in normalized for marker in ("fgcqr", "sgbl"))
        )

    def test_fresh_process_build_loads_no_candidate_only_runtime(self) -> None:
        program = """
import json
import importlib
import sys
import types
from pathlib import Path
root = Path.cwd()
for name, path in (
    ('recursive_horizons', root / 'src' / 'recursive_horizons'),
    ('recursive_horizons.fgc', root / 'src' / 'recursive_horizons' / 'fgc'),
    ('recursive_horizons.fgc.evolution', root / 'src' / 'recursive_horizons' / 'fgc' / 'evolution'),
):
    package = types.ModuleType(name)
    package.__package__ = name
    package.__path__ = [str(path)]
    sys.modules[name] = package
build_static_gr0_shells = importlib.import_module(
    'recursive_horizons.fgc.evolution.proto19_gr0_static_factory'
).build_static_gr0_shells
build_static_gr0_shells(Path.cwd())
forbidden = (
    'recursive_horizons.fgc.evolution.initial_state_bridge',
    'recursive_horizons.fgc.evolution.nonlinear_source',
    'recursive_horizons.fgc.evolution.static_initial_admission',
    'scripts.run_fgc_gr0_calibration_v13',
)
print(json.dumps({
    'forbidden_loaded': sorted(name for name in forbidden if name in sys.modules),
    'shared_definitions_loaded': all(name in sys.modules for name in (
        'recursive_horizons.fgc.evolution.vectorized_source',
        'recursive_horizons.fgc.initial_data_family',
    )),
}, sort_keys=True))
"""
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(ROOT / "src")
        completed = subprocess.run(
            [sys.executable, "-c", program],
            cwd=ROOT,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(
            json.loads(completed.stdout),
            {"forbidden_loaded": [], "shared_definitions_loaded": True},
        )


class Proto19GR0StaticFactoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.direct = factory.build_static_gr0_shells(ROOT)
        # This import is the historical oracle only.  The isolated test above
        # proves that production construction does not load this path.
        from scripts.run_fgc_gr0_calibration_v13 import (
            _build_frozen_member_shells,
        )

        cls.legacy = {
            key: Proto14RunMember.from_proto7(member)
            for key, member in _build_frozen_member_shells().items()
        }

    def test_exact_six_member_legacy_runtime_oracle(self) -> None:
        self.assertEqual(tuple(self.direct), MEMBER_KEYS)
        self.assertEqual(tuple(self.legacy), MEMBER_KEYS)
        for key in MEMBER_KEYS:
            with self.subTest(member=key):
                direct = self.direct[key]
                legacy = self.legacy[key]
                self.assertIs(type(direct), Proto14RunMember)
                self.assertEqual(_template_identity(direct), _template_identity(legacy))
                for actual, expected in zip(
                    _state_arrays(direct.state),
                    _state_arrays(legacy.state),
                    strict=True,
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())

                direct_before = array_content_sha256(*_state_arrays(direct.state))
                legacy_before = array_content_sha256(*_state_arrays(legacy.state))
                direct_rhs = direct.operator(0.0, direct.state)
                legacy_rhs = legacy.operator(0.0, legacy.state)
                self.assertEqual(
                    direct_before,
                    array_content_sha256(*_state_arrays(direct.state)),
                )
                self.assertEqual(
                    legacy_before,
                    array_content_sha256(*_state_arrays(legacy.state)),
                )
                for actual, expected in zip(
                    (direct_rhs.du, direct_rhs.dp, direct_rhs.dq),
                    (legacy_rhs.du, legacy_rhs.dp, legacy_rhs.dq),
                    strict=True,
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())
                self.assertEqual(direct_rhs.diagnostics, legacy_rhs.diagnostics)

                direct_projected = direct.projector(0.0, direct.state)
                legacy_projected = legacy.projector(0.0, legacy.state)
                for actual, expected in zip(
                    _state_arrays(direct_projected),
                    _state_arrays(legacy_projected),
                    strict=True,
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())
                for actual, expected in zip(
                    _tracer_arrays(direct.tracers),
                    _tracer_arrays(legacy.tracers),
                    strict=True,
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())
                direct_preview = direct.tracers.preview_advance(
                    old_state=direct.state,
                    new_state=direct.state,
                    coordinates=direct.initial.grid.coordinates,
                    step_size=2.0**-20,
                )
                legacy_preview = legacy.tracers.preview_advance(
                    old_state=legacy.state,
                    new_state=legacy.state,
                    coordinates=legacy.initial.grid.coordinates,
                    step_size=2.0**-20,
                )
                for actual, expected in zip(
                    direct_preview, legacy_preview, strict=True
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())

                self.assertEqual(
                    direct.transaction.thresholds, legacy.transaction.thresholds
                )
                self.assertEqual(
                    direct.transaction.causal_state,
                    legacy.transaction.causal_state,
                )
                self.assertEqual(
                    direct.transaction.boundary_geometry,
                    legacy.transaction.boundary_geometry,
                )
                self.assertEqual(direct.transaction.state, legacy.transaction.state)
                self.assertEqual(direct.temporal_ledger, legacy.temporal_ledger)

    def test_tracer_uses_descriptor_compatible_identity_without_module_import(
        self,
    ) -> None:
        self.assertEqual(
            factory.NormalFlowTracers.__module__,
            "scripts.run_fgc_gr0_calibration",
        )
        for member in self.direct.values():
            self.assertIs(type(member.tracers), factory.NormalFlowTracers)

    def test_captured_bytes_are_filesystem_independent_and_equivalent(self) -> None:
        with patch.object(
            factory,
            "_read_static_inputs",
            side_effect=AssertionError("captured construction touched the filesystem"),
        ):
            captured = factory.build_static_gr0_shells(
                ROOT / "path-that-does-not-exist",
                static_input_bytes=_captured_static_inputs(),
            )
        self.assertEqual(tuple(captured), MEMBER_KEYS)
        for key in MEMBER_KEYS:
            with self.subTest(member=key):
                self.assertEqual(
                    _template_identity(captured[key]),
                    _template_identity(self.direct[key]),
                )
                for actual, expected in zip(
                    _state_arrays(captured[key].state),
                    _state_arrays(self.direct[key].state),
                    strict=True,
                ):
                    self.assertEqual(actual.tobytes(), expected.tobytes())

    def test_dormant_candidate_entrypoints_are_not_invoked(self) -> None:
        from recursive_horizons.fgc import initial_data_family
        from recursive_horizons.fgc.evolution import vectorized_source

        refusal = AssertionError("dormant candidate entrypoint was invoked")
        with (
            patch.object(vectorized_source, "_coerce_action", side_effect=refusal),
            patch.object(vectorized_source, "batch_ref1_residual", side_effect=refusal),
            patch.object(
                vectorized_source, "solve_grid_accelerations", side_effect=refusal
            ),
            patch.object(
                initial_data_family, "solve_initial_data", side_effect=refusal
            ),
        ):
            captured = factory.build_static_gr0_shells(
                ROOT / "path-that-does-not-exist",
                static_input_bytes=_captured_static_inputs(),
            )
            member = captured["RK4-2049"]
            before = array_content_sha256(*_state_arrays(member.state))
            rhs = member.operator(0.0, member.state)
            self.assertTrue(rhs.diagnostics["source_raw_gate_passed"])
            self.assertEqual(before, array_content_sha256(*_state_arrays(member.state)))

    def test_captured_byte_inventory_and_identity_fail_closed(self) -> None:
        captured = _captured_static_inputs()
        missing = dict(captured)
        missing.pop(factory.STATIC_INPUT_PATHS[-1])
        extra = {**captured, "results/not-authorized.json": b"{}"}
        mutated = dict(captured)
        mutated[factory.STATIC_INPUT_PATHS[0]] += b"\n"
        cases = (
            ("missing", missing, "captured static GR-0 input inventory differs"),
            ("extra", extra, "captured static GR-0 input inventory differs"),
            ("mutated", mutated, "static GR-0 input identity differs"),
        )
        for label, value, message in cases:
            with (
                self.subTest(case=label),
                self.assertRaisesRegex(
                    factory.Proto19GR0StaticFactoryError,
                    message,
                ),
            ):
                factory.build_static_gr0_shells(
                    ROOT / "path-that-does-not-exist",
                    static_input_bytes=value,
                )

    def test_templates_match_the_persisted_generation8_descriptors(self) -> None:
        checkpoint = json.loads(GENERATION8_CHECKPOINT.read_bytes())
        self.assertEqual(checkpoint["generation"], 8)
        self.assertEqual(set(checkpoint["members"]), set(MEMBER_KEYS))
        for key in MEMBER_KEYS:
            with self.subTest(member=key):
                descriptor_sha256 = checkpoint["members"][key]["descriptor_sha256"]
                descriptor_path = ROOT / (
                    "runs/fgc-2-sf1/proto17/calibration/states/"
                    f"{descriptor_sha256}.json"
                )
                descriptor = json.loads(descriptor_path.read_bytes())
                persisted = descriptor["metadata"]
                current = _template_identity(self.direct[key])
                self.assertEqual(current, persisted["template"])
                self.assertEqual(_digest(current), persisted["template_sha256"])
                self.assertEqual(
                    current["tracer_type"],
                    "scripts.run_fgc_gr0_calibration.NormalFlowTracers",
                )

    def test_every_static_input_mutation_fails_before_construction(self) -> None:
        for relative in factory.STATIC_INPUT_PATHS:
            with self.subTest(input=relative), TemporaryDirectory() as directory:
                root = Path(directory)
                _copy_static_inputs(root)
                target = root / relative
                target.write_bytes(target.read_bytes() + b"\n")
                with self.assertRaisesRegex(
                    factory.Proto19GR0StaticFactoryError,
                    "static GR-0 input identity differs",
                ):
                    factory.build_static_gr0_shells(root)

    def test_filesystem_reader_rejects_symlink_hardlink_and_nonregular_leaf(
        self,
    ) -> None:
        relative = factory.STATIC_INPUT_PATHS[0]
        cases = ("symlink", "hardlink", "directory")
        for case in cases:
            with self.subTest(case=case), TemporaryDirectory() as directory:
                root = Path(directory)
                _copy_static_inputs(root)
                target = root / relative
                if case == "symlink":
                    exact_copy = root / "exact-symlink-target"
                    shutil.copy2(target, exact_copy)
                    target.unlink()
                    target.symlink_to(exact_copy)
                    message = "unavailable or unsafe"
                elif case == "hardlink":
                    exact_copy = root / "exact-hardlink-target"
                    shutil.copy2(target, exact_copy)
                    target.unlink()
                    os.link(exact_copy, target)
                    self.assertEqual(target.stat().st_nlink, 2)
                    message = "single-link regular file"
                else:
                    target.unlink()
                    target.mkdir()
                    message = "single-link regular file"
                with self.assertRaisesRegex(
                    factory.Proto19GR0StaticFactoryError,
                    message,
                ):
                    factory.build_static_gr0_shells(root)

    def test_filesystem_reader_rejects_lexical_path_escape(self) -> None:
        with TemporaryDirectory() as directory:
            with self.assertRaisesRegex(
                factory.Proto19GR0StaticFactoryError,
                "path escapes its root",
            ):
                factory._read_regular_file_no_follow(
                    Path(directory), "../outside-static-input"
                )

    def test_constructed_state_drift_fails_closed(self) -> None:
        original = factory.project_gr0_reference_balanced_state

        def changed(initial, *, spatial_order):
            state = original(initial, spatial_order=spatial_order)
            u = state.u.copy()
            u[-1, -1] = np.nextafter(u[-1, -1], np.inf)
            return type(state)(u, state.p, state.q)

        with patch.object(factory, "project_gr0_reference_balanced_state", changed):
            with self.assertRaisesRegex(
                factory.Proto19GR0StaticFactoryError,
                "constructed GR-0 state differs",
            ):
                factory.build_static_gr0_shells(ROOT)


if __name__ == "__main__":
    unittest.main()
