from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from types import CodeType, FunctionType
import sys
import unittest
from unittest.mock import patch

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import gr0_calibration  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_imp1_bridge import (  # noqa: E402
    synthetic_callable_binding_digest,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    UniformRadialGrid,
)
from recursive_horizons.fgc.evolution.pro20_source_binding import (  # noqa: E402
    SCHEMA,
    Pro20GR0SourceBinding,
    Pro20SourceBindingError,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    Proto12GR0EvolutionOperator,
)


def _cell(value: object):
    def inner():
        return value
    assert inner.__closure__ is not None
    return inner.__closure__[0]


def _state() -> EvolutionState:
    grid = UniformRadialGrid(0.0, 8.0, 9)
    u = np.zeros((9, 6), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _projector(reference: EvolutionState, rows: int = 2) -> FunctionType:
    codes = tuple(
        item for item in gr0_calibration.make_gr0_center_boundary_projector.__code__.co_consts
        if type(item) is CodeType and item.co_name == "projector"
    )
    assert len(codes) == 1
    values = {"fixed_outer_rows": rows, "reference": reference}
    return FunctionType(
        codes[0],
        gr0_calibration.__dict__,
        name="projector",
        closure=tuple(_cell(values[name]) for name in codes[0].co_freevars),
    )


def _operator() -> Proto12GR0EvolutionOperator:
    return Proto12GR0EvolutionOperator(
        UniformRadialGrid(0.0, 8.0, 9),
        spatial_order=4,
        ko_dissipation=1.0 / 64.0,
        raw_tolerance=1.0e-12,
        kinetic_condition_maximum=1.0e8,
        maximum_refinement_iterations=16,
        point_batch_size=9,
    )


def _binding() -> Pro20GR0SourceBinding:
    return Pro20GR0SourceBinding(_operator(), _projector(_state()))


class Pro20SourceBindingTests(unittest.TestCase):
    def test_capture_is_outcome_blind_and_includes_omitted_solver_controls(self) -> None:
        binding = _binding()
        record = binding.as_mapping()
        self.assertEqual(record["schema"], SCHEMA)
        self.assertEqual(record["operator"]["maximum_refinement_iterations"], 16)
        self.assertEqual(record["operator"]["point_batch_size"], 9)
        self.assertEqual(record["operator"]["spatial_order"], 4)
        self.assertEqual(record["projector"]["fixed_outer_rows"], 2)
        for label in (
            "physical_source_authenticated",
            "source_manifest_authenticated",
            "historical_origin_authenticated",
            "campaign_execution_authorized",
        ):
            self.assertIs(record[label], False)

    def test_construction_and_validation_do_not_call_the_source(self) -> None:
        with patch(
            "recursive_horizons.fgc.evolution.pro20_source_binding._evaluate_rhs",
            side_effect=AssertionError("source called"),
        ):
            binding = _binding()
            binding.validate()
            binding.hlt17_synthetic_callable_binding()

    def test_returned_mapping_and_bytes_do_not_mutate_the_binding(self) -> None:
        binding = _binding()
        expected = binding.captured_bytes
        record = binding.as_mapping()
        record["operator"]["point_batch_size"] = 1
        self.assertEqual(binding.captured_bytes, expected)
        self.assertNotEqual(binding.as_mapping(), record)
        with self.assertRaises(FrozenInstanceError):
            binding._captured_bytes = b"forged"  # type: ignore[misc]
        with self.assertRaises((TypeError, ValueError)):
            replace(binding, _captured_bytes=b"forged")

    def test_every_mutable_operator_setting_is_rechecked(self) -> None:
        changes = {
            "ko_dissipation": 0.5,
            "raw_tolerance": 0.5,
            "kinetic_condition_maximum": 2.0,
            "maximum_refinement_iterations": 3,
            "point_batch_size": 3,
            "grid": UniformRadialGrid(0.0, 9.0, 9),
            "derivative": gr0_calibration.SBPFirstDerivative(
                UniformRadialGrid(0.0, 8.0, 9), 2
            ),
        }
        for name, value in changes.items():
            with self.subTest(name=name):
                binding = _binding()
                setattr(binding.operator, name, value)
                with self.assertRaisesRegex(Pro20SourceBindingError, "changed"):
                    binding.validate()

    def test_missing_extra_bool_and_oversized_operator_configuration_refuse(self) -> None:
        operator = _operator()
        del operator.point_batch_size
        with self.assertRaisesRegex(Pro20SourceBindingError, "missing or extra"):
            Pro20GR0SourceBinding(operator, _projector(_state()))
        operator = _operator()
        operator.foreign_backend = object()
        with self.assertRaisesRegex(Pro20SourceBindingError, "missing or extra"):
            Pro20GR0SourceBinding(operator, _projector(_state()))
        operator = _operator()
        operator.maximum_refinement_iterations = True
        with self.assertRaisesRegex(Pro20SourceBindingError, "integer"):
            Pro20GR0SourceBinding(operator, _projector(_state()))
        operator = _operator()
        operator.grid = UniformRadialGrid(0.0, 8.0, 16386)
        with self.assertRaisesRegex(Pro20SourceBindingError, "point_count"):
            Pro20GR0SourceBinding(operator, _projector(_state()))

    def test_projector_closure_reference_and_global_identity_are_bound(self) -> None:
        binding = _binding()
        assert binding.projector.__closure__ is not None
        cells = dict(zip(binding.projector.__code__.co_freevars, binding.projector.__closure__, strict=True))
        cells["fixed_outer_rows"].cell_contents = 3
        with self.assertRaisesRegex(Pro20SourceBindingError, "changed"):
            binding.validate()
        foreign_globals = dict(gr0_calibration.__dict__)
        foreign = FunctionType(
            binding.projector.__code__,
            foreign_globals,
            name="projector",
            closure=binding.projector.__closure__,
        )
        with self.assertRaisesRegex(Pro20SourceBindingError, "global owner"):
            Pro20GR0SourceBinding(_operator(), foreign)

    def test_projector_reference_array_drift_refuses(self) -> None:
        binding = _binding()
        assert binding.projector.__closure__ is not None
        cells = dict(zip(binding.projector.__code__.co_freevars, binding.projector.__closure__, strict=True))
        reference = cells["reference"].cell_contents
        object.__setattr__(reference, "u", np.ones_like(reference.u))
        with self.assertRaisesRegex(Pro20SourceBindingError, "changed"):
            binding.validate()

    def test_hlt17_seam_records_local_identity_without_source_authority(self) -> None:
        binding = _binding()
        material = binding.hlt17_synthetic_callable_binding()
        self.assertEqual(material["configuration_sha256"], binding.configuration_sha256)
        self.assertIs(material["physical_source_authenticated"], False)
        digest = synthetic_callable_binding_digest(binding.rhs, binding.project)
        self.assertEqual(len(digest), 64)
        self.assertNotEqual(digest, binding.configuration_sha256)

    def test_rhs_forwarding_is_guarded_and_input_state_is_not_changed(self) -> None:
        binding = _binding()
        state = _state()
        before = (state.u.copy(), state.p.copy(), state.q.copy())
        zeros = np.zeros_like(state.u)
        expected = EvolutionRHS(zeros, zeros, zeros, {"synthetic_source_control": True})
        calls: list[tuple[object, object, object]] = []

        def fake(operator, time, current):
            calls.append((operator, time, current))
            return expected

        with patch(
            "recursive_horizons.fgc.evolution.pro20_source_binding._evaluate_rhs",
            side_effect=fake,
        ):
            result = binding.rhs(0.25, state)
        self.assertIs(result, expected)
        self.assertEqual(calls, [(binding.operator, 0.25, state)])
        for observed, original in zip((state.u, state.p, state.q), before, strict=True):
            np.testing.assert_array_equal(observed, original)

    def test_mutation_during_rhs_is_refused_after_delegation(self) -> None:
        binding = _binding()
        state = _state()
        zeros = np.zeros_like(state.u)

        def mutate(operator, _time, _state):
            operator.point_batch_size = 3
            return EvolutionRHS(zeros, zeros, zeros, {})

        with patch(
            "recursive_horizons.fgc.evolution.pro20_source_binding._evaluate_rhs",
            side_effect=mutate,
        ), self.assertRaisesRegex(Pro20SourceBindingError, "changed"):
            binding.rhs(0.25, state)

    def test_projector_forwarding_matches_preserved_function_and_is_pure(self) -> None:
        binding = _binding()
        state = _state()
        altered = EvolutionState(state.u + 0.125, state.p, state.q)
        before = (altered.u.copy(), altered.p.copy(), altered.q.copy())
        expected = binding.projector(0.25, altered)
        observed = binding.project(0.25, altered)
        for left, right in zip((observed.u, observed.p, observed.q),
                               (expected.u, expected.p, expected.q), strict=True):
            np.testing.assert_array_equal(left, right)
        for current, original in zip((altered.u, altered.p, altered.q), before, strict=True):
            np.testing.assert_array_equal(current, original)

    def test_foreign_input_shape_refuses_before_source_call(self) -> None:
        binding = _binding()
        foreign = EvolutionState(np.zeros((10, 6)), np.zeros((10, 6)), np.zeros((10, 6)))
        with patch(
            "recursive_horizons.fgc.evolution.pro20_source_binding._evaluate_rhs",
            side_effect=AssertionError("source called"),
        ), self.assertRaisesRegex(Pro20SourceBindingError, "shape"):
            binding.rhs(0.25, foreign)


if __name__ == "__main__":
    unittest.main()
