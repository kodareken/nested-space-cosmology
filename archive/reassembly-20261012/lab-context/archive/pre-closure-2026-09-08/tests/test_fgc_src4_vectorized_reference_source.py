from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import struct
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    exact_affine_system,
    exact_gaussian_solve,
    exact_ref1_residual,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    _inverse_and_derivative,
    _reference_covariant_connection_difference,
    _reference_covariant_ricci,
    diagnose_gr0_reference_balanced_accelerations,
    gr0_reference_balanced_ref1_residual_batch,
)
from recursive_horizons.fgc.evolution.gr0_direct_source import (  # noqa: E402
    _metric_jet_from_adm,
)
from recursive_horizons.fgc.evolution.src4_vectorized_reference_source import (  # noqa: E402
    _vectorized_reference_covariant_connection_difference,
    _vectorized_reference_covariant_ricci,
    diagnose_gr0_vectorized_reference_accelerations,
    gr0_vectorized_reference_ref1_residual_batch,
)


POINT_FIXTURE = json.loads(
    (REPOSITORY / "configs/fgc/fgc-1-src2-pref11-point0.json").read_text(
        encoding="utf-8"
    )
)
CONTROL_FIXTURE = json.loads(
    (REPOSITORY / "configs/fgc/fgc-1-src3-controls.json").read_text(
        encoding="utf-8"
    )
)
RAW_FIXTURE = (
    REPOSITORY
    / "runs/fgc-2-sf1/src2/affine-wall-capture/affine-wall-fixture.npz"
)


def _hex_row(values: list[str]) -> np.ndarray:
    return np.asarray([float.fromhex(value) for value in values], dtype=np.float64)


def _compact_lower_jet() -> tuple[np.ndarray, ...]:
    lower = POINT_FIXTURE["lower_jet"]
    return tuple(_hex_row(lower[name]) for name in ("u", "p", "q", "p_r", "q_r"))


def _ordered_float_bits(value: float) -> int:
    signed = struct.unpack(">q", struct.pack(">d", value))[0]
    return 0x8000000000000000 - signed if signed < 0 else signed


def _ulp_distance(left: float, right: float) -> int:
    if left == right:
        return 0
    return abs(_ordered_float_bits(left) - _ordered_float_bits(right))


class FGCSRC4VectorizedReferenceSourceTests(unittest.TestCase):
    def test_exact_spherical_minkowski_remains_bitwise_zero(self) -> None:
        radii = np.asarray(
            [
                float.fromhex(value)
                for value in CONTROL_FIXTURE["exact_reference_radii_binary64_hex"]
            ],
            dtype=np.float64,
        )
        points = radii.size
        u = np.zeros((points, 6), dtype=np.float64)
        p = np.zeros_like(u)
        q = np.zeros_like(u)
        p_r = np.zeros_like(u)
        q_r = np.zeros_like(u)
        u[:, 0] = 1.0
        u[:, 2] = 1.0
        u[:, 3] = radii
        q[:, 3] = 1.0
        result = gr0_vectorized_reference_ref1_residual_batch(
            u,
            p,
            q,
            np.zeros((1, points, 6), dtype=np.float64),
            p_r,
            q_r,
            radii,
        )
        for value in (
            result.full_residual,
            result.unredefined_metric_residual,
            result.scalar_residual,
            result.gauge_constraint,
            result.hamiltonian_constraint,
            result.momentum_constraint,
            result.ricci_scalar,
            result.ricci_squared,
        ):
            self.assertTrue(np.array_equal(value, np.zeros_like(value)))
        solved = diagnose_gr0_vectorized_reference_accelerations(
            u, p, q, p_r, q_r, radii
        )
        self.assertTrue(solved.raw_gate_passed)
        self.assertEqual(solved.residual_infinity, 0.0)
        self.assertTrue(
            np.array_equal(solved.accelerations, np.zeros_like(solved.accelerations))
        )

    def test_vector_contractions_match_the_independent_src3_index_loops(self) -> None:
        controls = CONTROL_FIXTURE["independent_nontrivial_controls"]
        lower = {
            name: np.concatenate(
                [
                    _hex_row(control["lower_jet"][name])[None, :]
                    for control in controls
                ],
                axis=0,
            )
            for name in ("u", "p", "q", "p_r", "q_r")
        }
        radii = np.asarray(
            [float.fromhex(control["radius_binary64_hex"]) for control in controls],
            dtype=np.float64,
        )
        seeds = np.zeros((7, len(controls), 6), dtype=np.float64)
        for field in range(6):
            seeds[field + 1, :, field] = 1.0
        metric, first, second = _metric_jet_from_adm(
            lower["u"],
            lower["p"],
            lower["q"],
            seeds,
            lower["p_r"],
            lower["q_r"],
        )
        inverse, inverse_derivative = _inverse_and_derivative(metric, first)
        scalar = _reference_covariant_connection_difference(
            metric,
            first,
            second,
            inverse,
            inverse_derivative,
            radii,
        )
        vector = _vectorized_reference_covariant_connection_difference(
            metric,
            first,
            second,
            inverse,
            inverse_derivative,
            radii,
        )
        self.assertTrue(np.array_equal(scalar[0], vector[0]))
        self.assertTrue(np.array_equal(scalar[1], vector[1]))
        self.assertLessEqual(float(np.max(np.abs(scalar[2] - vector[2]))), 1.0e-15)
        self.assertLessEqual(float(np.max(np.abs(scalar[3] - vector[3]))), 2.0e-14)
        scalar_ricci = _reference_covariant_ricci(scalar[0], scalar[2], scalar[3])
        vector_ricci = _vectorized_reference_covariant_ricci(
            vector[0], vector[2], vector[3]
        )
        scale = max(1.0, float(np.max(np.abs(scalar_ricci))))
        self.assertLessEqual(
            float(np.max(np.abs(scalar_ricci - vector_ricci))) / scale,
            64.0 * np.finfo(np.float64).eps,
        )

    def test_pref11_point_matches_exact_oracle_and_src3_root(self) -> None:
        lower = _compact_lower_jet()
        radius = float.fromhex(POINT_FIXTURE["coordinate_radius"])
        exact_constant, exact_matrix = exact_affine_system(*lower, radius)
        exact_root = exact_gaussian_solve(
            exact_matrix,
            tuple(-value for value in exact_constant),
        )
        rounded = np.asarray([float(value) for value in exact_root], dtype=np.float64)
        solved = diagnose_gr0_vectorized_reference_accelerations(
            *(value[None, :] for value in lower),
            np.asarray((radius,), dtype=np.float64),
        )
        src3 = diagnose_gr0_reference_balanced_accelerations(
            *(value[None, :] for value in lower),
            np.asarray((radius,), dtype=np.float64),
        )
        exact_residual = exact_ref1_residual(*lower, radius, solved.accelerations[0])
        self.assertTrue(solved.raw_gate_passed)
        self.assertLess(solved.residual_infinity, 1.0e-24)
        self.assertLess(
            float(max((abs(value) for value in exact_residual), default=Fraction(0))),
            1.0e-24,
        )
        self.assertLessEqual(
            max(
                _ulp_distance(left, right)
                for left, right in zip(solved.accelerations[0], rounded, strict=True)
            ),
            16,
        )
        self.assertLessEqual(
            float(np.max(np.abs(solved.accelerations - src3.accelerations))),
            2.0e-15,
        )

    def test_independent_controls_match_exact_oracle_and_src3(self) -> None:
        for control in CONTROL_FIXTURE["independent_nontrivial_controls"]:
            lower = tuple(
                _hex_row(control["lower_jet"][name])
                for name in ("u", "p", "q", "p_r", "q_r")
            )
            radius = float.fromhex(control["radius_binary64_hex"])
            seeds = np.zeros((7, 1, 6), dtype=np.float64)
            for field in range(6):
                seeds[field + 1, 0, field] = 1.0
            vector = gr0_vectorized_reference_ref1_residual_batch(
                *(value[None, :] for value in lower[:3]),
                seeds,
                *(value[None, :] for value in lower[3:]),
                np.asarray((radius,), dtype=np.float64),
            ).full_residual[:, 0]
            scalar = gr0_reference_balanced_ref1_residual_batch(
                *(value[None, :] for value in lower[:3]),
                seeds,
                *(value[None, :] for value in lower[3:]),
                np.asarray((radius,), dtype=np.float64),
            ).full_residual[:, 0]
            exact = np.asarray(
                [
                    [
                        float(value)
                        for value in exact_ref1_residual(
                            *lower, radius, acceleration
                        )
                    ]
                    for acceleration in seeds[:, 0]
                ],
                dtype=np.float64,
            )
            scale = max(1.0, float(np.max(np.abs(exact))))
            self.assertLessEqual(
                float(np.max(np.abs(vector - exact))) / scale,
                64.0 * np.finfo(np.float64).eps,
            )
            self.assertLessEqual(
                float(np.max(np.abs(vector - scalar))) / scale,
                64.0 * np.finfo(np.float64).eps,
            )
            solved = diagnose_gr0_vectorized_reference_accelerations(
                *(value[None, :] for value in lower),
                np.asarray((radius,), dtype=np.float64),
            )
            self.assertTrue(solved.raw_gate_passed)

    @unittest.skipUnless(RAW_FIXTURE.is_file(), "optional CAP1 raw fixture absent")
    def test_captured_full_grid_preserves_the_strict_src3_gate(self) -> None:
        with np.load(RAW_FIXTURE, allow_pickle=False) as raw:
            arguments = tuple(
                raw[name] for name in ("u", "p", "q", "p_r", "q_r", "radii")
            )
            src3 = diagnose_gr0_reference_balanced_accelerations(*arguments)
            src4 = diagnose_gr0_vectorized_reference_accelerations(*arguments)
        self.assertTrue(src3.raw_gate_passed)
        self.assertTrue(src4.raw_gate_passed)
        self.assertLess(src4.residual_infinity, 1.0e-12)
        self.assertLessEqual(
            float(np.max(np.abs(src4.accelerations - src3.accelerations))),
            2.0e-15,
        )
        self.assertLessEqual(
            abs(
                src4.kinetic_condition_infinity_maximum
                - src3.kinetic_condition_infinity_maximum
            ),
            2.0e-9,
        )

    def test_one_bit_mutation_and_invalid_domains_remain_visible(self) -> None:
        lower = list(_compact_lower_jet())
        radius = float.fromhex(POINT_FIXTURE["coordinate_radius"])
        baseline = diagnose_gr0_vectorized_reference_accelerations(
            *(value[None, :] for value in lower),
            np.asarray((radius,), dtype=np.float64),
        )
        mutated_u = lower[0].copy()
        mutated_u[0] = np.nextafter(mutated_u[0], np.inf)
        changed = diagnose_gr0_vectorized_reference_accelerations(
            mutated_u[None, :],
            *(value[None, :] for value in lower[1:]),
            np.asarray((radius,), dtype=np.float64),
        )
        self.assertFalse(np.array_equal(changed.accelerations, baseline.accelerations))
        arrays = tuple(value[None, :] for value in lower)
        with self.assertRaisesRegex(ValueError, "radii must be positive"):
            diagnose_gr0_vectorized_reference_accelerations(
                *arrays,
                np.asarray((0.0,), dtype=np.float64),
            )
        with self.assertRaisesRegex(ValueError, "condition limit"):
            diagnose_gr0_vectorized_reference_accelerations(
                *arrays,
                np.asarray((radius,), dtype=np.float64),
                condition_number_maximum=1.0,
            )
        with self.assertRaisesRegex(ValueError, "1<tilde<hat"):
            gr0_vectorized_reference_ref1_residual_batch(
                *arrays[:3],
                np.zeros((1, 1, 6), dtype=np.float64),
                *arrays[3:],
                np.asarray((radius,), dtype=np.float64),
                tilde_normal_factor=4.0,
                hat_normal_factor=4.0,
            )


if __name__ == "__main__":
    unittest.main()
