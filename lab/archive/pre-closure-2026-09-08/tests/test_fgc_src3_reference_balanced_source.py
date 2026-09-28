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

from recursive_horizons.fgc.evolution.gr0_direct_source import (  # noqa: E402
    gr0_ref1_residual_batch,
)
from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    exact_affine_system,
    exact_gaussian_solve,
    exact_ref1_residual,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    diagnose_gr0_reference_balanced_accelerations,
    gr0_reference_balanced_ref1_residual_batch,
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


class FGCSRC3ReferenceBalancedSourceTests(unittest.TestCase):
    def test_exact_spherical_minkowski_is_bitwise_zero(self) -> None:
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
        residual = gr0_reference_balanced_ref1_residual_batch(
            u,
            p,
            q,
            np.zeros((1, points, 6), dtype=np.float64),
            p_r,
            q_r,
            radii,
        )
        for value in (
            residual.full_residual,
            residual.unredefined_metric_residual,
            residual.scalar_residual,
            residual.gauge_constraint,
            residual.hamiltonian_constraint,
            residual.momentum_constraint,
            residual.ricci_scalar,
            residual.ricci_squared,
        ):
            self.assertTrue(np.array_equal(value, np.zeros_like(value)))
        solved = diagnose_gr0_reference_balanced_accelerations(
            u, p, q, p_r, q_r, radii
        )
        self.assertTrue(solved.raw_gate_passed)
        self.assertEqual(solved.residual_infinity, 0.0)
        self.assertTrue(
            np.array_equal(solved.accelerations, np.zeros_like(solved.accelerations))
        )

    def test_pref11_point_converges_to_the_exact_oracle(self) -> None:
        lower = _compact_lower_jet()
        radius = float.fromhex(POINT_FIXTURE["coordinate_radius"])
        exact_constant, exact_matrix = exact_affine_system(*lower, radius)
        exact_root = exact_gaussian_solve(
            exact_matrix, tuple(-value for value in exact_constant)
        )
        rounded_exact_root = np.asarray(
            [float(value) for value in exact_root], dtype=np.float64
        )
        solved = diagnose_gr0_reference_balanced_accelerations(
            *(value[None, :] for value in lower),
            np.asarray((radius,), dtype=np.float64),
        )
        exact_residual = exact_ref1_residual(
            *lower, radius, solved.accelerations[0]
        )
        self.assertTrue(solved.raw_gate_passed)
        self.assertLess(solved.residual_infinity, 1.0e-24)
        self.assertLess(
            float(max((abs(value) for value in exact_residual), default=Fraction(0))),
            1.0e-24,
        )
        self.assertLessEqual(
            max(
                _ulp_distance(left, right)
                for left, right in zip(
                    solved.accelerations[0], rounded_exact_root, strict=True
                )
            ),
            8,
        )
        legacy = gr0_ref1_residual_batch(
            *(value[None, :] for value in lower[:3]),
            solved.accelerations[None, :, :],
            *(value[None, :] for value in lower[3:]),
            np.asarray((radius,), dtype=np.float64),
        ).full_residual[0, 0]
        self.assertGreater(float(np.max(np.abs(legacy))), 1.0e-12)

    def test_two_independent_nontrivial_controls_match_exact_zero_unit_seeds(self) -> None:
        epsilon_bound = 16.0 * np.finfo(np.float64).eps
        for control in CONTROL_FIXTURE["independent_nontrivial_controls"]:
            lower = tuple(
                _hex_row(control["lower_jet"][name])
                for name in ("u", "p", "q", "p_r", "q_r")
            )
            radius = float.fromhex(control["radius_binary64_hex"])
            seeds = np.zeros((7, 1, 6), dtype=np.float64)
            for field in range(6):
                seeds[field + 1, 0, field] = 1.0
            stable = gr0_reference_balanced_ref1_residual_batch(
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
            normalized_error = float(np.max(np.abs(stable - exact))) / max(
                1.0, float(np.max(np.abs(exact)))
            )
            self.assertLessEqual(normalized_error, epsilon_bound)
            solved = diagnose_gr0_reference_balanced_accelerations(
                *(value[None, :] for value in lower),
                np.asarray((radius,), dtype=np.float64),
            )
            self.assertTrue(solved.raw_gate_passed)

    @unittest.skipUnless(RAW_FIXTURE.is_file(), "optional CAP1 raw fixture absent")
    def test_captured_full_grid_clears_only_the_unchanged_raw_source_gate(self) -> None:
        with np.load(RAW_FIXTURE, allow_pickle=False) as raw:
            solved = diagnose_gr0_reference_balanced_accelerations(
                *(raw[name] for name in ("u", "p", "q", "p_r", "q_r", "radii"))
            )
        self.assertTrue(solved.raw_gate_passed)
        self.assertEqual(solved.residual_infinity, 3.296668493746324e-14)
        self.assertEqual(solved.kinetic_condition_infinity_maximum, 2859227.196975944)
        self.assertEqual(solved.iterations[-1].maximum_residual_point_index, 872)
        self.assertEqual(solved.iterations[-1].maximum_residual_row_index, 3)

    def test_one_bit_input_mutation_remains_visible(self) -> None:
        lower = list(_compact_lower_jet())
        radius = float.fromhex(POINT_FIXTURE["coordinate_radius"])
        baseline = diagnose_gr0_reference_balanced_accelerations(
            *(value[None, :] for value in lower), np.asarray((radius,))
        )
        mutated_u = lower[0].copy()
        mutated_u[0] = np.nextafter(mutated_u[0], np.inf)
        changed = diagnose_gr0_reference_balanced_accelerations(
            mutated_u[None, :],
            *(value[None, :] for value in lower[1:]),
            np.asarray((radius,)),
        )
        self.assertFalse(np.array_equal(changed.accelerations, baseline.accelerations))

    def test_invalid_domain_and_condition_inputs_fail_closed(self) -> None:
        lower = _compact_lower_jet()
        radius = float.fromhex(POINT_FIXTURE["coordinate_radius"])
        arrays = tuple(value[None, :] for value in lower)
        with self.assertRaisesRegex(ValueError, "radii must be positive"):
            diagnose_gr0_reference_balanced_accelerations(
                *arrays, np.asarray((0.0,))
            )
        with self.assertRaisesRegex(ValueError, "condition limit"):
            diagnose_gr0_reference_balanced_accelerations(
                *arrays,
                np.asarray((radius,)),
                condition_number_maximum=1.0,
            )
        with self.assertRaisesRegex(ValueError, "1<tilde<hat"):
            gr0_reference_balanced_ref1_residual_batch(
                *arrays[:3],
                np.zeros((1, 1, 6)),
                *arrays[3:],
                np.asarray((radius,)),
                tilde_normal_factor=4.0,
                hat_normal_factor=4.0,
            )


if __name__ == "__main__":
    unittest.main()
