from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.exact_interval import interval  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_radial_boundary import (  # noqa: E402
    frozen_radial_modal_boundary_certificate,
)


def synthetic_uhyp1():
    negative = (interval(-2, -1), interval(-3, -2))
    positive = (interval(1, 2), interval(2, 3))

    def modes(count: int):
        values = []
        for index in range(count):
            source = negative[index % 2] if index < count // 2 else positive[index % 2]
            values.append({"speed_box": source})
        return values

    return {
        "classification": "exact_rational_nonzero_compact_radial_REF1_branch_strong_hyperbolicity_certificate",
        "domain": {
            "all_30_parameter_axes_nonzero": True,
            "all_6_acceleration_axes_nonzero": True,
        },
        "modes": {
            "tilde": modes(4),
            "hat": modes(4),
            "physical_chi": modes(2),
            "regulator": modes(2),
        },
        "eigenframe": {
            "all_enclosed_frames_invertible": True,
            "real_smooth_radial_eigenframe": True,
            "neumann_inverse": {"dimension": 12, "rho_infinity": Q(1, 10)},
        },
        "radial_symmetrizer": {
            "definition": "H=V^-T V^-1",
            "HA_symmetric": True,
            "euclidean_coercivity_lower_bound": Q(1, 20),
            "euclidean_coercivity_upper_bound": Q(30),
        },
        "nonclaims": {"IBVP_proven": False, "evolution_authorized": False},
    }


class RadialBoundaryTests(unittest.TestCase):
    def test_exact_modal_boundary_split_and_flux(self) -> None:
        out = frozen_radial_modal_boundary_certificate(synthetic_uhyp1())
        characteristic = out["characteristics"]
        self.assertEqual(len(characteristic["outer_incoming_indices"]), 6)
        self.assertEqual(len(characteristic["inner_incoming_indices"]), 6)
        self.assertEqual(
            characteristic["outer_incoming_indices"],
            characteristic["inner_outgoing_indices"],
        )
        self.assertEqual(
            characteristic["inner_incoming_indices"],
            characteristic["outer_outgoing_indices"],
        )
        self.assertEqual(out["boundary_operators"]["outer"]["rank"], 6)
        self.assertEqual(out["boundary_operators"]["inner"]["rank"], 6)
        self.assertEqual(
            out["frozen_boundary_stability"]["normalized_modal_Lopatinski_determinant"],
            Q(1),
        )
        self.assertTrue(out["energy_flux"]["homogeneous_frozen_energy_nonincreasing"])
        self.assertTrue(out["uniform_frozen_radial_main_system_boundary_dissipation_proven"])
        self.assertTrue(all(value is False for value in out["nonclaims"].values()))

    def test_speed_frame_symmetrizer_and_nonclaim_mutations_fail_closed(self) -> None:
        crossing = synthetic_uhyp1()
        crossing["modes"]["hat"][0]["speed_box"] = interval(-1, 1)
        with self.assertRaisesRegex(ValueError, "crosses zero"):
            frozen_radial_modal_boundary_certificate(crossing)

        frame = synthetic_uhyp1()
        frame["eigenframe"]["all_enclosed_frames_invertible"] = False
        with self.assertRaisesRegex(ValueError, "eigenframe"):
            frozen_radial_modal_boundary_certificate(frame)

        symmetrizer = synthetic_uhyp1()
        symmetrizer["radial_symmetrizer"]["euclidean_coercivity_lower_bound"] = Q(0)
        with self.assertRaisesRegex(ValueError, "coercivity"):
            frozen_radial_modal_boundary_certificate(symmetrizer)

        promoted = synthetic_uhyp1()
        promoted["nonclaims"]["IBVP_proven"] = True
        with self.assertRaisesRegex(ValueError, "promoted"):
            frozen_radial_modal_boundary_certificate(promoted)

    def test_wrong_count_and_boundary_rank_mutations_fail_closed(self) -> None:
        missing = synthetic_uhyp1()
        missing["modes"]["regulator"].pop()
        with self.assertRaisesRegex(ValueError, "twelve"):
            frozen_radial_modal_boundary_certificate(missing)

        def duplicate_selector(indices, dimension):
            rows = tuple(
                tuple(Q(int(column == indices[0])) for column in range(dimension))
                for _ in indices
            )
            return rows

        with patch(
            "recursive_horizons.fgc.modified_harmonic_radial_boundary._selector",
            side_effect=duplicate_selector,
        ):
            with self.assertRaisesRegex(ValueError, "rank"):
                frozen_radial_modal_boundary_certificate(synthetic_uhyp1())


if __name__ == "__main__":
    unittest.main()
