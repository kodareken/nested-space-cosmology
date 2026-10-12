from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hyp1_bnd1_md1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical_json_text,
    _rehydrate_uhyp1,
    load_canonical_result,
    load_config,
    record,
)
from scripts.reproduce_fgc_hyp1_dom3_uhyp1 import (  # noqa: E402
    load_canonical_result as load_uhyp1_result,
)
from recursive_horizons.fgc.modified_harmonic_radial_boundary import (  # noqa: E402
    frozen_radial_modal_boundary_certificate,
)


def _contains_float(value) -> bool:
    if isinstance(value, float):
        return True
    if isinstance(value, dict):
        return any(_contains_float(item) for item in value.values())
    if isinstance(value, list):
        return any(_contains_float(item) for item in value)
    return False


class BND1MD1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_real_predecessor_boundary_certificate_and_scope(self) -> None:
        boundary = self.payload["uniform_frozen_boundary_certificate"]
        characteristics = boundary["characteristics"]
        self.assertEqual(len(characteristics["inner_incoming_indices"]), 6)
        self.assertEqual(len(characteristics["outer_incoming_indices"]), 6)
        self.assertEqual(boundary["boundary_operators"]["inner"]["rank"], 6)
        self.assertEqual(boundary["boundary_operators"]["outer"]["rank"], 6)
        self.assertTrue(boundary["energy_flux"]["homogeneous_frozen_energy_nonincreasing"])
        self.assertEqual(
            boundary["frozen_boundary_stability"]["normalized_modal_Lopatinski_determinant"],
            "1",
        )
        compatibility = self.payload["kinematic_reduction_compatibility"]
        self.assertEqual(
            compatibility["independent_reduction_constraint_boundary_conditions_required"],
            0,
        )
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))
        self.assertFalse(_contains_float(self.payload))

    def test_canonical_record_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))
        self.assertEqual(self.payload, load_canonical_result(DEFAULT_OUTPUT))

    def test_config_result_and_speed_crossing_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bad_config = root / "bad.toml"
            bad_config.write_text(
                source.replace(
                    "constraint_preserving_ACT1_IBVP_proven = false",
                    "constraint_preserving_ACT1_IBVP_proven = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(bad_config)
            bad_json = root / "bad.json"
            bad_json.write_text(
                _canonical_json_text(self.payload).rstrip(), encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(bad_json)

        serialized = deepcopy(load_uhyp1_result())
        speed = serialized["uniform_radial_certificate"]["modes"]["hat"][0]["speed_box"]
        speed["lower"], speed["upper"] = "-1", "1"
        with self.assertRaisesRegex(ValueError, "crosses zero"):
            frozen_radial_modal_boundary_certificate(_rehydrate_uhyp1(serialized))


if __name__ == "__main__":
    unittest.main()
