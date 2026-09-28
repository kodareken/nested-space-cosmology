from __future__ import annotations

from dataclasses import replace
from fractions import Fraction as Q
import json
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY))
sys.path.insert(0, str(REPOSITORY / "src"))

from scripts.reproduce_fgc_hyp1_con1_comp1 import load_config as load_comp1_config  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_cauchy import (  # noqa: E402
    conditional_gauge_cauchy_certificate,
)
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)


class ModifiedHarmonicCauchyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        config = load_comp1_config()
        cls.state = activated_compatible_state(config.qift1_flat_fixture)["state"]
        cls.comp1 = json.loads(
            (REPOSITORY / "results/fgc-1-hyp1-con1-comp1.json").read_text()
        )
        cls.metric = json.loads(
            (REPOSITORY / "results/fgc-1-hyp1-con2-mprop1.json").read_text()
        )["controls"]["activated_comp1_deterministic_third_jet"][
            "metric_derived_gauge_subsidiary"
        ]
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=Q(1, 2)
        )

    def certificate(self, **changes):
        arguments = {
            "state": self.state,
            "reference": self.reference,
            "coordinate_radius": Q(4),
            "tilde_normal_factor": Q(4),
            "hat_normal_factor": Q(9),
            "metric_subsidiary_record": self.metric,
            "comp1_constraint_record": self.comp1,
        }
        arguments.update(changes)
        return conditional_gauge_cauchy_certificate(**arguments)

    def test_exact_point_witness_and_conditional_theorem(self) -> None:
        output = self.certificate()
        self.assertTrue(output["subsidiary_operator"]["normally_hyperbolic_at_COMP1"])
        self.assertTrue(output["compatible_point_witness"]["physical_H_zero"])
        self.assertTrue(output["compatible_point_witness"]["physical_M_zero"])
        self.assertNotEqual(
            output["compatible_point_witness"]["normal_gauge_extension_map"]["determinant"],
            0,
        )
        self.assertTrue(
            output["conditional_theorem"][
                "conditional_boundary_free_Cauchy_uniqueness_statement_derived"
            ]
        )
        self.assertTrue(all(value is False for value in output["nonclaims"].values()))

    def test_missing_scalar_shell_and_nonclaim_promotion_fail_closed(self) -> None:
        metric = dict(self.metric)
        metric["both_scalar_equations_required"] = False
        with self.assertRaisesRegex(ValueError, "premise"):
            self.certificate(metric_subsidiary_record=metric)

        metric = dict(self.metric)
        metric["nonclaims"] = dict(metric["nonclaims"])
        metric["nonclaims"]["initial_boundary_value_problem_proven"] = True
        with self.assertRaisesRegex(ValueError, "promoted"):
            self.certificate(metric_subsidiary_record=metric)

    def test_nonpositive_F_and_wrong_auxiliary_order_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "effective Planck coefficient F must be positive"):
            replace(self.state, beta=Q(-(2**36)))
        with self.assertRaisesRegex(ValueError, "1 < tilde"):
            self.certificate(hat_normal_factor=Q(4))


if __name__ == "__main__":
    unittest.main()
