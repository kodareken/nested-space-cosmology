import unittest

import numpy as np

from recursive_horizons.nsc_causal_common import (
    BoundaryDomain,
    CausalCommonFunctional,
    GaugeHistory,
    LinkHistory,
    MetricHistory,
    SpectralInducedSource,
)
from recursive_horizons.nsc_influence import ground_covariance


def control_problem():
    times = np.array([0.0, 0.13, 0.31])
    parent = np.repeat(np.array([[[-0.8]]], complex), len(times), axis=0)
    child = np.repeat(np.array([[[0.9]]], complex), len(times), axis=0)
    links = np.repeat(np.array([[[0.23 + 0.04j]]], complex), len(times), axis=0)
    full = np.array([[-0.8, 0.23 + 0.04j], [0.23 - 0.04j, 0.9]])
    _, _, covariance = ground_covariance(full)
    vertices = {
        "N": np.array([[0.4, 0.1], [0.1, -0.2]]),
        "beta": np.array([[0.0, 0.2j], [-0.2j, 0.0]]),
        "q": np.diag([0.3, -0.1]),
        "r": np.array([[0.0, 0.07], [0.07, 0.0]]),
    }
    metric = MetricHistory(times, parent, child, vertices)
    gauge = GaugeHistory(np.zeros_like(parent), np.zeros_like(child))
    link = LinkHistory(links)
    induced = SpectralInducedSource(
        forces={"N": 0.01, "beta": -0.02, "q": 0.03, "r": -0.04},
        coefficient_owner="finite common spectral control",
        unresolved_terms=("full continuum coefficient match",),
    )
    domain = BoundaryDomain(
        name="two-mode control",
        parent_dimension=1,
        child_dimension=1,
        spin_frame="common canonical frame",
        units="hbar=1",
        stress_factors={
            "rho": ("N", 2.0), "T_01": ("beta", 3.0),
            "p_parallel": ("q", -4.0), "p_perp": ("r", -5.0),
        },
        power_weights={"N": -0.5, "beta": 0.25, "q": 0.75},
    )
    return metric, gauge, link, covariance, induced, domain


class CausalCommonTests(unittest.TestCase):
    def test_equal_history_is_normalized_and_energy_is_balanced(self):
        result = CausalCommonFunctional().evaluate(*control_problem())
        self.assertLess(result.ward_residuals["equal_history_action"], 2e-15)
        self.assertLess(result.ward_residuals["boundary_schur_equal_history_action"], 2e-15)
        self.assertLess(abs(result.ward_residuals["energy_work_balance"]), 2e-14)
        self.assertLess(abs(result.ward_residuals["covariance_trace"]), 2e-14)
        self.assertLess(result.ward_residuals["unitarity"], 2e-14)

    def test_induced_source_is_counted_once_and_response_is_causal(self):
        result = CausalCommonFunctional().evaluate(*control_problem())
        for name in ("N", "beta", "q", "r"):
            self.assertAlmostEqual(
                result.forces[name],
                result.matter_forces[name] + result.induced_forces[name],
            )
            self.assertEqual(result.retarded_response[name][0], 0.0)
        self.assertAlmostEqual(result.stress["rho"], 2 * result.forces["N"])
        self.assertEqual(result.unresolved_terms, ("full continuum coefficient match",))

    def test_dimension_mismatch_is_rejected(self):
        metric, gauge, link, covariance, induced, _ = control_problem()
        wrong = BoundaryDomain("wrong", 2, 1, "frame", "hbar=1")
        with self.assertRaises(ValueError):
            CausalCommonFunctional().evaluate(metric, gauge, link, covariance, induced, wrong)


if __name__ == "__main__":
    unittest.main()
