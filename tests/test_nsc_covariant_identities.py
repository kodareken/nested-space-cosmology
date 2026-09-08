import json
import unittest

import sympy as sp

from recursive_horizons.nsc_covariant_identities import (
    dirac_heat_coefficients, exact_checks, flat_measure_identities,
    lichnerowicz_identities, paired_operator_identities,
)


class CovariantIdentitiesTests(unittest.TestCase):
    @staticmethod
    def expression(value):
        x = sp.symbols("x", real=True)
        names = {name: sp.Function(name)
                 for name in ("N", "q", "r", "f", "u0", "u1", "u2", "u3")}
        names["x"] = x
        return sp.sympify(value, locals=names), x, names

    def test_lapse_spin_connection_is_required_by_flat_measure(self):
        result = flat_measure_identities()
        self.assertEqual(result["hilbert_weight_residual"], "0")
        self.assertEqual(result["flat_radial_operator_residual"], "0")
        defect, x, names = self.expression(result["omitted_lapse_connection_defect"])
        value = defect.subs({names["N"](x): sp.exp(x), names["q"](x): 1,
                             names["f"](x): 1}).doit().subs(x, 0)
        self.assertEqual(value, sp.I / 2)

    def test_full_frequency_square_detects_a_lapse_gradient(self):
        result = paired_operator_identities()
        for name in ("physical_chirality_residual", "frequency_full_square_residual",
                     "ultrastatic_square_residual", "local_weyl_sandwich_residual",
                     "metric_tangent_residual", "uniform_lapse_frequency_residual"):
            self.assertEqual(result[name], [["0"] for _ in range(4)], name)
        values = []
        for row in result["omitted_lapse_square_defect"]:
            expression, x, names = self.expression(row[0])
            substitutions = {names["N"](x): sp.exp(x), names["q"](x): 1,
                             sp.Symbol("omega"): 2}
            substitutions.update({names[f"u{j}"](x): int(j == 0) for j in range(4)})
            values.append(expression.subs(substitutions).doit().subs(x, 0))
        self.assertEqual(values, [0, 2 * sp.I, 0, 0])

    def test_coordinate_curvature_and_lichnerowicz_recover_round_sphere(self):
        result = lichnerowicz_identities()
        self.assertEqual(result["full_unseparated_lichnerowicz_residual"],
                         [["0"] * 4 for _ in range(4)])
        self.assertEqual(result["scalar_curvature_warped_product_residual"], "0")
        expression, x, names = self.expression(result["scalar_curvature"])
        radius = sp.symbols("radius", positive=True)
        cylinder = expression.subs({names["N"](x): 1, names["q"](x): 1,
                                    names["r"](x): radius}).doit()
        self.assertEqual(sp.simplify(cylinder), 2 / radius**2)

    def test_clifford_trace_fixes_dirac_heat_sign_and_multiplicity(self):
        result = dirac_heat_coefficients()
        self.assertEqual(result["spin_curvature_riemann_ratio"], "-1/2")
        self.assertEqual(result["spin_curvature_contraction_residual"], "0")
        self.assertEqual(result["a0"], "4")
        self.assertEqual(result["a2"], "-R/3")
        scalar, ricci2, riemann2, box = sp.symbols("R Ricci2 Riemann2 BoxR")
        expected = (5 * scalar**2 - 8 * ricci2 - 7 * riemann2 - 12 * box) / 360
        self.assertEqual(sp.simplify(sp.sympify(result["a4"]) - expected), 0)
        self.assertEqual(result["cylinder_a2"], "-2/(3*r**2)")
        self.assertEqual(result["cylinder_a4"], "-1/(15*r**4)")

    def test_exact_packet_is_json_serializable_and_reproducible(self):
        first = json.dumps(exact_checks(), sort_keys=True)
        exact_checks.cache_clear()
        self.assertEqual(json.dumps(exact_checks(), sort_keys=True), first)


if __name__ == "__main__":
    unittest.main()
