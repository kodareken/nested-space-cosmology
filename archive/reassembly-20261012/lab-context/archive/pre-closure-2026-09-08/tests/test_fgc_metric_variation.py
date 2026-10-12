from __future__ import annotations

from fractions import Fraction
import unittest

from scripts.reproduce_fgc_action import DEFAULT_CONFIG as ACTION_CONFIG
from scripts.reproduce_fgc_action import load_config as load_action_config
from recursive_horizons.fgc import (
    EXPANDED_GB_COEFFICIENTS,
    branch_metric_source_identity,
    compact_gauss_bonnet_metric_source,
    contract_expanded_gauss_bonnet_basis,
    expanded_gauss_bonnet_basis,
    flrw_compact_gb_energy_density,
    flrw_lapse_variation_gb_t00,
    flrw_reduced_gb_lagrangian,
    flrw_reduced_gb_lapse_derivative,
    generate_algebraic_curvature_fixtures,
    metric_variation_certificate,
    required_var1_nonclaims,
    validate_algebraic_curvature_fixture,
)


class FGCMetricVariationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.models = load_action_config(ACTION_CONFIG).models
        self.fixtures = generate_algebraic_curvature_fixtures(20260821, 8, 3)

    def test_generator_contract_and_fixture_hashes_are_frozen(self) -> None:
        self.assertEqual(
            [fixture.fixture_id for fixture in self.fixtures],
            [
                "lcg64-00-draw-01",
                "lcg64-01-draw-02",
                "lcg64-02-draw-03",
                "lcg64-03-draw-04",
                "lcg64-04-draw-05",
                "lcg64-05-draw-08",
                "lcg64-06-draw-09",
                "lcg64-07-draw-10",
            ],
        )
        self.assertEqual(
            [fixture.source_sha256() for fixture in self.fixtures],
            [
                "6551b8447a137c622cb2aa057ef3a7ae470c9a1911465e42d65d2188ffe94006",
                "d17d4611c389098558595e84c82485d9d2853904c243f1956e918a6f1ee8e925",
                "8261a44a5b2aebc2e4d75e4d3873ad5dbc4289611b6b57cc59da5a4349571b5f",
                "4687abf3ceaa7c6508e740b2b3213237cc93aa4706b461d96d6be94c43dcff0a",
                "a27edf21354514cf6931ef8be6a7c2a6a66151329be0fb08c0b3c955064e3780",
                "c0bfd8f9e4b28e2b2b1d301806ec25b76aebe15a75c89cd008c3d43fdda6ff98",
                "0d80982995f7dd66d51d7e17a8f374ead579eef73239d8680ef2af8c05e782aa",
                "3bf066766a6c5cd0eb0b1ac486d913d46a51e448df3fb2eef74b484b5690ce76",
            ],
        )
        duplicate = generate_algebraic_curvature_fixtures(20260821, 8, 3)
        self.assertEqual(self.fixtures, duplicate)

    def test_every_exact_curvature_and_noether_contraction_gate_passes(self) -> None:
        for fixture in self.fixtures:
            with self.subTest(fixture=fixture.fixture_id):
                validation = validate_algebraic_curvature_fixture(fixture)
                for key, value in validation.items():
                    if key.endswith("_residual"):
                        self.assertEqual(value, "0", key)
                    else:
                        self.assertIs(value, True, key)

    def test_compact_and_expanded_routes_match_and_mutation_is_detected(self) -> None:
        fixture = self.fixtures[0]
        compact = compact_gauss_bonnet_metric_source(
            fixture.metric_diagonal,
            fixture.riemann_lower,
            fixture.hessian_scalar_lower,
        )
        basis = expanded_gauss_bonnet_basis(
            fixture.metric_diagonal,
            fixture.riemann_lower,
            fixture.hessian_scalar_lower,
        )
        expanded = contract_expanded_gauss_bonnet_basis(basis)
        self.assertEqual(compact, expanded)

        mutated = dict(EXPANDED_GB_COEFFICIENTS)
        mutated["riemann_hessian"] = Fraction(-7)
        self.assertNotEqual(
            compact,
            contract_expanded_gauss_bonnet_basis(basis, mutated),
        )

    def test_all_action_branches_have_required_exact_source_coverage(self) -> None:
        coverage = {
            "GR-0": {"gb": False, "nonminimal": False},
            "SGB-L": {"gb": False, "nonminimal": False},
            "FGC-QR": {"gb": False, "nonminimal": False},
        }
        for fixture in self.fixtures:
            for model in self.models:
                record = branch_metric_source_identity(model, fixture)
                self.assertEqual(
                    record["gb_compact_expanded_max_exact_residual"], "0"
                )
                self.assertEqual(
                    record["nonminimal_compact_expanded_max_exact_residual"],
                    "0",
                )
                coverage[model.model_id.value]["gb"] |= record[
                    "gb_source_nonzero"
                ]
                coverage[model.model_id.value]["nonminimal"] |= record[
                    "nonminimal_source_nonzero"
                ]
        self.assertEqual(
            coverage,
            {
                "GR-0": {"gb": False, "nonminimal": False},
                "SGB-L": {"gb": True, "nonminimal": False},
                "FGC-QR": {"gb": True, "nonminimal": True},
            },
        )

    def test_independent_flrw_lapse_variation_fixes_sign_and_factor(self) -> None:
        inputs = {
            "lapse": Fraction(1),
            "scale_factor": Fraction(7, 5),
            "hubble": Fraction(3, 2),
            "hubble_derivative": Fraction(-2, 3),
            "f_time_derivative": Fraction(5, 7),
            "f_second_time_derivative": Fraction(-11, 13),
        }
        compact = flrw_compact_gb_energy_density(
            inputs["hubble"],
            inputs["hubble_derivative"],
            inputs["f_time_derivative"],
            inputs["f_second_time_derivative"],
        )
        scale_factor_velocity = (
            inputs["lapse"] * inputs["scale_factor"] * inputs["hubble"]
        )
        reduced_lagrangian = flrw_reduced_gb_lagrangian(
            inputs["lapse"],
            scale_factor_velocity,
            inputs["f_time_derivative"],
        )
        lapse_derivative = flrw_reduced_gb_lapse_derivative(
            inputs["lapse"],
            scale_factor_velocity,
            inputs["f_time_derivative"],
        )
        lapse = flrw_lapse_variation_gb_t00(
            inputs["lapse"],
            inputs["scale_factor"],
            scale_factor_velocity,
            inputs["f_time_derivative"],
        )
        self.assertEqual(reduced_lagrangian, Fraction(-1323, 25))
        self.assertEqual(lapse_derivative, Fraction(3969, 25))
        self.assertEqual(flrw_reduced_gb_lagrangian(2, 3, 5), Fraction(-135))
        self.assertEqual(
            flrw_reduced_gb_lapse_derivative(2, 3, 5), Fraction(405, 2)
        )
        self.assertEqual(compact, Fraction(-405, 7))
        self.assertEqual(compact, lapse)

    def test_certificate_passes_only_var1_and_retains_physical_nonclaims(self) -> None:
        certificate = metric_variation_certificate(self.models, self.fixtures)
        self.assertEqual(certificate["aggregate"]["maximum_exact_residual"], "0")
        self.assertTrue(certificate["aggregate"]["all_required_coverage"])
        self.assertTrue(all(certificate["verified_exact_checks"].values()))
        self.assertEqual(certificate["nonclaims"], required_var1_nonclaims())
        self.assertTrue(
            all(value is False for value in certificate["nonclaims"].values())
        )

    def test_invalid_inputs_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            generate_algebraic_curvature_fixtures(1, 0, 3)
        with self.assertRaises(ValueError):
            generate_algebraic_curvature_fixtures(1, 1, 0)
        fixture = self.fixtures[0]
        bad_hessian = tuple(
            tuple(
                fixture.hessian_scalar_lower[a][b]
                + (Fraction(1) if (a, b) == (0, 1) else Fraction(0))
                for b in range(4)
            )
            for a in range(4)
        )
        with self.assertRaisesRegex(ValueError, "symmetric"):
            compact_gauss_bonnet_metric_source(
                fixture.metric_diagonal, fixture.riemann_lower, bad_hessian
            )
        basis = expanded_gauss_bonnet_basis(
            fixture.metric_diagonal,
            fixture.riemann_lower,
            fixture.hessian_scalar_lower,
        )
        with self.assertRaisesRegex(ValueError, "keys differ"):
            contract_expanded_gauss_bonnet_basis(
                basis, {"riemann_hessian": Fraction(-8)}
            )
        with self.assertRaisesRegex(ValueError, "lapse must be positive"):
            flrw_reduced_gb_lapse_derivative(0, 1, 1)
        with self.assertRaisesRegex(ValueError, "scale_factor must be positive"):
            flrw_lapse_variation_gb_t00(1, 0, 1, 1)


if __name__ == "__main__":
    unittest.main()
