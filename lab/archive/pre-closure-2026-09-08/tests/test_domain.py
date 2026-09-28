from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.core import G  # noqa: E402
from recursive_horizons.domain import (  # noqa: E402
    DomainScalarSpec,
    dust_scaled_density,
    exchange_ledger,
    instantaneous_split,
    integer_damped_period_benchmark,
    kottler_balance_radius,
    kottler_radial_acceleration,
    monomial_rapid_oscillation,
    quadratic_cycle_average,
    underdamped_spectator_state,
)


class DomainScalarFluidTests(unittest.TestCase):
    def setUp(self) -> None:
        self.spec = DomainScalarSpec(vacuum_energy=9.0, mass=2.0, decay_rate=0.25)

    def test_immutable_validated_spec(self) -> None:
        self.assertEqual(self.spec.vacuum_energy, 9.0)
        with self.assertRaises((AttributeError, TypeError)):
            self.spec.mass = 3.0  # type: ignore[misc]
        for values in ((-1.0, 1.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, -0.1)):
            with self.subTest(values=values):
                with self.assertRaises(ValueError):
                    DomainScalarSpec(*values)

    def test_instantaneous_split_and_cycle_average(self) -> None:
        state = instantaneous_split(self.spec, displacement=3.0, velocity=4.0)
        self.assertAlmostEqual(state["rho_vacuum"], 9.0)
        self.assertAlmostEqual(state["p_vacuum"], -9.0)
        self.assertEqual(state["w_vacuum"], -1.0)
        self.assertAlmostEqual(state["rho_oscillatory"], 26.0)
        self.assertAlmostEqual(state["p_oscillatory"], -10.0)
        self.assertAlmostEqual(state["rho_total"], 35.0)
        self.assertAlmostEqual(state["p_total"], -19.0)
        self.assertAlmostEqual(state["w_total"], -19.0 / 35.0)

        averaged = quadratic_cycle_average(self.spec, amplitude=3.0)
        self.assertAlmostEqual(averaged["rho_oscillatory"], 18.0)
        self.assertEqual(averaged["w_oscillatory"], 0.0)
        self.assertAlmostEqual(averaged["w_total"], -9.0 / 27.0)

    def test_monomial_and_dust_scaling(self) -> None:
        quadratic = monomial_rapid_oscillation(2.0)
        self.assertEqual(quadratic["w_average"], 0.0)
        self.assertEqual(quadratic["density_scale_exponent"], 3.0)
        quartic = monomial_rapid_oscillation(4.0)
        self.assertAlmostEqual(quartic["w_average"], 1.0 / 3.0)
        self.assertEqual(quartic["density_scale_exponent"], 4.0)
        self.assertAlmostEqual(dust_scaled_density(80.0, 2.0, 4.0), 10.0)


class SpectatorAndExchangeTests(unittest.TestCase):
    def test_integer_period_endpoint_identity(self) -> None:
        spec = DomainScalarSpec(vacuum_energy=1.0, mass=10.0, decay_rate=0.4)
        benchmark = integer_damped_period_benchmark(spec, amplitude=2.0, hubble=0.2, periods=9)
        self.assertAlmostEqual(
            benchmark["compensated_rho_oscillatory"],
            benchmark["initial_rho_oscillatory"],
            places=10,
        )
        self.assertAlmostEqual(benchmark["endpoint_residual"], 0.0, places=10)

    def test_exact_spectator_initial_conditions_and_exchange_ledger(self) -> None:
        spec = DomainScalarSpec(vacuum_energy=2.0, mass=5.0, decay_rate=0.7)
        state = underdamped_spectator_state(spec, amplitude=3.0, hubble=0.3, time=0.0)
        self.assertAlmostEqual(state["displacement"], 3.0)
        self.assertAlmostEqual(state["velocity"], 0.0)
        self.assertAlmostEqual(state["scale_factor"], 1.0)
        ledger = exchange_ledger(spec, displacement=2.0, velocity=3.0, hubble=0.4, radiation_density=7.0)
        self.assertAlmostEqual(ledger["transfer_q"], 0.7 * 9.0)
        self.assertAlmostEqual(ledger["total_continuity_residual"], 0.0, places=13)

    def test_invalid_inputs_and_json_safe_outputs(self) -> None:
        spec = DomainScalarSpec(vacuum_energy=0.0, mass=2.0)
        invalid = (True, "1", complex(1.0, 0.0), math.inf, math.nan)
        for value in invalid:
            with self.subTest(value=repr(value)):
                with self.assertRaises(ValueError):
                    instantaneous_split(spec, value, 0.0)
        with self.assertRaises(ValueError):
            instantaneous_split(
                DomainScalarSpec(vacuum_energy=0.0, mass=1.0e308),
                1.0e308,
                0.0,
            )
        for value in (0.0, -1.0, True, math.inf, math.nan):
            with self.subTest(exponent=repr(value)):
                with self.assertRaises(ValueError):
                    monomial_rapid_oscillation(value)
        with self.assertRaises(ValueError):
            underdamped_spectator_state(spec, 1.0, 2.0, 0.0)
        with self.assertRaises(ValueError):
            integer_damped_period_benchmark(spec, 1.0, 0.1, True)
        outputs = [
            instantaneous_split(spec, 0.0, 0.0),
            quadratic_cycle_average(spec, 0.0),
            monomial_rapid_oscillation(2.0),
            integer_damped_period_benchmark(spec, 1.0, 0.1, 2),
            exchange_ledger(spec, 1.0, 1.0, 0.1, 1.0),
        ]
        json.dumps(outputs, allow_nan=False)
        for output in outputs:
            self.assertTrue(
                all(
                    value is None or math.isfinite(value)
                    for value in output.values()
                )
            )
        self.assertIsNone(instantaneous_split(spec, 0.0, 0.0)["w_vacuum"])
        self.assertIsNone(quadratic_cycle_average(spec, 0.0)["w_vacuum"])

    def test_endpoint_rejects_antidamping_and_out_of_range_compensation(self) -> None:
        spec = DomainScalarSpec(vacuum_energy=0.0, mass=10.0, decay_rate=0.0)
        with self.assertRaises(ValueError):
            underdamped_spectator_state(
                spec,
                amplitude=1.0,
                hubble=-1.0,
                time=1.0,
            )
        with self.assertRaises(ValueError):
            integer_damped_period_benchmark(
                spec,
                amplitude=1.0,
                hubble=1.0,
                periods=1_000,
            )


class KottlerBenchmarkTests(unittest.TestCase):
    def test_balance_radius_zeros_acceleration_and_changes_sign(self) -> None:
        mass_kg = 2.0e30
        lambda_m2 = 1.0e-52
        radius = kottler_balance_radius(mass_kg, lambda_m2)
        term_scale = G * mass_kg / radius**2
        residual = kottler_radial_acceleration(mass_kg, radius, lambda_m2)
        self.assertLess(abs(residual) / term_scale, 1.0e-12)
        self.assertLess(
            kottler_radial_acceleration(mass_kg, radius / 2.0, lambda_m2),
            0.0,
        )
        self.assertGreater(
            kottler_radial_acceleration(mass_kg, radius * 2.0, lambda_m2),
            0.0,
        )

    def test_kottler_requires_strict_positive_finite_inputs(self) -> None:
        for bad in (0.0, -1.0, True, math.inf, math.nan):
            with self.subTest(bad=repr(bad)):
                with self.assertRaises(ValueError):
                    kottler_balance_radius(bad, 1.0e-52)
                with self.assertRaises(ValueError):
                    kottler_radial_acceleration(1.0e30, bad, 1.0e-52)
                with self.assertRaises(ValueError):
                    kottler_radial_acceleration(1.0e30, 1.0e10, 1.0e-52, bad, 1.0)
                with self.assertRaises(ValueError):
                    kottler_balance_radius(1.0e30, 1.0e-52, 1.0, bad)


if __name__ == "__main__":
    unittest.main()
