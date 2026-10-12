from __future__ import annotations

from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.def1_geometry_error import (  # noqa: E402
    INPUT_NAMES,
    JET_COMPONENTS,
    METRIC_FIELDS,
    QGeometryPremiseError,
    QGeometryResourceError,
    evaluate_q_box,
    q_input_error_bound,
    q_sensitivity_enclosures,
)
from recursive_horizons.fgc.exact_interval import Interval  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    Jet2,
    SphericalState,
    _hessian,
    direct_4d_curvature,
)


def _minkowski():
    inputs = {name: Q(0) for name in INPUT_NAMES}
    inputs.update({"alpha.value": Q(1), "lambda.value": Q(1),
                   "R.value": Q(2), "R.dr": Q(1), "k.t": Q(1), "k.r": Q(1)})
    return inputs


def _de_sitter():
    # ds^2=-dt^2+exp(-2t)(dr^2+r^2 dOmega^2), t=0, r=2.
    inputs = _minkowski()
    inputs.update({"lambda.dt": -Q(1), "lambda.dtt": Q(1),
                   "R.dt": -Q(2), "R.dtt": Q(2), "R.dtr": -Q(1)})
    return inputs


def _direct(inputs):
    # Independent direct-4D curvature route; no producer evaluation helper.
    jets = {field: Jet2(**{part: inputs[f"{field}.{part}"] for part in JET_COMPONENTS})
            for field in METRIC_FIELDS}
    a, v, lam, radius = (jets[field] for field in METRIC_FIELDS)
    state = SphericalState(
        h_tt=-(a * a) + lam * lam * v * v,
        h_tr=lam * lam * v, h_rr=lam * lam,
        areal_radius=radius, phi=Jet2(0), chi=Jet2(0),
        planck_mass=Q(1), mu=Q(1), g4=Q(1),
        alpha=Q(0), beta=Q(0), eta=Q(0), branch="GR-0",
    )
    curvature, inverse, connection = direct_4d_curvature(state)
    ricci = [[sum(inverse[c][d] * curvature[c][a][d][b]
                  for c, d in product(range(4), repeat=2)) for b in range(4)]
             for a in range(4)]
    tangent = (inputs["k.t"], inputs["k.r"])
    contraction = sum(ricci[a][b] * tangent[a] * tangent[b]
                      for a, b in product(range(2), repeat=2))
    theta = 2 * (radius.dt * tangent[0] + radius.dr * tangent[1]) / radius.value
    hessian = _hessian(radius, connection)
    shortcut = -2 * sum(hessian[a][b] * tangent[a] * tangent[b]
                        for a, b in product(range(2), repeat=2)) / radius.value
    return {"q": -theta**2 / 2 - contraction, "theta": theta,
            "ricci_kk": contraction, "null_hessian_shortcut": shortcut}


class QGeometryControlTests(unittest.TestCase):
    def test_minkowski_and_trapped_contracting_de_sitter_are_not_defocusing(self):
        flat = evaluate_q_box(_minkowski())
        self.assertEqual(flat["q"], -Q(1, 2))
        self.assertEqual(flat["theta"], 1)
        self.assertEqual(flat["ricci_kk"], 0)
        self.assertEqual(flat["metric_null_residual"], 0)
        collapsing = _de_sitter()
        outgoing = evaluate_q_box(collapsing)
        ingoing = evaluate_q_box(collapsing | {"k.r": -Q(1)})
        self.assertEqual(outgoing["theta"], -1)
        self.assertEqual(ingoing["theta"], -3)
        self.assertEqual(outgoing["q"], -Q(1, 2))
        self.assertEqual(ingoing["q"], -Q(9, 2))

    def test_full_ricci_retains_the_off_null_base_curvature_term(self):
        nonnull = _de_sitter() | {"k.r": Q(2)}
        observed = evaluate_q_box(nonnull)
        direct = _direct(nonnull)
        self.assertEqual(observed["metric_null_residual"], 3)
        self.assertEqual(observed["ricci_kk"], 9)
        self.assertEqual(observed["q"], -9)
        self.assertEqual(direct["null_hessian_shortcut"], 6)
        self.assertEqual(direct["ricci_kk"] - direct["null_hessian_shortcut"], 3)

    def test_positive_coordinate_gtt_does_not_remove_a_timelike_dt_foliation(self):
        # Flat space with areal radius R=r+2t: g_tt=3 while g^{tt}=-1.
        shifted = _minkowski() | {"shift.value": Q(2), "R.dt": Q(2), "k.r": -Q(1)}
        observed = evaluate_q_box(shifted)
        self.assertEqual(observed["metric_null_residual"], 0)
        self.assertEqual(observed["theta"], 1)
        self.assertEqual(observed["q"], -Q(1, 2))
        self.assertEqual(observed["q"], _direct(shifted)["q"])

    def test_generic_nonzero_shift_two_jets_match_direct_four_dimensions(self):
        for index in (1, 2, 3):
            with self.subTest(index=index):
                inputs = {name: Q((place + index) % 7 - 3, 11 + place)
                          for place, name in enumerate(INPUT_NAMES)}
                inputs.update({"alpha.value": Q(2), "shift.value": Q(index, 7),
                               "lambda.value": Q(3, 2), "R.value": Q(4),
                               "k.t": Q(3, 2), "k.r": -Q(index, 5)})
                expected, observed = _direct(inputs), evaluate_q_box(inputs)
                for name in ("theta", "ricci_kk", "q"):
                    self.assertEqual(observed[name], expected[name])

    def test_positive_q_point_does_not_supply_trappedness_or_a_physical_gate(self):
        marginal = _minkowski() | {"R.dr": Q(0), "R.dtt": Q(1)}
        observed = evaluate_q_box(marginal)
        self.assertEqual(observed["theta"], 0)
        self.assertEqual(observed["q"], 1)
        self.assertEqual(set(observed), {"theta", "q", "ricci_kk", "metric_null_residual"})


class QGeometrySensitivityTests(unittest.TestCase):
    def test_exact_singleton_derivatives_have_the_complete_input_inventory(self):
        derivatives = q_sensitivity_enclosures(_minkowski())
        self.assertEqual(tuple(derivatives), INPUT_NAMES)
        self.assertEqual(derivatives["k.r"], -1)
        self.assertEqual(derivatives["k.t"], 0)
        self.assertEqual(derivatives["R.dtt"], 1)

    def test_linear_injected_curvature_error_has_an_exact_sharp_bound(self):
        nominal = _minkowski()
        radii = {name: Q(0) for name in INPUT_NAMES}
        epsilon = radii["R.dtt"] = Q(1, 64)
        bound = q_input_error_bound(nominal, radii)
        self.assertEqual(bound, epsilon)
        original = _direct(nominal)["q"]
        for error in (-epsilon, epsilon / 2, epsilon):
            actual = _direct(nominal | {"R.dtt": error})["q"]
            self.assertLessEqual(abs(actual - original), bound)
        # The same input uncertainty is not reduced by a favorable Q value.
        positive = nominal | {"R.dtt": Q(1)}
        self.assertGreater(_direct(positive)["q"], 0)
        self.assertEqual(q_input_error_bound(positive, radii), bound)

    def test_whole_box_derivatives_cover_nonlinear_and_correlated_errors(self):
        nominal = _minkowski()
        radii = {name: Q(0) for name in INPUT_NAMES}
        radii["R.dtt"], radii["k.r"] = Q(1, 50), Q(1, 4)
        box = {name: Interval(nominal[name] - radii[name], nominal[name] + radii[name])
               for name in INPUT_NAMES}
        derivative = q_sensitivity_enclosures(box)["k.r"]
        self.assertEqual(derivative, Interval(-Q(5, 4), -Q(3, 4)))
        bound = q_input_error_bound(nominal, radii)
        self.assertEqual(bound, Q(1, 50) + Q(5, 16))
        original = _direct(nominal)["q"]
        enclosure = evaluate_q_box(box)["q"]
        for sign_r, sign_k in product((-1, 1), repeat=2):
            varied = nominal | {"R.dtt": sign_r * radii["R.dtt"],
                                "k.r": 1 + sign_k * radii["k.r"]}
            actual = _direct(varied)["q"]
            self.assertLessEqual(abs(actual - original), bound)
            self.assertLessEqual(enclosure.lower, actual)
            self.assertLessEqual(actual, enclosure.upper)


class QGeometryInputTests(unittest.TestCase):
    def test_missing_and_extra_inputs_are_not_filled_in(self):
        valid = _minkowski()
        for name in INPUT_NAMES:
            missing = {key: value for key, value in valid.items() if key != name}
            with self.subTest(name=name), self.assertRaises(QGeometryPremiseError):
                evaluate_q_box(missing)
        with self.assertRaises(QGeometryPremiseError):
            evaluate_q_box(valid | {"measured_q": Q(1)})
        with self.assertRaises(QGeometryPremiseError):
            q_input_error_bound(valid, {})

    def test_float_bool_and_nonfinite_aliases_refuse_in_every_slot(self):
        valid = _minkowski()
        for name, bad in product(INPUT_NAMES, (True, 1.0, float("inf"), float("nan"), "1")):
            with self.subTest(name=name, bad=bad), self.assertRaises(TypeError):
                evaluate_q_box(valid | {name: bad})

    def test_whole_box_domain_and_input_rational_limits_are_fail_closed(self):
        valid = _minkowski()
        for name in ("alpha.value", "lambda.value", "R.value", "k.t"):
            with self.subTest(name=name), self.assertRaises(QGeometryPremiseError):
                evaluate_q_box(valid | {name: Interval(-1, 1)})
        with self.assertRaises(QGeometryResourceError):
            evaluate_q_box(valid | {"R.dtt": 1 << 5000})
        radii = {name: Q(0) for name in INPUT_NAMES}
        with self.assertRaises(QGeometryPremiseError):
            q_input_error_bound(valid, radii | {"R.dtt": -Q(1, 2)})
        with self.assertRaises(TypeError):
            q_input_error_bound(valid, radii | {"R.dtt": True})
        with self.assertRaises(QGeometryPremiseError):
            q_input_error_bound(valid, radii | {"R.value": Q(2)})


if __name__ == "__main__":
    unittest.main()
