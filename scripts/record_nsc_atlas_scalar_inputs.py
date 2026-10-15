#!/usr/bin/env python3
"""Create or authenticate finite scalar lift and cumulative-forcing inputs.

The comparison includes its preparation in one cumulative compute budget.
Authentication checks saved bytes and algebra without a physical solve.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from record_nsc_atlas_scalar_boundary import (
    SOURCES as BOUNDARY_SOURCES, BASE as BOUNDARY_BASE,
    TIMES, PROFILE_COUNTS, digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-scalar-inputs-comparison-v1"
STATUS = "FORMAL_CUBIC_SCALAR_INPUT_REFINEMENT"
SCOPE = "formal cubic scalar lift and forcing inputs; no evolved correction, central waveform or finite-amplitude remainder"
BASE = {**BOUNDARY_BASE, "field_nodes": 129}
VARIATIONS = {"lift": ("boundary_knots", 9, 17), "forcing": ("forcing_radial_order", 2, 4)}
WITHHELD = tuple((2*i+1)/32 for i in range(16))
DOMAIN_KEYS = {"U0", "Ua", "Ub", "Va", "Vb"}
LIFT_SCHEMA = "nsc-atlas-scalar-lift-v1"
LIFT_SCOPE = "finite Hermite boundary and exact polynomial lift; withheld sensitivity is not a uniform error enclosure"
JET_DEFECTS = ("potential_identity_defect", "potential_u_identity_defect", "mass_null_identity_defect",
               "arrival_rate_defect", "wall_coordinate_rate_defect")
MARGINS = ("no_return_margin", "quiet_centre_time_margin", "initial_time_margin",
           "last_profile_margin", "first_profile_margin", "active_radius_lower")
KNOT_KEYS = ("t", "U", "V", "b", "b_U", "b_V", "U_V")
JET_KEYS = (*KNOT_KEYS, "x", "r_U", "c1_r", "J_minus", "J_minus_u", "J_minus_s", "q_minus_s",
            "q_plus_s", "K1_wall", "wall_profile", "integral_Ku_f", "positive_integral", "negative_integral",
            "ray_nodes", *JET_DEFECTS, *MARGINS)
LIFT_DEFECTS = ("b_defect", "scaled_slope_defect", "lift_derivative_defect", "wall_transport_defect")
VALIDATION_KEYS = ("t", "U", "V", "b_reference", "b_interpolated", "b_defect", "b_U_reference",
                   "b_U_interpolated", "scaled_slope_defect", "lift_derivative_defect", "wall_transport_defect", "B")
FORCING_ARRAYS = ("U", "K0", "F", "G", "mu_E", "zeta_E")
FORCING_COUNTS = ("radial_panels", "radial_nodes", "field_points", "field_batches")
FORCING_SCALARS = ("t", "V", "U_wall", "g1", "g1_V", "source_mu_E", "source_zeta_E",
    "source_H_E", "source_H_E_V", "source_chi_r", "source_chi", "source_r_U", "source_P0_r",
    "source_p_defect", "source_p_V_defect", "wall_chi_U", "chi_U_lower_limit_bound",
    "wall_integral_L_plus", "wall_mass_defect")
TARGET_RECORD = "results/development/nsc-atlas-scalar-boundary-v1/field/boundary.json"
TARGET_MANIFEST = "results/development/nsc-atlas-scalar-boundary-v1/field/manifest.json"
SOURCES = tuple(dict.fromkeys((*BOUNDARY_SOURCES,
    "docs/nsc-atlas-scalar-forcing.md", "src/recursive_horizons/nsc_atlas_scalar_lift.py",
    "src/recursive_horizons/nsc_atlas_scalar_forcing.py",
    "tests/test_nsc_atlas_scalar_lift.py", "tests/test_nsc_atlas_scalar_forcing.py",
    "scripts/record_nsc_atlas_scalar_inputs.py", "tests/test_record_nsc_atlas_scalar_inputs.py",
    TARGET_RECORD, TARGET_MANIFEST)))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def interpolation_target():
    """Diagnostic target from the nine saved field-refinement samples, not a bound."""
    raw = (ROOT/TARGET_RECORD).read_bytes()
    manifest = json.loads((ROOT/TARGET_MANIFEST).read_bytes())
    if manifest["outputs"]["boundary.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("scalar input target reference output mismatch")
    reference = json.loads(raw)
    rows = [solution["boundary_samples"] for solution in reference["solutions"]]
    if (reference["comparison"] != "field" or len(rows) != 2
            or [setting["field_nodes"] for setting in reference["settings"]] != [129, 257]
            or any([row["t"] for row in values] != list(TIMES) for values in rows)):
        raise ValueError("scalar input target requires the original nine-cut field comparison")
    def b(row):
        r_U = -.64*row["V"]*math.exp(-row["x"]/.8)/row["x"]
        return r_U*row["c1_r"]
    sensitivity = max(abs(b(left)-b(right)) for left, right in zip(*rows))
    return {"source": TARGET_RECORD, "source_sha256": digest(raw), "reference_cuts": len(TIMES),
        "sampled_field_b_sensitivity": sensitivity, "fraction": .1, "floor": 1e-11,
        "absolute_b_target": max(1e-11, .1*sensitivity), "is_uniform_bound": False}


def _prepare():
    from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, solve_prepared_field
    from recursive_horizons.nsc_atlas_wedge_field import WedgeField
    from recursive_horizons.nsc_atlas_wall_evolution import solve_wall_evolution
    from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
    shared = PreparedFieldBudget()
    metrics = BulkBudget(started=shared.started, external=shared)
    geometry = source_geometry()
    field = solve_prepared_field(BASE["field_nodes"], BASE["geometry_knots"], budget=shared)
    wedge = WedgeField(field)
    wall_settings = {key.removeprefix("wall_"): value for key, value in BASE.items() if key.startswith("wall_")}
    wall = solve_wall_evolution(wedge, geometry, budget=metrics, source_scale=1., **wall_settings)
    cap = CapLeadingField(field, **{key: BASE[key] for key in ("derivative_order", "small_radius_order", "centre_radius")})
    return shared, metrics, geometry, wedge, wall, cap


def _compare_lift(wall, cap, metrics):
    from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary
    from recursive_horizons.nsc_atlas_scalar_lift import BoundaryJet, BoundaryLift
    boundary = ScalarBoundary(cap, wall, radial_order=BASE["radial_order"], ray_order=BASE["ray_order"], budget=metrics)
    jets = BoundaryJet(boundary)
    fits = [BoundaryLift(jets, knots=count, validation_times=WITHHELD) for count in (9, 17)]
    solutions = []
    for fit in fits:
        saved = fit.summary()
        solutions.append({"schema": saved["schema"], "scope": saved["scope"], "knot_count": saved["knot_count"],
            "knot_samples": [{key: row[key] for key in KNOT_KEYS} for row in saved["knot_samples"]],
            "validation_samples": [{**{key: row[key] for key in VALIDATION_KEYS if key != "B"},
                                    "B": float(fit.B(row["U"]))} for row in saved["validation_samples"]],
            "maximum_absolute_defects": saved["maximum_absolute_defects"], "terminal_lift": saved["terminal_lift"]})
    direct = [{key: jets.cache[t][key] for key in JET_KEYS} for t in sorted(jets.cache)]
    target = interpolation_target()
    diagnostics = {"direct_samples": direct, "unique_boundary_jet_evaluations": len(direct),
        "ray_nodes": jets.ray_nodes, "interpolation_variable": "U_W", "interpolation_method": "cubic Hermite",
        "maximum_absolute_jet_defects": {key: max(abs(row[key]) for row in direct) for key in JET_DEFECTS},
        "minimum_domain_margins": {key: min(row[key] for row in direct) for key in MARGINS},
        "interpolation_target": target,
        "sampled_b_target_met": [solution["maximum_absolute_defects"]["b_defect"] <= target["absolute_b_target"]
                                  for solution in solutions],
        "maximum_lift_difference": max(abs(solutions[0]["terminal_lift"]-solutions[1]["terminal_lift"]),
            *(abs(a["B"]-b["B"]) for a, b in zip(solutions[0]["validation_samples"], solutions[1]["validation_samples"])))}
    return solutions, diagnostics


def _forcing_differences(solutions):
    result = {}
    for key in (*FORCING_ARRAYS[1:], "g1", "g1_V"):
        left, right = (np.asarray([row[key] for row in solution["samples"]]) for solution in solutions)
        result[key] = float(np.max(abs(left-right)))
    return result


def _compare_forcing(wedge, geometry, wall, metrics):
    from recursive_horizons.nsc_atlas_scalar_forcing import ScalarForcing
    wall_mass = {t: wall._wall(t)["integral_L_plus"] for t in TIMES}
    solutions = []
    for order in (2, 4):
        forcing = ScalarForcing(wedge, geometry, wedge.grid, radial_order=order,
            field_order=BASE["wall_field_order"], source_order=BASE["wall_source_order"], source_scale=1., budget=metrics)
        rows = []
        for t in TIMES:
            row = forcing.slice(t)
            row.update(wall_integral_L_plus=wall_mass[t], wall_mass_defect=float(row["mu_E"][-1])-wall_mass[t])
            rows.append(_plain(row))
        solutions.append({"schema": "nsc-atlas-scalar-forcing-v1", "settings": _plain(forcing.settings),
            "samples": rows, "slice_evaluations": forcing.slice_evaluations,
            "field_points": forcing.field_points, "field_batches": forcing.field_batches,
            "maximum_absolute_defects": {key: max(abs(row[key]) for row in rows)
                for key in ("wall_mass_defect", "source_p_defect", "source_p_V_defect", "wall_chi_U", "chi_U_lower_limit_bound")}})
    return solutions, {"requested_U_grid": _plain(wedge.grid), "maximum_absolute_differences": _forcing_differences(solutions)}


def compare(kind):
    if kind not in VARIATIONS:
        raise ValueError("scalar input comparison is not admitted")
    shared, metrics, geometry, wedge, wall, cap = _prepare()
    varied, low, high = VARIATIONS[kind]
    solutions, diagnostics = (_compare_lift(wall, cap, metrics) if kind == "lift"
                              else _compare_forcing(wedge, geometry, wall, metrics))
    shared.check()
    metrics.check()
    return _plain({"schema": SCHEMA, "comparison": kind, "status": STATUS, "scope": SCOPE,
        "source_scale": 1., "runtime": runtime_versions(), "comparison_base": BASE,
        "varied_quantity": varied, "settings": [{**BASE, varied: value} for value in (low, high)],
        "domain": {key: float(geometry[key]) for key in DOMAIN_KEYS},
        "solutions": solutions, "diagnostics": diagnostics,
        "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict(cap.evaluation_counts),
            "maximum_inverse_bracket_width": cap.max_inverse_bracket_width},
        "execution": {"cpu_seconds": time.process_time()-shared.started, "cpu_seconds_max": shared.cpu_seconds,
            "scalar_kernel_values": shared.kernel_values, "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "rhs_calls": metrics.rhs_calls, "field_solves": 1, "wall_solves": 1}})


def _real(value, label, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not math.isfinite(value) or (nonnegative and value < 0)):
        raise ValueError("invalid scalar input "+label)


def _count(value, label, maximum=None):
    if type(value) is not int or value < 0 or (maximum is not None and value > maximum):
        raise ValueError("invalid scalar input "+label)


def _equal(left, right):
    return math.isclose(left, right, rel_tol=1e-12, abs_tol=1e-13)


def _validate_lift(record):
    from scipy.interpolate import CubicHermiteSpline
    d, diagnostics = record["domain"], record["diagnostics"]
    direct = diagnostics["direct_samples"]
    if (not isinstance(direct, list) or len(direct) != 33
            or [row["t"] for row in direct] != [i/32 for i in range(33)]
            or type(diagnostics["unique_boundary_jet_evaluations"]) is not int
            or diagnostics["unique_boundary_jet_evaluations"] != 33
            or diagnostics["interpolation_variable"] != "U_W" or diagnostics["interpolation_method"] != "cubic Hermite"):
        raise ValueError("scalar lift query or interpolation inventory mismatch")
    by_time = {row["t"]: row for row in direct}
    for row in direct:
        if set(row) != set(JET_KEYS):
            raise ValueError("scalar lift direct jet inventory mismatch")
        for key, value in row.items():
            if key == "ray_nodes":
                _count(value, key, 17*BASE["ray_order"])
                if value % BASE["ray_order"]:
                    raise ValueError("scalar lift ray count mismatch")
            else:
                _real(value, "jet "+key)
        if (not d["Ua"]-1e-13 <= row["U"] <= d["Ub"]+1e-13
                or not _equal(row["V"], d["Va"]+(d["Vb"]-d["Va"])*row["t"])
                or row["U_V"] <= 0 or not .45 <= row["x"] <= .55 or row["r_U"] <= 0
                or row["positive_integral"] < 0 or row["negative_integral"] > 0
                or any(row[key] < -1e-13 for key in MARGINS)):
            raise ValueError("scalar lift jet exits the declared domain")
        identities = {"b": row["r_U"]*row["c1_r"], "b_V": row["b_U"]*row["U_V"],
            "J_minus": row["positive_integral"]+row["negative_integral"],
            "J_minus_u": -row["integral_Ku_f"],
            "J_minus_s": -row["q_plus_s"]*row["K1_wall"]*row["wall_profile"]-row["q_minus_s"]*row["integral_Ku_f"]}
        if any(not _equal(row[key], value) for key, value in identities.items()):
            raise ValueError("scalar lift saved derivative identity disagrees")
    if (not _equal(direct[0]["U"], d["Ua"]) or not _equal(direct[-1]["U"], d["Ub"])
            or any(b["U"] <= a["U"] for a, b in zip(direct[:-1], direct[1:]))
            or any(abs(direct[0][key]) > 1e-24 for key in ("b", "b_U", "b_V", "J_minus", "J_minus_s"))):
        raise ValueError("scalar lift ordering or quiet anchor lost")
    _count(diagnostics["ray_nodes"], "ray nodes")
    if diagnostics["ray_nodes"] != sum(row["ray_nodes"] for row in direct):
        raise ValueError("scalar lift cumulative ray count mismatch")
    for key, expected in (("maximum_absolute_jet_defects", {k: max(abs(row[k]) for row in direct) for k in JET_DEFECTS}),
                          ("minimum_domain_margins", {k: min(row[k] for row in direct) for k in MARGINS})):
        if diagnostics[key] != expected:
            raise ValueError("scalar lift saved diagnostics disagree")
    for solution, count in zip(record["solutions"], (9, 17)):
        if (solution["schema"] != LIFT_SCHEMA or solution["scope"] != LIFT_SCOPE
                or type(solution["knot_count"]) is not int or solution["knot_count"] != count):
            raise ValueError("scalar lift solution settings or scope mismatch")
        knots = solution["knot_samples"]
        expected = [{key: by_time[i/(count-1)][key] for key in KNOT_KEYS} for i in range(count)]
        if knots != expected:
            raise ValueError("scalar lift fitted knots are not the declared direct subset")
        U = np.array([row["U"] for row in knots])
        U[[0, -1]] = [d["Ua"], d["Ub"]]
        polynomial = CubicHermiteSpline(U, [row["b"] for row in knots], [row["b_U"] for row in knots], extrapolate=False)
        primitive = polynomial.antiderivative()
        anchor = float(primitive(d["Ua"]))
        withheld = solution["validation_samples"]
        if not isinstance(withheld, list) or len(withheld) != 16 or [row["t"] for row in withheld] != list(WITHHELD):
            raise ValueError("scalar lift withheld inventory mismatch")
        for row in withheld:
            if set(row) != set(VALIDATION_KEYS):
                raise ValueError("scalar lift withheld sample inventory mismatch")
            for value in row.values():
                _real(value, "withheld value")
            truth = by_time[row["t"]]
            index = min(int(np.searchsorted(U, row["U"], side="right")-1), count-2)
            width = U[index+1]-U[index]
            fitted, slope = float(polynomial(row["U"])), float(polynomial(row["U"], 1))
            expected = {"U": truth["U"], "V": truth["V"], "b_reference": truth["b"], "b_U_reference": truth["b_U"],
                "b_interpolated": fitted, "b_U_interpolated": slope, "b_defect": fitted-truth["b"],
                "scaled_slope_defect": width*(slope-truth["b_U"]), "B": float(primitive(row["U"]))-anchor,
                "lift_derivative_defect": float(primitive(row["U"], 1))-fitted,
                "wall_transport_defect": (float(primitive(row["U"], 1))-fitted)*truth["U_V"]}
            if any(not _equal(row[key], value) for key, value in expected.items()):
                raise ValueError("scalar lift saved interpolation or lift identity disagrees")
        _real(solution["terminal_lift"], "terminal lift")
        if not _equal(solution["terminal_lift"], float(primitive(d["Ub"]))-anchor):
            raise ValueError("scalar lift terminal antiderivative disagrees")
        if solution["maximum_absolute_defects"] != {key: max(abs(row[key]) for row in withheld) for key in LIFT_DEFECTS}:
            raise ValueError("scalar lift withheld maxima disagree")
    target = interpolation_target()
    if diagnostics["interpolation_target"] != target:
        raise ValueError("scalar lift diagnostic target reference disagrees")
    if diagnostics["sampled_b_target_met"] != [solution["maximum_absolute_defects"]["b_defect"] <= target["absolute_b_target"]
                                               for solution in record["solutions"]]:
        raise ValueError("scalar lift sampled target status disagrees")
    if any(type(value) is not bool for value in diagnostics["sampled_b_target_met"]):
        raise ValueError("scalar lift target status must be boolean")
    left, right = record["solutions"]
    expected = max(abs(left["terminal_lift"]-right["terminal_lift"]),
        *(abs(a["B"]-b["B"]) for a, b in zip(left["validation_samples"], right["validation_samples"])))
    if diagnostics["maximum_lift_difference"] != expected:
        raise ValueError("scalar lift saved comparison disagrees")


def _validate_forcing(record):
    d = record["domain"]
    grid = record["diagnostics"]["requested_U_grid"]
    if not isinstance(grid, list) or len(grid) != BASE["field_nodes"]:
        raise ValueError("scalar forcing grid inventory mismatch")
    for value in grid:
        _real(value, "grid point")
    if (grid[0] != d["U0"] or grid[-1] != d["Ub"] or any(b <= a for a, b in zip(grid[:-1], grid[1:]))
            or not np.allclose(grid, np.linspace(d["U0"], d["Ub"], BASE["field_nodes"]), rtol=0, atol=1e-15)):
        raise ValueError("scalar forcing does not use the original field grid")
    for solution, order in zip(record["solutions"], (2, 4)):
        expected = dict(radial_order=order, field_order=32, source_order=64, source_scale=1., grid_nodes=BASE["field_nodes"])
        if solution["schema"] != "nsc-atlas-scalar-forcing-v1" or solution["settings"] != expected:
            raise ValueError("scalar forcing solution settings mismatch")
        if any(type(solution["settings"][key]) is not int for key in expected if key != "source_scale"):
            raise ValueError("scalar forcing integer settings mismatch")
        _real(solution["settings"]["source_scale"], "source scale")
        rows = solution["samples"]
        if not isinstance(rows, list) or len(rows) != len(TIMES) or [row["t"] for row in rows] != list(TIMES):
            raise ValueError("scalar forcing time-cut inventory mismatch")
        for row in rows:
            if set(row) != set(FORCING_ARRAYS+FORCING_SCALARS+FORCING_COUNTS):
                raise ValueError("scalar forcing sample inventory mismatch")
            for key in FORCING_ARRAYS:
                if not isinstance(row[key], list) or len(row[key]) != len(grid):
                    raise ValueError("scalar forcing array shape mismatch")
                for value in row[key]:
                    _real(value, "array "+key)
            for key in FORCING_SCALARS:
                _real(row[key], key)
            for key in FORCING_COUNTS:
                _count(row[key], key, 513*5)
            clipped = [min(value, row["U_wall"]) for value in grid]
            clipped[-1] = row["U_wall"]
            if (not d["Ua"]-1e-13 <= row["U_wall"] <= d["Ub"]+1e-13
                    or not _equal(row["V"], d["Va"]+(d["Vb"]-d["Va"])*row["t"])
                    or row["U"] != clipped
                    or any(value <= 0 for value in row["K0"]) or row["source_r_U"] <= 0
                    or row["chi_U_lower_limit_bound"] < 0 or row["wall_chi_U"] != 0 or row["F"][-1] != 0):
                raise ValueError("scalar forcing exits its retained domain or wall condition")
            cuts = np.unique([d["U0"], *[u for u in grid if d["U0"] < u < row["U_wall"]], *clipped, row["U_wall"]])
            panels = len(cuts)-1
            nodes = np.polynomial.legendre.leggauss(order)[0]
            gauss = cuts[:-1, None]+(np.diff(cuts)/2)[:, None]*(1+nodes)
            points = len(np.unique(np.r_[cuts, gauss.ravel()]))
            if (row["radial_panels"] != panels or row["radial_nodes"] != order*panels
                    or row["field_points"] != points
                    or row["field_batches"] != (row["field_points"]+512)//513):
                raise ValueError("scalar forcing quadrature accounting mismatch")
            g1V = -(row["source_H_E_V"]+row["source_H_E"]*row["source_P0_r"]/2)*row["source_chi_r"]
            g1V += row["source_H_E"]*row["K0"][0]*row["source_chi"]/row["source_r_U"]
            if (not _equal(row["g1"], -row["source_H_E"]*row["source_chi_r"])
                    or not _equal(row["g1_V"], g1V)
                    or row["mu_E"][0] != row["source_mu_E"] or row["zeta_E"][0] != row["source_zeta_E"]
                    or not _equal(row["wall_mass_defect"], row["mu_E"][-1]-row["wall_integral_L_plus"])):
                raise ValueError("scalar forcing source or wall identity disagrees")
            for key in FORCING_ARRAYS[1:]:
                wall_value = row[key][-1]
                if any(value != wall_value for u, value in zip(row["U"], row[key]) if u == row["U_wall"]):
                    raise ValueError("scalar forcing clipped continuation disagrees")
        if any(abs(rows[0][key]) > 1e-24 for key in ("g1", "g1_V", "wall_integral_L_plus")) or any(
                abs(value) > 1e-24 for key in ("F", "G", "mu_E", "zeta_E") for value in rows[0][key]):
            raise ValueError("scalar forcing quiet anchor lost")
        for key, expected in (("slice_evaluations", len(TIMES)), ("field_points", sum(row["field_points"] for row in rows)),
                              ("field_batches", sum(row["field_batches"] for row in rows))):
            if type(solution[key]) is not int or solution[key] != expected:
                raise ValueError("scalar forcing cumulative count mismatch")
        expected = {key: max(abs(row[key]) for row in rows)
                    for key in ("wall_mass_defect", "source_p_defect", "source_p_V_defect", "wall_chi_U", "chi_U_lower_limit_bound")}
        if solution["maximum_absolute_defects"] != expected:
            raise ValueError("scalar forcing diagnostic maxima disagree")
    for left, right in zip(*(solution["samples"] for solution in record["solutions"])):
        if any(left[key] != right[key] for key in ("U", "V", "K0", "U_wall", "wall_integral_L_plus")):
            raise ValueError("scalar forcing comparison coordinates differ")
    if record["diagnostics"]["maximum_absolute_differences"] != _forcing_differences(record["solutions"]):
        raise ValueError("scalar forcing saved differences disagree")


def validate_record(record, kind):
    try:
        if kind not in VARIATIONS or not isinstance(record, dict):
            raise ValueError("scalar input declaration mismatch")
        expected = {"schema": SCHEMA, "comparison": kind, "status": STATUS, "scope": SCOPE}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError("scalar input declaration mismatch")
        if (record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("scalar input source or non-enclosure mismatch")
        varied, low, high = VARIATIONS[kind]
        settings = [{**BASE, varied: value} for value in (low, high)]
        if record["varied_quantity"] != varied or record["settings"] != settings or record["comparison_base"] != BASE:
            raise ValueError("scalar input settings mismatch")
        for setting in record["settings"]:
            for key, value in setting.items():
                if key in ("wall_tolerance", "centre_radius"):
                    _real(value, key)
                elif type(value) is not int:
                    raise ValueError("scalar input integer setting mismatch")
        runtime = record["runtime"]
        if (not isinstance(runtime, dict) or set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("scalar input runtime mismatch")
        domain = record["domain"]
        if set(domain) != DOMAIN_KEYS:
            raise ValueError("scalar input domain inventory mismatch")
        for value in domain.values():
            _real(value, "domain")
        if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0 and domain["Va"] < domain["Vb"] < 0
                and domain["U0"]*domain["Va"] < 1):
            raise ValueError("scalar input domain exits the retained wedge")
        if not isinstance(record["solutions"], list) or len(record["solutions"]) != 2:
            raise ValueError("scalar input comparison requires two solutions")
        (_validate_lift if kind == "lift" else _validate_forcing)(record)
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("scalar input profile inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        counts = profile["counts"]
        if not counts["second_points"] <= counts["wall_points"] <= counts["profile_requested"] or counts["cache_hits"] > counts["profile_requested"]:
            raise ValueError("scalar input profile count mismatch")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(domain["Va"]), abs(domain["Vb"])):
            raise ValueError("scalar input inverse bracket exceeds its numerical tolerance")
        cost = record["execution"]
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("scalar input budget ceilings changed")
        _real(cost["cpu_seconds"], "CPU time", True)
        if cost["cpu_seconds"] > 120:
            raise ValueError("scalar input CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("scalar input metric budget exceeded")
        _count(cost["rhs_calls"], "RHS calls", 12000)
        if any(type(cost[key]) is not int or cost[key] != 1 for key in ("field_solves", "wall_solves")):
            raise ValueError("scalar input solve count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete scalar input record") from error
    return True


def create(directory, kind):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    if kind not in VARIATIONS:
        raise ValueError("scalar input comparison is not admitted")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare(kind)
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record, kind)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings, "comparison": kind,
        "outputs": {"inputs.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"inputs.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    diagnostics = record["diagnostics"]
    summary = ({"withheld_maxima": [solution["maximum_absolute_defects"] for solution in record["solutions"]],
        "sampled_b_target_met": diagnostics["sampled_b_target_met"],
        "absolute_b_target": diagnostics["interpolation_target"]["absolute_b_target"],
        "maximum_lift_difference": diagnostics["maximum_lift_difference"]} if kind == "lift" else
        {"maximum_absolute_differences": diagnostics["maximum_absolute_differences"],
         "maximum_absolute_defects": [solution["maximum_absolute_defects"] for solution in record["solutions"]]})
    return {"created": str(directory), "comparison": kind,
            "execution": record["execution"], "summary": summary}


def check(directory):
    directory = Path(directory)
    manifest_path = directory/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    if (manifest_path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings() or manifest.get("comparison") not in VARIATIONS):
        raise ValueError("scalar input source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "inputs.json"}
            or set(manifest.get("outputs", {})) != {"inputs.json"}):
        raise ValueError("scalar input inventory mismatch")
    path = directory/"inputs.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["inputs.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("scalar input output mismatch")
    validate_record(json.loads(raw), manifest["comparison"])
    return {"checked": str(directory), "comparison": manifest["comparison"], "physical_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", choices=VARIATIONS, required=False, default="lift")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output, args.comparison), sort_keys=True))


if __name__ == "__main__":
    main()
