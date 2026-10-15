#!/usr/bin/env python3
"""Compare profile-moment orders with the original radial boundary and mass.

All preparation and independent reference queries share one finite budget.
Authentication checks saved bytes and identities without physical replay.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from record_nsc_atlas_scalar_inputs import (
    SOURCES as INPUTS_SOURCES, BASE, TIMES, JET_KEYS, JET_DEFECTS, MARGINS,
    PROFILE_COUNTS, _prepare, _plain, _real, _count, _equal, digest, encoded, runtime_versions,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-profile-moment-comparison-v1"
SCOPE = "formal pre-arrival profile-moment mass and boundary inputs; no refined boundary interpolant, evolved correction or uniform remainder"
STATUS = "FORMAL_PROFILE_MOMENT_REFINEMENT"
ORDERS = (4, 8)
FRACTIONS = (1/3, 2/3, 1.)
QUANTITIES = ("J_minus", "J_minus_s", "b", "b_U")
SAMPLE_KEYS = (*JET_KEYS, "q_minus", "q_plus")
DOMAIN_KEYS = {"U0", "Ua", "Ub", "Va", "Vb", "qa", "qb", "Ta", "Tb", "qmin", "xa", "xb"}
SOURCES = tuple(dict.fromkeys((*INPUTS_SOURCES,
    "docs/nsc-atlas-profile-moment.md", "src/recursive_horizons/nsc_atlas_profile_moment.py",
    "tests/test_nsc_atlas_profile_moment.py", "scripts/record_nsc_atlas_profile_moment.py",
    "tests/test_record_nsc_atlas_profile_moment.py")))


def source_bindings():
    return {name: digest((ROOT/name).read_bytes()) for name in SOURCES}


def _differences(left, right):
    return {key: max(abs(a[key]-b[key]) for a, b in zip(left, right)) for key in QUANTITIES}


def diagnostics(reference, solutions):
    return {"boundary_order_differences": _differences(solutions[0]["boundary_samples"], solutions[1]["boundary_samples"]),
        "boundary_radial_differences": [_differences(solution["boundary_samples"], reference["boundary_samples"]) for solution in solutions],
        "mass_order_difference": max(abs(a["nu_E"]-b["nu_E"]) for a, b in zip(solutions[0]["mass_samples"], solutions[1]["mass_samples"])),
        "mass_radial_differences": [max(abs(a["nu_E"]-b["nu_E"]) for a, b in zip(solution["mass_samples"], reference["mass_samples"])) for solution in solutions]}


def compare():
    from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary
    from recursive_horizons.nsc_atlas_scalar_lift import BoundaryJet
    from recursive_horizons.nsc_atlas_cap_response import radial_metric
    from recursive_horizons.nsc_atlas_profile_moment import ProfileMoment, MomentBoundaryJet
    shared, metrics, _, _, wall, cap = _prepare()
    boundary = ScalarBoundary(cap, wall, radial_order=BASE["radial_order"], ray_order=BASE["ray_order"], budget=metrics)
    original = BoundaryJet(boundary)
    reference_rows = [original.sample(t) for t in TIMES]
    points, reference_mass = [], []
    for row in reference_rows[1:]:
        u, qW = row["q_minus"], row["q_plus"]
        for fraction in FRACTIONS:
            q = cap.domain["qa"]+fraction*(qW-cap.domain["qa"])
            T, r = (u+q)/2, math.tanh((q-u)/2)
            point = dict(t=row["t"], fraction=fraction, q_minus=u, q=q, T=T, r=r)
            points.append(point)
            result = radial_metric(cap, T, r, radial_order=BASE["radial_order"], budget=metrics)
            reference_mass.append({**point, "nu_E": float(result["nu_E"])})
    reference = {"method": "original radial metric and joint boundary jet", "radial_order": BASE["radial_order"],
        "ray_order": BASE["ray_order"], "boundary_samples": [{key: row[key] for key in SAMPLE_KEYS} for row in reference_rows],
        "boundary_evaluations": len(original.cache), "ray_nodes": original.ray_nodes,
        "mass_samples": reference_mass, "mass_evaluations": len(points)}
    solutions = []
    for order in ORDERS:
        moment = ProfileMoment(cap, order=order, budget=metrics)
        jet = MomentBoundaryJet(boundary, moment_order=order, moment=moment)
        rows = [jet.sample(t) for t in TIMES]
        masses = []
        for point in points:
            result = moment.mass(point["T"], point["r"])
            masses.append({**point, **{key: float(result[key]) for key in ("nu_E", "f", "C")}})
        solutions.append({"moment_order": order, "boundary_samples": [{key: row[key] for key in SAMPLE_KEYS} for row in rows],
            "boundary_evaluations": len(jet.cache), "ray_nodes": jet.ray_nodes,
            "mass_samples": masses, "mass_evaluations": len(masses), "moment": moment.summary()})
    shared.check()
    metrics.check()
    return _plain({"schema": SCHEMA, "status": STATUS, "scope": SCOPE, "runtime": runtime_versions(),
        "comparison_base": BASE, "varied_quantity": "moment_order", "settings": [{**BASE, "moment_order": order} for order in ORDERS],
        "source_scale": 1., "domain": dict(cap.domain), "reference": reference, "solutions": solutions,
        "diagnostics": diagnostics(reference, solutions), "refinement_is_rigorous_error_enclosure": False,
        "profile_evaluation": {"counts": dict(cap.evaluation_counts), "maximum_inverse_bracket_width": cap.max_inverse_bracket_width},
        "execution": {"cpu_seconds": time.process_time()-shared.started, "cpu_seconds_max": shared.cpu_seconds,
            "scalar_kernel_values": shared.kernel_values, "scalar_kernel_values_max": shared.kernel_values_max,
            "metric_values": metrics.metric_values, "source_quadrature_nodes": metrics.source_nodes,
            "rhs_calls": metrics.rhs_calls, "field_solves": 1, "wall_solves": 1, "moment_evaluators": 2}})


def _validate_samples(solution, domain, moment=False):
    rows = solution["boundary_samples"]
    if not isinstance(rows, list) or len(rows) != len(TIMES) or [row["t"] for row in rows] != list(TIMES):
        raise ValueError("profile moment boundary cut inventory mismatch")
    for row in rows:
        if set(row) != set(SAMPLE_KEYS):
            raise ValueError("profile moment boundary sample inventory mismatch")
        for key, value in row.items():
            if key == "ray_nodes":
                _count(value, key, 17*BASE["ray_order"])
                if value % BASE["ray_order"]:
                    raise ValueError("profile moment ray count mismatch")
            else:
                _real(value, key)
        if (not domain["Ua"]-1e-13 <= row["U"] <= domain["Ub"]+1e-13
                or not _equal(row["V"], domain["Va"]+(domain["Vb"]-domain["Va"])*row["t"])
                or row["U_V"] <= 0 or row["r_U"] <= 0 or not .45 <= row["x"] <= .55
                or not domain["qa"]-1e-13 <= row["q_plus"] <= domain["qb"]+1e-13
                or not row["q_minus"] < domain["qa"] or any(row[key] < -1e-13 for key in MARGINS)
                or row["positive_integral"] < 0 or row["negative_integral"] > 0):
            raise ValueError("profile moment boundary leaves its pre-arrival domain")
        identities = {"b": row["r_U"]*row["c1_r"], "b_V": row["b_U"]*row["U_V"],
            "q_minus": row["q_plus"]-2*math.atanh(row["x"]),
            "J_minus": row["positive_integral"]+row["negative_integral"], "J_minus_u": -row["integral_Ku_f"],
            "J_minus_s": -row["q_plus_s"]*row["K1_wall"]*row["wall_profile"]-row["q_minus_s"]*row["integral_Ku_f"]}
        identities.update(no_return_margin=domain["qa"]-row["q_minus"],
            quiet_centre_time_margin=domain["qa"]-(row["q_minus"]+row["q_plus"])/2,
            initial_time_margin=(row["q_minus"]+domain["qa"])/2-domain["Ta"],
            last_profile_margin=domain["qb"]-row["q_plus"], first_profile_margin=row["q_minus"]-domain["qmin"],
            active_radius_lower=math.tanh((domain["qa"]-row["q_minus"])/2))
        if any(not _equal(row[key], value) for key, value in identities.items()):
            raise ValueError("profile moment boundary derivative identity disagrees")
    if (not _equal(rows[0]["q_plus"], domain["qa"]) or not _equal(rows[-1]["q_plus"], domain["qb"])
            or any(b["U"] <= a["U"] or b["q_plus"] <= a["q_plus"] for a, b in zip(rows[:-1], rows[1:]))
            or any(abs(rows[0][key]) > 1e-24 for key in (*QUANTITIES, "b_V"))):
        raise ValueError("profile moment ordered cuts or quiet anchor lost")
    if type(solution["boundary_evaluations"]) is not int or solution["boundary_evaluations"] != len(TIMES):
        raise ValueError("profile moment boundary count mismatch")
    _count(solution["ray_nodes"], "total ray nodes")
    if solution["ray_nodes"] != sum(row["ray_nodes"] for row in rows):
        raise ValueError("profile moment cumulative ray count mismatch")
    masses = solution["mass_samples"]
    if (not isinstance(masses, list) or len(masses) != 24
            or type(solution["mass_evaluations"]) is not int or solution["mass_evaluations"] != 24):
        raise ValueError("profile moment mass sample count mismatch")
    keys = {"t", "fraction", "q_minus", "q", "T", "r", "nu_E"} | ({"f", "C"} if moment else set())
    for point, (row, fraction) in zip(masses, ((row, fraction) for row in rows[1:] for fraction in FRACTIONS)):
        if set(point) != keys:
            raise ValueError("profile moment mass sample inventory mismatch")
        for value in point.values():
            _real(value, "mass sample")
        q = domain["qa"]+fraction*(row["q_plus"]-domain["qa"])
        u = row["q_minus"]
        if (point["t"] != row["t"] or point["fraction"] != fraction or point["q_minus"] != u
                or not _equal(point["q"], q) or not _equal(point["T"], (u+q)/2)
                or not _equal(point["r"], math.tanh((q-u)/2))
                or not domain["Ta"]-1e-13 <= point["T"] <= domain["Tb"]+1e-13
                or point["T"] >= domain["qa"] or not 0 < point["r"] < 1
                or (not moment and point["nu_E"] < 0)):
            raise ValueError("profile moment mass point leaves its declared ray")
        if moment and (point["C"] < 0 or not _equal(point["nu_E"], point["C"]-(1-point["r"]**2)*point["f"]**2/(4*point["r"]))):
            raise ValueError("profile moment mass identity disagrees")


def _validate_moment_summary(summary, order, domain, ray_nodes):
    counts = ("cells", "prefix_nodes", "partial_nodes", "profile_points", "profile_batches",
              "value_queries", "cache_hits", "cached_values", "mass_points")
    if (set(summary) != {"schema", "domain", "settings", "terminal_C", *counts}
            or summary["schema"] != "nsc-atlas-profile-moment-v1" or summary["domain"] != domain
            or summary["settings"] != {"moment_order": order}
            or type(summary["settings"]["moment_order"]) is not int):
        raise ValueError("profile moment summary declarations mismatch")
    for key in counts:
        _count(summary[key], key, 30_000_000)
    _real(summary["terminal_C"], "terminal moment", True)
    if (not 1 <= summary["cells"] <= 1024 or summary["prefix_nodes"] != order*summary["cells"]
            or summary["partial_nodes"] % order or summary["partial_nodes"] > order*summary["value_queries"]
            or summary["mass_points"] != ray_nodes+len(TIMES)+24
            or summary["value_queries"] != summary["mass_points"]
            or summary["cached_values"] > min(32768, summary["value_queries"])
            or summary["cache_hits"] > summary["value_queries"]
            or summary["profile_points"] > summary["prefix_nodes"]+summary["partial_nodes"]+summary["mass_points"]
            or not (summary["profile_points"]+512)//513 <= summary["profile_batches"] <= summary["profile_points"]):
        raise ValueError("profile moment cumulative accounting mismatch")


def validate_record(record):
    try:
        if any(record.get(key) != value for key, value in {"schema": SCHEMA, "status": STATUS, "scope": SCOPE}.items()):
            raise ValueError("profile moment scope or schema mismatch")
        if (record.get("refinement_is_rigorous_error_enclosure") is not False
                or record.get("source_scale") != 1. or isinstance(record.get("source_scale"), bool)):
            raise ValueError("profile moment source or non-enclosure mismatch")
        if (record["comparison_base"] != BASE or record["varied_quantity"] != "moment_order"
                or record["settings"] != [{**BASE, "moment_order": order} for order in ORDERS]):
            raise ValueError("profile moment settings mismatch")
        for setting in record["settings"]:
            for key, value in setting.items():
                if key in ("wall_tolerance", "centre_radius"):
                    _real(value, key)
                elif type(value) is not int:
                    raise ValueError("profile moment integer setting mismatch")
        runtime = record["runtime"]
        if (set(runtime) != {"python", "numpy", "scipy"}
                or any(not isinstance(v, str) or not v or not v[0].isdigit() or "." not in v for v in runtime.values())):
            raise ValueError("profile moment runtime mismatch")
        domain = record["domain"]
        if set(domain) != DOMAIN_KEYS:
            raise ValueError("profile moment domain inventory mismatch")
        for value in domain.values():
            _real(value, "domain")
        if not (domain["U0"] < domain["Ua"] < domain["Ub"] < 0 and domain["Va"] < domain["Vb"] < 0
                and domain["U0"]*domain["Va"] < 1 and domain["Ta"] < domain["Tb"] < domain["qa"] < domain["qb"]
                and _equal(domain["xa"], .45) and _equal(domain["xb"], .55)
                and _equal(domain["Ta"], domain["qa"]-math.atanh(domain["xa"]))
                and _equal(domain["Tb"], domain["qb"]-math.atanh(domain["xb"]))
                and _equal(domain["qmin"], domain["Ta"]-math.atanh(domain["xa"]))):
            raise ValueError("profile moment domain exits the original pre-arrival cap")
        reference = record["reference"]
        if (reference["method"] != "original radial metric and joint boundary jet"
                or type(reference["radial_order"]) is not int or reference["radial_order"] != 16
                or type(reference["ray_order"]) is not int or reference["ray_order"] != 16):
            raise ValueError("profile moment independent reference settings mismatch")
        _validate_samples(reference, domain)
        solutions = record["solutions"]
        if not isinstance(solutions, list) or len(solutions) != 2:
            raise ValueError("profile moment comparison requires two solutions")
        for solution, order in zip(solutions, ORDERS):
            if type(solution["moment_order"]) is not int or solution["moment_order"] != order:
                raise ValueError("profile moment order mismatch")
            _validate_samples(solution, domain, moment=True)
            for row, ref in zip(solution["boundary_samples"], reference["boundary_samples"]):
                if any(row[key] != ref[key] for key in ("t", "U", "V", "x", "q_minus", "q_plus", "U_V", "r_U")):
                    raise ValueError("profile moment reference coordinates differ")
            _validate_moment_summary(solution["moment"], order, domain, solution["ray_nodes"])
        if record["diagnostics"] != diagnostics(reference, solutions):
            raise ValueError("profile moment saved differences disagree")
        for key in ("boundary_order_differences",):
            for value in record["diagnostics"][key].values():
                _real(value, "comparison difference", True)
        for row in record["diagnostics"]["boundary_radial_differences"]:
            for value in row.values():
                _real(value, "reference difference", True)
        _real(record["diagnostics"]["mass_order_difference"], "mass difference", True)
        for value in record["diagnostics"]["mass_radial_differences"]:
            _real(value, "reference mass difference", True)
        profile = record["profile_evaluation"]
        if set(profile["counts"]) != PROFILE_COUNTS:
            raise ValueError("profile moment profile count inventory mismatch")
        for key, value in profile["counts"].items():
            _count(value, key)
        counts = profile["counts"]
        if counts["second_points"] != 0 or not counts["wall_points"] <= counts["profile_requested"] or counts["cache_hits"] > counts["profile_requested"]:
            raise ValueError("profile moment uses undeclared second derivatives or inconsistent counts")
        _real(profile["maximum_inverse_bracket_width"], "inverse bracket", True)
        if profile["maximum_inverse_bracket_width"] > 4*np.finfo(float).eps*max(abs(domain["Va"]), abs(domain["Vb"])):
            raise ValueError("profile moment inverse bracket exceeds its numerical tolerance")
        cost = record["execution"]
        _real(cost["cpu_seconds_max"], "CPU ceiling")
        _count(cost["scalar_kernel_values_max"], "kernel ceiling")
        if cost["cpu_seconds_max"] != 120 or cost["scalar_kernel_values_max"] != 30_000_000:
            raise ValueError("profile moment budget ceilings changed")
        _real(cost["cpu_seconds"], "CPU time", True)
        if cost["cpu_seconds"] > 120:
            raise ValueError("profile moment CPU budget exceeded")
        for key in ("scalar_kernel_values", "metric_values", "source_quadrature_nodes"):
            _count(cost[key], key, 30_000_000)
        if cost["metric_values"]+cost["source_quadrature_nodes"] > 30_000_000:
            raise ValueError("profile moment metric budget exceeded")
        _count(cost["rhs_calls"], "RHS calls", 12000)
        for key, value in (("field_solves", 1), ("wall_solves", 1), ("moment_evaluators", 2)):
            if type(cost[key]) is not int or cost[key] != value:
                raise ValueError("profile moment evaluator count mismatch")
    except (KeyError, TypeError, IndexError) as error:
        raise ValueError("incomplete profile moment record") from error
    return True


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a fresh directory or --check")
    bindings, runtime = source_bindings(), runtime_versions()
    record = compare()
    if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
        raise RuntimeError("calculation inputs or runtime changed during comparison")
    validate_record(record)
    raw = encoded(record)
    manifest = {"schema": SCHEMA+"-manifest", "source_bindings": bindings,
        "outputs": {"moment.json": {"sha256": digest(raw), "bytes": len(raw)}}}
    directory.mkdir(parents=True, exist_ok=False)
    (directory/"moment.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return {"created": str(directory), "execution": record["execution"], "diagnostics": record["diagnostics"]}


def check(directory):
    directory = Path(directory)
    manifest_path = directory/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    if (manifest_path.is_symlink() or manifest.get("schema") != SCHEMA+"-manifest"
            or manifest.get("source_bindings") != source_bindings()):
        raise ValueError("profile moment source binding or schema mismatch")
    if ({p.name for p in directory.iterdir()} != {"manifest.json", "moment.json"}
            or set(manifest.get("outputs", {})) != {"moment.json"}):
        raise ValueError("profile moment inventory mismatch")
    path = directory/"moment.json"
    raw = path.read_bytes()
    if path.is_symlink() or manifest["outputs"]["moment.json"] != {"sha256": digest(raw), "bytes": len(raw)}:
        raise ValueError("profile moment output mismatch")
    validate_record(json.loads(raw))
    return {"checked": str(directory), "physical_solves": 0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
