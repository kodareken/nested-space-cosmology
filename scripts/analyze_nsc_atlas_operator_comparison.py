"""Same-grid saved-state comparison of two exterior departure operators.

The baseline's bound arrays are referenced, never copied or evolved again.
Diagnostic derivatives remain independent of either evolution gradient.
"""
from __future__ import annotations

import copy
import numpy as np

import analyze_nsc_atlas_spatial_refinement as spatial

inherited = spatial.inherited
BASELINE_PATH = "results/development/nsc-atlas-spatial-refinement-v1/spatial.json"


def compare_constraints(left, right, budget, *, time_sampling=False):
    """Keep same-state time sampling distinct from equal-width method changes."""
    if type(time_sampling) is not bool:
        raise ValueError("time_sampling must be boolean")
    result = spatial.compare_constraints(left, right, budget,
                                        same_spatial_stencil=time_sampling)
    result.pop("same_spatial_stencil")
    result.update(comparison_kind="same-state time sampling" if time_sampling else "same-grid operator change",
        identical_pure_U_samples_required=time_sampling,
        same_spatial_diagnostic_widths=True,
        interpretation="signed sampled difference; no convergence order or error enclosure")
    return result


def diagnostics(record, baseline64, base513, cache1025, profile, budget):
    """New zero/finite64-cell windows versus the frozen64-cell baseline.

    base513 supplies the original saved-cut helper inputs; cache1025 supplies
    the exact evolution keys. Cap comparisons use TOTAL physical integrals
    directly, since both methods have the same normalized radial cells.
    No source, cache, old stage, old trajectory or restricted32 calculation
    occurs here. The one analytic vacuum wall query is shared between cases.
    """
    budget.check()
    initial = budget.evaluations
    if len(record["cases"]) != 2 or len(baseline64["cases"]) != 2:
        raise ValueError("ordered zero/finite cases required")
    coordinates = spatial.array(cache1025["x_nodes"], (1025,), "cache coordinates", budget)
    x = coordinates[::8]
    gauge = inherited.subset(base513["gauge"]["rows"], slice(None, None, 4))
    spatial.same(gauge["x"], x, "base/cache saved cuts")
    for key in ("x", "W", "W_x", "W_xx"):
        spatial.same(gauge[key], np.asarray(cache1025["gauge_rows"][key])[::8], "base/cache gauge "+key)
    W = np.asarray(gauge["W"])
    vacuum = spatial.vacuum_slice(W, x, budget=budget)
    histories, transports, tables, old_tables, caps, controls, wall_cases = [], [], [], [], [], [], []
    for ci, label in enumerate(("zero", "finite")):
        case, old = record["cases"][ci], baseline64["cases"][ci]
        for current in (case, old):
            if (current["amplitude"], current["cells"], current["steps"]) != (label, 64, 512):
                raise ValueError("require the same64-cell512-step cases")
            spatial.same([s["x"] for s in current["snapshots"]], x, "saved cuts")
        if case["alpha"] != old["alpha"] or case["alpha"] != base513["initial_data"][label]["alpha"]:
            raise ValueError("amplitude must remain fixed")
        rows = base513["initial_data"][label]["snapshot_rows"]
        source = inherited.subset(rows, slice(None, None, 4))
        cached = inherited.subset(cache1025["source_rows"][label], slice(None, None, 8))
        for key in ("x", "r", "r_x", "r_xx", "log_C", "log_C_x", "c", "c_x"):
            spatial.same(source[key], cached[key], "base/cache source "+key)
        table = inherited.tables(case, rows, budget)
        old_table = inherited.tables(old, rows, budget)
        tables.append(table)
        old_tables.append(old_table)
        history, grids, method = {}, {}, {}
        old_grids = baseline64["diagnostics"]["constraint_transport"][ci]["grids"]
        for nt, stride in ((129, 1), (65, 2)):
            raw, departed = inherited._histories(base513, case, budget, stride)
            history["raw"+str(nt)] = inherited.compact_history(raw)
            history["departure"+str(nt)] = inherited.compact_history(departed)
            src = {key: np.asarray(value)[::stride] for key, value in source.items()}
            grid = spatial.constraint_grid(departed, x[::stride], np.linspace(0, 1, 65),
                W[::stride], table["exterior"]["r"][::stride], src, case["alpha"],
                base513["geometry"]["support_x"], budget)
            grids[str(nt)] = grid
            method[str(nt)] = compare_constraints(grid, old_grids["native64_"+str(nt)], budget)
        sampling = compare_constraints(grids["65"],
            spatial.select_grid(grids["129"], rows=slice(1, None, 2)), budget, time_sampling=True)
        histories.append(history)
        transports.append(dict(amplitude=label, grids=grids, sampling=sampling,
            method_difference=method,
            baseline_grid_paths={str(nt): f"diagnostics.constraint_transport[{ci}].grids.native64_{nt}"
                                 for nt in (129, 65)},
            baseline_history_path=f"diagnostics.histories[{ci}].native64"))
        caps.append(inherited.cap_history(case, profile, budget))
        quiet = x <= base513["geometry"]["support_x"][0]
        channels = {group+"."+key: inherited.statistics(np.asarray(table[group][key])[quiet])
                    for group, key in (("exterior", "c"), ("wall", "c"), ("cap", "Ip"), ("cap", "Iq"))}
        budget.check(sum(np.asarray(table[group][key])[quiet].size for group, key in
                        (("exterior", "c"), ("wall", "c"), ("cap", "Ip"), ("cap", "Iq"))))
        controls.append(dict(quiet_scalar_channels=channels,
            minimum_upper_label_margin=min(caps[-1]["arrays"]["upper_label_margin"]),
            final_wall_kinematic_defect=case["final"]["wall_kinematic_defect"],
            final_canonical_wall_flux_defect=case["final"]["canonical_wall_flux_defect"],
            final_conformal_wall_flux_defect=case["final"]["conformal_wall_flux_defect"],
            constructed_identities_are_accuracy_tests=False))
        Z, wall = spatial.wall_departures(case, source, vacuum, budget)
        old_wall = baseline64["diagnostics"]["wall_curvature"]["cases"][ci]
        if old_wall["amplitude"] != label:
            raise ValueError("baseline wall case order mismatch")
        estimates, wall_method, widths = {}, {}, {}
        for points in (3, 5):
            for stride in (1, 2):
                key = f"p{points}_s{stride}"
                estimate = spatial.wall_estimate(Z, wall, vacuum, W, case["alpha"],
                    points=points, stride=stride, budget=budget)
                prior = old_wall["estimates"]["fine64"][key]
                for name in ("points", "stride", "eta_spacing", "first_indices", "second_indices"):
                    spatial.same(estimate[name], prior[name], "same wall stencil "+name)
                estimates[key] = estimate
                wall_method[key] = spatial.wall_comparison(estimate, prior, x, budget)
            widths[str(points)] = spatial.wall_comparison(estimates[f"p{points}_s1"],
                                                          estimates[f"p{points}_s2"], x, budget)
        wall_cases.append(dict(amplitude=label, estimates=estimates, method_difference=wall_method,
            native_minus_wide=widths,
            baseline_estimates_path=f"diagnostics.wall_curvature.cases[{ci}].estimates.fine64"))
    paired = inherited._difference(tables[1], tables[0])
    old_paired = inherited._difference(old_tables[1], old_tables[0])
    oracle = copy.deepcopy(record["cases"][0])
    for row in oracle["snapshots"]:
        row["state"]["cap"] = row["physical"]["cap"]
    result = dict(baseline_reference=dict(path=BASELINE_PATH,
            arrays="authenticated baseline arrays are referenced, not copied"),
        histories=histories, constraint_transport=transports, cap_histories=caps, controls=controls,
        zero_cap_pullback=inherited.cap_pullback(oracle, budget),
        finite_minus_zero64=inherited.comparison(tables[1], tables[0], x, budget),
        method_zero=inherited.comparison(tables[0], old_tables[0], x, budget),
        method_finite=inherited.comparison(tables[1], old_tables[1], x, budget),
        paired_method=inherited.comparison(paired, old_paired, x, budget),
        wall_curvature=dict(x=x, W=W, cases=wall_cases, analytic_vacuum_points=len(x),
            birth_transport_evaluations=0, scope="independent saved curvature at identical diagnostic widths"),
        comparison_coordinates="same x, exterior eta and cap normalized cells; physical cap radii may differ",
        cap_comparison="direct TOTAL physical dr integrals; no cell summation, interpolation or raw-departure substitution",
        source_constraint_is_zero=False, spatial_order_established=False, uniform_error_enclosure=False,
        old_evolutions=0, restricted32_diagnostics=0)
    budget.check()
    result["diagnostic_evaluations"] = budget.evaluations-initial
    return inherited.plain(result)
