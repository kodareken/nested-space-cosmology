#!/usr/bin/env python3
"""Evaluate a fixed IV preparation through the wall to a finite interior signal."""
from __future__ import annotations

import argparse
from copy import deepcopy
from importlib.metadata import version
from io import BytesIO
import json
from pathlib import Path
import platform

import numpy as np

import derive_nsc_cosmological_atlas_v7 as previous
from recursive_horizons import nsc_atlas_prepared_response as response_model

first = previous.first
REPO = previous.REPO
SCHEMA = "nsc-cosmological-atlas-v8"
DEFAULT = REPO / "results/development/nsc-cosmological-atlas-v8"
OWNER = "docs/nsc-atlas-prepared-response.md"
EVIDENCE = "results/development/nsc-cosmological-atlas-v8/prepared-response.json"


def build_record():
    record, programme, prepared = previous.build_record()
    response = response_model.prepared_response_record()
    scenario = deepcopy(record["scenarios"][-1])
    scenario.update(id="finite_iv_prepared_response",
        title="An independently prepared exterior signal at the interior centre",
        solution_status="finite source-to-observer waveform evaluated with refinement and a positive post-input tail",
        field_record="prepared-response.json")
    scenario["assumptions"] = [item for item in scenario["assumptions"]
        if "no global parent source or waveform evaluation" not in item]
    scenario["assumptions"] += [
        "fixed classical linear probe and the same v7 characteristic input and zero cap at R/ell=0.45",
        "finite numerical refinement is a diagnostic, not a continuum error enclosure",
        "positive post-input field tail only within the declared window; later evolution remains uncomputed",
    ]
    for region in scenario["regions"]:
        if region["id"] == "C":
            region["state"] = "zero cap at wall R/ell=0.45 followed by the evaluated response to independent IV data within the recorded proper-time window"
        elif region["id"] == "W":
            region["state"] = "evaluated trace from the unchanged independent IV preparation; no freely adjusted wall pulse"
    record["sources"]["finite_iv_response"] = {
        "path": OWNER, "use": "fixed prepared scalar transfer, numerical refinement and analytic finite-tail sign",
        "limit": "prescribed IV data; no remote source generation, quantum state, permanent memory or backreaction"}
    for link in scenario["links"]:
        link["source_refs"] = sorted(set(link["source_refs"] + ["finite_iv_response"]))
        if link["id"] == "preparation-IV-W":
            link.update(evidence="prepared-response.json evaluates the retained-IV Volterra operator for the fixed g(V)",
                status="conditional_prepared_field_transfer",
                missing=["physical generation of the independent IV data", "later returns and backreaction"])
        elif link["id"] == "scalar-W-C":
            link.update(evidence="prepared-response.json evaluates phi_C[T(tau)]=I(tau)/Tprime(tau) for the exterior-derived wall trace",
                status="conditional_prepared_observer_response",
                missing=["detector dynamics and energy flux", "evolution beyond the finite clock window"])
    preparation_claim = scenario["connection_records"][0]
    preparation_claim.update(id="finite-IV-prepared-response",
        missing=["physical generation of the independent IV data", "global quantum state and later evolution"],
        waveform_evaluated=True)
    preparation_claim["owner_refs"].append(OWNER)
    preparation_claim["evidence_refs"].append(EVIDENCE)
    preparation_claim["independent_check_refs"].append("tests/test_nsc_atlas_prepared_response.py")
    for node in programme["calculation_graph"]:
        if node["id"] in {"wall-trace", "observer-response"}:
            node.update(status="evaluated_for_fixed_finite_IV_preparation",
                evaluation_scope=scenario["id"],
                missing=["remote physical source generation", "later clock-domain evolution",
                         "a second compatible preparation for observer inference"],
                decision="Reuse the recorded finite response; admit further work only to compare compatible sources or observables.")
            node["owner_refs"] = sorted(set(node["owner_refs"] + [OWNER]))
            node["evidence_refs"] = sorted(set(node["evidence_refs"] + [EVIDENCE]))
    programme["completed_calculation"] = {
        "admission_ref": "docs/prior-art-and-open-claim.md#prepared-scalar-response-admission-12-october-2026", "scope": scenario["id"],
        "future_calculations_require_admission": True,
        "physical_source_generation_solved": False, "observer_equivalence_established": False}
    record["programme"] = programme
    record["scenarios"].append(scenario)
    record["assessments"][scenario["id"]] = {
        "record": "prepared-response.json", "transfer_operator_certified": True,
        "wall_waveform_evaluated": True, "observer_waveform_evaluated": True,
        "right_parent_access_changed": False, "theory_or_preparation_retuned": False,
        "positive_tail_in_finite_window": True, "permanent_memory_established": False,
        "global_preparation_or_quantum_state_solved": False,
        "numerical_refinement_is_rigorous_enclosure": False,
        "v7_analytic_operator_norm_bound": prepared["operator"]["analytic_norm_bound"]}
    record["admission"].update(category="fixed exterior-data-to-observer response and finite post-input tail",
        owner=OWNER, missing_arrow="a second physically compatible response or independently specified observer discriminator")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    record["predecessor_manifest"] = "../nsc-cosmological-atlas-v7/manifest.json"
    return record, programme, prepared, response


def source_bindings():
    bindings = previous.source_bindings()
    paths = {Path(__file__), Path(response_model.__file__), REPO / OWNER,
        REPO / "tests/test_nsc_atlas_prepared_response.py",
        REPO / "tests/test_nsc_cosmological_atlas_v8.py"}
    # The canonical paper consumes this evidence and must not bind its producer.
    bindings.update({path.relative_to(REPO).as_posix(): first._sha(path) for path in paths})
    return dict(sorted(bindings.items()))


def waveform_figure(response):
    """Three fixed-input views; the tail panel is a wall trace, not memory."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    samples = response["samples"]
    s = np.asarray(samples["source_coordinate"])
    T = np.asarray(samples["arrival_t_C_over_ell"])
    tail = s >= .75
    controls = response["analytic_tail"]["independent_controls"]
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans",
                                "font.size": 13, "axes.titlesize": 15}):
        fig, axes = plt.subplots(3, 1, figsize=(8.4, 10.2), layout="constrained")
        ax = axes[0]
        ax.plot(s, samples["source_chi_over_epsilon"], color="#b26825", lw=2.2,
                label=r"Independent input $g/\varepsilon$")
        ax.plot(s, samples["wall_chi_over_epsilon"], color="#246485", lw=2.0,
                label=r"Computed wall $\chi_W/\varepsilon$")
        ax.axvspan(.75, 1., color="#edf3e8", zorder=0)
        ax.set(xlabel=r"Exterior null label $s=(V-V_a)/(V_b-V_a)$",
               ylabel="Reduced field / input amplitude", title="Independent data determine the wall trace")
        ax.legend(fontsize=11.5, loc="upper right")
        ax = axes[1]
        ax.plot(T, samples["ell_phi_C_over_epsilon"], color="#246485", lw=2.2,
                label="Full exterior propagation")
        ax.plot(T, samples["zero_kernel_phi_C_over_epsilon"], color="#b26825", lw=1.6,
                linestyle="--", label="Kernel omitted: diagnostic only")
        ax.axvspan(float(np.interp(.75, s, T)), T[-1], color="#edf3e8", zorder=0)
        ax.set(xlabel=r"Central proper time $t_C/\ell$", ylabel=r"Observed field $\ell\phi_C/\varepsilon$",
               title="Interior response on its admitted clock interval")
        ax.legend(fontsize=11.5, loc="upper right")
        ax = axes[2]
        scale = 1e6
        J = np.asarray(samples["wall_dchi_dV_over_epsilon"])
        ax.plot(s[tail], scale*J[tail], color="#246485", lw=2.2, label="Computed post-input trace")
        for index, row in enumerate(controls):
            x, lo, hi = row["source_coordinate"], scale*row["lower_from_B"], scale*row["upper_from_B"]
            ax.plot([x, x], [lo, hi], color="#8b5425", lw=2.4,
                    label=r"Independent $(1-q_*)B\leq J\leq B$" if index == 0 else None)
            ax.plot([x], [scale*row["J"]], "o", color="#246485", ms=5)
            ax.plot([x-.003, x+.003], [lo, lo], color="#8b5425", lw=1.5)
            ax.plot([x-.003, x+.003], [hi, hi], color="#8b5425", lw=1.5)
        floor = scale*response["analytic_tail"]["uniform_lower_J_over_epsilon"]
        ax.axhline(floor, color="#607b3d", ls=":", lw=2, label="Analytic uniform positive lower bound")
        ax.axhline(0, color="#777777", lw=.7)
        ax.set(xlim=(.745, 1.008), xlabel=r"Post-input null label $s$ (input ends at $3/4$)",
               ylabel=r"$10^6 J/\varepsilon$,  $J=\partial_V\chi_W$",
               title="Magnified finite tail, with independent bounds")
        ax.legend(fontsize=10.5, loc="best")
        for ax in axes:
            ax.grid(alpha=.18)
            ax.ticklabel_format(style="plain", axis="y", useOffset=False)
        fig.suptitle("Exterior preparation → wall → interior observer", fontsize=17)
        fig.supxlabel("The same fixed geometry, scalar law and v7 data are used throughout.\nA finite field tail does not establish permanent memory or a detector record.", fontsize=11)
        outputs = {}
        for fmt in ("svg", "png"):
            stream = BytesIO()
            metadata = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
            fig.savefig(stream, format=fmt, dpi=180, metadata=metadata, facecolor="white")
            outputs[f"prepared-response.{fmt}"] = stream.getvalue()
        plt.close(fig)
    return outputs


def narrative(prepared, response):
    q = prepared["operator"]["analytic_norm_bound"]
    tail = response["analytic_tail"]
    checks, metrics = response["checks"], response["metrics"]
    rows = "\n".join(f"| {row['source_coordinate']:.4g} | {row['lower_from_B']:.10g} | {row['J']:.10g} | {row['upper_from_B']:.10g} | {row['B_quadrature_difference']:.3g} |"
                     for row in tail["independent_controls"])
    return f"""# An independently prepared signal reaches the interior observer

![Source, wall, interior response and finite tail](prepared-response.svg)

The v7 exterior datum now determines a computed wall characteristic and a
central classical waveform. Geometry, coupling, input profile, zero initial
cap at R/ell=0.45 and the wall interval [0.45,0.55] remain fixed. The designated
right parent remains causally excluded. Source generation and quantum-state
selection remain open. All ten earlier scenarios are preserved.

The existing Volterra norm bound is q*={q:.12f}<1. Product integration and
the borrowed DOP853 method evaluate that operator on the retained domain;
the exterior kernel is never evaluated inside the excised cap. The observer
map is the existing phi_C[T(tau)]=I(tau)/Tprime(tau).

The finite response has sampled extrema {metrics['negative_sampled_peak']:.9g}
and {metrics['positive_sampled_peak']:.9g} in ell phi_C/epsilon. Omitting the
kernel changes the waveform by at most {metrics['scattering_max_absolute_change']:.9g}
on the recorded samples. That omission is a diagnostic, not a second physical
spacetime or an observer-equivalence result.

After the input ends, the fixed kernel and data imply a positive wall derivative
J on the remaining window. The analytic floor is
J/epsilon >= {tail['uniform_lower_J_over_epsilon']:.12g}. Independent integration
gives the more informative brackets below; its 32/64-order differences are
numerical diagnostics, not rigorous enclosures.

| Normalized null label s | (1-q*) B / epsilon | Computed J / epsilon | B / epsilon | B 32/64 difference |
|---:|---:|---:|---:|---:|
{rows}

Field-grid maximum waveform differences are
{checks['field_grid_refinement_max_differences']}; the 129/257 geometry-knot
comparison differs by {checks['geometry_129_257_max_phi_difference']:.4g}.
The independent observer-area identity differs by
{checks['observer_area_identity_difference']:.4g}. Refinement is not a rigorous
continuum error bound. Complete tolerances and finite budgets are recorded in
[prepared-response.json](prepared-response.json).

The positive finite-window field tail does not establish permanent memory,
a detector record, energy density or late-time asymptotics. The upper observer
time is a domain exit: no zero tail or continuation to t_C=0 is inferred.
The next question is whether another independently specified, compatible
preparation can be distinguished by this observer; that comparison is open.

[atlas.json](atlas.json) and [programme.json](programme.json) scope the completed
connection and keep the other routes open. [preparation-access.json](preparation-access.json)
retains the v7 preparation certificate. The scientific owner is
[prepared response](../../../{OWNER}); the living
[paper](../../../paper/nested-space-cosmology.md) remains the manuscript.
"""


def payloads():
    atlas, programme, prepared, response = build_record()
    outputs = {"atlas.json": first.atlas.canonical_bytes(atlas),
        "programme.json": first.atlas.canonical_bytes(programme),
        "preparation-access.json": first.atlas.canonical_bytes(prepared),
        "prepared-response.json": first.atlas.canonical_bytes(response),
        "narrative.md": narrative(prepared, response).encode("utf-8")}
    outputs.update(waveform_figure(response))
    return outputs


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only atlas: use --check or a fresh versioned directory")
    before = source_bindings()
    outputs = payloads()
    if source_bindings() != before:
        raise RuntimeError("source changed during atlas generation")
    manifest = {"schema": SCHEMA, "source_bindings": before,
        "input_sha256": first.atlas.hashlib.sha256(outputs["atlas.json"]).hexdigest(),
        "outputs": {name: {"sha256": first.atlas.hashlib.sha256(data).hexdigest(), "bytes": len(data)} for name, data in outputs.items()},
        "runtime": {"python": platform.python_version(), **{name: version(name) for name in ("numpy", "scipy", "sympy", "matplotlib")}},
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12},
        "scope": "fixed finite IV source-to-observer response and positive field tail; no permanent memory"}
    directory.mkdir(parents=True, exist_ok=False)
    for name, data in outputs.items():
        with (directory / name).open("xb") as handle:
            handle.write(data)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(first.atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "permanent_memory_established": False}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != SCHEMA or manifest["source_bindings"] != source_bindings():
        raise ValueError("atlas source binding mismatch")
    expected = payloads()
    if set(manifest["outputs"]) != set(expected):
        raise ValueError("atlas output inventory mismatch")
    for name, replay in expected.items():
        path, binding = directory / name, manifest["outputs"][name]
        if first._sha(path) != binding["sha256"] or path.stat().st_size != binding["bytes"]:
            raise ValueError("artifact binding mismatch: " + name)
        equal = first._close(json.loads(path.read_bytes()), json.loads(replay)) if name.endswith(".json") else path.read_bytes() == replay
        if not equal:
            raise ValueError("atlas replay mismatch: " + name)
    if manifest["input_sha256"] != first._sha(directory / "atlas.json"):
        raise ValueError("canonical input binding mismatch")
    return {"checked": directory.name, "outputs": len(expected), "sources": len(manifest["source_bindings"])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
