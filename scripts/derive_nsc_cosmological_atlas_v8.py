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



# Fixed-preparation control inherited from the local source-response producer.
# It used at most four solves / 300 CPU seconds, 7601 columns, alpha=.001,
# delta-alpha=.0001, 60/120 time steps and a 1e-11 temporal indicator.
from hashlib import sha256
from types import SimpleNamespace


def _prepared_chart_coordinates(rho):
    from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
    return chart_coordinates(rho)


def _prepared_coefficients_from_records(*records):
    from derive_nsc_evolved_incoming_constraints import coefficients_from_records
    return coefficients_from_records(*records)


_prepared_ROOT = Path(__file__).resolve().parents[1]

_prepared_OUTPUT = 'results/development/nsc-local-prepared-response.json'

_prepared_ARTIFACT = 'results/development/artifacts/nsc-local-prepared-response.npz'

_prepared_OWNED = ('scripts/derive_nsc_cosmological_atlas_v8.py', 'docs/nsc-atlas-prepared-response.md')

_prepared_INPUTS = ('results/development/nsc-local-reference-accuracy.json', 'results/development/nsc-incoming-source-update-v5.json', 'results/development/nsc-incoming-surface-coefficients.json', 'results/development/nsc-incoming-surface-regular-branch.json', 'results/development/nsc-incoming-fixed-transfer.json', 'src/recursive_horizons/nsc_exact_phase_prepared_state.py', 'src/recursive_horizons/nsc_evolved_incoming_state.py', 'src/recursive_horizons/nsc_evolved_incoming_constraints.py', 'src/recursive_horizons/nsc_local_incoming_family.py', 'src/recursive_horizons/nsc_compatible_history_geometry.py', 'src/recursive_horizons/nsc_transmitting_history_jets.py', 'src/recursive_horizons/nsc_transmitting_history_modes.py', 'scripts/derive_nsc_evolved_incoming_constraints.py')

_prepared_ALPHA = 0.001

_prepared_STEP = 0.0001

_prepared_CAP = 300.0

_prepared_CASES = (('center', _prepared_ALPHA, 60, True), ('plus', _prepared_ALPHA + _prepared_STEP, 60, False), ('minus', _prepared_ALPHA - _prepared_STEP, 60, False), ('time_fine', _prepared_ALPHA, 120, True))

def _prepared_digest(path):
    return sha256((_prepared_ROOT / path).read_bytes()).hexdigest()

def _prepared_signatures():
    return {p: _prepared_digest(p) for p in (*_prepared_OWNED, *_prepared_INPUTS)}

def _prepared_norm(x):
    return float(np.max(abs(x), initial=0.0))

def _prepared_family(alpha):
    from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
    from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
    center = _prepared_chart_coordinates(1.0)[1] + 0.15
    w = LocalAxialFunction((0.0, 1.0), center)
    zero = LocalAxialFunction((0.0,), center)
    return CompatibleIncomingMetric((alpha,), (CompatibleRadiusDirection(w, zero, 0.007, 0.03),))

def _prepared_load_inputs():
    records = [json.loads((_prepared_ROOT / p).read_text()) for p in _prepared_INPUTS[:5]]
    p = records[0]['payload']
    if _prepared_digest(p['path']) != p['sha256']:
        raise ValueError('refined initial column artifact changed')
    with np.load(_prepared_ROOT / p['path'], allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    channel = records[4]['actual_source_bindings']['channel']
    if channel['compact_mass'] != float(arrays['source/mass'].item()) or channel['angular_eigenvalue'] != float(arrays['source/angular'].item()):
        raise ValueError('source/channel operator changed')
    return (records, arrays, channel)

def _prepared_evaluate(arrays, name, alpha, coeff, baseline, parameters):
    from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter, surface_geometry_response, compatible_history_slots

    def get(key):
        return arrays[name + '/' + key]
    times = get('times')
    z = times + _prepared_chart_coordinates(1.0)[1]
    F, Fz, dF, dFz = (get(k) for k in ('F', 'Fz', 'dF', 'dFz'))
    weights = arrays['source/column_weights']
    C = arrays['source/covariance']
    E = arrays['source/energies_repeated']
    current = source_column_matter(F * weights, Fz * weights, C, dF * weights, dFz * weights, **parameters)
    ref = arrays['reference_F0'][None] * np.exp(-1j * times[:, None, None] * E)
    refz = -1j * E * ref
    zero = np.zeros((0, *ref.shape), complex)
    reference = source_column_matter(ref * weights, refz * weights, C, zero, zero, **parameters)
    slots, ds = compatible_history_slots(_prepared_family(alpha), z, len(dF))
    geom = surface_geometry_response(slots, ds, coeff)
    value = baseline + geom['action_gradient_change'] + current['action_gradient'] - reference['action_gradient']
    tangent = geom['action_gradient_tangent'] + current['action_gradient_tangent']
    return {'residual': value, 'tangent': tangent, 'matter': current['action_gradient'], 'matter_tangent': current['action_gradient_tangent'], 'matter_correction': current['action_gradient'] - reference['action_gradient']}

def _prepared_analyze(arrays, meta):
    rows = {}
    for name, alpha, steps, _ in _prepared_CASES:
        if name not in meta['completed']:
            continue
        got = _prepared_evaluate(arrays, name, alpha, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        mask = np.arange(steps + 1) >= 2 * steps // 3
        rows[name] = {'sampled_local_constraint_approximant_max': np.max(abs(got['residual'][mask]), axis=0).tolist(), 'sampled_local_matter_correction_max': np.max(abs(got['matter_correction'][mask]), axis=0).tolist(), 'phase_residual': meta['cases'][name]['phase_residual'], 'flux_residual': meta['cases'][name]['flux_residual']}
    checks = {}
    if all((n in rows for n in ('center', 'plus', 'minus'))):
        center = _prepared_evaluate(arrays, 'center', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        plus = _prepared_evaluate(arrays, 'plus', _prepared_ALPHA + _prepared_STEP, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        minus = _prepared_evaluate(arrays, 'minus', _prepared_ALPHA - _prepared_STEP, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        fd = (plus['residual'] - minus['residual']) / (2 * _prepared_STEP)
        checks['complete_local_constraint_derivative'] = _prepared_norm((fd - center['tangent'][0])[40:])
        checks['local_matter_derivative'] = _prepared_norm(((plus['matter'] - minus['matter']) / (2 * _prepared_STEP) - center['matter_tangent'][0])[40:])
        checks['dropped_deltaC_effect'] = _prepared_norm(center['matter_tangent'][0, 40:])
        checks['same_fixed_preparation'] = len({meta['cases'][n]['preparation'] for n in ('center', 'plus', 'minus')}) == 1
        checks['derivative_passed'] = checks['complete_local_constraint_derivative'] <= 3e-08 and checks['same_fixed_preparation']
    if 'time_fine' in rows and 'center' in rows:
        coarse = _prepared_evaluate(arrays, 'center', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        fine = _prepared_evaluate(arrays, 'time_fine', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        checks['local_time_change'] = np.max(abs(fine['matter'][::2][40:] - coarse['matter'][40:]), axis=0).tolist()
        checks['local_tangent_time_change'] = np.max(abs(fine['matter_tangent'][0, ::2][40:] - coarse['matter_tangent'][0, 40:]), axis=0).tolist()
        checks['time_indicator_target'] = 1e-11
        checks['time_indicator_passed'] = max(checks['local_time_change']) <= 1e-11
    return {'per_case': rows, 'checks': checks}

def _prepared_record(arrays, meta, payload):
    measurements = _prepared_analyze(arrays, meta)
    return {'schema': 'NSC-LOCAL-PREPARED-RESPONSE-v1', 'accountable_author': 'Douglas Ek', 'status': 'OPEN: source-fixed nonzero-history control; full-source local gate uncertified', 'control': {'amplitude': _prepared_ALPHA, 'finite_difference_step': _prepared_STEP, 'runs_max': 4, 'profile': 'w=alpha*Chebyshev_T1*flat axial cutoff, U=0; normal .007/.03', 'family': '14_1', 'points': 7601, 'local_interval': 'S(1)+[.12,.18]', 'derivative_tolerance': 3e-08, 'temporal_indicator_tolerance': 1e-11}, 'measurements': measurements, 'runtime': meta['runtime'], 'scope': {'physical_local_gate': 'OPEN', 'full_source_error': None, 'field_and_derivative_error': None, 'between_node_error': None, 'all_channel_coverage': False, 'source_changed': False, 'stress_drift_subtracted': False, 'frozen_C0_imposed': False, 'metric_timestep': False}, 'source_hashes': {p: meta['signature'][p] for p in _prepared_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _prepared_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'call prepared_response.check() from scripts/derive_nsc_cosmological_atlas_v8.py'}

def _prepared_check():
    saved = json.loads((_prepared_ROOT / _prepared_OUTPUT).read_text())
    _prepared_load_inputs()
    if _prepared_digest(_prepared_ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('control payload changed')
    with np.load(_prepared_ROOT / _prepared_ARTIFACT, allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != _prepared_signatures():
        raise ValueError('control source/input changed')
    if _prepared_record(arrays, meta, saved['payload']) != saved:
        raise ValueError('array-only replay differs')
    return saved

prepared_response = SimpleNamespace(
    ALPHA=_prepared_ALPHA, STEP=_prepared_STEP, CAP=_prepared_CAP, CASES=_prepared_CASES,
    INPUTS=_prepared_INPUTS, OWNED=_prepared_OWNED, OUTPUT=_prepared_OUTPUT,
    ARTIFACT=_prepared_ARTIFACT, family=_prepared_family, load_inputs=_prepared_load_inputs,
    evaluate=_prepared_evaluate, analyze=_prepared_analyze, check=_prepared_check,
    digest=_prepared_digest, coefficients_from_records=_prepared_coefficients_from_records,
    chart_coordinates=_prepared_chart_coordinates)

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
