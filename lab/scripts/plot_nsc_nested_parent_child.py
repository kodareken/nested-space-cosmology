#!/usr/bin/env python3
"""Replay one scientific figure from final nested-pair records, never evolve.

The primary witness remains nf256; nf512 confirms the named parent-to-child
comparisons. Physical-window probabilities are not fixed-mode occupations.
All control effects subtract their own initial offsets at common coordinate T.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".build/nsc-nested-parent-child-figure/mpl"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

WITNESS = LAB / "results/development/nsc-nested-parent-child-v1.json"
CONFIRMATION = LAB / "results/development/nsc-nested-parent-child-confirmation-v2.json"
OUTPUT = LAB / "results/development/nsc-nested-parent-child-figure.png"
MANIFEST = OUTPUT.with_suffix(".json")
PRIMARY = "nf256_baseline_dt0.0005"
FINE = "nf512_baseline_dt0.0005"
EFFECTS = {
    "parent_to_child_probability": ("parent_only", "child_probability"),
    "child_to_parent_annulus_probability": ("child_only", "parent_annulus_probability"),
    "parent_to_child_radius": ("parent_only", "child_r_proper_mean"),
    "child_to_parent_annulus_radius": ("child_only", "parent_annulus_r_proper_mean"),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_record(path, schema):
    record = json.loads(path.read_text())
    if record.get("schema") != schema or digest(path.with_suffix(".npz")) != record["payload_sha256"]:
        raise ValueError(f"record/payload binding failed: {path.name}")
    return record


def times(case):
    result = np.array([row["time"] for row in case["rows"]], dtype=float)
    if len(result) < 2 or result[0] != 0 or abs(result[-1] - .3) > 1e-12 or np.any(np.diff(result) <= 0):
        raise ValueError("figure requires complete ordered witness frames on [0,.3]")
    return result


def values(case, observable):
    if observable == "child_probability":
        return np.array([row["windows"]["child"]["probability"] for row in case["rows"]])
    if observable == "parent_annulus_probability":
        return np.array([row["windows"]["left_annulus"]["probability"] +
                         row["windows"]["right_annulus"]["probability"] for row in case["rows"]])
    return np.array([row["metrics"][observable] for row in case["rows"]])


def effect(cases, nf, control, observable, cap=.0005):
    baseline, changed = cases[f"nf{nf}_baseline_dt{cap:g}"], cases[f"nf{nf}_{control}_dt{cap:g}"]
    if not np.allclose(times(baseline), times(changed), rtol=0, atol=1e-12):
        raise ValueError("control frames use different coordinate clocks")
    first, second = values(baseline, observable), values(changed, observable)
    return (second - second[0]) - (first - first[0])


def initial_density(path, case_name, weights):
    """AP interpolation of SAVED T0 columns, with correct half-density weight."""
    with np.load(path, allow_pickle=False) as payload:
        phi0 = np.array(payload[case_name + "_phi0"][0])
        phi1 = np.array(payload[case_name + "_phi1"][0])
    nf = phi0.shape[0]
    if phi0.shape != (nf, 6) or phi1.shape != (nf, 6):
        raise ValueError("initial figure source must be the saved rank-six columns")
    length, nq = 8., 4 * nf
    coarse_x, x = np.arange(nf) * length / nf, np.arange(nq) * length / nq
    modes = np.fft.fftfreq(nf) * nf + .5
    phase = np.exp(-1j * np.pi * coarse_x / length)[:, None]
    synthesis = np.exp(2j * np.pi * x[:, None] * modes[None, :] / length)
    interpolate = lambda columns: synthesis @ (np.fft.fft(columns * phase, axis=0) / nf) / np.sqrt(length / nf)
    fine0, fine1 = interpolate(phi0), interpolate(phi1)
    density = np.sum((np.abs(fine0) ** 2 + np.abs(fine1) ** 2) * np.array(weights)[None, :], axis=1)
    mass = float(length / nq * np.sum(density))
    if abs(mass - np.sum(weights)) > 1e-10:
        raise ValueError("initial probability-density normalization failed")
    return {"x": x.tolist(), "density": density.tolist(), "integral": mass,
            "fermion_count": nf, "quadrature_count": nq,
            "angular_multiplicity_applied": False}


def dataset():
    witness = load_record(WITNESS, "NSC-NESTED-PARENT-CHILD-v1")
    confirmation = load_record(CONFIRMATION, "NSC-NESTED-PARENT-CHILD-CONFIRMATION-v2")
    if confirmation.get("stage") != "complete_confirmation" or confirmation.get("status") != "MEASURED_CONFIRMED_NESTED_PAIR":
        raise ValueError("only the FINAL confirmed v2 record is a figure input")
    cases = witness["results"]
    fine_cases = confirmation["results"]
    coordinate_times = times(cases[PRIMARY])
    if not np.allclose(coordinate_times, times(fine_cases[FINE]), rtol=0, atol=1e-12):
        raise ValueError("witness and confirmation frame clocks do not coincide")
    curves, indicators = {}, {}
    for name, (control, observable) in EFFECTS.items():
        primary = effect(cases, 256, control, observable)
        coarse = effect(cases, 128, control, observable)
        temporal = effect(cases, 256, control, observable, .001)
        curves[name] = primary.tolist()
        peak = float(np.max(np.abs(primary)))
        time_gap = float(np.max(np.abs(temporal - primary)))
        comparison = effect(fine_cases, 512, control, observable) if control == "parent_only" else coarse
        gap = float(np.max(np.abs(comparison - primary)))
        indicators[name] = {"primary_absolute_peak": peak, "time_indicator": time_gap,
                            "time_fraction_of_primary_effect": time_gap / peak,
                            "space_indicator": gap, "space_fraction_of_primary_effect": gap / peak,
                            "space_comparison_nf": [256, 512] if control == "parent_only" else [128, 256],
                            "observable": observable, "control": control}
        if control == "parent_only":
            curves[name + "_nf512_confirmation"] = comparison.tolist()
    baseline = cases[PRIMARY]
    fine_baseline = fine_cases[FINE]
    clocks = np.array(baseline["summary"]["normal_clocks"])
    fine_clocks = np.array(fine_baseline["summary"]["normal_clocks"])
    if clocks.shape != (len(coordinate_times), 3) or fine_clocks.shape != clocks.shape:
        raise ValueError("saved original normal-clock curves changed shape")
    density = initial_density(WITNESS.with_suffix(".npz"), PRIMARY, baseline["source"]["occupations"])
    fine_density = initial_density(CONFIRMATION.with_suffix(".npz"), FINE, fine_baseline["source"]["occupations"])
    data = {"times": coordinate_times.tolist(), "initial_nf256": density,
            "initial_nf512": fine_density, "effects": curves,
            "baseline_nf256_proper_length_ratio": values(baseline, "proper_length_ratio").tolist(),
            "baseline_nf512_proper_length_ratio": values(fine_baseline, "proper_length_ratio").tolist(),
            "baseline_nf256_normal_clocks": clocks.tolist(),
            "baseline_nf512_normal_clocks": fine_clocks.tolist()}
    bindings = {str(path.relative_to(ROOT)): digest(path) for path in
                (Path(__file__), WITNESS, WITNESS.with_suffix(".npz"),
                 CONFIRMATION, CONFIRMATION.with_suffix(".npz"))}
    return data, indicators, bindings, witness, confirmation


def render(data):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.titlesize": 11, "axes.labelsize": 10,
                         "legend.fontsize": 8.5, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 160})
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    blue, red, dark = "#2166ac", "#b2182b", "#222222"
    top = axes[0, 0]
    top.axvspan(0, 4, color="#dbe9f5", label=r"Parent $I_P=[0,4]$")
    top.axvspan(1, 3, color="#f5dfbc", label=r"Child $I_C=[1,3]$")
    for key, style, label in (("initial_nf256", "-", "Saved initial source, nf256"),
                              ("initial_nf512", "--", "Saved initial source, nf512")):
        row = data[key]
        top.plot(row["x"], row["density"], style, color=dark, lw=1.25, label=label)
    top.set(xlim=(0, 8), ylim=(0, None), title="(a) Separated initial field on period 8",
            xlabel="Radial coordinate x (dimensionless)",
            ylabel="Canonical probability density (dimensionless)")
    top.legend(loc="upper right", framealpha=.95)
    time = np.array(data["times"])
    effects = data["effects"]
    for axis, stem, title, ylabel in (
            (axes[0, 1], "probability", "(b) Two-way physical state response", "Probability effect (dimensionless)"),
            (axes[1, 0], "radius", "(c) Two-way physical geometry response", "Proper-mean radius effect (reference units)")):
        parent = "parent_to_child_" + stem
        child = "child_to_parent_annulus_" + stem
        axis.plot(time, effects[parent], color=blue, lw=1.7, label="Parent → child, nf256 witness")
        axis.plot(time, effects[parent + "_nf512_confirmation"], color=blue, lw=1.25, ls="--",
                  label="Parent → child, nf512 confirmation")
        axis.plot(time, effects[child], color=red, lw=1.7, label="Child → parent annulus, nf256")
        axis.axhline(0, color=".65", lw=.65)
        axis.set(xlim=(0, time[-1]), title=title, xlabel="Common coordinate T (dimensionless)", ylabel=ylabel)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
        axis.legend(loc="best", framealpha=.95)
    bottom = axes[1, 1]
    bottom.plot(time, data["baseline_nf256_proper_length_ratio"], color=dark, lw=1.7,
                label=r"$\ell_P/\ell_C$, nf256")
    bottom.plot(time, data["baseline_nf512_proper_length_ratio"], color=dark, lw=1.2, ls="--",
                label=r"$\ell_P/\ell_C$, nf512")
    bottom.set(xlim=(0, time[-1]), title="(d) Measured scale ratio and original normal clocks",
               xlabel="Common coordinate T (dimensionless)", ylabel="Proper-length ratio (dimensionless)")
    clock_axis = bottom.twinx()
    clock_axis.spines["right"].set_visible(True)
    clocks = np.array(data["baseline_nf256_normal_clocks"])
    clock_axis.plot(time, clocks[:, 0], color="#4d9221", lw=1.35, label=r"$\tau_{x=1}$, nf256")
    clock_axis.plot(time, clocks[:, 1], color="#762a83", lw=1.35, ls=":", label=r"$\tau_{x=2}$, nf256")
    clock_axis.set_ylabel("Normal proper clock τ (reference units)")
    handles, labels = bottom.get_legend_handles_labels()
    second_handles, second_labels = clock_axis.get_legend_handles_labels()
    bottom.legend(handles + second_handles, labels + second_labels, loc="upper left", framealpha=.95)
    for axis in axes.flat:
        axis.grid(alpha=.2, lw=.5)
    fig.suptitle("One finite nested parent–child pair: physical domains and two-way response", fontsize=14)
    fig.supxlabel("Control effects: (control − its T0 value) − (baseline − its T0 value).\n"
                  "Physical window probabilities are distinct from fixed-mode occupations; clocks are displayed, not matched.", fontsize=9)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=160, metadata={"Title": "Finite nested parent-child pair",
                "Description": "Saved nf256 witness and nf512 confirmation; no re-evolution",
                "Software": "Matplotlib " + matplotlib.__version__})
    plt.close(fig)
    return buffer.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="replay curves and PNG bytes without writing")
    options = parser.parse_args()
    data, indicators, bindings, witness, confirmation = dataset()
    image = render(data)
    image_hash = hashlib.sha256(image).hexdigest()
    result = {"schema": "NSC-NESTED-PARENT-CHILD-FIGURE-v1", "bindings": bindings,
              "figure_sha256": image_hash, "figure_bytes": len(image),
              "matplotlib_version": matplotlib.__version__, "numpy_version": np.__version__,
              "primary_case": PRIMARY, "confirmation_case": FINE,
              "witness_record_status": witness["status"], "confirmation_record_status": confirmation["status"],
              "witness_physical_cases_complete": True, "trajectory_re_evolved": False,
              "coordinate_clock": "same global T, all plotted times [0,.3] dimensionless",
              "proper_clock_locations": [1., 2.], "proper_clock_matching_claimed": False,
              "parent_domain": [0., 4.], "child_domain": [1., 3.],
              "parent_annulus": [[0., 1.], [3., 4.]], "fixed_modal_observer_plotted": False,
              "initial_probability_definition": "sum_a c_a (|phi0_a|^2+|phi1_a|^2)/dx; AP interpolation; no angular M",
              "effect_definition": "(control-control[0])-(baseline-baseline[0])",
              "radius_definition": "spatial proper-length-weighted mean of full reconstructed areal radius",
              "indicator_scope": "named physical effects; numerical comparisons, not a continuum certificate",
              "effect_indicators_recomputed": indicators, "curve_count": 12, "data": data}
    if options.check:
        stored = json.loads(MANIFEST.read_text())
        if stored != result or digest(OUTPUT) != image_hash:
            raise ValueError("bound curve/figure replay differs")
        print("PASS: 4 panels, 12 bound curves, 4 physical-effect comparisons replay exactly")
    else:
        if OUTPUT.exists() or MANIFEST.exists():
            raise FileExistsError("figure outputs already exist; use --check for read-only replay")
        OUTPUT.write_bytes(image)
        MANIFEST.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print("SAVED", OUTPUT, "bytes", len(image), "curves", result["curve_count"])


if __name__ == "__main__":
    main()
