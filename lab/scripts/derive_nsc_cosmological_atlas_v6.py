#!/usr/bin/env python3
"""Record the conditional conformal-field response in the interior atlas."""
from __future__ import annotations

import argparse
from copy import deepcopy
from importlib.metadata import version
from io import BytesIO
import json
from pathlib import Path
import platform

import numpy as np

import derive_nsc_cosmological_atlas_v5 as previous
from recursive_horizons import nsc_atlas_field_transfer as field

first = previous.first
LAB, REPO = previous.LAB, previous.REPO
SCHEMA = "nsc-cosmological-atlas-v6"
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v6"


def build_record():
    record, _, _ = previous.build_record()
    response = field.field_transfer_record()
    record["sources"].update({
        "scalar_matching": {"path": "lab/docs/nsc-atlas-field-matching.md", "use": "borrowed distributional scalar matching and conformal reduction"},
        "scalar_state_transfer": {"path": "lab/docs/nsc-atlas-field-transfer.md", "use": "specified classical incoming state and exact conditional interior transfer"},
    })
    scenario = deepcopy(record["scenarios"][-1])
    scenario.update(id="early_conformal_field_response", title="An early scalar imprint at the interior centre",
        solution_status="conditional field response for prescribed incoming characteristic data",
        field_record="field-transfer.json")
    scenario.pop("access_record", None)
    scenario["theory"].update(test_field="classical massless scalar, angular mode zero", curvature_coupling=1/6,
        added_surface_interaction=False, backreaction=False)
    scenario["assumptions"] += ["zero earlier interior classical data", "compact exterior-side incoming characteristic prescribed in wall proper time",
        "global exterior initial state and source preparation not solved", "fixed incoming history, including zero outside its pulse"]
    for region in scenario["regions"]:
        region["state"] = {"C": "zero earlier classical interior data; response to the declared incident trace",
            "W": "one prescribed incoming classical characteristic; outgoing component follows from the interior",
            "P": "global exterior field state not supplied"}[region["id"]]
    scenario["links"] = [first.atlas._link("scalar-W-C", "W", "C", "causal", ["scalar_matching", "scalar_state_transfer"],
        "phi_C(T(tau_W)) = I(tau_W) / Tprime(tau_W); finite waveform in field-transfer.json",
        status="conditional_field_transfer", missing=["global exterior preparation and incoming-state selection", "backreaction and detector dynamics"],
        delay=None, interval=None, clock=None)]
    record["scenarios"].append(scenario)
    record["assessments"]["early_conformal_field_response"] = {
        "record": "field-transfer.json", "incident_state_is_input": True,
        "remote_parent_source_solved": False, "quantum_vacuum_selected": False,
        "field_signal_computed": True, "energy_density_or_dark_component_computed": False}
    record["admission"].update(category="classical field-state and wall-to-centre transfer",
        missing_arrow="derive admissible exterior incident data, or compare a field with a supported bulk tail")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    record["predecessor_manifest"] = "../nsc-cosmological-atlas-v5/manifest.json"
    return record, response


def source_bindings():
    bindings = previous.source_bindings()
    paths = [Path(__file__), Path(field.__file__), LAB / "tests/test_nsc_atlas_field_transfer.py",
        LAB / "tests/test_nsc_cosmological_atlas_v6.py", LAB / "docs/nsc-atlas-field-matching.md",
        LAB / "docs/nsc-atlas-field-transfer.md", LAB / "docs/nsc-cosmological-atlas-v6.md",
        REPO / "paper/atlas-working-draft-v6.md"]
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def render_map():
    text, panel = first._text, first._panel
    out = first._base(1440, 1360, "What an early signal leaves inside",
        "Atlas v6 | Fixed geometry + explicit classical state + borrowed field transmission")
    out += [panel(45, 130, 1350, 118, "#153329", "#70ac94"),
        text(75, 174, "A causal path now carries a calculated field waveform.", 27, "#a5e6c9", "600"),
        text(75, 216, "The signal is an input at the early wall. Its remote preparation remains an open connection.", 20)]
    for x, title, lines in ((45, "1  PREPARED INPUT", ["Initially zero scalar field in cap", "Compact incoming scalar pulse", "Profile fixed in wall proper time"]),
        (505, "2  MATCHED WALL", ["Classical conformal scalar", "Curvature jump included", "One incoming characteristic"]),
        (965, "3  CENTRAL OBSERVER", ["Exact arrival and amplitude map", "A finite, distorted transient", "Zero after pulse through t_C=0"])):
        out += [panel(x, 300, 430, 210, "#152c40", "#567f96"), text(x+20, 343, title, 22, "#a1d9f3", "600")]
        out += [text(x+20, 391+39*i, line, 18) for i, line in enumerate(lines)]
    out += [text(480, 415, ">", 26, "#92dfc4"), text(940, 415, ">", 26, "#92dfc4"),
        panel(45, 555, 1350, 156, "#2a2637", "#9c82b7"),
        text(75, 602, "THE TRANSFER LAW", 23, "#d0b5f0", "600"),
        text(75, 649, "central field at arrival = incoming characteristic / proper-time stretch", 26),
        text(75, 686, "The clock map and amplitude factor come from the same fixed wall trajectory.", 18),
        text(55, 771, "THE COMMON THREAD THROUGH THE PAPERS", 24, "#a1d9f3", "600")]
    for x, title, line in ((45, "Foundation", "Boundary, spectrum and inherited state"),
        (505, "Finite companion", "Measured transfer, feedback and memory"),
        (965, "Cosmological atlas", "Geometric access and observer response")):
        out += [panel(x, 804, 430, 110), text(x+20, 844, title, 23, "#c9b2ec", "600"), text(x+20, 882, line, 17)]
    out += [text(55, 957, "Shared question: how do boundaries and surroundings shape a measurement inside a region?", 21),
        text(55, 993, "Local round spheres supply areal radius and radial propagation; the many-body parent may be nonspherical.", 18),
        panel(45, 1040, 1350, 197, "#302c26", "#b49363"),
        text(75, 1086, "OPEN CONNECTIONS AND ALTERNATIVE ROUTES", 23, "#edbf8d", "600"),
        text(75, 1130, "Exterior preparation and scattering  |  A minimally coupled field with a tail  |  Quantum state and backreaction", 18),
        text(75, 1172, "Multiple daughters  |  Communicating mouths  |  The glued parent's global future", 20),
        text(75, 1210, "A transient here does not rule out memory in another field, state or geometry.", 18),
        text(55, 1300, "Conditional classical response. No global parent signal, persistent dark component or cosmic age is inferred.", 18, "#b7c9d5"),
        "</svg>"]
    return "\n".join(out)+"\n"


def narrative():
    return """# From a causal connection to an interior field measurement

![Field-conditioned map](map.svg)

The new connection carries a classical massless conformal scalar in its
spherical mode. The fixed wall's curvature distribution is included. It
gives continuity of the reduced field R phi and its normal derivative;
it does not erase smooth Schwarzschild scattering outside the wall.

The cap begins with zero classical field and time derivative. A single
exterior-side incoming characteristic is prescribed on the wall. Its compact
profile is fixed in wall proper time around the early R/ell=0.5 event.
This is an initial-boundary state, not a derived remote source or quantum
vacuum. The right-parent causal exclusion in atlas v5 remains unchanged.

![Input, clock map and central response](field-response.svg)

The exact conditional transfer is phi_C[T(tau_W)]=I(tau_W)/Tprime(tau_W).
It predicts a finite transient and a proper-time distortion of the supplied
profile. The central field returns to zero after this pulse within the stated
observation window. Squared-amplitude integrals are waveform diagnostics,
not a stress-energy or dark-matter calculation. Exterior tails have not been
computed; zero incoming trace outside the pulse is part of the input.

The common thread is a local observer's response to regional boundaries,
geometry and state. The foundation proposes inheritance; the finite companion
computes transport and memory under its own assumptions; the atlas supplies
physical geometry and causal access. The companion's unresolved mode-space
complement has not been identified with a cosmological parent. Round local
spheres are compatible with a nonspherical multi-hole assembly.

[field-transfer.json](field-transfer.json) records the state, response and
checks; [atlas.json](atlas.json) preserves all earlier scenarios;
[manifest.json](manifest.json) binds their sources. The
[paper v6](../../../../paper/atlas-working-draft-v6.md) develops the argument.
The minimally coupled field with a bulk tail, exterior-state construction,
multiple-daughter and communicating-mouth routes remain open.
"""


def response_figure(response):
    """Compare specified input and received field in their own proper clocks."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    s, support, peaks = response["samples"], response["support"], response["peaks"]
    wall = np.asarray(s["tau_wall_over_ell"])
    centre = np.asarray(s["arrival_t_C_over_ell"])
    incoming = np.asarray(s["ell_I_over_epsilon"])
    received = np.asarray(s["ell_phi_C_over_epsilon"])
    stretch = np.asarray(s["dT_dtau"])
    tau0 = response["incoming_state"]["wall_clock_centre_over_ell"]
    t0 = float(centre[np.argmin(abs(wall-tau0))])
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 10}):
        fig, axes = plt.subplots(2, 2, figsize=(14, 9), layout="constrained")
        ax = axes[0, 0]
        ax.plot(wall-tau0, incoming, color="#99662a", linewidth=2.2)
        ax.axhline(0, color="#8a939b", linewidth=.8)
        ax.fill_between(wall-tau0, incoming, alpha=.13, color="#99662a")
        ax.set(xlabel=r"Wall proper time $(\tau_W-\tau_*)/\ell$", ylabel=r"Incoming characteristic $\ell I/\varepsilon$",
            title="Declared input on the exterior wall side")
        ax.text(.03, .05, r"$I=\varepsilon B'(z)/(0.03\ell)$; compact initial-boundary data", transform=ax.transAxes, fontsize=10)
        ax = axes[0, 1]
        ax.plot(centre, received, color="#236e98", linewidth=2.2)
        ax.axhline(0, color="#8a939b", linewidth=.8)
        ax.fill_between(centre, received, alpha=.13, color="#236e98")
        ax.axvline(t0, color="#777777", linestyle=":", linewidth=1)
        ax.set(xlabel=r"Central proper time $t_C/\ell$", ylabel=r"Measured field $\ell\phi_C/\varepsilon$",
            title="Predicted response at the regular centre")
        ratio = response["metrics"]["negative_to_positive_peak_magnitude"]
        ax.text(.03, .05, f"Two unequal lobes; magnitude ratio {ratio:.4f}", transform=ax.transAxes, fontsize=10)
        ax = axes[1, 0]
        ax.plot(wall-tau0, stretch, color="#85529e", linewidth=2.2)
        ax.set(xlabel=r"Wall proper time $(\tau_W-\tau_*)/\ell$", ylabel=r"Arrival stretch $dt_C/d\tau_W$",
            title="One clock map controls duration and amplitude")
        mean_stretch = response["metrics"]["support_stretch_ratio"]
        ax.text(.03, .06, f"Support lasts {mean_stretch:.4f} times the wall pulse duration", transform=ax.transAxes, fontsize=10)
        ax = axes[1, 1]
        start, end = support["centre_arrival_t_over_ell"]
        xx = np.concatenate(([start-.10], centre, [end+.03, 0.]))
        yy = np.concatenate(([0.], received, [0., 0.]))
        ax.plot(xx, yy, color="#236e98", linewidth=2)
        ax.axvspan(end, 0, color="#d5e9dc", alpha=.7)
        ax.axhline(0, color="#77877f", linewidth=.8)
        ax.axvline(0, color="#505050", linestyle="--", linewidth=1)
        ax.set(xlabel=r"Central proper time $t_C/\ell$", ylabel=r"Measured field $\ell\phi_C/\varepsilon$",
            title="A finite transient, with no late field in this window")
        ax.text(end+.025, .70*max(received), "Zero after pulse\nthrough the chosen slice", fontsize=10)
        ax.text(.49, .07, "No detector, energy flux\nor backreaction inferred", transform=ax.transAxes, fontsize=9)
        for ax in axes.flat:
            ax.grid(alpha=.16)
        fig.suptitle("An early scalar pulse reaches the interior observer", fontsize=18)
        fig.supxlabel("Fixed conformal field and wall. Incident history is supplied, including zero outside the pulse; global exterior preparation remains open.", fontsize=10)
        out = {}
        for fmt in ("svg", "png"):
            stream = BytesIO()
            metadata = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
            fig.savefig(stream, format=fmt, dpi=155, metadata=metadata)
            out[f"field-response.{fmt}"] = stream.getvalue()
        plt.close(fig)
    return out


def payloads():
    record, response = build_record()
    svg = render_map()
    outputs = {"atlas.json": first.atlas.canonical_bytes(record),
        "field-transfer.json": first.atlas.canonical_bytes(response),
        "map.svg": svg.encode(), "map.png": first.png_preview(svg), "narrative.md": narrative().encode()}
    outputs.update(response_figure(response))
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
        "outputs": {n: {"sha256": first.atlas.hashlib.sha256(b).hexdigest(), "bytes": len(b)} for n, b in outputs.items()},
        "runtime": {"python": platform.python_version(), **{n: version(n) for n in ("numpy", "scipy", "sympy", "matplotlib")}},
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12},
        "scope": "conditional classical field transfer on the isolated reference; incident history is input"}
    directory.mkdir(parents=True, exist_ok=False)
    for n, b in outputs.items():
        with (directory/n).open("xb") as handle:
            handle.write(b)
    with (directory/"manifest.json").open("xb") as handle:
        handle.write(first.atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "remote_source_solved": False}


def check(directory):
    directory = Path(directory)
    m = json.loads((directory/"manifest.json").read_bytes())
    if m.get("schema") != SCHEMA or m["source_bindings"] != source_bindings():
        raise ValueError("atlas source binding mismatch")
    expected = payloads()
    if set(expected) != set(m["outputs"]):
        raise ValueError("atlas output inventory mismatch")
    for n, replay in expected.items():
        path, binding = directory/n, m["outputs"][n]
        if first._sha(path) != binding["sha256"] or path.stat().st_size != binding["bytes"]:
            raise ValueError(f"artifact binding mismatch: {n}")
        equal = first._close(json.loads(path.read_bytes()), json.loads(replay)) if n.endswith(".json") else path.read_bytes() == replay
        if not equal:
            raise ValueError(f"atlas replay mismatch: {n}")
    if m["input_sha256"] != first._sha(directory/"atlas.json"):
        raise ValueError("canonical input binding mismatch")
    return {"checked": directory.name, "outputs": len(expected), "sources": len(m["source_bindings"])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
