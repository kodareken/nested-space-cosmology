#!/usr/bin/env python3
"""Build a source-bound causal atlas around the fixed interior observer."""
from __future__ import annotations

import argparse
from copy import deepcopy
from io import BytesIO
import json
import math
from pathlib import Path
import platform
from importlib.metadata import version

import numpy as np

import derive_nsc_cosmological_atlas_v4 as previous
from recursive_horizons import nsc_atlas_causal_reference as reference
from recursive_horizons import nsc_atlas_causal_observers as observers

first = previous.first
LAB, REPO = previous.LAB, previous.REPO
SCHEMA = "nsc-cosmological-atlas-v5"
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v5"


def observer_access_record(ref, obs):
    events = obs["named_events"]
    threshold = obs["parent_escape_threshold"]
    early_c = events["C_emits_to_early_wall"]["t_flat_over_ell"]
    reception = events["C_receives_early_wall"]["t_flat_over_ell"]
    arrival = ref["events"]["early_wall_parent_arrival"]
    rows = [
        {"id": "early-C-W", "from": "C", "to": "W", "emission": f"C proper time t/ell={early_c:.10f}",
         "verdict": "allowed geometric ray", "reception": "W at R/ell=0.5", "basis": "outgoing de Sitter null intercept and retained-history coverage"},
        {"id": "early-W-C", "from": "W", "to": "C", "emission": "W at R/ell=0.5",
         "verdict": "allowed geometric ray", "reception": f"C proper time t/ell={reception:.10f}", "basis": "ingoing eta+r constant and negative"},
        {"id": "early-C-P", "from": "C", "to": "P", "emission": f"C proper time t/ell={early_c:.10f}, via W at R/ell=0.5",
         "verdict": "allowed geometric path", "reception": f"P proper time/ell={arrival['parent_proper_time_over_ell']:.10f} on its own clock",
         "basis": "outgoing cap ray followed by constant U<0 from white-hole sector to right exterior; transmission law unspecified"},
        {"id": "historical-C-W0", "from": "C", "to": "W", "emission": "C proper time t/ell=-log(2.5)",
         "verdict": "allowed geometric ray", "reception": "W0 at the chosen initial slice",
         "basis": "eta-r=-2.5ell; emission is within the retained matched history"},
        {"id": "C0-W", "from": "C", "to": "W", "emission": "C0 or later, t>=0",
         "verdict": "excluded in reference", "reception": "no future wall event",
         "basis": "remaining outgoing conformal radius<=ell; future wall comoving radius>=1.5ell"},
        {"id": "W0-C", "from": "W", "to": "C", "emission": "W0 or later, R>=1.5ell",
         "verdict": "excluded in reference", "reception": "no finite future central event",
         "basis": "inward null intercept eta+r>0 lies beyond future de Sitter boundary"},
        {"id": "C0-P", "from": "C", "to": "P", "emission": "C0 or later, t>=0",
         "verdict": "excluded in reference", "reception": "no future right-parent event",
         "basis": "no future exit through the enclosing wall; late wall also has U>0"},
        {"id": "P-C", "from": "P", "to": "C", "emission": "any event on the right static parent worldline",
         "verdict": "excluded in reference", "reception": "no central event via this attachment",
         "basis": "parent V>0; wall V<0; future-directed causal V cannot decrease"},
        {"id": "last-C-P", "from": "C", "to": "P", "emission": f"threshold t/ell={threshold['centre_emission_t_over_ell']:.10f}",
         "verdict": "boundary, not finite reception", "reception": "horizon generator; equality is excluded",
         "basis": "outgoing cap ray meets W at R=2M and continues as U=0; finite parent access requires earlier emission"},
        {"id": "siblings", "from": "A interior", "to": "B or descendants", "emission": "unspecified",
         "verdict": "unresolved", "reception": "no global many-body future supplied",
         "basis": "initial-data gluing alone gives no inter-body causal completion"},
    ]
    for row in rows:
        row["scope"] = "glued multi-body geometry" if row["id"] == "siblings" else "isolated fixed spherical reference"
    return {"schema": "nsc-atlas-observer-access-v1", "observers": deepcopy(obs["observers"]),
        "clock_conventions": deepcopy(obs["conventions"]), "rows": rows,
        "centre_escape_threshold": deepcopy(threshold),
        "initial_fields_or_messages_prescribed": False, "signal_amplitudes_calculated": False,
        "field_boundary_condition": "a causal geometric path does not select a material transmission or quantum state",
        "many_body_scope": "only the prior v4 smooth-vacuum protected-domain statement is transferred; no global reference causal verdict is transferred"}


def build_record():
    record, _, _, _ = previous.build_record()
    ref, obs = reference.causal_reference_record(), observers.causal_observers_record()
    access = observer_access_record(ref, obs)
    record["sources"].update({
        "reference_causality": {"path": "docs/nsc-atlas-causal-classification.md", "use": "known signed-wall completion and flat exterior limit"},
        "interior_observers": {"path": "docs/nsc-atlas-causal-observers.md", "use": "event-specific interior null access"},
    })
    region, link = first.atlas._region, first.atlas._link
    regions = [region(identifier, title, "observer_domain", {"chart": chart, "worldline": worldline}, source,
        "classical reference; signal state not prescribed", {"kind": clock}, {"H_ell": 1., "k_ell": 2., "M_over_ell": .4})
        for identifier, title, chart, worldline, source, clock in (
            ("C", "Central interior observer", "expanding de Sitter planar", "r_comoving=0", "positive vacuum energy", "central proper time t_flat"),
            ("W", "Moving wall observer", "matched worldsheet", "fixed signed expanding shell", "fixed positive tension", "wall proper time tau_wall"),
            ("P", "Right parent observer", "right Schwarzschild exterior I", "R=4M", "vacuum", "sqrt(1-2M/R) times right Schwarzschild T"))]
    links = []
    for a, b, status, evidence in (
        ("C", "W", "reference_event_conditioned", "historical null witness; no emission at or after C0 reaches future wall"),
        ("W", "C", "reference_event_conditioned", "early wall witness; no wall emission at or after W0 reaches centre"),
        ("C", "P", "reference_event_conditioned", "early path via white-hole wall; strict escape threshold; no escape for emissions at or after C0"),
        ("P", "C", "reference_excluded", "future V monotonicity forbids reaching the matched wall")):
        links.append(link(f"reference-{a}-{b}", a, b, "causal", ["reference_causality", "interior_observers"], evidence,
            status=status, missing=["field transmission and emitted state for an operational signal"],
            delay=None, interval=None, clock=None))
    record["scenarios"].append({"id": "interior_observer_reference", "title": "Interior observer in the isolated causal reference",
        "theory": {"family": "fixed de Sitter/Schwarzschild positive-tension thin wall", "couplings": {"G": 1, "c": 1, "H_ell": 1, "k_ell": 2}},
        "assumptions": ["isolated spherical reference, not a many-body evolution", "published singular-past monotonic completion",
            "event-specific geometric access; no generic inter-region delay", "independent proper clocks with explicit origins"],
        "regions": regions, "links": links, "max_descent_depth": 0,
        "solution_status": "borrowed reference completion with checked null paths and analytic exclusions",
        "access_record": "observer-access.json"})
    record["assessments"].update({"interior_observer_reference": {"record": "causal-map.json", "status": ref["status"]},
        "observer_access": {"record": "observer-access.json", "global_many_body_future_solved": False}})
    record["admission"].update(category="reference causal completion and interior observer access",
        missing_arrow="field/state data carrying an early imprint, or a separately justified alternative communicating attachment")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    return record, {"reference": ref, "interior": obs}, access


def source_bindings():
    bindings = previous.source_bindings()
    paths = [Path(__file__), Path(reference.__file__), Path(observers.__file__),
        LAB / "tests/test_nsc_atlas_causal_reference.py", LAB / "tests/test_nsc_atlas_causal_observers.py",
        LAB / "tests/test_nsc_cosmological_atlas_v5.py", LAB / "docs/nsc-atlas-causal-classification.md",
        LAB / "docs/nsc-atlas-causal-observers.md", LAB / "docs/nsc-cosmological-atlas-v5.md",
        REPO / "paper/atlas-working-draft-v5.md"]
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def access_markdown(access):
    lines = ["# Observer access in the fixed isolated reference", "",
        "C is the central observer; W is the wall; P is the right static parent observer.",
        "All times use the named observer's clock. Allowed paths are geometric test-light possibilities; field transmission and signal amplitudes are unspecified.",
        "Earlier emissions may arrive after the initial slice; exclusions concern the emission events listed below.", "",
        "| From → to | Emission | Verdict / reception | Mathematical basis |",
        "| --- | --- | --- | --- |"]
    for row in access["rows"]:
        lines.append(f"| {row['from']} → {row['to']} | {row['emission']} | {row['verdict']}; {row['reception']} | {row['basis']} |")
    lines += ["", "The B/descendant row remains unresolved. The other rows belong only to the isolated spherical reference.",
        "Its singular past is not a regular formation history. The glued parent's global future is not supplied.", ""]
    return "\n".join(lines)


def render_map(access):
    tcrit = access["centre_escape_threshold"]["centre_emission_t_over_ell"]
    text, panel = first._text, first._panel
    out = first._base(1440, 1240, "What can an interior observer know?",
        "Atlas v5 | Event-specific causal access in one borrowed reference history")
    out += [panel(45, 130, 1350, 118, "#153329", "#70ac94"),
        text(75, 173, "The same connection has different causal meanings at different times.", 25, "#a5e6c9", "600"),
        text(75, 216, "Fixed source and initial slice. The isolated reference has a singular past; regular formation remains open.", 19)]
    for x, title, color, lines in (
        (45, "EARLIER REFERENCE HISTORY", "#92dfc4", ["Interior centre -> wall: allowed null paths", "Early wall -> centre: allowed null paths",
            "Early centre -> parent: a geometric escape path", f"Last central escape threshold: t / ell = {tcrit:.6f}", "The threshold itself has no finite parent reception."]),
        (745, "EMISSIONS AT / AFTER THE SLICE", "#edbf8d", ["Centre -> future wall: excluded in this branch", "Wall -> future centre: excluded in this branch",
            "Centre -> right parent: excluded in this branch", "Right parent -> centre: excluded throughout reference", "These are observer/event statements, not generic sibling rules."])):
        out += [panel(x, 285, 650, 305, "#152c40", "#567f96"), text(x+25, 332, title, 22, color, "600")]
        for i, line in enumerate(lines):
            out.append(text(x+25, 376+43*i, line, 17 if i > 2 else 19))
    out += [panel(45, 635, 1350, 180, "#2a2637", "#9c82b7"),
        text(75, 680, "TWO DIFFERENT HORIZONS", 24, "#d0b5f0", "600"),
        text(75, 722, "The right-parent escape horizon crosses the wall at R = 2M and continues through the cap.", 21),
        text(75, 763, "The central observer's de Sitter horizon limits incoming light. They are different null surfaces.", 21),
        text(55, 880, "WHAT THIS UNLOCKS", 25, "#a1d9f3", "600"),
        text(55, 922, "An early imprint needs a specified field, state and wall transmission law before its local response can be calculated.", 19),
        text(55, 965, "Earlier signals may arrive later. These exclusions concern new emissions; other candidate geometries stay open.", 19)]
    for x, title, subtitle in ((45, "Multiple daughters", "One-hole branching: distinct source assumptions"),
        (505, "Communicating mouths", "Separate charged, semiclassical route"),
        (965, "Shared many-body future", "Still undetermined beyond certified domains")):
        out += [panel(x, 1010, 430, 112), text(x+20, 1050, title, 21, "#c9b2ec", "600"), text(x+20, 1090, subtitle, 16)]
    out += [text(55, 1181, "Spatial attachment, descent, causal influence and correlation remain distinct. No inherited constants or cosmic ages inferred.", 17, "#b7c9d5"), "</svg>"]
    return "\n".join(out)+"\n"


def spacetime_figure(causal):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import Polygon
    ref, obs = causal["reference"], causal["interior"]
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 10}):
        fig, axes = plt.subplots(1, 2, figsize=(16, 9), layout="constrained", gridspec_kw={"width_ratios": [1, 1.45]})
        ax = axes[0]
        samples, e = obs["trajectory"]["samples"], obs["named_events"]
        eta, r = np.asarray(samples["eta_over_ell"]), np.asarray(samples["r_wall_over_ell"])
        ax.fill_betweenx(eta, 0, r, color="#e1f0e8", alpha=.85)
        ax.plot(r, eta, color="#a3682b", linewidth=2.8, label="W: matched wall (finite samples)")
        ax.plot([0, 0], [eta[0], 0], color="#246184", linewidth=3, label="C: central proper-time observer")
        ax.plot([0, 1.5], [-1, -1], color="#505050", linewidth=1.8, linestyle="--", label="v4 initial slice")
        # The two horizons are distinct exact radial null surfaces.
        wh = e["W_Schwarzschild_horizon"]
        eh = obs["parent_escape_threshold"]["eta_minus_r_over_ell"]
        ax.plot([0, wh["r_wall_over_ell"]], [eh, wh["eta_over_ell"]], color="#85529e", linewidth=2, linestyle="-.", label="Escape horizon to right parent")
        dh = e["W_dS_horizon"]
        ax.plot([0, dh["r_wall_over_ell"]], [0, dh["eta_over_ell"]], color="#666666", linewidth=1.5, linestyle=":", label="C's incoming de Sitter horizon")
        early = e["W_early"]
        ce, cr = e["C_emits_to_early_wall"], e["C_receives_early_wall"]
        ax.annotate("", xy=(early["r_wall_over_ell"], early["eta_over_ell"]), xytext=(0, ce["eta_over_ell"]),
                    arrowprops={"arrowstyle": "->", "color": "#267b56", "linewidth": 2})
        ax.annotate("", xy=(0, cr["eta_over_ell"]), xytext=(early["r_wall_over_ell"], early["eta_over_ell"]),
                    arrowprops={"arrowstyle": "->", "color": "#267b56", "linewidth": 2})
        for name, label, offset in (("W_early", "W early: R=0.5 ell", (10, -19)),
                ("W_Schwarzschild_horizon", "Parent horizon\nR=2M", (10, -20)),
                ("W_dS_horizon", "R=ell", (12, -5)), ("W0", "W0: R=1.5 ell", (15, 5))):
            p = e[name]
            ax.scatter([p["r_wall_over_ell"]], [p["eta_over_ell"]], s=22, color="#a3682b", zorder=5)
            ax.annotate(label, (p["r_wall_over_ell"], p["eta_over_ell"]), xytext=offset, textcoords="offset points", fontsize=9)
        ax.scatter([0, 0], [-1, eh], color="#246184", s=30, zorder=5)
        ax.annotate("C0", (0, -1), xytext=(8, -15), textcoords="offset points", fontsize=10)
        ax.annotate("Escape threshold\n(no finite P reception)", (0, eh), xytext=(8, 10), textcoords="offset points", fontsize=9, color="#85529e")
        ax.axhline(0, color="#444444", linewidth=1.5)
        ax.text(.12, .08, "Future boundary: infinite C proper time", fontsize=9)
        ax.annotate("Wall continues expanding", (r[-1], eta[-1]), xytext=(.2, -.3),
                    fontsize=9, arrowprops={"arrowstyle": "->", "color": "#777777"})
        ax.text(.12, eta[0]-.25, "Earlier singular endpoint below sample window", fontsize=8, color="#7d4444")
        ax.set(xlim=(-.16, 2.9), ylim=(eta[0]-.48, .3), xlabel=r"Interior comoving radius $r/\ell$",
               ylabel=r"Interior conformal time $\eta/\ell$", title="Interior chart: light follows 45-degree lines")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.085), fontsize=8, framealpha=.92)
        ax.grid(alpha=.12)
        ax = axes[1]
        p = math.pi/4
        outline = [(-2*p, 0), (-p, -p), (p, -p), (2*p, 0), (p, p), (-p, p)]
        ax.add_patch(Polygon(outline, facecolor="#edf1f5", edgecolor="#6c8091", linewidth=1.5))
        wall = ref["wall_samples"]
        X, T = np.asarray(wall["X"]), np.asarray(wall["T"])
        radii = np.asarray(wall["R_over_ell"])
        ih = int(np.flatnonzero(radii == .8)[0])
        idh = int(np.flatnonzero(radii == 1.)[0])
        # All changes to the retained Schwarzschild region lie left of this
        # timelike wall; shaded wedge is the replaced Schwarzschild side.
        left_boundary = np.abs(T)-math.pi/2
        ax.fill_betweenx(T, left_boundary, X, color="#f6e7d8", alpha=.9)
        ax.plot(X, T, color="#a3682b", linewidth=2.8, label="Same W, Schwarzschild chart")
        ax.scatter(X[[0, -1]], T[[0, -1]], facecolors="white", edgecolors="#a3682b", s=30, zorder=6)
        ax.annotate("R=6 ell; curve continues", (X[-1], T[-1]), xytext=(-7, 28), textcoords="offset points", fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": "#a3682b"})
        ax.annotate("R=0.02 ell sample", (X[0], T[0]), xytext=(-93, -34), textcoords="offset points", fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": "#a3682b"})
        ax.plot([-p, p], [p, -p], color="#8a6aa5", linewidth=1.1, linestyle=":")
        # U=0: parent event-horizon generator, continuous from the cap panel.
        ax.plot([X[ih], p], [T[ih], p], color="#85529e", linestyle="-.", linewidth=1.8, label="U=0 horizon to parent (retained)")
        # Continue C's incoming horizon through W(R=ell), at constant V.
        vh = wall["V"][idh]
        uh = np.linspace(1/vh, wall["U"][idh], 100)
        ax.plot((np.arctan(vh)-np.arctan(uh))/2, (np.arctan(vh)+np.arctan(uh))/2,
                color="#666666", linewidth=1.5, linestyle=":", label="C's incoming horizon (same V)")
        ax.plot([-p, p], [p, p], color="#8b4141", linewidth=2.8)
        ax.plot([-p, p], [-p, -p], color="#a4a4a4", linewidth=1.2, linestyle="--")
        ax.text(-.55, p+.025, "future R=0 singular boundary", fontsize=9)
        ax.text(-.32, -p-.15, "Past R=0 chart boundary\n(left part excised)", fontsize=8)
        initial = ref["initial_slice"]
        ax.plot(initial["X"], initial["T"], color="#505050", linewidth=1.6, linestyle="--", label="v4 initial slice")
        tt = np.linspace(-7, 7, 161)
        f = math.sqrt(math.exp(2))
        U, V = -f*np.exp(-tt/1.6), f*np.exp(tt/1.6)
        ax.plot((np.arctan(V)-np.arctan(U))/2, (np.arctan(V)+np.arctan(U))/2,
                color="#246184", linewidth=2.4, label="P: right static observer, R=4M")
        w, arrival = ref["events"]["early_wall"], ref["events"]["early_wall_parent_arrival"]
        ax.annotate("", xy=(arrival["X"], arrival["T"]), xytext=(w["X"], w["T"]),
            arrowprops={"arrowstyle": "->", "color": "#267b56", "linewidth": 2})
        now = ref["events"]["cut_wall"]
        for event, label, xytext in ((w, "W early", (-47, 11)), (now, "W0", (-29, 9)), (arrival, "early ray reaches P", (8, 9))):
            ax.scatter([event["X"]], [event["T"]], s=22, color="#a3682b" if event != arrival else "#267b56", zorder=5)
            ax.annotate(label, (event["X"], event["T"]), xytext=xytext, textcoords="offset points", fontsize=9)
        for x, y, label in ((.96, -.25, "I: right parent"), (-.42, .3, "III: left exterior"), (-.24, -.48, "IV: white hole"), (.05, .58, "II: black hole")):
            ax.text(x, y, label, fontsize=8.5, ha="center", color="#536575")
        ax.text(-1.42, -.05, "Replaced\nby cap", fontsize=9, color="#9a6d42")
        ax.set(xlim=(-1.7, 1.7), ylim=(-1.02, 1.02), xlabel="Compact Kruskal space X", ylabel="Compact Kruskal time T",
               title="Parent chart: the spatial bridge is not traversable")
        ax.set_aspect("equal", adjustable="box")
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.15), fontsize=8, ncol=1)
        fig.suptitle("One interior, three observers, two distinct horizons", fontsize=18)
        fig.supxlabel("Matched reference charts share W events, not coordinate clocks. Finite wall curves; endpoint classification is analytic. The glued many-body future remains open.", fontsize=9)
        out = {}
        for fmt in ("svg", "png"):
            stream = BytesIO()
            meta = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
            fig.savefig(stream, format=fmt, dpi=155, metadata=meta)
            out[f"spacetime.{fmt}"] = stream.getvalue()
        plt.close(fig)
    return out


def narrative(causal, access):
    threshold = access["centre_escape_threshold"]["centre_emission_t_over_ell"]
    return f"""# The geometry seen by an interior observer

![Causal map](map.svg)

The fixed wall belongs to a published monotonic expanding reference branch.
Its past reaches a Schwarzschild singular boundary; its late interior expands.
This is a causal completion of the isolated spherical body, not a regular
formation mechanism or a computed evolution of the glued many-body parent.

![Matched spacetime charts](spacetime.svg)

The centre C, wall W and right static parent P carry different proper clocks.
The panels share wall events while retaining their own coordinates. The
initial slice matches v4. Horizon curves stop at their matching wall events;
the dotted central horizon continues across both panels. The faint dashed
past Schwarzschild boundary is uncut chart context, partly excised by the cap.
Green rays illustrate earlier geometric access:
light can travel from C to the early wall, from the early wall back to C,
and from the early wall through the white-hole sector to P. A field's actual
transmission and amplitude require additional matter/state data.

The central escape threshold to the right parent is
t_C/ell={threshold:.12f}. Earlier emissions within the certified history can
escape; equality follows a horizon generator and does not reach P at finite
proper time. The parent's horizon continues through the cap and is not a
constant-radius sphere there. It differs from C's incoming de Sitter horizon.

Emissions at or after the chosen slice cannot travel from C to the future wall,
from the wall to C, or from C to P in this fixed branch. Earlier wall emissions
can reach C arbitrarily late; these exclusions do not erase incoming history.
P cannot reach the wall or C anywhere in this reference history. These are
specific observer and emission-domain conclusions, not a generic rule for all
black holes, siblings or nested geometries.

[Observer access table](observer-access.md) records every verdict and basis.
[causal-map.json](causal-map.json), [observer-access.json](observer-access.json),
[atlas.json](atlas.json) and [manifest.json](manifest.json) retain exact inputs
and provenance. The [paper v5](../../../../paper/atlas-working-draft-v5.md)
connects these results to the prior spatial atlas.

The shared parent's future, second cosmological child and descendants remain
unresolved. Multiple-daughter, communicating-mouth and exterior-response
models remain separate routes. The next response question can concern a
specified early field imprint or a different attachment with its own causal
channel. No late-parent driving, inherited constant or cosmic age is inferred
from this reference. Known junction and causality mathematics is borrowed;
the new work closes a project-level observer map.
"""


def payloads():
    record, causal, access = build_record()
    svg = render_map(access)
    outputs = {"atlas.json": first.atlas.canonical_bytes(record), "causal-map.json": first.atlas.canonical_bytes(causal),
        "observer-access.json": first.atlas.canonical_bytes(access), "observer-access.md": access_markdown(access).encode(),
        "map.svg": svg.encode(), "map.png": first.png_preview(svg), "narrative.md": narrative(causal, access).encode()}
    outputs.update(spacetime_figure(causal))
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
        "scope": "isolated reference history and event-specific causal geometry, not global many-body evolution"}
    directory.mkdir(parents=True, exist_ok=False)
    for n, b in outputs.items():
        with (directory / n).open("xb") as handle:
            handle.write(b)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(first.atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "new_physical_mechanism": False}


def check(directory):
    directory = Path(directory)
    m = json.loads((directory / "manifest.json").read_bytes())
    if m.get("schema") != SCHEMA or m["source_bindings"] != source_bindings():
        raise ValueError("atlas source binding mismatch")
    expected = payloads()
    if set(expected) != set(m["outputs"]):
        raise ValueError("atlas output inventory mismatch")
    for n, replay in expected.items():
        path, binding = directory / n, m["outputs"][n]
        if first._sha(path) != binding["sha256"] or path.stat().st_size != binding["bytes"]:
            raise ValueError(f"artifact binding mismatch: {n}")
        equal = first._close(json.loads(path.read_bytes()), json.loads(replay)) if n.endswith(".json") else path.read_bytes() == replay
        if not equal:
            raise ValueError(f"atlas replay mismatch: {n}")
    if m["input_sha256"] != first._sha(directory / "atlas.json"):
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
