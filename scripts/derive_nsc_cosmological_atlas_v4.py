#!/usr/bin/env python3
"""Create or replay the first sourced interior connection and observer domain."""
from __future__ import annotations

import argparse
from copy import deepcopy
from io import BytesIO
import json
from pathlib import Path
import platform
from importlib.metadata import version

import numpy as np

import derive_nsc_cosmological_atlas_v3 as previous
from recursive_horizons import nsc_atlas_connection_slice as slices
from recursive_horizons import nsc_atlas_connection_observer as observers

first = previous.first
LAB, REPO = previous.LAB, previous.REPO
SCHEMA = "nsc-cosmological-atlas-v4"
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v4"


def parent_assembly_record(connection):
    observers.require_compatible_slice(connection)
    return {
        "schema": "nsc-atlas-parent-assembly-v1",
        "status": "END_ONLY_GLUING_WITH_VERIFIED_SOURCED_BODY",
        "source": {"url": "https://arxiv.org/html/1004.1355v1",
                   "read": "Theorem1, following multiple-end discussion and proof support of deformations",
                   "owner": "docs/nsc-atlas-connection-gluing.md"},
        "admission": {"borrowed": "smooth vacuum AE-end gluing plus locality of constraints",
            "missing": "how to preserve a verified distributional body while using smooth gluing inputs",
            "decides": "the theorem acts only on its smooth vacuum end; the body is reattached on an identical open overlap",
            "unlocks": "a shared initial-data parent containing the first sourced interior",
            "finite_domain": "two selected ends and prescribed retained regions",
            "budget": "theorem application only, no new elliptic solve",
            "stop": "no smooth vacuum end or no identical open overlap"},
        "bodies": [
            {"id": "A", "mass_over_ell": .4, "source": "verified fixed-tension wall and de Sitter cap",
             "record": "connection-slice.json", "selected_end": "right Schwarzschild AE end",
             "retained_right_R_cut_over_ell": 4., "open_overlap_beyond_cut_required": True,
             "other_AE_ends": 0},
            {"id": "B", "mass_over_ell": .6, "source": "complete time-symmetric vacuum Schwarzschild bridge",
             "selected_end": "one chosen AE end", "other_AE_ends": 1,
             "cosmological_child_attached": False}],
        "hypotheses": {"selected_inputs_smooth_vacuum_AE": True,
            "ADM_four_momenta": "(M_A,0,0,0), (M_B,0,0,0): future timelike",
            "end_deformations_only": True,
            "source_remainder_class": "piecewise smooth h,K satisfying distributional constraints at A wall",
            "open_overlap": "identical h,K, hence no extra jump or matter layer at reattachment",
            "whole_body_smooth_corollary_used": False},
        "topology": {"designated_common_outer_AE_ends": 1, "retained_other_AE_ends": 1,
            "total_AE_ends": 2, "second_cosmological_child_established": False},
        "separation": {"domain": "sufficiently large under the published existence theorem",
            "numeric_threshold": None, "BL_d12_certified": False},
        "conclusions": {"shared_initial_data_existence": True,
            "explicit_connector_metric": False, "unchanged_BL_parent": False,
            "one_sourced_interior_connected": True, "formation_history": False,
            "global_shell_evolution": False, "perturbation_stability": False,
            "sibling_signalling": False, "new_to_world_claim": False},
    }


def build_record():
    record, _, _ = previous.build_record()
    connection = slices.connection_slice_record()
    assembly = parent_assembly_record(connection)
    observer = observers.causal_window_record(connection)
    record["sources"].update({
        "connected_slice": {"path": "docs/nsc-atlas-connection-slice.md", "use": "explicit sourced initial data"},
        "end_only_gluing": deepcopy(assembly["source"]),
        "vacuum_observer": {"path": "docs/nsc-atlas-connection-observer.md", "use": "finite vacuum domain of dependence"},
    })
    region, link = first.atlas._region, first.atlas._link
    regions = [
        region("parent", "Common outer parent", "parent",
            {"chart": "many-body existence family", "globally_spherical": False},
            "vacuum connector", "classical initial data", {"kind": "initial slice"}, {}),
        region("bridge_a", "A: vacuum bridge and retained exterior", "interface",
            {"chart": "explicit intrinsic spherical slice and exact right Schwarzschild end", "symmetry_orbits": "S2"},
            "vacuum", "fixed mass M_A/ell=0.4",
            {"kind": "local static parent observer", "proper_window_over_ell": [0, observer["window"]["observer_proper_duration_over_ell"]],
             "upper_endpoint_included": False}, {"M_over_ell": .4}),
        region("interior_a", "A: regular expanding de Sitter cap", "interior",
            {"chart": "finite flat spatial ball", "symmetry_orbits": "S2", "closed_S3_asserted": False},
            "positive vacuum energy; fixed-tension boundary", "classical cap at R_wall/ell=1.5",
            {"kind": "local de Sitter slice and wall proper time; no shared age"}, {"H_ell": 1., "k_ell": 2.}),
        region("bridge_b", "B: unequal vacuum body", "interface",
            {"chart": "complete Schwarzschild bridge; other end retained", "symmetry_orbits": "S2"},
            "vacuum", "time-symmetric", {"kind": "initial slice"}, {"M_over_ell": .6}),
    ]
    links = [
        link("parent-a", "parent", "bridge_a", "attachment", ["end_only_gluing"],
            "Published end-only gluing with identical open vacuum overlap",
            status="justified_initial_data_existence", missing=["explicit connector metric and numerical separation threshold", "global evolution"]),
        link("a-interior", "bridge_a", "interior_a", "attachment", ["connected_slice"],
            "Explicit common slice with bulk constraints and wall delta constraints checked",
            status="constructed_sourced_initial_slice", missing=["formation history", "global evolution and stability"]),
        link("parent-b", "parent", "bridge_b", "attachment", ["end_only_gluing"],
            "Unequal smooth vacuum end under the same gluing theorem",
            status="justified_initial_data_existence", missing=["cosmological child attachment to B"]),
    ]
    record["scenarios"].append({"id": "first_connected_interior", "title": "First sourced interior in a common parent initial-data family",
        "theory": {"family": "Einstein constraints with fixed positive vacuum energy and pure tension in A; vacuum elsewhere",
                   "couplings": {"G": 1, "c": 1, "H_ell": 1, "k_ell": 2}},
        "assumptions": ["fixed source and shell state from v3", "sufficiently large separation",
            "end-only gluing; exact open overlaps", "B other AE end is retained, not a formed child"],
        "regions": regions, "links": links, "max_descent_depth": 0,
        "solution_status": "explicit first body plus justified common-parent initial-data existence; spacetime formation open"})
    record["assessments"].update({
        "first_connected_interior": {"record": "connection-slice.json", "status": connection["status"]},
        "end_only_parent_assembly": {"record": "parent-assembly.json", "status": assembly["status"]},
        "finite_observer_degeneracy": {"record": "causal-window.json", "status": observer["status"]},
    })
    record["historical_scenario_status"] = {
        "moving_vacuum_wall": "v3 local evidence retained; first_connected_interior supplies one initial slice successor",
        "many_body_spherical_collars": "v3 open-interface record retained; end-only_parent_assembly supplies the narrower application",
        "other_candidates": "independent assumptions and unresolved links unchanged"}
    record["admission"].update(category="first sourced initial-data connection and finite observer equivalence",
        missing_arrow="global causal completion and discriminating exterior response beyond the protected vacuum domain")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    return record, connection, assembly, observer


def source_bindings():
    bindings = previous.source_bindings()
    paths = [Path(__file__), Path(slices.__file__), Path(observers.__file__),
        LAB / "tests/test_nsc_atlas_connection_slice.py", LAB / "tests/test_nsc_atlas_connection_observer.py",
        LAB / "tests/test_nsc_cosmological_atlas_v4.py", LAB / "docs/nsc-atlas-connection-slice.md",
        LAB / "docs/nsc-atlas-connection-observer.md", LAB / "docs/nsc-atlas-connection-gluing.md",
        LAB / "docs/nsc-cosmological-atlas-v4.md", REPO / "paper/atlas-working-draft-v4.md"]
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def render_map(observer):
    text, panel = first._text, first._panel
    out = first._base(1440, 1460, "The first interior-to-parent connection",
        "Atlas v4 | A verified sourced body, a borrowed parent construction and a finite observer domain")
    out += [panel(45, 135, 1350, 250, "#112d42", "#6397b1"),
        text(78, 181, "COMMON OUTER PARENT", 25, "#a2e0f4", "600"),
        text(78, 219, "Two selected vacuum ends joined by the published many-body existence theorem.", 21),
        text(78, 254, "No global spherical symmetry; no numerical separation threshold or explicit connector metric.", 18, "#b9cedc"),
        text(78, 313, "A: sourced body   M / ell = 0.4", 23, "#92dfc4", "600"),
        text(810, 313, "B: vacuum body   M / ell = 0.6", 23, "#efca89", "600"),
        text(78, 350, "End-only gluing leaves the verified wall untouched.", 18),
        text(810, 350, "Its other end stays; child still unresolved.", 18)]
    out += [panel(45, 425, 1350, 335, "#163730", "#6fb59a"),
        text(78, 468, "A: AN EXPLICIT SPACELIKE CONNECTION", 24, "#92dfc4", "600")]
    centers = [200, 525, 850, 1175]
    for i, (cx, title, sub, radius) in enumerate(zip(centers,
            ["Regular cap", "Moving wall", "Vacuum bridge", "Parent end"],
            ["flat de Sitter ball", "fixed positive tension", "minimum R / ell = 0.8", "exact Schwarzschild"],
            [51, 43, 26, 60])):
        if i:
            out.append(f'<path d="M {centers[i-1]+65} 560 L {cx-65} 560" stroke="#86c5ad" stroke-width="3"/>')
        out.append(f'<circle cx="{cx}" cy="560" r="{radius}" fill="#0b2421" stroke="#92dfc4" stroke-width="3"/>')
        out += [text(cx-110, 652, title, 23, "#edf6f1", "600"), text(cx-110, 686, sub, 17, "#b9dace")]
    out += [text(78, 729, "Bulk constraints + shell energy and momentum jumps checked. Diagrammatic circles are not an embedding.", 18),
        panel(45, 800, 650, 253, "#162f49", "#70a8d0"),
        text(72, 843, "WHAT A PARENT OBSERVER CAN INFER", 21, "#a2d9ff", "600"),
        text(72, 883, "Retain 2M < R < 10M; observe at R = 4M.", 20),
        text(72, 923, f"Same local vacuum metric for 0 <= tau / ell < {observer['window']['observer_proper_duration_over_ell']:.6f}.", 20),
        text(72, 964, "Different distant assemblies are indistinguishable", 19),
        text(72, 997, "inside this domain. This does not predict an arrival.", 19),
        panel(735, 800, 660, 253, "#352d23", "#c6a16c"),
        text(762, 843, "WHAT REMAINS A RESEARCH QUESTION", 21, "#f0ca8c", "600"),
        text(762, 884, "Formation, causal completion and stability", 21),
        text(762, 924, "Second child, descendants and state inheritance", 20),
        text(762, 964, "A distinguishing observable beyond this window", 19),
        text(762, 1004, "No new descent or sibling signal is asserted.", 19),
        text(55, 1110, "OTHER ROUTES STAY OPEN", 25, "#a2e0f4", "600")]
    for x, title, lines in ((45, "Multiple daughters", ["One-hole branching geometry", "formation microphysics unresolved"]),
            (505, "Communicating mouths", ["Charged semiclassical constructions", "separate states and assumptions"]),
            (965, "Exterior response", ["Finite memory and equivalence maps", "physical coupling still required"])):
        out += [panel(x, 1143, 430, 156), text(x+22, 1182, title, 22, "#c9aeea", "600"),
                text(x+22, 1222, lines[0], 17), text(x+22, 1258, lines[1], 17, "#afbfca")]
    out += [text(55, 1353, "Attachment / descent / causal influence / correlation are different relations.", 22, "#a2e0f4"),
        text(55, 1394, "An initial-data connection is established in one model. It is not a formation history or a cosmic matter budget.", 18, "#b8cbd8"), "</svg>"]
    return "\n".join(out)+"\n"


def scientific_figure(connection, observer):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 11}):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5.9), layout="constrained")
        ax = axes[0]
        for sample, label, color in ((connection["cap"], "Regular de Sitter cap", "#267858"),
            (connection["transition"]["samples"], "Vacuum slicing transition", "#b16c24"),
            (connection["bridge"]["samples"], "Schwarzschild bridge + parent end", "#245c89")):
            ax.plot(sample["l_over_ell"], sample["R_over_ell"], label=label, color=color, linewidth=2.8)
        wall, bridge = connection["wall"], connection["bridge"]
        ax.scatter([wall["l_over_ell"], bridge["throat_l_over_ell"]], [wall["R_over_ell"], bridge["throat_R_over_ell"]], color="#7d343e", zorder=4)
        ax.annotate("Wall: fixed tension", (1.5, 1.5), xytext=(.2, 2.8), arrowprops={"arrowstyle": "->", "color": "#555555"})
        ax.annotate("Regular minimum sphere", (bridge["throat_l_over_ell"], .8), xytext=(3.5, .22), arrowprops={"arrowstyle": "->", "color": "#555555"}, fontsize=10)
        ax.set(xlabel=r"Proper radial distance $l/\ell$ from the centre", ylabel=r"Areal radius $R/\ell$",
               title="A single connected initial slice", ylim=(0, 6.4))
        ax.legend(loc="upper left", fontsize=9)
        ax.text(.03, .58, "Intrinsic profile, not a Euclidean embedding\nThe exact vacuum end continues to infinity", transform=ax.transAxes, fontsize=9, color="#444444")
        ax.grid(alpha=.18)
        ax = axes[1]
        M = .4
        xo, xc = observers.tortoise(1.6, M), observers.tortoise(4., M)
        bound = observer["window"]["coordinate_duration_over_ell"]
        xmin, xmax, ymax = xo-bound-.5, xc+.4, bound*1.2
        t = np.linspace(0, ymax, 101)
        boundary = xc-t
        ax.fill_betweenx(t, xmin, boundary, color="#e1eee8", label="Protected vacuum domain")
        ax.fill_betweenx(t, boundary, xmax, color="#f5e8d6", label="Far data may matter")
        ax.plot(boundary, t, color="#a16f25", linewidth=2)
        ax.plot([xo, xo], [0, bound], color="#245c89", linewidth=3, label="Static parent observer")
        ax.plot([xo, xo], [bound, ymax], color="#245c89", linewidth=1.5, linestyle="--")
        ax.scatter([xo], [bound], facecolors="white", edgecolors="#245c89", s=65, zorder=4)
        ax.plot([xmin, xc], [0, 0], color="#267858", linewidth=4)
        ax.annotate("Open certificate endpoint", (xo, bound), xytext=(xmin+.3, ymax-.17), fontsize=9,
                    arrowprops={"arrowstyle": "->", "color": "#555555"})
        ax.text(xmin+.2, .22, "Initial vacuum region extends left to r* = -infinity", fontsize=9)
        ax.set(xlim=(xmin, xmax), ylim=(-.12, ymax), xlabel=r"Schwarzschild tortoise coordinate $r_*/\ell$",
               ylabel=r"Parent Schwarzschild time $t/\ell$", title="A finite domain of indistinguishable exteriors")
        ax.legend(loc="upper right", fontsize=8)
        fig.suptitle("From a sourced interior to a parent observer", fontsize=17)
        fig.supxlabel("Initial-data construction and local vacuum uniqueness; formation, global shell evolution and sibling signalling remain open.", fontsize=10)
        outputs = {}
        for fmt in ("svg", "png"):
            stream = BytesIO()
            metadata = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
            fig.savefig(stream, format=fmt, dpi=145, metadata=metadata)
            outputs[f"connection-slice.{fmt}"] = stream.getvalue()
        plt.close(fig)
    return outputs


def narrative(connection, observer):
    return f"""# A first interior connection, and what a parent observer can know

![The current map](map.svg)

The first sourced body now has one connected initial slice: a regular expanding
de Sitter cap, its fixed positive-tension wall, a Schwarzschild vacuum transition,
a regular minimum sphere, and an exact asymptotically Euclidean end. The wall
source H ell=1, k ell=2 and state M/ell=0.4 are unchanged from v3. The transition
changes the slicing of vacuum Schwarzschild, not its matter content.

Both bulk constraints and the wall's distributional energy/momentum constraints
pass. An independent metric-derivative calculation checks the vacuum residuals.
The cap is a flat spatial ball, not an assumed closed S3 cosmos. The spatial
bridge is not a demonstrated signal channel or a black-hole formation history.

![Intrinsic radial profile and protected observer domain](connection-slice.svg)

The many-body theorem can operate on the smooth vacuum end alone. Its unchanged
open overlap reattaches the sourced body without imposing a global smoothness
condition on the wall. Pairing this body with an unequal ordinary Schwarzschild
bridge gives a justified initial-data existence family: a common outer parent
and B's retained other end. The connector metric and separation threshold have
not been calculated. The old explicit Brill-Lindquist data remain a separate
candidate and are not replaced by a claimed spherical version of themselves.

A static observer at R/ell=1.6 is protected by the retained right vacuum initial
region 0.8<R/ell<4. The open bound is t/ell <
{observer['window']['coordinate_duration_over_ell']:.12f}, or local proper time
tau_obs/ell < {observer['window']['observer_proper_duration_over_ell']:.12f}.
Different far assemblies that preserve those initial data give the same vacuum
geometry inside this domain. This is an explicit local degeneracy, not a
predicted arrival or a cross-universe age. Independent field/state data can
still differ. The inward optical endpoint lies at minus infinity.

The next question is the sourced body's causal completion and stability, or an
observable capable of breaking this finite equivalence. Known completions and
obstructions should be imported before a new calculation is admitted. The
other branching, communicating-mouth and exterior-response routes remain open.

[Slice record](connection-slice.json), [parent assembly](parent-assembly.json),
[observer window](causal-window.json), [atlas](atlas.json), and
[manifest](manifest.json) retain inputs, scope and evidence. The
[working paper v4](../../../../paper/atlas-working-draft-v4.md) explains the
equations and the preserved alternatives. This closes a project interface in
one model; no new-to-world cosmological discovery is claimed.
"""


def payloads():
    record, connection, assembly, observer = build_record()
    svg = render_map(observer)
    outputs = {"atlas.json": first.atlas.canonical_bytes(record),
        "connection-slice.json": first.atlas.canonical_bytes(connection),
        "parent-assembly.json": first.atlas.canonical_bytes(assembly),
        "causal-window.json": first.atlas.canonical_bytes(observer),
        "map.svg": svg.encode(), "map.png": first.png_preview(svg),
        "narrative.md": narrative(connection, observer).encode()}
    outputs.update(scientific_figure(connection, observer))
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
        "outputs": {name: {"sha256": first.atlas.hashlib.sha256(blob).hexdigest(), "bytes": len(blob)} for name, blob in outputs.items()},
        "runtime": {"python": platform.python_version(), **{name: version(name) for name in ("numpy", "scipy", "sympy", "matplotlib")}},
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12},
        "scope": "explicit sourced initial slice, end-only parent existence and finite vacuum observer equivalence"}
    directory.mkdir(parents=True, exist_ok=False)
    for name, blob in outputs.items():
        with (directory / name).open("xb") as handle:
            handle.write(blob)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(first.atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "new_physical_mechanism": False}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA or manifest["source_bindings"] != source_bindings():
        raise ValueError("atlas source binding mismatch")
    expected = payloads()
    if set(expected) != set(manifest["outputs"]):
        raise ValueError("atlas output inventory mismatch")
    for name, replay in expected.items():
        path, binding = directory / name, manifest["outputs"][name]
        if first._sha(path) != binding["sha256"] or path.stat().st_size != binding["bytes"]:
            raise ValueError(f"artifact binding mismatch: {name}")
        equal = first._close(json.loads(path.read_bytes()), json.loads(replay)) if name.endswith(".json") else path.read_bytes() == replay
        if not equal:
            raise ValueError(f"atlas replay mismatch: {name}")
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
