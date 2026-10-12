#!/usr/bin/env python3
"""Record a finite exterior preparation operator and the atlas research dependencies."""
from __future__ import annotations

import argparse
from copy import deepcopy
from importlib.metadata import version
from io import BytesIO
import json
from pathlib import Path
import platform
import textwrap

import numpy as np

import derive_nsc_cosmological_atlas_v6 as previous
from recursive_horizons import nsc_atlas_field_preparation as preparation

first = previous.first
REPO = previous.REPO
SCHEMA = "nsc-cosmological-atlas-v7"
DEFAULT = REPO / "results/development/nsc-cosmological-atlas-v7"


def build_record():
    record, _ = previous.build_record()
    prepared = preparation.preparation_record()
    programme = first.atlas.programme_record()
    scenario = deepcopy(record["scenarios"][-1])
    scenario.update(id="finite_iv_field_preparation", title="Independent characteristic data in a retained IV wedge",
        solution_status="finite transfer operator certified; transmitted waveform not evaluated",
        preparation_record="preparation-access.json")
    scenario.pop("field_record", None)
    scenario["assumptions"] = [
        "isolated spherical reference with the published singular-past monotonic completion",
        "event-specific geometric access; independent proper clocks with explicit origins",
        "independent compact C3 characteristic data g(V) on U=U0 in retained region IV",
        "zero interior cap data at the initial wall radius R/ell=0.45 and compatible zero bottom characteristic",
        "finite wall-radius interval 0.45 <= R/ell <= 0.55 with a regular-centre no-return margin",
        "classical linear probe on the fixed isolated reference; no global parent source or waveform evaluation",
    ]
    for region in scenario["regions"]:
        if region["id"] == "C":
            region["state"] = "zero classical cap data at wall R/ell=0.45; subsequent response defined by the finite composed operator, not evaluated"
        elif region["id"] == "W":
            region["state"] = "trace uniquely determined by independent retained-IV characteristic data; trace and response not evaluated"
    scenario["links"] = [first.atlas._link(
        "scalar-W-C", "W", "C", "causal", ["scalar_matching", "scalar_state_transfer"],
        "Borrowed wall-to-centre transfer composed with the finite IV preparation operator; neither new wall trace nor observer waveform is evaluated",
        status="conditional_field_transfer_operator",
        missing=["evaluation of the composed operator for a separately admitted observable", "backreaction and detector dynamics"],
        delay=None, interval=None, clock=None)]
    scenario["regions"].append(first.atlas._region(
        "S_IV", "Finite IV preparation characteristic", "preparation_surface",
        {"chart": "retained Schwarzschild region IV in Kruskal U,V",
         "symmetry_orbits": "S2", "domain": deepcopy(prepared["geometry"]),
         "embedding": "finite characteristic boundary, not another child universe"},
        "specified initial-boundary data for the borrowed scalar operator",
        "independently prescribed compact C3 g(V); zero compatible bottom characteristic",
        {"kind": "Kruskal characteristic labels; not an observer proper clock"},
        {"field_amplitude": "epsilon is a test-field input"}))
    record["sources"].update({"finite_iv_preparation": {
        "path": "docs/nsc-atlas-field-preparation.md",
        "use": "finite retained-IV characteristic problem and bounded transfer operator",
        "limit": "independent boundary data, not an astrophysical preparation or solved waveform"}})
    scenario["links"].append(first.atlas._link(
        "preparation-IV-W", "S_IV", "W", "causal", ["finite_iv_preparation"],
        "A bounded finite characteristic operator delivers a wall trace; its waveform is not evaluated",
        status="conditional_preparation_operator",
        missing=["numerical wall trace and composed observer waveform", "physical production of the independent IV data"],
        delay=None, interval=None, clock=None))
    owner = ["docs/nsc-atlas-field-preparation.md"]
    evidence = ["results/development/nsc-cosmological-atlas-v7/preparation-access.json"]
    common = {"domain": "Fixed linear conformal scalar; finite retained IV wedge and early wall window R/ell in [0.45,0.55]",
        "assumptions": ["fixed geometry and scalar law", "independent compatible characteristic input", "no centre return within the admitted wedge"],
        "owner_refs": owner, "evidence_refs": evidence}
    scenario["connection_records"] = [
        {**deepcopy(common), "id": "finite-IV-preparation", "kind": "preparation", "status": "established_in_domain",
         "origin": "prescribed_at_boundary", "remote_preparation_established": False,
         "prescribed_boundary": "U=U0 outside the wall in retained IV; g(V) is not chosen from a target wall trace",
         "missing": ["astrophysical production of the characteristic data", "evaluated transmitted waveform"],
         "region_refs": ["S_IV", "W", "C"], "link_refs": ["preparation-IV-W", "scalar-W-C"],
         "prerequisite_ids": ["boundary-classification", "wall-transmission"], "unlock_ids": ["wall-trace", "observer-response"],
         "independent_check_refs": ["tests/test_nsc_atlas_field_preparation.py"]},
        {**deepcopy(common), "id": "open-physical-inheritance", "kind": "transition_inheritance", "status": "unresolved",
         "map_ref": None, "transferred_quantities": [], "physical_inheritance_established": False,
         "missing": ["a physical child-formation and dimensionless-law transition map"],
         "region_refs": ["C", "W"], "link_refs": [], "prerequisite_ids": ["boundary-classification"],
         "unlock_ids": ["regional-scales", "field-orientation", "branching-recurrence"]},
        {**deepcopy(common), "id": "open-observer-equivalence", "kind": "observer_equivalence", "status": "unresolved",
         "observer_region": "C", "observables": ["central scalar field within the finite arrival interval"],
         "alternatives": ["the recorded compact IV input", "a distinct compatible IV input or exterior model"],
         "equivalence_established": False,
         "discriminator": "Compare independent observer measurements after both traces are physically specified",
         "missing": ["a second fixed compatible candidate", "a computed response comparison or exact degeneracy"],
         "region_refs": ["C", "S_IV"], "link_refs": [], "prerequisite_ids": ["observer-response"],
         "unlock_ids": ["observer-inference"]},
    ]
    # Enrich the reusable programme with this finite result without declaring
    # the global preparation question closed or admitting a new simulation.
    for node in programme["calculation_graph"]:
        if node["id"] == "parent-preparation":
            node.update(status="finite_IV_operator_certified",
                missing=["physical generation of the independent IV data", "extension beyond this finite IV wedge"],
                decision="Reuse the finite transfer-operator proof; do not mistake it for preparation by the excluded right parent.")
            node["owner_refs"] = sorted(set(node["owner_refs"] + owner))
            node["evidence_refs"] = sorted(set(node["evidence_refs"] + evidence))
        elif node["id"] == "wall-trace":
            node.update(status="operator_available_waveform_not_evaluated",
                missing=["evaluation of the justified operator on the independent g(V), if admitted for a named new observable"],
                decision="Use the certified operator for a later admitted wall trace; no waveform is computed in this bundle.")
    record["programme"] = programme
    record["scenarios"].append(scenario)
    record["assessments"][scenario["id"]] = {
        "record": "preparation-access.json", "transfer_operator_certified": True,
        "independent_IV_data_prescribed": True, "right_parent_access_changed": False,
        "wall_waveform_evaluated": False, "observer_waveform_evaluated": False,
        "global_preparation_or_quantum_state_solved": False,
    }
    record["admission"].update(category="finite characteristic preparation and transfer-operator existence",
        missing_arrow="evaluate an admitted observable from the independent IV data, or supply a different physically compatible preparation")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    record["predecessor_manifest"] = "../nsc-cosmological-atlas-v6/manifest.json"
    return record, programme, prepared


def source_bindings():
    bindings = previous.source_bindings()
    paths = {Path(__file__), Path(preparation.__file__),
        REPO / "tests/test_nsc_atlas_field_preparation.py", REPO / "tests/test_nsc_atlas_programme.py",
        REPO / "tests/test_nsc_cosmological_atlas_v7.py", REPO / "docs/nsc-atlas-field-preparation.md",
        REPO / "docs/nsc-cosmological-atlas-v7.md"}
    # Owners are part of the reusable inventory. The living paper is a consumer
    # and is deliberately excluded, preventing a figure/manuscript hash cycle.
    for question in first.atlas.programme_record()["questions"]:
        paths.update(REPO / ref.split("#", 1)[0] for ref in question["owner_refs"] + question["evidence_refs"])
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def _figure_bytes(figure, stem):
    result = {}
    for fmt in ("svg", "png"):
        stream = BytesIO()
        metadata = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
        figure.savefig(stream, format=fmt, dpi=170, metadata=metadata, facecolor="white")
        result[f"{stem}.{fmt}"] = stream.getvalue()
    return result


def preparation_figure(record):
    """Plot actual finite-domain coordinates and independent input; no response."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    geometry, boundary = record["geometry"], record["boundary_samples"]
    wall, source, bottom = boundary["wall"], boundary["source"], boundary["bottom"]
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 14}):
        fig, (ax, input_ax) = plt.subplots(1, 2, figsize=(13, 7.8), layout="constrained", gridspec_kw={"width_ratios": [1.35, 1]})
        ax.fill([geometry["U0"], *wall["U"], geometry["U0"]],
                [geometry["Va"], *wall["V"], geometry["Vb"]], color="#dcecf2", label="Retained IV domain")
        ax.plot(source["U"], source["V"], color="#b26825", linewidth=3, label=r"Initial data: $U=U_0$")
        ax.plot(wall["U"], wall["V"], color="#246485", linewidth=2.5, label=r"Matched wall: $U_W(V)$")
        ax.plot(bottom["U"], bottom["V"], color="#555555", linestyle="--", linewidth=1.8, label=r"Bottom: $\chi=0$")
        ax.plot([geometry["U0"], geometry["Ub"]], [geometry["Vb"]] * 2, color="#8399a5", linestyle=":", linewidth=1.2)
        ax.annotate("Future null\ndirections", xy=(.25, .77), xycoords="axes fraction", fontsize=13, color="#183b4b")
        for end in ((.33, .69), (.14, .9)):
            ax.annotate("", xy=end, xytext=(.14, .69), xycoords="axes fraction",
                        arrowprops={"arrowstyle": "-|>", "color": "#183b4b", "lw": 1.3})
        ax.text(.12, .60, r"Retained side: $U\leq U_W(V)$", transform=ax.transAxes,
                fontsize=13, color="#183b4b")
        ax.set(xlabel="Kruskal U\n(characteristic label)", ylabel="Kruskal V (characteristic label)",
               title="Finite characteristic domain\noutside the wall")
        ax.ticklabel_format(style="plain", useOffset=False)
        ax.legend(fontsize=12, loc="lower right", framealpha=.94)
        ax.grid(alpha=.17)
        input_ax.plot(source["V"], source["chi_over_epsilon"], color="#b26825", linewidth=2.4)
        input_ax.fill_between(source["V"], source["chi_over_epsilon"], color="#ead4ba", alpha=.6)
        input_ax.axhline(0, color="#777777", linewidth=.8)
        input_ax.set(xlabel="Kruskal V on the independent\ninitial characteristic", ylabel=r"Prescribed $g(V)/\varepsilon$",
                     title="Independent input\nnot selected from a wall trace")
        input_ax.grid(alpha=.17)
        q = record["operator"]["analytic_norm_bound"]
        input_ax.text(.04, .05, f"Operator bound: q* = {q:.8f} < 1\nNo transmitted or observed\nwaveform evaluated",
                      transform=input_ax.transAxes, fontsize=12, color="#183b4b",
                      bbox={"facecolor": "white", "edgecolor": "#cad8df", "alpha": .95})
        fig.suptitle("Independent exterior data reach the early wall\nthrough a certified finite operator", fontsize=18)
        fig.supxlabel("Region IV data are prescribed; right-parent signalling is excluded.\nRemote source and quantum state remain unspecified.", fontsize=12)
        result = _figure_bytes(fig, "preparation-domain")
        plt.close(fig)
    return result


def dependency_figure(programme):
    """Research prerequisites are not physical arrows between universes."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
    from matplotlib.path import Path as FigurePath

    positions = {"geometry-baseline": (1.2, 6), "boundary-classification": (1.2, 5),
        "wall-transmission": (-.6, 4), "parent-preparation": (1.2, 4), "wall-trace": (1.2, 3),
        "observer-response": (1.2, 2), "observer-inference": (1.2, 1), "inheritance-map": (3, 4),
        "regional-scales": (4.6, 3), "field-orientation": (3, 3), "clock-entropy": (3, 2), "branching-recurrence": (3, 1)}
    labels = {"geometry-baseline": "Borrowed region geometry", "boundary-classification": "Boundary and source compatibility",
        "wall-transmission": "Borrowed wall transmission", "parent-preparation": "Independent IV preparation",
        "wall-trace": "Prepared wall trace", "observer-response": "Interior measurement",
        "observer-inference": "Distinguish or equate exteriors", "inheritance-map": "Physical inheritance map",
        "regional-scales": "Dimensionless regional laws", "field-orientation": "Charge and field orientation",
        "clock-entropy": "Proper clocks and entropy", "branching-recurrence": "Branching and recycling"}
    colors = {"reusable_result": "#deebd6", "known_equations_missing_data": "#dcecf6",
              "known_equations_missing_compatibility": "#e8e0f0", "missing_physical_mechanism": "#f3e3d5"}
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 10}):
        fig, ax = plt.subplots(figsize=(13, 10.8))
        fig.subplots_adjust(left=.02, right=.98, bottom=.23, top=.93)
        boxes = {}
        for row in programme["calculation_graph"]:
            x, y = positions[row["id"]]
            box = FancyBboxPatch((x-.66, y-.35), 1.32, .70, boxstyle="round,pad=.03",
                facecolor=colors[row["gap_kind"]], edgecolor="#82949f", linewidth=.9, zorder=2)
            boxes[row["id"]] = box
            ax.add_patch(box)
            ax.text(x, y, textwrap.fill(labels[row["id"]], 18), ha="center", va="center", fontsize=16, color="#183441", zorder=3)
        for row in programme["calculation_graph"]:
            for previous_id in row["prerequisite_ids"]:
                settings = {"arrowstyle": "-|>", "mutation_scale": 13, "color": "#657581", "linewidth": 1.3, "zorder": 1}
                if previous_id == "inheritance-map" and row["id"] in {"clock-entropy", "branching-recurrence"}:
                    # Separate long prerequisites from the intermediate nodes;
                    # crossing a node would imply an unrecorded dependency.
                    y = positions[row["id"]][1]
                    route_x = 3.85 if row["id"] == "clock-entropy" else 5.5
                    route = FigurePath([(3.72, 4), (route_x, 4), (route_x, y), (3.72, y)])
                    arrow = FancyArrowPatch(path=route, **settings)
                else:
                    arrow = FancyArrowPatch(positions[previous_id], positions[row["id"]],
                        patchA=boxes[previous_id], patchB=boxes[row["id"]], shrinkA=3, shrinkB=3, **settings)
                ax.add_patch(arrow)
        ax.set(xlim=(-1.4, 5.7), ylim=(.4, 6.6))
        ax.axis("off")
        legend = [Patch(facecolor=colors[key], edgecolor="#82949f", label=textwrap.fill(label, 35)) for key, label in (
            ("reusable_result", "Reuse a known result"), ("known_equations_missing_data", "Known equations; data/trace gap"),
            ("known_equations_missing_compatibility", "Compatibility or identifiability gap"),
            ("missing_physical_mechanism", "Missing physical mechanism"))]
        fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=13, frameon=False, bbox_to_anchor=(.5, .075))
        fig.suptitle("Which missing connection unlocks the next question?", fontsize=18)
        fig.supxlabel("Research prerequisites, not causal paths. Candidate routes remain separate.\nEach new simulation requires its own calculation admission.", fontsize=12, y=.018)
        result = _figure_bytes(fig, "research-dependencies")
        plt.close(fig)
    return result


def narrative(prepared):
    q = prepared["operator"]["analytic_norm_bound"]
    return f"""# An independent preparation and a finite transfer operator

![Finite retained-IV preparation and independent input](preparation-domain.svg)

The new entry specifies compact C3 data g(V) on a finite retained-IV
characteristic U=U0, with a compatible zero bottom characteristic. The datum
is chosen independently of a desired wall trace. This is a linear classical
initial-boundary problem on the existing isolated geometry; the right parent
observer's causal exclusion is unchanged.

The interior cap is initialized with zero data at wall radius R/ell=0.45.
This is a separate finite preparation from the v6 pulse and its earlier cap.
Connecting both preparations requires an additional no-earlier-incoming-history
assumption; this bundle does not assert that connection.

Borrowed characteristic equations and the matched wall give a finite transfer
operator with analytic contraction bound **q*={q:.12f}<1**. This establishes
existence and uniqueness in the declared finite domain. It does not evaluate
a new transmitted waveform, prepare an astrophysical source, select a quantum
vacuum or establish cosmological recurrence. The previously computed v6 pulse
remains its own conditional wall-input result, with all earlier routes retained.

| Observer or surface | Access in this record | Missing connection or limit |
|---|---|---|
| Independent IV characteristic | Supplies g(V) to the retained finite wedge | Physical production of the prescribed data is not modelled |
| Wall W | Receives a uniquely defined trace through the certified operator | The transmitted waveform has not been evaluated |
| Interior centre C | Existing early wall-to-centre transfer can be composed in its clock domain | No new central waveform or detector response is claimed |
| Right static parent P | No P-to-W or P-to-C signal in this isolated reference | Spatial connection does not reopen the excluded causal route |
| Other children or exterior models | Their original candidate branches remain available | No new sibling signalling or physical inheritance map is added |

The recorded wall interval is R/ell in [0.45,0.55]. The lower radius bound,
finite proper-time arrival interval and regular-centre no-return margin belong
to [preparation-access.json](preparation-access.json). Outside that domain the
argument must be re-established. A characteristic label is not an observer's age.

![Research prerequisite graph](research-dependencies.svg)

The eight questions are spheres, inside/outside, quantum fields, barriers,
regional scales, charge, time and branching. The diagram shows research
prerequisites; it does not add physical connections. The useful next calculation
must close a named gap or distinguish alternatives. Recomputing a known identity
without new data, compatibility or an observable does not advance this map.

[programme.json](programme.json) records those dependencies, reading boundaries
and stopping rules. [atlas.json](atlas.json) retains all nine prior scenarios
and appends this finite preparation. [manifest.json](manifest.json) authenticates
the source and output bytes. The scientific owner is
[field preparation](../../../docs/nsc-atlas-field-preparation.md); the one living
[paper](../../../paper/nested-space-cosmology.md) remains the manuscript.
"""


def payloads():
    atlas, programme, prepared = build_record()
    outputs = {"atlas.json": first.atlas.canonical_bytes(atlas),
        "programme.json": first.atlas.canonical_bytes(programme),
        "preparation-access.json": first.atlas.canonical_bytes(prepared),
        "narrative.md": narrative(prepared).encode("utf-8")}
    outputs.update(preparation_figure(prepared))
    outputs.update(dependency_figure(programme))
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
        "scope": "finite retained-IV preparation and exact transfer-operator existence; no evaluated new waveform"}
    directory.mkdir(parents=True, exist_ok=False)
    for name, data in outputs.items():
        with (directory / name).open("xb") as handle:
            handle.write(data)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(first.atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "waveform_evaluated": False}


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
