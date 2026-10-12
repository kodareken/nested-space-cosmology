#!/usr/bin/env python3
"""Create or replay a finite atlas of moving walls and retained spherical collars."""
from __future__ import annotations

import argparse
from copy import deepcopy
from io import BytesIO
import json
from pathlib import Path
import platform
from importlib.metadata import version

import derive_nsc_cosmological_atlas_v2 as previous
from recursive_horizons import nsc_atlas_moving_boundaries as moving

first = previous.first
LAB, REPO = previous.LAB, previous.REPO
SCHEMA = "nsc-cosmological-atlas-v3"
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v3"


def parent_gluing_record():
    """A literature contract, not a numerically constructed gluing metric."""
    return {
        "schema": "nsc-atlas-parent-gluing-v1",
        "status": "IMPORTED_EXISTENCE_FAMILY",
        "source": {"url": "https://arxiv.org/html/0909.1101v1",
                   "read": "Theorem 1.1, ensuing interior-preservation discussion and proof obstruction balance",
                   "owner": "docs/nsc-atlas-parent-gluing.md"},
        "number_of_bodies": 2,
        "body_masses": {"m_A": "positive", "m_B": "positive, need not equal m_A"},
        "preserved": "chosen compact Schwarzschild body regions, including selected throat collars",
        "constraints": "K_ij=0 and R[h]=0 on the vacuum construction",
        "parent": "one shared outer asymptotically Schwarzschild end; body inner ends may remain",
        "separation": {"condition": "sufficiently large, depending on the chosen bodies",
                       "numeric_threshold": None, "BL_d12_certified": False},
        "outer_mass": "arbitrarily close to m_A+m_B; not stipulated exactly equal",
        "construction": "deform distant transition regions and balance mass/centre obstruction moments",
        "globally_spherical": False,
        "explicit_metric_constructed_here": False,
        "unchanged_brill_lindquist": False,
        "moving_wall_attached": False,
        "formation_established": False,
        "non_time_symmetric_successor": {
            "url": "https://arxiv.org/html/1004.1355v1",
            "read": "Theorem 1 and Corollary 1.1, introduction and regularity assumptions",
            "use": "preserve selected smooth (h,K) bodies with vacuum AE ends and timelike ADM four-momenta",
            "missing": "an ideal delta-layer is not automatically in the theorem's smooth-data class"},
        "missing_conditions": [
            "select a spacelike Cauchy slice through the moving wall and retained Schwarzschild side",
            "match induced spatial metric, extrinsic curvature and source constraints on overlap",
            "supply smooth source data or justify a gluing theorem for the distributional shell",
            "specify retained collar extent and its finite domain of dependence",
            "verify theorem hypotheses before choosing a numerical separation"],
        "admission": {"kind": "imported theorem, no new elliptic solve",
            "borrowed": "0909.1101v1 Theorem 1.1 and 1004.1355v1 Theorem 1 / Corollary 1.1",
            "missing": "wall-body initial data satisfying the gluing hypotheses",
            "decides": "exact local spherical geometry is available in a separate many-body family",
            "unlocks": "a Cauchy-interface check, not a new gluing-theorem derivation",
            "finite_domain": "two bodies with specified compact retained regions",
            "budget": "literature import only; no PDE solve or inferred numerical threshold",
            "stop": "leave the wall attachment unresolved until its full data are supplied"},
    }


def build_record():
    record, p, j = previous.build_record()
    m, g = moving.moving_boundary_record(), parent_gluing_record()
    record["sources"].update({
        "moving_wall": {"path": "docs/nsc-atlas-moving-boundaries.md",
                        "use": "signed local matching and fixed-source conservation screen"},
        "many_body_gluing": deepcopy(g["source"]),
    })
    region, link = first.atlas._region, first.atlas._link
    wall = m["vacuum_wall_control"]
    wall_regions = []
    for key, title, metric, source in (
        ("schwarzschild", "Selected Schwarzschild side", "f_plus=1-2M/R", "vacuum"),
        ("de_sitter", "Selected de Sitter side", "f_minus=1-H^2 R^2", "positive vacuum energy")):
        wall_regions.append(region(key, title, "local_patch",
            {"chart": "signed spherical shell patch", "metric": metric,
             "embedding": "no global embedding supplied", "symmetry_orbits": "S2",
             "wall_radius_domain_over_ell": deepcopy(wall["domain"]["R_over_ell"])},
            source, "classical source; no quantum state supplied",
            {"kind": "shell proper time on the recorded finite segment only",
             "duration_over_ell": wall["domain"]["proper_duration_over_ell"],
             "static_coordinate": wall["clock"]["de_sitter_static_t" if key == "de_sitter" else "schwarzschild_patch_t"]},
            {"H_ell": wall["theory"]["H_ell"]} if key == "de_sitter" else
            {"M_over_ell": wall["theory"]["M_over_ell"]}))
    record["scenarios"].append({
        "id": "moving_vacuum_wall", "title": "A separate moving vacuum wall",
        "theory": {"family": "de Sitter / Schwarzschild with fixed positive pure tension",
                   "couplings": {"G": 1, "c": 1, "H_times_ell": 1, "k_times_ell": 2}},
        "assumptions": ["M/ell=0.4", "expanding finite shell segment",
                        "signed beta_minus positive and beta_plus negative",
                        "no parent embedding, formation or stability asserted"],
        "regions": wall_regions,
        "links": [link("vacuum-wall", "schwarzschild", "de_sitter", "attachment",
            ["moving_wall"], "Unsquared angular and temporal Israel equations with fixed pure tension",
            status="imported_local_dynamics",
            missing=["global Cauchy slice", "shared parent embedding", "formation history"])],
        "max_descent_depth": 0,
        "solution_status": "imported finite local junction, not established offspring",
    })
    parent_regions = [region("glued_parent", "Many-body parent family", "parent",
        {"chart": "existence family, no explicit metric supplied here", "globally_spherical": False},
        "vacuum constraints in the gluing region", "classical initial data",
        {"kind": "initial slice only"}, {})]
    parent_links = []
    for key in ("a", "b"):
        parent_regions.append(region(f"collar_{key}", f"Retained spherical collar {key.upper()}", "interface",
            {"chart": "exact Schwarzschild spatial region preserved by the theorem",
             "symmetry_orbits": "S2", "explicit_embedding_supplied": False},
            "vacuum body data", "time-symmetric (h,K=0)", {"kind": "no shared observer clock"},
            {"mass": f"m_{key.upper()}>0; unequal allowed"}, contained_in="glued_parent"))
        parent_links.append(link(f"retained-{key}", "glued_parent", f"collar_{key}", "attachment",
            ["many_body_gluing"], "Theorem preserves chosen compact body regions in a common outer parent",
            status="imported_existence_family",
            missing=["explicit construction and separation domain", "moving wall Cauchy data"]))
    record["scenarios"].append({
        "id": "many_body_spherical_collars", "title": "A parent family with exact spherical collars",
        "theory": {"family": "time-symmetric vacuum Einstein constraints", "couplings": {"G": 1, "c": 1}},
        "assumptions": ["two smooth positive-mass AE Schwarzschild bodies", "sufficiently large separation",
                        "transition geometry may change; original BL metric not preserved"],
        "regions": parent_regions, "links": parent_links, "max_descent_depth": 0,
        "solution_status": "borrowed existence family; moving-wall attachment unresolved",
    })
    record["assessments"].update({
        "moving_boundary": {"record": "moving-boundary.json", "new_physics_result": False,
                            "embedded_in_shared_parent": False},
        "spherical_collar_parent": {"record": "parent-gluing.json", "status": g["status"],
                                    "wall_attached": False},
    })
    record["admission"].update(category="source conservation and imported parent-family compatibility",
        missing_arrow="complete wall-body Cauchy data and a justified source regularity class")
    record["schema"] = first.atlas.SCHEMA
    first.atlas.validate_atlas(record)
    record["schema"], record["predecessor"] = SCHEMA, previous.SCHEMA
    return record, m, g


def source_bindings():
    bindings = previous.source_bindings()
    paths = [Path(__file__), Path(moving.__file__),
        LAB / "tests/test_nsc_atlas_moving_boundaries.py",
        LAB / "tests/test_nsc_cosmological_atlas_v3.py",
        LAB / "docs/nsc-atlas-moving-boundaries.md",
        LAB / "docs/nsc-atlas-parent-gluing.md",
        LAB / "docs/nsc-cosmological-atlas-v3.md",
        REPO / "paper/atlas-working-draft-v3.md"]
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def render_map():
    text, panel = first._text, first._panel
    out = first._base(1400, 1460, "Which boundaries can become part of the map?",
        "Atlas v3 | Two parent families, a source constraint and one moving local junction")
    for x, title, lines in (
        (50, "EXPLICIT SHARED PARENT", ["Brill-Lindquist: two unequal punctures",
            "One vacuum slice, with angular dependence.", "Local spheres are generally distorted.",
            "The v2 metric and evidence remain unchanged."]),
        (725, "PARENT FROM A PUBLISHED THEOREM", ["Two retained Schwarzschild body regions",
            "A many-body theorem supplies a shared parent.", "Exact spherical collars; sufficiently far apart.",
            "No explicit separation threshold or new metric."])):
        out += [panel(x, 135, 625, 300, "#102c41", "#47758e"), text(x+25, 177, title, 21, "#86d6f0", "600")]
        for i, line in enumerate(lines):
            out.append(text(x+25, 218+34*i, line, 18 if i else 20, "#e5eef3"))
        for cx, rad in ((x+215, 26), (x+380, 35)):
            out.append(f'<circle cx="{cx}" cy="393" r="{rad}" fill="#020810" stroke="#e6be7a" stroke-width="2"/>')
    out += [panel(50, 475, 1300, 152, "#302d28", "#c69f62"),
        text(80, 516, "MISSING INTERFACE: PUT THE WALL AND PARENT ON ONE SPACELIKE SLICE", 21, "#f4ce8c", "600"),
        text(80, 553, "Match the spatial metric, extrinsic curvature, retained side and source constraints.", 20),
        text(80, 590, "Smooth-data gluing does not automatically include an ideal delta-layer. A shared mass is insufficient.", 18, "#c7c9c9")]
    for cx in (360, 1040):
        out.append(f'<path d="M {cx} 437 L {cx} 472" stroke="#e6be7a" stroke-width="2" stroke-dasharray="6 6"/>')
    out += [panel(50, 669, 625, 270, "#382633", "#bf7881"),
        text(75, 711, "RADIATION + PURE FIXED TENSION", 21, "#f0b1b8", "600"),
        text(75, 751, "Vacuum outside: no compensating bulk flux.", 19),
        text(75, 789, "Conservation forces zero relative fluid speed.", 19),
        text(75, 827, "Motion alone does not evade the v2 obstruction.", 18),
        text(75, 869, "Change the source or bulk assumptions to proceed.", 17, "#e0b6bf"),
        text(75, 907, "Scoped exclusion; other child models remain open.", 17, "#c7c9c9"),
        panel(725, 669, 625, 270, "#153630", "#61aa91"),
        text(750, 711, "DE SITTER + PURE FIXED TENSION", 21, "#8fe0bc", "600"),
        text(750, 751, "Imported expanding shell with a fixed source.", 19),
        text(750, 789, "Signed normals: beta-minus > 0, beta-plus < 0.", 18),
        text(750, 827, "Local proper time and both junction stresses checked.", 17),
        text(750, 869, "A separate spherical patch; parent embedding open.", 17, "#a9cfc1"),
        text(750, 907, "No formation history or stability result assigned.", 17, "#c7c9c9"),
        text(60, 987, "THE FINITE BRANCHING QUESTION REMAINS", 22, "#86d6f0", "600")]
    for x, label in ((80, "Possible child A"), (880, "Possible child B")):
        out += [panel(x, 1015, 425, 95), text(x+28, 1055, label, 23, "#e5eef3", "600"),
                text(x+28, 1090, "Formation and state transfer unresolved", 17, "#e6be7a")]
    out += ['<path d="M 505 1060 L 877 1060" stroke="#ab8bda" stroke-width="2" stroke-dasharray="3 8"/>',
        text(560, 1040, "possible correlation", 17, "#c9aaf0"),
        text(550, 1095, "no sibling signal path", 17, "#adc0cd"),
        '<path d="M 290 1111 L 290 1160" stroke="#e6be7a" stroke-width="2" stroke-dasharray="6 6" marker-end="url(#arrow)"/>',
        panel(80, 1170, 425, 83), text(105, 1205, "Possible grandchild", 22, "#86d6f0", "600"),
        text(105, 1236, "Another transition; depth cutoff two", 17),
        text(570, 1175, "One-hole multiple daughters", 20, "#86d6f0"),
        text(570, 1210, "Communicating mouths / indistinguishable exteriors", 19),
        text(570, 1245, "Independent candidates retain their own sources and domains.", 17, "#adc0cd"),
        text(60, 1320, "Descent / attachment / causal influence / correlation remain distinct relations.", 21, "#86d6f0"),
        text(60, 1360, "Known mathematics changes which connection to test next. It is not a new cosmological discovery.", 18, "#adc0cd"),
        text(60, 1400, "Schematic candidate map. Sizes are not distances; no global spacetime embedding is drawn.", 18, "#adc0cd"), "</svg>"]
    return "\n".join(out) + "\n"


def scientific_figure(m):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt
    samples = m["vacuum_wall_control"]["samples"]
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans", "font.size": 11}):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), layout="constrained")
        axes[0].plot(samples["tau_over_ell"], samples["R_over_ell"], color="#175778", linewidth=2.5)
        axes[0].set(xlabel=r"Shell proper time $\tau/\ell$", ylabel=r"Areal radius $R/\ell$",
                    title="A finite expanding segment")
        radii = samples["R_over_ell"]
        for key, color, label in (("beta_minus", "#257256", r"$\beta_-$: de Sitter side"),
                                 ("beta_plus", "#ac4d37", r"$\beta_+$: Schwarzschild side")):
            axes[1].plot(radii, samples[key], color=color, label=label, linewidth=2.5)
        axes[1].axhline(0, color="#777777", linewidth=.8)
        axes[1].set(xlabel=r"Areal radius $R/\ell$", ylabel="Signed normal derivative of radius",
                    title="Orientation survives squaring")
        axes[1].legend(fontsize=10)
        for ax in axes:
            ax.grid(alpha=.2)
        fig.suptitle(r"Published vacuum-wall dynamics: $H\ell=1$, $4\pi\sigma\ell=2$, $M/\ell=0.4$", fontsize=16)
        fig.supxlabel("Imported local junction. No parent embedding, formation history, common cosmological clock or stability result.", fontsize=10)
        outputs = {}
        for fmt in ("svg", "png"):
            stream = BytesIO()
            metadata = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
            fig.savefig(stream, format=fmt, dpi=145, metadata=metadata)
            outputs[f"moving-boundary.{fmt}"] = stream.getvalue()
        plt.close(fig)
    return outputs


def narrative():
    return """# A boundary needs a source; a shared parent needs compatible data

![Successor candidate map](map.svg)

The explicit Brill-Lindquist parent remains our sampled baseline. A published
many-body gluing theorem supplies a different family with two unequal, exactly
spherical Schwarzschild body regions at sufficiently large separation. The
theorem is imported; its unknown separation bound is not replaced by the
baseline's coordinate distance.

For the radiation-filled candidate, a fixed pure-tension wall has no intrinsic
stress divergence. With vacuum outside, conservation forces its speed relative
to the radiation to vanish. Motion alone therefore does not evade the previous
comoving-wall obstruction. This is a scoped application of known conservation.

![Finite signed moving-wall control](moving-boundary.svg)

A de Sitter interior has rho+p=0 and avoids that particular flux obstruction.
Its published constant-tension junction supplies this finite expanding segment.
The negative Schwarzschild-side normal derivative is retained: replacing it by
a positive square root would describe the wrong attachment. Proper time belongs
to this wall only. Local matching does not prove offspring formation or stability.

The next interface is complete initial data: a spacelike slice through the wall,
the correct Schwarzschild side, matching spatial metric and extrinsic curvature,
and a source regularity class covered by the chosen gluing result. The smooth
non-time-symmetric many-body theorem is promising, but does not automatically
include a distributional Israel layer. We have not attached the wall to either
shared parent.

[atlas.json](atlas.json), [moving-boundary.json](moving-boundary.json),
[parent-gluing.json](parent-gluing.json) and [manifest.json](manifest.json)
record the scenarios, borrowed mathematics, finite checks and provenance.
The [working manuscript v3](../../../../paper/atlas-working-draft-v3.md)
develops the mathematical narrative. The three independent v1 alternatives,
unequal potential children and the finite grandchild remain visible; no new
causal path or local dark-sector consequence is assigned.
"""


def payloads():
    record, m, g = build_record()
    svg = render_map()
    outputs = {"atlas.json": first.atlas.canonical_bytes(record),
        "moving-boundary.json": first.atlas.canonical_bytes(m),
        "parent-gluing.json": first.atlas.canonical_bytes(g),
        "map.svg": svg.encode(), "map.png": first.png_preview(svg),
        "narrative.md": narrative().encode()}
    outputs.update(scientific_figure(m))
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
        "outputs": {name: {"sha256": first.atlas.hashlib.sha256(blob).hexdigest(), "bytes": len(blob)}
                    for name, blob in outputs.items()},
        "runtime": {"python": platform.python_version(),
                    **{name: version(name) for name in ("numpy", "scipy", "sympy", "matplotlib")}},
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12},
        "scope": "imported local dynamics and parent existence; shared wall-parent spacetime unresolved"}
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
