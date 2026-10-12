#!/usr/bin/env python3
"""Create a bounded atlas once, or verify its portable source/output bindings."""
from __future__ import annotations

import argparse
import ast
from html import escape
from io import BytesIO
from importlib.metadata import version
import json
from pathlib import Path
import platform
import re
import xml.etree.ElementTree as ET

from recursive_horizons import nsc_cosmological_atlas as atlas

LAB = Path(__file__).resolve().parents[1]
REPO = LAB
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v1"


def _sha(path):
    return atlas.hashlib.sha256(path.read_bytes()).hexdigest()


def source_bindings():
    """Hash transitive local Python imports plus owned test and explanation."""
    package = LAB / "src/recursive_horizons"
    pending = [package / "nsc_cosmological_atlas.py"]
    found = {Path(__file__).resolve(), LAB / "tests/test_nsc_cosmological_atlas.py",
             LAB / "docs/nsc-cosmological-atlas.md"}
    cited_documents = {REPO / source["path"] for source in atlas.SOURCES.values() if "path" in source}
    while pending:
        path = pending.pop()
        if path in found:
            continue
        found.add(path)
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.ImportFrom) and node.level == 1 and node.module:
                dependency = package / (node.module.replace(".", "/") + ".py")
                if dependency.is_file():
                    pending.append(dependency)
    found.update(cited_documents)
    return {p.relative_to(REPO).as_posix(): _sha(p) for p in sorted(found)}


def _text(x, y, content, size=18, color="#dce8f3", weight="normal"):
    return (f'<text x="{x}" y="{y}" font-family="Segoe UI,Arial,sans-serif" '
            f'font-size="{size}" fill="{color}" font-weight="{weight}">{escape(content)}</text>')


def _base(width, height, title, subtitle):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-label="{escape(title)}">',
            f'<title>{escape(title)}</title>',
            '<defs><radialGradient id="sea"><stop stop-color="#123a54"/><stop offset="1" stop-color="#071424"/></radialGradient>'
            '<marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
            '<path d="M 0 0 L 10 5 L 0 10 z" fill="#86d6f0"/></marker></defs>',
            f'<rect width="{width}" height="{height}" fill="#07111e"/>',
            _text(55, 64, title, 34, weight="600"), _text(55, 100, subtitle, 17, "#9fb7c9")]


def _panel(x, y, w, h, color="#18394d", stroke="#36576b"):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="18" fill="{color}" stroke="{stroke}"/>'


def render_default(atlas_record):
    """An explanatory patch diagram: positions/radii are not a spatial embedding."""
    scenario = atlas_record["scenarios"][0]
    regions = {r["id"]: r for r in scenario["regions"]}
    out = _base(1400, 1160, "A cosmos seen from its possible exterior",
                "Finite candidate atlas | Two distinct black-hole interfaces | No completed global solution")
    out += [_panel(50, 130, 1300, 410, "url(#sea)", "#47758e"),
            _text(85, 177, "PARENT REGION", 20, "#86d6f0", "600"),
            _text(85, 210, "Its full geometry and state remain unknown.", 22),
            _text(85, 243, "Separate local spherical charts can sit off centre in a nonspherical parent.", 18, "#adc0cd")]
    for key, x, r in (("a", 410, 56), ("b", 990, 70)):
        mass = regions[f"hole_{key}"]["dimensionless_parameters"]["mu_over_reference_length"]
        out += [f'<circle cx="{x}" cy="355" r="{r+12}" fill="none" stroke="#d9ae67" stroke-width="2" opacity=".7"/>',
                f'<circle cx="{x}" cy="355" r="{r}" fill="#020810" stroke="#f4ce8c" stroke-width="3"/>',
                _text(x-35, 361, f"BH {key.upper()}", 23, "#f4ce8c", "600"),
                _text(x-115, 459, f"Local chart: f = 1 - 2mu/R", 17),
                _text(x-115, 487, f"Illustrative mu/L = {mass:g}; R > 2mu", 16, "#adc0cd")]
    out += [_text(495, 296, "No exact two-hole parent metric has been supplied.", 17, "#adc0cd")]
    for key, x, center in (("a", 60, 410), ("b", 750, 990)):
        out += [f'<path d="M {center} 540 L {center} 599" stroke="#86d6f0" stroke-width="2.5" '
                'stroke-dasharray="8 7" marker-end="url(#arrow)"/>',
                _text(center+15, 575, "descent + unresolved neck", 15, "#86d6f0"),
                _panel(x, 605, 590, 230),
                _text(x+28, 645, f"POTENTIAL CHILD {key.upper()}", 21, "#86d6f0", "600"),
                _text(x+28, 684, "Closed FLRW is one candidate chart, not a result.", 17),
                _text(x+28, 718, "S3 describes a possible spatial slice; S2 a symmetry orbit.", 16, "#adc0cd"),
                _text(x+28, 752, "Its source, state and clock can differ from its sibling.", 17),
                _text(x+28, 792, "Needed: matching + physical layer + transition state", 16, "#e6be7a")]
    out += ['<path d="M 650 740 C 687 695 714 695 750 740" fill="none" stroke="#b497dc" '
            'stroke-width="2" stroke-dasharray="3 6"/>',
            _text(558, 870, "Shared history may correlate children; no sibling signal route is supplied.", 18, "#cbb6e7"),
            '<path d="M 280 835 L 280 920" stroke="#86d6f0" stroke-width="2" stroke-dasharray="8 7" marker-end="url(#arrow)"/>',
            _panel(60, 925, 430, 100, "#172838", "#466178"),
            _text(83, 961, "POTENTIAL GRANDCHILD", 18, "#86d6f0", "600"),
            _text(83, 993, "Formation interface missing. Depth stops at 2.", 15),
            _text(558, 925, "FOUR RELATIONS, FOUR DIFFERENT QUESTIONS", 18, "#86d6f0", "600"),
            _text(558, 963, "Descent: origin   /   Attachment: geometric matching", 18),
            _text(558, 999, "Causal: allowed influence   /   Correlation: shared state", 18),
            _text(55, 1085, "Circles and panel sizes are schematic. This is a patch diagram, not an embedding of one universe inside a 3D ball.", 17, "#adc0cd"),
            _text(55, 1118, "No dark-sector fractions are fitted. The next value is a justified connection from an exterior to a local observer.", 17, "#adc0cd"), '</svg>']
    return "\n".join(out) + "\n"


def render_candidates():
    out = _base(1400, 910, "Three different ways the exterior could matter",
                "Compare mathematical possibilities before committing to a source theory")
    cards = [
        ("01", "One hole, several daughters", "Null-shell construction", "hep-th/9511136v1",
         ["R0^3 = 2 m a^2", "f_A f_B = f_C f_D at a collision", "l > a: disconnected bubble boundaries"],
         ["Near singularity; equal R0; a << R0 << 2m.", "Borrowed kinematics, fixed limiting curvature.",
          "Nucleation and layer microphysics remain open."]),
        ("02", "Three connected mouths", "Conditional causal network", "2012.07821v2",
         ["E_small < |E_binding|", "G m_small < r_e", "m_small < N_f Q / (8 r_e)"],
         ["Magnetic, near-extremal, many-species domain.", "Not the one-field neutral-child construction."]),
        ("03", "Different, yet locally identical", "Response-equivalence control", "Known finite Schur algebra",
         ["Visible exterior spectrum: {2}", "Hidden spectrum: {3} or {5}", "Same local resolvent; hidden link = 0"],
         ["A matrix control, not a spacetime solution.", "Local data need not identify an entire exterior."]),
    ]
    for index, (number, title, label, source, formulas, limits) in enumerate(cards):
        x = 50 + index * 450
        out += [_panel(x, 140, 415, 665), _text(x+25, 190, number, 30, "#d9ae67", "600"),
                _text(x+25, 231, title, 22, weight="600"), _text(x+25, 263, label, 17, "#86d6f0")]
        if index == 0:
            out += [f'<circle cx="{x+207}" cy="355" r="67" fill="#061320" stroke="#e6be7a" stroke-width="2"/>']
            for offset in (-30, 30):
                out.append(f'<circle cx="{x+207+offset}" cy="358" r="24" fill="#17445c" stroke="#86d6f0"/>')
        elif index == 1:
            out += [f'<path d="M {x+90} 330 Q {x+200} 380 {x+323} 330 M {x+206} 358 L {x+206} 405" fill="none" stroke="#86d6f0" stroke-width="9"/>']
            for dx, y in ((90, 330), (323, 330), (206, 405)):
                out.append(f'<circle cx="{x+dx}" cy="{y}" r="19" fill="#061320" stroke="#e6be7a" stroke-width="2"/>')
        else:
            out += [f'<path d="M {x+88} 355 L {x+210} 355" stroke="#86d6f0" stroke-width="3"/>']
            for dx, color in ((88, "#86d6f0"), (210, "#86d6f0"), (331, "#b497dc")):
                out.append(f'<circle cx="{x+dx}" cy="355" r="25" fill="#061320" stroke="{color}" stroke-width="2"/>')
            out += [_text(x+310, 402, "hidden", 16, "#b497dc")]
        for j, formula in enumerate(formulas):
            out.append(_text(x+25, 485+j*42, formula, 18))
        for j, line in enumerate(limits):
            out.append(_text(x+25, 639+j*32, line, 15, "#e6be7a"))
        out += [_text(x+25, 755, source, 17, "#adc0cd")]
    out += [_text(55, 858, "A published candidate stays inside its source domain. An untested bridge remains visible as an open connection.", 19, "#adc0cd"), '</svg>']
    return "\n".join(out) + "\n"


def narrative(record):
    return """# A cosmos seen from its possible exterior

![Finite exterior-first patch map](map.svg)

Start outside the local cosmos. Two black holes in a parent region may be
interfaces to two children. The pictured parent is not assumed globally
spherical: spherical symmetry belongs to separate local patches. A child's
possible closed spatial slice is S3; its round horizon sections are S2.
The drawing gives neither a global spatial embedding nor a solved spacetime.

The children may inherit a common law while differing in state and history.
Descent, geometric attachment, causal influence and correlation are four
independent relations. Shared ancestry alone gives no signal route. The
grandchild marks a finite second generation; its formation interface is open.

## Borrow before extending

![Three candidate families](candidates.svg)

The one-hole multiple-daughter construction supplies null-shell matching and
conditional boundary collisions. Its l>a disconnection condition uses a
near-singularity approximation, equal nucleation radius, and a<<R0<<2m.
It differs from two separate black holes
hosting different children. The connected-mouth construction supplies actual
conditional traversability, within its charged, many-species semiclassical
domain. Neither is imported as the project's completed geometry.

The response-equivalence control asks a different question: can a local
observer distinguish all exteriors? An exactly uncoupled mode changes the
exterior spectrum without changing the local resolvent. This known algebra
is a control against overclaiming what a local measurement can reconstruct.

## The next useful mathematical connections

1. Supply a common parent geometry and match each local interface, carrying
   induced metric, orientation, source and conservation together.
2. Supply the state and clock conversions that turn an exterior response
   into a local observer prediction. Do not turn correlation into signaling.
3. Separate inherited dimensionless law, environmental response and chosen
   state. Hold each candidate theory fixed when perturbing its state.
4. Then ask which local dark-sector, matter-antimatter or time observables
   can distinguish these exterior candidates. No abundance fit starts this map.

## What the saved records establish

`atlas.json` records every region, connection, source, assumption and gap.
`controls.json` checks imported junction and response algebra. The separately
labeled `integration.json` passes an explicitly supplied four-mode matrix to
the existing retained-region reducer. Its matrix coefficients are software
fixtures; no physical Hamiltonian, formation probability or relativistic
light cone is derived from them. Numerical discrepancies are not error bounds.

`manifest.json` binds portable source paths, the canonical input, and each
output. The producer refuses an existing output directory; `--check` verifies
the existing bundle without rewriting it. Full source reading boundaries
and versioned arXiv URLs are included in `atlas.json`.
"""


def png_preview(svg):
    """Raster preview of our small SVG vocabulary using the existing matplotlib.

    Preview uses a flat equivalent of the decorative radial gradient. All
    geometry/text come from the SVG, so the preview cannot drift in content.
    This is deliberately not a general SVG interpreter.
    """
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.patches import Circle, FancyBboxPatch, PathPatch
    from matplotlib.path import Path as MplPath

    root = ET.fromstring(svg)
    width, height = float(root.attrib["width"]), float(root.attrib["height"])
    fig = Figure(figsize=(width / 100, height / 100), dpi=100, facecolor="#07111e")
    FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.set_axis_off()
    for element in root:
        kind = element.tag.split("}")[-1]
        data = element.attrib
        fill = data.get("fill", "none")
        fill = "#123a54" if fill.startswith("url(") else fill
        stroke = data.get("stroke", "none")
        common = {"facecolor": fill, "edgecolor": stroke,
                  "linewidth": float(data.get("stroke-width", 1)) * .72,
                  "alpha": float(data.get("opacity", 1))}
        if "stroke-dasharray" in data:
            common["linestyle"] = (0, tuple(float(x) * .72 for x in data["stroke-dasharray"].replace(",", " ").split()))
        if kind == "rect":
            ax.add_patch(FancyBboxPatch((float(data.get("x", 0)), float(data.get("y", 0))),
                float(data["width"]), float(data["height"]),
                boxstyle=f'round,pad=0,rounding_size={data.get("rx", 0)}', **common))
        elif kind == "circle":
            ax.add_patch(Circle((float(data["cx"]), float(data["cy"])), float(data["r"]), **common))
        elif kind == "text":
            ax.text(float(data["x"]), float(data["y"]), element.text or "",
                    fontsize=float(data.get("font-size", 18)) * .72, color=fill,
                    family="DejaVu Sans", weight="bold" if data.get("font-weight") == "600" else "normal",
                    va="baseline", ha="left")
        elif kind == "path":
            tokens = re.findall(r"[MLCQz]|[-+]?(?:\d*\.\d+|\d+)", data["d"])
            vertices, codes, cursor = [], [], 0
            while cursor < len(tokens):
                command = tokens[cursor]
                cursor += 1
                count = {"M": 1, "L": 1, "Q": 2, "C": 3, "z": 0}[command]
                code = {"M": MplPath.MOVETO, "L": MplPath.LINETO,
                        "Q": MplPath.CURVE3, "C": MplPath.CURVE4, "z": MplPath.CLOSEPOLY}[command]
                if not count:
                    vertices.append(vertices[0])
                    codes.append(code)
                for _ in range(count):
                    vertices.append((float(tokens[cursor]), float(tokens[cursor+1])))
                    codes.append(code)
                    cursor += 2
            ax.add_patch(PathPatch(MplPath(vertices, codes), **common))
            if data.get("marker-end") and len(vertices) >= 2:
                ax.annotate("", xy=vertices[-1], xytext=vertices[-2],
                            arrowprops={"arrowstyle": "-|>", "color": stroke, "lw": 0,
                                        "mutation_scale": 15})
    stream = BytesIO()
    fig.savefig(stream, format="png", dpi=100, metadata={"Software": "NSC finite atlas SVG preview"})
    return stream.getvalue()


def payloads():
    record = atlas.build_atlas()
    map_svg, candidates_svg = render_default(record), render_candidates()
    return {
        "atlas.json": atlas.canonical_bytes(record),
        "controls.json": atlas.canonical_bytes(atlas.compatibility_controls()),
        "integration.json": atlas.canonical_bytes(atlas.integration_control()),
        "map.svg": map_svg.encode("utf-8"),
        "candidates.svg": candidates_svg.encode("utf-8"),
        "map.png": png_preview(map_svg),
        "candidates.png": png_preview(candidates_svg),
        "narrative.md": narrative(record).encode("utf-8"),
    }


def _close(actual, expected):
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        return actual == expected
    if isinstance(expected, (int, float)):
        return isinstance(actual, (int, float)) and atlas.math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12)
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(_close(a, b) for a, b in zip(actual, expected))
    return isinstance(actual, dict) and actual.keys() == expected.keys() and all(_close(actual[k], expected[k]) for k in expected)


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only atlas output already exists; use --check or a new version")
    before = source_bindings()
    outputs = payloads()
    if source_bindings() != before:
        raise RuntimeError("source changed during atlas production")
    manifest = {"schema": atlas.SCHEMA, "source_bindings": before,
        "input_sha256": atlas.hashlib.sha256(outputs["atlas.json"]).hexdigest(),
        "outputs": {name: {"sha256": atlas.hashlib.sha256(blob).hexdigest(), "bytes": len(blob)} for name, blob in outputs.items()},
        "runtime": {"python": platform.python_version(),
                    **{name: version(name) for name in ("numpy", "scipy", "sympy", "matplotlib")}},
        "provenance": "repository-relative byte hashes; no historical producer revision claimed",
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12}}
    # Atomic directory claim: a second producer cannot open or overwrite it.
    directory.mkdir(parents=True, exist_ok=False)
    for name, blob in outputs.items():
        with (directory / name).open("xb") as handle:
            handle.write(blob)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(atlas.canonical_bytes(manifest))
    return {"created": directory.name, "outputs": len(outputs), "new_physics_result": False}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != atlas.SCHEMA or manifest["source_bindings"] != source_bindings():
        raise ValueError("atlas source binding mismatch")
    expected = payloads()
    if set(manifest["outputs"]) != set(expected):
        raise ValueError("atlas output inventory mismatch")
    for name, replay in expected.items():
        path = directory / name
        binding = manifest["outputs"][name]
        if _sha(path) != binding["sha256"] or path.stat().st_size != binding["bytes"]:
            raise ValueError(f"atlas artifact binding mismatch: {name}")
        equal = _close(json.loads(path.read_bytes()), json.loads(replay)) if name.endswith(".json") else path.read_bytes() == replay
        if not equal:
            raise ValueError(f"atlas replay mismatch: {name}")
    if manifest["input_sha256"] != _sha(directory / "atlas.json"):
        raise ValueError("atlas canonical input binding mismatch")
    return {"checked": directory.name, "sources": len(manifest["source_bindings"]), "outputs": len(expected)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
