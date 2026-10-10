#!/usr/bin/env python3
"""Create or verify the shared-parent atlas and finite fixed-source screen."""
from __future__ import annotations

import argparse
from copy import deepcopy
from io import BytesIO
import json
from pathlib import Path
import platform
from importlib.metadata import version

import numpy as np

import derive_nsc_cosmological_atlas as first
from recursive_horizons import nsc_atlas_parent_initial_data as parent
from recursive_horizons import nsc_atlas_junction_sources as junction

LAB, REPO = first.LAB, first.REPO
DEFAULT = LAB / "results/development/nsc-cosmological-atlas-v2"
SCHEMA = "nsc-cosmological-atlas-v2"


def build_record():
    """Refine one spatial baseline; preserve independent scenario assumptions."""
    record = first.atlas.build_atlas()
    p, j = parent.parent_record(), junction.junction_record()
    record["sources"]["bl_initial_data"] = deepcopy(p["source"])
    record["sources"]["fixed_wall_screen"] = {
        "path": "lab/docs/nsc-atlas-junction-sources.md",
        "use": "Finite source compatibility in a separate spherical problem",
        "limit": "Not a junction embedded in the Brill-Lindquist parent",
    }
    scenario = record["scenarios"][0]
    scenario["title"] = "Two unequal holes in a shared spatial parent baseline"
    scenario["theory"] = {
        "family": "vacuum Einstein initial data; child source unselected",
        "couplings": {"G": 1, "c": 1},
        "rule": "units and initial data do not determine a formation law",
    }
    scenario["solution_status"] = "imported spatial initial data; physical child interfaces unresolved"
    regions = {r["id"]: r for r in scenario["regions"]}
    regions["parent"].update(
        title="Brill-Lindquist spatial parent baseline",
        geometry={"chart": "one time-symmetric spatial slice",
                  "metric": "h_ij=psi^4 delta_ij; K_ij=0",
                  "spatial_topology": "R3 minus two punctures; three asymptotically flat ends",
                  "symmetry_orbits": "axial symmetry only, not globally spherical",
                  "globally_spherical": False, "embedding": "none; conformal coordinates"},
        source="vacuum Einstein constraint solution, away from puncture coordinates",
        state="classical initial data; no quantum state supplied",
        clock={"kind": "initial slice only; lapse and interpatch clocks not supplied"},
        dimensionless_parameters=deepcopy(p["dimensionless_parameters"]),
        validity="spatial initial data only; no global spacetime evolution supplied",
    )
    for index, key in enumerate(("a", "b")):
        hole = regions[f"hole_{key}"]
        hole.update(
            title=f"Puncture neighborhood {key.upper()}",
            geometry={"chart": "restriction of the common Brill-Lindquist slice",
                      "metric": "h_ij=psi^4 delta_ij with both puncture terms",
                      "symmetry_orbits": "not exact S2 symmetry at finite separation",
                      "spatial_topology": "neighborhood of an additional asymptotically flat end",
                      "coordinate_position": p["inputs"]["positions"][index],
                      "embedding": "none; coordinate placement on an initial slice",
                      "exact_schwarzschild_patch": False},
            source="shared vacuum initial data; companion contribution retained",
            clock={"kind": "no static Killing clock supplied for the binary spacetime"},
            dimensionless_parameters={
                "bare_puncture_mass_over_L": p["inputs"]["bare_masses"][index],
                "end_ADM_mass_over_L": p["adm_masses"]["puncture_ends"][index]},
            validity="puncture coordinate excluded; horizons not located",
        )
        edge = next(e for e in scenario["links"] if e["id"] == f"parent-hole-{key}")
        edge.update(source_refs=["bl_initial_data"], status="imported_spatial_restriction",
                    evidence="Both neighborhoods use one published vacuum spatial metric",
                    missing=["spacetime development", "physical matching to a child interface"])
    scenario["assumptions"] = [
        "time-symmetric, nonspinning, vacuum parent initial slice",
        "asymptotic ends are not offspring or established communicating regions",
        "unequal initial puncture parameters; finite separation",
        "possible children and grandchild are still hypotheses, cutoff two",
    ]
    record["assessments"] = {
        "shared_parent": {"record": "parent.json", "status": "IMPORTED_INITIAL_DATA",
                          "spacetime_development_supplied": False},
        "fixed_source_junction": {"record": "junction.json", "status": j["status"],
            "scope": "separate comoving closed-radiation/Schwarzschild/constant-wall candidate",
            "embedded_in_shared_parent": False, "all_child_attachments_excluded": False},
    }
    record["admission"].update(
        new_physics_result=False,
        category="imported parent geometry and scoped source-compatibility application",
        missing_arrow="source-supported child dynamics and nonspherical parent matching",
    )
    first.atlas.validate_atlas(record)
    record["schema"] = SCHEMA
    record["predecessor"] = "nsc-cosmological-atlas-v1"
    return record, p, j


def source_bindings():
    bindings = first.source_bindings()
    paths = [Path(__file__), Path(parent.__file__), Path(junction.__file__),
        LAB / "tests/test_nsc_atlas_parent_initial_data.py",
        LAB / "tests/test_nsc_atlas_junction_sources.py",
        LAB / "tests/test_nsc_cosmological_atlas_v2.py",
        LAB / "docs/nsc-atlas-parent-initial-data.md",
        LAB / "docs/nsc-atlas-junction-sources.md",
        LAB / "docs/nsc-cosmological-atlas-v2.md",
        LAB / "docs/child-neck/balanced-layers-note.md",
        REPO / "paper/atlas-working-draft.md",
    ]
    bindings.update({p.relative_to(REPO).as_posix(): first._sha(p) for p in paths})
    return dict(sorted(bindings.items()))


def render_map(record, p, j):
    text, panel = first._text, first._panel
    out = first._base(1400, 1230, "A shared parent, and one attachment decided",
        "Atlas v2 | Published spatial geometry + a separate finite source test")
    out += [panel(50, 130, 1300, 345, "#102c41", "#47758e"),
        text(80, 175, "SHARED PARENT INITIAL SLICE", 23, "#86d6f0", "600"),
        text(80, 210, "Brill-Lindquist: h = psi^4 delta, K = 0; both holes contribute to the same geometry.", 19),
        text(80, 242, "Vacuum constraints satisfied away from punctures. No spacetime evolution or horizon map supplied.", 17, "#adc0cd")]
    for i, x in enumerate((365, 970)):
        out += [f'<circle cx="{x}" cy="335" r="44" fill="#020810" stroke="#f4ce8c" stroke-width="2"/>',
            text(x-16, 343, "AB"[i], 24, "#f4ce8c", "600"),
            text(x-180, 407, f'Bare parameter {p["inputs"]["bare_masses"][i]:g}; end ADM {p["adm_masses"]["puncture_ends"][i]:g}', 18),
            text(x-180, 439, "Asymptotic end; offspring not established", 16, "#e6be7a")]
    out += [text(580, 337, f'coordinate d/L = {p["inputs"]["coordinate_separation"]:g}', 18, "#adc0cd")]
    for name, x, center in (("A", 50, 365), ("B", 735, 970)):
        out += [f'<path d="M {center} 475 L {center} 555" stroke="#e6be7a" stroke-width="2" stroke-dasharray="8 7" marker-end="url(#arrow)"/>',
            text(center+15, 515, "formation + matching open", 16, "#e6be7a"),
            panel(x, 565, 615, 160), text(x+25, 608, f"POTENTIAL CHILD {name}", 22, "#86d6f0", "600"),
            text(x+25, 646, "Its geometry, source and state require an attachment law.", 18),
            text(x+25, 684, "No causal path or cross-universe age is assigned.", 17, "#adc0cd")]
    out += ['<path d="M 300 725 L 300 770" stroke="#e6be7a" stroke-width="2" stroke-dasharray="8 7" marker-end="url(#arrow)"/>',
        panel(65, 778, 580, 90), text(90, 814, "POTENTIAL GRANDCHILD", 20, "#86d6f0", "600"),
        text(90, 845, "A further formation interface is still missing.", 17),
        text(755, 801, "Finite descent stops at generation two.", 18),
        text(755, 838, "Shared history does not supply sibling signaling.", 17, "#adc0cd"),
        panel(50, 900, 1300, 172, "#382633", "#bf7881"),
        text(80, 942, "SEPARATE SOURCE SCREEN: THIS COMBINATION IS INCOMPATIBLE", 21, "#f0b1b8", "600"),
        text(80, 980, "Comoving radiation-filled child + fixed Schwarzschild mass + constant-tension wall", 20),
        text(80, 1016, "Excluded on an expanding interval by the junction equation. Moving walls and other interiors remain open.", 17),
        text(80, 1048, "This spherical test has not been embedded in the two-hole parent above.", 16, "#e6be7a"),
        text(55, 1120, "Relations remain separate: descent / geometric attachment / causal influence / correlation.", 19, "#86d6f0"),
        text(55, 1157, "One-hole multiple daughters, communicating mouths and indistinguishable exteriors retain their own assumptions.", 17, "#adc0cd"),
        text(55, 1194, "Patch sizes are schematic. The separate parent-slice figure plots conformal coordinates, not a physical embedding.", 16, "#adc0cd"), "</svg>"]
    return "\n".join(out) + "\n"


def scientific_figures(p, j):
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    outputs = {}
    with matplotlib.rc_context({"svg.hashsalt": SCHEMA, "font.family": "DejaVu Sans",
        "font.size": 12, "axes.titlesize": 17, "figure.facecolor": "white"}):
        def save(fig, name):
            for fmt in ("svg", "png"):
                stream = BytesIO()
                meta = {"Date": None, "Creator": SCHEMA} if fmt == "svg" else {"Software": SCHEMA}
                fig.savefig(stream, format=fmt, dpi=145, metadata=meta)
                outputs[f"{name}.{fmt}"] = stream.getvalue()
            plt.close(fig)

        z, x = np.meshgrid(np.linspace(-10, 10, 181), np.linspace(-5, 5, 101))
        points = np.stack((x, np.zeros_like(x), z), axis=-1)
        positions = np.asarray(p["inputs"]["positions"])
        distance = np.linalg.norm(points[..., None, :] - positions, axis=-1)
        mask = np.any(distance < .28, axis=-1)
        psi = np.full(x.shape, np.nan)
        psi[~mask] = parent.conformal_factor(points[~mask], p["inputs"]["bare_masses"], positions)
        fig, ax = plt.subplots(figsize=(13, 7), layout="constrained")
        field = ax.contourf(z, x, np.log10(psi), levels=np.linspace(0, .75, 31), cmap="magma", extend="max")
        ax.contour(z, x, psi, levels=[1.12, 1.2, 1.35, 1.6, 2.0], colors="white", linewidths=.65, alpha=.8)
        for i, pos in enumerate(positions):
            ax.plot(pos[2], pos[0], "x", color="#70ddfa", markersize=10, markeredgewidth=2)
            ax.annotate(f'Puncture {"AB"[i]}\n$m/L={p["inputs"]["bare_masses"][i]:g}$',
                        (pos[2], pos[0]), xytext=(12, 22), textcoords="offset points", color="white", fontsize=12)
        ax.set(xlabel="Conformal coordinate z / L", ylabel="Conformal coordinate x / L", aspect="equal",
               title="Two unequal holes: one Brill-Lindquist spatial metric")
        fig.colorbar(field, ax=ax, label=r"$\log_{10}\psi$, with $h_{ij}=\psi^4\delta_{ij}$")
        fig.supxlabel("y = 0 slice; small coordinate disks masked. Crosses are punctures, not measured horizons. No time evolution shown.", fontsize=11)
        save(fig, "parent-slice")

        rows = j["samples"]
        tau = [r["tau_over_A"] for r in rows]
        fig, axes = plt.subplots(1, 2, figsize=(13, 5.7), layout="constrained")
        for ax, key, wallkey, label in (
                (axes[0], "A_sigma_required", "A_sigma_wall", r"$A\sigma$"),
                (axes[1], "A_pressure_required", "A_pressure_wall", r"$AP$")):
            ax.plot(tau, [r[key] for r in rows], color="#175778", linewidth=2.5, label="Required by prescribed geometry")
            ax.plot(tau, [r[wallkey] for r in rows], color="#b54b31", linestyle="--", linewidth=2, label="Prescribed constant-tension wall")
            ax.set(xlabel=r"Shell proper time $\tau/A$", ylabel=label)
            ax.grid(alpha=.2)
            ax.legend(fontsize=9)
        fig.suptitle("One fixed wall cannot supply this comoving radiation junction", fontsize=17)
        fig.supxlabel("Finite illustration of the algebraic obstruction; no fit. This spherical test is separate from the two-hole initial slice.", fontsize=10)
        save(fig, "junction")
    return outputs


def narrative(record, p, j):
    return f"""# A shared parent, and one attachment decided

![Updated candidate map](map.svg)

The parent is now a published spatial initial-data construction containing two
unequal punctures. Additional ends are part of its topology; their formation
as children and communication are not established. The same metric contains
both holes, so companion-induced angular dependence is retained.

![Conformal initial slice](parent-slice.svg)

The bare parameters are {p['inputs']['bare_masses']}; the outer ADM mass is
{p['adm_masses']['outer_end']:g} and the two end masses are
{p['adm_masses']['puncture_ends']}, in the chosen geometric unit. These are
imported initial-data quantities, not inherited cosmological constants.
The plot is a conformal-coordinate section, not a spatial embedding.

## A source question with an answer

![Required junction stress and fixed wall](junction.svg)

The source verdict is **{j['status']}**, within its declared assumptions.
A closed radiation FLRW region, a comoving shell, a fixed Schwarzschild mass
and a constant-tension wall cannot satisfy the angular junction throughout an
open interval of changing radius. With k=4 pi sigma and b=cos(chi), it would
require k^2 R^4 - 2 b k R^3 + 2 M R - A^2 sin(chi)^4 = 0 identically.
Its nonzero constant coefficient excludes that identity for A>0 and
0<chi<pi. The independent temporal-curvature and conservation checks support
the implementation; the polynomial argument supplies the interval result.

This is an application of established matching theory, without a priority
claim. The spherical test is not an attachment to the Brill-Lindquist slice.
Moving boundaries, other interiors and other physically specified sources
remain separate possibilities. The earlier adjusted composite layers still
lack a derived fill and exchange law.

## Alternatives remain visible

![Independent candidate families](candidates.svg)

The one-hole multiple-daughter, communicating-mouth and response-equivalence
branches preserve their v1 assumptions. No shared ancestry or initial spatial
connection acquires a signal delay. The existing response adapter still needs
a physically supplied Hamiltonian and observer clock before a local prediction.

## Recorded inputs and replay

[atlas.json](atlas.json), [parent.json](parent.json) and [junction.json](junction.json)
contain the scenario map, source/domain cards, fixed inputs, admission records,
sampled outputs and checks. [manifest.json](manifest.json) binds sources and
outputs. The producer is creation-only; --check replays without rewriting.
The [working manuscript](../../../../paper/atlas-working-draft.md) develops
the argument and distinguishes imported mathematics, this scoped exclusion,
and open physical connections. No dark-sector fraction, antimatter asymmetry,
cross-universe age or constant-inheritance law has been predicted.
"""


def payloads():
    record, p, j = build_record()
    map_svg = render_map(record, p, j)
    candidates_svg = first.render_candidates()
    outputs = {"atlas.json": first.atlas.canonical_bytes(record),
        "parent.json": first.atlas.canonical_bytes(p),
        "junction.json": first.atlas.canonical_bytes(j),
        "map.svg": map_svg.encode(), "map.png": first.png_preview(map_svg),
        "candidates.svg": candidates_svg.encode(), "candidates.png": first.png_preview(candidates_svg),
        "narrative.md": narrative(record, p, j).encode()}
    outputs.update(scientific_figures(p, j))
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
        "scope": "imported initial data and a scoped incompatibility; no global cosmology or novelty claim"}
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
        path = directory / name
        binding = manifest["outputs"][name]
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
