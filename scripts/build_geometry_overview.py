"""Render the living paper's overview from stored atlas records only.

No field integration or geometric solve occurs here. The JSON companion binds
the selected record excerpts and their hashes to explicit diagram semantics.
Run from any directory with the repository's scientific Python interpreter.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "motif": "results/development/nsc-cosmological-atlas-v1/atlas.json",
    "slice": "results/development/nsc-cosmological-atlas-v4/connection-slice.json",
    "assembly": "results/development/nsc-cosmological-atlas-v4/parent-assembly.json",
    "access": "results/development/nsc-cosmological-atlas-v5/observer-access.json",
    "field": "results/development/nsc-cosmological-atlas-v6/field-transfer.json",
    "preparation": "results/development/nsc-cosmological-atlas-v7/preparation-access.json",
    "programme": "results/development/nsc-cosmological-atlas-v7/programme.json",
}
COLORS = {
    "ink": "#183346", "muted": "#536977", "panel": "#f1f5f7",
    "attachment": "#2879a0", "descent": "#b27920", "causal": "#21816c",
    "correlation": "#865a98", "white": "#ffffff", "line": "#c9d5dc",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_evidence() -> tuple[dict, dict]:
    records, provenance = {}, {}
    for name, relative in INPUTS.items():
        raw = (ROOT / relative).read_bytes()
        records[name] = json.loads(raw)
        provenance[name] = {"path": relative, "sha256": sha(raw)}
    if records["slice"]["status"] != "CONNECTED_SOURCED_INITIAL_SLICE":
        raise ValueError("The sourced-slice record no longer supports this overview.")
    if records["field"]["status"] != "CONDITIONAL_EARLY_CLASSICAL_SCALAR_RESPONSE":
        raise ValueError("The field record no longer supports this overview.")
    return records, provenance


def render(output: Path) -> dict:
    records, provenance = load_evidence()
    motif = next(s for s in records["motif"]["scenarios"] if s["id"] == "separate_holes")
    wanted = {
        "parent-hole-a", "parent-hole-b", "descent-a", "descent-b", "neck-a",
        "neck-b", "siblings", "grandchild-descent", "grandchild-neck",
    }
    links = [link for link in motif["links"] if link["id"] in wanted]
    if {link["id"] for link in links} != wanted:
        raise ValueError("The original motif's selected connections have changed.")
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 11,
        "svg.hashsalt": "nsc-geometry-overview-2026-10-12",
        "axes.unicode_minus": False,
    })
    fig, ax = plt.subplots(figsize=(16, 10))
    fig.patch.set_facecolor(COLORS["white"])
    ax.set(xlim=(0, 16), ylim=(0, 10))
    ax.axis("off")

    def label(x, y, text, size=11, color="ink", weight="normal", ha="center", **kw):
        ax.text(x, y, text, fontsize=size, color=COLORS.get(color, color),
                fontweight=weight, ha=ha, va="center", linespacing=1.4, **kw)

    def arrow(start, end, kind, dashed=False, both=False, rad=0.0):
        ax.add_patch(FancyArrowPatch(
            start, end, arrowstyle="<->" if both else "-|>",
            mutation_scale=14, lw=1.8, color=COLORS[kind],
            linestyle=(0, (4, 3)) if dashed else "-",
            connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0,
        ))

    def circle(x, y, radius, text, sub=None, color="attachment"):
        ax.add_patch(Circle((x, y), radius, fc=COLORS["white"],
                            ec=COLORS[color], lw=1.8))
        label(x, y + (0.12 if sub else 0), text, size=11, weight="bold")
        if sub:
            label(x, y - 0.23, sub, size=9, color="muted")

    def box(x, y, w, h, text, color="attachment", size=10):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                                   boxstyle="round,pad=0.08,rounding_size=0.12",
                                   fc="white", ec=COLORS[color], lw=1.6))
        label(x, y, text, size=size)

    label(0.55, 9.58, "NESTED-SPACE COSMOLOGY", 12, "muted", "bold", ha="left")
    label(0.55, 9.10, "From possible regions to a measured interior", 23,
          "ink", "bold", ha="left")
    for x0 in (0.35, 8.15):
        ax.add_patch(FancyBboxPatch((x0, 1.35), 7.5, 7.13,
                                   boxstyle="round,pad=0.02,rounding_size=0.18",
                                   fc=COLORS["panel"], ec="none"))
    label(0.65, 8.08, "1  A finite branching question", 15, weight="bold", ha="left")
    label(8.45, 8.08, "2  One checked connection", 15, weight="bold", ha="left")

    ax.add_patch(Ellipse((4.05, 6.78), 6.3, 1.83,
                         fc="white", ec=COLORS["line"], lw=1.6))
    label(4.05, 7.42, "Common parent • not globally spherical", 11, weight="bold")
    circle(2.50, 6.60, 0.37, "A")
    circle(5.45, 6.60, 0.46, "B")
    label(4.05, 6.04, "Two off-center black-hole boundaries", 10, color="muted")
    # Parent-to-hole attachment remains a missing common-metric connection.
    arrow((3.65, 7.05), (2.81, 6.85), "attachment", dashed=True)
    arrow((4.44, 7.05), (5.04, 6.90), "attachment", dashed=True)
    circle(2.50, 4.69, 0.71, "Child A?", "candidate")
    circle(5.45, 4.69, 0.62, "Child B?", "may differ")
    for x, top in ((2.50, 6.17), (5.45, 6.10)):
        arrow((x-0.13, top), (x-0.13, 5.44), "descent", dashed=True)
        arrow((x+0.13, top), (x+0.13, 5.44), "attachment", dashed=True)
    arrow((3.27, 4.69), (4.75, 4.69), "correlation", dashed=True, both=True)
    label(4.02, 4.24, "shared ancestry?\nno signal assumed", 9, color="correlation")
    circle(2.50, 2.67, 0.60, "Grandchild?", "one more level")
    arrow((2.37, 3.91), (2.37, 3.32), "descent", dashed=True)
    arrow((2.63, 3.91), (2.63, 3.32), "attachment", dashed=True)
    label(5.35, 2.96, "Descent is still open.\nMultiple daughters and\ncommunicating mouths retain\ntheir own assumptions.",
          11, color="muted")
    label(4.05, 1.73, "Dashed links: a question, not a formation history", 10, color="muted")

    label(12.0, 7.50, "Sourced initial slice • fixed wall parameters", 10, color="muted")
    for x, text, w in ((9.17, "Regular\ncap", 1.14), (10.96, "Material\nwall", 1.14),
                       (12.75, "Vacuum\nbridge", 1.14), (14.54, "Right\nend", 1.14)):
        box(x, 6.76, w, 0.83, text)
    for x in (9.82, 11.61, 13.40):
        arrow((x, 6.76), (x+0.47, 6.76), "attachment")
    label(12.0, 5.95, "End-only gluing preserves this body in a\nsufficiently separated initial-data family.",
          10, color="muted")
    ax.plot([8.6, 15.15], [5.47, 5.47], color=COLORS["line"], lw=1)
    label(12.0, 5.10, "Isolated reference • event-specific observer access", 11, weight="bold")
    for x, name, sub in ((9.17, "C", "center"), (11.84, "W", "wall"), (14.54, "P", "right parent")):
        circle(x, 4.37, 0.34, name, color="causal")
        label(x, 3.81, sub, 10, color="muted")
    arrow((9.60, 4.37), (11.41, 4.37), "causal", both=True)
    arrow((12.27, 4.37), (14.11, 4.37), "causal")
    label(10.49, 4.75, "early paths", 9, color="causal")
    label(13.18, 4.75, "early escape", 9, color="causal")
    label(12.0, 3.30, "After the anchor: no new cross-observer emission\nreaches another of these three worldlines.", 10, color="muted")
    box(12.0, 2.34, 6.13, 0.78,
        "Wall-input calibration → transient at C\nFinite IV data → admitted transfer operator",
        color="causal", size=10)
    label(12.0, 1.73, "IV waveform and remote generation remain open", 10, color="muted")

    for x, kind, name in ((0.66, "attachment", "Attachment"), (4.32, "descent", "Descent"),
                          (7.83, "causal", "Causal influence"), (12.00, "correlation", "Correlation")):
        arrow((x, 0.90), (x+0.56, 0.90), kind)
        label(x+0.71, 0.90, name, 11, ha="left")
    label(8.0, 0.30, "Circles and spacing are schematic. Local spherical charts do not establish global topology.",
          10, color="muted")
    fig.subplots_adjust(left=0.015, right=0.985, top=0.985, bottom=0.015)
    output.mkdir(parents=True, exist_ok=True)
    for suffix in ("svg", "png"):
        fig.savefig(output / f"geometry-overview.{suffix}", dpi=190,
                    metadata={"Date": None} if suffix == "svg" else {"Software": "NSC stored-evidence overview"})
    plt.close(fig)
    manifest = {
        "schema": "nsc.geometry-overview.v1",
        "purpose": "Living-manuscript narrative figure; no new scientific solve",
        "producer": "scripts/build_geometry_overview.py",
        "producer_sha256": sha(Path(__file__).read_bytes()),
        "inputs": provenance,
        "left_panel": {
            "scenario": "separate_holes",
            "region_ids": [r["id"] for r in motif["regions"]],
            "connections": links,
            "status": "Open motif; v1 links are not promoted by later separate models",
        },
        "right_panel": {
            "slice_status": records["slice"]["status"],
            "assembly_scope": "Theorem-backed sufficiently separated initial-data family, not solved global future",
            "observer_access_rows": records["access"]["rows"],
            "field_status": records["field"]["status"],
            "incoming_state": records["field"]["incoming_state"],
            "field_metrics": records["field"]["metrics"],
            "preparation": {
                "geometry": records["preparation"]["geometry"],
                "operator": records["preparation"]["operator"],
                "waveform_evaluated": False,
                "remote_generation_established": False,
            },
            "displayed_connections": [
                {"from": "cap", "to": "wall", "kind": "attachment", "status": "checked initial slice"},
                {"from": "wall", "to": "bridge", "kind": "attachment", "status": "checked initial slice"},
                {"from": "bridge", "to": "right end", "kind": "attachment", "status": "checked initial slice"},
                {"from": "C", "to": "W", "kind": "causal", "status": "selected early emissions only"},
                {"from": "W", "to": "C", "kind": "causal", "status": "selected early emissions only"},
                {"from": "W", "to": "P", "kind": "causal", "status": "early white-hole escape only"},
                {"from": "prescribed incoming trace", "to": "C", "kind": "causal",
                 "status": "conditional classical scalar response; exterior origin uncomputed"},
                {"from": "independent finite IV data", "to": "wall and centre", "kind": "causal",
                 "status": "finite conditional transfer operator; waveform not evaluated"},
            ],
        },
        "limits": [
            "No spatial embedding, global topology or relative physical scales are depicted.",
            "No parent observer prepares the early scalar trace in the isolated reference.",
            "The admitted IV data are distinct from the old prescribed wall trace and have no modelled remote source.",
            "Early reception after the anchor is compatible with exclusions on new emissions.",
            "No second cosmological child, regular formation, stability or dark component is established.",
        ],
    }
    manifest["outputs"] = {
        f"geometry-overview.{suffix}": sha((output / f"geometry-overview.{suffix}").read_bytes())
        for suffix in ("svg", "png")
    }
    (output / "geometry-overview.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "paper/figures")
    args = parser.parse_args()
    record = render(args.output_dir.resolve())
    print(json.dumps({"output_dir": str(args.output_dir), "inputs": len(record["inputs"]),
                      "new_physics_calculations": 0}))


if __name__ == "__main__":
    main()
