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
    "programme": "results/development/nsc-cosmological-atlas-v8/programme.json",
    "prepared_response": "results/development/nsc-cosmological-atlas-v8/prepared-response.json",
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
    response = records["prepared_response"]
    if (response["status"] != "FINITE_PREPARED_CLASSICAL_RESPONSE_WITH_POSITIVE_TAIL"
            or not response["analytic_tail"]["positive_after_input_in_declared_window"]):
        raise ValueError("The prepared response no longer supports the finite-tail view.")
    if (response["inputs"] != records["preparation"]["inputs"]
            or response["preparation"] != records["preparation"]["preparation"]):
        raise ValueError("The prepared response must retain the certified geometry and input.")
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
    # The manuscript caps a figure at 530 pt tall. At 11.2 in native height,
    # the minimum 13 pt label remains 8.54 pt in the rendered paper.
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 13,
        "svg.hashsalt": "nsc-geometry-overview-2026-10-12",
        "axes.unicode_minus": False,
    })
    fig, ax = plt.subplots(figsize=(9, 11.2))
    fig.patch.set_facecolor(COLORS["white"])
    ax.set(xlim=(0, 9), ylim=(0, 11.2))
    ax.axis("off")

    def label(x, y, text, size=13, color="ink", weight="normal", ha="center", **kw):
        ax.text(x, y, text, fontsize=size, color=COLORS.get(color, color),
                fontweight=weight, ha=ha, va="center", linespacing=1.3, **kw)

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
        label(x, y + (0.12 if sub else 0), text, size=13, weight="bold")
        if sub:
            label(x, y - 0.17, sub, size=13, color="muted")

    def box(x, y, w, h, text, color="attachment", size=13):
        ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                                   boxstyle="round,pad=0.06,rounding_size=0.10",
                                   fc="white", ec=COLORS[color], lw=1.6))
        label(x, y, text, size=size)

    label(0.45, 10.91, "NESTED-SPACE COSMOLOGY", 13, "muted", "bold", ha="left")
    label(0.45, 10.51, "From possible regions to an interior observer", 19,
          "ink", "bold", ha="left")
    for y0, height in ((6.04, 4.10), (1.04, 4.78)):
        ax.add_patch(FancyBboxPatch((0.25, y0), 8.5, height,
                                   boxstyle="round,pad=0.02,rounding_size=0.18",
                                   fc=COLORS["panel"], ec="none"))
    label(0.52, 9.79, "1  A finite branching question", 16, weight="bold", ha="left")

    ax.add_patch(Ellipse((4.5, 8.96), 7.5, 1.15,
                         fc="white", ec=COLORS["line"], lw=1.6))
    label(4.5, 9.29, "Common parent • not globally spherical", 13, weight="bold")
    circle(2.70, 8.84, 0.27, "A")
    circle(6.55, 8.84, 0.33, "B")
    label(4.5, 8.83, "Two off-center\nblack-hole boundaries", 13, color="muted")
    # Parent-to-hole attachment remains a missing common-metric connection.
    arrow((3.46, 9.14), (3.03, 8.99), "attachment", dashed=True)
    arrow((5.57, 9.14), (6.15, 9.02), "attachment", dashed=True)
    circle(2.70, 7.81, 0.56, "Child A?", "candidate")
    circle(6.55, 7.81, 0.51, "Child B?", "may differ")
    for x, top, bottom in ((2.70, 8.52, 8.41), (6.55, 8.45, 8.37)):
        arrow((x-0.13, top), (x-0.13, bottom), "descent", dashed=True)
        arrow((x+0.13, top), (x+0.13, bottom), "attachment", dashed=True)
    arrow((3.34, 7.81), (5.96, 7.81), "correlation", dashed=True, both=True)
    label(4.63, 7.31, "shared ancestry?\nno signal assumed", 13, color="correlation")
    circle(1.44, 6.97, 0.65, "Grandchild?")
    arrow((2.18, 7.54), (1.95, 7.38), "descent", dashed=True)
    arrow((2.37, 7.34), (2.14, 7.18), "attachment", dashed=True)
    label(6.13, 6.79, "Descent remains open.\nMultiple daughters and communicating\nmouths retain their own assumptions.",
          13, color="muted")
    label(4.5, 6.16, "Dashed links: questions, not a formation history", 13, color="muted")

    label(0.52, 5.49, "2  One checked connection", 16, weight="bold", ha="left")
    label(4.5, 5.11, "Sourced initial slice • fixed wall parameters", 13, color="muted")
    for x, text in ((1.35, "Regular\ncap"), (3.45, "Material\nwall"),
                    (5.55, "Vacuum\nbridge"), (7.65, "Right\nend")):
        box(x, 4.61, 1.44, 0.61, text)
    for x in (2.16, 4.26, 6.36):
        arrow((x, 4.61), (x+0.48, 4.61), "attachment")
    label(4.5, 3.94, "End-only gluing preserves this body in a\nsufficiently separated initial-data family.",
          13, color="muted")
    ax.plot([0.55, 8.45], [3.65, 3.65], color=COLORS["line"], lw=1)
    label(4.5, 3.39, "Isolated reference • event-specific observer access", 13, weight="bold")
    for x, name, sub in ((1.35, "C", "center"), (4.50, "W", "wall"), (7.65, "P", "right parent")):
        circle(x, 2.72, 0.23, name, color="causal")
        label(x, 2.32, sub, 13, color="muted")
    arrow((1.67, 2.72), (4.18, 2.72), "causal", both=True)
    arrow((4.82, 2.72), (7.33, 2.72), "causal")
    label(2.93, 3.00, "early paths", 13, color="causal")
    label(6.08, 3.00, "early escape", 13, color="causal")
    label(4.5, 1.92, "After the anchor: no new cross-observer emission\nreaches another of these three worldlines.", 13, color="muted")
    box(4.5, 1.26, 7.58, 0.59,
        "Wall-input calibration → transient at C\nFinite IV data → calculated pulse + positive finite tail",
        color="causal", size=13)
    label(4.5, 0.72, "Remote preparation and later evolution remain open", 13, color="muted")

    for x, kind, name in ((0.39, "attachment", "Attachment"), (2.76, "descent", "Descent"),
                          (4.62, "causal", "Causal influence"), (7.18, "correlation", "Correlation")):
        arrow((x, 0.37), (x+0.35, 0.37), kind)
        label(x+0.45, 0.37, name, 13, ha="left")
    label(4.5, 0.13, "Schematic circles and spacing; local spherical charts do not fix global topology.",
          13, color="muted")
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
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
        "layout": {"orientation": "portrait", "native_inches": [9, 11.2],
                   "minimum_font_points": 13, "paper_height_cap_points": 530,
                   "left_panel_display": "upper", "right_panel_display": "lower"},
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
                "waveform_evaluated_in_v7_certificate": False,
                "remote_generation_established": False,
            },
            "prepared_response": {
                "status": records["prepared_response"]["status"],
                "domain": records["prepared_response"]["domain"],
                "metrics": records["prepared_response"]["metrics"],
                "analytic_tail": records["prepared_response"]["analytic_tail"],
                "waveform_evaluated": True,
                "permanent_memory_established": False,
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
                 "status": "calculated pulse and positive post-input tail within the finite admitted window"},
            ],
        },
        "limits": [
            "No spatial embedding, global topology or relative physical scales are depicted.",
            "No parent observer prepares the early scalar trace in the isolated reference.",
            "The admitted IV data are distinct from the old prescribed wall trace and have no modelled remote source.",
            "The positive prepared-field tail is confined to the admitted early window; later evolution and permanent memory remain open.",
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
