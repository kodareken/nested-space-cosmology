#!/usr/bin/env python3
"""Render the canonical Markdown manuscript as a publication-form PDF."""

from __future__ import annotations

import argparse
import hashlib
import html
from io import BytesIO
import json
import math
import os
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE = REPOSITORY / "paper" / "nested-space-cosmology.md"
OUTPUT = REPOSITORY / "paper" / "nested-space-cosmology.pdf"
METADATA = REPOSITORY / "paper" / "metadata.json"
BUILD_MANIFEST = REPOSITORY / "paper" / "build-manifest.json"
BUILD_DIRECTORY = REPOSITORY / "paper" / ".build"


def scientific_build_inputs(source: Path) -> list[dict[str, str]]:
    """Authenticate the renderer and every dataset used by a manuscript figure."""
    paths = ["scripts/build_paper.py"]
    text = source.read_text(encoding="utf-8")
    for marker, path in (
        ("geometric-gap", "results/nsc-3-geometric-chain.json"),
        ("vacuum-work", "results/nsc-6-vacuum-work.json"),
        ("compact-source", "results/development/compact-casimir.json"),
        ("source-matching", "results/development/warped-source.json"),
        ("source-matching", "results/development/compact-matching.json"),
        ("parent-source", "results/development/unruh-state.json"),
        ("parent-source", "results/development/state-regulator.json"),
    ):
        if f"<!-- nsc-figure:{marker} -->" in text:
            paths.append(path)
    for path in (
        "results/development/scale-binding.json",
        "results/development/recursive-source-binding.json",
        "results/development/charged-ctp-neck-source.json",
        "results/development/charged-compact-ctp-completion.json",
        "results/development/nsc-background-projection.json",
        "results/development/nsc-child-metric-backreaction.json",
        "results/development/nsc-constraint-complete-neck.json",
    ):
        if path in text:
            paths.append(path)
    return [{"path": path, "sha256": hashlib.sha256((REPOSITORY/path).read_bytes()).hexdigest()}
            for path in sorted(paths)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="move the generated PDF to a recoverable local build archive",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="build twice without touching tracked files and require identical bytes",
    )
    return parser.parse_args()


def clean_output(path: Path) -> int:
    if not path.exists():
        print(f"nothing to clean: {path}")
        return 0
    trashit = shutil.which("trashit")
    if trashit:
        import subprocess

        subprocess.run([trashit, str(path)], check=True)
        print(f"moved generated PDF to the system Trash: {path}")
        return 0
    recovery = BUILD_DIRECTORY / "removed"
    recovery.mkdir(parents=True, exist_ok=True)
    destination = recovery / f"{path.name}.{time.time_ns()}"
    path.replace(destination)
    print(f"moved generated PDF to recoverable build storage: {destination}")
    return 0


def load_metadata() -> dict[str, object]:
    value = json.loads(METADATA.read_text(encoding="utf-8"))
    required = {
        "author",
        "date",
        "date_iso",
        "project",
        "source_date_epoch",
        "status",
        "subtitle",
        "title",
        "version",
    }
    missing = sorted(required - set(value))
    if missing:
        raise RuntimeError(f"paper metadata is missing: {', '.join(missing)}")
    return value


def compact_source_plot():
    """Plot authenticated stored results; never launch a scientific generator."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    record = json.loads((REPOSITORY / "results/development/compact-casimir.json").read_text())
    figure = Figure(figsize=(6.4, 2.8), dpi=220)
    FigureCanvasAgg(figure)
    left, right = figure.subplots(1, 2)
    phases = record["holonomy"]["cutoff_controls"][-1]["phases"]
    left.plot([row["phase"] for row in phases],
              [row["difference_from_periodic"] for row in phases],
              "o-", color="#2357A6", markersize=4, linewidth=1.2)
    left.axvline(.5, color="#526779", linewidth=.7, linestyle=":")
    left.set_xlabel("effective spatial phase (turns)")
    left.set_ylabel(r"$\Delta E\,r_0$")
    left.set_title("Full free-fermion phase potential", fontsize=9)
    left.set_xticks([0, .25, .5, .75, 1])
    left.grid(True, alpha=.25, linewidth=.5)
    cells = record["domain_family"]["observed"]
    right.plot([row["axial_length"] for row in cells],
               [1000 * row["neck_null"] for row in cells],
               "o-", color="#8C4264", markersize=4, linewidth=1.2)
    right.axhline(0, color="#526779", linewidth=.7)
    right.set_xlabel(r"axial circumference $L_x/r_0$")
    right.set_ylabel(r"$10^3(\rho+p_r)_{\rm int}\,r_0^4$")
    right.set_title("Compact interaction at the neck", fontsize=9)
    right.set_xticks([4, 6, 8, 12])
    right.grid(True, alpha=.25, linewidth=.5)
    for axis in (left, right):
        axis.tick_params(labelsize=7)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
    figure.tight_layout(pad=.6)
    return figure


def source_matching_plot():
    """Two distinct source projections, read only from authenticated records."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    warped = json.loads((REPOSITORY / "results/development/warped-source.json").read_text())
    matched = json.loads((REPOSITORY / "results/development/compact-matching.json").read_text())
    figure = Figure(figsize=(6.4, 2.9), dpi=220)
    FigureCanvasAgg(figure)
    left, right = figure.subplots(1, 2)
    reference = warped["flux_cases"]["1"]
    values = [reference["flat_upper_comparison"], reference["response"]["potential"]]
    left.bar([0, 1], values, color=["#2357A6", "#21766C"], width=.58)
    left.set_xticks([0, 1], ["Unwarped", "Gaussian warp"])
    left.set_ylim(0, 2.7)
    left.set_title("Homogeneous spherical potential", fontsize=9)
    left.set_ylabel(r"$V_2$ (reference length$^{-2}$)")
    for x, y in enumerate(values):
        left.text(x, y+.05, f"{y:.6f}", ha="center", fontsize=7)
    left.text(.5, 2.53, f"change: {100*(values[1]/values[0]-1):.2f}%",
              ha="center", fontsize=8, color="#334155")
    rows = matched["matched_coefficients"]
    cutoffs = [row["matching_cutoff"] for row in rows]
    complement = [row["V_Dirac"] for row in rows]
    light = [nu**4/(4*math.pi)**2 for nu in cutoffs]
    total = [a+b for a, b in zip(complement, light)]
    right.plot(cutoffs, complement, "o-", color="#2357A6", label="Complement", linewidth=1.2, markersize=3)
    right.plot(cutoffs, light, "s-", color="#A15A38", label="Retained light", linewidth=1.2, markersize=3)
    right.plot(cutoffs, total, "--", color="#21766C", label="Sum", linewidth=1.3)
    right.set_title("Matching changes the partition", fontsize=9)
    right.set_xlabel(r"matching cutoff $\nu$ (length$^{-1}$)")
    right.set_ylabel(r"$V_D$ (reference length$^{-4}$)")
    right.set_ylim(0, .205)
    right.legend(fontsize=6.5, loc="center right", frameon=False)
    for axis in (left, right):
        axis.grid(True, axis="y", alpha=.2, linewidth=.5)
        axis.set_axisbelow(True)
        axis.tick_params(labelsize=7)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
    figure.tight_layout(pad=.7)
    return figure


def parent_source_plot():
    """The curved canonical source and a separate flat regulator control."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    parent = json.loads((REPOSITORY / "results/development/unruh-state.json").read_text())
    regulator = json.loads((REPOSITORY / "results/development/state-regulator.json").read_text())
    neck = parent["neck_source_budget"]["candidate_tensor"]
    figure = Figure(figsize=(6.4, 3.0), dpi=220)
    FigureCanvasAgg(figure)
    left, right = figure.subplots(1, 2)
    fields = ("rho", "p_parallel", "p_sphere", "radial_null")
    values = [neck[key] for key in fields]
    left.bar(range(4), values, color=["#2357A6", "#21766C", "#21766C", "#8C4264"], width=.6)
    left.set_xticks(range(4), [r"$\rho$", r"$p_\parallel$", r"$p_\perp$", r"$\rho+p_\parallel$"])
    left.set_ylim(-.067, .006)
    left.set_ylabel(r"canonical source $\times L_{\rm throat}^{4}$")
    left.set_title("Curved neck: parent-state candidate", fontsize=8.5)
    for x, value in enumerate(values):
        left.text(x, value-.0022, f"{value:.4f}", ha="center", va="top", fontsize=6.5)
    selected = [row["images"] for row in regulator["thermal_source_rows"]
                if row["images"]["temperature"] in (.25, .5)]
    for j, row in enumerate(selected):
        right.bar(j-.17, row["canonical"]["rho"], width=.32, color="#2357A6",
                  label="Canonical state" if j==0 else None)
        right.bar(j+.17, row["finite_endpoint"]["rho"], width=.32, color="#A15A38",
                  label="Raw finite endpoint" if j==0 else None)
    right.set_xticks([0, 1], ["0.25", "0.50"])
    right.set_xlabel(r"thermal control $T/\nu$")
    right.set_ylabel(r"thermal density $\rho/\nu^4$")
    right.set_title("Flat equilibrium: matching is required", fontsize=8.5)
    right.set_ylim(-.01,.09)
    right.legend(fontsize=6.5, frameon=False, loc="upper left")
    for axis in (left, right):
        axis.axhline(0, color="#526779", linewidth=.7)
        axis.set_axisbelow(True)
        axis.grid(True, axis="y", alpha=.2, linewidth=.5)
        axis.tick_params(labelsize=7)
        axis.xaxis.label.set_size(8)
        axis.yaxis.label.set_size(8)
    figure.tight_layout(pad=.7)
    return figure


def build_pdf(
    source: Path, output: Path, metadata: dict[str, object]
) -> tuple[int, int, str]:
    try:
        import matplotlib
        import reportlab
        from matplotlib.font_manager import FontProperties
        from matplotlib.mathtext import math_to_image
        from reportlab.lib import colors
        from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import mm
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.platypus import (
            BaseDocTemplate,
            Frame,
            Image,
            KeepTogether,
            ListFlowable,
            ListItem,
            NextPageTemplate,
            PageBreak,
            PageTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
            XPreformatted,
        )
        from reportlab.platypus.tableofcontents import TableOfContents
    except ImportError as exc:  # pragma: no cover - exercised only without dependency
        raise SystemExit(
            "ReportLab is required. Run: python3 -m pip install -r "
            "requirements-paper.txt"
        ) from exc
    if getattr(reportlab, "Version", None) != "4.4.9":
        raise SystemExit(
            "ReportLab 4.4.9 is required for the canonical PDF; observed "
            f"{getattr(reportlab, 'Version', 'unknown')}"
        )

    page_width, page_height = A4
    left_margin = 19 * mm
    right_margin = 19 * mm
    top_margin = 21 * mm
    bottom_margin = 18 * mm
    body_width = page_width - left_margin - right_margin

    # Always use Matplotlib's pinned bundled DejaVu family.  It contains the
    # Greek glyphs used by the paper and is identical on macOS and Linux.
    font_directory = Path(matplotlib.get_data_path()) / "fonts" / "ttf"
    regular_font, bold_font, italic_font, bold_italic_font = (
        font_directory / "DejaVuSans.ttf",
        font_directory / "DejaVuSans-Bold.ttf",
        font_directory / "DejaVuSans-Oblique.ttf",
        font_directory / "DejaVuSans-BoldOblique.ttf",
    )
    if not all(
        path.is_file()
        for path in (regular_font, bold_font, italic_font, bold_italic_font)
    ):
        raise RuntimeError("Matplotlib's bundled DejaVu Sans fonts are absent")
    mono_font = font_directory / "DejaVuSansMono.ttf"
    if not mono_font.is_file():
        raise RuntimeError("Matplotlib's bundled DejaVu Sans Mono font is absent")
    pdfmetrics.registerFont(TTFont("RHRegular", regular_font))
    pdfmetrics.registerFont(TTFont("RHBold", bold_font))
    pdfmetrics.registerFont(TTFont("RHItalic", italic_font))
    pdfmetrics.registerFont(TTFont("RHBoldItalic", bold_italic_font))
    pdfmetrics.registerFont(TTFont("RHMono", mono_font))
    pdfmetrics.registerFontFamily(
        "RHRegular",
        normal="RHRegular",
        bold="RHBold",
        italic="RHItalic",
        boldItalic="RHBoldItalic",
    )

    ink = colors.HexColor("#122033")
    navy = colors.HexColor("#0A1730")
    blue = colors.HexColor("#2357A6")
    cyan = colors.HexColor("#64C9E8")
    muted = colors.HexColor("#536274")
    pale = colors.HexColor("#EFF4F8")
    line = colors.HexColor("#D3DCE6")
    white = colors.white

    sample = getSampleStyleSheet()
    styles: dict[str, ParagraphStyle] = {
        "Body": ParagraphStyle(
            "RHBody",
            parent=sample["BodyText"],
            fontName="RHRegular",
            fontSize=9.25,
            leading=13.35,
            textColor=ink,
            alignment=TA_JUSTIFY,
            spaceAfter=6.5,
            allowWidows=0,
            allowOrphans=0,
        ),
        "Heading1": ParagraphStyle(
            "RHHeading1",
            parent=sample["Heading1"],
            fontName="RHBold",
            fontSize=19,
            leading=23,
            textColor=navy,
            spaceBefore=14,
            spaceAfter=8,
            keepWithNext=True,
        ),
        "Heading2": ParagraphStyle(
            "RHHeading2",
            parent=sample["Heading2"],
            fontName="RHBold",
            fontSize=13.2,
            leading=16.5,
            textColor=blue,
            spaceBefore=10,
            spaceAfter=5,
            keepWithNext=True,
        ),
        "Heading3": ParagraphStyle(
            "RHHeading3",
            parent=sample["Heading3"],
            fontName="RHBold",
            fontSize=10.5,
            leading=13.5,
            textColor=ink,
            spaceBefore=8,
            spaceAfter=4,
            keepWithNext=True,
        ),
        "Bullet": ParagraphStyle(
            "RHBullet",
            parent=sample["BodyText"],
            fontName="RHRegular",
            fontSize=9.1,
            leading=12.8,
            textColor=ink,
            leftIndent=13,
            firstLineIndent=0,
            spaceAfter=3.5,
        ),
        "Numbered": ParagraphStyle(
            "RHNumbered",
            parent=sample["BodyText"],
            fontName="RHRegular",
            fontSize=8.9,
            leading=12.4,
            textColor=ink,
            leftIndent=30,
            firstLineIndent=0,
            bulletIndent=0,
            bulletFontName="RHBold",
            bulletFontSize=8.2,
            bulletColor=blue,
            spaceAfter=4.2,
        ),
        "Quote": ParagraphStyle(
            "RHQuote",
            parent=sample["BodyText"],
            fontName="RHItalic",
            fontSize=9.3,
            leading=13.4,
            textColor=ink,
            spaceAfter=0,
        ),
        "Code": ParagraphStyle(
            "RHCode",
            parent=sample["Code"],
            fontName="RHMono",
            fontSize=7.6,
            leading=10.2,
            textColor=colors.HexColor("#EAF2FA"),
            leftIndent=4,
            rightIndent=4,
            spaceAfter=0,
        ),
        "Table": ParagraphStyle(
            "RHTable",
            parent=sample["BodyText"],
            fontName="RHRegular",
            fontSize=7.4,
            leading=10.0,
            textColor=ink,
        ),
        "TableHead": ParagraphStyle(
            "RHTableHead",
            parent=sample["BodyText"],
            fontName="RHBold",
            fontSize=7.5,
            leading=10.1,
            textColor=white,
        ),
        "CoverTitle": ParagraphStyle(
            "RHCoverTitle",
            parent=sample["Title"],
            fontName="RHBold",
            fontSize=34,
            leading=39,
            textColor=white,
            alignment=TA_LEFT,
            spaceAfter=12,
        ),
        "CoverDeck": ParagraphStyle(
            "RHCoverDeck",
            parent=sample["Heading2"],
            fontName="RHRegular",
            fontSize=17.5,
            leading=23,
            textColor=cyan,
            alignment=TA_LEFT,
            spaceAfter=24,
        ),
        "CoverSubtitle": ParagraphStyle(
            "RHCoverSubtitle",
            parent=sample["BodyText"],
            fontName="RHItalic",
            fontSize=14,
            leading=19,
            textColor=white,
            alignment=TA_LEFT,
            spaceAfter=18,
        ),
        "CoverMeta": ParagraphStyle(
            "RHCoverMeta",
            parent=sample["BodyText"],
            fontName="RHRegular",
            fontSize=9.5,
            leading=14,
            textColor=colors.HexColor("#B8C8DB"),
            alignment=TA_LEFT,
        ),
        "TocTitle": ParagraphStyle(
            "RHTocTitle",
            parent=sample["Heading1"],
            fontName="RHBold",
            fontSize=24,
            leading=28,
            textColor=navy,
            spaceAfter=18,
        ),
    }

    token_pattern = re.compile(r"@@RH_TOKEN_(\d+)@@")

    greek = {
        "Pi": "Π",
        "Psi": "Ψ",
        "psi": "ψ",
        "xi": "ξ",
        "Xi": "Ξ",
        "Sigma": "Σ",
        "nu": "ν",
        "chi": "χ",
        "alpha": "α",
        "beta": "β",
        "gamma": "γ",
        "Gamma": "Γ",
        "delta": "δ",
        "Delta": "Δ",
        "epsilon": "ε",
        "varepsilon": "ε",
        "eta": "η",
        "theta": "θ",
        "Theta": "Θ",
        "kappa": "κ",
        "Lambda": "Λ",
        "lambda": "λ",
        "mu": "μ",
        "Omega": "Ω",
        "omega": "ω",
        "Phi": "Φ",
        "phi": "φ",
        "varphi": "φ",
        "pi": "π",
        "rho": "ρ",
        "sigma": "σ",
        "tau": "τ",
        "zeta": "ζ",
    }

    def latex_to_plain(value: str) -> str:
        """Make short inline formulae legible in ReportLab paragraphs."""

        value = value.replace(r"\!", "").replace(r"\,", " ")
        value = value.replace(r"\;", " ").replace(r"\quad", " ")
        value = value.replace(r"\qquad", "  ")
        value = value.replace(r"\dagger", "†").replace(r"\star", "⋆")
        value = value.replace(r"\to", "→").replace(r"\mapsto", "↦")
        value = re.sub(r"\\(?:geq|ge)(?![A-Za-z])", "≥", value)
        value = re.sub(r"\\(?:leq|le)(?![A-Za-z])", "≤", value)
        value = value.replace(r"\neq", "≠").replace(r"\approx", "≈")
        value = re.sub(r"\\in(?![A-Za-z])", "∈", value).replace(r"\times", "×")
        value = re.sub(r"\\int(?![A-Za-z])", "∫", value)
        value = re.sub(r"\\oint(?![A-Za-z])", "∮", value)
        value = value.replace(r"\ell", "ℓ").replace(r"\odot", "⊙")
        value = value.replace(r"\ddots", "⋱").replace(r"\cdots", "⋯")
        value = value.replace(r"\otimes", "⊗").replace(r"\cdot", "·")
        value = value.replace(r"\partial", "∂").replace(r"\infty", "∞")
        value = value.replace(r"\hbar", "ℏ").replace(r"\pm", "±").replace(r"\mp", "∓")
        value = value.replace(r"\Box", "□")
        value = value.replace(r"\parallel", "∥").replace(r"\perp", "⊥")
        value = value.replace(r"\bar\psi", "ψ̄")
        value = re.sub(r"\\bar\s+([A-Za-z])", lambda match: match.group(1) + "̄", value)
        value = value.replace(r"\sqrt", "√")
        value = re.sub(
            r"\\frac\{([^{}]+)\}\{([^{}]+)\}", r"(\1)/(\2)", value
        )
        for command in ("mathrm", "mathbf", "mathbb", "mathcal", "operatorname", "text"):
            value = re.sub(rf"\\{command}\{{([^{{}}]+)\}}", r"\1", value)
            value = re.sub(rf"\\{command}\s+([A-Za-z])", r"\1", value)
        for name, symbol in greek.items():
            value = value.replace(f"\\{name}", symbol)
        value = re.sub(r"\\(log|ln|exp|cos|sin|tan|det|Tr)(?![A-Za-z])", r"\1", value)
        value = re.sub(r"_\{([^{}]+)\}", r"_(\1)", value)
        value = re.sub(r"\^\{([^{}]+)\}", r"^(\1)", value)
        value = value.replace("{", "").replace("}", "")
        value = value.replace(r"\left", "").replace(r"\right", "")
        value = value.replace(r"\langle", "⟨").replace(r"\rangle", "⟩")
        return re.sub(r"\s+", " ", value).strip()

    def inline_math_markup(value: str) -> str:
        plain = html.escape(latex_to_plain(value), quote=False)
        def script(match: re.Match[str]) -> str:
            tag = "sub" if match.group(1) == "_" else "super"
            return f"<{tag}>{match.group(2)}</{tag}>"
        plain = re.sub(r"([_^])\(([^()]*)\)", script, plain)
        return re.sub(r"([_^])([\w*⋆†±+∥⊥-])", script, plain)

    def inline_markup(value: str) -> str:
        tokens: list[str] = []

        def reserve(markup: str) -> str:
            tokens.append(markup)
            return f"@@RH_TOKEN_{len(tokens) - 1}@@"

        def link_replacement(match: re.Match[str]) -> str:
            label = html.escape(match.group(1).strip("`"), quote=False)
            url = html.escape(match.group(2), quote=True)
            return reserve(f'<link href="{url}" color="#2357A6"><u>{label}</u></link>')

        value = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link_replacement, value)
        value = re.sub(
            r"`([^`]+)`",
            lambda match: reserve(
                '<font name="Courier" color="#173D6B">'
                + html.escape(match.group(1), quote=False)
                + "</font>"
            ),
            value,
        )
        value = re.sub(
            r"\\\((.+?)\\\)",
            lambda match: reserve(
                '<font name="RHItalic">'
                + inline_math_markup(match.group(1))
                + "</font>"
            ),
            value,
        )
        # GitHub Markdown uses single-dollar delimiters for inline mathematics.
        # Keep support for the legacy parenthesized form above so historical
        # source notes can still be rendered, but make the canonical manuscript
        # portable between GitHub and this PDF builder.
        value = re.sub(
            r"(?<!\$)\$(?!\$)(.+?)(?<!\$)\$(?!\$)",
            lambda match: reserve(
                '<font name="RHItalic">'
                + inline_math_markup(match.group(1))
                + "</font>"
            ),
            value,
        )
        value = html.escape(value, quote=False)
        value = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", value)
        value = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<i>\1</i>", value)
        value = value.replace("  ", " ")
        return token_pattern.sub(lambda match: tokens[int(match.group(1))], value)

    def table_flowable(lines: list[str]) -> Table:
        raw_rows = [[cell.strip() for cell in row.strip().strip("|").split("|")] for row in lines]
        if len(raw_rows) > 1 and all(
            re.fullmatch(r":?-{3,}:?", cell.replace(" ", ""))
            for cell in raw_rows[1]
        ):
            raw_rows.pop(1)
        columns = max(len(row) for row in raw_rows)
        lengths = []
        for column in range(columns):
            longest = max(
                (len(row[column]) if column < len(row) else 0) for row in raw_rows
            )
            lengths.append(max(8, min(longest, 42)))
        total = sum(lengths)
        widths = [body_width * length / total for length in lengths]
        if columns >= 3:
            widths = [min(width, body_width * 0.46) for width in widths]
            scale = body_width / sum(widths)
            widths = [width * scale for width in widths]

        rows = []
        for row_index, row in enumerate(raw_rows):
            style = styles["TableHead"] if row_index == 0 else styles["Table"]
            padded = row + [""] * (columns - len(row))
            rows.append([Paragraph(inline_markup(cell), style) for cell in padded])

        table = Table(rows, colWidths=widths, repeatRows=1, hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), blue),
                    ("BACKGROUND", (0, 1), (-1, -1), colors.white),
                    ("BOX", (0, 0), (-1, -1), 0.45, line),
                    ("INNERGRID", (0, 0), (-1, -1), 0.3, line),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 5),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                    ("TOPPADDING", (0, 0), (-1, -1), 4),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, pale]),
                ]
            )
        )
        return table

    def quote_flowable(lines: list[str]) -> Table:
        # A Table can split between rows but ReportLab cannot split one tall
        # Paragraph inside a single table row.  Keep short source blockquotes
        # visually continuous while giving long status/audit notes safe page
        # break points at the manuscript's own wrapped-line boundaries.
        source_lines = [line.removeprefix(">").strip() for line in lines]
        chunks = [
            " ".join(source_lines[index : index + 4])
            for index in range(0, len(source_lines), 4)
        ]
        quote = Table(
            [[Paragraph(inline_markup(chunk), styles["Quote"])] for chunk in chunks],
            colWidths=[body_width - 6],
            splitByRow=1,
        )
        quote.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), pale),
                    ("LINEBEFORE", (0, 0), (0, -1), 3, cyan),
                    ("LEFTPADDING", (0, 0), (-1, -1), 10),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                    ("TOPPADDING", (0, 0), (-1, -1), 0),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
                    ("TOPPADDING", (0, 0), (-1, 0), 8),
                    ("BOTTOMPADDING", (0, -1), (-1, -1), 8),
                ]
            )
        )
        return quote

    def code_flowable(code: str) -> Table:
        block = XPreformatted(html.escape(code.rstrip(), quote=False), styles["Code"])
        table = Table([[block]], colWidths=[body_width], hAlign="LEFT")
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), navy),
                    ("BOX", (0, 0), (-1, -1), 0.4, colors.HexColor("#20395B")),
                    ("LEFTPADDING", (0, 0), (-1, -1), 7),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 7),
                    ("TOPPADDING", (0, 0), (-1, -1), 6),
                    ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
                ]
            )
        )
        return table

    def matrix_plaintext(source: str) -> str:
        before, _, remainder = source.partition(r"\begin{pmatrix}")
        body, _, _ = remainder.partition(r"\end{pmatrix}")
        heading = latex_to_plain(before).strip()
        rows = []
        for row in re.split(r"\\\\", body):
            cells = [latex_to_plain(cell.strip()) or " " for cell in row.split("&")]
            if any(cell.strip() for cell in cells):
                rows.append(cells)
        widths = [
            max(len(row[index]) if index < len(row) else 0 for row in rows)
            for index in range(max(map(len, rows), default=0))
        ]
        rendered = []
        for row in rows:
            padded = row + [""] * (len(widths) - len(row))
            rendered.append("[ " + "   ".join(value.ljust(widths[index]) for index, value in enumerate(padded)) + " ]")
        return "\n".join(([heading] if heading else []) + rendered)

    def equation_flowable(source: str) -> object:
        expression = " ".join(line.strip() for line in source.splitlines()).strip()
        if r"\begin{pmatrix}" in expression:
            return code_flowable(matrix_plaintext(source))
        if expression.startswith(r"\boxed{") and expression.endswith("}"):
            expression = expression[len(r"\boxed{") : -1].strip()
        expression = re.sub(
            r"\\(mathbb|mathcal|mathrm|mathbf)\s+([A-Za-z])",
            r"\\\1{\2}",
            expression,
        )
        segments = [
            segment.strip().rstrip(",")
            for segment in re.split(r"\\qquad", expression)
            if segment.strip().rstrip(",")
        ]
        rows = []
        for segment in segments:
            buffer = BytesIO()
            math_to_image(
                f"${segment}$",
                buffer,
                prop=FontProperties(size=13),
                dpi=220,
                format="png",
                color="#122033",
            )
            buffer.seek(0)
            rendered = Image(buffer)
            rendered._nsc_buffer = buffer
            rendered.drawWidth *= 72.0 / 220.0
            rendered.drawHeight *= 72.0 / 220.0
            maximum = body_width - 28
            if rendered.drawWidth > maximum:
                scale = maximum / rendered.drawWidth
                rendered.drawWidth *= scale
                rendered.drawHeight *= scale
            rendered.hAlign = "CENTER"
            rows.append([rendered])
        table = Table(rows, colWidths=[body_width], hAlign="CENTER")
        table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, -1), white),
                    ("BOX", (0, 0), (-1, -1), 0.45, line),
                    ("LEFTPADDING", (0, 0), (-1, -1), 10),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                    ("TOPPADDING", (0, 0), (-1, 0), 8),
                    ("BOTTOMPADDING", (0, -1), (-1, -1), 8),
                ]
            )
        )
        return KeepTogether([Spacer(1, 3), table, Spacer(1, 7)])

    def geometric_gap_figure() -> object:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        record = json.loads(
            (REPOSITORY / "results" / "nsc-3-geometric-chain.json").read_text(
                encoding="utf-8"
            )
        )
        families = record["gap_families"]
        radii = [row["radius"] for row in families]
        edges = [row["continuum"]["band_edge"] for row in families]
        figure = Figure(figsize=(6.4, 2.6), dpi=220)
        FigureCanvasAgg(figure)
        left, right = figure.subplots(1, 2)
        left.plot(radii, edges, "o-", color="#2357A6", markersize=5, linewidth=1.2)
        left.set_xlabel("motif radius R")
        left.set_ylabel("continuum first band edge")
        left.set_title("Gap versus declared R")
        left.grid(True, alpha=0.3, linewidth=0.4)
        for family, colour in zip(families, ("#2357A6", "#0A1730", "#64C9E8")):
            intervals = [row["intervals"] for row in family["lattice"]]
            errors = family["absolute_errors"]
            right.loglog(
                intervals,
                errors,
                "o-",
                color=colour,
                markersize=4,
                linewidth=1.1,
                label=f"R={family['radius']:g}",
            )
        right.set_xlabel("intervals per motif")
        right.set_xticks([32, 64, 128], labels=["32", "64", "128"])
        right.minorticks_off()
        right.set_ylabel("absolute lattice error")
        right.set_title("Second-order lattice agreement")
        right.grid(True, which="both", alpha=0.3, linewidth=0.4)
        right.legend(frameon=False, fontsize=7)
        figure.tight_layout(pad=0.4)
        buffer = BytesIO()
        figure.savefig(buffer, format="png", dpi=220)
        buffer.seek(0)
        rendered = Image(buffer)
        rendered._nsc_buffer = buffer
        rendered.drawWidth = body_width
        rendered.drawHeight = body_width * 2.6 / 6.4
        caption = Paragraph(
            inline_markup(
                "**Figure.** Continuum first band edges and lattice errors from "
                "[`nsc-3-geometric-chain.json`](../results/nsc-3-geometric-chain.json). "
                "R is a declared input."
            ),
            styles["Body"],
        )
        return KeepTogether([Spacer(1, 4), rendered, caption, Spacer(1, 6)])

    def vacuum_work_figure() -> object:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure

        record = json.loads((REPOSITORY / "results/nsc-6-vacuum-work.json").read_text())
        figure = Figure(figsize=(6.4, 2.8), dpi=220)
        FigureCanvasAgg(figure)
        left, right = figure.subplots(1, 2)
        grids = record["spatial_refinement"]
        left.plot([r["points"] for r in grids],
                  [r["full_angular_through_kappa8_energy_coefficient"] for r in grids],
                  "o-", color="#2357A6", linewidth=1.2)
        left.set_xlabel("axial grid points")
        left.set_ylabel("energy / pulse amplitude squared")
        left.set_title("Angular sum through kappa = 8")
        left.grid(True, alpha=.25)
        rows = [r for r in record["finite_real_time_controls"] if r["epsilon"] > 0]
        x = [r["epsilon"]**2 for r in rows]
        right.plot(x, [r["external_work"] for r in rows], "o", color="#2357A6", label="supplied work")
        right.plot(x, [r["excitation_energy"] for r in rows], "+", color="#0A1730", label="excitation energy")
        coefficient = grids[0]["response"]["energy_coefficient"]
        right.plot([0, max(x)], [0, coefficient*max(x)], "--", color="#64C9E8", label="leading spectrum")
        right.set_xlabel("pulse amplitude squared")
        right.set_ylabel("energy in throat units")
        right.set_title("One channel; 32 axial points")
        right.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
        right.grid(True, alpha=.25)
        right.legend(frameon=False, fontsize=6.5)
        figure.tight_layout(pad=.6)
        buffer = BytesIO()
        figure.savefig(buffer, format="png", dpi=220)
        buffer.seek(0)
        rendered = Image(buffer)
        rendered._nsc_buffer = buffer
        rendered.drawWidth = body_width
        rendered.drawHeight = body_width*2.8/6.4
        caption = Paragraph(inline_markup(
            "**Figure.** Response to a prescribed geometry pulse. Left: spatial refinement of the angular sum. "
            "Right: independent finite evolution accounts for supplied work in one channel. "
            "The pulse is an input; this is not a self-sourced cosmology. "
            "[Vacuum-work evidence](../results/nsc-6-vacuum-work.json)."), styles["Body"])
        return KeepTogether([Spacer(1,4),rendered,caption,Spacer(1,6)])

    def compact_source_figure() -> object:
        buffer = BytesIO()
        figure = compact_source_plot()
        figure.savefig(buffer, format="png", dpi=220)
        buffer.seek(0)
        rendered = Image(buffer)
        rendered._nsc_buffer = buffer
        rendered.drawWidth = body_width
        rendered.drawHeight = body_width*2.8/6.4
        caption = Paragraph(inline_markup(
            "**Figure.** Left: full free-KK holonomy energy differences on the ultrastatic "
            "closed cell; the effective AP phase minimizes this one-loop real potential. "
            "Right: only the compact endpoint interaction's effective null source, at that phase. "
            "The transverse separation is 2 in g4 neck units. Increasing the axial circle "
            "changes the global geometry and the source sign. The remaining bulk source and "
            "transmitting geometry are still required. "
            "[Source and state record](../results/development/compact-casimir.json)."), styles["Body"])
        return KeepTogether([Spacer(1,4),rendered,caption,Spacer(1,6)])

    def source_matching_figure() -> object:
        buffer = BytesIO()
        source_matching_plot().savefig(buffer, format="png", dpi=220)
        buffer.seek(0)
        rendered = Image(buffer)
        rendered._nsc_buffer = buffer
        rendered.drawWidth = body_width
        rendered.drawHeight = body_width*2.9/6.4
        caption = Paragraph(inline_markup(
            "**Figure.** Left: two-dimensional-area potential on R2 times the unit sphere, "
            "with magnetic flux label 1. Right: the four-dimensional local vacuum coefficient "
            "split between a retained light field and its complement; their sum is unchanged. "
            "Both use cutoff 2 and transverse length 2. The panels have different dimensions "
            "and are different projections of the same free determinant. They do not plot "
            "cosmological evolution. "
            "[Warped source](../results/development/warped-source.json); "
            "[light-field matching](../results/development/compact-matching.json)."), styles["Body"])
        return KeepTogether([Spacer(1,4), rendered, caption, Spacer(1,6)])

    def parent_source_figure() -> object:
        buffer = BytesIO()
        parent_source_plot().savefig(buffer, format="png", dpi=220)
        buffer.seek(0)
        rendered = Image(buffer)
        rendered._nsc_buffer = buffer
        rendered.drawWidth = body_width
        rendered.drawHeight = body_width*3.0/6.4
        caption = Paragraph(inline_markup(
            "**Figure.** Left: the canonical massless Dirac source at the imposed curved neck. "
            "The last bar is the mean of the two radial null contractions; nonzero flux "
            "splits them slightly, with both remaining negative. Density, anisotropy and flux "
            "leave independent geometric residuals. Right: a separate flat thermal control "
            "shows the finite-endpoint/state conversion requirement. The panels use distinct "
            "geometries and normalizations. They do not depict cosmic evolution or assign a "
            "temperature to the child. [Parent source](../results/development/unruh-state.json); "
            "[state/regulator control](../results/development/state-regulator.json)."), styles["Body"])
        return KeepTogether([Spacer(1,4), rendered, caption, Spacer(1,6)])

    def parse_markdown(markdown: str) -> list[object]:
        lines = markdown.splitlines()
        # Cover content is built separately. Begin at the technical abstract.
        start = next(
            (index for index, line_value in enumerate(lines) if line_value == "## Abstract"),
            0,
        )
        lines = lines[start:]
        story: list[object] = []
        index = 0
        paragraph: list[str] = []

        def flush_paragraph() -> None:
            if paragraph:
                story.append(
                    Paragraph(inline_markup(" ".join(part.strip() for part in paragraph)), styles["Body"])
                )
                paragraph.clear()

        while index < len(lines):
            line_value = lines[index]
            stripped = line_value.strip()
            if not stripped:
                flush_paragraph()
                index += 1
                continue

            if stripped == "<!-- nsc-figure:geometric-gap -->":
                flush_paragraph()
                story.append(geometric_gap_figure())
                index += 1
                continue

            if stripped == "<!-- nsc-figure:vacuum-work -->":
                flush_paragraph()
                story.append(vacuum_work_figure())
                index += 1
                continue

            if stripped == "<!-- nsc-figure:compact-source -->":
                flush_paragraph()
                story.append(compact_source_figure())
                index += 1
                continue

            if stripped == "<!-- nsc-figure:source-matching -->":
                flush_paragraph()
                story.append(source_matching_figure())
                index += 1
                continue

            if stripped == "<!-- nsc-figure:parent-source -->":
                flush_paragraph()
                story.append(parent_source_figure())
                index += 1
                continue

            if stripped.startswith("<!--") and stripped.endswith("-->"):
                flush_paragraph()
                index += 1
                continue

            if stripped.startswith("```"):
                flush_paragraph()
                index += 1
                code_lines: list[str] = []
                while index < len(lines) and not lines[index].strip().startswith("```"):
                    code_lines.append(lines[index])
                    index += 1
                index += 1
                story.extend([code_flowable("\n".join(code_lines)), Spacer(1, 6)])
                continue

            if stripped == "$$":
                flush_paragraph()
                index += 1
                equation_lines: list[str] = []
                while index < len(lines) and lines[index].strip() != "$$":
                    equation_lines.append(lines[index])
                    index += 1
                if index >= len(lines):
                    raise RuntimeError("unclosed display-math block in manuscript")
                index += 1
                story.append(equation_flowable("\n".join(equation_lines)))
                continue

            if stripped.startswith("|"):
                flush_paragraph()
                table_lines: list[str] = []
                while index < len(lines) and lines[index].strip().startswith("|"):
                    table_lines.append(lines[index])
                    index += 1
                table = table_flowable(table_lines)
                if len(table_lines) <= 6:
                    story.append(KeepTogether([table, Spacer(1, 8)]))
                else:
                    story.extend([table, Spacer(1, 8)])
                continue

            if stripped.startswith(">"):
                flush_paragraph()
                quote_lines: list[str] = []
                while index < len(lines) and lines[index].strip().startswith(">"):
                    quote_lines.append(lines[index].strip())
                    index += 1
                story.extend([quote_flowable(quote_lines), Spacer(1, 8)])
                continue

            heading = re.match(r"^(#{1,3})\s+(.+)$", stripped)
            if heading:
                flush_paragraph()
                level = len(heading.group(1))
                story.append(
                    Paragraph(inline_markup(heading.group(2)), styles[f"Heading{level}"])
                )
                index += 1
                continue

            if re.match(r"^[-*]\s+", stripped):
                flush_paragraph()
                items: list[ListItem] = []
                while index < len(lines) and re.match(r"^\s*[-*]\s+", lines[index]):
                    item_text = re.sub(r"^\s*[-*]\s+", "", lines[index]).strip()
                    items.append(ListItem(Paragraph(inline_markup(item_text), styles["Bullet"])))
                    index += 1
                story.append(
                    ListFlowable(
                        items,
                        bulletType="bullet",
                        bulletFontName="RHRegular",
                        bulletFontSize=7,
                        leftIndent=15,
                        bulletColor=blue,
                        spaceAfter=5,
                    )
                )
                continue

            numbered = re.match(r"^(\d+)\.\s+(.+)$", stripped)
            if numbered:
                flush_paragraph()
                while index < len(lines):
                    match = re.match(r"^\s*(\d+)\.\s+(.+)$", lines[index])
                    if not match:
                        break
                    story.append(
                        Paragraph(
                            inline_markup(match.group(2)),
                            styles["Numbered"],
                            bulletText=f"{match.group(1)}.",
                        )
                    )
                    index += 1
                continue

            paragraph.append(line_value)
            index += 1

        flush_paragraph()
        return story

    class NestedSpaceDocTemplate(BaseDocTemplate):
        def afterFlowable(self, flowable: object) -> None:  # noqa: N802 - ReportLab API
            if isinstance(flowable, Paragraph):
                if flowable.style.name == "RHHeading1":
                    self.notify("TOCEntry", (0, flowable.getPlainText(), self.page))
                elif flowable.style.name == "RHHeading2":
                    self.notify("TOCEntry", (1, flowable.getPlainText(), self.page))

    def draw_cover(canvas: object, document: object) -> None:
        canvas.saveState()
        canvas.setTitle(str(metadata["title"]))
        canvas.setAuthor(str(metadata["author"]))
        canvas.setSubject(
            "One inherited spectrum connecting matter, geometry, dark response, and nested space"
        )
        canvas.setKeywords(
            "nested-space cosmology, frequency spectrum, black holes, antimatter, dark matter, dark energy, spectral geometry"
        )
        canvas.setFillColor(navy)
        canvas.rect(0, 0, page_width, page_height, fill=1, stroke=0)
        canvas.setStrokeColor(colors.HexColor("#16345D"))
        canvas.setLineWidth(1.0)
        for radius in (58, 95, 142, 205, 282):
            canvas.circle(page_width * 0.89, page_height * 0.18, radius, fill=0, stroke=1)
        canvas.setStrokeColor(cyan)
        canvas.setLineWidth(2.0)
        canvas.line(left_margin, page_height - 39 * mm, left_margin + 34 * mm, page_height - 39 * mm)
        canvas.restoreState()

    def draw_body(canvas: object, document: object) -> None:
        canvas.saveState()
        canvas.setFont("RHRegular", 6.8)
        canvas.setFillColor(muted)
        canvas.drawString(left_margin, page_height - 11 * mm, str(metadata["project"]).upper())
        right = f"WORKING PREPRINT  ·  VERSION {metadata['version']}"
        canvas.drawRightString(page_width - right_margin, page_height - 11 * mm, right)
        canvas.setStrokeColor(line)
        canvas.setLineWidth(0.5)
        canvas.line(left_margin, page_height - 13 * mm, page_width - right_margin, page_height - 13 * mm)
        canvas.line(left_margin, 12 * mm, page_width - right_margin, 12 * mm)
        canvas.setFont("RHRegular", 7.2)
        canvas.drawString(left_margin, 8.2 * mm, "Douglas Ek · Equations and evidence")
        canvas.drawRightString(page_width - right_margin, 8.2 * mm, str(document.page))
        canvas.restoreState()

    cover_frame = Frame(
        left_margin,
        25 * mm,
        body_width,
        page_height - 50 * mm,
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
        id="cover-frame",
    )
    body_frame = Frame(
        left_margin,
        bottom_margin,
        body_width,
        page_height - top_margin - bottom_margin,
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
        id="body-frame",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    BUILD_DIRECTORY.mkdir(parents=True, exist_ok=True)
    temporary = BUILD_DIRECTORY / f"{output.name}.tmp"
    if temporary.exists():
        temporary.unlink()

    os.environ["SOURCE_DATE_EPOCH"] = str(metadata["source_date_epoch"])
    os.environ["TZ"] = "UTC"
    if hasattr(time, "tzset"):
        time.tzset()
    document = NestedSpaceDocTemplate(
        str(temporary),
        pagesize=A4,
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=top_margin,
        bottomMargin=bottom_margin,
        title=str(metadata["title"]),
        author=str(metadata["author"]),
        invariant=1,
    )
    document.addPageTemplates(
        [
            PageTemplate(id="Cover", frames=[cover_frame], onPage=draw_cover),
            PageTemplate(id="Body", frames=[body_frame], onPage=draw_body),
        ]
    )

    toc = TableOfContents()
    toc.levelStyles = [
        ParagraphStyle(
            "RHTocLevel1",
            fontName="RHBold",
            fontSize=9.2,
            leading=13,
            textColor=ink,
            leftIndent=0,
            firstLineIndent=0,
            spaceBefore=3,
        ),
        ParagraphStyle(
            "RHTocLevel2",
            fontName="RHRegular",
            fontSize=8.1,
            leading=11.5,
            textColor=muted,
            leftIndent=14,
            firstLineIndent=0,
            spaceBefore=1,
        ),
    ]

    cover_story: list[object] = [
        Spacer(1, 31 * mm),
        Paragraph("NESTED-SPACE<br/>COSMOLOGY", styles["CoverTitle"]),
        Paragraph(
            str(metadata["subtitle"]),
            styles["CoverDeck"],
        ),
        Paragraph(
            "What if one inherited spectrum is enough?",
            styles["CoverSubtitle"],
        ),
        Spacer(1, 18 * mm),
        Paragraph(
            f"AUTHOR<br/>{str(metadata['author']).upper()}<br/><br/>"
            f"{str(metadata['project']).upper()}<br/>"
            f"VERSION {metadata['version']} · {str(metadata['date']).upper()}<br/><br/>"
            f"STATUS: {metadata['status']}",
            styles["CoverMeta"],
        ),
        NextPageTemplate("Body"),
        PageBreak(),
        Paragraph("Contents", styles["TocTitle"]),
        toc,
        PageBreak(),
    ]
    manuscript = source.read_text(encoding="utf-8")
    story = cover_story + parse_markdown(manuscript)
    document.multiBuild(story)

    if not temporary.is_file() or temporary.stat().st_size < 20_000:
        raise RuntimeError("PDF build did not produce a plausible output file")
    with temporary.open("rb") as handle:
        if handle.read(5) != b"%PDF-":
            raise RuntimeError("generated file does not have a PDF signature")
    temporary.replace(output)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    return int(document.page), output.stat().st_size, digest


def write_build_manifest(
    source: Path,
    output: Path,
    metadata: dict[str, object],
    pages: int,
    size: int,
    digest: str,
) -> None:
    record = {
        "schema": "NSC-PAPER-BUILD-v1",
        "source": str(source.relative_to(REPOSITORY)),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "metadata": str(METADATA.relative_to(REPOSITORY)),
        "metadata_sha256": hashlib.sha256(METADATA.read_bytes()).hexdigest(),
        "output": str(output.relative_to(REPOSITORY)),
        "pdf_sha256": digest,
        "pages": pages,
        "bytes": size,
        "version": metadata["version"],
        "scientific_build_inputs": scientific_build_inputs(source),
    }
    BUILD_MANIFEST.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def check_build(source: Path, output: Path, metadata: dict[str, object]) -> int:
    if not output.is_file() or not BUILD_MANIFEST.is_file():
        print("tracked PDF or build manifest is absent", file=sys.stderr)
        return 1
    with tempfile.TemporaryDirectory(prefix="nsc-paper-check-") as directory:
        first = Path(directory) / "first.pdf"
        second = Path(directory) / "second.pdf"
        first_result = build_pdf(source, first, metadata)
        second_result = build_pdf(source, second, metadata)
        if first.read_bytes() != second.read_bytes():
            print("two identical paper builds produced different bytes", file=sys.stderr)
            return 1
        if first.read_bytes() != output.read_bytes():
            print("tracked PDF does not match a fresh deterministic build", file=sys.stderr)
            return 1
        manifest = json.loads(BUILD_MANIFEST.read_text(encoding="utf-8"))
        expected = {
            "scientific_build_inputs": scientific_build_inputs(source),
            "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "metadata_sha256": hashlib.sha256(METADATA.read_bytes()).hexdigest(),
            "pdf_sha256": first_result[2],
            "pages": first_result[0],
            "bytes": first_result[1],
        }
        for key, value in expected.items():
            if manifest.get(key) != value:
                print(f"paper manifest mismatch for {key}", file=sys.stderr)
                return 1
        if first_result != second_result:
            print("paper build metadata is nondeterministic", file=sys.stderr)
            return 1
    print(f"paper check passed: {output} ({expected['pdf_sha256']})")
    return 0


def main() -> int:
    args = parse_args()
    output = args.output.resolve()
    if args.clean:
        return clean_output(output)
    source = args.source.resolve()
    if not source.is_file():
        print(f"source manuscript is absent: {source}", file=sys.stderr)
        return 1
    metadata = load_metadata()
    if args.check:
        return check_build(source, output, metadata)
    pages, size, digest = build_pdf(source, output, metadata)
    write_build_manifest(source, output, metadata, pages, size, digest)
    print(f"wrote {output}")
    print(f"pages: {pages}")
    print(f"bytes: {size}")
    print(f"sha256: {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
