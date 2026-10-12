#!/usr/bin/env python3
"""Render the canonical Markdown manuscript as a publication-form PDF."""

from __future__ import annotations

import argparse
import hashlib
import html
import os
import re
import shutil
import sys
import time
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE = REPOSITORY / "paper" / "recursive-horizons.md"
OUTPUT = REPOSITORY / "paper" / "recursive-horizons.pdf"
BUILD_DIRECTORY = REPOSITORY / "paper" / ".build"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--clean",
        action="store_true",
        help="move the generated PDF to a recoverable local build archive",
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


def build_pdf(source: Path, output: Path) -> tuple[int, int, str]:
    try:
        import reportlab
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

    # Prefer a host font with complete Greek coverage.  Linux GitHub runners
    # normally provide DejaVu Sans; macOS provides Arial.  ReportLab's bundled
    # Bitstream Vera remains a final fallback for environments without either.
    font_directory = Path(reportlab.__file__).resolve().parent / "fonts"
    font_candidates = [
        (
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf"),
            Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf"),
        ),
        (
            Path("/System/Library/Fonts/Supplemental/Arial.ttf"),
            Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf"),
            Path("/System/Library/Fonts/Supplemental/Arial Italic.ttf"),
            Path("/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf"),
        ),
        (
            font_directory / "Vera.ttf",
            font_directory / "VeraBd.ttf",
            font_directory / "VeraIt.ttf",
            font_directory / "VeraBI.ttf",
        ),
    ]
    regular_font, bold_font, italic_font, bold_italic_font = next(
        candidates for candidates in font_candidates if all(path.is_file() for path in candidates)
    )
    pdfmetrics.registerFont(TTFont("RHRegular", regular_font))
    pdfmetrics.registerFont(TTFont("RHBold", bold_font))
    pdfmetrics.registerFont(TTFont("RHItalic", italic_font))
    pdfmetrics.registerFont(TTFont("RHBoldItalic", bold_italic_font))
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
            fontName="Courier",
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

    def inline_markup(value: str) -> str:
        tokens: list[str] = []

        def reserve(markup: str) -> str:
            tokens.append(markup)
            return f"@@RH_TOKEN_{len(tokens) - 1}@@"

        def link_replacement(match: re.Match[str]) -> str:
            label = html.escape(match.group(1), quote=False)
            url = html.escape(match.group(2), quote=True)
            return reserve(f'<link href="{url}" color="#2357A6"><u>{label}</u></link>')

        value = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", link_replacement, value)
        value = re.sub(
            r"`([^`]+)`",
            lambda match: reserve(
                '<font name="Courier" color="#173D6B">'
                + html.escape(match.group(1), quote=False)
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

    def parse_markdown(markdown: str) -> list[object]:
        lines = markdown.splitlines()
        # Cover content is built separately. Begin at the status statement.
        start = next(
            (index for index, line_value in enumerate(lines) if line_value.startswith("> **Status statement.")),
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

            if stripped.startswith("|"):
                flush_paragraph()
                table_lines: list[str] = []
                while index < len(lines) and lines[index].strip().startswith("|"):
                    table_lines.append(lines[index])
                    index += 1
                story.append(KeepTogether([table_flowable(table_lines), Spacer(1, 8)]))
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
                if heading.group(2) == "References":
                    story.append(PageBreak())
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

    class RecursiveHorizonsDocTemplate(BaseDocTemplate):
        def afterFlowable(self, flowable: object) -> None:  # noqa: N802 - ReportLab API
            if isinstance(flowable, Paragraph):
                if flowable.style.name == "RHHeading1":
                    self.notify("TOCEntry", (0, flowable.getPlainText(), self.page))
                elif flowable.style.name == "RHHeading2":
                    self.notify("TOCEntry", (1, flowable.getPlainText(), self.page))

    def draw_cover(canvas: object, document: object) -> None:
        canvas.saveState()
        canvas.setTitle("Recursive Horizons")
        canvas.setAuthor("Douglas Ek & ChatGPT 5.6 Sol")
        canvas.setSubject(
            "A conjectural thermodynamic research programme for black-hole-to-cosmology transitions"
        )
        canvas.setKeywords(
            "black holes, cosmology, Finite Gradient Closure, Nested Gradient Spectrum, horizon thermodynamics, Einstein-Cartan, falsifiability"
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
        canvas.drawString(left_margin, page_height - 11 * mm, "RECURSIVE HORIZONS")
        right = "RESEARCH PROGRAMME  ·  VERSION 0.11.0"
        canvas.drawRightString(page_width - right_margin, page_height - 11 * mm, right)
        canvas.setStrokeColor(line)
        canvas.setLineWidth(0.5)
        canvas.line(left_margin, page_height - 13 * mm, page_width - right_margin, page_height - 13 * mm)
        canvas.line(left_margin, 12 * mm, page_width - right_margin, 12 * mm)
        canvas.setFont("RHRegular", 7.2)
        canvas.drawString(left_margin, 8.2 * mm, "Conjectural · not peer reviewed")
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

    # Fix the PDF timestamps at the declared release date so identical inputs
    # produce identical bytes on repeated builds with the same font set.
    os.environ["SOURCE_DATE_EPOCH"] = "1787270400"  # 2026-08-21 00:00:00 UTC
    document = RecursiveHorizonsDocTemplate(
        str(temporary),
        pagesize=A4,
        leftMargin=left_margin,
        rightMargin=right_margin,
        topMargin=top_margin,
        bottomMargin=bottom_margin,
        title="Recursive Horizons",
        author="Douglas Ek & ChatGPT 5.6 Sol",
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
        Paragraph("RECURSIVE<br/>HORIZONS", styles["CoverTitle"]),
        Paragraph(
            "A Thermodynamic Research Programme for<br/>Black-Hole-to-Cosmology Transitions",
            styles["CoverDeck"],
        ),
        Paragraph(
            "Can Einstein and Hawking Finally Rest Easy Together? Here Is What Must Be Proven",
            styles["CoverSubtitle"],
        ),
        Spacer(1, 18 * mm),
        Paragraph(
            "AUTHORED AND PRODUCED BY<br/>DOUGLAS EK &amp; CHATGPT 5.6 SOL<br/><br/>"
            "BLACKHOLES-INFINITY PROJECT<br/>VERSION 0.11.0 · 21 AUGUST 2026<br/><br/>"
            "STATUS: CONJECTURAL RESEARCH PROGRAMME · NOT PEER REVIEWED",
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

    if not temporary.is_file() or temporary.stat().st_size < 50_000:
        raise RuntimeError("PDF build did not produce a plausible output file")
    with temporary.open("rb") as handle:
        if handle.read(5) != b"%PDF-":
            raise RuntimeError("generated file does not have a PDF signature")
    temporary.replace(output)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()
    return int(document.page), output.stat().st_size, digest


def main() -> int:
    args = parse_args()
    output = args.output.resolve()
    if args.clean:
        return clean_output(output)
    source = args.source.resolve()
    if not source.is_file():
        print(f"source manuscript is absent: {source}", file=sys.stderr)
        return 1
    pages, size, digest = build_pdf(source, output)
    print(f"wrote {output}")
    print(f"pages: {pages}")
    print(f"bytes: {size}")
    print(f"sha256: {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
