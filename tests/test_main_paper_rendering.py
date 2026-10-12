"""The living paper renders and authenticates local scientific illustrations."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("main_paper_renderer", ROOT / "scripts/build_paper.py")
paper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(paper)


@pytest.fixture
def illustrated_source(tmp_path, monkeypatch):
    from PIL import Image

    monkeypatch.setattr(paper, "REPOSITORY", tmp_path)
    (tmp_path / "paper/figures").mkdir(parents=True)
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts/build_paper.py").write_bytes(Path(paper.__file__).read_bytes())
    Image.new("RGB", (400, 200), "navy").save(tmp_path / "paper/figures/a map.png")
    (tmp_path / "results").mkdir()
    (tmp_path / "results/response.json").write_text('{"domain": "finite"}')
    source = tmp_path / "paper/main.md"
    source.write_text(
        "# Nested-Space Cosmology\n\n## Abstract\n\nA finite observer response.\n\n"
        "![Observer map](<figures/a map.png>)\n\n*Figure 1. A local chart.*\n\n"
        "$$\n\\phi_C(T)=I/T^{\\prime}\n$$\n\n"
        "| Observer | Access |\n|---|---|\n| Centre | Conditional early pulse |\n\n"
        "[Evidence](../results/response.json)\n", encoding="utf-8")
    return source


def test_image_and_evidence_bytes_are_bound(illustrated_source):
    inputs = {row["path"]: row["sha256"] for row in paper.scientific_build_inputs(illustrated_source)}
    for path in ("paper/figures/a map.png", "results/response.json", "scripts/build_paper.py"):
        assert inputs[path] == hashlib.sha256((paper.REPOSITORY / path).read_bytes()).hexdigest()
    original = inputs["paper/figures/a map.png"]
    (paper.REPOSITORY / "paper/figures/a map.png").write_bytes(b"replaced")
    changed = {row["path"]: row["sha256"] for row in paper.scientific_build_inputs(illustrated_source)}
    assert changed["paper/figures/a map.png"] != original


@pytest.mark.parametrize("target", ["https://example.org/map.png", "../../escape.png", "figures/missing.png"])
def test_missing_remote_or_escaping_image_is_rejected(illustrated_source, target):
    with pytest.raises(ValueError):
        paper.local_image_path(illustrated_source, target)


def test_illustrated_math_and_table_render_identically_twice(illustrated_source, tmp_path):
    metadata = json.loads((ROOT / "paper/metadata.json").read_text(encoding="utf-8"))
    first, second = tmp_path / "first.pdf", tmp_path / "second.pdf"
    result_a = paper.build_pdf(illustrated_source, first, metadata)
    result_b = paper.build_pdf(illustrated_source, second, metadata)
    assert result_a == result_b
    assert first.read_bytes() == second.read_bytes()
    # ReportLab emits a PDF image XObject for the map and for display math.
    assert first.read_bytes().count(b"/Subtype /Image") >= 2

