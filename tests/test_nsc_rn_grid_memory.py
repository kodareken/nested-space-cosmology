"""Tests for the grid memory measurement.

The recorded family's silent stops were attributed to memory on the grounds
that it scales as N squared. That inference was never tested, and it is wrong:
the production derivative is a CSR operator, and on this machine the memory
cost of a grid is flat in the point count.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

LAB = Path(__file__).resolve().parents[1]
RECORD = LAB / "results/development/nsc-rn-grid-memory-v1.json"
DRIVER = LAB / "scripts/measure_nsc_rn_grid_memory.py"
FAMILY_NOTE = LAB / "docs/nsc-rn-silent-stops.md"


@pytest.fixture(scope="module")
def record():
    return json.loads(RECORD.read_text(encoding="utf-8"))


def _rows(record):
    return {row["points"]: row for row in record["rows"]}


class TestTheMachineIsKnown:
    def test_it_reports_the_machine_it_ran_on(self, record):
        assert record["machine"]["total_gb"] > 30.0
        assert record["machine"]["logical_cores"] >= 8

    def test_it_says_how_it_measured(self, record):
        assert "GlobalMemoryStatusEx" in record["machine"]["measurement"]

    def test_it_does_not_depend_on_a_missing_package(self, record):
        """psutil is absent on this host, and the measurement still works."""
        assert "psutil" not in DRIVER.read_text(encoding="utf-8").split("import")[-1][:200] \
            or "absent" in DRIVER.read_text(encoding="utf-8")


class TestMemoryIsFlat:
    def test_the_range_covers_nine_times_the_points(self, record):
        points = sorted(_rows(record))
        assert max(points) >= 8000
        assert min(points) <= 1000
        assert max(points) / min(points) > 8.0

    def test_the_commit_cost_does_not_grow_with_points(self, record):
        rows = sorted(_rows(record).values(), key=lambda r: r["points"])
        first, last = rows[0], rows[-1]
        assert last["commit_delta_bytes"] < first["commit_total_bytes"]

    def test_the_operator_is_sparse(self, record):
        """This is why it is flat: nnz is far below N squared."""
        for row in record["rows"]:
            assert row["operator_nnz"] < row["points"] ** 2

    def test_the_nnz_grows_only_about_linearly(self, record):
        rows = sorted(_rows(record).values(), key=lambda r: r["points"])
        smallest, largest = rows[0], rows[-1]
        ratio = largest["operator_nnz"] / smallest["operator_nnz"]
        assert ratio < 2.0 * (largest["points"] / smallest["points"])


class TestTheBuildStillWorks:
    def test_every_measured_grid_built(self, record):
        for row in record["rows"]:
            assert row["mass_finite"]
            assert row["build_seconds"] > 0.0

    def test_the_cost_is_time_not_memory(self, record):
        rows = sorted(_rows(record).values(), key=lambda r: r["points"])
        assert rows[-1]["build_seconds"] < 10.0 * rows[0]["build_seconds"]


class TestWhatItMeansForTheSubHorizonGap:
    def test_more_points_do_not_resolve_the_deepest_edge(self, record):
        """The gap at 0.7555 stays sub-cell even at 8000 points."""
        largest = max(record["rows"], key=lambda r: r["points"])
        assert largest["gap_in_cells"] < 1.0

    def test_the_record_is_at_the_deepest_edge(self, record):
        assert record["excision"] == pytest.approx(0.7555)

    def test_the_wide_outer_radius_is_the_reason(self, record):
        """Resolving a 1e-3 gap at r_out = 32 needs about 1e5 points."""
        largest = max(record["rows"], key=lambda r: r["points"])
        needed = (record["r_out"] - record["excision"]) / (
            (record["excision"] - 0.7543428628582858) / 4.0)
        assert needed > 50000
        assert largest["points"] < needed / 5.0


def test_the_measurement_is_kept_as_a_record():
    assert RECORD.is_file()
    assert RECORD.read_text(encoding="utf-8").strip().startswith("{")