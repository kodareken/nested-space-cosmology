"""Census and established-control locks. This file does not solve a row."""
from pathlib import Path
import resource
import runpy
import pytest
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_subgap_source_window import (
    EXPECTED_POSITIVE_ROWS, PANEL_ROWS, _cpu_limit, settings_for, subgap_catalogue)
ROOT = Path(__file__).resolve().parents[1]


def test_established_controls_cover_the_54_and_do_not_stick_the_cpu_limit():
    assert sum(len(rows) for _panel, rows in PANEL_ROWS) == EXPECTED_POSITIVE_ROWS == 54
    assert [panel for panel, _rows in PANEL_ROWS] == [
        'group14/low16_1', 'group14/low32_1', 'subgap8/14_1']
    assert PANEL_ROWS[0][1] == tuple(range(1, 15))
    assert PANEL_ROWS[1][1] == tuple(range(0, 32))
    assert PANEL_ROWS[2][1] == tuple(range(0, 8))
    standard = settings_for('group14/low16_1', 1)
    wide = settings_for('subgap8/14_1', 7)
    assert standard == {
        'transport_profile': 'row15-mixed-transport', 'tube': '0.00001', 'defect_subdivisions': 4}
    assert wide == {
        'transport_profile': 'near-threshold-wider-tube', 'tube': '0.001', 'defect_subdivisions': 16}
    before = resource.getrlimit(resource.RLIMIT_CPU)
    with _cpu_limit(30):
        assert resource.getrlimit(resource.RLIMIT_CPU)[1] == before[1]
    assert resource.getrlimit(resource.RLIMIT_CPU) == before
    source = (ROOT / 'scripts/derive_nsc_subgap_source_window.py').read_text()
    assert '--check' in source and 'solve_jost' not in source


def test_catalogue_is_the_uncovered_subgap_rows():
    archive = RetainedUpstreamArchive(ROOT)
    rows = subgap_catalogue(archive)
    assert [(row.panel, row.row) for row in rows] == [
        (panel, index) for panel, indexes in PANEL_ROWS for index in indexes]
    assert rows[0].energy.hex() == '0x1.c60a99e906500p-8'
    assert rows[-1].energy.hex() == '0x1.8f38f9b035a88p+0'
    assert all(row.negative_energy == -row.energy for row in rows)
    assert all(row.negative_angular == -row.angular for row in rows)
    assert all(row.positive_source_digest != row.negative_source_digest for row in rows)
    driver = runpy.run_path(str(ROOT / 'scripts/derive_nsc_subgap_source_window.py'))
    with pytest.raises(SystemExit):
        driver['main'](['--check', '--row-budget', '1'])
