"""Physical units, saved source/state scope, and creation-only authenticated diagnostics."""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import json
from pathlib import Path

import numpy as np
import pytest

import assess_nsc_discovery_eft_domain as assess


def test_bertotti_robinson_has_zero_scalar_and_weyl_but_nonzero_physical_curvature():
    factor = abs(assess.coupling.locked_coefficients()["C_W"])/assess.coupling.locked_coefficients()["A"]
    values = []
    for radius in (1., 2.):
        Q = np.array([.5])
        package = assess.tidal.curvature_from_local_jets(
            Q, np.array([radius]), 0.*Q, 0.*Q, -Q**3, 0.*Q, 0.*Q,
            0.*Q, 0.*Q, 0.*Q, 0.*Q, 0.*Q)
        result = assess.physical_curvature(package, factor, 0)
        eps = result["cw_over_a_times_physical_scale"]
        assert eps["abs_R4"] == pytest.approx(0., abs=1e-14)
        assert eps["sqrt_Weyl2"] == pytest.approx(0., abs=1e-14)
        assert eps["Riemann_components"] == pytest.approx(factor/radius**2)
        assert eps["sqrt_abs_K"] == pytest.approx(factor*np.sqrt(8)/radius**2)
        assert result["physical_values"]["Ricci2"]["max_abs"] == pytest.approx(4/radius**4)
        values.append(eps["Riemann_components"])
    assert values[0] == pytest.approx(4*values[1])
    assert assess.units_and_scope()["units"]["Q_L_chi_R_h"] == "dimensionless"
    assert assess.units_and_scope()["raw_R_h_used_as_physical_curvature"] is False


@pytest.fixture(scope="module")
def report():
    if not (assess.DEFAULT_DYNAMIC/"manifest.json").exists() or not (assess.DEFAULT_BALANCED/"manifest.json").exists():
        pytest.skip("local saved production datasets are unavailable")
    return assess.build_report()


def test_actual_stations_preserve_chart_stop_and_do_not_relabel_occupied_modes(report):
    assert len(report["rows"]) == 24
    assert report["abs_CW_over_A"] == pytest.approx(.01875)
    assert report["extra_flat_pole_mass"] == pytest.approx(np.sqrt(80/3))
    uniform = next(row for row in report["rows"] if row["preparation_case"] == "uniform" and row["actual_time"] == 0.)
    Q = uniform["Q"]["min"]
    expected = np.sqrt(np.repeat(np.array([np.pi/8, 3*np.pi/8, 5*np.pi/8]), 2)**2+Q**2)
    assert np.max(abs(np.asarray(uniform["frequencies"]["coordinate_positive_levels"])-expected)) < 1e-11
    assert uniform["frequencies"]["source_temporal_RMS_over_pole_range"][1] == pytest.approx(.17967048665, abs=1e-9)
    lower = next(row for row in report["rows"] if row["preparation_case"] == "chi_lower" and row["requested_station"] == 3.)
    assert lower["status"] == "chart_exit" and lower["actual_time"] < 1.5
    assert lower["requested_station_reached"] is False
    for row in report["rows"]:
        assert row["state_unchanged"] and row["chart_positive"]
        assert row["frequencies"]["evolved_source_occupies_six_instantaneous_modes_claimed"] is False
        assert row["curvature"]["chi_substituted"] is False
    final = next(row for row in report["rows"] if row["preparation_case"] == "uniform" and row["actual_time"] == 3.)
    assert final["curvature"]["cw_over_a_times_physical_scale"]["normal_tidal_components"] > 4e6
    assert final["curvature"]["cw_over_a_times_physical_scale"]["sqrt_Weyl2"] < 1
    assert report["input_and_source_hashes_unchanged"] is True


def test_postprocessor_has_no_preparation_or_evolution_call(monkeypatch):
    directory = assess.DEFAULT_DYNAMIC
    if not (directory/"manifest.json").exists():
        pytest.skip("saved stations unavailable")
    def forbidden(*args, **kwargs):
        raise AssertionError("postprocessing must not prepare or advance a state")
    monkeypatch.setattr(assess.nested, "initial_state", forbidden)
    monkeypatch.setattr(assess.nested, "rk4_step", forbidden)
    monkeypatch.setattr(assess.episode, "advance_case", forbidden)
    row = assess.measure_station(directory, "nf128_pattern_dt0.0005", 2, .01875, np.sqrt(80/3))
    assert row["actual_time"] == 1. and row["state_unchanged"]


def test_default_is_readonly_and_write_requires_commit(tmp_path, report, monkeypatch, capsys):
    monkeypatch.setattr(assess, "build_report", lambda *args: report)
    assert assess.main([]) == 0
    assert json.loads(capsys.readouterr().out)["record_written"] is False
    assert not list(tmp_path.iterdir())
    with pytest.raises(SystemExit):
        assess.main(["--write", str(tmp_path/"not-created.json")])
    with pytest.raises(SystemExit):
        assess.main(["--producer-commit", "0"*40])
    assert not list(tmp_path.iterdir())


def test_exact_output_integrity_closure_and_numeric_replay(tmp_path, report, monkeypatch):
    # This unit fixture checks dispatch/integrity. It is not a published Git provenance proof.
    calls = []
    def authenticate(commit, hashes):
        calls.append((commit, dict(hashes)))
        return {"authenticated": True, "test_fixture": True}
    monkeypatch.setattr(assess, "authenticate_sources", authenticate)
    output = tmp_path/"test-only-record.json"
    stored = assess.write_report(report, output, "0"*40)
    assert calls[0] == ("0"*40, report["measurement_source_hashes"])
    assert len(calls) == 1+len(report["saved_producer_bindings"])
    assert output.stat().st_mode & 0o222 == 0
    assert assess.check_report(output)["rows_checked"] == 24
    with pytest.raises(FileExistsError, match="exact"):
        assess.write_report(report, output, "0"*40)
    output.chmod(0o644)
    tampered = dict(stored)
    tampered["rows"] = []
    output.write_text(json.dumps(tampered)); output.chmod(0o444)
    with pytest.raises(ValueError, match="digest"):
        assess.check_report(output)


def test_provenance_dispatch_and_authentication_failure_creates_nothing(tmp_path, report, monkeypatch):
    calls = []
    def blob(root, path, expected, *, commit):
        calls.append((root, path, expected, commit))
        return b"test-dispatch-only"
    monkeypatch.setattr(assess.provenance, "resolve_pinned_source_bytes", blob)
    result = assess.authenticate_sources("1"*40, {"lab/scripts/fixture.py": "2"*64})
    assert result["authenticated"] and calls == [(assess.REPO, "lab/scripts/fixture.py", "2"*64, "1"*40)]
    def failure(*args, **kwargs): raise RuntimeError("frozen source unavailable")
    monkeypatch.setattr(assess, "authenticate_sources", failure)
    output = tmp_path/"uncreated.json"
    with pytest.raises(RuntimeError, match="unavailable"):
        assess.write_report(report, output, "0"*40)
    assert not output.exists()
