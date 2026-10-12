"""Bounded source-free control, with original discrete rates and no campaign input."""
import json

import numpy as np
import pytest

import derive_nsc_discovery_vacuum_control as cli
from recursive_horizons import nsc_discovery_vacuum_control as control


@pytest.fixture(scope="module")
def measured():
    return control.run_control()


def test_locked_action_selects_radius_and_canonical_momenta():
    grid = control.control_grid(32)
    values = control.exact_values(grid, 3.0)
    assert values["radius_square"] == pytest.approx(1.0)
    assert values["Q0"] == pytest.approx(0.2513249333)
    assert grid.gauge == "conformal"
    assert np.array_equal(grid.fine.occupations, np.zeros(6))
    coefficients = control.coupling.locked_coefficients()
    assert coefficients == grid.fine.coefficients
    D = values["D"]
    p = control.action.momenta(
        values["Q"], values["Q"], values["Q_t"], 0.0, values["radius"],
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, grid.fine.A, grid.fine.C_W,
    )
    assert p == pytest.approx((0.0, values["p_r"], values["p_chi"]))
    assert values["p_chi"] == pytest.approx(4 * values["alpha"] * D)
    assert values["Q_t"] < 0
    assert D**2 + values["Q"]**2 == pytest.approx(values["Q0"]**2)
    assert control.action.feedback_V(values["radius"], 0.0, grid.fine.A,
                                    grid.fine.C_W, grid.fine.C_F, grid.fine.flux) == pytest.approx(0.0)


def test_exact_solution_satisfies_original_rates_constraints_and_tides(measured):
    report, _arrays = measured
    assert report["cpu_seconds"] < 30.0
    assert report["source_hashes_unchanged"] and report["input_hashes_unchanged"]
    for rows in report["exact"].values():
        assert [row["coordinate_time"] for row in rows] == list(control.EXACT_TIMES)
        for row in rows:
            assert max(row["rate_residuals"].values()) < 1e-7
            assert row["constraints"]["full_hamilton_max"] < 1e-7
            assert row["constraints"]["full_momentum_max"] < 1e-7
            assert abs(row["energy"]) < 1e-7
            assert row["Pi_max"] < 1e-10
            assert max(v["max_error_from_exact"] for v in row["metric_invariants"].values()) < 1e-6
            assert row["metric_invariants"]["R_0101"]["mean"] == pytest.approx(-1.0, abs=1e-7)
            assert row["metric_invariants"]["K"]["mean"] == pytest.approx(8.0, abs=1e-6)


def test_zero_covariance_is_car_admissible_without_packet_gram_identity(measured):
    report, arrays = measured
    assert report["source"]["occupations"] == [0.0] * 6
    assert report["source"]["occupied_rank6_Gram_identity_required"] is False
    assert report["source"]["magnetic_and_induced_action_retained"] is True
    for name, values in arrays.items():
        if name.endswith(("_phi0", "_phi1")):
            assert np.array_equal(values, np.zeros_like(values))
    for rows in report["exact"].values():
        for row in rows:
            assert row["covariance_max"] == 0.0
            assert all(value == 0 for value in row["covariance_eigenvalues"])
            assert row["source_force_max"] == row["fieldwork_power"] == 0.0


def test_short_independent_rk4_matches_formula_and_records_nodal_pi(measured):
    report, arrays = measured
    for nf, case in report["numerical"].items():
        assert case["steps"] == 300
        assert [row["coordinate_time"] for row in case["stations"]] == list(control.RUN_TIMES)
        late = case["stations"][-1]
        assert max(late["state_errors"].values()) < 1e-6
        assert abs(late["energy"]) < 1e-7
        assert late["constraints"]["full_hamilton_max"] < 1e-6
        assert max(v["max_error_from_exact"] for v in late["metric_invariants"].values()) < 1e-5
        for name in ("Q", "r", "chi"):
            nodal = arrays[f"nf{nf}_station2_p_{name}"]
            canonical = arrays[f"nf{nf}_station2_canonical_pi_{name}"]
            assert np.array_equal(canonical, arrays[f"nf{nf}_dx_g"] * nodal)
        early = case["stations"][0]
        assert late["proper_lengths"]["period"] < early["proper_lengths"]["period"]
        assert late["r"] == pytest.approx(early["r"], abs=1e-6)
    assert report["scope"]["same_initial_geometry_as_discovery_family"] is False
    assert report["scope"]["family_instability_decided"] is False


def test_preview_writes_nothing_and_checker_does_not_evolve(measured, tmp_path, monkeypatch):
    report, arrays = measured
    monkeypatch.setattr(control, "run_control", lambda: (report, arrays))
    assert cli.main([]) == 0
    assert not list(tmp_path.iterdir())
    destination = tmp_path / "control-v1"
    monkeypatch.setattr(control, "OUTPUT", destination)
    assert cli.main(["--write", str(destination)]) == 0
    with pytest.raises(FileExistsError):
        cli.main(["--write", str(destination)])
    monkeypatch.setattr(control, "run_control", lambda: pytest.fail("checker evolved"))
    assert cli.main(["--check", str(destination)]) == 0
    saved = json.loads((destination / "record.json").read_text())
    assert saved["payload_bytes"] < control.MAX_PAYLOAD_BYTES
    with (destination / "payload.npz").open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="payload hash changed"):
        cli.main(["--check", str(destination)])


def test_output_scope_is_bounded():
    with pytest.raises(ValueError, match="output must"):
        control.output_directory("/tmp/not-an-owned-control-output")
    with pytest.raises(ValueError, match="nf=32 or nf=64"):
        control.control_grid(128)
