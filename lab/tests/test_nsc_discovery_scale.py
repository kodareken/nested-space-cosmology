"""Initial source-selected scale: owned constraint, not a dynamical solve."""
import json
import subprocess
import sys

import numpy as np
import pytest

from recursive_horizons import nsc_discovery_scale as scale
from recursive_horizons import nsc_spherical_coupling as coupling


@pytest.fixture(scope="module")
def report():
    return scale.derivation_report()


def test_uniform_source_analytic_control(report):
    analytic = report["analytic"]
    assert analytic["uniform_residual_max"] < 1e-9
    assert analytic["uniform_formula_gap"] < 1e-12
    assert analytic["uniform_clock_rate"] == pytest.approx(
        analytic["uniform_radius"] * scale.owned_couplings()["Q"])
    assert analytic["uniform_proper_ratio"] == pytest.approx(analytic["uniform_radius"])
    comparison = analytic["variable_owner_comparison"]
    assert comparison["identity_gap"] < 1e-9
    assert comparison["owner_residual_max"] == pytest.approx(comparison["contrast_max"])
    assert comparison["owner_residual_max"] > 1e-4
    assert analytic["one_percent_pass_declared"] is False


def test_integrated_identity_and_max_principle(report):
    analytic = report["analytic"]
    assert analytic["integrated_identity_gap"] < 1e-9
    assert analytic["manufactured_G_min"] > 0.0
    assert analytic["max_principle_holds"] is True
    assert analytic["max_principle"]["lower_applicable"] is True
    shape = analytic["shape"]
    assert shape["same_flat_radius"] is True
    assert shape["contrast_increased"] is True
    assert shape["cases"][1]["conditional_absolute_gap"] > shape["cases"][0]["conditional_absolute_gap"]


def test_owned_projected_radius_matches_derived_weak_form(report):
    projection = report["projection"]
    assert projection["algebraic_owner_gap"] < 1e-8
    assert projection["projected_pull_gap"] < 1e-8
    assert projection["defect_explains_y_gap"] < 1e-8
    assert projection["owned_coefficient"] == 4.0
    assert projection["reviewer_unit_coefficient_gap"] > 1e-4
    assert projection["reviewer_unit_coefficient_gap"] > 100.0 * projection["algebraic_owner_gap"]


def test_coordinate_scaling_leaves_areal_radius_weight_zero(report):
    coordinate = report["coordinate"]
    assert coordinate["areal_radius_weight"] == pytest.approx(0.0, abs=1e-12)
    assert coordinate["pulled_r_flat2"] == pytest.approx(coordinate["base_r_flat2"])
    assert coordinate["period_only_changes_radius"] is True
    pullback = coordinate["pullback"]
    assert pullback["areal_radius_weight"] == pytest.approx(0.0, abs=1e-12)
    assert pullback["areal_radius_unchanged"] is True
    assert pullback["Q_weight"] == pytest.approx(1.0)
    assert pullback["proper_ratio_weight"] == pytest.approx(0.0)
    assert pullback["clock_rate_weight"] == pytest.approx(1.0)


def test_heldout_phase_frequency_gram_current_and_car(report):
    held = report["heldout"]
    assert held["saved_layouts_used_as_held_out"] is False
    assert held["initial_solver_called"] is False
    assert held["production_record_written"] is False
    members = held["members"]
    assert [member["phase_sign"] for member in members] == [1.0, -1.0]
    integrals = []
    for member in members:
        assert member["held_out"] is True
        assert member["role"] == "held_out_initial_preparation"
        assert member["carrier_frequency"] == scale.HELD_OUT_CARRIER_NU
        assert member["carrier_frequency"] != coupling.CARRIER_K
        assert member["gram_max"] < 1e-9
        assert member["column_gram_max"] < 1e-9
        assert member["child_frame_repaired"] is False
        assert member["child_source_equals_reference"] is False
        assert member["source_equals_reference"] is False
        assert member["mean_current_deleted"] is False
        assert member["shift_residual_mean_subtracted"] is False
        assert member["occupations"] == coupling.OCCUPATIONS.tolist()
        assert member["rho_independent_of_r_max"] == 0.0
        assert member["source_positivity"] is True
        assert member["S1_real_max"] < 1e-7
        assert member["prediction"]["solver_called"] is False
        assert member["prediction"]["one_percent_pass_declared"] is False
        assert member["prediction"]["hand_set_hubble"] is False
        assert member["prediction"]["hand_set_magnetic_field"] is False
        assert member["prediction"]["probe_radius_used_as_prediction"] is False
        comparison = member["owner_comparison"]
        assert comparison["relative_identity_gap"] < 1e-9
        assert comparison["owner_residual_max"] > 1.0
        assert member["prediction"]["clock_rate"] == pytest.approx(
            member["prediction"]["radius"] * scale.owned_couplings()["Q"])
        assert member["prediction"]["proper_ratio"] == pytest.approx(member["prediction"]["radius"])
        for column in member["column_signs"]:
            assert column["K_sum"] > 0.0
            assert column["S1_abs_max"] < 1e-7
            if column["chirality"] == "plus":
                assert column["Pmom_sum"] > 0.0
            else:
                assert column["Pmom_sum"] < 0.0
        integrals.append(member["rho_integral"])
    assert integrals[0] == pytest.approx(integrals[1], abs=1e-9)
    assert held["phase_rho_integral_gap"] < 1e-9
    control = held["frequency_control"]
    assert control["held_out"] is False
    assert control["role"] == "override_not_held_out"
    assert control["carrier_frequency"] == coupling.CARRIER_K
    assert held["frequency_changes_kinetic_invariant"] is True
    assert abs(held["kinetic_sum_gap"]) / abs(members[0]["K_sum"]) > 1e-3


def test_saved_layouts_are_not_held_out():
    for layout in ("original", "separated"):
        assert scale.preparation_role(layout) == "previously_saved_layout"
        assert scale.preparation_role(
            layout, carrier_frequency=scale.HELD_OUT_CARRIER_NU, chiral_family=True,
        ) == "previously_saved_layout"
    assert scale.preparation_role(
        "override", carrier_frequency=scale.HELD_OUT_CARRIER_NU, chiral_family=True,
    ) == "held_out_initial_preparation"
    assert scale.preparation_role(
        "override", carrier_frequency=coupling.CARRIER_K, chiral_family=True,
    ) == "override_not_held_out"


def test_prediction_uses_measured_invariants_and_states_applicability(report):
    member = report["heldout"]["members"][0]
    measured = member["measured"]
    constants = scale.owned_couplings()
    expected = scale.uniform_radius_square(
        measured["rmag2"], measured["rho_integral"], measured["A"], measured["Q"], measured["length"],
    )
    assert member["prediction"]["r_flat2"] == pytest.approx(expected)
    assert expected == pytest.approx(
        constants["rmag2"] + measured["rho_integral"] / (
            8.0 * np.pi * constants["A"] * constants["Q"] * measured["length"]
        )
    )
    assert member["prediction"]["conditional_bound"]
    assert member["prediction"]["exact_for_uniform_positive_G"] is False
    assert member["G_contrast_max"] > 0.0
    assert member["prediction"]["conditional_relative_gap"] == pytest.approx(
        member["G_contrast_max"] / member["G_mean"])
    assert "one_percent" not in json.dumps(member["prediction"]["conditional_bound"])
    assert member["owner_comparison"]["relative_identity_gap"] < 1e-9
    assert member["owner_comparison"]["owner_residual_max"] == pytest.approx(
        member["owner_comparison"]["contrast_max"], rel=1e-9)


def test_multiplicity_and_source_signs_match_owned_forces(report):
    measured = report["heldout"]["members"][0]["measured"]
    assert measured["multiplicity"] == 4 * measured["kappa"]
    assert measured["multiplicity_applied_once"] is True
    assert measured["ambient_shift_changes_force_L"] is False
    current = measured["current"]
    assert measured["current_mean"] == pytest.approx(float(np.mean(current)))
    assert measured["shift_residual_mean_subtracted"] is False
    # The retained mean is the occupation-weighted chiral sum, not a forced zero.
    assert measured["current_mean"] == pytest.approx(float(np.mean(current)), abs=0.0)


def test_script_prints_summary_and_refuses_a_record():
    script = scale.__file__.replace(
        "src/recursive_horizons/nsc_discovery_scale.py",
        "scripts/derive_nsc_discovery_scale.py",
    )
    completed = subprocess.run(
        [sys.executable, script],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    payload = json.loads(completed.stdout)
    assert payload["initial_solver_called"] is False
    assert payload["production_record_written"] is False
    assert payload["uniform_residual_max"] < 1e-9
    assert payload["algebraic_owner_gap"] < 1e-8
    assert payload["areal_radius_weight"] == pytest.approx(0.0, abs=1e-12)
    assert payload["held_out_source_positivity"] is True
    refused = subprocess.run(
        [sys.executable, script, "--solve", "--output", "unused.json"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert refused.returncode != 0
