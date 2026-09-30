"""Direct-order conformal ADM source, against the unchanged historical owner.

The historical finite-grid product rule is a different ordering. These
checks do not mark its recorded results false. Naive C=I/2 is the
one-block spin identity, not the recorded 6-mode state.
"""
from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_adm_source import adm_hamiltonian
from recursive_horizons.nsc_conformal_adm_source import (
    FRAME_MAPPING,
    PRODUCT_RULE_NAME,
    SPIN_IDENTITY_NOTE,
    VACUUM_BRANCH_DEFINITION,
    VACUUM_BRANCH_OWNER,
    block_source,
    conformal_source,
    direct_hamiltonian,
    fixed_gaussian_covariance,
    ordering_comparison,
    physical_forces_from_conformal,
    prescribed_work_ledger,
    sector_multiplicity,
    spin_identity_covariance,
)
from recursive_horizons.nsc_covariant_operator import cylinder_metric, smooth_metric
from recursive_horizons.nsc_influence import ground_covariance


LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
HISTORICAL_ADM = (
    "0dba93a1090728c0c26fca2f04ac13ae74201fbb66b1aee09155a7146e5d2d09"
)
HISTORICAL_PATHS = (
    LAB / "src/recursive_horizons/nsc_adm_source.py",
    ROOT / "src/recursive_horizons/nsc_adm_source.py",
)


def _metric():
    return smooth_metric(12, general=True)


def _shift(metric):
    angle = 2 * np.pi * metric.x / metric.length
    return 0.08 * np.cos(angle)


def _state(metric):
    return fixed_gaussian_covariance(2 * metric.points)


def _energy(metric, shift, state, kappa):
    return float(np.trace(state @ direct_hamiltonian(metric, shift, kappa)).real)


def _central(metric, shift, state, kappa, field, node, step):
    samples = []
    for sign in (-1.0, 1.0):
        if field == "L":
            radius = metric.sphere_radius
            length = metric.lapse / radius
            radial = metric.radial_scale / radius
            length = length.copy()
            length[node] += sign * step
            changed = replace(
                metric,
                lapse=length * radius,
                radial_scale=radial * radius,
            )
            beta = shift
        elif field == "Q":
            radius = metric.sphere_radius
            length = metric.lapse / radius
            radial = metric.radial_scale / radius
            radial = radial.copy()
            radial[node] += sign * step
            changed = replace(
                metric,
                lapse=length * radius,
                radial_scale=radial * radius,
            )
            beta = shift
        elif field == "beta":
            changed = metric
            beta = np.array(shift, dtype=float, copy=True)
            beta[node] += sign * step
        elif field == "N":
            lapse = metric.lapse.copy()
            lapse[node] += sign * step
            changed = replace(metric, lapse=lapse)
            beta = shift
        elif field == "q":
            radial = metric.radial_scale.copy()
            radial[node] += sign * step
            changed = replace(metric, radial_scale=radial)
            beta = shift
        elif field == "r":
            radius = metric.sphere_radius.copy()
            radius[node] += sign * step
            changed = replace(metric, sphere_radius=radius)
            beta = shift
        else:
            raise AssertionError(field)
        samples.append(_energy(changed, beta, state, kappa))
    return (samples[1] - samples[0]) / (2 * step)


def test_historical_adm_source_bytes_are_unchanged():
    for path in HISTORICAL_PATHS:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == HISTORICAL_ADM


def test_direct_ordering_is_hermitian_and_matches_constants():
    metric = _metric()
    shift = _shift(metric)
    hamiltonian = direct_hamiltonian(metric, shift, 1)
    assert hamiltonian.shape == (2 * metric.points, 2 * metric.points)
    assert np.max(np.abs(hamiltonian - hamiltonian.conj().T)) < 1e-12
    flat = cylinder_metric(points=12, length=4.0, radius=1.3)
    zero = np.zeros(flat.points)
    assert np.max(np.abs(direct_hamiltonian(flat, zero, 2) - adm_hamiltonian(flat, zero, 2))) < 1e-12


def test_multiplicity_is_four_kappa_once_and_not_a_negative_copy():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    source = conformal_source(metric, shift, state, 3)
    assert sector_multiplicity(3) == 12
    assert source["multiplicity"] == 12
    assert source["energy"] == pytest.approx(12 * source["sector_energy"])
    positive = _energy(metric, shift, state, 3)
    negative = _energy(metric, shift, state, -3)
    assert positive != pytest.approx(negative, rel=0, abs=1e-8)
    doubled = 12 * (positive + negative)
    assert source["energy"] != pytest.approx(doubled, rel=0, abs=1e-6)
    with pytest.raises(ValueError):
        conformal_source(metric, shift, state, -3)
    with pytest.raises(ValueError):
        conformal_source(metric, shift, state, 0)
    with pytest.raises(ValueError):
        sector_multiplicity(True)


def test_finite_difference_gradients_are_nontrivial_and_match():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    kappa = 1
    block = block_source(metric, shift, state, kappa)
    coupled = conformal_source(metric, shift, state, kappa)
    assert coupled["energy"] == pytest.approx(4 * block["energy"])
    assert coupled["hamiltonian_is_one_block"] is True
    assert np.max(np.abs(coupled["hamiltonian"] - direct_hamiltonian(metric, shift, kappa))) < 1e-15
    assert coupled["energy"] == pytest.approx(
        4 * float(np.trace(state @ coupled["hamiltonian"]).real)
    )
    for name in ("L", "Q", "beta", "N", "q", "r"):
        assert np.linalg.norm(block["nodal"][name]) > 1e-3
        for node in (0, 5, 11):
            step = 1e-6
            difference = _central(metric, shift, state, kappa, name, node, step)
            assert block["nodal"][name][node] == pytest.approx(difference, rel=1e-6, abs=1e-7)
            assert coupled["nodal"][name][node] == pytest.approx(4 * difference, rel=1e-6, abs=4e-7)
    chain = physical_forces_from_conformal(
        block["nodal"]["L"], block["nodal"]["Q"], block["nodal"]["beta"],
        metric.lapse / metric.sphere_radius,
        metric.radial_scale / metric.sphere_radius,
        metric.sphere_radius,
    )
    for name in ("N", "q", "r", "beta"):
        assert np.max(np.abs(chain[name] - block["nodal"][name])) < 1e-12
    assert np.max(np.abs(block["mass_radius_force_disagreement"])) < 1e-12


def test_fixed_lq_radius_invariance_and_massless_ward():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    comparison = ordering_comparison(metric, shift, state, 1)
    assert comparison["ordering_name"] == PRODUCT_RULE_NAME
    assert comparison["historical_results_marked_false"] is False
    assert comparison["frame_mapping"] == FRAME_MAPPING
    assert comparison["same_physical_covariance"] is True
    assert comparison["direct_fixed_LQ_hamiltonian_max"] < 1e-12
    assert abs(comparison["direct_fixed_LQ_trace_change"]) < 1e-12
    assert comparison["direct_ward_l2"] < 1e-12
    assert comparison["product_rule_fixed_LQ_hamiltonian_frobenius"] > 1e-3
    assert abs(comparison["product_rule_fixed_LQ_trace_change"]) > 1e-5
    assert comparison["product_rule_ward_l2"] > 1e-4
    assert comparison["physical_radius_force_l2"] > 1e-3
    assert comparison["radius_force_disagreement_l2"] < 1e-12
    assert comparison["shift_force_disagreement_l2"] < 1e-12
    assert comparison["lapse_force_disagreement_l2"] > 1e-4

    angle = 2 * np.pi * metric.x / metric.length
    probe = np.cos(angle)
    base = _energy(metric, shift, state, 1)

    def _weyl(epsilon):
        factor = np.exp(epsilon * probe)
        changed = replace(
            metric,
            lapse=metric.lapse * factor,
            radial_scale=metric.radial_scale * factor,
            sphere_radius=metric.sphere_radius * factor,
        )
        return _energy(changed, shift, state, 1)

    direct_slope = (_weyl(1e-6) - _weyl(-1e-6)) / 2e-6
    assert direct_slope == pytest.approx(0, abs=1e-8)

    def _old_weyl(epsilon):
        factor = np.exp(epsilon * probe)
        changed = replace(
            metric,
            lapse=metric.lapse * factor,
            radial_scale=metric.radial_scale * factor,
            sphere_radius=metric.sphere_radius * factor,
        )
        return float(np.trace(state @ adm_hamiltonian(changed, shift, 1)).real)

    old_slope = (_old_weyl(1e-6) - _old_weyl(-1e-6)) / 2e-6
    assert abs(old_slope) > 1e-5

    massless = block_source(metric, shift, state, 0)
    assert np.linalg.norm(massless["nodal"]["r"]) < 1e-12
    assert np.linalg.norm(massless["massless_ward_residual"]) < 1e-12
    scale = np.exp(1e-6 * probe)
    scaled = replace(
        metric,
        lapse=metric.lapse * scale,
        radial_scale=metric.radial_scale * scale,
    )
    massless_scale = _energy(scaled, shift, state, 0) - _energy(metric, shift, state, 0)
    assert abs(massless_scale) < 1e-9
    radius = metric.sphere_radius.copy()
    radius[4] += 1e-5
    moved = replace(metric, sphere_radius=radius)
    assert abs(_energy(moved, shift, state, 0) - _energy(metric, shift, state, 0)) < 1e-12


def test_naive_c0_spin_identity_raw_zero_is_not_renormalized_source():
    metric = _metric()
    shift = _shift(metric)
    identity = spin_identity_covariance(metric.points)
    source = conformal_source(metric, shift, identity, 1)
    assert abs(source["sector_energy"]) < 1e-12
    assert abs(source["energy"]) < 1e-12
    for values in source["nodal"].values():
        assert np.max(np.abs(values)) < 1e-12
    assert np.max(np.abs(source["massless_ward_residual"])) < 1e-12
    assert source["limits"]["naive_spin_identity_means_zero_renormalized_source"] is False
    assert "not" in SPIN_IDENTITY_NOTE or "does not" in SPIN_IDENTITY_NOTE
    assert source["spin_identity_note"] == SPIN_IDENTITY_NOTE
    assert source["limits"]["six_mode_embedding"] is False


def test_relative_difference_is_explicit_and_vacuum_is_not_silent():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    sea = ground_covariance(direct_hamiltonian(metric, shift, 1))[2]
    naked = conformal_source(metric, shift, state, 2)
    sea_source = conformal_source(metric, shift, sea, 2)
    difference = state - sea
    named = conformal_source(
        metric, shift, difference, 2,
        state_interface="relative_difference",
        relative_name="gaussian_minus_canonical_sea",
    )
    assert named["vacuum_subtracted"] is False
    assert named["relative_name"] == "gaussian_minus_canonical_sea"
    assert named["energy"] == pytest.approx(naked["energy"] - sea_source["energy"])
    for name in named["nodal"]:
        assert np.max(np.abs(
            named["nodal"][name] - (naked["nodal"][name] - sea_source["nodal"][name])
        )) < 1e-10
    assert naked["energy"] != pytest.approx(named["energy"], rel=0, abs=1e-6)
    assert naked["limits"]["vacuum_branch_accounted_in_this_module"] is False
    assert naked["limits"]["vacuum_branch_owner"] == VACUUM_BRANCH_OWNER
    assert naked["limits"]["vacuum_branch_definition"] == VACUUM_BRANCH_DEFINITION
    assert naked["limits"]["absolute_renormalized_stress"] is False
    assert naked["limits"]["full_source_closure"] is False
    assert naked["limits"]["frame_mapping"] == "missing"
    with pytest.raises(ValueError):
        conformal_source(
            metric, shift, difference, 1, state_interface="relative_difference",
        )
    with pytest.raises(ValueError):
        conformal_source(
            metric, shift, difference, 1, state_interface="covariance",
        )


def test_prescribed_geometry_work_ledger_is_not_autonomous():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    ledger = prescribed_work_ledger(metric, shift, state, 1)
    assert ledger["autonomous_metric_evolution"] is False
    assert ledger["geometry_role"] == "prescribed_control"
    assert ledger["vacuum_branch_included"] is False
    assert ledger["multiplicity"] == 4
    assert ledger["unitary_residual_max"] < 1e-9
    assert ledger["work_residual_max"] < 1e-9
    assert abs(ledger["ledger_closure"]) < 1e-9
    assert abs(ledger["total_work"]) > 1e-4
    assert ledger["energy_change"] == pytest.approx(ledger["total_work"], abs=1e-9)
    assert abs(ledger["linear_residual"]) < 1e-8
    assert ledger["force_alignment_of_prescribed_step"] > -0.95
    assert ledger["final_covariance_min"] > -1e-9
    assert ledger["final_covariance_max"] < 1 + 1e-9


def test_same_covariance_product_rule_magnitude_is_reported():
    metric = _metric()
    shift = _shift(metric)
    state = _state(metric)
    comparison = ordering_comparison(metric, shift, state, 1, amplitude=0.12)
    assert comparison["source_interfaces"] == ["covariance", "relative_difference"]
    assert comparison["historical_results_marked_false"] is False
    assert comparison["physical_radius_force_l2"] == pytest.approx(
        comparison["historical_radius_force_l2"], abs=1e-12,
    )
    # One-block control on the fixed Gaussian covariance, M=12 smooth cell.
    # The product-rule trace change is not F_r, and neither number retracts
    # the historical factorization.
    assert comparison["direct_ward_l2"] == pytest.approx(0, abs=1e-12)
    assert comparison["product_rule_ward_l2"] == pytest.approx(0.010872496849216692, rel=1e-9, abs=1e-15)
    assert comparison["physical_radius_force_l2"] == pytest.approx(0.03308026372835154, rel=1e-9, abs=1e-15)
    assert comparison["product_rule_fixed_LQ_trace_change"] == pytest.approx(
        0.00014262357291895732, rel=1e-9, abs=1e-15,
    )
    assert comparison["product_rule_fixed_LQ_hamiltonian_frobenius"] == pytest.approx(
        0.0533295845381964, rel=1e-9, abs=1e-15,
    )
    assert comparison["lapse_force_disagreement_l2"] == pytest.approx(
        0.021077411660050126, rel=1e-9, abs=1e-15,
    )


def test_sea_fixed_lq_defect_is_not_the_radius_force():
    """Canonical sea of this block, not a recorded 6-mode state.

    On the M=12 smooth cell the product-rule fixed-(L, Q) trace change is
    about -0.0461 while ||F_r||_2 is about 0.650. That is the ordering
    defect next to the radius force, not a false historical result.
    """
    metric = _metric()
    shift = _shift(metric)
    sea = ground_covariance(direct_hamiltonian(metric, shift, 1))[2]
    comparison = ordering_comparison(metric, shift, sea, 1, amplitude=0.12)
    assert comparison["same_physical_covariance"] is True
    assert comparison["historical_results_marked_false"] is False
    assert comparison["direct_fixed_LQ_trace_change"] == pytest.approx(0, abs=1e-12)
    assert comparison["direct_ward_l2"] == pytest.approx(0, abs=1e-12)
    assert comparison["product_rule_fixed_LQ_trace_change"] == pytest.approx(
        -0.0461283850516527, rel=1e-9, abs=1e-12,
    )
    assert comparison["product_rule_ward_l2"] == pytest.approx(
        0.13566530208932837, rel=1e-9, abs=1e-12,
    )
    assert comparison["physical_radius_force_l2"] == pytest.approx(
        0.6496334743012157, rel=1e-9, abs=1e-12,
    )
    assert comparison["radius_force_disagreement_l2"] < 1e-12
