"""Independent jets, horizon roots and the read-only RN calibration."""
import numpy as np
import pytest

from recursive_horizons import nsc_rn_reference as rn


A = rn.AUDITED_A
C_F = rn.AUDITED_C_F
FLUX = rn.AUDITED_FLUX


def _metric_at(mass, charge, radius, theta):
    beta = np.sqrt(2.0 * mass / radius - charge / radius ** 2)
    factor = 1.0 - beta ** 2
    metric = np.zeros((4, 4))
    metric[0, 0] = factor
    metric[0, 1] = metric[1, 0] = -beta
    metric[1, 1] = -1.0
    metric[2, 2] = -radius ** 2
    metric[3, 3] = -radius ** 2 * np.sin(theta) ** 2
    return metric


def _metric_derivatives(mass, charge, radius, theta):
    beta = np.sqrt(2.0 * mass / radius - charge / radius ** 2)
    beta_r = (-mass / radius ** 2 + charge / radius ** 3) / beta
    factor_r = 2.0 * mass / radius ** 2 - 2.0 * charge / radius ** 3
    derivatives = [np.zeros((4, 4)) for _ in range(4)]
    derivatives[1][0, 0] = factor_r
    derivatives[1][0, 1] = derivatives[1][1, 0] = -beta_r
    derivatives[1][2, 2] = -2.0 * radius
    derivatives[1][3, 3] = -2.0 * radius * np.sin(theta) ** 2
    derivatives[2][3, 3] = -2.0 * radius ** 2 * np.sin(theta) * np.cos(theta)
    return derivatives


def _christoffel(mass, charge, radius, theta):
    metric = _metric_at(mass, charge, radius, theta)
    inverse = np.linalg.inv(metric)
    derivatives = _metric_derivatives(mass, charge, radius, theta)
    gamma = np.zeros((4, 4, 4))
    for lam in range(4):
        for mu in range(4):
            for nu in range(4):
                accumulator = 0.0
                for sig in range(4):
                    accumulator += inverse[lam, sig] * (
                        derivatives[mu][nu, sig] + derivatives[nu][mu, sig] - derivatives[sig][mu, nu]
                    )
                gamma[lam, mu, nu] = 0.5 * accumulator
    return gamma


def _riemann(mass, charge, radius, theta, step=1e-6):
    """Finite-difference Riemann of the ingoing chart, independent of the closed form."""
    gamma = _christoffel(mass, charge, radius, theta)
    partial = [np.zeros((4, 4, 4)) for _ in range(4)]
    partial[1] = (
        _christoffel(mass, charge, radius + step, theta) - _christoffel(mass, charge, radius - step, theta)
    ) / (2.0 * step)
    partial[2] = (
        _christoffel(mass, charge, radius, theta + step) - _christoffel(mass, charge, radius, theta - step)
    ) / (2.0 * step)
    riemann = np.zeros((4, 4, 4, 4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    quadratic = 0.0
                    for e in range(4):
                        quadratic += gamma[a, c, e] * gamma[e, d, b] - gamma[a, d, e] * gamma[e, c, b]
                    riemann[a, b, c, d] = partial[c][a, d, b] - partial[d][a, c, b] + quadratic
    return riemann


def _independent_invariants(mass, charge, radius, theta=0.7):
    gamma_metric = _metric_at(mass, charge, radius, theta)
    inverse = np.linalg.inv(gamma_metric)
    riemann = _riemann(mass, charge, radius, theta)
    ricci = np.zeros((4, 4))
    for b in range(4):
        for d in range(4):
            ricci[b, d] = sum(riemann[a, b, a, d] for a in range(4))
    scalar = sum(inverse[b, d] * ricci[b, d] for b in range(4) for d in range(4))
    ricci_up = np.zeros((4, 4))
    for b in range(4):
        for d in range(4):
            ricci_up[b, d] = sum(inverse[b, s] * inverse[d, t] * ricci[s, t] for s in range(4) for t in range(4))
    ricci2 = sum(ricci[b, d] * ricci_up[b, d] for b in range(4) for d in range(4))
    lowered = np.zeros_like(riemann)
    for e in range(4):
        lowered[e] = sum(gamma_metric[e, a] * riemann[a] for a in range(4))
    kretschmann = 0.0
    for e in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    kretschmann += lowered[e, b, c, d] * sum(
                        inverse[e, ee] * inverse[b, bb] * inverse[c, cc] * inverse[d, dd]
                        * lowered[ee, bb, cc, dd]
                        for ee in range(4) for bb in range(4) for cc in range(4) for dd in range(4)
                    )
    beta = np.sqrt(2.0 * mass / radius - charge / radius ** 2)
    frame = [
        np.array([1.0, -beta, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0 / radius, 0.0]),
        np.array([0.0, 0.0, 0.0, 1.0 / (radius * np.sin(theta))]),
    ]

    def component(i, j, k, l):
        total = 0.0
        for e in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        total += lowered[e, b, c, d] * frame[i][e] * frame[j][b] * frame[k][c] * frame[l][d]
        return total

    return {
        "R4": scalar,
        "Ricci2": ricci2,
        "K": kretschmann,
        "R_0101": component(0, 1, 0, 1),
        "R_0202": component(0, 2, 0, 2),
    }


def _bisect_root(function, left, right):
    flo, fhi = function(left), function(right)
    if flo * fhi > 0.0:
        raise AssertionError("horizon bracket does not change sign")
    for _ in range(80):
        mid = 0.5 * (left + right)
        if function(left) * function(mid) <= 0.0:
            right = mid
        else:
            left = mid
    return 0.5 * (left + right)


def test_action_dictionary_and_vacuum_lapse():
    reference = rn.RNReference.from_action(A, C_F, FLUX, 2.0)
    assert reference.G_N == pytest.approx(1.0 / (16.0 * np.pi * A))
    assert reference.magnetic_r2 == pytest.approx(C_F * FLUX ** 2 / (4.0 * A))
    assert rn.ACTION["weyl_pole"] is False
    assert rn.ACTION["dirac"]["multiplicity"] == 4
    assert rn.ACTION["dirac"]["multiplicity_counted_once"] is True
    assert rn.ACTION["k_is_mass"] is False
    jets = reference.pg_metric_jets(6.0)
    assert jets["N"] == pytest.approx(1.0)
    assert jets["N_r"] == pytest.approx(0.0)
    assert jets["coupled_lapse_kept_at_one"] is False
    assert jets["g_tt"] == pytest.approx(reference.f(6.0))
    assert jets["g_tr"] == pytest.approx(-reference.beta_rn(6.0))
    step = 1e-6
    radius = 6.0
    numeric = (reference.beta_rn(radius + step) - reference.beta_rn(radius - step)) / (2.0 * step)
    assert reference.beta_derivative(radius) == pytest.approx(numeric, rel=1e-8)


def test_independent_tetrad_jets_at_fixed_charge():
    charge = C_F * FLUX ** 2 / (4.0 * A)
    for mass, radius in ((2.0, 8.0), (1.3, 6.0)):
        reference = rn.RNReference.from_action(A, C_F, FLUX, mass)
        assert reference.V4 == 0.0
        assert reference.magnetic_r2 == pytest.approx(charge)
        sample = rn.sourcefree_static_rn(radius, mass)
        assert sample["field_state"] is None
        assert sample["campaign"] is None
        assert sample["P2"] == pytest.approx(charge)
        assert sample["br_radius_rescaled"] is False
        assert sample["charged_mass"] == pytest.approx(mass, rel=1e-12)
        closed = reference.invariants(radius)
        tetrad = reference.tetrad_riemann(radius)
        independent = _independent_invariants(mass, charge, radius)
        assert independent["R4"] == pytest.approx(0.0, abs=1e-8)
        assert independent["Ricci2"] == pytest.approx(closed["Ricci2"], rel=1e-7)
        assert independent["K"] == pytest.approx(closed["K"], rel=1e-7)
        assert independent["R_0101"] == pytest.approx(tetrad["R_0101"], rel=1e-6)
        assert independent["R_0202"] == pytest.approx(tetrad["R_0202"], rel=1e-6)
    throat_mass = float(np.sqrt(charge))
    throat = rn.RNReference.from_action(A, C_F, FLUX, throat_mass)
    values = throat.invariants(throat_mass)
    assert values["Ricci2"] == pytest.approx(4.0 / throat_mass ** 4)
    assert values["K"] == pytest.approx(8.0 / throat_mass ** 4)


def test_horizon_roots_at_fixed_charge():
    charge = C_F * FLUX ** 2 / (4.0 * A)
    for mass in (2.0, 1.3):
        reference = rn.RNReference.from_action(A, C_F, FLUX, mass)
        assert reference.magnetic_r2 == pytest.approx(charge)
        horizons = reference.horizons

        def factor(radius, mass=mass, charge=charge):
            return 1.0 - 2.0 * mass / radius + charge / radius ** 2

        outer = _bisect_root(factor, 0.5 * (horizons["r_minus"] + horizons["r_plus"]), 4.0 * mass)
        inner = _bisect_root(factor, reference.shift_domain_lower * 1.0000001, outer * 0.999)
        assert horizons["r_plus"] == pytest.approx(outer, rel=1e-10)
        assert horizons["r_minus"] == pytest.approx(inner, rel=1e-10)
        assert reference.f(horizons["r_plus"]) == pytest.approx(0.0, abs=1e-12)
    assert rn.RNReference.from_action(A, C_F, FLUX, 2.0).horizons["r_plus"] != pytest.approx(
        rn.RNReference.from_action(A, C_F, FLUX, 1.3).horizons["r_plus"]
    )
    assert rn.RNReference.from_action(A, C_F, FLUX, 0.4).horizons["black_hole"] is False


def test_excision_is_pure_outflow_and_domains_differ():
    reference = rn.RNReference.from_action(A, C_F, FLUX, 2.0)
    horizons = reference.horizons
    excision = 0.5 * (horizons["r_minus"] + horizons["r_plus"])
    outgoing, ingoing = reference.characteristic_speeds(excision)
    assert outgoing < 0.0 and ingoing < 0.0
    report = rn.placement(
        reference,
        parent=(8.0, 30.0),
        child=(horizons["r_plus"] + 0.05, horizons["r_plus"] + 1.0),
        excision=excision,
        outer="absorbing",
    )
    assert report["supported"] is True
    assert report["periodic"] is False
    same = rn.placement(
        reference, parent=(8.0, 30.0), child=(8.0, 30.0), excision=excision, outer="absorbing",
    )
    assert same["supported"] is False
    periodic = rn.placement(
        reference,
        parent=(8.0, 30.0),
        child=(horizons["r_plus"] + 0.05, horizons["r_plus"] + 1.0),
        excision=excision,
        outer="periodic",
    )
    assert periodic["supported"] is False


def test_manufactured_quadrature_does_not_claim_positive_energy():
    shell = rn.manufactured_quadrature_control(
        center=6.0, width=0.4, r_inner=4.0, r_outer=9.0, points=21,
        occupations=np.array([0.25, 0.5]),
    )
    assert shell["phi"].shape == (2, 21, 2)
    assert np.allclose(shell["sigma2"], np.diag([1.0, -1.0]))
    assert np.allclose(shell["gram"], np.eye(2))
    assert shell["weights"].sum() == pytest.approx(5.0)
    assert shell["positivity_unestablished"] is True
    assert shell["prepared_positive_energy_source"] is False
    assert shell["exterior_killing_frequency"] is None
    assert shell["angular_is_4d_mass"] is False


def test_decode_ignores_forged_k_and_gradient_uses_r_to_the_fourth():
    count = 4
    spacing = 0.2
    frame = np.eye(count)
    frame[0, 1] = 0.3
    nodal = np.array([0.2, -0.4, 0.15, 0.05])
    stored = spacing * (frame.T @ nodal)
    record = {
        "momentum_representation": "canonical_pi",
        "dx_g": spacing,
        "common_k": 0.02,
        "nf": count + 1,
    }
    arrays = {"W": frame, "pi_Q": stored, "pi_r": stored.copy()}
    decoded = rn.decode_canonical_pi(record, arrays)
    assert np.allclose(decoded["p_Q"], frame @ (stored / spacing))
    forged = dict(record, common_k=80.0)
    again = rn.decode_canonical_pi(forged, arrays)
    assert np.allclose(again["p_Q"], decoded["p_Q"])
    assert np.allclose(again["p_r"], decoded["p_r"])
    missing = dict(record)
    del missing["momentum_representation"]
    with pytest.raises(ValueError):
        rn.decode_canonical_pi(missing, arrays)
    radius = 2.5
    lapse = 1.7
    slope = 0.15
    momentum = 0.4
    coefficient = -8.0 * np.pi * A
    rate = lapse * momentum / (2.0 * coefficient * radius)
    expected = (rate ** 2 - slope ** 2) / (radius ** 2 * lapse ** 2)
    assert rn.conformal_grad_r_squared(lapse, radius, slope, momentum, A) == pytest.approx(expected)
    stale = momentum ** 2 / (4.0 * coefficient ** 2 * radius ** 2) - slope ** 2 / (radius ** 2 * lapse ** 2)
    assert expected != pytest.approx(stale)
    charge = C_F * FLUX ** 2 / (4.0 * A)
    charged = rn.charged_mass_from_gradient(radius, expected, charge)
    bare = rn.misner_sharp(radius, expected)
    assert charged - bare == pytest.approx(charge / (2.0 * radius))
    assert charged != pytest.approx(0.02)


def test_historical_decode_and_unresolved_curvature():
    evidence = rn.assess_historical()
    assert evidence["authenticated"] is True
    assert evidence["final_record_written"] is False
    assert evidence["coupled_calibration"] is False
    assert evidence["audited_coefficients"]["V4"] == 0.0
    strong = evidence["campaigns"]["strong-empty-k-v1"]["plus"]
    assert strong["n"] == 127
    assert strong["dx_g"] == pytest.approx(8.0 / 127.0)
    assert strong["curvature"]["physical_ground_truth"] is False
    assert strong["curvature"]["status"] == "unresolved"
    assert strong["mass"]["fitted_to_kinetic_anchor"] is False
    assert strong["mass"]["kinetic_anchor_from_decoded_p"] == pytest.approx(strong["declared_k"], rel=1e-12)
    assert strong["mass"]["charged_mass"]["mean"] != pytest.approx(strong["declared_k"])
    assert strong["stored_constraints"]["raw_C_max"] < 1e-8
    weak = evidence["campaigns"]["old-weak-control"]
    assert weak["mass"]["charged_mass"]["mean"] != pytest.approx(weak["declared_k"])
    minus = evidence["campaigns"]["strong-empty-k-v1"]["minus"]
    assert minus["mass"]["charged_mass"]["mean"] == pytest.approx(strong["mass"]["charged_mass"]["mean"], rel=1e-12)
    confirmation = evidence["campaigns"]["cut-confirmation-v1"]
    assert confirmation["complete"] is True
    assert confirmation["momentum_decoded"] is False
    assert confirmation["curvature"]["physical_ground_truth"] is False
    vacuum = evidence["campaigns"]["vacuum-br-control"]
    assert vacuum["used_as_rn_test"] is False
    assert vacuum["authenticated"] is True
    assert any("0.757" in item and "0.167" in item for item in evidence["unresolved_jets"])
    assert any("k=0.02" in item for item in evidence["unresolved_jets"])
    assert not rn.ASSESSMENT_OUTPUT.exists()


def test_evidence_writer_refuses_every_path_except_the_new_record(tmp_path):
    payload = {
        "schema": rn.ASSESSMENT_SCHEMA,
        "input_hashes": {"results/development/example.npz": "abc"},
    }
    with pytest.raises(RuntimeError):
        rn.write_assessment(payload, tmp_path / "assessment.json", producer_commit=None)
    with pytest.raises(RuntimeError):
        rn.write_assessment(payload, tmp_path / "assessment.json", producer_commit=rn.CUT_CONFIRMATION_COMMIT)
    sealed = rn.SEALED_INPUT_DIRECTORIES[0] / "forged.json"
    with pytest.raises(RuntimeError):
        rn.write_assessment(payload, sealed, producer_commit=rn.CUT_CONFIRMATION_COMMIT)
    with pytest.raises(RuntimeError):
        rn.write_assessment(
            {"schema": rn.ASSESSMENT_SCHEMA},
            rn.ASSESSMENT_OUTPUT,
            producer_commit=rn.CUT_CONFIRMATION_COMMIT,
        )
    limit = rn.MAX_ASSESSMENT_BYTES
    rn.MAX_ASSESSMENT_BYTES = 32
    try:
        with pytest.raises(RuntimeError):
            rn.write_assessment(payload, rn.ASSESSMENT_OUTPUT, producer_commit=rn.CUT_CONFIRMATION_COMMIT)
    finally:
        rn.MAX_ASSESSMENT_BYTES = limit
    assert not rn.ASSESSMENT_OUTPUT.exists()
    assert not sealed.exists()
