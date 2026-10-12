"""Synthetic algorithm controls for the local history Newton CONTROL owner.

Fixtures are labelled synthetic evaluators. They are not NSC source
campaigns, frozen physical certificates, or production field evolution.
"""
import numpy as np
import pytest

from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots
from recursive_horizons.nsc_local_history_newton import (
    ALLOWED_COEFFICIENT_COUNTS, DeclaredJetBoundary, GeometryBoundEvaluatorCache,
    HistoryCollocation, HistoryEvaluatorReceipt, LocalHistoryNewtonSettings,
    NUMERICAL_CONTROL_TOLERANCE, STOP_BUDGET, STOP_CONVERGED, STOP_LINE_SEARCH,
    STOP_RANK, STOP_VERIFICATION, bind_history_evaluator_receipt, solve_local_history,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


SYNTHETIC_SOURCE = 'synthetic-evaluator-fixture:coupled-nonlinear'


def star_coefficients(n=8):
    coeff = np.zeros((2, n))
    coeff[0, :4] = [0.012, -0.006, 0.003, -0.001]
    coeff[1, :4] = [0.020, 0.008, -0.004, 0.002]
    return coeff


def make_family(coeff, n=None):
    coeff = np.asarray(coeff, float)
    if coeff.ndim == 1:
        n = 8 if n is None else n
        coeff = coeff.reshape(2, n)
    return LocalIncomingFamily(coeff)


def family_history_identity(family):
    coeff = np.asarray(family.coefficients, float)
    return (
        tuple(coeff[0].tolist()), tuple(coeff[1].tolist()),
        float(family.center), float(family.axial_inner), float(family.axial_outer),
        float(family.normal_inner), float(family.normal_outer),
    )


def synthetic_contract(family, source_identity):
    return {
        'source_identity': source_identity,
        'history_identity': family_history_identity(family),
        'full_retarded_state_derivative': True,
        'constant_frozen_matter_used_as_physical_state': False,
        'reference_state_tangent_subtracted': False,
    }


def declared_jet(family, point=None):
    point = family.center if point is None else float(point)
    w, _ = family.functions
    return DeclaredJetBoundary(point, w(point, 0), w(point, 1), w(point, 2))


def _map(slots, kappa):
    w, wz, wzz, wzzz, U, Uz = slots.T
    n_component = U + wzz + w * U + 0.25 * w * w + kappa * w[0] * U
    beta = Uz + wzzz + w * wz
    return np.stack((n_component, beta), axis=-1)


def _tangent(slots, tangents, kappa):
    w, wz, _, _, U, _ = slots.T
    dw, dwz, dwzz, dwzzz, dU, dUz = (tangents[..., k] for k in range(6))
    dw_left = tangents[:, 0, 0][:, None]
    dN = dU + dwzz + dw * U + w * dU + 0.5 * w * dw + kappa * (dw_left * U + w[0] * dU)
    dbeta = dUz + dwzzz + dw * wz + w * dwz
    return np.stack((dN, dbeta), axis=-1)


def _geom_tangent(tangents):
    dN = tangents[..., 4] + tangents[..., 2]
    dbeta = tangents[..., 5] + tangents[..., 3]
    return np.stack((dN, dbeta), axis=-1)


class CoupledNonlinearSyntheticEvaluator:
    """SYNTHETIC coupled nonlinear residual. Not a physical NSC source."""

    def __init__(self, star, kappa=0.2, source_identity=SYNTHETIC_SOURCE):
        self.star = star
        self.kappa = kappa
        self.source_identity = source_identity
        self.calls = 0
        self.geometries = []

    def evaluate(self, family, z=None):
        self.calls += 1
        self.geometries.append(tuple(np.asarray(family.coefficients, float).ravel()))
        z = np.asarray(family.collocation_nodes(len(family.coefficients[0])), float) if z is None else np.asarray(z, float)
        n = len(family.coefficients[0])
        ndir = 2 * n
        slots, tangents = compatible_history_slots(family.metric(), z, ndir)
        star_slots, _ = compatible_history_slots(self.star.metric(), z, ndir)
        return {
            'constraint_order': ('N', 'beta'),
            'z': z,
            'action_gradient': _map(slots, self.kappa) - _map(star_slots, self.kappa),
            'history_jacobian': _tangent(slots, tangents, self.kappa),
            'local_reference_gradient_tangent': _geom_tangent(tangents),
            'incoming_slots': slots,
            'incoming_slot_tangents': tangents,
            'error_budget': None,
            'physical_constraint_status': 'OPEN',
            'physical_EXISTENCE_certificate': False,
            'fixture': 'synthetic-evaluator-fixture:coupled-nonlinear',
            **synthetic_contract(family, self.source_identity),
        }


class ZeroJacobianSyntheticEvaluator(CoupledNonlinearSyntheticEvaluator):
    """SYNTHETIC rank-deficient Jacobian. Not a non-existence certificate."""

    def evaluate(self, family, z=None):
        raw = super().evaluate(family, z)
        raw = dict(raw)
        raw['action_gradient'] = np.ones_like(raw['action_gradient'])
        raw['history_jacobian'] = np.zeros_like(raw['history_jacobian'])
        return raw


class InterpolatingSyntheticEvaluator:
    """SYNTHETIC nodal interpolation of (w,U). Used to allow n=16 without the
    mixed-derivative collocation ill-conditioning of the coupled fixture.
    """

    def __init__(self, star, source_identity='synthetic-evaluator-fixture:interpolating'):
        self.star = star
        self.source_identity = source_identity
        self.calls = 0

    def evaluate(self, family, z=None):
        self.calls += 1
        z = np.asarray(family.collocation_nodes(len(family.coefficients[0])), float) if z is None else np.asarray(z, float)
        n = len(family.coefficients[0])
        ndir = 2 * n
        slots, tangents = compatible_history_slots(family.metric(), z, ndir)
        star_slots, _ = compatible_history_slots(self.star.metric(), z, ndir)
        gradient = np.stack((slots[:, 0] - star_slots[:, 0], slots[:, 4] - star_slots[:, 4]), axis=-1)
        jacobian = np.stack((tangents[..., 0], tangents[..., 4]), axis=-1)
        return {
            'constraint_order': ('N', 'beta'),
            'z': z,
            'action_gradient': gradient,
            'history_jacobian': jacobian,
            'incoming_slots': slots,
            'incoming_slot_tangents': tangents,
            'error_budget': None,
            'physical_constraint_status': 'OPEN',
            'fixture': 'synthetic-evaluator-fixture:interpolating',
            **synthetic_contract(family, self.source_identity),
        }


class AmplitudeTargetSyntheticEvaluator:
    """SYNTHETIC radius-wall driver. Not physical matching data.

    N = w - target, beta = U, with a full retarded slot Jacobian. The
    target amplitude is outside the positive-radius ball.
    """

    def __init__(self, target_w0=80.0, source_identity='synthetic-evaluator-fixture:radius-wall'):
        self.target_w0 = target_w0
        self.source_identity = source_identity
        self.calls = 0

    def evaluate(self, family, z=None):
        self.calls += 1
        z = np.asarray(family.collocation_nodes(len(family.coefficients[0])), float) if z is None else np.asarray(z, float)
        n = len(family.coefficients[0])
        ndir = 2 * n
        slots, tangents = compatible_history_slots(family.metric(), z, ndir)
        gradient = np.stack((slots[:, 0] - self.target_w0, slots[:, 4]), axis=-1)
        jacobian = np.stack((tangents[..., 0], tangents[..., 4]), axis=-1)
        return {
            'constraint_order': ('N', 'beta'),
            'z': z,
            'action_gradient': gradient,
            'history_jacobian': jacobian,
            'incoming_slots': slots,
            'incoming_slot_tangents': tangents,
            'error_budget': None,
            'physical_constraint_status': 'OPEN',
            'fixture': 'synthetic-evaluator-fixture:radius-wall',
            **synthetic_contract(family, self.source_identity),
        }


def _assert_open(result):
    assert result['physical_constraint_status'] == 'OPEN'
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['physical_NONEXISTENCE_certificate'] is False
    assert result['stop_reason_is_nonexistence'] is False
    assert result['optimization_failure_is_nonexistence'] is False
    assert result['certified_constraint_residual'] is None
    assert result['full_source_error_bound'] is None
    assert result['missing_error_budget'] is None
    assert result['verification']['between_node_certified'] is False
    assert result['verification']['automatic_resolution_increase'] is False
    assert result['scope']['physical_constraint_status'] == 'OPEN'
    if result['numerical_control_converged']:
        assert result['stop_reason'] == STOP_CONVERGED
        assert result['control_converged'] is True
        assert result['solve_all_nodes_pass'] is True
        assert result['verification_nodes_pass'] is True
    for name, bound in result['error_budget'].items():
        assert bound is None, name


def test_known_coupled_nonlinear_solution_converges_numerical_control_open():
    star = make_family(star_coefficients())
    start = star_coefficients()
    start[0, 0] += 0.004
    start[1, 1] -= 0.003
    start[0, 2] += 0.002
    evaluator = CoupledNonlinearSyntheticEvaluator(star)
    result = solve_local_history(
        evaluator, make_family(start), boundary=declared_jet(star),
        collocation=HistoryCollocation(8, 'tau'))
    _assert_open(result)
    assert result['stop_reason'] in (STOP_CONVERGED, STOP_VERIFICATION)
    assert result['numerical_control_converged'] is (result['stop_reason'] == STOP_CONVERGED)
    assert result['control_converged'] is True
    assert result['control_residual_max'] <= NUMERICAL_CONTROL_TOLERANCE
    assert result['newton_jacobian_source'] == 'full_retarded_history_jacobian'
    assert result['geometric_jacobian_used_as_newton_matrix'] is False
    assert result['full_retarded_jacobian_consumed'] is True
    assert result['full_minus_geometric_jacobian_max'] > 1e-8
    np.testing.assert_allclose(result['coefficients'], star_coefficients(), atol=1e-6, rtol=0)
    assert result['boundary_control']['physical_matching_data'] is False
    assert result['collocation']['layout'] == 'tau'
    assert result['collocation']['equations'] == {'N': 8, 'beta': 5, 'boundary': 3, 'unknowns': 16}
    assert result['verification']['nodes_separate_from_solve_nodes'] is True
    assert not np.array_equal(result['z'], result['verification']['z'])
    assert result['verification']['between_node_certified'] is False
    assert result['history_jacobian'].shape == (16, 8, 2)
    assert evaluator.calls >= 2


def test_rectangular_layout_recovers_known_solution_and_keeps_all_beta_nodes():
    star = make_family(star_coefficients())
    start = star_coefficients()
    start[0, 1] += 0.003
    evaluator = CoupledNonlinearSyntheticEvaluator(star)
    result = solve_local_history(
        evaluator, make_family(start), boundary=declared_jet(star),
        collocation=HistoryCollocation(8, 'rectangular'))
    _assert_open(result)
    assert result['stop_reason'] == STOP_CONVERGED
    assert result['collocation']['layout'] == 'rectangular'
    assert result['collocation']['equations'] == {'N': 8, 'beta': 8, 'boundary': 3, 'unknowns': 16}
    assert result['linear_algebra']['n_equations'] == 19
    np.testing.assert_allclose(result['coefficients'], star_coefficients(), atol=1e-6, rtol=0)
    assert 'retarded' in result['collocation']['rectangular_justification']


def test_dropped_or_missing_retarded_jacobian_rejected():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    slots, tangents = compatible_history_slots(star.metric(), z, 16)
    base = CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z)

    missing = dict(base)
    missing.pop('history_jacobian')
    with pytest.raises(ValueError, match='history_jacobian'):
        bind_history_evaluator_receipt(missing, star, expected_nodes=z)

    none = dict(base)
    none['history_jacobian'] = None
    with pytest.raises(ValueError, match='history_jacobian'):
        bind_history_evaluator_receipt(none, star, expected_nodes=z)

    dropped = dict(base)
    dropped['history_jacobian'] = np.zeros((len(z), 2, 6))
    with pytest.raises(ValueError, match='ndir, nz, 2'):
        bind_history_evaluator_receipt(dropped, star, expected_nodes=z)

    frozen = dict(base)
    frozen['incoming_C0'] = np.zeros((2, 2))
    with pytest.raises(ValueError, match='frozen C0'):
        bind_history_evaluator_receipt(frozen, star, expected_nodes=z)
    assert slots.shape == (len(z), 6) and tangents.shape[0] == 16


def test_source_identity_change_rejected():
    star = make_family(star_coefficients())
    start = star_coefficients()
    start[0, 0] += 0.004
    evaluator = CoupledNonlinearSyntheticEvaluator(star)

    class FlipSource:
        def __init__(self):
            self.inner = evaluator
            self.calls = 0

        def evaluate(self, family, z=None):
            self.calls += 1
            raw = dict(self.inner.evaluate(family, z))
            raw['source_identity'] = 'synthetic-source-A' if self.calls == 1 else 'synthetic-source-B'
            return raw

    with pytest.raises(ValueError, match='source identity'):
        solve_local_history(
            FlipSource(), make_family(start), boundary=declared_jet(star),
            collocation=HistoryCollocation(8, 'tau'))


def test_radius_violating_step_is_damped():
    start = np.zeros((2, 8))
    start[0, 0] = 32.6
    family = make_family(start)
    assert family.radius_lower_bound() > 0
    evaluator = AmplitudeTargetSyntheticEvaluator(target_w0=80.0)
    boundary = DeclaredJetBoundary(family.center, 80.0, 0.0, 0.0)
    result = solve_local_history(
        evaluator, family, boundary=boundary,
        collocation=HistoryCollocation(8, 'rectangular'),
        settings=LocalHistoryNewtonSettings(
            max_iterations=8, max_step_norm=8.0, min_line_search_scale=1 / 64))
    _assert_open(result)
    assert result['stop_reason'] in (STOP_LINE_SEARCH, STOP_BUDGET, STOP_RANK)
    assert result['stop_reason'] != STOP_CONVERGED
    assert any(step['reason'] == 'nonpositive_radius' for step in result['rejected_steps'])
    assert result['radius_lower_bound'] > 0


def test_rank_deficiency_reported_open():
    star = make_family(star_coefficients())
    evaluator = ZeroJacobianSyntheticEvaluator(star)
    result = solve_local_history(
        evaluator, star, boundary=declared_jet(star),
        collocation=HistoryCollocation(8, 'tau'))
    _assert_open(result)
    assert result['stop_reason'] == STOP_RANK
    assert result['numerical_control_converged'] is False
    assert result['linear_algebra']['rank_or_conditioning'] is True
    assert result['linear_algebra']['rank'] < result['linear_algebra']['n_unknowns']


def test_iteration_budget_stop_open():
    star = make_family(star_coefficients())
    start = star_coefficients()
    start[0, 0] += 0.004
    evaluator = CoupledNonlinearSyntheticEvaluator(star)
    result = solve_local_history(
        evaluator, make_family(start), boundary=declared_jet(star),
        collocation=HistoryCollocation(8, 'tau'),
        settings=LocalHistoryNewtonSettings(max_iterations=0))
    _assert_open(result)
    assert result['stop_reason'] == STOP_BUDGET
    assert result['numerical_control_converged'] is False
    assert result['iterations'] == 0


def test_mismatched_nodes_or_history_rejected():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    evaluator = CoupledNonlinearSyntheticEvaluator(star)
    raw = evaluator.evaluate(star, z)

    shifted = dict(raw)
    shifted['z'] = z + 0.01
    with pytest.raises(ValueError, match='mismatched incoming nodes'):
        bind_history_evaluator_receipt(shifted, star, expected_nodes=z)

    other = make_family(star_coefficients() + 0.002)
    other_slots, _ = compatible_history_slots(other.metric(), z, 16)
    mismatched = dict(raw)
    mismatched['incoming_slots'] = other_slots
    with pytest.raises(ValueError, match='mismatched history'):
        bind_history_evaluator_receipt(mismatched, star, expected_nodes=z)

    class WrongNodes:
        def evaluate(self, family, z=None):
            payload = dict(evaluator.evaluate(family, z))
            payload['z'] = np.linspace(family.interval[0], family.interval[1], 9)
            payload['action_gradient'] = np.zeros((9, 2))
            payload['history_jacobian'] = np.zeros((16, 9, 2))
            payload.pop('incoming_slots', None)
            payload.pop('incoming_slot_tangents', None)
            return payload

    with pytest.raises(ValueError, match='mismatched incoming nodes'):
        solve_local_history(
            WrongNodes(), star, boundary=declared_jet(star),
            collocation=HistoryCollocation(8, 'tau'))


def test_boundary_required_and_recorded_as_solver_control():
    star = make_family(star_coefficients())
    with pytest.raises(TypeError, match='declared w'):
        solve_local_history(
            CoupledNonlinearSyntheticEvaluator(star), star, boundary=None)
    result = solve_local_history(
        CoupledNonlinearSyntheticEvaluator(star), star, boundary=declared_jet(star))
    _assert_open(result)
    assert result['boundary_control']['role'] == 'solver_control'
    assert result['boundary_control']['physical_matching_data'] is False
    assert result['stop_reason'] == STOP_CONVERGED


def test_caller_n16_allowed_no_automatic_resolution():
    star = make_family(star_coefficients(16))
    start = star_coefficients(16)
    start[0, 0] += 0.002
    result = solve_local_history(
        InterpolatingSyntheticEvaluator(star), make_family(start),
        boundary=declared_jet(star), collocation=HistoryCollocation(16, 'rectangular'))
    _assert_open(result)
    assert result['stop_reason'] == STOP_CONVERGED
    assert result['collocation']['coefficient_count'] == 16
    assert result['collocation']['automatic_resolution'] is False
    assert result['history_jacobian'].shape == (32, 16, 2)
    with pytest.raises(ValueError, match='not raised automatically'):
        solve_local_history(
            CoupledNonlinearSyntheticEvaluator(star), make_family(star_coefficients()),
            boundary=declared_jet(make_family(star_coefficients())),
            collocation=HistoryCollocation(16, 'tau'))
    with pytest.raises(ValueError, match='8, 16 or 32'):
        HistoryCollocation(64, 'tau')
    assert ALLOWED_COEFFICIENT_COUNTS == (8, 16, 32)


def test_geometry_cache_is_bound_to_coefficients_and_nodes():
    star = make_family(star_coefficients())
    evaluator = CoupledNonlinearSyntheticEvaluator(star)
    cache = GeometryBoundEvaluatorCache(evaluator)
    z = star.collocation_nodes(8)
    cache.evaluate(star, z)
    cache.evaluate(star, z)
    assert evaluator.calls == 1
    perturbed = make_family(star_coefficients() + 0.001)
    cache.evaluate(perturbed, z)
    assert evaluator.calls == 2
    cache.evaluate(star, star.collocation_nodes(17))
    assert evaluator.calls == 3


def test_receipt_rejects_physical_existence_claim():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    raw = CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z)
    raw = dict(raw)
    raw['physical_EXISTENCE_certificate'] = True
    with pytest.raises(ValueError, match='EXISTENCE'):
        bind_history_evaluator_receipt(raw, star, expected_nodes=z)
    typed = bind_history_evaluator_receipt(
        CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z), star, expected_nodes=z)
    assert isinstance(typed, HistoryEvaluatorReceipt)
    assert typed.physical_constraint_status == 'OPEN'
    assert typed.missing_error_budget is None
    assert typed.full_retarded_state_derivative is True


class HeldOutBetaSyntheticEvaluator:
    """SYNTHETIC tau residual that vanishes on selected beta rows only."""

    def __init__(self, source_identity='synthetic-evaluator-fixture:held-out-beta'):
        self.source_identity = source_identity
        self.calls = 0

    def evaluate(self, family, z=None):
        self.calls += 1
        z = np.asarray(family.collocation_nodes(len(family.coefficients[0])), float) if z is None else np.asarray(z, float)
        n = len(family.coefficients[0])
        ndir = 2 * n
        slots, tangents = compatible_history_slots(family.metric(), z, ndir)
        gradient = np.zeros((len(z), 2))
        held = np.array([0, n - 2, n - 1], dtype=int)
        gradient[held, 1] = 1e-6
        jacobian = np.stack((tangents[..., 0], tangents[..., 4]), axis=-1)
        return {
            'constraint_order': ('N', 'beta'),
            'z': z,
            'action_gradient': gradient,
            'history_jacobian': jacobian,
            'incoming_slots': slots,
            'incoming_slot_tangents': tangents,
            'error_budget': None,
            'physical_constraint_status': 'OPEN',
            'fixture': 'synthetic-evaluator-fixture:held-out-beta',
            **synthetic_contract(family, self.source_identity),
        }


def test_prep_digest_change_rejected_when_source_covariance_unchanged():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    base = CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z)

    def with_batches(preparation):
        raw = dict(base)
        raw.pop('source_identity', None)
        raw['batch_records'] = {
            '14_1': {'source_digest': 'csrc-fixed', 'preparation_digest': preparation},
            '14_-1': {'source_digest': 'csrc-fixed', 'preparation_digest': preparation + '-neg'},
        }
        return raw

    first = bind_history_evaluator_receipt(with_batches('A_up-rho-mass-v1'), star, expected_nodes=z)
    with pytest.raises(ValueError, match='source identity'):
        bind_history_evaluator_receipt(
            with_batches('A_up-rho-mass-v2'), star, expected_nodes=z,
            expected_source_identity=first.source_identity)
    missing = dict(base)
    missing.pop('source_identity', None)
    missing['batch_records'] = {'14_1': {'source_digest': 'csrc-fixed'}}
    with pytest.raises(ValueError, match='partial missing batch bindings'):
        bind_history_evaluator_receipt(missing, star, expected_nodes=z)
    mutable = {'source_digest': 'csrc-fixed', 'preparation_digest': 'A_up-v1'}
    aliased = dict(base)
    aliased['source_identity'] = mutable
    bound = bind_history_evaluator_receipt(aliased, star, expected_nodes=z)
    mutable['preparation_digest'] = 'A_up-mutated'
    mutated = dict(base)
    mutated['source_identity'] = mutable
    with pytest.raises(ValueError, match='source identity'):
        bind_history_evaluator_receipt(
            mutated, star, expected_nodes=z, expected_source_identity=bound.source_identity)


def test_missing_history_proof_rejected():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    raw = CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z)
    missing = dict(raw)
    missing.pop('history_identity')
    with pytest.raises(ValueError, match='history identity'):
        bind_history_evaluator_receipt(missing, star, expected_nodes=z)
    undeclared = dict(raw)
    undeclared.pop('full_retarded_state_derivative')
    with pytest.raises(ValueError, match='full_retarded_state_derivative'):
        bind_history_evaluator_receipt(undeclared, star, expected_nodes=z)
    frozen = dict(raw)
    frozen['scope'] = {'constant_frozen_matter_used_as_physical_state': True}
    with pytest.raises(ValueError, match='frozen matter'):
        bind_history_evaluator_receipt(frozen, star, expected_nodes=z)
    removed = dict(raw)
    removed['reference_state_tangent_subtracted'] = True
    with pytest.raises(ValueError, match='removed state derivative'):
        bind_history_evaluator_receipt(removed, star, expected_nodes=z)


def test_cache_copies_immutable_receipt_before_buffer_mutation():
    star = make_family(star_coefficients())
    z = np.array(star.collocation_nodes(8), dtype=float, copy=True)
    holder = {}

    class Buffer:
        def evaluate(self, family, z=None):
            raw = CoupledNonlinearSyntheticEvaluator(star).evaluate(family, z)
            holder['raw'] = raw
            return raw

    cache = GeometryBoundEvaluatorCache(Buffer())
    first = cache.evaluate(star, z)
    original = np.array(first.action_gradient, copy=True)
    holder['raw']['action_gradient'][:] = 123.0
    holder['raw']['source_identity'] = 'mutated-after-cache'
    second = cache.evaluate(star, z)
    np.testing.assert_array_equal(second.action_gradient, original)
    assert second.source_identity == first.source_identity
    assert second.source_identity != 'mutated-after-cache'


def test_initial_negative_radius_never_calls_evaluator():
    coeff = np.zeros((2, 8))
    coeff[0, 0] = 100.0
    family = make_family(coeff)
    assert family.radius_lower_bound() < 0
    evaluator = CoupledNonlinearSyntheticEvaluator(make_family(star_coefficients()))
    with pytest.raises(ValueError, match='positive before evaluator call'):
        solve_local_history(
            evaluator, family, boundary=declared_jet(make_family(star_coefficients())),
            collocation=HistoryCollocation(8, 'tau'))
    assert evaluator.calls == 0


def test_tau_held_out_beta_is_not_numerical_control_converged():
    family = make_family(star_coefficients())
    evaluator = HeldOutBetaSyntheticEvaluator()
    result = solve_local_history(
        evaluator, family, boundary=declared_jet(family),
        collocation=HistoryCollocation(8, 'tau'))
    _assert_open(result)
    assert result['stop_reason'] == STOP_VERIFICATION
    assert result['stop_reason_phrase'] == 'resolution/verification'
    assert result['control_converged'] is True
    assert result['solve_all_nodes_pass'] is False
    assert result['numerical_control_converged'] is False
    assert result['verification']['between_node_certified'] is False
    assert result['collocation']['automatic_resolution'] is False


def test_requested_nodes_must_match_exactly():
    star = make_family(star_coefficients())
    z = star.collocation_nodes(8)
    raw = CoupledNonlinearSyntheticEvaluator(star).evaluate(star, z)
    shifted = dict(raw)
    shifted['z'] = np.array(z, dtype=float, copy=True) + 3e-12
    with pytest.raises(ValueError, match='mismatched incoming nodes'):
        bind_history_evaluator_receipt(shifted, star, expected_nodes=z)


def test_armijo_and_line_search_scale_bounds():
    LocalHistoryNewtonSettings(armijo=0.0, min_line_search_scale=1.0)
    with pytest.raises(ValueError, match='armijo'):
        LocalHistoryNewtonSettings(armijo=1.0)
    with pytest.raises(ValueError, match='armijo'):
        LocalHistoryNewtonSettings(armijo=-1e-8)
    with pytest.raises(ValueError, match='min_line_search_scale'):
        LocalHistoryNewtonSettings(min_line_search_scale=0.0)
    with pytest.raises(ValueError, match='min_line_search_scale'):
        LocalHistoryNewtonSettings(min_line_search_scale=1.1)


def test_analytic_profile_digest_is_supported_and_coefficients_alone_are_not_a_binding():
    from recursive_horizons.nsc_ks_profile_identity import profile_identity
    family=make_family(star_coefficients());z=family.collocation_nodes(8)
    raw=CoupledNonlinearSyntheticEvaluator(family).evaluate(family,z)
    raw.pop('history_identity');raw['profile_identity']=profile_identity(family)
    receipt=bind_history_evaluator_receipt(raw,family,expected_nodes=z)
    assert receipt.family is family
    raw.pop('profile_identity');raw['history_identity']=np.array(family.coefficients,copy=True)
    with pytest.raises(ValueError,match='history identity'):
        bind_history_evaluator_receipt(raw,family,expected_nodes=z)
