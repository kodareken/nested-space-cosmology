"""Independent equation, current, tail and coherent sewing controls."""
import numpy as np
import pytest
from flint import acb, arb, ctx
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_metric_horizon_frame import (
    _profile, horizon_q, metric_horizon_frame, reflection_from_phase,
)
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_unruh_state import ParentDirac

H = 1.9006916054701435
E = .0013248831260437577
M = np.pi/2
L = np.sqrt(5.)


def midpoint(frame):
    return np.array([[complex(v) for v in row] for row in frame])


@pytest.mark.parametrize('interior', [False, True])
@pytest.mark.parametrize('sign', [-1, 1])
def test_full_metric_frame_preserves_exact_dirac_current(interior, sign):
    f = metric_horizon_frame(E, M, sign*L, H)
    with ctx.workprec(192):
        a, tail = f.evaluate('0.00003', interior=interior)
        metric = (1, 1) if interior else (1, -1)
        for j in range(2):
            for k in range(2):
                product = sum((a[i][j].conjugate()*metric[i]*a[i][k]
                               for i in range(2)), acb(0))
                assert product.contains(metric[j] if j == k else 0)
        assert tail < arb('4e-21')


@pytest.mark.parametrize('interior', [False, True])
def test_independent_ode_integration_matches_frame(interior):
    # Direct closed-form metric, standard DOP853, no coefficient recurrence
    # in the comparison equation. This is a numerical sign/chart control;
    # it is not used to supply the Cauchy tail certificate.
    f = metric_horizon_frame(.7, M, L, H)
    q = float(f.q)
    a0 = midpoint(f.evaluate('0.00003', interior=interior)[0])
    def rhs(y, flat):
        d = np.exp(y)
        u = -d if interior else d
        W = 3*(np.pi-q-u)+1.5*np.sin(2*(q+u))-np.sin(q+u)**2
        h = -W/u
        p = np.array([[0, -M/np.sin(q+u)-1j*L],
                      [-M/np.sin(q+u)+1j*L, 0]])
        K = 1j*.7/h*np.diag([1.,-1.]) + (-1j if interior else 1)*np.sqrt(d/h)*p
        return (K@flat.reshape(2,2)).ravel()
    run = solve_ivp(rhs, (np.log(3e-5),np.log(2e-4)),a0.ravel(),
                    method='DOP853',rtol=2e-12,atol=2e-14,max_step=.03)
    expected = midpoint(f.evaluate('0.0002', interior=interior)[0])
    assert run.success
    assert np.max(np.abs(run.y[:,-1].reshape(2,2)-expected)) < 2e-9


def test_exterior_series_generator_matches_original_radial_riccati():
    # Check the ORIGINAL rho-chart Riccati, including its Jacobian. Equating
    # d/dlog(q-qh) and d/dlog(rho-h) without that Jacobian gives a false gap.
    f = metric_horizon_frame(.4,M,L,H,order=30)
    with ctx.workprec(192):
        u = arb('0.00003')
        t = u.sqrt()
        F = midpoint(f.evaluate(u)[0])
        df = np.zeros((2,2),complex)
        for col,sign in enumerate((1,-1)):
            lam = acb(0,sign*f.energy/f.slope)
            phase = (lam*u.log()).exp()
            for row in range(2):
                df[row,col] = complex(phase*sum(
                    ((lam+arb(n)/2)*c[row][col]*t**n
                     for n,c in enumerate(f.coefficients)),acb(0)))
        K = df@np.linalg.inv(F)
        rho = float(-(f.q+u).cos()/(f.q+u).sin())
    background = ParentDirac(H,.23832579963401956)
    A = background.A_from_offset(rho-H)
    v1,v2 = L*np.sqrt(A)/np.sqrt(1+rho*rho),M*np.sqrt(A)
    for z in (.3+.7j,-.8+.1j):
        source_rhs = (1j*(v1+1j*v2)-.8j*z+1j*(v1-1j*v2)*z*z)/A
        jacobian = (1+rho*rho)*float(u)
        frame_rhs = K[1,0]+K[1,1]*z-z*(K[0,0]+K[0,1]*z)
        assert abs(frame_rhs-jacobian*source_rhs) < 2e-10
        wrong = (rho-H)*source_rhs
        assert abs(frame_rhs-wrong) > 1e-6


def test_tail_bounds_the_discarded_coefficients_without_refinement_assumption():
    coarse = metric_horizon_frame(E, M, L, H, order=6)
    fine = metric_horizon_frame(E, M, L, H, order=30)
    with ctx.workprec(192):
        c, tail = coarse.evaluate('0.0001')
        f, fine_tail = fine.evaluate('0.0001')
        assert tail > fine_tail*10**15
        for row in range(2):
            for col in range(2):
                assert c[row][col].contains(f[row][col])


def test_zero_channel_has_exact_constant_frame_and_zero_tail():
    f = metric_horizon_frame(0, 0, 0, H)
    a, tail = f.evaluate('0.001', interior=True)
    assert tail.is_zero()
    assert a == [[acb(1), acb(0)], [acb(0), acb(1)]]


def test_reflection_map_retains_phase_uncertainty_and_unit_current():
    with ctx.workprec(192):
        free = metric_horizon_frame(0,0,0,H)
        theta = arb('.4 +/- 0.0001')
        R, distance, tail = reflection_from_phase(free,arb(H)+arb('.0001'),theta)
        assert R.contains(acb(0,theta).exp())
        assert distance > 0 and tail.is_zero()
        f = metric_horizon_frame(E,M,L,H)
        R, _, tail = reflection_from_phase(f,arb(H)+arb('.0001'),theta)
        assert abs(R).contains(1)
        assert tail < arb('1e-21')
        with pytest.raises(ValueError):
            reflection_from_phase(f,arb(H)-arb('.0001'),theta)


def test_common_phase_sewing_preserves_full_coherent_covariance():
    # This is the actual three-port source, including its off-diagonal pair.
    # Transforming only the frame or only R fails; transforming both is gauge.
    source = source_covariance(.2, .23832579963401956, 1., M)
    c = .37
    phase = np.diag(np.exp([1j*c,-1j*c]))
    R, T = .6*np.exp(.2j), .64
    sewing = lambda r: np.array([[0,1,0],[r,0,np.sqrt(T)]], complex)
    srcphase = np.diag(np.exp([1j*c,1j*c,-1j*c]))
    assert abs(source[0,1]) > .01
    assert np.allclose(phase@sewing(np.exp(2j*c)*R),sewing(R)@srcphase,atol=2e-16)
    original = sewing(R)@source@sewing(R).conj().T
    updated = phase@sewing(np.exp(2j*c)*R)
    assert np.allclose(updated@source@updated.conj().T,original,atol=2e-16)
    wrong = phase@sewing(R)
    assert np.linalg.norm(wrong@source@wrong.conj().T-original) > .01


@pytest.mark.parametrize('kwargs', [{'order':True},{'order':1},{'bits':40},
                                  {'analytic_radius':'0.9'}])
def test_unproved_analytic_domains_rejected(kwargs):
    with pytest.raises((ValueError,ArithmeticError)):
        metric_horizon_frame(E,M,L,H,**kwargs)


def test_root_and_evaluation_domains_and_context_restoration():
    before = ctx.prec, ctx.cap
    f = metric_horizon_frame(E,M,L,H)
    with ctx.workprec(192):
        assert _profile(f.q).contains(0)
        assert f.slope > arb('.47')
    for bad in (0, -.01, .011):
        with pytest.raises(ValueError):
            f.evaluate(bad)
    with pytest.raises(ValueError):
        horizon_q(2.1)
    assert (ctx.prec,ctx.cap) == before
