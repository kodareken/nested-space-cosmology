"""Independent cap constraints and actual conservative background transport."""
import time

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import (
    background_cap_metric, background_cap_profile, background_cap_rhs,
)


@pytest.fixture(scope="module")
def budget():
    b = Budget(maximum_evaluations=100000, cpu_seconds=30)
    yield b
    b.check()
    print("background-cap evaluations:", b.evaluations,
          "process CPU:", time.process_time()-b.started)


def profile(scale=0.):
    x = np.array([0., .3, .6, 1.3])
    return CharacteristicProfile(x, scale*x**3, 3*scale*x*x, 6*scale*x)


def state(n=4, R=.5, f=0., Pi=0.):
    y = np.linspace(0, 1, n+1)
    r = R*y
    return y, {"Ip": Pi*np.diff(r*r)/2, "Iq": f*np.diff(r)}


@pytest.mark.parametrize("alpha", [0., .08])
def test_quiet_vacuum_exact(alpha, budget):
    y, s = state()
    m = background_cap_metric(y, .5, s, T=0., profile=profile(), alpha=alpha, budget=budget)
    np.testing.assert_array_equal(m["m_faces"], m["r_faces"]**3/2)
    np.testing.assert_array_equal(m["delta_faces"], np.zeros(5))
    for v in (-.1, .1):
        rhs = background_cap_rhs(m, T_x=.8, radius_x=v, incoming_w_plus=0., budget=budget)
        for rate in rhs["rates"].values():
            np.testing.assert_array_equal(rate, np.zeros(4))
        assert rhs["wall_chi_x"] == rhs["wall_c_x"] == 0


def test_constant_field_closed_mass_clock_and_reaction(budget):
    f, alpha, R = .3, .2, .5
    y, s = state(R=R, f=f)
    m = background_cap_metric(y, R, s, T=0., profile=profile(), alpha=alpha, budget=budget)
    z = np.sqrt(alpha/6)*f
    P, U = 3*np.cosh(z)**4, 2*f*np.sinh(z)/z*np.cosh(z)**3
    r, rn = m["r_faces"], m["r_nodes"]
    np.testing.assert_allclose(m["m_faces"], P*r**3/6, atol=2e-16)
    np.testing.assert_allclose(m["delta_faces"], -np.log(np.cosh(z)), atol=2e-16)
    np.testing.assert_allclose(m["A_faces"], 1-P*r*r/3, atol=2e-16)
    rhs = background_cap_rhs(m, T_x=1., radius_x=0., incoming_w_plus=0., budget=budget)
    expected = rn*(U-2*P*f/3)/np.cosh(z)
    np.testing.assert_allclose(rhs["reaction_nodes"], expected, atol=2e-16)
    np.testing.assert_allclose(rhs["reaction_integral"],
                               (U-2*P*f/3)*np.diff(r*r)/2/np.cosh(z), atol=2e-16)


def test_nonzero_pi_against_independent_integrating_factor(budget):
    f, pi, alpha, R = .2, .6, .12, .55
    y, s = state(R=R, f=f, Pi=pi)
    errors = []
    z = np.sqrt(alpha/6)*f
    Pdiff = (3*np.cosh(z)**4-3)/alpha
    for order in (4, 8):
        m = background_cap_metric(y, R, s, T=.1, profile=profile(), alpha=alpha,
                                  quadrature_order=order, budget=budget)
        expected = []
        for r in m["r_faces"]:
            def integrand(v):
                budget.check(1)
                return np.exp(alpha*pi*pi*v*v/4)*v*v*(Pdiff/2+(1-v*v)*pi*pi/4)
            expected.append(np.exp(-alpha*pi*pi*r*r/4)*quad(integrand,0,r,epsabs=1e-13)[0])
        errors.append(np.max(abs(m["n_faces"]-expected)))
        np.testing.assert_allclose(m["delta_faces"],
            -np.log(np.cosh(z))+alpha*pi*pi*m["r_faces"]**2/4, atol=3e-16)
    assert errors[1] < 2e-14
    assert errors[0] < 2e-14  # Both quadratures already reach roundoff here.


def test_total_field_energy_cross_terms_against_radial_ode(budget):
    # Explicit incoming cubic profile and constant departure: the ODE oracle
    # independently evaluates the analytic total field, not the reconstruction.
    R, T, alpha, amp, fdep, pidep = .5, .12, .08, .06, .15, .2
    y, s = state(R=R, f=fdep, Pi=pidep)
    m = background_cap_metric(y, R, s, T=T, profile=profile(amp), alpha=alpha,
                              quadrature_order=8, budget=budget)
    def ode(r, values):
        budget.check(1)
        if r == 0:
            return np.zeros(2)
        xi, A0 = np.arctanh(r), 1-r*r
        zp, zm = max(T+xi,0), max(T-xi,0)
        chi = amp*(zp**3-zm**3)
        p = 3*amp*(zp*zp-zm*zm)/A0
        q = 3*amp*(zp*zp+zm*zm)/A0
        f, Pi, Phi = chi/r+fdep, p/r+pidep, q/r-chi/r**2
        zz = np.sqrt(alpha/6)*f
        Pdiff = 3*np.sinh(zz)**2*(np.cosh(zz)**2+1)/alpha
        W = Pi*Pi+Phi*Phi
        a = alpha*r*W/2
        return [a, r*r*(Pdiff/2+(1-r*r)*W/4)-a*values[1]]
    # Independently split the quiet/active reflected boundary.
    last = np.zeros(2)
    for a,b in zip([0.,np.tanh(T)],[np.tanh(T),R]):
        sol = solve_ivp(ode,(a,b),last,rtol=2e-11,atol=2e-13)
        assert sol.success
        last = sol.y[:,-1]
    assert m["n_faces"][-1] == pytest.approx(last[1],abs=3e-11)
    f0 = 6*amp*T*T+fdep
    assert m["delta_faces"][-1] == pytest.approx(last[0]-np.log(np.cosh(np.sqrt(alpha/6)*f0)),abs=3e-11)
    assert m["center"]["Pi"] == pytest.approx(12*amp*T+pidep)


def test_exact_background_advances_with_moving_cells_without_spreading(budget):
    # Eight real Euler steps on zero departures: only analytic physical
    # cell integrals change, and an interior pulse is not advected numerically.
    prof = profile(.09)
    y,s = state(n=3,R=.45)
    dt = .002
    for k in range(8):
        T, R = .1+k*dt, .45+.1*k*dt
        m = background_cap_metric(y,R,s,T=T,profile=prof,alpha=0.,budget=budget)
        bg = m["background_faces"]
        wp = (bg["p"][-1]+bg["q"][-1]-m["wall"]["f"])/R
        rhs = background_cap_rhs(m,T_x=1.,radius_x=.1,incoming_w_plus=wp,budget=budget)
        assert max(np.max(abs(v)) for v in rhs["rates"].values()) < 2e-15
        s = {key:s[key]+dt*rhs["rates"][key] for key in s}
        r = R*y
        exact = .09*(np.maximum(T+np.arctanh(r),0)**3-np.maximum(T-np.arctanh(r),0)**3)
        np.testing.assert_allclose(m["chi_faces"],exact,atol=3e-16)
        assert abs(rhs["constructed_wall_flux_defect"]) < 1e-16


def test_incoming_discrepancy_is_free_and_wall_chain_uses_background(budget):
    prof = profile(.05)
    R,Tx,Rx = .5,.8,.1
    y,s = state(R=R)
    m = background_cap_metric(y,R,s,T=.1,profile=prof,alpha=.04,budget=budget)
    bg=m["background_faces"]
    wp=(bg["p"][-1]+bg["q"][-1]-m["wall"]["f"])/R
    a=background_cap_rhs(m,T_x=Tx,radius_x=Rx,incoming_w_plus=wp,budget=budget)
    b=background_cap_rhs(m,T_x=Tx,radius_x=Rx,incoming_w_plus=wp+.03,budget=budget)
    assert b["incoming_departure"]-a["incoming_departure"] == pytest.approx(R*.03)
    expected = (Tx*m["wall"]["B"]+Rx)*R*.03/2
    assert b["rates"]["Iq"][-1]-a["rates"]["Iq"][-1] == pytest.approx(expected)
    np.testing.assert_array_equal(b["rates"]["Iq"][:-1],a["rates"]["Iq"][:-1])
    assert abs(b["constructed_wall_flux_defect"]) < 1e-15
    assert b["wall_chi_x"] == pytest.approx(b["baseline_wall_chi_x"]+sum(b["rates"]["Iq"]))
    # Independent canonical/conformal directional chain rule.
    h=1e-6
    chi=m["wall"]["chi"]
    def conv(t):
        rr=R+t*Rx
        return rr*np.sinh(np.sqrt(.04/6)*(chi+t*b["wall_chi_x"])/rr)/np.sqrt(.04/6)
    assert b["wall_c_x"] == pytest.approx((conv(h)-conv(-h))/(2*h),abs=2e-11)


def test_nonzero_departure_transport_refines_against_distinct_exact_wave(budget):
    T,R,Tx,Rx,amp,eps=.1,.5,1.,.1,.05,.2
    errors=[]
    for n in (8,16):
        y=np.linspace(0,1,n+1)
        r=R*y
        xp,xm=np.maximum(T+np.arctanh(r),0),np.maximum(T-np.arctanh(r),0)
        s={"Ip":eps*amp*np.diff(xp**3+xm**3),"Iq":eps*amp*np.diff(xp**3-xm**3)}
        m=background_cap_metric(y,R,s,T=T,profile=profile(amp),alpha=0.,budget=budget)
        total_p=3*(1+eps)*amp*(xp*xp-xm*xm)/(1-r*r)
        total_q=3*(1+eps)*amp*(xp*xp+xm*xm)/(1-r*r)
        wp=(total_p[-1]+total_q[-1]-m["wall"]["f"])/R
        rhs=background_cap_rhs(m,T_x=Tx,radius_x=Rx,incoming_w_plus=wp,budget=budget)
        expected_p=np.diff(3*eps*amp*((xp*xp+xm*xm)*Tx
                                     +(xp*xp-xm*xm)*Rx*y/(1-r*r)))
        expected_q=np.diff(3*eps*amp*((xp*xp-xm*xm)*Tx
                                     +(xp*xp+xm*xm)*Rx*y/(1-r*r)))
        errors.append(sum(np.sum(abs(rhs["rates"][key]-expected))
                          for key,expected in (("Ip",expected_p),("Iq",expected_q))))
        assert np.max(abs(rhs["rates"]["Iq"])) > 1e-4
    assert errors[1] < errors[0]/1.5


def test_profile_and_panel_bounds(budget):
    y,s=state()
    m=background_cap_metric(y,.5,s,T=.1,profile=profile(.03),alpha=.02,budget=budget)
    assert all(v in m["panel_cuts"] for v in m["r_faces"])
    assert np.tanh(.1) in m["panel_cuts"]
    row=background_cap_profile(m,np.array([0.,.2,.5]),budget=budget)
    assert row["chi"][0] == row["p"][0] == row["Phi"][0] == 0
    assert row["chi"][-1] == pytest.approx(m["wall"]["chi"])
    with pytest.raises(ValueError,match="inside"):
        background_cap_profile(m,np.array([.6]),budget=budget)
    with pytest.raises(ValueError,match="shared"):
        background_cap_rhs(m,T_x=1.,radius_x=0.,incoming_w_plus=0.,budget=Budget())
    with pytest.raises(ValueError,match="timelike"):
        background_cap_rhs(m,T_x=1.,radius_x=2.,incoming_w_plus=0.,budget=budget)


def test_narrow_pulse_energy_quadrature_against_separate_radial_integration(budget):
    # The pulse is narrower than a cap cell. Its exact knot splits, not grid
    # smoothing or a cell mean, must resolve the gradient-squared source.
    T,R,alpha=.1,.5,1e-4
    knots=np.array([0.,.25,.27,.29,1.])
    prof=CharacteristicProfile(knots,[0.,0.,.2,0.,0.],np.zeros(5),np.zeros(5))
    y,s=state(n=3,R=R)
    def Fj(v):
        if .25 <= v <= .27:
            z=(v-.25)/.02
            return .2*(10*z**3-15*z**4+6*z**5),10*(30*z*z-60*z**3+30*z**4)
        if .27 < v <= .29:
            z=(v-.27)/.02
            return .2*(1-10*z**3+15*z**4-6*z**5),-10*(30*z*z-60*z**3+30*z**4)
        return 0.,0.
    def rhs(r, v):
        budget.check(1)
        if r==0:
            return np.zeros(2)
        F,Fp=Fj(T+np.arctanh(r))
        f=F/r
        Pi=Fp/((1-r*r)*r)
        Phi=Pi-F/r**2
        z=np.sqrt(alpha/6)*f
        pdiff=3*np.sinh(z)**2*(np.cosh(z)**2+1)/alpha
        a=alpha*r*(Pi*Pi+Phi*Phi)/2
        b=r*r*(pdiff/2+(1-r*r)*(Pi*Pi+Phi*Phi)/4)
        return [a,b-a*v[1]]
    last=np.zeros(2)
    cuts=[0.,*np.tanh(np.array([.25,.27,.29])-T),R]
    for a,b in zip(cuts[:-1],cuts[1:]):
        sol=solve_ivp(rhs,(a,b),last,rtol=1e-10,atol=1e-12)
        assert sol.success
        last=sol.y[:,-1]
    errors=[]
    gamma_disagreements=[]
    for order in (4,8,16):
        m=background_cap_metric(y,R,s,T=T,profile=prof,alpha=alpha,
                                quadrature_order=order,budget=budget)
        errors.append(abs(m["n_faces"][-1]-last[1]))
        gamma_disagreements.append(np.max(abs(m["gamma_interpolant_defect_nodes"])))
        assert np.all(m["n_nodes"]>=-1e-12)
        assert np.all(m["delta_nodes"]>=-1e-12)
    assert errors[1] < errors[0]/100
    assert errors[2] < 1e-7
    assert gamma_disagreements[1] < gamma_disagreements[0]/10
    assert gamma_disagreements[2] < gamma_disagreements[1]/10
    print("narrow-pulse mass errors:",errors,"gamma consistency:",gamma_disagreements)
