"""Independent curvature, orientation and detector-clock controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_nested_spacetime import (
    central_inward_flash, regular_kottler_chart, static_slice_frame,
    central_radar_clock, invert_central_radar_clock, radar_junction_branches,
    radar_cap_input_bounds,
    radar_clock_metric,
)


def test_regular_metric_einstein_equations_from_warped_curvature():
    # Generic exact symbols independently differentiate the metric and r;
    # no target host parameters or numerical wall application is evaluated.
    u,w,C,L=s.symbols('u w C L', nonzero=True, real=True)
    k=(1-3*L*C*C)/(2*C);r=C+k*u*w
    F=w*w/r*(1/C+L*(r-C))
    Fu,Fw,Fww=s.diff(F,u),s.diff(F,w),s.diff(F,w,2)
    hu=Fw*k*w/2-(Fu+F*Fw)*k*u/2
    huw=k-Fw*k*u/2
    radial_norm=2*k*k*u*w-F*k*k*u*u
    f=1-(C-L*C**3)/r-L*r*r
    for residual in (radial_norm-f, Fww*F/2-2*hu/r-3*L*F,
                     Fww/2-2*huw/r-3*L, 1-f-2*r*huw-3*L*r*r):
        assert s.factor(residual)==0


@pytest.mark.parametrize('C,H', [(2.,0.),(1.,.25),(3.,.25)])
def test_chart_jets_and_mass_invariant(C,H):
    # Includes a negative-surface-gravity cosmological horizon.
    mass=C*(1-H*H*C*C)/2
    u,w=s.symbols('u w');k=(1-3*H*H*C*C)/(2*C)
    r=C+k*u*w;F=w*w/r*(1/C+H*H*(r-C))
    row=regular_kottler_chart(.2,-.3,mass=mass,H=H,horizon=C)
    for name,expression in [('F',F),('F_u',s.diff(F,u)),
                            ('F_w',s.diff(F,w)),('F_ww',s.diff(F,w,2))]:
        expected=float(expression.subs({u:.2,w:-.3}))
        assert float(row[name])==pytest.approx(expected,rel=2e-13,abs=1e-14)
    norm=2*row['R_u']*row['R_w']-row['F']*row['R_w']**2
    assert norm==pytest.approx(1-2*mass/row['R']-(H*row['R'])**2,abs=3e-15)
    axis=regular_kottler_chart(np.array([0.,.1,0.]),np.array([.2,0.,0.]),
                               mass=mass,H=H,horizon=C)
    assert np.all(axis['R']==C)
    assert np.all(np.isfinite(axis['F']))
    assert axis['radial_metric_determinant']==-1


@pytest.mark.parametrize('orientation',[-1,1])
def test_oriented_frame_preserves_one_killing_time(orientation):
    q=static_slice_frame(3.,.2,mass=1.,H=0.,orientation=orientation)
    f=q['f'];metric=np.diag([-f,1/f])
    normal=np.array([q['N_T'],q['N_R']]);tangent=np.array([q['s_T'],q['s_R']])
    assert normal@metric@normal==pytest.approx(-1,abs=1e-15)
    assert tangent@metric@tangent==pytest.approx(1,abs=1e-15)
    assert normal@metric@tangent==pytest.approx(0,abs=1e-15)
    assert np.linalg.det(np.column_stack([normal,tangent]))==pytest.approx(1)
    assert math.copysign(1,q['N_T'])==orientation
    assert q['outgoing_l']==pytest.approx(q['T_l']-q['Rprime']/f)


def test_flash_flat_limit_and_independent_doppler_factor():
    v=.3;g=1/math.sqrt(1-v*v)
    flat=central_inward_flash(radius=2.,H=0.,beta_minus=g,Rdot=g*v)
    assert flat['central_proper_time']==2.
    assert flat['one_plus_redshift']==pytest.approx(math.sqrt((1+v)/(1-v)))
    cap=central_inward_flash(radius=.5,H=.5,beta_minus=1.,Rdot=.25)
    assert cap['central_proper_time']==pytest.approx(2*math.log(4/3))
    assert cap['one_plus_redshift']==pytest.approx(4/3)


def test_flash_and_chart_length_scaling():
    flash=central_inward_flash(radius=.5,H=.5,beta_minus=1.,Rdot=.25)
    scaled=central_inward_flash(radius=3.5,H=.5/7,beta_minus=1.,Rdot=.25)
    assert scaled['central_proper_time']==pytest.approx(7*flash['central_proper_time'])
    assert scaled['one_plus_redshift']==flash['one_plus_redshift']
    base=regular_kottler_chart(.2,-.3,mass=1.,H=0.,horizon=2.)
    big=regular_kottler_chart(1.4,-2.1,mass=7.,H=0.,horizon=14.)
    assert big['R']==pytest.approx(7*base['R'])
    assert big['F']==pytest.approx(base['F'])
    assert big['F_ww']==pytest.approx(base['F_ww']/49)


def test_invalid_chart_and_flash_domains_stay_explicit():
    with pytest.raises(ValueError,match='root'):
        regular_kottler_chart(0.,0.,mass=1.,H=0.,horizon=3.)
    with pytest.raises(ValueError,match='simple'):
        regular_kottler_chart(0.,0.,mass=1/3,H=1/math.sqrt(3),horizon=1.)
    with pytest.raises(ValueError,match='R>0'):
        regular_kottler_chart(10.,-10.,mass=1.,H=0.,horizon=2.)
    with pytest.raises(ValueError,match='f>0'):
        static_slice_frame(2.,0.,mass=1.,H=0.,orientation=1)
    with pytest.raises(ValueError,match=r'H\*R<1'):
        central_inward_flash(radius=1.,H=1.,beta_minus=1.,Rdot=1.)
    with pytest.raises(ValueError,match='factors'):
        central_inward_flash(radius=.5,H=.5,beta_minus=2.,Rdot=.25)


def _clock(R=.6, H=.4, beta=1.2, acceleration=.3, direction=1):
    speed = direction*math.sqrt(beta*beta-1+(H*R)**2)
    return central_radar_clock(radius=R,H=H,beta_minus=beta,
                               Rdot=speed,Rddot=acceleration)


def test_radar_recovers_metric_from_two_independent_null_clock_derivatives():
    # Direct normalization in radar coordinates works even when F(D) is not
    # de Sitter and the angular radius function is left arbitrary.
    F, D_tau = math.exp(-.4), .7
    T_tau = math.sqrt(1/F+D_tau*D_tau)
    receive, send = T_tau+D_tau, T_tau-D_tau
    assert (receive/send)/receive**2 == pytest.approx(F,rel=3e-16)
    for direction in (-1,1):
        clock = _clock(direction=direction)
        inferred = invert_central_radar_clock(**{key:clock[key] for key in
            ('radar_distance','one_plus_redshift','echo_stretch','log_redshift_rate')})
        assert inferred['radius'] == pytest.approx(.6,rel=2e-15)
        assert inferred['H'] == pytest.approx(.4,rel=2e-15)
        assert inferred['Rdot']*direction > 0
        assert inferred['beta_minus'] == pytest.approx(1.2)


def test_radar_metric_readout_survives_an_evolving_nonvacuum_geometry():
    # Exact worldline in ds2=exp(2*a*T)*(-dT2+dD2), D=v*T.
    # Its factor C>1 is outside the admitted quiet de Sitter inverse,
    # yet proper-clock readout remains valid without Einstein equations.
    a,v=.3,.2
    for tau in (.1,.4,.7):
        T=math.log1p(a*tau/math.sqrt(1-v*v))/a
        rate=math.exp(-a*T)/math.sqrt(1-v*v)
        Z=(1+v)*rate
        K=(1+v)/(1-v)
        measured=radar_clock_metric(one_plus_redshift=Z,echo_stretch=K)
        assert measured['radial_metric_factor']==pytest.approx(math.exp(2*a*T),rel=3e-16)
        assert measured['centre_clock_rate']==pytest.approx(rate,rel=3e-16)
        assert measured['radar_distance_rate']==pytest.approx(v*rate,rel=8e-16)
        assert measured['radar_velocity']==pytest.approx(v,rel=8e-16)
        with pytest.raises(ValueError,match='incompatible'):
            invert_central_radar_clock(radar_distance=v*T,one_plus_redshift=Z,
                                       echo_stretch=K,log_redshift_rate=0.)


def test_radar_drift_from_differentiated_null_clocks_symbolically():
    R,H,a,c=s.symbols('R H a c', positive=True)
    beta=a/R**2+c*R
    v=s.sqrt(beta**2-1+H**2*R**2)
    Z=(beta+v)/(1-H**2*R**2)
    drift=v*s.diff(Z,R)/Z**2
    assert s.simplify(drift-s.diff(beta,R)/Z-H**2*R)==0


def test_radar_flat_limit_and_no_positive_tension_hypothesis():
    clock=_clock(H=0.)
    row=invert_central_radar_clock(**{key:clock[key] for key in
        ('radar_distance','one_plus_redshift','echo_stretch','log_redshift_rate')})
    # H itself has a square-root conditioning singularity at zero; H^2
    # retains the clock roundoff. Do not erase a small real curvature.
    assert row['H']**2*.6**2 < 8*np.finfo(float).eps
    assert row['radius']==pytest.approx(.6)
    assert row['Lambda_cap']*.6**2 < 24*np.finfo(float).eps
    flat=invert_central_radar_clock(radar_distance=.6,one_plus_redshift=2.,
                                    echo_stretch=4.,log_redshift_rate=.1)
    assert flat['H']==0 and flat['radius']==.6
    assert radar_junction_branches(dict(H=1.,mass_over_tension=.1,
        junction_radius_coefficient=.1),host_H=2.)==[]


def test_radar_junction_inverse_keeps_two_admissible_tensions():
    # H_host^2-H_cap^2=.75, c=1: k=.5 and 1.5 have the SAME
    # a,c and therefore every cap-side clock history. Do not invent uniqueness.
    inferred=dict(H=.5,mass_over_tension=.2,junction_radius_coefficient=1.)
    branches=radar_junction_branches(inferred,host_H=1.)
    assert [row['tension'] for row in branches]==[.5,1.5]
    assert [row['mass'] for row in branches]==pytest.approx([.1,.3])


def test_radar_inverse_recovers_junction_from_existing_signed_bulk_geometry():
    from recursive_horizons.nsc_atlas_nested_slice import tension_wall
    state=tension_wall(.6,0.,.4,.4,0.,2.)
    clock=central_radar_clock(radius=.6,H=.4,beta_minus=state['beta_minus'],
                             Rdot=state['Rdot'],Rddot=state['Rddot'])
    row=invert_central_radar_clock(**{key:clock[key] for key in
        ('radar_distance','one_plus_redshift','echo_stretch','log_redshift_rate')})
    branches=radar_junction_branches(row,host_H=0.)
    assert len(branches)==1
    assert branches[0]['tension']==pytest.approx(2.,rel=3e-15)
    assert branches[0]['mass']==pytest.approx(.4,rel=3e-15)


def test_radar_input_box_contains_independent_clock_perturbations():
    clock=_clock()
    keys=('radar_distance','one_plus_redshift','echo_stretch')
    bounds={key:(clock[key]*(1-1e-4),clock[key]*(1+1e-4)) for key in keys}
    enclosure=radar_cap_input_bounds(**bounds)
    rng=np.random.default_rng(407)
    for _ in range(150):
        query={key:rng.uniform(*bounds[key]) for key in keys}
        query['log_redshift_rate']=clock['log_redshift_rate']
        row=invert_central_radar_clock(**query)
        for key in ('radius','Lambda_cap'):
            assert enclosure[key][0] <= row[key] <= enclosure[key][1]


def test_radar_inverse_length_and_clock_scaling():
    clock=_clock()
    query={key:clock[key] for key in
        ('radar_distance','one_plus_redshift','echo_stretch','log_redshift_rate')}
    small=invert_central_radar_clock(**query)
    query['radar_distance']*=7
    query['log_redshift_rate']/=7
    large=invert_central_radar_clock(**query)
    for key,factor in (('radius',7),('H',1/7),('Lambda_cap',1/49),
                       ('mass_over_tension',49),('junction_radius_coefficient',1/7)):
        assert large[key]==pytest.approx(small[key]*factor,rel=3e-15)


def test_radar_rejects_outside_horizon_and_incompatible_clocks():
    with pytest.raises(ValueError,match='static cap'):
        central_radar_clock(radius=1.,H=1.,beta_minus=1.,Rdot=1.,Rddot=0.)
    with pytest.raises(ValueError,match='disagree'):
        central_radar_clock(radius=.5,H=1.,beta_minus=1.,Rdot=0.,Rddot=0.)
    with pytest.raises(ValueError,match='incompatible'):
        invert_central_radar_clock(radar_distance=1.,one_plus_redshift=1.,
                                  echo_stretch=2.,log_redshift_rate=0.)
    with pytest.raises(ValueError,match='entire input box'):
        radar_cap_input_bounds(radar_distance=(1.,2.),
            one_plus_redshift=(1.,1.1),echo_stretch=(2.,3.))


@pytest.mark.parametrize('R0,mass,host',[(.5,.4,0.),(.8,.01,1.)])
def test_finite_clock_prediction_against_proper_time_and_changed_material(R0,mass,host):
    from check_nsc_atlas_radar_inverse import independent_wall_history, predict_reception
    from recursive_horizons.nsc_atlas_nested_slice import tension_wall
    state=tension_wall(R0,0.,1.,mass,host,2.)
    clock=central_radar_clock(radius=R0,H=1.,beta_minus=state['beta_minus'],
                             Rdot=state['Rdot'],Rddot=state['Rddot'])
    row=invert_central_radar_clock(**{key:clock[key] for key in
        ('radar_distance','one_plus_redshift','echo_stretch','log_redshift_rate')})
    predicted=predict_reception(row,.04,64)
    direct=independent_wall_history(R0,mass,host,.04)
    alternative=independent_wall_history(R0,mass,host,.04,change_law=.5)
    assert direct['radius']==pytest.approx(predicted['radius'],abs=3e-13)
    assert direct['one_plus_redshift']==pytest.approx(predicted['one_plus_redshift'],abs=3e-12)
    assert abs(direct['temporal_drift_check'])<1e-8
    assert alternative['positive_wall_surface_density']
    assert abs(alternative['one_plus_redshift']-predicted['one_plus_redshift'])>1e-7


def test_nested_clock_prediction_against_eighty_digit_quadrature():
    import mpmath as mp
    from check_nsc_atlas_radar_inverse import predict_reception
    with mp.workdps(80):
        r0=mp.mpf(4)/5
        beta=lambda r:mp.mpf(1)/200/r**2+r
        velocity=lambda r:mp.sqrt(beta(r)**2-1+r*r)
        Z=lambda r:(beta(r)+velocity(r))/(1-r*r)
        elapsed=lambda r:mp.quad(lambda x:Z(x)/velocity(x),[r0,r])
        endpoint=mp.findroot(lambda r:elapsed(r)-mp.mpf(1)/25,
                             (r0,r0+mp.mpf(1)/100))
        row=dict(radius=.8,H=1.,mass_over_tension=.005,junction_radius_coefficient=1.)
        result=predict_reception(row,.04,64)
        assert abs(mp.mpf(result['radius'])-endpoint)<mp.mpf('2e-15')
        assert abs(mp.mpf(result['one_plus_redshift'])-Z(endpoint))<mp.mpf('3e-14')
