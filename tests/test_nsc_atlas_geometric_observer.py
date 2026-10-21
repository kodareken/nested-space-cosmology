import math
import numpy as np
import pytest
from scipy.integrate import quad

from recursive_horizons.nsc_atlas_geometric_observer import (
    WallGeometry,ScalarWallObserver,regular_wall_geometry,regular_bulk_Ruu,
    equivalent_one_wall_example,infer_mixing_and_charge)
from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux


@pytest.fixture(scope='module')
def reference():
    o=ScalarWallObserver()
    sol,_=o.prepared.solve(tolerance=1e-9)
    return o,sol


def test_one_sided_geometry_recovers_regular_field_and_tidal_projection(reference):
    o,sol=reference
    truth=o.prepared.scalar.minimum(o.prepared.child_flux)+sol.y[2,-1]
    for eta in (0.,.5,1.):
        row=o.reconstruct(regular_wall_geometry(o.prepared,sol,eta),measured_Ruu=regular_bulk_Ruu(sol))
        assert row['scalar_value']==pytest.approx(truth,abs=2e-9)
        assert row['normal_magnitude']==pytest.approx(abs(sol.y[3,-1]),abs=2e-8)
        assert row['signed_normal'] is None  # geometry alone does not choose the sign at Pi=0
        assert abs(row['tidal_residual'])<2e-8
        assert row['inferred_vacuum']==pytest.approx(.125,abs=2e-8)


def test_mass_can_grow_with_no_influx_due_to_independent_pressure_work(reference):
    o,sol=reference
    a=o.reconstruct(regular_wall_geometry(o.prepared,sol,0.))
    b=o.reconstruct(regular_wall_geometry(o.prepared,sol,1.))
    R0=sol.y[0,-1]
    U=o.prepared.child_vacuum+o.prepared.scalar.omega**2*sol.y[2,-1]**2/2
    p=sol.y[3,-1]**2/2-U
    work=quad(lambda t:-o.prepared.scalar.kappa*(R0*math.cosh(t/R0))**2*math.sinh(t/R0)
              *(p+o.prepared.child_vacuum)/2,0,R0,epsabs=1e-13)[0]
    assert a['inward_flux']==b['inward_flux']==0
    growth=b['vacuum_subtracted_mass']-a['vacuum_subtracted_mass']
    assert growth>0
    assert growth==pytest.approx(work,abs=1e-9)


def synthetic_jet(o,Pi=-.02,Phi=-.1):
    # Local algebraic recovery fixture, not a claim of a bulk evolution.
    m=o.prepared
    phi=m.scalar.minimum(m.child_flux)-.01
    U=m.scalar.potential(phi,m.child_flux)
    p=(Pi*Pi+Phi*Phi)/2-U
    P=m.scalar.wall_pressure(phi,m.parent_flux,m.charge)
    R,v,B=2.,.1,.3
    A=P/m.tension+1.5*m.k-2*B
    beta=B*R
    mass=R*(1+v*v-beta*beta)/2
    d=-mass/R**2-m.scalar.kappa*R*p/2+A*beta
    Bd=v/R*(A-B)-m.scalar.kappa*Pi*Phi/2
    Ad=o.force_scale*Pi/m.tension-2*Bd
    return WallGeometry(R,v,d,A,B,Ad,Bd),phi,m.scalar.kappa*(Pi*Pi-U)


def test_nonzero_synthetic_flux_recovers_the_signed_normal_gradient():
    o=ScalarWallObserver()
    for Pi,Phi in ((-.02,-.1),(.03,-.15),(-.04,.2)):
        geo,phi,Ruu=synthetic_jet(o,Pi,Phi)
        row=o.reconstruct(geo,measured_Ruu=Ruu)
        assert row['scalar_value']==pytest.approx(phi)
        assert row['scalar_time_derivative']==pytest.approx(Pi)
        assert row['signed_normal']==pytest.approx(Phi)
        assert row['inward_flux']==pytest.approx(Pi*Phi)
        assert abs(row['scalar_stress_residual'])<1e-12
        assert abs(row['tidal_residual'])<1e-12


def test_changed_flux_geometry_is_detected_even_with_unchanged_pressure_history():
    o=ScalarWallObserver()
    geo,_,_=synthetic_jet(o)
    bad=WallGeometry(geo.radius,geo.radius_dot,geo.radius_ddot,geo.A,geo.B,
                     geo.A_dot-.02,geo.B_dot+.01)
    row=o.reconstruct(bad)
    assert abs(row['scalar_stress_residual'])>1e-5


def test_independent_tidal_check_rejects_a_gauss_only_false_positive():
    o=ScalarWallObserver()
    minkowski_tube=WallGeometry(1.,0.,0.,0.,1.,0.,0.)
    row=o.reconstruct(minkowski_tube,measured_Ruu=0.)
    assert row['scalar_stress_residual']==0
    assert row['normal_square']>0
    assert row['tidal_residual']>.1


def sample(o,P,Pi):
    chi=(P-o.ground_pressure)/o.force_scale
    U=o.prepared.child_vacuum+o.prepared.scalar.omega**2*chi*chi/2
    return dict(pressure=P,pressure_dot=o.force_scale*Pi,Ruu=o.prepared.scalar.kappa*(Pi*Pi-U))


def test_two_calibration_samples_predict_a_withheld_tidal_projection():
    o=ScalarWallObserver()
    samples=[sample(o,p,pi) for p,pi in ((.2,.03),(.22,.04),(.24,.02))]
    fit=o.calibrate_tidal_pair(*samples[:2])
    assert fit['inferred_child_flux']==pytest.approx(1.)
    assert fit['vacuum']==pytest.approx(.125)
    assert abs(fit['fixed_action_vacuum_residual'])<1e-12
    predicted=o.predict_tidal(.24,samples[2]['pressure_dot'],fit)
    assert predicted==pytest.approx(samples[2]['Ruu'],abs=1e-12)
    assert o.divided_difference(samples)==pytest.approx(-8/3,abs=1e-11)
    samples[2]['Ruu']+=.001
    assert abs(o.divided_difference(samples)+8/3)>1.


def test_vacuum_offset_cannot_be_silently_retuned_after_calibration():
    o=ScalarWallObserver()
    pair=[sample(o,.2,.03),sample(o,.22,.04)]
    for s in pair:s['Ruu']+=.01
    fit=o.calibrate_tidal_pair(*pair)
    assert fit['ground_pressure']==pytest.approx(.25)
    assert fit['fixed_action_vacuum_residual']==pytest.approx(-.01)


def test_degenerate_pressure_samples_and_zero_scalar_coupling_are_explicit():
    o=ScalarWallObserver()
    s=sample(o,.2,.03)
    with pytest.raises(ValueError):o.calibrate_tidal_pair(s,s)
    pure=ScalarWallObserver(ScalarPreparation(scalar=ScalarFlux(epsilon=1.)))
    with pytest.raises(ValueError):pure.reconstruct(WallGeometry(1.,0.,0.,0.,1.,0.,0.))


def test_negative_normal_square_is_not_hidden_by_a_square_root():
    row=ScalarWallObserver().reconstruct(WallGeometry(1.,0.,10.,0.,1.,0.,0.))
    assert row['normal_square']<0
    assert row['normal_magnitude'] is None


def test_one_wall_couplings_are_degenerate_but_next_sector_is_not():
    first,second,shift=equivalent_one_wall_example()
    for phi in (-.08,-.05,.02):
        for qa,qb in ((first.child_flux,second.child_flux),(first.parent_flux,second.parent_flux)):
            assert second.scalar.potential(phi+shift,qb)==pytest.approx(first.scalar.potential(phi,qa),abs=1e-12)
        assert second.scalar.wall_pressure(phi+shift,second.parent_flux,second.charge)==pytest.approx(
            first.scalar.wall_pressure(phi,first.parent_flux,first.charge),abs=1e-12)
    # These are precisely the inputs of the regular-cap displacement/gravity equations.
    for name in ('child_vacuum','parent_vacuum','shift','k'):
        assert getattr(first,name)==pytest.approx(getattr(second,name),abs=1e-12)
    assert first.scalar.omega==second.scalar.omega
    assert first.scalar.vacuum(first.child_flux-first.charge)==pytest.approx(1/32)
    assert second.scalar.vacuum(second.child_flux-second.charge)==pytest.approx(49/32)


def test_additional_ground_sector_resolves_mixing_and_charge_under_shared_law():
    first,second,_=equivalent_one_wall_example()
    for model in (first,second):
        o=ScalarWallObserver(model)
        values=[model.scalar.kappa*model.scalar.vacuum(q) for q in
                (model.parent_flux,model.child_flux,model.child_flux-model.charge)]
        D=values[0]-2*values[1]+values[2]
        result=infer_mixing_and_charge(o.tidal_quadratic_coefficient,D,kappa=model.scalar.kappa)
        assert result['epsilon']==pytest.approx(model.scalar.epsilon)
        assert result['charge_magnitude']==pytest.approx(model.charge)
    with pytest.raises(ValueError):infer_mixing_and_charge(8/3,-.1)
