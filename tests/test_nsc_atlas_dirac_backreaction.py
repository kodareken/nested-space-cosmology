"""Independent charged-family and conserved-mass reconstruction controls."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_mean import RadarPointEvaluator
from recursive_horizons.nsc_atlas_dirac_gravity import null_stress_from_orthonormal
from recursive_horizons.nsc_atlas_dirac_backreaction import (
    MeanBackreaction,ef_lapse_v,mass_one_form)


class MaxwellControl:
    """Derivative in equal squared charge, ONLY a finite-domain exact control."""
    def __init__(self,chart,*,region='nested',strength=1.):
        self.chart,self.region,self.strength = chart,region,strength
        self.points = RadarPointEvaluator(chart,region=region)
        self.breakpoints = np.array([])
        self.start = 0.

    def evaluate(self,U,V,*,region=None):
        fields = self.points(U,V,region=region)
        rho = self.strength/fields['radius']**4
        pr,pt,j = -rho,rho,np.zeros_like(rho)
        null = null_stress_from_orthonormal(C=fields['conformal_factor'],rho=rho,
            radial_pressure=pr,outward_energy_flux=j,angular_pressure=pt)
        return dict(fields=fields,null=null,orthonormal=dict(
            energy_density=rho,radial_pressure=pr,angular_pressure=pt,
            outward_energy_flux=j))


def exact_charged_clock(chart,U,V,charge_squared,*,region='nested'):
    """Full nonlinear Einstein-Maxwell quiet family, separate from mean solver."""
    e = charge_squared
    beta_m = lambda r:.4/(2*r*r)+3*r/4
    beta_p = lambda r:beta_m(r)-2*r
    velocity = lambda r:np.sqrt(beta_m(r)**2-(1-r*r+e/r**2))
    wall = solve_ivp(lambda t,y:velocity(y),(0.,chart.duration),
        np.array([.45],dtype=complex),method='DOP853',rtol=2e-12,atol=2e-14,
        max_step=chart.duration/12,dense_output=True)
    assert wall.success
    cap = region=='cap' or U>=V
    beta = beta_m if cap else beta_p
    f = (lambda r:1-r*r+e/r**2) if cap else (lambda r:1-.8/r+e/r**2)
    fp = (lambda r:-2*r-2*e/r**3) if cap else (lambda r:.8/r**2-2*e/r**3)
    lapse = lambda u:1/(beta(wall.sol(u)[0])+velocity(wall.sol(u)[0]))
    r = wall.sol(V)[0]
    rv = (velocity(r)+beta(r))/2
    ray = solve_ivp(lambda u,y:np.array([-lapse(u)*f(y[0])/2,
            -lapse(u)*fp(y[0])*y[1]/2]),(V,U),np.array([r,rv]),method='DOP853',
        rtol=2e-12,atol=2e-14,max_step=chart.duration/20)
    assert ray.success
    radius = ray.y[0,-1]
    C = 2*lapse(U)*ray.y[1,-1]
    return radius,C


@pytest.fixture(scope='module')
def reference():
    chart = WallRadarChart()
    source = MaxwellControl(chart)
    model = MeanBackreaction(source,chart.duration,
        initial_masses=(-1/(2*.45),-1/(2*.45)),wall_rtol=2e-10,wall_atol=2e-12)
    return chart,source,model


def test_nonvacuum_mass_one_form_and_EF_constraint(reference):
    chart,source,model = reference
    u,v = chart.maps(.475)[0],chart.maps(.525)[0]
    for a,b in ((u,v),(v,u)):
        row = source.evaluate(a,b)
        fields = row['fields'];r = fields['radius']
        ru = (fields['radius_t']-fields['radius_x'])/2
        rv = (fields['radius_t']+fields['radius_x'])/2
        dm = mass_one_form(fields,row['null'])
        assert dm['u'] == pytest.approx(ru/(2*r*r),rel=2e-14)
        assert dm['v'] == pytest.approx(rv/(2*r*r),rel=2e-14)
        assert ef_lapse_v(fields,row['null']) == 0.
        _,mass,_ = model.ef_fields(a,b,order=8)
        assert mass == pytest.approx(-1/(2*r),abs=2e-9)


def test_full_charged_metric_and_clock_derivative(reference):
    chart,source,model = reference
    u,v = chart.maps(.475)[0],chart.maps(.525)[0]
    for a,b in ((u,v),(v,u)):
        measured = model.clock_response(a,b,order=8,ray_rtol=2e-10,ray_atol=2e-12)
        exact_r,exact_C = exact_charged_clock(chart,a,b,1e-6j)
        background = source.points(a,b)
        assert exact_r.real == pytest.approx(background['radius'],abs=3e-10)
        assert exact_C.real == pytest.approx(background['conformal_factor'],abs=3e-9)
        assert measured['radius_change'] == pytest.approx(exact_r.imag/1e-6,abs=3e-7)
        assert measured['log_radar_C_change'] == pytest.approx(
            exact_C.imag/(1e-6*background['conformal_factor']),abs=8e-7)
        assert abs(measured['mass_path_gap']) < 3e-9
        assert abs(measured['wall_clock_log_C_at_corner']) < 2e-14
    # The signed exterior beta crosses zero without a singular reconstruction.
    tstar = float(chart.maps(.16**(1/3))[0])
    wall = model.clock_response(tstar,tstar)
    assert abs(wall['log_radar_C_change']) < 2e-14


def test_smooth_cap_control_has_shared_observer_motion_not_a_membrane(reference):
    chart,_,nested = reference
    source = MaxwellControl(chart,region='cap')
    control = MeanBackreaction(source,chart.duration,
        initial_masses=(-1/(2*.45),-1/(2*.45)),shared_worldtube=nested,
        wall_rtol=2e-10,wall_atol=2e-12)
    u,v = chart.maps(.475)[0],chart.maps(.525)[0]
    measured = control.clock_response(u,v,ray_rtol=2e-10,ray_atol=2e-12)
    exact_r,exact_C = exact_charged_clock(chart,u,v,1e-6j,region='cap')
    background = source.points(u,v)
    assert measured['radius_change'] == pytest.approx(exact_r.imag/1e-6,abs=3e-7)
    assert measured['log_radar_C_change'] == pytest.approx(
        exact_C.imag/(1e-6*background['conformal_factor']),abs=8e-7)
    assert abs(measured['mass_path_gap']) < 3e-9
    assert control.wall_state(v)['displacement'] == nested.wall_state(v)['displacement']


def test_zero_source_has_no_mass_clock_or_radius_response():
    chart = WallRadarChart()
    source = MaxwellControl(chart,strength=0.)
    u,v = chart.maps(.475)[0],chart.maps(.525)[0]
    model = MeanBackreaction(source,v)
    result = model.clock_response(u,v)
    for key in ('radius_change','log_radar_C_change','frequency_fractional_change',
                'received_over_wall_frequency_change','mass','mass_path_gap'):
        assert result[key] == 0.


def test_nonzero_lapse_source_against_exact_static_fluid_metric():
    # A second exact family exercises psi_r, absent in Maxwell's T_rr=0.
    # rho=e/kappa, p_r=0, F=1-(1+e/3)r² and
    # psi=-e*log(F)/(4*(1+e/3)); p_t follows static conservation.
    # The comparison worldtube is driven along the original areal trajectory.
    chart = WallRadarChart()
    class FixedWorldtube:
        def wall_state(self,t):
            return dict(displacement=0.,displacement_rate=0.)
    class FluidControl(MaxwellControl):
        def evaluate(self,U,V,*,region=None):
            fields = self.points(U,V,region='cap')
            r,N,f = fields['radius'],fields['lapse'],fields['static_f']
            boost = fields['radius_t']/(N*np.sqrt(f))
            gamma = np.sqrt(1+boost*boost)
            rho,pr,j = gamma*gamma,boost*boost,-gamma*boost
            pt = -r*r/(2*f)
            null = null_stress_from_orthonormal(C=fields['conformal_factor'],
                rho=rho,radial_pressure=pr,outward_energy_flux=j,angular_pressure=pt)
            return dict(fields=fields,null=null,orthonormal=dict(
                energy_density=rho,radial_pressure=pr,angular_pressure=pt,outward_energy_flux=j))
    source = FluidControl(chart,region='cap')
    model = MeanBackreaction(source,chart.duration,initial_masses=(.45**3/6,.45**3/6),
        shared_worldtube=FixedWorldtube(),wall_rtol=2e-10,wall_atol=2e-12)
    U,V = chart.maps(.475)[0],chart.maps(.525)[0]
    measured = model.clock_response(U,V,ray_rtol=2e-10,ray_atol=2e-12)
    row = source.evaluate(U,V);r = row['fields']['radius']
    assert ef_lapse_v(row['fields'],row['null'])>0
    assert measured['mass'] == pytest.approx(r**3/6,abs=2e-9)
    assert abs(measured['mass_path_gap'])<3e-9
    beta = lambda r:.4/(2*r*r)+3*r/4
    velocity = lambda r:np.sqrt(beta(r)**2-(1-r*r))
    wall = solve_ivp(lambda t,y:velocity(y),(0.,chart.duration),[.45],
        method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
    e=1e-6j
    F = lambda r:1-(1+e/3)*r*r
    Fp = lambda r:-2*(1+e/3)*r
    psi = lambda r:-e*np.log(F(r))/(4*(1+e/3))
    psip = lambda r:e*r/(2*F(r))
    def a(u,r):
        R = wall.sol(u)[0];v = velocity(R)
        b = np.sqrt(v*v+F(R))
        return np.exp(psi(r)-psi(R))/(b+v)
    R = wall.sol(V)[0];v = velocity(R);b=np.sqrt(v*v+F(R))
    ray = solve_ivp(lambda u,y:np.array([-a(u,y[0])*F(y[0])/2,
            -a(u,y[0])*(Fp(y[0])+psip(y[0])*F(y[0]))*y[1]/2]),
        (V,U),np.array([R,(v+b)/2],dtype=complex),method='DOP853',
        rtol=2e-12,atol=2e-14,max_step=chart.duration/20)
    assert ray.success
    exact_r = ray.y[0,-1]
    exact_C = 2*a(U,exact_r)*ray.y[1,-1]
    assert measured['radius_change'] == pytest.approx(exact_r.imag/1e-6,abs=3e-7)
    assert measured['log_radar_C_change'] == pytest.approx(
        exact_C.imag/(1e-6*row['fields']['conformal_factor']),abs=8e-7)
