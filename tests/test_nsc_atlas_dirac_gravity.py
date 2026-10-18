"""Independent curvature, Bianchi, wall-normal and gauge controls."""
import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart
from recursive_horizons.nsc_atlas_dirac_gravity import (
    constraint_transport, einstein_from_null_jets, linearized_null_rhs,
    null_stress_from_orthonormal, quiet_null_jets, quiet_wall_jets,
    radar_observer_variation, wall_metric_variation)


def test_full_four_dimensional_einstein_contraction():
    u,v,theta,phi = sp.symbols('u v theta phi',real=True)
    coords = (u,v,theta,phi)
    at0 = {u:0,v:0,theta:sp.pi/2,phi:0}
    L = sp.log(sp.Rational(3,2))+2*u/7-v/5+u*u/17+u*v/11-v*v/23
    r = sp.Rational(4,5)+u/3+v/4+u*u/13+u*v/19+v*v/29
    C = sp.exp(L)
    metric = sp.zeros(4)
    metric[0,1] = metric[1,0] = -C/2
    metric[2,2],metric[3,3] = r*r,r*r*sp.sin(theta)**2
    g = metric.subs(at0);inv = g.inv()
    dg = [metric.diff(x).subs(at0) for x in coords]
    ddg = [[metric.diff(x,y).subs(at0) for y in coords] for x in coords]
    di = [-inv@d@inv for d in dg]
    gamma = [[[sum(inv[a,d]*(dg[b][d,c]+dg[c][d,b]-dg[d][b,c])/2
        for d in range(4)) for c in range(4)] for b in range(4)] for a in range(4)]
    dgamma = [[[[sum((di[m][a,d]*(dg[b][d,c]+dg[c][d,b]-dg[d][b,c])
        +inv[a,d]*(ddg[m][b][d,c]+ddg[m][c][d,b]-ddg[m][d][b,c]))/2
        for d in range(4)) for c in range(4)] for b in range(4)] for a in range(4)] for m in range(4)]
    ricci = sp.Matrix(4,4,lambda b,c:sum(dgamma[a][a][b][c]-dgamma[c][a][b][a]
        +sum(gamma[a][a][d]*gamma[d][b][c]-gamma[a][c][d]*gamma[d][b][a]
             for d in range(4)) for a in range(4)))
    scalar = sum(inv[a,b]*ricci[a,b] for a in range(4) for b in range(4))
    tensor = ricci-g*scalar/2
    jets = dict(r=r,C=C,ru=r.diff(u),rv=r.diff(v),ruu=r.diff(u,2),
        ruv=r.diff(u,v),rvv=r.diff(v,2),Lu=L.diff(u),Lv=L.diff(v),Luv=L.diff(u,v))
    actual = einstein_from_null_jets(**{k:x.subs(at0) for k,x in jets.items()})
    assert actual['uu'] == tensor[0,0]
    assert actual['uv'] == tensor[0,1]
    assert actual['vv'] == tensor[1,1]
    assert actual['pt'] == tensor[2,2]/r.subs(at0)**2


def test_exact_constraint_transport_with_arbitrary_nonvacuum_stress():
    r,C,ru,rv,Lu,Lv,ruu,rvv,lam,k = sp.symbols('r C ru rv Lu Lv ruu rvv lam k')
    a,b,c,p,au,av,bu,bv,cu,cv = sp.symbols('a b c p au av bu bv cu cv')
    F = -ru*rv/r-C/(4*r)+lam*C*r/4+k*r*b/2
    G = 2*ru*rv/r**2+C/(2*r**2)-k*b-k*C*p/2
    # Differentiate F directly; mixed equations substitute only after differentiation.
    dU = lambda f:sp.diff(f,r)*ru+sp.diff(f,C)*C*Lu+sp.diff(f,ru)*ruu+sp.diff(f,rv)*F+sp.diff(f,b)*bu
    dV = lambda f:sp.diff(f,r)*rv+sp.diff(f,C)*C*Lv+sp.diff(f,rv)*rvv+sp.diff(f,ru)*F+sp.diff(f,b)*bv
    Eu,Ev = ruu-Lu*ru+k*r*a/2,rvv-Lv*rv+k*r*c/2
    Qu = av+bu+2*rv*a/r+(2*ru/r-Lu)*b+C*ru*p/r
    Qv = cu+bv+2*ru*c/r+(2*rv/r-Lv)*b+C*rv*p/r
    exact_u = dU(F)-G*ru-Lu*F+k*(rv*a+r*av)/2
    exact_v = dV(F)-G*rv-Lv*F+k*(ru*c+r*cu)/2
    actual = constraint_transport(r=r,ru=ru,rv=rv,constraint_u=Eu,constraint_v=Ev,
        conservation_u=Qu,conservation_v=Qv,gravity_coupling=k)
    assert sp.expand(exact_u-actual['constraint_u_v']) == 0
    assert sp.expand(exact_v-actual['constraint_v_u']) == 0


def test_actual_quiet_einstein_linearization_and_null_stress_units():
    chart = WallRadarChart()
    stress = null_stress_from_orthonormal(C=1.3,rho=.7,radial_pressure=.2,
        outward_energy_flux=-.11,angular_pressure=.25)
    assert -4*stress['uv']/1.3+2*stress['pt'] == pytest.approx(0.,abs=2e-16)
    for region in ('cap','exterior'):
        jets = quiet_null_jets(chart.from_wall_radii(.48,.52,region=region))
        raw = einstein_from_null_jets(**{k:jets[k] for k in
            ('r','C','ru','rv','ruu','ruv','rvv','Lu','Lv','Luv')})
        assert raw['uu'] == pytest.approx(0.,abs=1e-14)
        assert raw['uv'] == pytest.approx(jets['Lambda']*jets['C']/2,abs=2e-14)
        assert raw['pt'] == pytest.approx(-jets['Lambda'],abs=2e-14)
        s,su,sv,q,qu,qv = .17,-.11,.08,-.13,.19,-.23
        tangent = linearized_null_rhs(jets,s=s,su=su,sv=sv,q=q,qu=qu,qv=qv,
            stress=stress,gravity_coupling=.7)
        # Complex-step the independently checked full Einstein tensor, not the
        # displayed mixed RHS. The vacuum term varies with the actual metric.
        h = 1e-30j
        varied = {k:jets[k] for k in ('r','C','ru','rv','ruu','ruv','rvv','Lu','Lv','Luv')}
        for key,variation in dict(r=s,ru=su,rv=sv,Lu=qu,Lv=qv,**tangent).items():
            varied[key] = jets[key]+h*variation
        varied['C'] = jets['C']*np.exp(h*q)
        tensor = einstein_from_null_jets(**varied)
        equations = dict(uu=tensor['uu']-.7*h*stress['uu'],
            uv=tensor['uv']-jets['Lambda']*varied['C']/2-.7*h*stress['uv'],
            vv=tensor['vv']-.7*h*stress['vv'],
            pt=tensor['pt']+jets['Lambda']-.7*h*stress['pt'])
        assert max(abs(x.imag/h.imag) for x in equations.values()) < 4e-14


def test_moving_wall_variation_from_normal_and_acceleration():
    R=.5
    wall = quiet_wall_jets(R)
    s,sx,q,qx = [np.array(x) for x in ((.13,-.07),(-.21,.17),(.08,-.11),(.19,.23))]
    w,wd,wdd = .04,-.03,.09
    actual = wall_metric_variation(R,s=s,sx=sx,q=q,qx=qx,
        displacement=w,velocity=wd,acceleration=wdd)
    def exact(epsilon):
        x = epsilon*w
        L = wall['Lx']*x+wall['Lxx']*x*x/2+epsilon*(q+qx*x)
        Lt = 2*wall['beta_RR']*wall['v']*x+epsilon*.15
        Lx = wall['Lx']+wall['Lxx']*x+epsilon*qx
        r = R+wall['beta']*x+wall['rxx']*x*x/2+epsilon*(s+sx*x)
        rt = wall['v']+wall['beta_R']*wall['v']*x+epsilon*.12
        rx = wall['beta']+wall['rxx']*x+epsilon*sx
        time,angle,clock = [],[],[]
        for i in range(2):
            C = np.exp(L[i])
            tangent = np.array([1.,epsilon*wd])
            metric = np.diag([-C,C])
            proper = np.sqrt(-tangent@metric@tangent)
            normal = np.array([epsilon*wd,1.])/proper
            gamma = np.array([[[Lt,Lx[i]],[Lx[i],Lt]],
                              [[Lx[i],Lt],[Lt,Lx[i]]]])/2
            acceleration = (np.array([0.,epsilon*wdd])
                +np.einsum('abc,b,c->a',gamma,tangent,tangent))/proper**2
            time.append(normal@metric@acceleration)
            angle.append(normal@np.array([rt[i],rx[i]])/r[i])
            clock.append(proper)
        return np.array(time),np.array(angle),np.array(clock)
    step=2e-6
    before,after = exact(-step),exact(step)
    for index,key in enumerate(('Ktime','Kangle','proper_clock_rate')):
        assert (after[index]-before[index])/(2*step) == pytest.approx(actual[key],abs=2e-9)


def test_wall_Gauss_sign_and_residual_null_gauge():
    for R in (.45,.5,.55):
        b = quiet_wall_jets(R)
        mean = np.mean(b['Ktime']+2*b['Kangle'])
        assert -2*b['k_wall']*mean == pytest.approx(3.,abs=2e-14)
        assert b['Ktime'][1]-b['Ktime'][0] == pytest.approx(-2.)
        assert b['Kangle'][1]-b['Kangle'][0] == pytest.approx(-2.)
        # Active residual null reparametrization. a is the common time
        # relabeling; w is the coordinate wall displacement, not a new source.
        a,ad,w,wd,wdd = .17,-.09,-.04,.025,-.07
        beta,br,br2,v = b['beta'],b['beta_R'],b['beta_RR'],b['v']
        s = a*v-w*beta
        q = 2*ad-2*w*br
        sx = ad*beta-wd*v+a*br*v-w*b['rxx']
        qx = 2*ad*br+2*a*br2*v-w*b['Lxx']-2*wdd
        result = wall_metric_variation(R,s=s,sx=sx,q=q,qx=qx,
            displacement=w,velocity=wd,acceleration=wdd)
        for key in ('radius_continuity','metric_continuity','Israel_time',
                    'Israel_angle','normal_force'):
            assert abs(result[key]) < 6e-15
        assert result['surface_radius'] == pytest.approx(np.full(2,a*v),abs=2e-15)
        assert result['proper_clock_rate'] == pytest.approx(np.full(2,ad),abs=2e-15)


def test_proper_wall_radar_readout_removes_residual_coordinate_changes():
    chart = WallRadarChart()
    background = quiet_null_jets(chart.from_wall_radii(.475,.525))
    u,v = chart.maps(.475)[0],chart.maps(.525)[0]
    f = lambda t:.13+.07*t+.11*t*t
    g = lambda t:-.05+.09*t-.08*t*t
    fp = lambda t:.07+.22*t
    gp = lambda t:.09-.16*t
    wall = lambda t:dict(time_shift=(f(t)+g(t))/2,
        displacement=(f(t)-g(t))/2,clock_rate=(fp(t)+gp(t))/2,
        velocity=(fp(t)-gp(t))/2)
    s = f(u)*background['ru']+g(v)*background['rv']
    q = f(u)*background['Lu']+g(v)*background['Lv']+fp(u)+gp(v)
    result = radar_observer_variation(background,s=s,q=q,wall_u=wall(u),wall_v=wall(v))
    for value in result.values():
        assert abs(value) < 2e-16
    # With already normalized wall radar coordinates, no endpoint gauge
    # contribution remains: delta ln(omega_rec/omega_wall)=-delta ln C/2.
    zero_wall = dict(time_shift=0.,displacement=0.,clock_rate=0.,velocity=0.)
    measured = radar_observer_variation(background,s=.03,q=.04,
        wall_u=zero_wall,wall_v=zero_wall)
    assert measured['received_over_wall_frequency_fraction'] == pytest.approx(-.02)
