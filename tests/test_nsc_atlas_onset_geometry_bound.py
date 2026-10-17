"""Independent algebra and interval-integral controls, without a physical run."""
from fractions import Fraction as Q

from flint import arb, ctx
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_onset_geometry_bound as m

CONTROLS = ControlBudget(maximum_evaluations=30000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('onset geometry controls:', CONTROLS.budget.evaluations,
          'charged values;', CONTROLS.active_cpu_seconds, 'active CPU seconds')


def test_mass_fluxes_from_tangent_normal_projections(monkeypatch, budget):
    budget.check(3000)
    monkeypatch.setattr(m, '_ball', lambda q:s.Rational(q.numerator, q.denominator))
    z, v, bp, sx, w, wx = s.symbols('z v bp sx w wx', nonzero=True)
    p, j, bm = w/z, wx/sx, bp+2*z
    A = (j-v*p)/z
    def kinetic(beta):
        B = (j-beta*p)/z
        return z*z*(beta*A*B-v*(A*A+B*B)/2)/2
    lp, ld, qp, qd = m.mass_terms(z, v, bp, sx, w, wx)
    # Substitute the independent reference normalization only after projection.
    def reduce_reference(expr):
        return s.factor(s.rem(s.together(expr).as_numer_denom()[0],
                             v*v-bp*bp+1-s.Rational(4,5)/z,v))
    assert reduce_reference(kinetic(bp)-lp) == 0
    assert s.factor(kinetic(bp)-kinetic(bm)-z*z*v*p*p/2-ld) == 0
    Bp, Bm = (j-bp*p)/z, (j-bm*p)/z
    qplus = s.Rational(2, 5)*p*p/12+z*z*p*(v*A-bp*Bp)/6
    qminus = z*z*p*(v*A-bm*Bm)/6-z**3*p*p/12
    assert reduce_reference(qplus-qp) == 0
    assert s.factor(reduce_reference(qplus-qminus-qd).subs(bp,1/(5*z*z)-5*z/4)) == 0


def test_radius_variation_and_frame_cancellation(budget):
    budget.check(3000)
    z,v,vz,bp,sx,w,j,ip,idiff,e = s.symbols('z v vz bp sx w j ip idiff e', nonzero=True)
    cp, cd = 1/(6*z)-1/(10*z*z), z/4-1/(10*z*z)
    mu, D = ip+(v-bp)*w*j/6+cp*w*w, idiff+z*w*j/3+cd*w*w
    rho = e+w*w/(12*z)
    rho_s = vz*rho+(mu/z+bp*D/(2*z*z))/v
    e_s = rho_s-w*j/(6*z)+v*w*w/(12*z*z)
    independent = sx*(e_s/v-e*vz/v)
    c = dict(z=z,v=v,bp=bp,sx=sx,
             Ce=vz/(12*z*v)+(cp/z+bp*cd/(2*z*z))/(v*v)+1/(12*z*z))
    assert s.factor(independent-m.radius_rate(c,ip,idiff,w)) == 0
    assert s.diff(m.radius_rate(c,ip,idiff,w),j) == 0


def test_cap_mass_lapse_density_and_clock_cancellations(budget):
    budget.check(3000)
    r,z,qc,p = s.symbols('r z qc p', nonzero=True)
    nu_r = r*r*p*p/2+r*r*(1-r*r)*qc/4
    delta_r = r*qc/2
    independent = 2*nu_r/z-(1-z*z)*delta_r
    factored = r*qc/(2*z)*(r-z)*(1-r*r-r*z-z*z)+r*r*p*p/z
    assert s.factor(independent-factored) == 0
    v,bm,bmz,sx,w,j,ip,e,rho,nu,delta = s.symbols('v bm bmz sx w j ip e rho nu delta')
    p, az = w/z,1-z*z
    cd = z/4-1/(10*z*z)
    D = ip+z*w*j/3+cd*w*w
    dbm = bmz*rho+D/(2*z*z)
    B = (j-bm*p)/z
    A1 = -2*z*e-2*nu/z
    direct = sx*((dbm-z*p*B/6-bm*p*p/12)/az-bm*delta/az-bm*A1/az**2)
    c = dict(z=z,bm=bm,bm_z=bmz,sx=sx,Az=az)
    reduced = m.cap_clock_rate(c,e,rho,ip,w,2*nu/z-az*delta)
    assert s.factor(direct-reduced) == 0
    assert s.diff(reduced,j) == 0


def test_eta_endpoint_and_history_cancellations(budget):
    budget.check(3000)
    z,d,v,bp,w,j,ip,idiff = s.symbols('z d v bp w j ip idiff', nonzero=True)
    cp,cd = 1/(6*z)-1/(10*z*z), z/4-1/(10*z*z)
    mu,D = ip+d*w*j/6+cp*w*w, idiff+z*w*j/3+cd*w*w
    assert s.factor(2*z*mu-d*D-(2*z*ip-d*idiff+(2*z*cp-d*cd)*w*w)) == 0
    zt,dt = s.symbols('zt dt',nonzero=True)
    assert s.factor(d*zt-z*dt-z*zt*(d/z-dt/zt)) == 0
    rho,vz,bpz = s.symbols('rho vz bpz')
    rho_s = vz*rho+(mu/z+bp*D/(2*z*z))/v
    dbp = bpz*rho+D/(2*z*z)
    target = (vz-bpz)*rho+(2*z*mu-d*D)/(2*z*z*v)
    assert s.factor((rho_s-dbp-target).subs(d,v-bp)) == 0


def test_source_logarithmic_derivative_and_collar_sign(budget):
    budget.check(1000)
    u = s.symbols('u',positive=True)
    logarithmic = -s.diff((1-u*u)**4,u)/(1-u*u)**4
    assert s.factor(logarithmic-8*u/(1-u*u)) == 0
    assert s.factor(8*u/(1-u*u)-s.Rational(16,3)-8*(2*u-1)*(u+2)/(3*(1-u*u))) == 0
    # u=-compact_coordinate is in [1/2,1), coordinate derivative >=1.


def test_moment_triangle_includes_nonconstant_exact_primitive(budget):
    budget.check(1000)
    with ctx.workprec(160):
        h = m._ball(Q(1,4))
        # density(t)=2+3t, first(0)=3, second(0)=5.
        density = arb(2).union(m._ball(Q(11,4)))
        F,G,endF,endG = m.moment_panel(arb(3),arb(5),density,h)
        for k in range(17):
            t = Q(k,64)
            f, g = 3+2*t+Q(3,2)*t*t, 5+3*t+t*t+Q(1,2)*t**3
            assert F.contains(m._ball(f)) and G.contains(m._ball(g))
        assert endF.contains(m._ball(f)) and endG.contains(m._ball(g))
        exact = m.moment_panel(arb(3),arb(5),arb(2),h)
        assert exact[2].contains(m._ball(Q(7,2)))
        assert exact[3].contains(m._ball(Q(93,16)))


def manufactured_gauge():
    def row(lo,hi=None):
        if hi is None:hi=lo
        def dy(q):
            q=Q(q); den=q.denominator
            assert den & (den-1) == 0
            return dict(mantissa=str(q.numerator),exponent=-(den.bit_length()-1))
        return dict(lower=dy(lo),upper=dy(hi))
    return dict(inputs=dict(left='7/16',right=str(Q(7,16)+Q(1,10000)),
                            radius_lower='15/32',radius_upper='15/32'),
        wall={k:row(v) for k,v in dict(s_x=1,v=1,beta_plus=Q(1,4),
            beta_minus=Q(5,4),T_x_plus_Z_x=4,W=Q(1,8)).items()},
        vacuum_wedge=dict(r=row(Q(7,16),Q(15,32)),r_U=row(Q(1,4)),K=row(2)))


def test_quiet_panel_and_hull_guards(budget):
    budget.check(2048)
    with ctx.workprec(160):
        gauge = manufactured_gauge()
        left,right=Q(gauge['inputs']['left']),Q(gauge['inputs']['right'])
        c=m.geometry_constants(gauge,left,right)
        quiet={k:arb(0) for k in m.STATE_KEYS}
        new=m._step(quiet,c,m._ball(left).union(m._ball(right)),m._ball(right-left),arb(0),arb(0))
        assert all(v==0 for v in new.values())
        gauge['wall']['s_x']['lower']['mantissa']='-1'
        with pytest.raises(ValueError,match='positive finite sx'):
            m.geometry_constants(gauge,left,right)


@pytest.mark.parametrize('panels',[0,129,True,1.5])
def test_partition_guard(panels,budget):
    with pytest.raises(ValueError,match='fixed partition'):
        m.onset_geometry_bound({}, {}, {}, budget=budget, panels=panels)


def test_budget_stops_before_source_or_geometry_work(budget):
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.onset_geometry_bound({}, {}, {}, budget=Spent())
