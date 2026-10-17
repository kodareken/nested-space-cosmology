"""Independent second-derivative and comparison controls, no original run."""
from fractions import Fraction as Q

from flint import arb,ctx
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_geometric_remainder_bound as m
from recursive_horizons.nsc_atlas_geometric_departure_bound import affine_image

CONTROLS=ControlBudget(maximum_evaluations=30000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        with ctx.workprec(160):yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('geometric remainder controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def contains(value,exact):
    exact=Q(str(s.simplify(exact)))
    assert m.dyadic(value.lower())<=exact<=m.dyadic(value.upper())


def test_jet_composition_against_polynomial_derivatives(budget):
    budget.check(1024)
    t=s.symbols('t');p=2+t+t*t;q=3-2*t+3*t*t
    pj=m.Jet(2,1,2);qj=m.Jet(3,-2,6)
    for actual,expected in [(pj*qj,p*q),(pj/qj,p/q),(pj**3,p**3),((-pj/2).exp(),s.exp(-p/2))]:
        # exp(-1) is evaluated by an independent Arb constant for this case.
        for jet,order in [(actual.value,0),(actual.first,1),(actual.second,2)]:
            target=s.diff(expected,t,order).subs(t,0)
            if not target.is_Rational:
                coefficient=s.simplify(target*s.E)
                assert jet.overlaps(m.ball(Q(str(coefficient)))*arb(-1).exp())
            else:contains(jet,target)


def test_even_product_keeps_one_sixth_and_second_field_jet(budget):
    budget.check(1024)
    alpha,ch,sh,q=m.even_product(Q(0),m.Jet(2,3,5))
    contains(q.value,0);contains(q.first,Q(2,3));contains(q.second,4)
    contains(ch.first,Q(1,3));contains(ch.second,Q(55,27))
    contains(sh.first,Q(1,9));contains(sh.second,Q(91,135))
    _,_,_,quiet=m.even_product(Q(1,10),m.Jet(2,0,0))
    contains(quiet.second,0)


def test_even_zero_limits_and_nonnegative_declared_domain(budget):
    budget.check(1024)
    _,ch,sh,_=m.even_product(Q(0),m.Jet(1,0,0))
    contains(ch.value,1);contains(ch.first,Q(1,12));contains(ch.second,Q(1,432))
    contains(sh.value,1);contains(sh.first,Q(1,36));contains(sh.second,Q(1,2160))
    alpha,ch,sh,q=m.even_product(Q(1,10),m.Jet(1))
    assert alpha.value.contains(m.ball(Q(1,10)))
    assert ch.value.contains(1) and sh.value.contains(1)
    # A known interior positive point from the first three positive terms.
    assert ch.value.contains(m.ball(1+Q(1,120)+Q(1,86400)))
    with pytest.raises(ValueError,match='nonnegative coupling'):
        m.even_product(Q(-1,10),m.Jet(1))


def test_bulk_jet_uses_independent_transverse_input(budget):
    budget.check(1024)
    t=s.symbols('t');expr=((Q(1,4)+t/7)*(Q(1,2)-t/9)+s.exp(t/11)/4)/(2+t/5)**2
    actual=m.bulk_kernel(m.Jet(2,Q(1,5)),m.Jet(Q(1,4),Q(1,7)),
                          m.Jet(Q(1,2),Q(-1,9)),m.Jet(0,Q(1,11)))
    for jet,order in [(actual.value,0),(actual.first,1),(actual.second,2)]:
        contains(jet,s.diff(expr,t,order).subs(t,0))
    changed=m.bulk_kernel(m.Jet(2),m.Jet(Q(1,4),1),m.Jet(Q(1,2)),m.Jet(0))
    contains(changed.first,Q(1,8));contains(changed.second,0)


def wall_inputs():
    values=dict(R=(Q(1,2),Q(1,17)),L=(0,Q(1,19)),qr=(Q(1,2),Q(-1,23)),
        X=(Q(1,16),Q(1,29)),Qx=(Q(1,8),Q(-1,31)),Em=(0,Q(1,37)),
        N=(Q(1,64),Q(1,41)),D=(Q(1,32),Q(1,43)))
    return values


def test_wall_jets_against_independent_symbolic_alpha_zero_map(budget):
    budget.check(4096)
    t=s.symbols('t');v={k:s.Rational(a)+t*s.Rational(b) for k,(a,b) in wall_inputs().items()}
    R,L,qr,X,Qx,Em,N,D=[v[k] for k in ('R','L','qr','X','Qx','Em','N','D')]
    f=X/R;q=t*f*f/6
    ch=1+q/2+q*q/24;sh=1+q/6+q*q/120
    Up=s.exp(-L/2)/2;Vp=2*s.exp(-L/2);Sw=2*ch**3
    A0=1-R*R;A=A0-2*t*N/R;mu=s.exp(t*D);H=2*Vp*qr+Sw*R;gap=A/H;bb=(2*Vp*qr+gap)/2
    K=-3*ch/4-t*f*sh/(2*R*ch)*(Vp*Qx-gap*Em/(2*A0)-f*bb)-2*bb/R-ch**3
    expected=dict(r_U=(Sw*R-gap)/(2*Up),L_U=4*((1+t)/7-2*K/Vp),T_x=(H+gap)/(2*A*mu*Vp))
    actual,_=m.wall_ports(**{k:m.Jet(a,b) for k,(a,b) in wall_inputs().items()},
                          Wx=Q(1,4),alpha_max=Q(0),qL=m.Jet(Q(1,7),Q(1,7)))
    for key in expected:
        for jet,order in [(actual[key].value,0),(actual[key].first,1),(actual[key].second,2)]:
            contains(jet,s.diff(expected[key],t,order).subs(t,0))


def test_signed_normal_trace_cancellation(budget):
    budget.check(1024)
    R,Up,Vp,qr,Qx,Em,A0,gap,Sw,root,u,f,alpha=s.symbols('R Up Vp qr Qx Em A0 gap Sw root u f alpha',nonzero=True)
    frame=u-root*f;J=Sw*R*(f-u/root);P=(gap*Em/A0+J)/(2*Up);pr=(Sw*R-gap)/(2*Up)
    cp=2*Vp*(root*Qx+frame*qr);cm=2*Up*(root*P+frame*pr);bb=Vp*qr+gap/2
    nc=(cp-cm)/2
    assert s.simplify(nc-u*bb-root*(Vp*Qx-gap*Em/(2*A0)-f*bb))==0
    old=-3*root/4-alpha*u*(nc-u*bb)/(2*R*root**2)-2*bb/R-Sw/2
    new=-3*root/4-alpha*u/(2*R*root)*(Vp*Qx-gap*Em/(2*A0)-f*bb)-2*bb/R-Sw/2
    assert s.simplify(old-new)==0


def test_affine_wall_lapse_slots_have_no_hessian(budget):
    budget.check(1024)
    fixed={k:m.Jet(a) for k,(a,_) in wall_inputs().items()}
    a,_=m.wall_ports(**fixed,Wx=Q(1,4),alpha_max=Q(0))
    b,_=m.wall_ports(**fixed,Wx=Q(1,4),alpha_max=Q(0),qL=m.Jet(2,3),Wxx=m.Jet(5,7))
    contains(b['L_U'].value-a['L_U'].value,-72)
    contains(b['L_U'].first-a['L_U'].first,-100)
    contains(b['L_U'].second-a['L_U'].second,0)


def test_homogeneous_image_reuses_existing_matrix(budget):
    budget.check(2048)
    c={k:arb(1) for k in ('dt','W','Wxhi','Wxlo','Upmin','Upmax','Vpmin','Vpmax','Amin','Hmin',
        'rx0hi','Sw','drt','dSw','u_wall','nc','bbar','mu_coefficient','N','vmax','pabs','k','Cmax','kinetic')}
    c.update(rlo=m.ball(Q(1,4)),rhi=m.ball(Q(1,2)),Rlo=m.ball(Q(1,2)),R0lo=m.ball(Q(1,2)),R0hi=m.ball(Q(1,2)))
    c['source']={k:arb(1) for k in ('r','log_C','log_C_x','r_times_r_x_minus_x')}
    zero={k:Q(0) for k in c['source']}
    nonlinear=dict(hessians={k:0 for k in ('r_U','L_U','T_x','K')},
                   terms={k:0 for k in ('radius_longitudinal','kinetic','radius_kernel_cross')})
    vector=(Q(2,3),Q(3,5),Q(5,7))
    before=affine_image(c,vector,0);after=m.remainder_image(c,zero,nonlinear,vector)
    assert {k:m._round_upper(v) for k,v in before.items()}=={k:m._round_upper(v) for k,v in after.items()}


def test_strict_sign_rejects_threshold_equality(budget):
    first={name:dict(lower=-2,upper=-1) for name in ('radius_fixed_source','clock_fixed_source')}
    result=m.endpoint_comparison(first,dict(radius=2,delta_T=0),Q(1,2),budget=budget)
    assert result['radius']['strict_on_full_interval'] is False
    assert result['radius']['strict_sign_alpha_max']==Q(1,4)
    assert result['radius']['normalized_upper_at_maximum']==0
    assert result['delta_T']['strict_on_full_interval'] is True


def test_positive_reciprocal_and_budget_guards(budget):
    budget.check(512)
    with pytest.raises(ValueError,match='positive reciprocal'):
        m.Jet(arb(0).union(arb(1))).inverse()
    class Spent:
        def check(self,count=0):raise RuntimeError('spent')
    with pytest.raises(RuntimeError,match='spent'):
        m.geometric_remainder_bound({}, {}, {}, {}, {}, {}, {}, {},budget=Spent())
