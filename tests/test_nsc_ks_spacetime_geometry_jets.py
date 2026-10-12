"""Clock, endpoint germs and genuine mixed high geometry derivatives."""
import numpy as np
import pytest
pytest.importorskip('flint')
from flint import arb,ctx
from recursive_horizons import nsc_ks_ball_geometry as G
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily
from recursive_horizons.nsc_ks_spacetime_geometry_jets import background_coordinate_time_jets,spacetime_geometry_jets
from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric,CompatibleRadiusDirection
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction


def family(alpha=.001):
    w=LocalAxialFunction((.2,.3,-.1,.07),1.2)
    U=LocalAxialFunction((.1,-.03,.05),1.2)
    return CompatibleIncomingMetric((alpha,),(CompatibleRadiusDirection(w,U,.007,.03),))


def test_coordinate_time_chain_rule_keeps_the_reference_clock():
    with ctx.workprec(160),G._JetWork(9):
        rho=arb(1.01);base=G.background_series(rho,9)
        jets=background_coordinate_time_jets(rho,9)
        assert jets['rho_increment'][0]==0
        assert jets['rho_increment'][1].overlaps(-base.a[0])
        assert jets['a'][1].overlaps(-base.a[1]*base.a[0])
        assert jets['r_reference'][1].overlaps(-rho/base.r[0]*base.a[0])
        # This nonconstant background cannot be relabelled as rho derivatives.
        assert not jets['r_reference'][1].overlaps(base.r[1])
        assert jets['a'][9].is_finite() and not jets['a'][9].contains(0)


def test_incoming_germ_and_mixed_derivative_factorials():
    f=family();z=arb(1.21)
    jet=spacetime_geometry_jets(f,arb(1),z,order=9)
    with ctx.workprec(160),G._JetWork(9):
        w,U=AnalyticRadiusFamily(f).profile_series(z,9)
        for j in range(7):
            expected=w[j]*__import__('math').factorial(j)
            actual=jet.derivative('r',1,j)-jet.derivative('r_reference',1,j)
            assert actual.overlaps(expected)
            expected_U=U[j]*__import__('math').factorial(j)
            actual_U=jet.derivative('r',3,j)-jet.derivative('r_reference',3,j)
            assert actual_U.overlaps(expected_U)
        assert jet.derivative('r',0,0).overlaps(jet.derivative('r_reference',0,0))
        assert (jet.derivative('r',2,2)-jet.derivative('r_reference',2,2)).contains(0)
    with pytest.raises(ValueError,match='retained'):
        jet.derivative('r',5,5)
    assert spacetime_geometry_jets(family(.002),1.,z).history_identity!=jet.history_identity


def test_transition_mixed_jets_match_independent_high_precision_first_derivatives():
    import mpmath as mp
    f=family();rho0=1.012;z=1.245
    jet=spacetime_geometry_jets(f,rho0,z,order=9)
    with mp.workdps(65):
        def axial(r):return mp.sqrt(3*((1+r*r)*(mp.pi/2-mp.atan(r))-r)-1)
        r0=mp.mpf(rho0);s0=-mp.quad(lambda r:1/axial(r),[1,r0])
        w=float(f.directions[0].w(z,0))*f.amplitudes[0]
        U=float(f.directions[0].U(z,0))*f.amplitudes[0]
        inner=mp.mpf(f.directions[0].inner_radius);outer=mp.mpf(f.directions[0].outer_radius)
        def radius(t):
            rho=mp.findroot(lambda r:-mp.quad(lambda x:1/axial(x),[r0,r])-t,r0)
            s=s0+t;u=(-s-inner)/(outer-inner)
            chi=1/(1+mp.exp(1/(1-u)-1/u))
            return mp.sqrt(1+rho*rho)+chi*(s*w+s**3*U/6)
        for n in (1,2):
            expected=float(mp.diff(radius,0,n))
            observed=float(jet.derivative('r',n,0))
            assert abs(observed-expected)<2e-12*max(1,abs(expected))
    assert jet.derivative('r',4,5).is_finite()
    assert not jet.derivative('r',4,5).contains(0)


def test_geometry_box_contains_mixed_derivatives_at_all_sampled_corners():
    f=family()
    with ctx.workprec(160):
        rho=arb(1.012,1e-6);z=arb(1.245,1e-6)
        enclosure=spacetime_geometry_jets(f,rho,z,order=9)
        for r in (1.0119995,1.0120005):
            for x in (1.2449995,1.2450005):
                point=spacetime_geometry_jets(f,r,x,order=9)
                for t,s in ((0,0),(1,1),(1,4),(4,5)):
                    assert enclosure.derivative('r',t,s).contains(point.derivative('r',t,s))
