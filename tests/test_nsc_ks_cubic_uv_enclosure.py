"""Directed triangular propagation, coordinate jets and original-history controls."""
import json
from pathlib import Path

import numpy as np
import pytest
from flint import acb,arb,arb_series,ctx

from recursive_horizons.nsc_ks_cubic_uv_enclosure import (
    CubicForcing,CubicGeometry,characteristic_distance,enclose_characteristic,
    taylor_volterra_cell,volterra_cell,
    _distance_integrand,
)
from recursive_horizons import nsc_ks_ball_geometry as G
from recursive_horizons.nsc_ks_cubic_uv_current import (
    geometry_jets,transport_massless_cubic_fields,cubic_pauli_densities,
)
from recursive_horizons.nsc_ks_current_uv_transport import HISTORY_RECORD,ARCHIVE_RHO_UP
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

ROOT=Path(__file__).resolve().parents[1]


def family():
    return LocalIncomingFamily(np.asarray(json.loads((ROOT/HISTORY_RECORD).read_text())['history']['coefficients']))


def test_phase_gradient_forcing_is_the_derivative_of_the_original_equation():
    import sympy as sp
    z=sp.symbols('z',real=True)
    a,ell=sp.symbols('a ell',real=True,nonzero=True)
    r=sp.Function('r',positive=True)(z)
    b=a*ell/(2*r)
    original=-ell**2/(2*r*r)
    assert sp.simplify(sp.diff(original,z)+4*b*sp.diff(b,z)/a**2)==0
    assert sp.simplify(sp.diff(original,z,2)+4*(sp.diff(b,z)**2+b*sp.diff(b,z,2))/a**2)==0


def test_distance_background_identity_and_interval_radius_preservation():
    with ctx.workprec(160):
        for rho in (arb(1),arb('1.015 +/- 0.00001'),arb(33)/32):
            direct=_distance_integrand(acb(rho),True)
            background=G.background_series(rho,0).inv_a2[0]
            assert direct.imag.contains(0)
            assert direct.real.overlaps(background)
            inside=3*((1+rho*rho)*(arb.pi()/2-rho.atan())-rho)-1
            assert inside>0
            assert (background*inside).contains(1)
        original=arb('1.2 +/- 0.0001')-characteristic_distance(arb('1.018 +/- 0.000002'))
        converted=G._binary_arb(original,'characteristic image')
        assert converted is original
        assert converted.rad()>0
        copied=arb(original)
        assert copied.contains(original) and original.contains(copied)
        assert arb(0,1).contains(-1) and arb(0,1).contains(1)
        assert arb('0.5','0.5').contains(0) and arb('0.5','0.5').contains(1)


def test_partial_cell_tube_encloses_nested_integral_in_decreasing_time():
    # qzz(t)=4+2t, J(t)=12t+3t^2, t in [-1,0].
    forcing=CubicForcing(arb(1),arb(2),arb(0),arb(3))
    out=volterra_cell((0,4,0),forcing,1,0)
    assert out[0].contains(-1) and out[1].contains(2)
    assert out[2].contains(-9)
    assert not arb(-12).contains(-9)  # freezing qzz omits the nested term
    dropped=volterra_cell((0,4,0),CubicForcing(arb(1),arb(2),arb(0),arb(0)),1,0)
    assert not dropped[2].contains(-9)


def test_taylor_cell_is_exact_for_a_polynomial_triangular_problem():
    # f2=2+t, k=3+t, qzz(0)=4. J'=12+10t+3.5t^2+0.5t^3.
    with ctx.workprec(160):
        point={k:arb_series(v,prec=5) for k,v in
               {'qz':[1],'qzz':[2,1],'current':[0],'coupling':[3,1]}.items()}
        ranges={k:arb_series(v,prec=5) for k,v in
                {'qz':[1],'qzz':[arb(1).union(2),1],
                 'current':[0],'coupling':[arb(2).union(3),1]}.items()}
        out=taylor_volterra_cell((0,4,0),point,ranges,1,0,4)
        assert out[0].contains(-1)
        assert out[1].contains(arb(5)/2)
        exact=-12+5-arb(7)/6+arb(1)/8
        assert out[2].contains(exact)
        assert out[2].rad() < arb('1e-40')


def test_coordinate_jets_match_existing_analytic_derivatives():
    f=family();model=CubicGeometry(f,np.sqrt(5))
    for rho,z in ((1.008,f.center+.012),(1.018,f.center-.041)):
        j=geometry_jets(f,rho,np.array([z]));a=j['a'];r=j['r'][0]
        b=a*np.sqrt(5)/(2*r)
        expected={'a':a,'a_rho':j['a_rho'],'b':b,
                  'b_z':-b*j['r_z'][0]/r,
                  'b_rho':b*(j['a_rho']/a-j['r_rho'][0]/r)}
        enclosed=model.jets(arb(rho),arb(z))
        for name,value in expected.items():
            # The binary owner is a diagnostic comparison; its roundoff does
            # not replace the ball owner's analytic enclosure.
            assert abs(float(enclosed[name])-value)<2e-12*(1+abs(value))


def test_characteristic_composition_matches_direct_uniform_forcing():
    f=family();model=CubicGeometry(f,np.sqrt(5))
    with ctx.workprec(160):
        rho=arb('1.018 +/- 0.000002');z=arb(f.center)
        for sign in (-1,1):
            direct=model.forcing(rho,z,sign)
            jets=model.forcing_series(rho,z,sign,4)
            for name,value in vars(direct).items():
                assert jets[name][0].overlaps(value)
            midpoint=arb('1.018');h=arb('0.0000001')
            plus=model.forcing(midpoint+h,z,sign)
            minus=model.forcing(midpoint-h,z,sign)
            for name in vars(direct):
                fd=(getattr(plus,name)-getattr(minus,name))/(2*h)
                assert jets[name][1].contains(fd)


def test_angular_pair_has_even_directed_forcing():
    f=family()
    plus,minus=(CubicGeometry(f,s*np.sqrt(5)) for s in (1,-1))
    with ctx.workprec(160):
        for sign in (-1,1):
            a=plus.forcing(arb('1.018 +/- 0.000001'),arb(f.center),sign)
            b=minus.forcing(arb('1.018 +/- 0.000001'),arb(f.center),sign)
            assert all(getattr(a,k).overlaps(getattr(b,k)) for k in vars(a))


@pytest.mark.parametrize('sign',[-1,1])
def test_original_nonzero_history_enclosure_contains_independent_transport(sign):
    f=family();model=CubicGeometry(f,np.sqrt(5))
    z=np.linspace(f.center-.0001,f.center+.0001,5)
    fields=transport_massless_cubic_fields(f,z,angular=np.sqrt(5),source_sign=sign,rho_steps=192)
    current=cubic_pauli_densities(geometry_jets(f,1.,z),fields,np.sqrt(5),sign)
    result=enclose_characteristic(model,arb(f.center),sign,ARCHIVE_RHO_UP,cells=256,order=8)
    assert result['J3'].contains(arb(float(current['J_I_3'][2])))
    reference=LocalIncomingFamily(np.zeros((2,32)))
    reference_fields=transport_massless_cubic_fields(reference,z,angular=np.sqrt(5),
        source_sign=sign,rho_steps=192)
    reference_current=cubic_pauli_densities(geometry_jets(reference,1.,z),
        reference_fields,np.sqrt(5),sign)
    delta_N=float(current['N_bracket_3'][2]-reference_current['N_bracket_3'][2])
    assert result['N_bracket3'].contains(arb(delta_N))
    assert result['J3'].rad() < arb('2e-6')
    assert result['higher_uv_remainder_bound'] is None
    assert result['physical_local_gate']=='OPEN'


def test_zero_history_encloses_zero_without_creating_a_physical_source_bound():
    result=enclose_characteristic(CubicGeometry(LocalIncomingFamily(np.zeros((2,8))),np.sqrt(5)),
        arb(1.2),-1,ARCHIVE_RHO_UP,cells=8,order=4)
    for key in ('q_z','q_zz','J3','N_bracket3'):
        assert result[key].contains(0)
        assert abs(result[key]).upper()<arb('1e-35')


def test_input_domains_and_precision_restoration():
    before=ctx.prec,ctx.cap
    model=CubicGeometry(family(),np.sqrt(5))
    assert characteristic_distance(arb(1)).is_zero()
    with pytest.raises(ValueError):
        model.forcing(arb(1.01),arb(1.2),0)
    with pytest.raises(ValueError):
        model.forcing_series(arb(1.01),arb(1.2),1,14)
    with pytest.raises(ValueError):
        enclose_characteristic(model,arb(1.2),1,arb('1.02 +/- .001'),cells=4)
    with pytest.raises(ValueError):
        enclose_characteristic(model,arb(1.2),1,ARCHIVE_RHO_UP,cells=True)
    with pytest.raises(ValueError,match='flat'):
        enclose_characteristic(model,arb(1.2),1,1.005,cells=4)
    assert (ctx.prec,ctx.cap)==before
