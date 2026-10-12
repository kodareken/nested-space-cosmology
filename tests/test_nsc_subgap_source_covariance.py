"""Current-preserving source covariance and directed defect controls."""
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp
from flint import acb,arb,arb_series,ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_massive_jost_phase_bound import _series_context
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame
from recursive_horizons.nsc_subgap_source_covariance import (
    BlochSource,metric_H,initial_bloch,capture_bloch,cell_defect,validate_bloch,
    original_covariance_error,
    original_covariance_distance_bounds,
)

ROOT=Path(__file__).resolve().parents[1]


def inputs():
    record=json.loads((ROOT/'results/development/nsc-metric-horizon-frame-v1.json').read_text())
    source=json.loads((ROOT/'results/development/nsc-ks-source-inventory.json').read_text())
    with np.load(ROOT/source['payload']['path'],allow_pickle=False) as data:
        config=json.loads(data['metadata_json'].tobytes())['config']
        F=data['group14/low16_1/amplitudes_at_one'][0].copy()
        C=data['group14/low16_1/covariance'][0].copy()
    ball=lambda v:restored_upper(v['lower']).union(restored_upper(v['upper']))
    with ctx.workprec(192):
        E=float.fromhex(record['energy_hex']);m=ball(record['mass']);ell=ball(record['angular_enclosure']) if 'angular_enclosure' in record else ball(record['angular'])
        frame=metric_horizon_frame(E,m,ell,config['horizon_rho'],bits=192)
        reflection=acb(ball(record['affine_reflection_direct_frame']['real']),
                       ball(record['affine_reflection_direct_frame']['imag']))
        start=float(ball(record['compact_distance']).log().mid())
        end=(frame.q-3*arb.pi()/4).log()
        initial=initial_bloch(frame,reflection,config['surface_gravity'],start)
        return BlochSource(frame.q,E,m,ell),initial,start,end,F,C


def test_pauli_commutator_has_the_proved_bloch_factor_and_sign():
    hx,hy,hz,nx,ny,nz=sp.symbols('hx hy hz nx ny nz',real=True)
    matrices=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))
    H=sum((x*s for x,s in zip((hx,hy,hz),matrices)),sp.zeros(2))
    Q=(sp.eye(2)+sum((x*s for x,s in zip((nx,ny,nz),matrices)),sp.zeros(2)))/2
    cross=(hy*nz-hz*ny,hz*nx-hx*nz,hx*ny-hy*nx)
    expected=sum((x*s for x,s in zip(cross,matrices)),sp.zeros(2))
    assert (-sp.I*(H*Q-Q*H)-expected).applyfunc(sp.expand)==sp.zeros(2)
    assert sp.expand(nx*cross[0]+ny*cross[1]+nz*cross[2])==0


def test_common_source_phase_and_interior_quarter_turn_preserve_coherence():
    x,y,f,c,t=sp.symbols('x y f c t',real=True)
    R=x+sp.I*y
    S=lambda r:sp.Matrix([[0,1,0],[r,0,t]])
    D=sp.diag(sp.exp(sp.I*c),sp.exp(-sp.I*c))
    Ds=sp.diag(sp.exp(sp.I*c),sp.exp(sp.I*c),sp.exp(-sp.I*c))
    assert (D*S(sp.exp(2*sp.I*c)*R)-S(R)*Ds).applyfunc(sp.simplify)==sp.zeros(2,3)
    coh,n=sp.symbols('coh n',real=True)
    C=sp.Matrix([[f,-sp.I*coh,0],[sp.I*coh,1-f,0],[0,0,n]])
    assert (Ds*C*Ds.H-C).applyfunc(sp.simplify)==sp.zeros(3)
    quarter=sp.diag(sp.exp(-sp.I*sp.pi/4),sp.exp(sp.I*sp.pi/4))
    sewn=(quarter*S(R).subs(t,0)*C*S(R).subs(t,0).H*quarter.H).applyfunc(sp.simplify)
    assert sp.simplify(sewn[0,1]-coh*sp.conjugate(R))==0
    assert sp.simplify(sewn[0,0]-(1-f))==0
    assert sp.simplify(sewn[1,1]-f*(x*x+y*y))==0


def test_stable_metric_tail_contains_closed_form_derivatives():
    with _series_context(512,6):
        q=arb('2.657264704679347')
        c=arb(3)/2*(2*q).sin()+arb(1)/2*(2*q).cos()
        s=arb(3)/2*(2*q).cos()-arb(1)/2*(2*q).sin()
        for center in ('-0.3','-0.1','0.001'):
            u=arb_series([arb(center),1],prec=7)
            closed=3-c*((2*u).cos()-1)/u-s*(2*u).sin()/u
            enclosed=metric_H(q,u,6,terms=16)
            for j in range(7):
                assert enclosed[j].contains(closed[j])


def test_actual_pg_restriction_preserves_the_working_covariance():
    from types import SimpleNamespace
    from recursive_horizons.nsc_common_ks_trace import restrict_resolved_modes
    from recursive_horizons.nsc_pg_massive_modes import MassivePGModeResolution
    from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain
    from recursive_horizons.nsc_lorentzian import geometry

    # frame() needs only the horizon position, not a new source preparation.
    resolution=object.__new__(MassivePGModeResolution)
    resolution.preparation=SimpleNamespace(horizon_rho=1.9006916054701435)
    energy=np.array([.0013248831260437577])
    columns=np.array([[1+.2j,.4j,-.3],[.1-.5j,.7,.2j]])
    source=np.array([[.4,.1j,0],[-.1j,.6,0],[0,0,.2]])
    expected=columns@source@columns.conj().T
    for rho in (1.,1.03):
        frame=resolution.frame(rho)[0]
        radius=np.sqrt(1+rho*rho)
        domain=TransmittingDiracSeamDomain(1.,1.,float(geometry(rho)[0]),radius)
        trace,_,_=domain.trace_map(np.ones(1))
        np.testing.assert_allclose(frame,radius*trace[0],rtol=3e-15,atol=3e-15)
        # Any real clock phase must disappear without changing source coherence.
        pg=(np.exp(.731j)*frame@columns)[None]
        canonical,_=restrict_resolved_modes(rho,energy,pg)
        actual=canonical[0]@source@canonical[0].conj().T
        np.testing.assert_allclose(actual,expected,rtol=4e-14,atol=4e-14)
        wrong=pg[0]@source@pg[0].conj().T
        assert np.linalg.norm(wrong-expected)>.01


class ConstantRotation:
    bits=160
    metric_terms=48
    def hamiltonian_series(self,y,order):
        return tuple(arb_series([v],prec=order+1) for v in (0,0,arb(1)/2))


def test_normalized_defect_is_not_multiplied_by_cell_width_twice():
    row=np.zeros(26);row[1]=.25;row[2]=row[5]=1
    error=cell_defect(ConstantRotation(),row)
    assert error.contains(arb(1)/4)
    actual=2*(arb(1)/8).sin()
    assert error>=actual
    assert error*arb('.25')<actual


def test_complete_source_covariance_bound_and_mutation_controls():
    model,initial,start,end,F,C=inputs()
    trace,_=capture_bloch(model,initial,start,float(end.mid()),max_step=.025)
    proof=validate_bloch(model,initial,trace,end)
    error=original_covariance_error(F,C,proof)
    assert error < arb('2e-11')
    assert proof['defect_integral'] < arb('2e-12')
    assert original_covariance_distance_bounds(F,np.diag(np.diag(C)),proof)['lower']>arb('.1')
    malformed=C.copy();malformed[0,1]+=.01
    with pytest.raises(ValueError,match='Hermitian'):
        original_covariance_error(F,malformed,proof)
    missing=trace.copy();missing[1,2]+=1e-5
    with pytest.raises(ValueError,match='join'):
        validate_bloch(model,initial,missing,end)
    assert proof['physical_local_gate']=='OPEN'


def test_point_initializer_encloses_unit_bloch_norm():
    _,initial,_,_,_,_=inputs()
    with ctx.workprec(192):
        assert sum((x*x for x in initial),arb(0)).contains(1)


def test_domains_and_flint_context_restored():
    before=ctx.prec,ctx.cap
    with pytest.raises(ValueError,match='subgap'):
        BlochSource(2.65,2,1,2)
    with _series_context(160,8):
        with pytest.raises(ValueError,match='tail order'):
            metric_H(arb(2.65),arb_series([-.01,1],prec=9),8,terms=8)
    model,initial,start,end,_,_=inputs()
    trace,_=capture_bloch(model,initial,start,float(end.mid()))
    with pytest.raises(ValueError):
        cell_defect(model,trace[0],degree=5)
    assert (ctx.prec,ctx.cap)==before
