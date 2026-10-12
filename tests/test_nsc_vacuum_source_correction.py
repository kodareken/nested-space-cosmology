"""Check the driven covariance correction against its original Bloch equation."""
from pathlib import Path

import numpy as np
import pytest
from flint import arb,arb_series,ctx

from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_massive_jost_phase_bound import _series_context
from recursive_horizons.nsc_subgap_source_covariance import _cross,original_covariance_distance_bounds
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import (
    VacuumCorrection,capture_correction,correction_cell_defect,validate_correction,
)

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def pilot():
    archive=RetainedUpstreamArchive(ROOT)
    entries=[b for b,_ in archive.family_entries((14,1)) if b.energy_sign>0]
    E,b,index=min(((float(E),b,i) for b in entries for i,E in enumerate(b.source.energies[::3])
                  if E>=16),key=lambda row:row[0])
    with ctx.workprec(192):
        expansion=VacuumSourceExpansion(archive.meta['config']['horizon_rho'],
            arb(b.mass).union(arb.pi()/2),arb(b.angular).union(arb(5).sqrt()),order=8)
        flow=VacuumCorrection(expansion,E)
        end=float(expansion.target_distance(b.rho_up).log().mid())
    trace,nfev=capture_correction(flow,-18.,end)
    proof=validate_correction(flow,trace,-18.,b.rho_up)
    return flow,trace,proof,b,index,nfev


def test_defect_jets_equal_the_original_equation_residual(pilot):
    flow,*_=pilot
    with _series_context(192,4):
        y=arb_series([arb('-1.3'),1],prec=5);delta=y.exp()
        _,ps,zs=flow.expansion.jets(delta[0],remainder_jet_order=4)
        increment=arb_series([0]+[delta[j] for j in range(1,5)],prec=5)
        def compose(p):
            value=arb_series([0],prec=5)
            for j in range(4,-1,-1):value=value*increment+p[j]
            return value
        E=flow.energy;M=flow.expansion.order
        n=[delta.sqrt()*sum((compose(ps[j][k])/E**j for j in range(1,M+1)),arb_series([0],prec=5))
           for k in (0,1)]
        n.append(sum((compose(zs[j])/E**j for j in range(M+1)),arb_series([0],prec=5)))
        dy=[arb_series([(j+1)*v[j+1] for j in range(4)],prec=4) for v in n]
        h=flow.hamiltonian_series(y,3);rotation=_cross(h,n)
        defect=flow.defect_series(y,3)
        for k in range(3):
            residual=dy[k]-2*rotation[k]
            for j in range(4):assert (residual[j]-defect[k][j]).contains(0)
        assert any(not (dy[k][0]-2*rotation[k][0]+defect[k][0]).contains(0) for k in (0,1))


class ConstantForcing:
    bits=160
    def hamiltonian_series(self,y,order):
        return tuple(arb_series([0],prec=order+1) for _ in range(3))
    def defect_series(self,y,order):
        return tuple(arb_series([v],prec=order+1) for v in (1,0,0))


def test_forcing_and_normalized_width_are_counted_once():
    row=np.zeros(26);row[1]=.25
    bound=correction_cell_defect(ConstantForcing(),row)
    assert bound==arb(1)/4
    assert bound*arb('.25')<arb(1)/4
    row[5]=-.25
    assert correction_cell_defect(ConstantForcing(),row)==0


def test_original_low_middle_covariance_error_is_resolved(pilot):
    flow,trace,proof,b,index,_=pilot
    assert proof['initial_error']>0
    assert proof['bloch_error']<arb('1e-13')
    assert np.linalg.norm(trace[-1,5:8])>1e-9
    sl=slice(3*index,3*index+3)
    distance=original_covariance_distance_bounds(b.initial_columns[:,sl],b.source.covariance[sl,sl],proof)
    assert distance['lower']>arb('8e-12')
    assert distance['upper']<arb('2e-11')
    assert proof['physical_local_gate']=='OPEN'


def test_join_and_initial_time_mutations_fail(pilot):
    flow,trace,_,b,_,_=pilot
    missing=np.delete(trace,len(trace)//2,axis=0)
    with pytest.raises(ValueError,match='join'):
        validate_correction(flow,missing,-18.,b.rho_up)
    with pytest.raises(ValueError,match='start'):
        validate_correction(flow,trace,-18.1,b.rho_up)
    bad=trace[0].copy();bad[1]=bad[0]
    with pytest.raises(ValueError):correction_cell_defect(flow,bad)


def test_domains_and_flint_context(pilot):
    flow,*_=pilot;before=ctx.prec,ctx.cap
    with pytest.raises(ValueError):VacuumCorrection(flow.expansion,None)
    with pytest.raises(ValueError):capture_correction(flow,0,-1)
    with pytest.raises(ValueError):correction_cell_defect(flow,np.zeros(26),degree=3)
    assert (ctx.prec,ctx.cap)==before


def test_replay_uses_the_saved_curve_and_rejects_source_mutation(pilot,monkeypatch):
    import runpy
    from copy import deepcopy
    from hashlib import sha256
    from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
    flow,trace,_,positive,index,nfev=pilot
    archive=RetainedUpstreamArchive(ROOT)
    negative=next(b for b,_ in archive.family_entries((14,1)) if b.energy_sign<0
        and b.original_panel==positive.original_panel and b.rows==positive.rows)
    sl=slice(3*index,3*index+3)
    arrays={'trace':trace,'positive_columns':positive.initial_columns[:,sl],
            'negative_columns':negative.initial_columns[:,sl],
            'positive_covariance':positive.source.covariance[sl,sl],
            'negative_covariance':negative.source.covariance[sl,sl]}
    raw=deterministic_npz_bytes(arrays)
    seed={'payload':{'sha256':sha256(raw).hexdigest()},'nfev':nfev}
    driver=runpy.run_path(str(ROOT/'scripts/derive_nsc_vacuum_source_correction.py'))
    def no_solve(*args,**kwargs):raise AssertionError('replay attempted another solve')
    monkeypatch.setitem(driver['calculate'].__globals__,'capture_correction',no_solve)
    result,replayed=driver['calculate']((seed,raw))
    assert replayed==raw
    assert result['initial_true_correction_assumed_zero'] is False
    assert result['source_quadrature_error_included'] is False
    assert result['physical_upstream_budget_component'] is None
    assert result['physical_local_gate']=='OPEN'
    changed=dict(arrays);changed['positive_columns']=arrays['positive_columns'].copy()
    changed['positive_columns'][0,1]+=.001
    altered=deterministic_npz_bytes(changed)
    forged=deepcopy(seed);forged['payload']['sha256']=sha256(altered).hexdigest()
    with pytest.raises(ValueError,match='original source'):
        driver['calculate']((forged,altered))
    with pytest.raises(ValueError,match='payload hash'):
        driver['calculate']((seed,altered))
