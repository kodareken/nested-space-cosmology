"""Common-geometry algebra and authenticated six-case replay; no new solves."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
SPEC=importlib.util.spec_from_file_location('upstream_compensator',ROOT/'scripts/derive_nsc_upstream_compensator.py')
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)
from recursive_horizons.nsc_upstream_compensator import bump,vertices,vertex_identities,fit_selected,combine,small_amplitude_range,PAULI


def test_bump_and_Weyl_vertices_are_actual_common_metric_directions():
    assert all(bump(rho)==0 for rho in (.9,1.,1.01,1.03))
    assert bump(1.02)==pytest.approx(1.)
    proof=vertex_identities()
    assert all(v=='0' for v in proof['raw_vertex_residuals']) and proof['Weyl_midpoint_residual']=='0'
    E=np.array([-.7,-.2,.2,.7]);v=vertices(1.02,E,1.3,2.)
    assert np.max(abs(v-v.transpose(0,2,1,4,3).conj()))==0
    # Opposite frequencies have zero Ebar kinetic insertion, rather than Ei.
    assert np.max(abs(v[1,0,3]))==0


def test_one_direct_real_fit_and_positive_chart_range():
    K=np.zeros((4,4,4,2,2),complex);target=np.array([.2,-.1,.3])
    K[0,3,3]=-1j*np.einsum('j,jab->ab',target,PAULI)
    K[1:,3,3]=-1j*PAULI
    c,checks=fit_selected(K)
    assert np.max(abs(c+target))<1e-15 and np.max(abs(combine(K,c)))<1e-15
    valid=small_amplitude_range(c);eta=Q(valid['strict_abs_eta_upper']['numerator'],valid['strict_abs_eta_upper']['denominator'])
    assert eta>0 and eta*Q(abs(float(c[0])))<Q(1,100)
    assert eta*Q(abs(float(c[1])))<Q(1,125)
    assert eta*(Q(abs(float(c[2])))+Q(3,100))<Q(1,2)
    bad=K.copy();bad[0,3,3]+=np.eye(2)
    with pytest.raises(ArithmeticError,match='anti-Hermitian'):fit_selected(bad)


def test_saved_replay_retains_all_source_fibers_and_fixed_coefficients(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('old or new scientific solve called during replay')
    monkeypatch.setattr(R,'operator_response',forbidden)
    monkeypatch.setattr(R.D,'short_response',forbidden)
    monkeypatch.setattr(R.P,'evolve_prepared_field_jets',forbidden)
    record=R.check()
    assert len(record['control']['cases_completed'])<=6
    assert record['scope']['global_C0_matching']=='OPEN'
    assert record['selected_compensation']['same_coefficients_all_cases']
    assert not record['control']['energy_weights_used']
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:
        c=a['coefficients'].copy()
        for case in record['control']['cases_completed']:
            name=case.split('/')[0];f=a['input/'+name+'/stationary'];C=a['input/'+name+'/source'];K=combine(a[case+'/K'],c)
            cov=f@C@f.swapaxes(-1,-2).conj();assert np.array_equal(cov,a['input/'+name+'/covariance'])
            for o,i in ((0,3),(3,3)):
                dAoi=K[o,i]@f[i];dAio=K[i,o]@f[o]
                full=dAoi@C[i]@f[i].conj().T+f[o]@C[o]@dAio.conj().T
                slab=K[o,i]@cov[i]+cov[o]@K[i,o].conj().T
                assert np.max(abs(full-slab))<3e-11
        for name in ('13_1','0_1'):
            if name+'/fine' in record['control']['cases_completed']:
                assert np.max(abs(a[name+'/fine/K'][0]))==0
    assert record['runtime']['CPU_seconds']<=60 or record['runtime']['cpu_budget_exceeded']
