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
from recursive_horizons.nsc_upstream_compensator import (
    bump,vertices,vertex_identities,fit_selected,combine,small_amplitude_range,PAULI,
    require_covariance_product)


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
    record=R.check(source_context='329ef809e466a86ed4a290ee32c19485f85add9d')
    assert len(record['control']['cases_completed'])<=6
    assert record['scope']['global_C0_matching']=='OPEN'
    assert record['selected_compensation']['same_coefficients_all_cases']
    assert not record['control']['energy_weights_used']
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:
        c=a['coefficients'].copy()
        for case in record['control']['cases_completed']:
            name=case.split('/')[0];f=a['input/'+name+'/stationary'];C=a['input/'+name+'/source'];K=combine(a[case+'/K'],c)
            cov=f@C@f.swapaxes(-1,-2).conj()
            require_covariance_product(cov,f,C)
            require_covariance_product(a['input/'+name+'/covariance'],f,C)
            for o,i in ((0,3),(3,3)):
                dAoi=K[o,i]@f[i];dAio=K[i,o]@f[o]
                full=dAoi@C[i]@f[i].conj().T+f[o]@C[o]@dAio.conj().T
                slab=K[o,i]@cov[i]+cov[o]@K[i,o].conj().T
                assert np.max(abs(full-slab))<3e-11
        for name in ('13_1','0_1'):
            if name+'/fine' in record['control']['cases_completed']:
                assert np.max(abs(a[name+'/fine/K'][0]))==0
    assert record['runtime']['CPU_seconds']<=60 or record['runtime']['cpu_budget_exceeded']


def test_fresh_signature_binds_retained_owners_and_history_is_explicit():
    from hashlib import sha256
    from recursive_horizons.nsc_frozen_sources import frozen_source_context
    current=R.signature()
    assert 'scripts/derive_nsc_retarded_response_resolution.py' not in current
    assert 'docs/nsc-upstream-compensator.md' not in current
    assert 'docs/nsc-preparation-domain-decision.md' in current
    assert 'docs/nsc-normalized-histories.md' in current
    assert 'docs/nsc-retarded-fixed-c0-certificate.md' in current
    for path,digest in current.items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==digest
    with frozen_source_context('329ef809e466a86ed4a290ee32c19485f85add9d'):
        original=R.signature()
    assert 'scripts/derive_nsc_retarded_response_resolution.py' in original
    assert 'docs/nsc-upstream-compensator.md' in original
    assert 'docs/nsc-retarded-response-resolution.md' in original


def covariance_operands():
    f=np.array([[[.1+.2j,.7-.3j,-.4+.6j],[-.2+.5j,.3+.1j,.8-.7j]]],dtype=np.complex128)
    source=np.array([[[1.,.3+.2j,-.1j],[.3-.2j,.7,.2+.4j],[.1j,.2-.4j,.9]]],dtype=np.complex128)
    return f,source,f@source@f.swapaxes(-1,-2).conj()


def test_covariance_checks_each_candidate_against_exact_operand_product():
    f,source,cov=covariance_operands()
    policy=require_covariance_product(cov,f,source)
    assert 0<=policy['maximum_exact_L1_residual']<=policy['maximum_operand_rounding_radius']
    assert policy['physical_error_bound'] is None
    bad=cov.copy();bad[0,0,0]+=1e-8
    # Equality of two identically corrupted products is not admission.
    assert np.array_equal(bad,bad.copy())
    for candidate in (bad,bad.copy()):
        with pytest.raises(ValueError,match='operand rounding envelope'):
            require_covariance_product(candidate,f,source)


def test_covariance_zero_mass_has_no_absolute_floor_and_underflow_is_scoped():
    f,source,cov=covariance_operands();source[:]=0;cov[:]=0
    assert require_covariance_product(cov,f,source)['maximum_operand_rounding_radius']==0
    bad=cov.copy();bad[0,0,0]=np.nextafter(0.,1.)
    with pytest.raises(ValueError,match='operand rounding envelope'):
        require_covariance_product(bad,f,source)
    # Nonzero operands can legitimately underflow; the explicit operation
    # count covers that case without loosening the exact-zero source case.
    f[:]=1e-200;source[:]=1e-200
    assert require_covariance_product(cov,f,source)['maximum_operand_rounding_radius']>0


@pytest.mark.parametrize('operand,mutation',[('covariance','dtype'),('stationary','layout'),
    ('source','nonfinite'),('covariance','list')])
def test_covariance_rejects_invalid_numeric_layouts(operand,mutation):
    f,source,cov=covariance_operands();values={'stationary':f,'source':source,'covariance':cov}
    if mutation=='dtype':values[operand]=values[operand].astype(np.complex64)
    elif mutation=='layout':values[operand]=values[operand][:,:,:2]
    elif mutation=='nonfinite':values[operand]=values[operand].copy();values[operand].flat[0]=np.nan
    else:values[operand]=values[operand].tolist()
    with pytest.raises(ValueError):require_covariance_product(**values)


def zero_family_inputs():
    data={'energies':np.array([-.7,-.2,.2,.7]),'mass':np.array(1.),'angular':np.array(0.),
          'source':np.zeros((4,3,3),complex),'projector':np.zeros((4,3,3),complex),
          'weights':np.ones(4),'stationary':np.zeros((4,2,3),complex),'covariance':np.zeros((4,2,2),complex)}
    arrays={'input/toy/'+key:value.copy() for key,value in data.items()}
    x=np.linspace(.8,1.2,801);phi=np.zeros((1602,12),complex)
    return arrays,data,x,phi


@pytest.mark.parametrize('mutation',['value','dtype','signed_zero'])
def test_raw_source_is_bit_exact_before_any_arithmetic_policy(mutation):
    arrays,data,x,phi=zero_family_inputs()
    if mutation=='value':data['source'][0,0,0]=np.nextafter(0.,1.)
    elif mutation=='dtype':data['source']=data['source'].astype(np.complex64)
    else:data['source'][0,0,0]=complex(-0.,0.)
    with pytest.raises(ValueError,match='saved original source changed: toy/source'):
        R.require_family_inputs(arrays,'toy',data,x,phi)


def test_stationary_map_and_covariance_cannot_self_enlarge_their_allowance():
    arrays,data,x,phi=zero_family_inputs()
    R.require_family_inputs(arrays,'toy',data,x,phi)
    arrays['input/toy/stationary'][0,0,0]=1e-8
    with pytest.raises(ValueError,match='archived.*map'):
        R.require_family_inputs(arrays,'toy',data,x,phi)


def test_co_corrupted_stationary_and_matching_covariance_fail_sealed_history_anchor():
    arrays,data,x,phi=zero_family_inputs()
    data['source'][:]=np.eye(3);arrays['input/toy/source']=data['source'].copy()
    data['stationary'][:]=1e-8
    data['covariance']=data['stationary']@data['source']@data['stationary'].swapaxes(-1,-2).conj()
    arrays['input/toy/stationary']=data['stationary'].copy()
    arrays['input/toy/covariance']=data['covariance'].copy()
    # The two maps agree, and both covariance products are arithmetically
    # correct, but neither map comes from the authenticated zero history.
    require_covariance_product(data['covariance'],data['stationary'],data['source'])
    with pytest.raises(ValueError,match='archived.*map'):
        R.require_family_inputs(arrays,'toy',data,x,phi)
