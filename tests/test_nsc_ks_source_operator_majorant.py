"""Continuous operator bounds checked against an independent matrix evolution."""
import numpy as np
import pytest
from flint import arb,ctx
from scipy.linalg import expm

from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_source_operator_majorant import (
    fundamental_difference_bounds,source_moment_insertion,history_source_insertion,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_transmitting_dirac_domain import S1,S2,S3


def test_characteristic_bounds_cover_independent_spatial_matrix_exponential():
    count=64;z=2*np.pi*np.arange(count)/count
    wave=np.fft.fftfreq(count,1/count)
    derivative=np.fft.ifft(1j*wave[:,None]*np.fft.fft(np.eye(count),axis=0),axis=0)
    mass,ell,energy,axial,radius,amplitude,duration=.3,1.,.2,1.,2.,.02,.1
    actual_radius=radius+amplitude*np.cos(z)
    generator=np.kron(S3,derivative)/axial**2-1j*energy/axial**2*np.kron(S3,np.eye(count))
    for i in range(2):
        for j in range(2):
            generator[i*count:(i+1)*count,j*count:(j+1)*count]+=np.diag(
                1j/axial*(-mass*S1[i,j]+ell/actual_radius*S2[i,j]))
    initial=np.zeros((2*count,2),complex);initial[:count,0]=1;initial[count:,1]=1
    changed=(expm(duration*generator)@initial).reshape(2,count,2).transpose(1,0,2)
    ref=expm(1j*duration/axial*(-mass*S1+ell/radius*S2-energy/axial*S3))
    reference=np.broadcast_to(ref,changed.shape)
    envelope_z=np.einsum('ij,jab->iab',derivative,changed)
    difference=changed-reference
    physical_z=envelope_z-1j*energy*changed
    reference_z=-1j*energy*reference
    with ctx.workprec(192):
        minimum=arb(radius)-arb(amplitude)
        B=arb(duration)/arb(axial)*(arb(mass)**2+(arb(ell)/minimum)**2).sqrt()
        M=arb(duration)*arb(ell)*arb(amplitude)/(arb(axial)*arb(radius)*minimum)
        Mz=arb(duration)*arb(ell)*arb(amplitude)/(arb(axial)*minimum**2)
        norms=fundamental_difference_bounds(B,M,Mz)
        assert np.max(np.linalg.norm(difference,axis=(1,2)))<float(norms['difference'])
        assert np.max(np.linalg.norm(envelope_z,axis=(1,2)))<float(norms['difference_envelope_z'])
        measured=np.max(np.linalg.norm(physical_z-reference_z,axis=(1,2)))
        assert measured<float(energy*norms['difference']+norms['difference_envelope_z'])
        assert measured>float(energy*norms['difference'])  # dropping Mz fails
        epsilon=.01;C=.4*np.eye(2);delta=epsilon*(S1+S2)/np.sqrt(2)
        zeros=np.zeros((0,*changed.shape),complex)
        def value(F,Fz,cov):
            return source_column_matter(F,Fz,cov,zeros,zeros,mass=mass,angular=ell,
                axial_scale=axial,radius=radius,multiplicity=1)['action_gradient']
        observed=(value(changed,physical_z,C+delta)-value(reference,reference_z,C+delta)
                  -value(changed,physical_z,C)+value(reference,reference_z,C))
        bound=source_moment_insertion(norms,epsilon,energy*epsilon,mass=mass,
            absolute_angular=ell,axial_lower=axial,radius_lower=radius,multiplicity=1)
        for i,name in enumerate(('N','beta')):
            assert np.max(abs(observed[:,i]))>1e-10
            assert np.max(abs(observed[:,i]))<float(restored_upper(bound[name]))


def test_zero_history_has_zero_source_difference_but_not_zero_baseline_error():
    zero=LocalIncomingFamily(np.zeros((2,8)))
    result=history_source_insertion(zero,1.03,1e-3,1e-2,mass=1,absolute_angular=2,multiplicity=12)
    assert restored_upper(result['N'])==0
    assert restored_upper(result['beta'])==0
    assert result['physical_upstream_budget_component'] is None
    assert result['physical_local_gate']=='OPEN'


def test_actual_history_bound_is_continuous_and_rejects_wrong_preparation_support():
    family,_=require_original_history()
    result=history_source_insertion(family,1.0300000000000002,1e-14,3e-13,
        mass=np.pi/2,absolute_angular=np.sqrt(5),multiplicity=12)
    assert result['whole_spatial_line'] and result['continuous_on_incoming_interval']
    assert not result['numerical_field_norms_used']
    assert restored_upper(result['N'])<arb('1e-14')
    with pytest.raises(ValueError,match='outside the history'):
        history_source_insertion(family,1.001,1e-14,3e-13,mass=1,absolute_angular=2,multiplicity=12)


def test_missing_bounds_are_not_zero_and_context_is_restored():
    before=ctx.prec,ctx.cap
    with pytest.raises(ValueError):fundamental_difference_bounds(1,None,1)
    with pytest.raises(ValueError):fundamental_difference_bounds(arb(0,'.1'),1,1)
    norms=fundamental_difference_bounds(.1,.01,.1)
    with pytest.raises(ValueError):
        source_moment_insertion(norms,None,1,mass=1,absolute_angular=2,
            axial_lower=1,radius_lower=1,multiplicity=1)
    assert (ctx.prec,ctx.cap)==before


def test_real_source_contribution_and_binding_mutations():
    import runpy
    from pathlib import Path
    from copy import deepcopy
    from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
    root=Path(__file__).resolve().parents[1]
    driver=runpy.run_path(str(root/'scripts/derive_nsc_ks_source_operator_majorant.py'))
    result=driver['calculate']()
    assert result['source_coverage']['covered_signed_rows']==772
    assert result['source_coverage']['remaining_signed_rows_in_family']==396
    action=result['continuous_action_error_upper']
    assert restored_upper(action['N'])<arb('5e-14')
    assert restored_upper(action['beta'])<arb('4e-14')
    assert action['continuous_on_incoming_interval']
    assert result['physical_upstream_budget_component'] is None
    records=[driver['read_bound'](driver[k]) for k in ('HIGH','LOW','MIDDLE')]
    rows=driver['selected_bounds'](*records)
    archive=RetainedUpstreamArchive(root);entries=archive.family_entries((14,1))
    with ctx.workprec(192):
        with pytest.raises(ValueError,match='duplicate'):
            driver['source_error_moments'](archive,entries,[rows[0],rows[0]])
        wrong=deepcopy(rows[0]);wrong['epsilon']=None
        with pytest.raises(ValueError,match='missing covariance'):
            driver['source_error_moments'](archive,entries,[wrong])
        wrong=deepcopy(rows[0]);wrong['energy_hex']=float(33).hex()
        with pytest.raises(ValueError,match='preparation'):
            driver['source_error_moments'](archive,entries,[wrong])
        wrong=deepcopy(rows[0]);wrong['preparation_digest']='0'*64
        with pytest.raises(ValueError,match='preparation'):
            driver['source_error_moments'](archive,entries,[wrong])
    wrong=deepcopy(records);wrong[0]['archived_comparison']['unweighted_covariance_errors']=False
    with pytest.raises(ValueError,match='unweighted'):
        driver['selected_bounds'](*wrong)
