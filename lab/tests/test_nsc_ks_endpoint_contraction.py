from fractions import Fraction

import numpy as np
import pytest
from flint import arb,ctx

from recursive_horizons.nsc_ks_endpoint_contraction import (
    endpoint_difference_norms,contract_signed_pair,complete_endpoint_inputs)
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_finite_matter_error import source_norm_upper
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter


@pytest.mark.parametrize('mode',[2,-4])
def test_exact_harmonic_retains_carrier_and_nyquist(mode):
    # These harmonic samples have exact binary real/imaginary components.
    samples = np.array([1,1j,-1,-1j]*2) if mode==2 else np.array([1,-1]*4)
    d = np.stack((samples,2*samples))[:,None,:]
    result=endpoint_difference_norms(d,[.5],[.7],Fraction(1),bits=120)
    # Compare to an independent, tighter enclosure, not a separately rounded
    # upper at the same precision (which need not be the smaller upper).
    with ctx.workprec(240):
        exact=arb(5).sqrt()/2
        actual_z=abs(2*arb.pi()*mode-arb(.7))*exact
        assert restored_upper(result['difference_norm']) >= exact.upper()
        assert restored_upper(result['difference_axial_norm']) >= actual_z.upper()
        assert float(restored_upper(result['difference_norm'])) < float(exact)*1.000001
        assert float(restored_upper(result['difference_axial_norm'])) < float(actual_z)*1.000001
        assert actual_z > arb(.7)*exact


def test_norm_is_not_just_the_maximum_at_grid_nodes():
    # Two modes agree at a point between nodes. Samples underresolve that peak;
    # the continuous Fourier bound must still enclose it.
    x=np.arange(8)/8
    values=1+np.exp(2j*np.pi*(x-1/16))
    d=np.stack((values,np.zeros(8)))[:,None,:]
    result=endpoint_difference_norms(d,[1.],[0.],1)
    assert np.max(abs(values)) < 2
    assert float(restored_upper(result['difference_norm'])) >= 2-1e-14
    assert not result['nodal_maximum_used_as_bound']


def test_signed_contraction_keeps_distinct_covariances_and_covers_perturbations():
    rng=np.random.default_rng(42)
    C=np.array([[.8,.05j],[-.05j,.6]])
    Cneg=np.eye(2)-C.conj()
    arrays=[rng.normal(size=(3,2,2))+1j*rng.normal(size=(3,2,2)) for _ in range(8)]
    A,D,Az,Dz=arrays[:4]
    ea,ed,eaz,edz=[v*1e-7 for v in arrays[4:]]
    norm=lambda v:np.nextafter(np.linalg.norm(v,axis=(1,2)).max(),np.inf)
    kwargs=dict(reference_norm=norm(A),difference_norm=norm(D),
        reference_axial_norm=norm(Az),difference_axial_norm=norm(Dz),
        reference_error=norm(ea),difference_error=norm(ed),
        reference_axial_error=norm(eaz),difference_axial_error=norm(edz),
        source_norm=source_norm_upper(C),mass=.3,absolute_angular=.8,
        axial_lower=.8,radius_lower=1.4,multiplicity=2.)
    result=contract_signed_pair(kwargs,Cneg)
    empty=np.zeros((0,*A.shape),complex)
    def gradient(F,Z,cov,angular):
        return source_column_matter(F,Z,cov,empty,empty,mass=.3,angular=angular,
            axial_scale=.8,radius=1.4,multiplicity=2.)['action_gradient']
    total=np.zeros((3,2))
    for cov,sign in [(C,1),(Cneg,-1)]:
        transform=lambda v:v if sign==1 else np.einsum('ab,zbs->zas',np.diag([1,-1]),v.conj())
        f,z,f0,z0=map(transform,(A+D,Az+Dz,A,Az))
        fp,zp,ap,azp=map(transform,(A+D+ea+ed,Az+Dz+eaz+edz,A+ea,Az+eaz))
        total += gradient(fp,zp,cov,sign*.8)-gradient(ap,azp,cov,sign*.8)-gradient(f,z,cov,sign*.8)+gradient(f0,z0,cov,sign*.8)
    assert np.all(np.max(abs(total),axis=0)<[float(restored_upper(result[k])) for k in ('N','beta')])
    assert result['covered_energy_sectors']==[1,-1]
    assert result['positive']['N'] != result['negative']['N']
    assert result['physical_local_gate']=='OPEN'


def test_invalid_endpoint_weights_rejected():
    with pytest.raises(ValueError):endpoint_difference_norms(np.zeros((2,1,8)),[0.],[1.],1.)


def test_contraction_geometry_denominators_round_down():
    from recursive_horizons.nsc_ks_current_field_cone import whole_cone_continuous_inputs
    from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
    family=LocalIncomingFamily(np.zeros((2,8)))
    inputs=whole_cone_continuous_inputs(family,mass=.1,angular=.8,rho_up=1.03,
        source_energies=[.7],bits=80)
    zero={'mantissa':'0','exponent':0}
    propagation={k:zero for k in ('weighted_reference_F_error_upper',
        'weighted_reference_F_z_error_upper','weighted_difference_F_error_upper',
        'weighted_difference_F_z_error_upper')}
    joined=complete_endpoint_inputs(inputs,np.zeros((2,1)),np.zeros((2,1,8)),
        [1.],[.7],Fraction(1),propagation,bits=80)
    for name in ('axial_lower','radius_lower'):
        packed=joined['matter'][name]['value']
        value=Fraction(int(packed['mantissa']))*Fraction(2)**packed['exponent']
        assert value <= Fraction(inputs['geometry'][name]['exact_rational'])
    assert joined['physical_rho1_source_error'] is None
