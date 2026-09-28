"""Independent finite-matrix/ODE controls; no authenticated field runs."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from scipy.sparse import bmat,csr_matrix,diags

from test_nsc_exact_phase_prepared_state import setup
from recursive_horizons.nsc_cf4_prepared_state import prepare_cf4_incoming
from recursive_horizons.nsc_evolved_incoming_state import AmplitudeOnlyMetric,CachedCompatibleIncomingMetric
from recursive_horizons.nsc_transmitting_history_jets import incoming_reference_columns


def test_constant_case_is_exact_augmented_flow_and_pde():
    owner,provider,family,times,phi,source=setup(amplitude=0.,ndir=0)
    result=prepare_cf4_incoming(owner,times,provider,phi,source)
    L,_=owner.generator(owner.reference);B=incoming_reference_columns(owner,phi)
    A=bmat([[L,csr_matrix(B)],[None,diags(-1j*source.energies)]]).toarray()
    initial=np.vstack((phi,np.eye(len(source.energies))))
    ref=np.array([expm(t*A)@initial for t in times])
    np.testing.assert_allclose(result.diagnostics['node_fields'],ref[:,:len(phi)],atol=3e-12,rtol=0)
    assert result.diagnostics['exact_harmonic_phase_residual']<3e-12
    result.require_history(provider,times,owner.x)
    altered=result.gauss_metrics.copy();altered[0,3,0]+=.01
    with pytest.raises(ValueError,match='Gauss profile'):
        replace(result,gauss_metrics=altered)


def test_nonzero_cf4_retarded_pde_tangent_matches_whole_family_difference():
    owner,provider,family,_,phi,source=setup()
    times=np.linspace(0.,.04,17)
    result=prepare_cf4_incoming(owner,times,provider,phi,source)
    h=1e-4
    op,pp,_,_,ip,sp=setup(amplitude=.003+h,ndir=0)
    om,pm,_,_,im,sm=setup(amplitude=.003-h,ndir=0)
    plus=prepare_cf4_incoming(op,times,pp,ip,sp)
    minus=prepare_cf4_incoming(om,times,pm,im,sm)
    np.testing.assert_allclose((plus.state.columns-minus.state.columns)/(2*h),result.state.column_tangents[0],atol=3e-8,rtol=0)
    np.testing.assert_allclose((plus.axial_columns-minus.axial_columns)/(2*h),result.axial_tangents[0],atol=3e-8,rtol=0)
    assert result.diagnostics['CF4_time_integrator'] is True
    assert result.fixed_preparation_digest==plus.fixed_preparation_digest==minus.fixed_preparation_digest
    assert result.diagnostics['metric_flux_tangent_residual']<3e-11


def test_cf4_converges_against_independent_adaptive_ode():
    owner,_,family,_,phi,source=setup(amplitude=.02,ndir=0)
    cache=CachedCompatibleIncomingMetric(owner.x,family)
    B=incoming_reference_columns(owner,phi)
    def rhs(t,y):
        L,_=owner.generator(cache.values(t,owner.x))
        return (L@y.reshape(phi.shape)+B*np.exp(-1j*source.energies*t)).ravel()
    reference=solve_ivp(rhs,(0.,.04),phi.ravel(),method='DOP853',rtol=2e-12,atol=2e-14,max_step=.0002)
    assert reference.success
    errors=[]
    for steps in (24,48,96):
        result=prepare_cf4_incoming(owner,np.linspace(0.,.04,steps+1),AmplitudeOnlyMetric(cache),phi,source)
        errors.append(float(np.max(abs(result.diagnostics['node_fields'][-1].ravel()-reference.y[:,-1]))))
    assert errors[2]<errors[1]/8, errors
    assert errors[1]<errors[0]/8, errors
