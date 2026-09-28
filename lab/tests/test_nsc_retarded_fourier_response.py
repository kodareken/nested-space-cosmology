"""Independent Gaussian-transform normalization and saved Fourier receipt."""
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
SPEC=importlib.util.spec_from_file_location('fourier_response',ROOT/'scripts/derive_nsc_retarded_fourier_response.py')
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)


def test_analytic_gaussian_phase_modulation_preserves_CAR_and_normalization():
    # Exact independent transform: integral exp(-z²)exp(i deltaE z) dz.
    E=np.array([-.3,.4]);shift=.7;z=np.linspace(-8.,8.,2049);t=z-shift
    trace=np.zeros((len(z),2,4),complex)
    for i,energy in enumerate(E):
        trace[:,:,2*i:2*i+2]=1j*np.exp(-z*z-1j*energy*z)[:,None,None]*np.eye(2)
    A=M.transform_modes(E,t,trace,shift)
    exact=1j*np.sqrt(np.pi)*np.exp(-(E[:,None]-E[None,:])**2/4)
    assert np.max(abs(A-exact[:,:,None,None]*np.eye(2)))<3e-12
    f=np.broadcast_to(np.eye(2),(2,2,2));C=np.array([np.diag([.2,.8]),np.diag([.7,.3])])
    expected=np.array([[exact[o,i]*(C[i]-C[o]) for i in range(2)] for o in range(2)])
    assert np.max(abs(M.covariance_pairs(A,f,C)-expected))<3e-12
    assert np.max(abs(M.covariance_pairs(A,f,f)))<3e-12
    # Tangents are not covariance matrices: the nonzero pair is not clipped.
    assert np.max(abs(expected))>.5


def test_invalid_time_or_source_shapes_are_rejected():
    E=np.array([.2]);trace=np.zeros((3,2,3))
    with pytest.raises(ValueError,match='ordered'):
        M.transform_modes(E,[0.,.2,.1],trace,0.)
    with pytest.raises(ValueError,match='source'):
        M.transform_modes(E,[0.,.1,.2],trace[:2],0.)
    with pytest.raises(ValueError,match='matching'):
        M.covariance_pairs(np.zeros((1,1,2,3)),np.zeros((1,2,3)),np.eye(2)[None])


def test_receipt_retains_full_pair_scope_and_authenticated_inputs():
    record=json.loads(M.OUTPUT.read_text())
    assert record['source_sha256']==M.signature()
    assert not record['convention']['energy_quadrature_weights_used']
    assert not record['convention']['numerical_trace_clipped']
    assert record['scope']['field_or_radial_solves']==0
    assert record['scope']['rigorous_continuum_response_error_bound'] is None
    assert record['scope']['exact_fixed_C0_exclusion_certificate']=='OPEN'
    assert not record['scope']['CAR_residual_subtracted']
    fine=record['per_case']['3201_256']
    assert fine['closed_input_tangent_column']==0
    assert fine['phase_origin_residual']<record['indicators']['algebra_tolerance']
