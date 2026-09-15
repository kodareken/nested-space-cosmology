"""Boundary orientation on the real shared surface, not two time endpoints."""
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_smooth_seam_variation import smooth_seam_match,reverse_normal
from derive_nsc_smooth_seam_variation import calculate,reference_jet


def test_actual_reference_seam_is_bound_to_the_same_local_density():
    result=calculate();matching=result['matching']
    assert max(matching['residuals'].values())<3e-11
    assert result['reference_metric_normal_jet']['r']==1.
    assert abs(matching['parent']['Weyl_scalar_F']-2)<3e-11
    assert abs(matching['parent']['weyl'][2]-.10900236876297109)<3e-11
    assert result['legacy_diagnostic']['Weyl_value']==93.54264532195464
    assert result['quantum_account']['Gamma_rest_assigned_zero'] is False
    assert matching['normal_vectors_PG']['parent'][0]>1
    assert np.max(abs(np.asarray(matching['normal_vectors_PG']['parent'])+
                      np.asarray(matching['normal_vectors_PG']['child'])))<3e-11


def test_shear_variation_reverses_but_its_conjugate_does_not():
    r=calculate();p=r['matching']['parent'];c=r['matching']['child']
    assert np.max(abs(np.asarray(p['weyl'])[:2]+np.asarray(c['weyl'])[:2]))<3e-11
    assert abs(p['weyl'][2]-c['weyl'][2])<3e-11
    assert abs(p['sigma']+c['sigma'])<3e-11
    assert abs(p['weyl'][2]+c['weyl'][2])>.2
    # Opposite normals are not an instruction to add a second minus to every
    # entry. The tangent-space pullback is part of the boundary comparison.
    assert np.max(abs(np.asarray(p['weyl'])+np.diag([1.,1.,-1.])@np.asarray(c['weyl'])))<3e-11


def test_smooth_gluing_identity_uses_the_common_jet_not_a_fitted_force():
    r=calculate();jet,_=reference_jet();jet=dict(jet)
    jet.update(N=1.1,Nd=.07,Ndd=-.02,ad=.6,add=1.2,addd=-.4,rd=.1,rdd=.5,rddd=.3)
    matched=smooth_seam_match(jet,r['locked_inputs'])
    assert max(matched['residuals'].values())<3e-11
    assert reverse_normal(reverse_normal(jet))==jet
    assert matched['quantum_boundary_remainder'] is None
    assert matched['bulk_metric_residual'] is None
