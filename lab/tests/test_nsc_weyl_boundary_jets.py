"""Same-action boundary-jet identities and the frozen nodal coefficients."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_weyl_boundary_jets import nodal_decomposition
from derive_nsc_weyl_boundary_jets import calculate


def test_all_frozen_nodal_coefficients_are_reconstructed():
    r=calculate()
    assert r['preserved_Weyl_value']==93.54264532195464
    assert r['residuals']['legacy_eight_coefficient_binding']==0
    assert r['residuals']['full_gradient_vs_locked']<3e-11
    assert r['residuals']['gradient_decomposition']<3e-11
    assert r['residuals']['stencil_SBP_defect_norm']>.1
    assert r['gate']['physical_two_sided_mismatch'] is None


def test_normal_shear_is_the_independent_Weyl_boundary_rate():
    r=calculate();proper=r['proper_boundary_coefficients']
    assert max(abs(v) for v in proper['proper_boundary_N'])<3e-11
    assert min(abs(v) for v in proper['proper_boundary_shear'])>.5
    assert r['residuals']['Weyl_normal_trace']<3e-11
    assert r['remaining_action_ownership']['Gamma_rest_assigned_zero'] is False


def test_fixed_values_do_not_fix_normal_derivative_variations():
    r=json.loads((ROOT/'results/development/nsc-landau-cauchy-isometry.json').read_text())
    with np.load(ROOT/r['payload']['path'],allow_pickle=False) as a:
        i=np.linspace(0,len(a['short_time'])-1,33,dtype=int)
        h=SimpleNamespace(time=a['short_time'][i],lapse=a['short_lapse'][i],
            a_parallel=a['short_a'][i],radius=a['short_r'][i])
    coefficient=json.loads((ROOT/'results/development/nsc-general-ks-local-history.json').read_text())['locked_inputs']['C_Weyl']
    d=nodal_decomposition(h,coefficient);D=d['derivative_matrix'];n=len(h.time)
    s=(h.time-h.time[0])/(h.time[-1]-h.time[0]);direction=s*(1-s)
    assert direction[0]==direction[-1]==0
    assert abs(d['canonical_boundary_component'][3]@direction)>.1
    # This projection restricts a test variation, not a quantum state.
    A=np.vstack((np.eye(n)[0],np.eye(n)[-1],D[0],D[-1]))
    fixed=direction-A.T@np.linalg.solve(A@A.T,A@direction)
    assert np.max(abs(A@fixed))<3e-11
    assert abs(d['canonical_boundary_component'][3]@fixed)<3e-11
