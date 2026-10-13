"""Covariance error transport in absolute Fourier coefficients and k moments.

The Fourier-diagonal Hermitian Dirac evolution is treated exactly as a
left/right unitary on each symbol coefficient. Only the spatially varying
potential enters the growth bound. No physical state or action is replaced.
"""
from flint import arb,ctx
from .nsc_bloch_trace_bound import _norm_upper,_positive_interval,local_stress_from_traces
from .nsc_ks_ball_trajectory import exact_upper,restored_upper


def propagate_covariance_moments(initial_zero,initial_first,residual_zero,residual_first,
                                 potential_zero_integral,potential_first_integral,*,bits=160):
    """M0 <= e^(2K0) A0; M1 <= e^(2K0)(A1+K1 A0).

    Mj=sum_l integral |k|^j ||D_l(k)||_1 dk/(2pi), j=0,1.
    Kj=integral sum_l |omega_l|^j ||delta V_l||_op |d rho|.
    A_j=initial moment + time-integrated true-defect moment. None is unknown.
    """
    with ctx.workprec(bits):
        i0,i1,r0,r1,v0,v1=[_norm_upper(v,n,bits=bits) for v,n in (
            (initial_zero,'initial D_up zeroth moment'),(initial_first,'initial D_up first moment'),
            (residual_zero,'time-integrated R zeroth moment'),(residual_first,'time-integrated R first moment'),
            (potential_zero_integral,'potential zeroth integral'),(potential_first_integral,'potential first integral'))]
        a0,a1=i0+r0,i1+r1;growth=(2*v0).exp()
        m0=growth*a0;m1=growth*(a1+v1*a0)
        return {'zeroth_moment_upper':exact_upper(m0),'first_moment_upper':exact_upper(m1),
            'growth_upper':exact_upper(growth),'moment_shift_coefficient':exact_upper(v1),
            'matrix_norm':'nuclear per2x2 Fourier coefficient; potential operator norm',
            'momentum_weights':('1','abs(k)'),'momentum_measure':'dk/(2pi)',
            'half_transfer_factor_retained':True,'Fourier_diagonal_evolution_unitary':True,
            'initial_D_up_retained':True,'time_integrated_R_retained':True,
            'physical_local_gate':'OPEN','source_law_changed':False,'new_action_term':False}


def pure_radius_potential_integrals(q_fourier_zero,q_fourier_first,*,absolute_angular,
                                    axial_lower,rho_up,rho_sigma=1,bits=160):
    """Uniform q=1/r_g-1/r_ref Fourier moments -> true potential integral bounds."""
    with ctx.workprec(bits):
        q0=_norm_upper(q_fourier_zero,'reciprocal zeroth Fourier moment',bits=bits)
        q1=_norm_upper(q_fourier_first,'reciprocal first Fourier moment',bits=bits)
        ell=_norm_upper(absolute_angular,'absolute angular',bits=bits)
        a=_positive_interval(axial_lower,'axial lower',bits=bits)
        up=_positive_interval(rho_up,'upstream coordinate',bits=bits)
        end=_positive_interval(rho_sigma,'incoming coordinate',bits=bits)
        if not up>=end:raise ValueError('decreasing preparation coordinate required')
        factor=(up-end)*ell/a
        return {'potential_zero_integral_upper':exact_upper(factor*q0),
            'potential_first_integral_upper':exact_upper(factor*q1),
            'homogeneous_mass_counted_as_growth':False,'potential_moment_measure':'sum of cell Fourier coefficients',
            'clock':'absolute d rho; H_rho potential=(-m S1+ell/r S2)/a',
            'physical_local_gate':'OPEN'}


def covariance_matter_error(moment_record,*,mass,absolute_angular,axial_lower,radius_lower,
                            multiplicity,continuum_embedding_established,bits=160):
    """Uniform coincidence/current bounds from M0/M1; no additional k powers."""
    if continuum_embedding_established is not True:
        raise ValueError('continuum local embedding must be established separately')
    with ctx.workprec(bits):
        stress=local_stress_from_traces(
            moment_record['zeroth_moment_upper'],moment_record['first_moment_upper'],mass=mass,
            absolute_angular=absolute_angular,axial_lower=axial_lower,radius_lower=radius_lower,
            multiplicity=multiplicity,bits=bits)
        return {**stress,'uniform_on_local_interval':True,'physical_local_gate':'OPEN',
            'covariance_coincidence_uses_zeroth_moment':True,'current_uses_first_moment':True,
            'numerical_operator_error_included':False,'source_accuracy_included':False}

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: radius_control
# Folded from scripts/derive_nsc_ks_radius_covariance_bounds.py at 329ef809; original receipts retain their producer identities.
import argparse as _radius_control_argparse
from fractions import Fraction as _radius_control_Q
from hashlib import sha256 as _radius_control_sha256
import json as _radius_control_json
from pathlib import Path as _radius_control_Path
import sys as _radius_control_sys
import numpy as _radius_control_np
from flint import arb as _radius_control_arb, ctx as _radius_control_ctx
from recursive_horizons.nsc_posix import ROOT as _radius_control_ROOT
_radius_control_sys.path.insert(0, str(_radius_control_ROOT / 'src'))
from recursive_horizons.nsc_ks_high_radius_derivatives import bound_pure_radius_difference as _radius_control_bound_pure_radius_difference, restore_packed as _radius_control_restore_packed
from recursive_horizons.nsc_fourier_covariance_moments import pure_radius_potential_integrals as _radius_control_pure_radius_potential_integrals
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper as _radius_control_exact_upper, restored_upper as _radius_control_restored_upper
_radius_control_OUTPUT = _radius_control_ROOT / 'results/development/nsc-ks-radius-covariance-bounds.json'
_radius_control_OWNED = ('scripts/derive_nsc_ks_radius_covariance_bounds.py', 'docs/nsc-ks-radius-covariance-bounds.md', 'src/recursive_horizons/nsc_ks_high_radius_derivatives.py', 'src/recursive_horizons/nsc_fourier_covariance_moments.py', 'src/recursive_horizons/nsc_bloch_trace_bound.py', 'src/recursive_horizons/nsc_ks_profile_fourier_bound.py')
_radius_control_RECORDS = ('results/development/nsc-ks-fine-profile-table.json', 'results/development/nsc-ks-radius-bound.json', 'results/development/nsc-ks-source-inventory.json', 'results/development/nsc-ks-fine-trajectory.json')
_radius_control_BITS = 160

def _radius_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_radius_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_radius_control_ROOT, relative)

def _radius_control_require_normal_support(geometry, history):
    sigma = _radius_control_Q(geometry['normal_support']['exact_rational'])
    directions = history['history_description']['directions']
    if not directions:
        raise ValueError('explicit normal windows required')
    for direction in directions:
        inner, outer = (_radius_control_Q(float.fromhex(direction[k])) for k in ('normal_inner', 'normal_outer'))
        if not 0 < inner < outer <= sigma:
            raise ValueError('history exceeds the certified normal support')

def _radius_control_compute():
    records = [_radius_control_json.loads((_radius_control_ROOT / p).read_text()) for p in _radius_control_RECORDS]
    for record in records:
        for category in ('source_hashes', 'input_hashes'):
            for p, h in record.get(category, {}).items():
                if _radius_control_digest(p) != h:
                    raise ValueError('input dependency changed: ' + p)
        if 'payload' in record:
            if _radius_control_digest(record['payload']['path']) != record['payload']['sha256']:
                raise ValueError('input payload changed')
    with _radius_control_np.load(_radius_control_ROOT / records[2]['payload']['path'], allow_pickle=False) as f:
        source = _radius_control_json.loads(f['metadata_json'].tobytes())
    with _radius_control_np.load(_radius_control_ROOT / records[3]['payload']['path'], allow_pickle=False) as f:
        history = _radius_control_json.loads(f['metadata_json'].tobytes())
    with _radius_control_ctx.workprec(_radius_control_BITS):
        geometry = _radius_control_bound_pure_radius_difference(root=_radius_control_ROOT, max_derivative=9, bits=_radius_control_BITS)
        if geometry['axial_profile_identity'] != history['axial_profile_identity']:
            raise ValueError('radius and evolved profile identities differ')
        _radius_control_require_normal_support(geometry, history)
        if history['rho_up'] > float(_radius_control_Q(records[1]['bounds']['rho_upper']['exact_rational'])):
            raise ValueError('history outside geometry slab')
        q0, q1 = (_radius_control_restore_packed(geometry['q_difference_fourier_l1'][j]) for j in (0, 1))
        rows = []
        growth_max = _radius_control_arb(0)
        shift_max = _radius_control_arb(0)
        for channel in source['channels']:
            potential = _radius_control_pure_radius_potential_integrals(q0, q1, absolute_angular=channel['angular_eigenvalue'], axial_lower=_radius_control_Q(records[1]['bounds']['axial_lower']['exact_rational']), rho_up=history['rho_up'], bits=_radius_control_BITS)
            K0 = _radius_control_restored_upper(potential['potential_zero_integral_upper'])
            K1 = _radius_control_restored_upper(potential['potential_first_integral_upper'])
            growth = (2 * K0).exp()
            shift = growth * K1
            growth_max = _radius_control_arb.max(growth_max, growth)
            shift_max = _radius_control_arb.max(shift_max, shift)
            rows.append({'group': channel['index'], 'angular_abs': channel['angular_eigenvalue'], 'compact_mass': channel['compact_mass'], 'actual_angular_signs': 1 if channel['angular_eigenvalue'] == 0 else 2, 'potential_integrals': potential, 'moment_comparison_diagonal_upper': _radius_control_exact_upper(growth), 'zeroth_to_first_moment_upper': _radius_control_exact_upper(shift), 'initial_moments': None, 'true_defect_moments': None})
        return {'schema': 'NSC-KS-RADIUS-COVARIANCE-BOUNDS-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: nonlinear spatial geometry and covariance-growth components; UV/local gate OPEN', 'geometry': geometry, 'groups': rows, 'actual_angular_families': sum((r['actual_angular_signs'] for r in rows)), 'maximum_moment_diagonal_upper': _radius_control_exact_upper(growth_max), 'maximum_moment_off_diagonal_upper': _radius_control_exact_upper(shift_max), 'history_identity': history['history_identity'], 'source_inventory_unchanged': True, 'source_hashes': {p: _radius_control_digest(p) for p in _radius_control_OWNED}, 'input_hashes': {p: _radius_control_digest(p) for p in _radius_control_RECORDS}, 'scope': {'physical_local_gate': 'OPEN', 'UV_error_bound': None, 'source_error_bound': None, 'actual_projector_derivative_bound': None, 'time_integrated_defect_bound': None, 'initial_C_up_minus_P_up_moments': None, 'new_state_or_action': False, 'new_field_or_source_runs': 0, 'metric_timestep': False}, 'reproducer': 'python scripts/derive_nsc_ks_radius_covariance_bounds.py --check'}

from types import SimpleNamespace as _radius_control_Namespace
radius_control = _radius_control_Namespace(
    BITS=_radius_control_BITS,
    OUTPUT=_radius_control_OUTPUT,
    OWNED=_radius_control_OWNED,
    Path=_radius_control_Path,
    Q=_radius_control_Q,
    RECORDS=_radius_control_RECORDS,
    ROOT=_radius_control_ROOT,
    arb=_radius_control_arb,
    argparse=_radius_control_argparse,
    bound_pure_radius_difference=_radius_control_bound_pure_radius_difference,
    compute=_radius_control_compute,
    ctx=_radius_control_ctx,
    digest=_radius_control_digest,
    exact_upper=_radius_control_exact_upper,
    json=_radius_control_json,
    np=_radius_control_np,
    pure_radius_potential_integrals=_radius_control_pure_radius_potential_integrals,
    require_normal_support=_radius_control_require_normal_support,
    restore_packed=_radius_control_restore_packed,
    restored_upper=_radius_control_restored_upper,
    sha256=_radius_control_sha256,
    sys=_radius_control_sys,
)
