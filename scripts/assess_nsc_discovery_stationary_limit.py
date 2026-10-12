#!/usr/bin/env python3
"""Independent bounded assessment of the sealed noncompact stationary query.

No optimization, eigenstate selection, dynamics, or prior evidence mutation.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name, '1')

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
INPUT = LAB / 'results/development/nsc-discovery-stationary-critical-v3.json'
OUTPUT = LAB / 'results/development/nsc-discovery-stationary-limit-v1.json'
PRODUCING_COMMIT = '697efe0155cdbf796c0174468c90a4c368bb7dbe'
SCHEMA = 'NSC-DISCOVERY-STATIONARY-LIMIT-ASSESSMENT-v1'
CPU_BUDGET = 10.
MAX_BYTES = 256 << 10
STATE_FIELDS = ('Q','r','chi','p_Q','p_r','p_chi','phi0','phi1')
BASE_WEIGHTS = np.array([.75,.75,.5,.5,.25,.25])
NUMERICAL_OWNERS = (
    'nsc_spherical_coupling.py', 'nsc_spherical_galerkin_coupling.py',
    'nsc_spherical_feedback_action.py', 'nsc_conformal_adm_source.py',
    'nsc_covariant_operator.py', 'nsc_regulated.py', 'nsc_nested_parent_child.py',
)


def cpu_now():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + child.ru_utime + child.ru_stime


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def repository_path(name):
    name = Path(name)
    if name.is_absolute() or '..' in name.parts:
        raise ValueError('bound references must be repository relative')
    return ROOT / name


def source_hashes():
    paths = [Path(__file__).resolve(), LAB/'tests/test_nsc_discovery_stationary_limit.py',
             LAB/'docs/nsc-discovery-stationary-limit.md']
    paths += [LAB/'src/recursive_horizons'/name for name in NUMERICAL_OWNERS]
    # The historical checker is watched; this assessment does not modify it.
    paths += [LAB/'scripts/derive_nsc_discovery_stationary.py',
              LAB/'src/recursive_horizons/nsc_discovery_stationary.py']
    return {str(path.relative_to(ROOT)): sha256(path) for path in paths}


def authenticate(require_current_numerical=True):
    """Immutable producer/input binding and current byte-authentic payload."""
    record = json.loads(INPUT.read_text())
    if (record.get('producing_commit') != PRODUCING_COMMIT
            or record.get('schema') != 'NSC-DISCOVERY-STATIONARY-CRITICAL-v3'
            or record.get('nf') != 128 or record.get('ng') != 127 or record.get('nq') != 512
            or record.get('status') != 'UNSATISFIED_LOCAL_CRITICAL_RELATIONS'):
        raise ValueError('unexpected sealed critical-v3 domain')
    if not record.get('source_hashes') or not record.get('input_hashes'):
        raise ValueError('sealed critical record lacks producer/input bindings')
    for group in ('source_hashes','input_hashes'):
        for name,digest in record[group].items():
            repository_path(name)
            data = subprocess.check_output(['git','show',PRODUCING_COMMIT+':'+name],cwd=ROOT)
            if hashlib.sha256(data).hexdigest() != digest:
                raise ValueError('immutable critical hash mismatch: '+name)
    if require_current_numerical:
        for name in NUMERICAL_OWNERS:
            path = LAB/'src/recursive_horizons'/name
            if sha256(path) != record['source_hashes'][str(path.relative_to(ROOT))]:
                raise ValueError('assessment requires the original numerical owner: '+name)
    payload_name = Path(record['payload'])
    if payload_name.name != str(payload_name):
        raise ValueError('critical payload must be beside its record')
    payload = INPUT.parent/payload_name
    if sha256(payload) != record['payload_sha256'] or payload.stat().st_size != record['payload_bytes']:
        raise ValueError('critical payload binding mismatch')
    with np.load(payload,allow_pickle=False) as archive:
        arrays = {name:archive[name].copy() for name in archive.files}
    if any(not np.isfinite(values).all() for values in arrays.values()):
        raise ValueError('nonfinite sealed critical array')
    return record,arrays,{
        'producing_commit':PRODUCING_COMMIT, 'record':str(INPUT.relative_to(ROOT)),
        'record_sha256':sha256(INPUT), 'payload':str(payload.relative_to(ROOT)),
        'payload_sha256':sha256(payload), 'source_hashes_authenticated':len(record['source_hashes']),
        'input_hashes_authenticated':len(record['input_hashes']),
        'envelope_and_payload_authenticated':True,
    }


def analytic_limits(coefficients, length=8., harmonic=8, multiplicity=4, weights=BASE_WEIGHTS):
    """Leading resolved-mode lemma, with the AP doublets and M counted once."""
    alpha = float(coupling.alpha_of(coefficients['C_W']))
    A = coefficients['A']
    k = 2*np.pi*harmonic/length
    eps = (2*np.pi/length)*np.array([.5,.5,1.5,1.5,2.5,2.5])
    energy = float(multiplicity*np.dot(weights,eps))
    auxiliary = float(12*alpha*length*k*k)
    return {
        'alpha':alpha,'period':float(length),'harmonic':int(harmonic),'k':float(k),
        'positive_AP_limit_eigenvalues':eps.tolist(),'base_occupations':np.asarray(weights).tolist(),
        'multiplicity_applied_once':multiplicity,'massless_base_field_energy':energy,
        'auxiliary_energy_limit':auxiliary,'source_eta_limit':auxiliary/energy,
        's_over_abs_a_limit':float(k/np.sqrt(6)),
        'r_times_abs_a_limit':float(3*np.sqrt(alpha/(np.pi*A))),
        'rQ_limit':float(k*np.sqrt(3*alpha/(2*np.pi*A))),
        'normal_den_over_s_square_limit':float(16*np.pi*A*length),
        'normal_den_over_a_square_limit':float(16*np.pi*A*length*k*k/6),
        'proper_period_length_limit':float(length*k*np.sqrt(3*alpha/(2*np.pi*A))),
        'finite_band_leading_orders':[0,harmonic,2*harmonic,3*harmonic],
        'new_spectral_trial_performed':False,
    }


def maximum(values):
    return float(np.max(np.abs(values)))


def dot_indicator(first,second,spacing=1.):
    """eps sum of absolute products: an arithmetic indicator, not a bound."""
    products = np.asarray(first)*np.asarray(second)
    return float(np.finfo(float).eps*spacing*np.sum(np.abs(products)))


def state_from_arrays(arrays,prefix=''):
    return coupling.CauchyState(*(arrays[prefix+name].copy() for name in STATE_FIELDS))


def actual_fields(grid,state):
    """Original discrete force on saved columns; never selects an eigensource."""
    rate,bundle = galerkin.compose_fine_hamiltonian(grid,state)
    constraints = coupling.constraint_residuals(bundle['fine_system'],bundle['fine_state'],bundle['source'])
    return rate,bundle,constraints


def stable_einstein_coefficient(grid,Q,R):
    """Isolate Einstein a_r^2 coefficient, without large auxiliary subtraction.

    This is an algebraic coefficient evaluation, not an alternative physical
    candidate: chi/columns/momenta vanish only in this diagnostic call.
    """
    zeros=np.zeros(grid.ng); columns=np.zeros((grid.nf,6),dtype=complex)
    state=coupling.CauchyState(Q.copy(),R.copy(),zeros.copy(),zeros.copy(),zeros.copy(),zeros.copy(),columns.copy(),columns.copy())
    rate,bundle,_=actual_fields(grid,state)
    fine=bundle['fine_state']; system=bundle['fine_system']
    magnetic=-4*np.pi*system.C_F*system.flux**2*fine.Q
    coefficient=rate.p_Q-galerkin.pull_geometry(grid,magnetic)
    source_gap=max(maximum(bundle['source'][name]) for name in ('force_Q','force_L','force_beta'))
    return coefficient,{'isolated_source_force_max':source_gap,'auxiliary_field_zero_for_coefficient_only':True,
                        'original_action_coefficients_changed':False,'candidate_columns_or_momenta_changed':False}


def strict_old_check(record,arrays):
    """Preserve the old check's failure; do not substitute a looser criterion."""
    dx=record['period']/record['ng']; T=arrays['radial_shape_tangents']; B=arrays['radius_pQ_remainder']
    E=record['reduction']['base_source_energy']; saved=arrays['eta_gradient']
    producer=dx*(T.T@B)/E
    checker=(dx*T.T)@B/E
    difference=maximum(checker-saved)
    expected='critical analytic gradient and fixed-eta force disagree'
    command=[sys.executable,str(ROOT/'scripts/lab.py'),'scripts/derive_nsc_discovery_stationary.py',
             '--check','--record',str(INPUT)]
    outcome=subprocess.run(command,cwd=ROOT,capture_output=True,text=True,check=False)
    error=outcome.stderr.strip()
    return {'old_checker_exit_code':outcome.returncode,'old_checker_stderr':error,
            'strict_original_numeric_replay_passed':outcome.returncode==0,
            'numeric_identity_unresolved':outcome.returncode!=0,
            'expected_failure_observed':outcome.returncode!=0 and expected in error,
            'historical_absolute_gradient_threshold':1e-11,
            'producer_order_gap_max':maximum(producer-saved),'checker_order_gap_max':difference,
            'eps_weighted_dot_indicators':[dot_indicator(T[:,j],B,dx)/E for j in range(T.shape[1])],
            'indicator_is_certified_bound':False,'old_checker_modified':False,
            'failure_message_not_replaced_with_success':True}


def saved_sequence(record,limits):
    """Selected improving rows only; scalar summaries, no new trial."""
    improving=[]; score=float('inf')
    for row in record['trials']:
        if row['status']!='MEASURED':continue
        value=float(np.dot(row['gradient'],row['gradient']))
        if value<score:improving.append(row);score=value
    selected=[]
    for threshold in (.17,.08,.03,.01,.003,.001,.0004):
        eligible=[row for row in improving if abs(row['coefficients'][0])<=threshold]
        if eligible and eligible[0]['index'] not in [row['index'] for row in selected]:selected.append(eligible[0])
    final=record['trials'][record['selected_trial_index']]
    if final['index'] not in [row['index'] for row in selected]:selected.append(final)
    rows=[]
    for row in selected:
        a=row['coefficients'][0]; reduction=row['reduction']; den=reduction['normal_radius_denominator']
        radius=np.sqrt(row['radius_amplitude_square'])
        rows.append({'saved_trial_index':row['index'],'a1':a,'a2_over_a1_square':row['coefficients'][1]/a**2,
                     'source_eta':row['source_eta'],'eta_minus_limit':row['source_eta']-limits['source_eta_limit'],
                     'r_times_abs_a':float(radius*abs(a)),'normal_denominator':den,
                     'normal_den_over_a_square':den/a**2,'gradient_max':row['gradient_max'],
                     'actual_pQ_max':row['raw_residual_max']['p_Q'],'actual_lapse_max':row['raw_residual_max']['lapse'],
                     'auxiliary_energy':reduction['auxiliary_term'],'base_source_energy':reduction['base_source_energy']})
    return rows


def assess(producer_commit=None):
    started=cpu_now(); before=source_hashes()
    own_commit=None
    if producer_commit is not None:
        own_commit=subprocess.check_output(['git','rev-parse','--verify',producer_commit+'^{commit}'],cwd=ROOT,text=True).strip()
        for name,digest in before.items():
            data=subprocess.check_output(['git','show',own_commit+':'+name],cwd=ROOT)
            if hashlib.sha256(data).hexdigest()!=digest:raise ValueError('assessment producer bytes do not match commit: '+name)
    record,arrays,binding=authenticate()
    inputs={binding['record']:binding['record_sha256'],binding['payload']:binding['payload_sha256']}
    limits=analytic_limits(record['locked_coefficients'],record['period'],record['preflight']['harmonic'],
                           record['multiplicity'],np.asarray(record['base_occupations']))
    with threadpool_limits(limits=1):
        strict=strict_old_check(record,arrays)
        grid=galerkin.build_grid(record['nf'],gauge='conformal')
        state=state_from_arrays(arrays)
        grid.fine=replace(grid.fine,occupations=np.asarray(record['occupations']))
        rate,bundle,constraints=actual_fields(grid,state)
        R=arrays['radial_shape']; Q=state.Q; fineR=galerkin.prolong_geometry(grid,R)
        stable,coefficient_scope=stable_einstein_coefficient(grid,Q,R)
        initial=state_from_arrays(arrays,'initial_'); initialR=arrays['initial_radial_shape']
        initial_grid=replace(grid,fine=replace(grid.fine,occupations=record['initial_source_eta']*BASE_WEIGHTS))
        unit=initial.copy();unit.r=initialR.copy();double=initial.copy();double.r=2*initialR
        unit_rate,_,_=actual_fields(initial_grid,unit);double_rate,_,_=actual_fields(initial_grid,double)
        initial_ordinary=(double_rate.p_Q-unit_rate.p_Q)/3
        initial_stable,_=stable_einstein_coefficient(initial_grid,initial.Q,initialR)
        # Saved candidate, at unit R, is used only for the affine remainder.
        candidate_unit=state.copy();candidate_unit.r=R.copy()
        candidate_unit_rate,_,_=actual_fields(grid,candidate_unit)
        stable_remainder=candidate_unit_rate.p_Q-stable
        ordinary=arrays['radius_pQ_coefficient'];old_remainder=arrays['radius_pQ_remainder']
        square=record['radius_amplitude_square'];E=record['reduction']['base_source_energy'];s=record['radial_scale']
        dx=grid.dx_g;pair=lambda first,second:float(dx*np.dot(first,second))
        denominator=pair(Q,stable);identity=float(16*np.pi*grid.fine.A*grid.dx_q*np.sum(bundle['fine_state'].Q**2*fineR**2))
        cosine=np.cos(2*np.pi*grid.xi_g[:,None]/grid.length*record['preflight']['harmonic']*np.arange(1,8))
        dQ0=s*cosine;mu=dx*(dQ0.T@stable)/denominator;tangent=dQ0-Q[:,None]*mu
        normal=pair(Q,rate.p_Q)
        conditioned_gradient=dx*(tangent.T@rate.p_Q)/E
        reduced_from_force=s*dx*(cosine.T@rate.p_Q)/E
        stable_envelope_gradient=dx*(tangent.T@stable_remainder)/E
        fine=bundle['fine_state'];a=arrays['shape_coefficients'][0]
        best_match={'a1':float(a),'s':s,'s_over_abs_a':s/abs(a),'r_mean':float(np.mean(state.r)),
                    'r_times_abs_a':float(np.mean(state.r)*abs(a)),'eta':record['source_eta'],
                    'eta_minus_predicted_limit':record['source_eta']-limits['source_eta_limit'],
                    'rQ_mean':float(np.mean(fine.r*fine.Q)),'rQ_min':float(np.min(fine.r*fine.Q)),
                    'rQ_max':float(np.max(fine.r*fine.Q)),'proper_period_length':float(grid.dx_q*np.sum(fine.r*fine.Q)),
                    'radial_shape_first_order_gap_max':maximum(R-(1-a*cosine[:,0]/3)),
                    'a_chi_plus_12_cos_max':maximum(a*state.chi+12*cosine[:,0]),
                    'auxiliary_energy':record['reduction']['auxiliary_term'],'base_source_energy':E}
        affine={'well_conditioned_initial_stable_vs_ordinary_max':maximum(initial_stable-initial_ordinary),
                'selected_stable_vs_stored_ordinary_max':maximum(stable-ordinary),
                'ordinary_affine_at_actual_radius_gap_max':maximum(rate.p_Q-(square*ordinary+old_remainder)),
                'stable_affine_at_actual_radius_gap_max':maximum(rate.p_Q-(square*stable+stable_remainder)),
                'ordinary_affine_terms_max':[maximum(square*ordinary),maximum(old_remainder)],
                'stable_affine_terms_max':[maximum(square*stable),maximum(stable_remainder)],
                'original_normal_denominator':pair(Q,ordinary),'stable_normal_denominator':denominator,
                'positive_normal_identity':identity,'stable_denominator_identity_defect':denominator-identity,
                'normal_den_over_s_square':denominator/s**2,'normal_den_over_a_square':denominator/a**2,
                'actual_normal_Q_force':normal,'actual_normal_Q_force_over_s':normal/s,
                'stored_ordinary_affine_normal_force':pair(Q,square*ordinary+old_remainder),
                'stable_affine_normal_force':pair(Q,square*stable+stable_remainder),
                'stable_normal_square_algebra_only':-pair(Q,stable_remainder)/denominator,
                'actual_candidate_replaced_with_stable_projection':False,
                'eps_normal_dot_indicator':dot_indicator(Q,rate.p_Q,dx),
                'eps_stable_denominator_dot_indicator':dot_indicator(Q,stable,dx),
                'eps_gradient_dot_indicators':[dot_indicator(tangent[:,j],rate.p_Q,dx)/E for j in range(7)],
                'indicators_are_certified_bounds':False,'coefficient_isolation':coefficient_scope}
        forcing={'actual_retained_pQ_max':maximum(rate.p_Q),'actual_fine_pQ_max':maximum(bundle['unprojected_rates'][3]),
                 'actual_retained_pQ_vs_saved_max_gap':maximum(rate.p_Q-arrays['residual_p_Q']),
                 'actual_projected_lapse_max':maximum(galerkin.pull_geometry(grid,constraints['hamilton'])),
                 'actual_full_lapse_max':constraints['hamilton_max'],'actual_full_shift_max':constraints['momentum_max'],
                 'recorded_eta_gradient':arrays['eta_gradient'].tolist(),
                 'conditioned_tangent_gradient_from_actual_force':conditioned_gradient.tolist(),
                 's_times_cosine_force_over_Ebase':reduced_from_force.tolist(),
                 'stable_envelope_gradient_at_unit_radius':stable_envelope_gradient.tolist(),
                 'force_cosine_coefficient_gain':float(2*E/(s*grid.length)),
                 'normal_pairing_term_in_conditioned_gradient':(-mu*normal/E).tolist(),
                 'physical_force_eta_held_fixed':True,'new_eigenstate_or_source_trial':False,
                 'full_force_satisfied':False}
    after=source_hashes()
    if after!=before or any(sha256(repository_path(name))!=digest for name,digest in inputs.items()):
        raise RuntimeError('assessment source/input bytes changed during consumption')
    elapsed=cpu_now()-started
    if elapsed>CPU_BUDGET:raise RuntimeError('stationary-limit assessment exceeded ten aggregate CPU seconds')
    return {'schema':SCHEMA,'status':'AUTHENTIC_NONCOMPACT_LIMIT_WITH_UNSATISFIED_FORCE_AND_UNRESOLVED_STRICT_REPLAY',
            'binding':binding,'strict_original_replay':strict,'analytic_limits':limits,
            'saved_sequence':saved_sequence(record,limits),'selected_saved_match':best_match,
            'einstein_affine_conditioning':affine,'actual_force':forcing,
            'scope':{'assessment_completed':True,'input_envelope_authenticated':True,
                     'original_strict_numeric_replay_passed':strict['strict_original_numeric_replay_passed'],
                     'noncompact_limit_explains_small_coordinate_gradient':True,'finite_branch_nonexistence_proven':False,
                     'stationary_solution_established':False,'holding_or_stability_claim':False,
                     'new_optimizer_eigenstate_source_or_trajectory':False,'old_evidence_or_checker_modified':False},
            'cpu_seconds':elapsed,'cpu_budget_seconds':CPU_BUDGET,'source_hashes':after,'input_hashes':inputs,
            'producing_commit':own_commit,'numpy_version':np.__version__}


def output_path(path=OUTPUT):
    path=Path(path).expanduser()
    if not path.is_absolute() and path.parts and path.parts[0]=='lab':path=ROOT/path
    path=path.resolve();temporary=Path('/tmp').resolve()
    if path!=OUTPUT.resolve() and temporary not in path.parents:
        raise ValueError('assessment output must use its new limit-v1 owner or a /tmp JSON path')
    if path.suffix!='.json':raise ValueError('assessment output must be JSON')
    return path


def write_record(report,path=OUTPUT):
    destination=output_path(path)
    if destination.exists():raise FileExistsError('refusing to replace '+str(destination))
    content=(json.dumps(report,indent=2,allow_nan=False)+'\n').encode()
    if len(content)>MAX_BYTES:raise ValueError('assessment record exceeds 256 KiB')
    destination.parent.mkdir(parents=True,exist_ok=True)
    with destination.open('xb') as handle:handle.write(content)
    return report


def check_record(path=OUTPUT):
    """Authenticate assessment bytes/bindings; retain the prior failure outcome."""
    destination=output_path(path);record=json.loads(destination.read_text())
    if record.get('schema')!=SCHEMA:raise ValueError('unexpected assessment schema')
    if set(record.get('source_hashes',{}))!=set(source_hashes()):raise ValueError('assessment producer set changed')
    critical,unused,authenticated=authenticate(require_current_numerical=False)
    if record.get('binding')!=authenticated:raise ValueError('assessment critical binding changed')
    expected_inputs={authenticated['record']:authenticated['record_sha256'],authenticated['payload']:authenticated['payload_sha256']}
    if record.get('input_hashes')!=expected_inputs:raise ValueError('assessment input set changed')
    for name,digest in record['source_hashes'].items():
        if record.get('producing_commit') is None:
            actual=sha256(repository_path(name))
        else:
            data=subprocess.check_output(['git','show',record['producing_commit']+':'+name],cwd=ROOT)
            actual=hashlib.sha256(data).hexdigest()
        if actual!=digest:raise ValueError('assessment producer binding changed: '+name)
    for name,digest in record['input_hashes'].items():
        if sha256(repository_path(name))!=digest:raise ValueError('assessment input binding changed: '+name)
    if record['binding']['producing_commit']!=PRODUCING_COMMIT:raise ValueError('unexpected critical producing commit')
    if (record['scope']['original_strict_numeric_replay_passed']!=record['strict_original_replay']['strict_original_numeric_replay_passed']
            or record['strict_original_replay']['old_checker_exit_code']==0
            or not record['strict_original_replay']['numeric_identity_unresolved']
            or record['actual_force']['full_force_satisfied']):
        raise ValueError('assessment erased the prior unsatisfied outcome')
    return {'assessment_record_authenticated':True,'schema':SCHEMA,
            'strict_original_numeric_replay_passed':False,'numeric_identity_unresolved':True,
            'full_force_satisfied':False,'bytes_written':0,'numerical_model_evaluated':False}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--record',type=Path,default=OUTPUT)
    parser.add_argument('--producer-commit')
    args=parser.parse_args(argv)
    if args.check:
        if args.write is not None or args.producer_commit:parser.error('--check cannot write or change a producer pin')
        report=check_record(args.record)
    else:
        if args.write is not None and output_path(args.write).exists():raise FileExistsError('refusing to replace '+str(args.write))
        report=assess(args.producer_commit)
        if args.write is not None:write_record(report,args.write)
    sys.stdout.write(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (OSError,ValueError,RuntimeError) as error:
        print(str(error),file=sys.stderr);raise SystemExit(1)
