#!/usr/bin/env python3
"""Compose authenticated full-energy PG covariances without replaying old runs.

--prepare consumes an explicit component-path manifest; it performs no mode
integration. --check recomputes every recorded numerical output from payloads.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from recursive_horizons.nsc_pg_retained_covariance import positive_panel,signed_covariance,covariance_residuals
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP as S
from recursive_horizons.nsc_pg_state_covariance import equal_time_ctp_blocks
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_pg_group14_covariance import assemble as assemble_group14
from derive_nsc_transmitting_boundary_binding import compare


def real_panel_quadrature(mass, upper=40., points=16, panels=None):
    """Real-frequency low panels exclude the separately retained (1,m) subgap."""
    from numpy.polynomial.legendre import leggauss
    if mass < 0 or upper <= 1 or points < 1:
        raise ValueError("nonnegative mass, upper>1 and positive quadrature order required")
    # Preserve the original low-window landmarks and exact panel membership.
    edges = sorted(set(min(float(v), upper) for v in
                       (0., .25, 1., mass, 4., 8., 12., 16., upper)))
    edges = [value for index, value in enumerate(edges)
             if index == 0 or value-edges[index-1] > 1e-12]
    allowed = list(zip(edges, edges[1:]))
    chosen = allowed if panels is None else list(panels)
    if any((left, right) not in allowed for left, right in chosen):
        raise ValueError("refinement must use already declared panel endpoints")
    x, weights = leggauss(points)
    energy, weight = [], []
    for left, right in chosen:
        if right <= left or (left >= 1. and right <= mass):
            continue
        energy.extend(left+(x+1)*(right-left)/2)
        weight.extend(weights*(right-left)/2)
    return np.asarray(energy), np.asarray(weight), edges


def middle_panel_quadrature(channel, left=40., right=160., points=16):
    """Eight-unit middle panels retain the validated group14 E=16 exception."""
    from numpy.polynomial.legendre import leggauss
    if left < 40 and not (channel == 14 and left == 16):
        raise ValueError("lower window requires a validated boundary-mode control")
    if right <= left or (right-left) % 8 or points < 1:
        raise ValueError("ordered eight-unit quadrature panels required")
    x, weight = leggauss(points)
    starts = np.arange(left, right, 8.)
    return (np.concatenate([value+4*(x+1) for value in starts]),
            np.tile(4*weight, len(starts)))


def complete_archived_real_rows(inverse_modes, accepted, archived_modes):
    """Complete rejected inverse rows only; accepted inverse bytes stay unchanged."""
    inverse = np.asarray(inverse_modes)
    accepted = np.asarray(accepted, dtype=bool)
    archived = np.asarray(archived_modes)
    if (inverse.shape != archived.shape or accepted.shape != (len(inverse),)
            or not np.isfinite(archived[~accepted]).all()):
        raise ValueError("complete finite archived fields and matching acceptance mask required")
    result = inverse.copy()
    result[~accepted] = archived[~accepted]
    if not np.array_equal(result[accepted], inverse[accepted]):
        raise ValueError("completion changed accepted inverse rows")
    if not np.isfinite(result).all():
        raise ValueError("completion left missing real rows")
    return result


def batched_real_input(energies, weights, channel, sign, preparation):
    """Physical horizon/infinity modes, not fitted seed covariance, feed low panels."""
    from recursive_horizons.nsc_pg_batched_modes import real_mode_family
    from recursive_horizons.nsc_pg_batch_projection import project_mode_family
    if sign not in (-1, 1) or (not channel["angular_eigenvalue"] and sign == -1):
        raise ValueError("an admissible signed retained family is required")
    energies, weights = np.asarray(energies), np.asarray(weights)
    if weights.shape != energies.shape or energies.ndim != 1:
        raise ValueError("matching real-frequency energies and quadrature weights required")
    mass, angular = channel["compact_mass"], sign*channel["angular_eigenvalue"]
    if np.any((energies > 1.) & (energies < mass)):
        raise ValueError("the retained subgap must not be counted in real panels")
    modes = real_mode_family(energies, mass, angular, preparation)
    projection = project_mode_family(energies, mass, angular,
                                     modes["mode_at_zero"], modes["mode_at_three"])
    return {"energy": energies, "weight": weights, "projection": projection,
            "source": modes["source"], "projector": modes["projector"],
            "current_residual": modes["current_residual"],
            "inner_unitarity_residual": modes["inner_unitarity_residual"],
            "seed_covariance_used_as_input": False}


def retained_family_projection(arrays, channel, sign):
    """Reuse authenticated A2 fields for a non-subgap signed packet projection."""
    from recursive_horizons.nsc_pg_batch_projection import project_mode_family
    labels = np.asarray(arrays["label"])
    selected = ((labels[:, 0] == channel) & (labels[:, 1] == sign)
                & ((labels[:, 5] < 1.) | (labels[:, 5] > labels[:, 3])))
    indices = np.flatnonzero(selected)
    if not len(indices):
        raise ValueError("empty declared projection family")
    rows = labels[indices]
    mass, angular = rows[0, 3:5]
    projection = project_mode_family(rows[:, 5], mass, angular,
                                     arrays["interior"][indices, 1],
                                     arrays["exterior"][indices, 1])
    return {"indices": indices, "labels": rows, "projection": projection,
            "source_covariance": arrays["source_covariance"][indices],
            "source_projector": arrays["source_projector"][indices]}


INPUT_METHODS = {
    "low_real": "same physical horizon/infinity preparation, batched or independent real mode equations; exclude compact (1,m) subgap",
    "middle": "eight-unit panels, normally E=40..160; validated group14 lower boundary E=16",
    "archived_high": "order16 Riccati fields fill rejected inverse rows while accepted rows and source covariance stay unchanged",
    "projection": "reuse authenticated A2 physical fields; keep signed family and original preparation",
    "scope": "finite panel coverage and order/refinement diagnostics; full differentiated subgap and infinite-tail source remain OPEN",
    "data": "https://doi.org/10.5281/zenodo.23090015",
}

OUTPUT='results/development/nsc-pg-retained-covariance.json'
INPUTS=('results/development/nsc-pg-group13-covariance.json',
 'results/development/nsc-pg-retained-subgap.json','results/development/nsc-pg-retained-tail.json',
 'results/development/nsc-pg-lll-preparation.json','results/development/nsc-pg-massive-mode-resolution.json',
 'results/development/nsc-mode-resolved-cauchy-state.json','results/development/nsc-compact-matched-restart.json',
 'results/development/nsc-transmitting-ctp-resolvent.json')
SOURCES=('src/recursive_horizons/nsc_pg_retained_covariance.py',
 'src/recursive_horizons/nsc_pg_batched_modes.py','src/recursive_horizons/nsc_pg_batch_projection.py',
 'scripts/derive_nsc_pg_retained_covariance.py','scripts/derive_nsc_pg_group14_covariance.py',
 'docs/nsc-pg-retained-join-audit.json',
 'tests/test_nsc_pg_retained_covariance.py','docs/nsc-pg-retained-covariance.md')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):
    with np.load(path,allow_pickle=False) as a:return {k:a[k].copy() for k in a.files}
def rec(name):return json.loads((ROOT/name).read_text())
def pairs(x):return np.stack((x.real,x.imag),axis=-1).tolist()
def unpairs(x):
    a=np.asarray(x);return a[...,0]+1j*a[...,1]
def authenticated_payload(record):
    p=ROOT/record['payload']['path']
    if sha(p)!=record['payload']['sha256']:raise ValueError('component payload changed: '+str(p))
    return read(p)
def unpack(a,prefix):return {k[len(prefix)+1:]:v for k,v in a.items() if k.startswith(prefix+'/')}

def prepare(manifest_path):
    manifest=json.loads(Path(manifest_path).read_text());arrays={};receipts=[]
    for item in manifest['components']:
        path=Path(item['path']);name=item['name'];data=read(path)
        for k,v in data.items():arrays[name+'/'+k]=v
        receipts.append({'name':name,'sha256':sha(path),'bytes':path.stat().st_size})
    metadata={'schema':'NSC-PG-RETAINED-COVARIANCE-INPUT-v1','components':receipts,
       'source_hashes':{s:sha(ROOT/s) for s in (
          'src/recursive_horizons/nsc_pg_batched_modes.py','src/recursive_horizons/nsc_pg_batch_projection.py',
          'scripts/derive_nsc_pg_retained_covariance.py')},
       'input_hashes':{s:sha(ROOT/s) for s in INPUTS},
       'seed_covariance_used_as_input':False,'start_convention':'matched_delta_q',
       'group14_input':manifest['group14_input'],
       'mode_controls':manifest.get('mode_controls',{}),
       'quadrature_scope':'finite numerical frequency panels plus separately authenticated subgap and infinite-tail integrals'}
    arrays['metadata_json']=np.frombuffer(json.dumps(metadata,sort_keys=True).encode(),np.uint8)
    raw=deterministic_npz_bytes(arrays);digest=hashlib.sha256(raw).hexdigest()
    path=ROOT/f'results/development/artifacts/nsc-pg-retained-covariance-input.{digest}.npz'
    path.write_bytes(raw);return path

def integrated(a,prefix,low=True):
    d=unpack(a,prefix)
    if not d:raise ValueError('missing physical component '+prefix)
    meta=json.loads(d['metadata_json'].tobytes())
    if low:
        if meta.get('seed_covariance_used_as_input') is not False:raise ValueError('seed input is forbidden')
        return positive_panel(d['projection'],d['weight'],d['source'],d['projector'])
    return positive_panel(d['projection'],d['weights'])

def assemble(a,metadata):
    tail=authenticated_payload(rec(INPUTS[2]));sub=authenticated_payload(rec(INPUTS[1]))
    channels=rec(INPUTS[5])['channels'];positive={};diagnostics={};bounds={}
    config=rec(INPUTS[6])['scattering_provenance']['config']
    low_deltas={};mid_deltas={};coarse_deltas={}
    for ch in channels[1:]:
        g=ch['index']
        if g in (13,14):continue
        signs=(1,-1) if ch['angular_eigenvalue'] else (1,)
        for s in signs:
            label=f'{g}_{s}';low=integrated(a,'low/'+label);mid=integrated(a,'mid/'+label,False)
            lm=json.loads(unpack(a,'low/'+label)['metadata_json'].tobytes())
            mm=json.loads(unpack(a,'mid/'+label)['metadata_json'].tobytes())
            ix=np.flatnonzero((tail['group']==g)&(tail['sign']==s))
            if len(ix)!=1:raise ValueError('missing signed tail')
            i=ix[0];L=float(tail['selected_lower'][i])
            if lm['channel']!=ch or lm['sign']!=s or lm['upper']!=40 or mm['left']!=40 or mm['right']!=L:
                raise ValueError('channel or low/middle/tail domain mismatch')
            E=unpack(a,'low/'+label)['energy'];m=ch['compact_mass']
            if np.any((E>1)&(E<m)):raise ValueError('subgap double counted')
            ld=unpack(a,'low/'+label)
            expected_C=np.array([source_covariance(e,config['surface_gravity'],config['omega'],m) for e in E])
            expected_P=np.array([np.diag([1.,1.,float(e>m)]) for e in E])
            source_residual=float(max(np.max(abs(ld['source']-expected_C)),np.max(abs(ld['projector']-expected_P))))
            if source_residual>3e-13:raise ValueError('low covariance differs from physical horizon/incoming source law')
            G=low[0]+mid[0]+tail['gram'][i];K=low[1]+mid[1]+tail['centered'][i]
            j=np.flatnonzero((sub['labels'][:,0]==g)&(sub['labels'][:,1]==s))
            if m:
                if len(j)!=1:raise ValueError('missing compact subgap')
                G+=sub['gram_positive'][j[0]]/(2*np.pi);K+=sub['centered_positive'][j[0]]/(2*np.pi)
            positive[g,s]=(G,K);bounds[g,s]=(float(tail['bound_C'][i]),float(tail['bound_G'][i]))
            diagnostics[g,s]={'source_law_residual':source_residual,'low_current':lm['current_residual'],'low_inner_unitarity':lm['inner_unitarity_residual'],
                               'low_nodes':len(E),'middle_nodes':len(unpack(a,'mid/'+label)['weights']),'tail_transition':L}
            for prefix,destination,base,islow in (('low_ref',low_deltas,low,True),('mid_ref',mid_deltas,mid,False),('low_coarse',coarse_deltas,low,True)):
                if unpack(a,prefix+'/'+label):
                    v=integrated(a,prefix+'/'+label,islow);destination[g,s]=(v[0]-base[0],v[1]-base[1])
    results={};checks={}
    for key in positive:
        g,s=key;partner=(g,-s) if channels[g]['angular_eigenvalue'] else key
        results[key]=signed_covariance(positive[key],positive[partner]);row=covariance_residuals(results[key]);row.update(diagnostics[key])
        row['tail_Fourier_normalization_bound_C']=bounds[key][0]+bounds[partner][0]
        row['tail_Fourier_normalization_bound_Gram']=bounds[key][1]+bounds[partner][1]
        for name,deltas in (('low_refinement',low_deltas),('middle_refinement',mid_deltas),('rejected_coarse_difference',coarse_deltas)):
            if key in deltas and partner in deltas:
                d=signed_covariance(deltas[key],deltas[partner])
                row[name+'_C']=float(np.linalg.norm(d['covariance'],2));row[name+'_Gram']=float(np.linalg.norm(d['car_gram'],2))
        checks[key]=row
    # Group13 is an immutable completed input. Group14 has independently
    # prepared low/middle windows and both angular signs, authenticated here.
    old=rec(INPUTS[0]);results[13,1]={k:unpairs(v) for k,v in old['computed'].items()}
    checks[13,1]=covariance_residuals(results[13,1]);checks[13,1]['imported_record']=INPUTS[0]
    p=ROOT/metadata['group14_input']['path']
    if sha(p)!=metadata['group14_input']['sha256']:raise ValueError('group14 input changed')
    d=read(p);d14={}
    for k,v in d.items():
        if '/' in k:
            head,tailkey=k.split('/',1);d14.setdefault(head,{})[tailkey]=v
        else:d14[k]=v
    r14,c14=assemble_group14(d14)
    for s in (1,-1):
        results[14,s]=r14[s];checks[14,s]=c14[s]
        checks[14,s]['tail_Fourier_normalization_bound_C']=float(d14[f'tail_bounds_{s}'][0]+d14[f'tail_bounds_{-s}'][0])
    lll=rec(INPUTS[3]);C=unpairs(lll['projection']['covariance_real_imag']);G=unpairs(lll['projection']['CAR_gram_real_imag'])
    results[0,1]={'covariance':C,'car_gram':G,**equal_time_ctp_blocks(C,G)}
    checks[0,1]=covariance_residuals(results[0,1]);checks[0,1]['imported_record']=INPUTS[3]
    if {k[0] for k in results}!=set(range(33)):raise ValueError('retained-state coverage incomplete')
    return results,checks

def make_record(path,a,metadata):
    for name,digest in metadata['source_hashes'].items():
        if sha(ROOT/name)!=digest:raise ValueError('new input-generation owner changed: '+name)
    for name,digest in metadata['input_hashes'].items():
        if sha(ROOT/name)!=digest:raise ValueError('locked input changed: '+name)
    results,checks=assemble(a,metadata);tol=3e-9
    control=unpack(a,'mode_control')
    if not control:raise ValueError('batched physical-mode control is missing')
    old_modes=authenticated_payload(rec(INPUTS[4]));indices=control['indices']
    mode_checks={
      'at_zero_vs_independent_mode_owner':float(np.max(abs(control['mode_at_zero']-old_modes['interior'][indices,1]))),
      'at_three_vs_independent_mode_owner':float(np.max(abs(control['mode_at_three']-old_modes['exterior'][indices,1]))),
      'current':float(control['current_residual'].item()),
      'inner_unitarity':float(control['inner_unitarity_residual'].item())}
    audit=rec('docs/nsc-pg-retained-join-audit.json');join_checks=[]
    for row in audit['cases']:
        if row['energy']<39:raise ValueError('join control must lie at the actual middle-band start')
        F=unpairs(row['exact_projection']);v=unpairs(row['approx_projection'])
        C=unpairs(row['source']);P=unpairs(row['projector'])
        join_checks.append({'group':row['group'],'sign':row['sign'],'energy':row['energy'],
            'centered':float(np.linalg.norm(F@(C-.5*P)@F.conj().T-v@np.diag([-.5,.5,-.5])@v.conj().T,2)),
            'gram':float(np.linalg.norm(F@P@F.conj().T-v@v.conj().T,2))})
    failures=[]
    for (g,s),v in checks.items():
        for name in ('CAR','Hermiticity','CTP_CAR','CTP_Keldysh','low_refinement_C','middle_refinement_C',
                     'low_refinement_Gram','middle_refinement_Gram','tail_Fourier_normalization_bound_C',
                     'tail_Fourier_normalization_bound_Gram','low_inner_unitarity','low_current'):
            if name in v and (not np.isfinite(v[name]) or v[name]>tol):failures.append([g,s,name,v[name]])
        if v['eigenvalue_min'] < -tol or v['eigenvalue_max'] > 1+tol:failures.append([g,s,'CAR_interval'])
    # Controls cover angular-only, first compact, mixed and maximum compact
    # angular families. They are numerical convergence evidence, not bounds.
    for g in (6,22,32):
        if any('low_refinement_C' not in checks[g,s] for s in (1,-1)):failures.append([g,'missing_low_refinement'])
    for g in (12,22,32):
        if any('middle_refinement_C' not in checks[g,s] for s in (1,-1)):failures.append([g,'missing_middle_refinement'])
    for name,value in mode_checks.items():
        if not np.isfinite(value) or value>tol:failures.append(['batched_mode',name,value])
    for row in join_checks:
        if max(row['centered'],row['gram'])>tol:failures.append(['middle_join',row])
    passed=not failures
    payload={'path':str(path.relative_to(ROOT)),'sha256':sha(path),'bytes':path.stat().st_size}
    return {'schema':'NSC-PG-RETAINED-COVARIANCE-v1','accountable_author':'Douglas Ek',
      'status':'PASS' if passed else 'OPEN','source_hashes':{s:sha(ROOT/s) for s in SOURCES},
      'input_hashes':metadata['input_hashes'],'payload':payload,'auxiliary_payloads':[metadata['group14_input']],
      'locked_inputs':rec(INPUTS[0])['locked_inputs'],
      'domain':{'group_indices':list(range(33)),'massive_groups':32,'surface':'fixed common PG Cauchy surface',
        'multiplicity_owner':'unchanged channel ledger in nsc-mode-resolved-cauchy-state.json; signed evaluations do not add field copies',
        'state':'C_PG=F(C_H direct_sum n_in)F dagger; physical matched horizon/infinity source modes',
        'angular_cross_block_Z':'zero from magnetic angular orthogonality; horizon-pair correlations retained inside C_H',
        'original_packet_coordinates':4,'massive_test_coordinates':8,'LLL_test_coordinates':7,
        'probe_complement':'Q=I-JJdagger retained by spectral_probe_bulk_blocks on the same global mode fields',
        'finite_probes_exhaust_bulk':False,'closed_packet_system':False,'start_convention':'matched_delta_q'},
      'computed':{f'{g}_{s}':{k:pairs(v) for k,v in r.items()} for (g,s),r in sorted(results.items())},
      'per_family':{f'{g}_{s}':v for (g,s),v in sorted(checks.items())},
      'batched_mode_control':{'archived_indices':indices.tolist(),'residuals':mode_checks,'tolerance':tol,
                             'scope':'largest mixed compact/angular group, negative angular sign; four existing frequencies'},
      'middle_join_controls':join_checks,
      'error_account':{'numerical_tolerance':tol,'low_refinement_groups':[6,22,32],
        'middle_refinement_groups':[12,22,32],'group13_controls':'imported immutable full-covariance record',
        'group14_controls':'independent low/middle refinement, both angular signs',
        'analytic_bound_scope':'Fourier endpoint-series and normalization remainder only',
        'mode_error_scope':'numerical mode/order/grid controls; not a uniform spectral error theorem',
        'omitted_source_log_norm_bound_above40':audit['source_correction_log_bound_above40'],
        'omitted_source_bound_rule':'log(6)-min(pi/kappa,2*pi/(Omega*kappa))*40; both energy signs conservatively included',
        'covariance_assembly':'half the independently integrated Gram plus centered integral; no identity completion'},
      'failures':failures,'gate':{'C1b_full_retained':'PASS' if passed else 'OPEN',
        'C1_memory':'imported PASS','transmitting_EndpointBranchJets':'OPEN: full unitary history derivatives and KS endpoint pullback',
        'Gamma_rest':'OPEN','two_sided_Weyl_mismatch':None,'Weyl_diagnostic':93.54264532195464,
        'stationarity':'OPEN','physical_Vc_selected':False,'metric_timestep':False,'stress_computed':False,
        'seed_covariance_used_as_input':False,'LLL_map_assigned_to_massive':False,'physical_history_selected':False,'Z3':'OUT OF SCOPE'},
      'comparison':{'fields':'all scientific fields','float_atol':3e-13,'float_rtol':3e-13,'exact':'labels, paths, hashes, schema, statuses'},
      'reproducer':'python3 scripts/derive_nsc_pg_retained_covariance.py --check'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--prepare',metavar='COMPONENT_MANIFEST');p.add_argument('--input');p.add_argument('--check',action='store_true');a=p.parse_args()
    if a.prepare:path=prepare(a.prepare)
    elif a.input:path=Path(a.input).resolve()
    elif a.check:path=ROOT/rec(OUTPUT)['payload']['path']
    else:p.error('choose --prepare, --input or --check')
    arrays=read(path);metadata=json.loads(arrays['metadata_json'].tobytes())
    if sha(path) not in path.name:raise ValueError('input filename/hash mismatch')
    actual=make_record(path,arrays,metadata)
    if a.check:compare(rec(OUTPUT),actual)
    else:(ROOT/OUTPUT).write_text(json.dumps(actual,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':actual['status'],'families':len(actual['per_family']),'failures':actual['failures'],
        'CAR_max':max(v['CAR'] for v in actual['per_family'].values()),'record':OUTPUT,'payload':actual['payload']},indent=2))

if __name__=='__main__':main()
