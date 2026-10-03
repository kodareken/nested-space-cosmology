#!/usr/bin/env python3
"""Read completed locked/held cut response; create one new readonly assessment."""
import argparse, hashlib, json, os, sys
from pathlib import Path
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name,'1')
import numpy as np
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_parent_cut_response as cut

LIMIT=64*1024*1024
COMMIT='7004a7915a9375c2d5a420853158c4cf19bcdca1'
STAGES=('prepare','prediction','measurement-a0.05-T2.25','measurement-a-0.05-T2.25')


def digest(path):return cut.parent._sha_file(path)
def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def authenticate(record):
    if record.get('producing_commit')!=COMMIT or not record.get('producers'):
        raise ValueError('stage does not bind the accepted frozen producer')
    for path,sha in record['producers'].items():
        cut.provenance.resolve_pinned_source_bytes(cut.ROOT,path,sha,commit=COMMIT)
    for path,sha in record['inputs'].items():
        if digest(path)!=sha:raise ValueError('baseline input hash changed: '+path)


def columns(pair,state,tangent=None,delta_tau=0.):
    fine=cut.leading.fine_state(pair,state);c=pair.weights;dx=pair.grid.dx_q
    norm=abs(fine.phi0)**2+abs(fine.phi1)**2
    integral=lambda a:np.array([cut.extent.real_interval_integral(pair.grid,a[:,j],pair.child_interval) for j in range(a.shape[1])])
    content=integral(norm*c/dx)
    if tangent is None:return content
    lifted=cut.leading.prolong(pair.grid,cut.leading.decode(pair,cut.response._as_state(tangent)))
    variation=norm*tangent.occupations+2*c*np.real(fine.phi0.conj()*lifted.phi0+fine.phi1.conj()*lifted.phi1)
    velocity=cut.leading.prolong(pair.grid,cut.leading.decode(pair,cut.response._as_state(cut.response.velocity_tangent(pair,state))))
    slope=integral(2*c*np.real(fine.phi0.conj()*velocity.phi0+fine.phi1.conj()*velocity.phi1)/dx)
    rate=cut.response.centre_clock(pair,state)
    if not np.isfinite(rate) or rate<=0:raise ValueError('nonpositive centre clock')
    return content,integral(variation/dx)-slope*delta_tau/rate


def effect(predicted,measured):
    residual=measured-predicted
    return dict(predicted_change=predicted,measured_change=measured,residual=residual,
                relative_to_predicted_effect=None if predicted==0 else abs(residual/predicted),
                relative_to_measured_effect=None if measured==0 else abs(residual/measured))


def report(directory=cut.OUTPUT):
    directory=Path(directory).resolve();hashes={};stages={};array_hashes={}
    owners=[m for n,m in sys.modules.items() if n.startswith('recursive_horizons.') and getattr(m,'__file__',None)]
    sources={str(Path(m.__file__).resolve()):digest(m.__file__) for m in owners};sources[str(Path(__file__).resolve())]=digest(__file__)
    for stem in STAGES:
        for suffix in ('.json','.npz'):
            path=directory/(stem+suffix)
            if path.stat().st_size>LIMIT:raise ValueError('stage must be <=64 MiB')
            hashes[str(path)]=digest(path)
        record,arrays=cut.response._load_stage(directory,stem);authenticate(record);stages[stem]=(record,arrays)
        hashes.update(record['inputs']);array_hashes[stem]={k:cut.parent._sha_array(v) for k,v in arrays.items()}
    prepared,pa=stages['prepare'];locked,la=stages['prediction']
    baseline,ba=cut.parent._load(Path(next(iter(prepared['inputs']))).parent)
    for path,sha in baseline['input_hashes'].items():hashes[str(cut.ROOT/path)]=sha
    if baseline!=prepared['baseline_record']:raise ValueError('baseline record differs from preparation')
    if not locked.get('forecast_locked') or locked['forecast']['stop']!='COMPLETE' or locked['prepare_json_sha256']!=hashes[str(directory/'prepare.json')]:raise ValueError('forecast lock/stage link invalid')
    if prepared['declared_held_amplitudes']!=[-.05,.05]:raise ValueError('held amplitudes differ')
    for n,sha in prepared['state_hashes'].items():
        if cut.parent._sha_array(pa[n])!=sha:raise ValueError('initial state array digest failed: '+n)
    if cut.parent._sha_array(pa['delta_c'])!=prepared['direction_hash']:raise ValueError('direction array digest failed')
    for n in ('Q','r','phi0','phi1','W','reference_columns','source_columns'):
        if not np.array_equal(pa[n],ba[n]):raise ValueError('baseline preparation changed: '+n)
    target=locked['forecast']['rows'][-1]
    if abs(target['time']-2.25)>1e-9:raise ValueError('forecast endpoint is not T2.25')
    holder=lambda a:cut.reconstruct(baseline,dict(a,pi_Q=a['p_Q'],pi_r=a['p_r']))
    arms={};held_columns=[]
    with threadpool_limits(limits=1),cut.parent.backend.fft_thread_limit(1):
        pair,state=holder(la);tangent=cut.response.ParentTangent(*(la['delta_'+n] for n in cut.leading.FIELDS),la['delta_c'])
        values,delta=cut.observables(pair,state,tangent)
        _,slope=cut.observables(pair,state,cut.response.velocity_tangent(pair,state))
        equal={n:cut.response.centre_clock_variation(delta[n],slope[n],target['delta_tau'],values['centre_clock_rate']) for n in values}
        content,column_delta=columns(pair,state,tangent,target['delta_tau'])
        if not np.isclose(sum(column_delta),equal['child_content'],rtol=1e-11,atol=1e-13):raise ValueError('column tangent sum differs from total')
        for n in values:
            if not np.isclose(values[n],target['values'][n],rtol=1e-10,atol=1e-12) or not np.isclose(equal[n],target['delta_equal_tau'][n],rtol=1e-10,atol=1e-12):raise ValueError('locked readout differs: '+n)
        for stem in STAGES[2:]:
            r,a=stages[stem];alpha=r['amplitude'];run=r['measurement']
            if run['stop']!='COMPLETE' or not r['comparison_valid'] or not r['matched_proper_clock'] or abs(run['tau']-target['tau'])>1e-9:raise ValueError('held arm has no matched completed clock')
            if r['prediction_json_sha256']!=hashes[str(directory/'prediction.json')]:raise ValueError('held prediction link differs')
            if not r['held_preparation']['converged'] or not r['held_preparation']['hard_anchor']['satisfied']:raise ValueError('held preparation unresolved')
            for n in ('W','reference_columns','source_columns'):
                if not np.array_equal(a[n],pa[n]):raise ValueError('held source/frame identity changed: '+n)
            if not np.array_equal(a['weights'],pa['weights']+alpha*pa['delta_c']):raise ValueError('held occupations differ')
            hp,hs=holder(a);measured,_=cut.observables(hp,hs);cc=columns(hp,hs);held_columns.append((alpha,cc))
            if not np.isclose(cc.sum(),measured['child_content'],rtol=1e-11,atol=1e-13):raise ValueError('held column sum differs')
            for n in measured:
                if not np.isclose(measured[n],r['values'][n],rtol=1e-10,atol=1e-12):raise ValueError('held array readout differs: '+n)
            arms[str(alpha)]={'time':run['time'],'tau':run['tau'],'clock_residual':run['tau']-target['tau'],
                'headline':{n:effect(alpha*equal[n],measured[n]-values[n]) for n in values},
                'columns':[effect(float(alpha*d),float(x-b)) for d,x,b in zip(column_delta,cc,content)],
                'column_norms':np.sum(abs(hs.phi0)**2+abs(hs.phi1)**2,axis=0).tolist(),
                'CAR':cut.car(hp,hp.weights,hs),'energy':cut.leading.energy(hp,hs),'CPU_seconds':r['run_CPU_seconds']}
        byalpha=dict(held_columns);h=prepared['alpha'];centered=(byalpha[h]-byalpha[-h])/(2*h);even=(byalpha[h]+byalpha[-h])/2-content
        total_cpu=locked['CPU_seconds']+sum(r['run_CPU_seconds'] for r,a in [stages[s] for s in STAGES[2:]])
        if total_cpu>cut.CPU_LIMIT:raise ValueError('completed batch exceeded CPU allowance')
        body={'schema':'NSC-DISCOVERY-PARENT-CUT-ASSESSMENT-v1','domain':{'nf':baseline['nf'],'k':baseline['k'],'cuts':list(cut.CUTS),'child_interval':list(pair.child_interval),'target_tau':target['tau'],
             'finite_band':True,'spatial_confirmation_measured':False,'historical_producing_commit':COMMIT},
             'baseline':{'initial_child_content':prepared['initial_row']['values']['child_content'],'values':values,'CAR':cut.car(pair,pair.weights,state),'energy':cut.leading.energy(pair,state),'column_content':content.tolist(),
                         'column_norms':np.sum(abs(state.phi0)**2+abs(state.phi1)**2,axis=0).tolist()},
             'locked_equal_tau_derivative':equal,'column_equal_tau_derivative':column_delta.tolist(),'held':arms,
             'column_centered_slopes':centered.tolist(),'column_centered_relative_residual':(abs((centered-column_delta)/column_delta)).tolist(),'column_even_remainder':even.tolist(),
             'CPU':{'sum_stage_seconds':total_cpu,'saved_aggregate_seconds':max(stages[s][0]['aggregate_CPU_seconds'] for s in STAGES[2:]),'allowance_seconds':cut.CPU_LIMIT},'input_sha256':hashes,'array_sha256':array_hashes,'measurement_source_sha256':sources,
             'historical_source_current_differences':[p for p,s in prepared['producers'].items() if digest(cut.ROOT/p)!=s],
             'evolved':False,'renewal_asserted':False,'column_tags':'original computational ancestry, not particles or energy',
             'frozen_column0_derivative':0.,'frozen_column0_scope':'unchanged generator and fixed c0; no frozen run requested',
             'limitations':['Initial held Q/r/Phi fixed by authenticated preparation owner; momenta re-solved.','Baseline child remains depleted; finite-band relation has no continuum/EFT or autonomous-renewal proof.'],
             'next_physical_question':'Does the original-child-column feedback survive spatial confirmation at the same proper clock?'}
    for path,sha in {**hashes,**sources}.items():
        if digest(path)!=sha:raise ValueError('input/evaluator changed during assessment: '+path)
    body['content_sha256']=hashlib.sha256(canonical(body).encode()).hexdigest();return body


def write_record(path,result,directory=cut.OUTPUT):
    path=Path(path).expanduser().resolve()
    if path.is_relative_to(Path(directory).resolve()):raise ValueError('output must stay outside sealed stage directory')
    payload=(json.dumps(result,indent=2,allow_nan=False)+'\n').encode()
    if len(payload)>LIMIT:raise ValueError('output exceeds 64 MiB')
    with path.open('xb') as stream:stream.write(payload)
    path.chmod(0o444)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path);args=parser.parse_args(argv)
    try:
        if args.output and args.output.exists():raise FileExistsError('refusing to overwrite assessment output')
        result=report()
        if args.output:write_record(args.output,result)
        print(json.dumps(result,indent=2,allow_nan=False));return 0
    except (ValueError,OSError,RuntimeError,KeyError) as error:
        print(type(error).__name__+': '+str(error),file=sys.stderr);return 2


if __name__=='__main__':raise SystemExit(main())
