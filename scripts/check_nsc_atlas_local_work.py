"""One fixed-state replay for the local fourth-order gravity constraint."""
import argparse,hashlib,json,time
from pathlib import Path
from threadpoolctl import threadpool_limits
from recursive_horizons.nsc_atlas_flow_evolution import PeriodicFlow
from recursive_horizons.nsc_atlas_local_work import summary
from recursive_horizons.nsc_atlas_rindler_bulk import checks,exact_radial_check
ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    out=(ROOT/parser.parse_args().output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    source=ROOT/'results/summaries/nsc-atlas-flow-evolution-v1.json'
    saved=json.loads(source.read_bytes());case=saved['cases'][-1]
    start=time.perf_counter();cpu=time.process_time();rows=[]
    with threadpool_limits(limits=1,user_api='blas'):
        f=PeriodicFlow(case['n'],workers=8);h=f.prepare();loss=0.
        rows.append(dict(t=0.,**summary(f,h)))
        for step in range(case['steps']):
            h,inc=f.step(h,case['dt']);loss+=inc
            if time.process_time()-cpu>900:raise RuntimeError('local-constraint CPU budget exhausted; no second campaign')
        d=f.diagnostics(h)
        replay={key:d[key]-case['final'][key] for key in ('energy','dissipation')}
        if any(abs(v)>1e-8 for v in replay.values()):raise ValueError('fixed-state replay differs from the retained result')
        rows.append(dict(t=case['endpoint'],**summary(f,h)))
    numeric=checks();exact=exact_radial_check()
    if not exact['all_coefficients_zero']:raise ValueError('exact bulk tensor control failed')
    result=dict(schema='nsc-atlas-local-work-v1',
        source='Compere--McFadden--Skenderis--Taylor 1103.3022v2, equations 3.2,4.2-4.17,6.1,6.5-6.8',
        source_formula='S=(u.grad p)/nu-u.Delta u+sigma_ij sigma_ij/2; sigma_ij=di uj+dj ui',
        average_identity='mean(S)=2D/(Vnu)=3 b2_prime',
        local_completion='w=b2_prime*y+grad chi, Delta chi=S-mean(S), mean(chi)=0',
        domain='Periodic leading profile; higher-order correction is affine-periodic on an open patch, not a fixed periodic vector field.',
        states=rows,retained_flow_replay_difference=replay,
        direct_bulk_tensor_controls=numeric,exact_radial_tensor_control=exact,
        tensor_scope='Published fourth-order metric independently contracted for an analytic decaying shear witness, including its affine correction. The full numerical vortex metric is not independently contracted here.',
        coordinate_scope='A time-dependent dilation of a flat boundary includes shift/lapse terms. Matching the average coefficient does not prove identity with the earlier shift-free expanding-boundary ansatz.',
        solved=['Mean obstruction and required affine expansion coefficient',
                'Local scalar constraint on the retained band-limited initial/final flow states',
                'Independent finite-order Ricci check of a horizon-regular local bulk witness'],
        remaining=['Fifth-order vector evolution of the correction and its transverse freedom',
                   'Global compact/parent attachment and inherited field-state compatibility',
                   'Controlled finite-amplitude remainder and child-region formation'],
        budget=dict(fixed_replay_only=True,new_spatial_refinement=False,FFT_workers=8,
            source_padding_factor=2,CPU_seconds_limit=900,full_arrays_saved=False),
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start)
    if result['cpu_seconds']>900:raise RuntimeError('local constraint exceeded declared compute budget')
    result['input_sha256']={source.relative_to(ROOT).as_posix():hashlib.sha256(source.read_bytes()).hexdigest()}
    names=('src/recursive_horizons/nsc_atlas_local_work.py','src/recursive_horizons/nsc_atlas_rindler_bulk.py',
           'src/recursive_horizons/nsc_atlas_flow_evolution.py','scripts/check_nsc_atlas_local_work.py')
    result['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(states=[{k:r[k] for k in ('t','mean','expected_mean','affine_rate','residual_rms','local_residual_relative_rms','periodic_correction_rms')} for r in rows],
        replay=replay,exact_bulk=exact,cpu=result['cpu_seconds'],wall=result['wall_seconds']),indent=2))


if __name__=='__main__':main()
