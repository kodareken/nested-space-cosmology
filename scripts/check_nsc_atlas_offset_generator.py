"""Retain the one-seed accumulated-offset construction and finite resolution."""
import argparse,hashlib,json,time
from pathlib import Path
from recursive_horizons.nsc_atlas_offset_generator import (
    calculation,winding_calculation,winding_transport_calculation,
    offset_connection_calculation,inherited_offset_calculation,
    intrinsic_compression_calculation,buckling_to_goldilocks_calculation,
    nonlocal_to_local_atlas_readout,parent_receiver_calculation,affine_receiving_birth,
    continuum_receiving_birth,hybrid_receiving_echo,winding_two_neck_calculation,
    split_waveform_calculation,split_onset_calculation,inherited_daughter_calculation)
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--winding-only',action='store_true',help='geometry-only calculation to stdout; no parent or files')
    group.add_argument('--transport-only',action='store_true',help='conservative winding transport to stdout; no parent or files')
    group.add_argument('--connection-only',action='store_true',help='retained isolated-offset comparison to stdout')
    group.add_argument('--inheritance-only',action='store_true',help='parent-derived compression and nesting identity to stdout')
    group.add_argument('--intrinsic-only',action='store_true',help='intrinsic deformation cost and constrained neck motion to stdout')
    group.add_argument('--buckling-only',action='store_true',help='helical buckling, damped +/0/- wave asymmetry and Goldilocks equilibrium to stdout')
    group.add_argument('--atlas-only',action='store_true',help='non-local to local multiverse atlas fixed point and translation dictionary to stdout')
    group.add_argument('--evolution-only',action='store_true',help='finite parent/receiver phase and geometry, with eight parallel controls; stdout only')
    group.add_argument('--birth-only',action='store_true',help='zero-radius affine continuation and seed comparison')
    group.add_argument('--continuum-birth-only',action='store_true',help='zero-radius characteristic continuation')
    group.add_argument('--last-echo-only',action='store_true',help='hybrid birth, future optical bounds and echoes')
    group.add_argument('--neck-split-only',action='store_true',help='winding-two equipotential neck and split')
    group.add_argument('--split-waveform-only',action='store_true',help='coupled neck motion and outgoing phase with six parallel controls')
    group.add_argument('--split-onset-only',action='store_true',help='exact whole-exterior response at neck onset')
    group.add_argument('--cascade-inheritance-only',action='store_true',help='inherited daughter shape, motion and next split direction')
    args=p.parse_args()
    if args.last_echo_only:
        if args.output is not None:p.error('echo calculations write to stdout; omit --output')
        print(json.dumps([hybrid_receiving_echo(samples=n, tolerance=tol, endpoint=end)
                          for n,tol,end in ((65,2e-9,16.),(129,2e-10,20.),(257,2e-11,24.))],indent=2))
        return
    if args.continuum_birth_only:
        if args.output is not None:p.error('birth calculations write to stdout; omit --output')
        print(json.dumps([continuum_receiving_birth(samples=n, emission_offset=.02)
                          for n in (65,129,257)],indent=2))
        return
    if args.birth_only:
        if args.output is not None:p.error('birth calculations write to stdout; omit --output')
        print(json.dumps([affine_receiving_birth(seed=seed) for seed in (1e-4,1e-5,1e-6)],indent=2))
        return
    for selected,run in ((args.winding_only,winding_calculation),
                         (args.transport_only,winding_transport_calculation),
                         (args.connection_only,offset_connection_calculation),
                         (args.inheritance_only,inherited_offset_calculation),
                         (args.intrinsic_only,intrinsic_compression_calculation),
                         (args.buckling_only,buckling_to_goldilocks_calculation),
                         (args.atlas_only,nonlocal_to_local_atlas_readout),
                         (args.neck_split_only,winding_two_neck_calculation),
                         (args.split_waveform_only,split_waveform_calculation),
                         (args.split_onset_only,split_onset_calculation),
                         (args.cascade_inheritance_only,inherited_daughter_calculation),
                         (args.evolution_only,parent_receiver_calculation)):
        if selected:
            if args.output is not None:p.error('geometry-only calculations write to stdout; omit --output')
            print(json.dumps(run(),indent=2,sort_keys=True));return
    if args.output is None:p.error('--output is required for the historical coverage calculation')
    out=(ROOT/args.output).resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh repository output required')
    start=time.process_time();r=calculation();r['cpu_seconds']=time.process_time()-start
    if r['cpu_seconds']>120:raise RuntimeError('finite geometry budget exceeded')
    names=('src/recursive_horizons/nsc_atlas_offset_generator.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
           'scripts/check_nsc_atlas_offset_generator.py')
    r['source_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(parent=r['inherited_parent'],cases=[dict(offset=c['offset'],coverage=[{k:x[k] for k in ('passes','ball_Hausdorff_distance_over_radius','S3_meridian_Hausdorff_distance_over_radius')} for x in c['coverage']]) for c in r['cases']],cpu=r['cpu_seconds']),indent=2))


if __name__=='__main__':main()
