"""Two bounded correction controls with the unchanged incoming occupation gap.

Each capture has a 30-second CPU cap. No original source or record is edited.
"""
from pathlib import Path
import time,json,argparse
from flint import arb,ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection,capture_correction,validate_correction
from recursive_horizons.nsc_subgap_source_covariance import original_covariance_distance_bounds
from recursive_horizons.nsc_source_occupation_enclosure import occupation_vacuum_distance
parser=argparse.ArgumentParser(description='Moderate-energy source diagnostic, not a certificate.')
parser.add_argument('--run',required=True,action='store_true');parser.parse_args()
archive=RetainedUpstreamArchive(Path.cwd()); config=archive.meta['config']
entries=[b for b,_ in archive.family_entries((14,1)) if b.energy_sign>0]
print(json.dumps(dict(scope='positive-row diagnostic only',input_hashes=archive.input_hashes,physical_local_gate='OPEN')),flush=True)
for threshold,order,start_y,step in ((2.,4,-40.,.025),(2.,8,-18.,.00625)):
    E,b,j=min(((float(E),b,j) for b in entries for j,E in enumerate(b.source.energies[::3]) if E>=threshold),key=lambda v:v[0])
    start=time.process_time()
    try:
        with ctx.workprec(192):
            expansion=VacuumSourceExpansion(config['horizon_rho'],arb(b.mass).union(arb.pi()/2),arb(b.angular).union(arb(5).sqrt()),order=order)
            flow=VacuumCorrection(expansion,E); end=float(expansion.target_distance(b.rho_up).log().mid())
            trace,nfev=capture_correction(flow,start_y,end,max_step=step,cpu_limit=30.)
            proof=validate_correction(flow,trace,start_y,b.rho_up)
            sl=slice(3*j,3*j+3)
            distance=original_covariance_distance_bounds(b.initial_columns[:,sl],b.source.covariance[sl,sl],proof)
            thermal=occupation_vacuum_distance(E,expansion.mass,config['surface_gravity'],config['omega'])['operator_distance_upper']
            print(json.dumps(dict(order=order,start_y=start_y,max_step=step,E=E,panel=b.original_panel,row=b.rows[0]+j,cells=len(trace),vacuum_error=float(proof['bloch_error']),thermal=float(thermal),covariance_upper=float((distance['upper']+thermal).upper()),cpu_seconds=time.process_time()-start)),flush=True)
    except Exception as exc: print(json.dumps(dict(E=E,error=str(exc),cpu_seconds=time.process_time()-start)),flush=True)
