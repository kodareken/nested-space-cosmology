"""Probe where the vacuum-correction source method remains useful.

No certificate or source payload is modified. Every case has a 30-second
capture CPU limit. Printed bounds are diagnostic, with no coverage promotion.
"""
from pathlib import Path
import time,json,argparse
from flint import arb,ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection,capture_correction,validate_correction
from recursive_horizons.nsc_subgap_source_covariance import original_covariance_distance_bounds
parser=argparse.ArgumentParser(description='Bounded preparation-method diagnostic; not a gate certificate.')
parser.add_argument('--run',action='store_true',required=True)
parser.parse_args()
archive=RetainedUpstreamArchive(Path.cwd()); config=archive.meta['config']
entries=[b for b,_ in archive.family_entries((14,1)) if b.energy_sign>0]
print(json.dumps(dict(kind="DIAGNOSTIC_ONLY", input_hashes=archive.input_hashes, physical_local_gate="OPEN")),flush=True)
for threshold,step in ((8.,.1),(4.,.1),(2.,.1),(4.,.025)):
    E,b,j=min(((float(E),b,j) for b in entries for j,E in enumerate(b.source.energies[::3]) if E>=threshold),key=lambda v:v[0])
    start=time.process_time()
    try:
        with ctx.workprec(192):
            expansion=VacuumSourceExpansion(config['horizon_rho'],arb(b.mass).union(arb.pi()/2),arb(b.angular).union(arb(5).sqrt()),order=8)
            flow=VacuumCorrection(expansion,E); end=float(expansion.target_distance(b.rho_up).log().mid())
            trace,nfev=capture_correction(flow,-18.,end,max_step=step,cpu_limit=30.)
            proof=validate_correction(flow,trace,-18.,b.rho_up)
            sl=slice(3*j,3*j+3)
            distance=original_covariance_distance_bounds(b.initial_columns[:,sl],b.source.covariance[sl,sl],proof)
            thermal=expansion.covariance_error(0,E,kappa=config['surface_gravity'],omega=config['omega'])['thermal_coherent_error']
            print(json.dumps(dict(max_step=step,E=E,panel=b.original_panel,row=b.rows[0]+j,cells=len(trace),vacuum_error=float(proof['bloch_error']),thermal=float(thermal),covariance_upper=float((distance['upper']+thermal).upper()),cpu_seconds=time.process_time()-start)),flush=True)
    except Exception as exc: print(json.dumps(dict(max_step=step,E=E,error=str(exc),cpu_seconds=time.process_time()-start)),flush=True)
