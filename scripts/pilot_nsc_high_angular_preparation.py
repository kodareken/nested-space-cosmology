"""Compare source-proof methods at three original high-angular rows.

Every correction capture has a 30-second CPU cap. No evidence is overwritten
and no source coverage or gate verdict is promoted by this diagnostic.
"""
from pathlib import Path
import time,json,argparse
from flint import arb,ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_vacuum_source_remainder import VacuumSourceExpansion
from recursive_horizons.nsc_vacuum_source_correction import VacuumCorrection,capture_correction,validate_correction
from recursive_horizons.nsc_subgap_source_covariance import original_covariance_distance_bounds
parser=argparse.ArgumentParser(description='Three original high-angular preparation controls; not a source certificate.')
parser.add_argument('--run',required=True,action='store_true');parser.parse_args()
archive=RetainedUpstreamArchive(Path.cwd());c=archive.meta['config']
print(json.dumps(dict(scope='diagnostic only',input_hashes=archive.input_hashes,physical_local_gate='OPEN')),flush=True)
for group in (12,22,32):
 start=time.process_time();channel=archive.meta['channels'][group]
 entries=[b for b,_ in archive.family_entries((group,1)) if b.energy_sign>0]
 E,b,j=min(((float(e),b,j) for b in entries for j,e in enumerate(b.source.energies[::3]) if e>=32),key=lambda v:v[0])
 try:
  with ctx.workprec(192):
   level=channel['angular_level'];m=arb(b.mass).union(channel['compact_level']*arb.pi()/2)
   ell=arb(b.angular).union(arb(level*(level+4)).sqrt());sl=slice(3*j,3*j+3)
   expansion=VacuumSourceExpansion(c['horizon_rho'],m,ell,order=8)
   flow=VacuumCorrection(expansion,E);end=float(expansion.target_distance(b.rho_up).log().mid())
   trace,nfev=capture_correction(flow,-18.,end,cpu_limit=30.)
   proof=validate_correction(flow,trace,-18.,b.rho_up)
   bound=original_covariance_distance_bounds(b.initial_columns[:,sl],b.source.covariance[sl,sl],proof)
   thermal=expansion.covariance_error(0,E,kappa=c['surface_gravity'],omega=c['omega'])['thermal_coherent_error']
   high=VacuumSourceExpansion(c['horizon_rho'],m,ell,order=16)
   C=high.remainder_constant(b.rho_up,cells=128)
   direct=high.compare_original(b.rho_up,E,b.initial_columns[:,sl],b.source.covariance[sl,sl],C,kappa=c['surface_gravity'],omega=c['omega'])
   print(json.dumps(dict(group=group,E=E,cells=len(trace),vacuum_bloch_error=float(proof['bloch_error']),correction_covariance_upper=float((bound['upper']+thermal).upper()),direct_covariance_upper=float(direct["original_source_operator_error"]),cpu_seconds=time.process_time()-start)),flush=True)
 except Exception as ex:print(json.dumps(dict(group=group,error=str(ex),cpu_seconds=time.process_time()-start)),flush=True)
