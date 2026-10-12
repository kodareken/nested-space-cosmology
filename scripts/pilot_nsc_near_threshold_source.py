"""Reproduce the wider-tube original subgap row control without writing records."""
from pathlib import Path
import time,json,argparse,numpy as np
from flint import arb,ctx
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_paired_horizon_preparation import PairedHorizonSeedMap
from recursive_horizons.nsc_massive_jost_modes import solve_jost
from recursive_horizons.nsc_massive_jost_transport_bound import original_dense_solution,phase_transport_segments
from recursive_horizons.nsc_massive_jost_mixed_transport import subgap_mixed_phase_transport_bound as subgap_phase_transport_bound
from recursive_horizons.nsc_metric_horizon_frame import metric_horizon_frame,reflection_from_phase
from recursive_horizons.nsc_subgap_source_covariance import BlochSource,initial_bloch,capture_bloch,validate_bloch,original_covariance_error
from recursive_horizons.nsc_subgap_upstream_covariance import negative_subgap_covariance_error
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
parser=argparse.ArgumentParser(description='One near-threshold source control; no production record.')
parser.add_argument('--run',required=True,action='store_true');parser.parse_args()
start=time.process_time(); archive=RetainedUpstreamArchive(Path.cwd()); c=archive.meta['config']
p=PairedHorizonSeedMap(c['horizon_rho'],c['surface_gravity'],c['omega'],c['horizon_offset'],c['scattering_tolerance'],c['outer_floor'])
bs=[b for b,_ in archive.family_entries((14,1)) if b.original_panel=='subgap8/14_1' and b.rows==(0,8)]
b=next(b for b in bs if b.energy_sign>0); neg=next(b for b in bs if b.energy_sign<0); j=7; E=float(b.source.energies[3*j]); sl=slice(3*j,3*j+3)
mode=solve_jost(p.background,E,b.mass,b.angular,order=8,radial_collar=1e-10,rtol=2e-13,atol=2e-15)
near=p.horizon_rho+1.01e-4
segments=phase_transport_segments(original_dense_solution(mode.run),p.horizon_rho,near)
print(json.dumps(dict(phase='captured',energy=E,cells=len(segments),cpu=time.process_time()-start)),flush=True)
with ctx.workprec(192):
 m=arb(b.mass).union(arb.pi()/2); ell=arb(b.angular).union(arb(5).sqrt())
 bound=subgap_phase_transport_bound(E,m,ell,p.background,mode=mode,inner_radius=near,degree=16,metric_terms=48,defect_subdivisions=16,tube="0.001")
 error=restored_upper(bound['phase_error_inner_upper']); last=segments[-1]
 rho=arb(p.horizon_rho)+arb(float(last.y_end)).exp(); theta=arb(float(last.theta_end))+arb(0,error)
 frame=metric_horizon_frame(E,m,ell,p.horizon_rho,bits=192)
 R,d,tail=reflection_from_phase(frame,rho,theta)
 y=-18.; init=initial_bloch(frame,R,c['surface_gravity'],y)
 model=BlochSource(frame.q,E,m,ell,bits=192); target=(frame.q-arb.pi()/2-arb(b.rho_up).atan()).log()
 trace,nfev=capture_bloch(model,init,y,float(target.mid()),max_step=.025)
 proof=validate_bloch(model,init,trace,target)
 ep=original_covariance_error(b.initial_columns[:,sl],b.source.covariance[sl,sl],proof)
 en=negative_subgap_covariance_error(neg.initial_columns[:,sl],neg.source.covariance[sl,sl],proof)
 print(json.dumps(dict(energy=E,phase_error=float(error),positive=float(ep),negative=float(en),cells=len(trace),cpu=time.process_time()-start,gate='OPEN')),flush=True)
