"""Bounded, state-pinned phase-vortex and actual momentum-curl check."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_s3_frame import point_values
from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry
from recursive_horizons.nsc_atlas_trapped_surfaces import round_surface,evaluate_surface

ROOT=Path(__file__).resolve().parents[1]
CORE=np.array([[0.,0.,1/math.sqrt(2),1/math.sqrt(2)],
               [0.,1/math.sqrt(2),1/math.sqrt(2),0.]])
NORMAL_AXES=((0,1),(0,3))
CORE_TANGENTS=np.array([[0.,0.,-1/math.sqrt(2),1/math.sqrt(2)],
                        [0.,-1/math.sqrt(2),1/math.sqrt(2),0.]])


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frame(X):
    x,y,z,w=np.asarray(X).T
    return np.stack((np.c_[w,z,-y,-x],np.c_[-z,w,x,-y],np.c_[y,-x,w,-z]),axis=1)


def transverse_loop(core_number,radius,count=128):
    angle=2*math.pi*np.arange(count)/count;X=np.broadcast_to(math.cos(radius)*CORE[core_number],(count,4)).copy()
    a,b=NORMAL_AXES[core_number]
    X[:,a]=math.sin(radius)*np.cos(angle);X[:,b]=math.sin(radius)*np.sin(angle)
    tangent=np.zeros_like(X);tangent[:,a]=-math.sin(radius)*np.sin(angle);tangent[:,b]=math.sin(radius)*np.cos(angle)
    return X,tangent


def winding(Z):
    modulus=abs(Z);steps=np.angle(np.roll(Z,-1)*np.conj(Z));value=float(sum(steps)/(2*math.pi))
    scale=float(max(modulus));floor=128*np.finfo(float).eps*max(1.,scale)
    return dict(sampled_winding=value,nearest_integer=int(round(value)),
        min_abs_Z=float(min(modulus)),max_abs_Z=scale,max_argument_step=float(max(abs(steps))),
        sampled_nonzero_margin=bool(min(modulus)>floor),
        phase_sampling_check=bool(max(abs(steps))<math.pi/2),
        scope=f'{len(Z)} sampled phase increments, not a certified continuous-loop nonvanishing bound.')


def normal_field_derivatives(provider,X):
    """Analytic E_i Pi from spectral first/ordered second scalar derivatives."""
    jet=provider._jets(np.asarray(X)[None,:]);g=jet['g'][0];dg=jet['dg'][0]
    inverse=np.linalg.inv(g);hi=np.linalg.inv(g[1:,1:]);N=1/math.sqrt(-inverse[0,0])
    beta=hi@g[0,1:];row=point_values(X,provider.system.basis.labels)
    C=provider.coefficients;D=provider.system.basis.generators;r=provider.r
    dp=np.array([row@(r*(d@C[22:24].T)) for d in D])
    d2q=np.array([[row@(r*r*(di@(dj@C[20:22].T))) for dj in D] for di in D])
    q=jet['q'][0];p=jet['p'][0];dq=jet['dq'][0]
    Pi=(p-beta@dq)/N;dPi=np.empty((3,2));dN=np.empty(3)
    for k in range(3):
        dhi=-hi@dg[k,1:,1:]@hi
        dbeta=dhi@g[0,1:]+hi@dg[k,0,1:]
        dinverse=-inverse@dg[k]@inverse;dN[k]=N**3*dinverse[0,0]/2
        dPi[k]=(dp[k]-dbeta@dq-beta@d2q[k])/N-Pi*dN[k]/N
    return q,Pi,dq,dPi


def core_readout(provider,number,out):
    X=CORE[number];h=out['h'];inverse=np.linalg.inv(h);K=out['K'];q=out['q'];Pi=out['Pi'];dq=out['dq']
    j_by_field=-dq*Pi[None,:];j=np.sum(j_by_field,axis=1)
    tangent=frame(X[None,:])[0]@CORE_TANGENTS[number]/provider.r
    tangent/=math.sqrt(tangent@h@tangent)
    q2,Pi2,dq2,dPi=normal_field_derivatives(provider,X)
    axes=np.eye(4)[list(NORMAL_AXES[number])]
    normal_vectors=axes@frame(X[None,:])[0].T/provider.r
    dnormal_chi=normal_vectors@dq2[:,1];dnormal_Pi=normal_vectors@dPi[:,1]
    jac=float(np.linalg.det(np.c_[dnormal_chi,dnormal_Pi]))
    Zmatrix=np.vstack((dnormal_chi,dnormal_Pi/math.sqrt(provider.system.parent.mass_squared)))
    singular_values=np.linalg.svd(Zmatrix,compute_uv=False)
    gram=normal_vectors@h@normal_vectors.T;proper_area=math.sqrt(np.linalg.det(gram))
    total_curl=np.zeros((3,3))
    for field in (0,1):total_curl+=np.outer(dq2[:,field],dPi[:,field])-np.outer(dPi[:,field],dq2[:,field])
    curl_normal=float(normal_vectors[0]@total_curl@normal_vectors[1])
    return dict(core=X.tolist(),normal_coordinates=['X'+str(i+1) for i in NORMAL_AXES[number]],
        chi=float(q[1]),Pi_chi=float(Pi[1]),canonical_S=float(q[0]),Pi_S=float(Pi[0]),
        matter_j_covector=j.tolist(),matter_j_norm=float(math.sqrt(max(0.,j@inverse@j))),
        chi_sector_j_covector=j_by_field[:,1].tolist(),volume_sector_j_covector=j_by_field[:,0].tolist(),
        declared_axial_tangent=CORE_TANGENTS[number].tolist(),matter_j_axial_component=float(j@tangent),
        normal_theta=float(-np.einsum('ij,ij',inverse,K)),K_eigenvalues=out['K_eigenvalues'].tolist(),
        future_tangential_capacity=float(out['future_tangential_capacity']),
        transverse_chi_Pi_jacobian=jac,
        transverse_Z_gradient_matrix=Zmatrix.tolist(),
        transverse_Z_gradient_singular_values=singular_values.tolist(),
        predicted_local_winding=int(np.sign(jac)),
        chi_sector_curl_per_proper_transverse_area=jac/proper_area,
        total_curl_per_proper_transverse_area=curl_normal/proper_area,
        normal_field_derivative_consistency_max=float(max(np.max(abs(q2-q)),np.max(abs(Pi2-Pi)),np.max(abs(dq2-dq)))),
        interpretation='A regular scalar-phase core. Momentum curl is not a fluid velocity or a proof of axial drainage or trapping.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--core-surfaces',action='store_true')
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    started=time.perf_counter();budget=180.;cases=[];failure=None
    names=['src/recursive_horizons/nsc_atlas_trapped_parent.py','src/recursive_horizons/nsc_atlas_s3_frame.py',
           'src/recursive_horizons/nsc_atlas_parent_gh.py','src/recursive_horizons/nsc_atlas_trapped_surfaces.py',
           'scripts/check_nsc_atlas_valve_vortex.py']
    before={name:sha(ROOT/name) for name in names};loop_specs=[];samples=[]
    for core in range(2):
        for radius in (.01,.02):
            X,tangent=transverse_loop(core,radius);start=sum(len(value) for value in samples)
            samples.append(X);loop_specs.append((core,radius,start,start+len(X),tangent))
    all_points=np.vstack(samples+[CORE]);core_start=len(all_points)-2
    try:
        with threadpool_limits(limits=4):
            for degree in (16,20,24):
                if time.perf_counter()-started>budget:raise TimeoutError('180-second vortex-check budget exceeded')
                base=ROOT/f'.build/valve-e1-L{degree}-extended-v1'
                record_path=base.with_suffix('.json');cache=base.with_suffix('.npz');history=Path(str(base)+'-history.npz')
                record=json.loads(record_path.read_bytes());cache_sha=sha(cache)
                sparse=GHSurfaceGeometry.from_state_cache(cache,record,expected_sha256=cache_sha,snapshot=3,chunk=256)
                t=record['snapshots'][3]['coordinate_t']
                dense=GHSurfaceGeometry.from_history(history,record,coordinate_t=t,chunk=256)
                if not math.isclose(sparse.proper_duration,.9,rel_tol=0,abs_tol=1e-10):raise ValueError('wrong physical snapshot')
                modulus_reference=math.sqrt(sparse.system.parent.mass_squared)
                sparse_out=sparse(all_points);dense_out=dense(all_points)
                loops=[]
                for core,radius,lo,hi,tangent in loop_specs:
                    sZ=sparse_out['q'][lo:hi,1]+1j*sparse_out['Pi'][lo:hi,1]/modulus_reference
                    dZ=dense_out['q'][lo:hi,1]+1j*dense_out['Pi'][lo:hi,1]/modulus_reference
                    velocity=np.einsum('nij,nj->ni',frame(all_points[lo:hi]),tangent)/sparse.r
                    chi_current=-sparse_out['Pi'][lo:hi,1,None]*sparse_out['dq'][lo:hi,:,1]
                    circulation=float(2*math.pi*np.mean(np.einsum('ni,ni->n',chi_current,velocity)))
                    chi=sparse_out['q'][lo:hi,1];Pi=sparse_out['Pi'][lo:hi,1]
                    polygon=float(.5*np.sum(chi*np.roll(Pi,-1)-Pi*np.roll(chi,-1)))
                    loops.append(dict(core_number=core+1,reference_angular_radius=radius,points=128,
                        sparse=winding(sZ),dense=winding(dZ),
                        dense_sparse_max_abs_Z_gap=float(max(abs(dZ-sZ))),
                        chi_momentum_circulation=circulation,
                        quadrature_plane_polygon_area=polygon,
                        circulation_minus_polygon=circulation-polygon,
                        phase_plane_identity='integral j_chi = 1/2 integral(chi dPi-Pi dchi); finite polygon area is a separate sampled quadrature.'))
                cores=[]
                for number in range(2):
                    out={name:value[core_start+number] for name,value in sparse_out.items()}
                    cores.append(core_readout(sparse,number,out))
                for item in loops:
                    area=4*math.pi*math.sin(item['reference_angular_radius']/2)**2
                    curl=item['chi_momentum_circulation']/area
                    core_curl=cores[item['core_number']-1]['transverse_chi_Pi_jacobian']
                    item.update(reference_angular_disk_area=area,
                        circulation_over_reference_angular_disk_area=curl,
                        relative_gap_to_core_transverse_curl=curl/core_curl-1)
                angular_checks=[]
                if degree==24:
                    for count in (256,512):
                        if time.perf_counter()-started>budget:raise TimeoutError('angular check exhausted the vortex budget')
                        refinement=[];spec=[]
                        for core in range(2):
                            for radius in (.01,.02):
                                X,_=transverse_loop(core,radius,count=count)
                                lo=sum(len(v) for v in refinement);refinement.append(X);spec.append((core,radius,lo,lo+count))
                        fields=sparse(np.vstack(refinement))
                        for core,radius,lo,hi in spec:
                            Z=fields['q'][lo:hi,1]+1j*fields['Pi'][lo:hi,1]/modulus_reference
                            angular_checks.append(dict(core_number=core+1,reference_angular_radius=radius,
                                points=count,result=winding(Z)))
                surface_rows=[]
                if args.core_surfaces:
                    for core in range(2):
                        for radius in np.linspace(.2,2.6,13):
                            if time.perf_counter()-started>budget:raise TimeoutError('core surface search exhausted the vortex budget')
                            evaluation=evaluate_surface(round_surface(CORE[core],float(radius)),sparse,
                                frame_scale=sparse.r,quadrature_order=(16,33))
                            surface_rows.append(dict(core_number=core+1,reference_radius=float(radius),**evaluation.summary()))
                cases.append(dict(degree=degree,epsilon=record['epsilon'],coordinate_t=t,
                    observer_proper_duration=sparse.proper_duration,producer_commit=record['producer_commit'],
                    producer_source_sha256=record['source_sha256'],record_sha256=sha(record_path),
                    caller_admitted_sparse_cache_sha256=cache_sha,producer_pinned_dense_history_sha256=record['dense_history']['sha256'],
                    sparse_provenance=sparse.provenance,dense_provenance=dense.provenance,
                    dense_sparse_physical_state_max_gap=float(max(abs(sparse.state[:-1]-dense.state[:-1]))),
                    mu_reference=modulus_reference,loops=loops,cores=cores,angular_refinement=angular_checks,
                    core_surface_search=surface_rows,core_surface_scope='Direct sampled round S2 surfaces about the phase cores; finite radii and angular sampling do not prove absence of other trapped surfaces.'))
                print(f'Vortex degree{degree}: windings {[row["sparse"]["nearest_integer"] for row in loops]}',flush=True)
    except (ValueError,TimeoutError,RuntimeError) as error:failure=str(error)
    after={name:sha(ROOT/name) for name in names}
    result=dict(schema='nsc-valve-vortex-check-v1',status='complete' if failure is None and len(cases)==3 else 'failed',
        failure=failure,elapsed_seconds=time.perf_counter()-started,budget_seconds=budget,
        source_sha256=before,analysis_sources_unchanged=before==after,
        current_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),cases=cases,
        analytic=dict(two_mode_zero_condition='If D=q1*p2-q2*p1!=0, Z=chi+iPi/mu=0 iff Y1=Y2=0.',
            zero_circles=['X1=X2=0','X1=X4=0'],
            fixed_point_symmetry='Both chi and slice-normal Pi have proper-flip character(--), so they vanish at the fixed cores at all nonlinear orders; this does not alone select unit winding.',
            momentum_curl='j=-Pi_chi dchi-Pi_S dS, hence dj=dchi wedge dPi_chi+dS wedge dPi_S.',
            two_mode_curl='dj_chi=D dY1 wedge dY2.',
            index='For positive mu and the declared positive normal-coordinate loops at the two selected cores, the nondegenerate two-mode index is sign(D). Nonlinear index is checked from actual fields.'),
        scope='Resolved-field sampled phase-vortex/momentum-curl check on fixed proper-age GH slices. No PDE, fluid velocity, new horizon, child formation or net global rotation is inferred.')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=result['status'],failure=failure,elapsed=result['elapsed_seconds'],
        sources_unchanged=result['analysis_sources_unchanged']),indent=2))
    if result['status']!='complete' or not result['analysis_sources_unchanged']:raise SystemExit(1)


if __name__=='__main__':main()
