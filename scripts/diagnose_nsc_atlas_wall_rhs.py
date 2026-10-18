"""Pointwise scalar-equation defect on one already archived P001 state."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits

import run_nsc_atlas_fixed_endpoint as endpoint
from diagnose_nsc_atlas_endpoint_radial import read_seed
from recursive_horizons.nsc_atlas_cap_pointwise_rhs import scalar_rhs,incoming_trace_defect
from recursive_horizons.nsc_atlas_cap_potential_pair import reference_potential_pair
from recursive_horizons import nsc_atlas_cap_fixed_reference as fixed
from recursive_horizons import nsc_atlas_cap_variational_reference as radial
from recursive_horizons.nsc_atlas_cap_panel_jets import history_piece_indices

ROOT=Path(__file__).resolve().parents[1]
REFERENCE='results/summaries/nsc-atlas-independent-cap-mass-lapse-v1.json'
CPU_LIMIT=60.
VALUE_LIMIT=40_000_000
BYTE_LIMIT=12000


def binding(path):
    raw=path.read_bytes()
    return dict(path=path.relative_to(ROOT).as_posix(),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())


def validate_recovery(expected,actual):
    comparisons={}
    for key,value in actual.items():
        target=float(expected[key]);scale=max(1.,abs(target),abs(value))
        gap=abs(value-target)
        comparisons[key]=dict(actual=value,reference=target,absolute_difference=gap)
        if not math.isfinite(value) or gap>1e-10*scale:
            raise ValueError('archived stage recovery differs: '+key)
    return comparisons


def projection_check(basis,state,cap,background,solutions,budget):
    """Independent strong-PDE load on the same captured radial metric."""
    Z=Fraction(cap['wall_coordinate']);T=cap['background']['time'];tx=cap['T_x']
    fields=fixed._Fields(basis,state['cap']['chi'],state['cap']['canonical_density'],budget)
    breaks,owners=fixed._background_partition(basis,Z,T,background)
    if len(solutions)!=len(owners):raise ValueError('captured radial panel count changed')
    q,w=leggauss(32);loads=np.zeros((2,basis.dimension));scales=np.zeros_like(loads)
    cT=cap['rates']['chi']/tx;dT=cap['rates']['canonical_density']/tx
    wall=cap['wall']
    for panel,(cell,solution) in enumerate(zip(owners,solutions)):
        lo,hi=breaks[panel:panel+2];half=(hi-lo)/2
        if solution.t[0]!=lo or solution.t[-1]!=hi:raise ValueError('captured radial span changed')
        points=lo+half*(q+1);metric=solution.sol(points)
        midpoint=lo+(hi-lo)/2
        zmid=fields.left[cell]+fields.width[cell]*(midpoint-cell)
        pieces=history_piece_indices(background.profile,[T-zmid,T+zmid],budget=budget)
        for j,point in enumerate(points):
            xi=float(point-cell);row=fields.evaluate(int(cell),xi,budget)
            v=basis.values(int(cell),xi,budget=budget);idx=v['indices'];z=row['z'];r=row['r']
            bg=background.reduced_jets(T,z,budget=budget,polynomial_indices=pieces)
            chi=state['cap']['chi'][idx];density=state['cap']['canonical_density'][idx]
            jets=dict(chi=row['chi']+bg['chi'],chi_z=row['chi_z']+bg['chi_z'],
                chi_zz=float(v['second']@chi)+bg['k_T'],k=row['k']+bg['k'],
                k_z=float(v['first']@density)+bg['k_z'])
            P,U=reference_potential_pair(jets['chi']/r,cap['alpha'],budget=budget)
            rhs=scalar_rhs(radius=r,**jets,mass=float(metric[0,j]),lapse=float(metric[1,j]),
                alpha=cap['alpha'],potential_difference=P,potential_derivative=U)
            actual=np.array([float(v['value']@cT[idx])+bg['k'],float(v['value']@dT[idx])+bg['k_T']])
            reference=np.array([rhs['chi_T'],rhs['k_T']])
            weights=half*fields.width[cell]*w[j]*v['value']
            loads[:,idx]+=(actual-reference)[:,None]*weights
            scales[:,idx]+=(abs(actual)+abs(reference))[:,None]*abs(weights)
            budget.check(64*basis.dimension)
    expected=np.zeros_like(loads)
    gram=basis.gram(Z,budget=budget);idx=gram.wall['indices'];e=gram.wall['value']
    expected[1,idx]=e*wall['B']*(wall['h_weak']-wall['h_trace'])
    differences=loads-expected
    scale=np.maximum(1.,np.maximum(np.max(abs(loads),axis=1),np.max(abs(expected),axis=1)))
    ratios=np.max(abs(differences),axis=1)/scale
    if np.max(ratios)>1e-8:raise ValueError('independent projected scalar equation differs from the natural boundary load')
    return dict(quadrature_order=32,panels=len(owners),points=32*len(owners),
        velocity_load=loads[0].tolist(),density_load=loads[1].tolist(),expected_density_boundary_load=expected[1].tolist(),
        maximum_absolute_load_differences=np.max(abs(differences),axis=1).tolist(),normalized_load_differences=ratios.tolist(),
        load_scales=scale.tolist(),integrated_absolute_operation_scales=np.max(scales,axis=1).tolist(),
        tolerance=1e-8,scope='Strong-PDE projection on the existing captured metric. Not a new radial solve or a physical error enclosure.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    if not args.output.parent.is_dir():raise ValueError('output parent must exist')
    commit,sources=endpoint.frozen_sources()
    test='tests/test_nsc_atlas_cap_pointwise_rhs.py'
    raw_test=(ROOT/test).read_bytes()
    if subprocess.check_output(['git','show',f'{commit}:{test}'],cwd=ROOT)!=raw_test:
        raise ValueError('freeze the equation controls before the calculation')
    sources[test]=hashlib.sha256(raw_test).hexdigest()
    budget=endpoint.Budget(maximum_evaluations=VALUE_LIMIT,cpu_seconds=CPU_LIMIT)
    reference=json.loads((ROOT/REFERENCE).read_bytes())
    old=reference['independent_mass_lapse']
    if old['outcome']!='passed':raise ValueError('the independent radial equation check must have passed')
    receipt_path=ROOT/reference['input_receipt']['path']
    if binding(receipt_path)!=reference['input_receipt']:raise ValueError('the archived receipt changed')
    receipt=json.loads(receipt_path.read_bytes());state,x=read_seed(receipt)
    if x!=.4370111883201181 or receipt['level']!=0:raise ValueError('the admitted archived point changed')
    report=dict(schema='nsc-atlas-wall-pointwise-rhs-v1',status='incomplete',source_context_commit=commit,
        source_graph_sha256=hashlib.sha256(json.dumps(sources,sort_keys=True).encode()).hexdigest(),
        source_files=len(sources),source_sha256={k:v for k,v in sources.items() if k in (
            'scripts/diagnose_nsc_atlas_wall_rhs.py','src/recursive_horizons/nsc_atlas_cap_pointwise_rhs.py',test)},
        independent_mass_check=binding(ROOT/REFERENCE),input_receipt=binding(receipt_path),stage_seed=receipt['failed_stage_seed'],
        x=x,level=0,scope='One archived coupled state. Local fixed-r scalar-equation defects; not the complete transmission derivative, an evolution or a closure theorem.',
        limits=dict(CPU_seconds=CPU_LIMIT,charged_values=VALUE_LIMIT,summary_bytes=BYTE_LIMIT),
        counts=dict(new_meshes=0,time_steps=0,coupled_stages=0,original_history_preparations=0))
    try:
        with threadpool_limits(limits=1):
            endpoint.authenticate(ROOT,budget)
            source_raw=(ROOT/endpoint.producer.BASE_INPUT).read_bytes()
            if hashlib.sha256(source_raw).hexdigest()!=receipt['input']['sha256']:
                raise ValueError('the physical source changed')
            base=json.loads(source_raw);geometry=base['geometry'];data=base['initial_data']['finite']
            rows={k:np.asarray(v,dtype=float) for k,v in data['exact_rows'].items()}
            spec=endpoint.level_spec(0);basis=endpoint.FixedCapBasis(endpoint.cap_knots(0),budget=budget)
            preparation={};report['counts']['original_history_preparations']=1
            background,offset=endpoint.prepare_history(budget,preparation,enrich_ordinary_pieces=True)
            for key in ('input_sha256','original_input_sha256'):
                if preparation[key]!=receipt['preparation'][key]:raise ValueError('the original history changed')
            report['prepared_input_sha256']=preparation['input_sha256']
            snapshot=endpoint.SourceSnapshot(geometry,alpha=data['alpha'],x_nodes=rows['x'],source_rows=rows,budget=budget)
            gauge=endpoint.ReferenceWallGauge(geometry,order=64,budget=budget)
            report['counts']['coupled_stages']=1
            original=radial.solve_ivp;solutions=[]
            def capture(*a,**kw):
                result=original(*a,**kw);solutions.append(result);return result
            radial.solve_ivp=capture
            try:
                stage=endpoint.coupled_fixed_stage(state,basis=basis,x=x,
                    eta=np.linspace(0.,1.,spec['exterior_cells']+1),source=snapshot.sample(x),gauge=gauge.sample(x),
                    alpha=data['alpha'],tube=endpoint.TUBE,budget=budget,background=background,background_clock_offset=offset,
                    rtol=2.5e-14,atol=1e-15,background_panel_sides=True,compiled_coefficients=False)
            finally:
                radial.solve_ivp=original
            cap=stage['connection']['cap'];wall=cap['wall'];R=cap['radius'];tx=cap['T_x']
            expected=dict(old['reference'],p=old['wall']['p'],h=old['wall']['h'])
            recovered=dict(full_mass=R**3/2+data['alpha']*cap['mass'],lapse=cap['lapse'],
                temporal_mass_residual=stage['connection']['residual']['physical_temporal_mass_actual'],p=wall['p_trace'],h=wall['h_trace'])
            report['saved_stage_recovery']=validate_recovery(expected,recovered)
            Z=Fraction(cap['wall_coordinate']);cell=basis.cells-1
            xi=float((Z-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2]))
            v=basis.values(cell,xi,budget=budget);idx=v['indices']
            bg=background.reduced_jets(cap['background']['time'],float(Z),budget=budget,
                polynomial_indices=cap['background']['wall_polynomial_indices'])
            chi=state['cap']['chi'][idx];density=state['cap']['canonical_density'][idx]
            cT=cap['rates']['chi'][idx]/tx;dT=cap['rates']['canonical_density'][idx]/tx
            jets=dict(chi=float(v['value']@chi)+bg['chi'],chi_z=float(v['first']@chi)+bg['chi_z'],
                chi_zz=float(v['second']@chi)+bg['k_T'],k=float(v['value']@density)+bg['k'],
                k_z=float(v['first']@density)+bg['k_z'])
            P,U=reference_potential_pair(jets['chi']/R,data['alpha'],budget=budget)
            continuum=scalar_rhs(radius=R,**jets,mass=cap['mass'],lapse=cap['lapse'],alpha=data['alpha'],
                potential_difference=P,potential_derivative=U)
            represented=dict(chi_T=float(v['value']@cT)+bg['k'],k_T=float(v['value']@dT)+bg['k_T'],
                chi_T_z=float(v['first']@cT)+bg['k_z'])
            defects=incoming_trace_defect(radius=R,continuum=continuum,
                **{'represented_'+k:v for k,v in represented.items()})
            identity=defects['chi_equation_defect']-wall['velocity_discrepancy']
            if abs(identity)>1e-10*max(1.,abs(wall['velocity_discrepancy'])):
                raise ValueError('the local kinematic defect does not recover the existing wall channel')
            report['projected_field_equations']=projection_check(basis,state,cap,background,solutions,budget)
            report.update(status='completed_archived_pointwise_rhs',point=dict(radius=R,wall_coordinate=float(Z),
                history_time=cap['background']['time'],T_x=tx,R_x=cap['radius_rate'],alpha=data['alpha']),
                continuum=continuum,represented=represented,defects=defects,kinematic_identity_defect=identity,
                wall_channels=dict(momentum_discrepancy=wall['momentum_discrepancy'],normal_discrepancy=wall['normal_discrepancy'],
                    incoming_discrepancy=wall['incoming_discrepancy'],conformal_flux=stage['connection']['residual']['conformal_wall_flux']),
                coupled_x_cap_defect=tx*defects['incoming_cap_rate_defect'])
            budget.check()
    except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
        report.update(status='stopped',reason=str(error))
    report['runtime_cpu_seconds']=time.process_time()-budget.started;report['charged_values']=budget.evaluations
    raw=json.dumps(endpoint.plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>BYTE_LIMIT:raise RuntimeError('pointwise RHS summary exceeds byte limit')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps({k:report[k] for k in ('status','runtime_cpu_seconds','charged_values')}))
    if report['status']!='completed_archived_pointwise_rhs':raise RuntimeError(report['reason'])


if __name__=='__main__':main()
