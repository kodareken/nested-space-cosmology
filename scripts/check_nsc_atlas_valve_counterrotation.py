"""One finite-amplitude counterrotating initial-data continuation and check."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from recursive_horizons.nsc_atlas_valve_counterrotation import CounterrotatingData,solve_counterrotating_initial_data

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def imported_sources():
    package=(ROOT/'src/recursive_horizons').resolve();paths={Path(__file__).resolve()}
    for module in list(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name:
            path=Path(name).resolve()
            if path.suffix=='.py' and path.is_relative_to(package):paths.add(path)
    return {path.relative_to(ROOT).as_posix():sha(path) for path in sorted(paths)}


def compare_prototype(data,path):
    if not path.exists():return dict(available=False,required_for_reproduction=False)
    old=json.loads(path.read_text());raw=ROOT/old['raw_cache']
    if sha(raw)!=old['raw_cache_sha256']:raise RuntimeError('private prototype cache changed')
    with np.load(raw,allow_pickle=False) as saved:
        if not np.array_equal(saved['labels'],np.asarray(data.labels)):
            raise RuntimeError('production/prototype basis differs')
        coefficient_difference=float(np.linalg.norm(data.coefficients-saved['coefficients']))
        tensor_difference=float(np.linalg.norm(data.tensor_coefficients-saved['York_tensor_coefficients']))
    new=data.metadata['independent_validation'];reference=old['validation']
    H_difference=abs(new['physical_Hamiltonian_RMS_over_H_squared']-reference['physical_Hamiltonian_RMS_over_H_squared'])
    energy_difference=abs(new['accounting']['total_proper_matter_energy']-reference['accounting']['total_proper_matter_energy'])
    passed=coefficient_difference<1e-10 and tensor_difference<1e-12 and H_difference<1e-8 and energy_difference<1e-7
    if not passed:raise RuntimeError('production counterrotation does not match its reviewed prototype')
    return dict(available=True,required_for_reproduction=False,record_sha256=sha(path),raw_cache_sha256=sha(raw),
                coefficient_L2_difference=coefficient_difference,tensor_L2_difference=tensor_difference,
                full_Hamiltonian_RMS_difference=H_difference,total_proper_matter_energy_difference=energy_difference,passed=passed)


def calculation(output,cache_directory,prototype):
    start=time.perf_counter();before=imported_sources()
    data=solve_counterrotating_initial_data()
    after=imported_sources()
    if before!=after:raise RuntimeError('counterrotation source changed during continuation')
    comparison=compare_prototype(data,prototype)
    directory=cache_directory.resolve()
    if not directory.is_relative_to(ROOT):raise ValueError('raw counterrotation cache must stay inside the repository')
    directory.mkdir(parents=True,exist_ok=True);raw=directory/'initial-data.npz'
    np.savez_compressed(raw,coefficients=data.coefficients,labels=np.asarray(data.labels),
        A=data.A,kick_amplitude=data.kick_amplitude,tensor_coefficients=data.tensor_coefficients,
        tensor_lrs=np.asarray([label[:3] for label in data.tensor_labels],dtype=int),
        tensor_kind=np.asarray([label[3] for label in data.tensor_labels]),
        vacuum=data.parent.vacuum,mass_squared=data.parent.mass_squared,offset=data.parent.offset,kappa=data.parent.kappa,
        metadata=json.dumps(dict(preparation_id='one-mode-counterrotation-v1',source_sha256=before,data_metadata=data.metadata)))
    provenance=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
        git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        role='Captured current working-tree solver sources; numeric state independently SHA-bound. No old source preset or historical result overwritten.')
    result=dict(schema='nsc-atlas-valve-counterrotation-v1',
        repository_root_relative_to_record=Path(os.path.relpath(ROOT,output.resolve().parent)).as_posix(),
        provenance=provenance,raw_cache=dict(path=raw.relative_to(ROOT).as_posix(),sha256=sha(raw)),
        preparation=dict(id='one-mode-counterrotation-v1',A=data.A,A_squared=data.A**2,kick_amplitude=data.kick_amplitude,
            chi='A*X1',momentum='Pi_chi=lambda*f*psi^-6; f=H^1.5 sqrt12/pi X2X4',
            S=0.,Pi_S=0.,trace_K=0.,freely_specified_conformal_TT=0.,
            K_reference='lambda*psi^-2*Abar; tensor(X) returns unit-lambda Abar in the original right frame.',
            status='Declared same-action benchmark preparation, not selected by inheritance or fitted to observations.'),
        metadata=data.metadata,prototype_comparison=comparison,
        verdict=dict(finite_amplitude_branch_completed=True,positive_finite_polynomial_bound=True,
            physical_constraints_are='Independent full residual estimates in metadata, not the tiny projected Newton residual.',
            evolution_computed=False,horizon_or_child_proved=False,SI_scale_selected=False,
            next='Use this explicit source/data API in a coupled formation test; legacy mixed-source initializers must not replacechi.'),
        compute=dict(elapsed_seconds=time.perf_counter()-start,BLAS_threads=6,raw_arrays_committed=False),
        reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_counterrotation.py --output results/summaries/nsc-atlas-valve-counterrotation-v1.json')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    # One real production check: strict authenticated loader plus unchanged
    # field/tensor/normal-momentum-density evaluations at original coordinates.
    loaded=CounterrotatingData.from_record(output)
    points=np.array([[1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.],
                     [0.,2**-.5,0.,2**-.5],[.5,.5,.5,.5]])
    checks=dict(authenticated_load=True,psi_max_difference=float(np.max(abs(data.values(points)-loaded.values(points)))),
                tensor_max_difference=float(np.max(abs(data.tensor(points)-loaded.tensor(points)))),
                f_max_difference=float(np.max(abs(data.f(points)-loaded.f(points)))),
                chi_max_difference=float(np.max(abs(data.chi(points)-loaded.chi(points)))))
    if any(v!=0 for key,v in checks.items() if key!='authenticated_load'):
        raise RuntimeError('authenticated state roundtrip changed the preparation')
    result['production_roundtrip_check']=checks
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-counterrotation-v1.json')
    parser.add_argument('--cache-directory',type=Path,default=ROOT/'.build/valve-counterrotation-v1')
    parser.add_argument('--compare-prototype',type=Path,default=ROOT/'.build/valve-counterrotating-continuation-j8.json')
    args=parser.parse_args();result=calculation(args.output,args.cache_directory,args.compare_prototype)
    print(json.dumps(dict(output=str(args.output),last_stage=result['metadata']['stages'][-1],
                         independent_constraints=result['metadata']['independent_validation'],
                         prototype=result['prototype_comparison'],roundtrip=result['production_roundtrip_check']),indent=2))


if __name__=='__main__':main()
