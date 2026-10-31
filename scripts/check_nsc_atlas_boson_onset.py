"""Independent local continuum onset from the inherited initial slice."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import PAIRS,ParentGH
from recursive_horizons.nsc_atlas_s3_frame import point_values
from record_nsc_atlas_gh_campaign import source_binding

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--history',type=Path,required=True)
    parser.add_argument('--background-record',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();output=(ROOT/args.output).resolve()
    if not output.is_relative_to(ROOT) or output.exists():raise ValueError('fresh repository output required')
    background_path=(ROOT/args.background_record).resolve();path=(ROOT/args.history).resolve()
    background=json.loads(background_path.read_bytes())
    if not source_binding(background) or hashlib.sha256(path.read_bytes()).hexdigest()!=background['dense_history']['sha256']:
        raise ValueError('verified dense background required')
    history=DenseGHHistory.load(path);m=history.metadata
    for name,expected in m['source_sha256'].items():
        if name.startswith('src/') and hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('current model differs from the inherited background')
    with threadpool_limits(limits=1):
        system=ParentGH(m['degree'],m['epsilon'],gauge=m['gauge'],backend=m['backend'],
            quadrature_order=m.get('quadrature_order'))
        state=history.evaluate(history.bounds[0][0]);C=system.unpack(state)
        X=state[system.dof:system.dof+4];X=X/np.linalg.norm(X)
        row=point_values(X,system.basis.labels)
        values=C@row
        g=np.empty((4,4))
        for k,(i,j) in enumerate(PAIRS):g[i,j]=g[j,i]=values[k]
        determinant=float(np.linalg.det(g[1:,1:]));lapse=1/math.sqrt(-np.linalg.inv(g)[0,0])
        if determinant<=0:raise ValueError('invalid initial spatial metric')
        s=system.sstar+values[20]/math.sqrt(system.K);chi=values[21]
        mixing=-12*400*math.exp(-12*s)*chi/math.sqrt(system.K)
        profiles=(math.sqrt(2)/math.pi*X[0],math.sqrt(12)/math.pi*X[1]*X[3])
        answers=[]
        for number,profile in enumerate(profiles,1):
            physical_profile=system.r**1.5*profile
            volume_coefficient=-mixing*physical_profile/(6*math.sqrt(determinant))
            log_coefficient=-7/math.sqrt(system.K)*volume_coefficient
            answers.append(dict(source_harmonic=number,physical_reference_harmonic_at_observer=physical_profile,
                local_volume_proper_cubic_coefficient=volume_coefficient,
                local_log_g_KK_proper_cubic_coefficient=log_coefficient,
                local_log_g_KK_coordinate_cubic_coefficient=lapse**3*log_coefficient))
    result=dict(schema='nsc-atlas-boson-onset-v1',background_history_sha256=background['dense_history']['sha256'],
        background_producer_commit=background['producer_commit'],background_degree=m['degree'],epsilon=m['epsilon'],
        initial_lapse=lapse,initial_volume_weight=math.sqrt(determinant),initial_mixing_mass_squared=mixing,
        observer_quaternion=X.tolist(),response=answers,
        derivation='Initially deltaS=normal(deltaS)=deltachi=0 and normal(deltachi)=f/sqrt(h). The differentiated field equation gives deltaS=-W_Schi*f*tau^3/(6sqrt(h))+O(tau^4). No projected Hamiltonian or chosen covariance is used.',
        scope='Local proper-age onset for a global unit canonical momentum impulse and fixed background geodesic. Background discretization remains; no primordial state, entropy, charge bias or SI scale is selected.',
        source_sha256={'scripts/check_nsc_atlas_boson_onset.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
