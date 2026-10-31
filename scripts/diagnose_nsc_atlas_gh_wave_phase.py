"""Diagnose stored resolved GH fields without repeating their evolution."""
import argparse
import json
from pathlib import Path

import numpy as np
from numba import set_num_threads
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_gh import ParentGH
from recursive_horizons.nsc_atlas_gh_wave_diagnostics import wave_phase_diagnostics

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--state-cache',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();source=(ROOT/a.state_cache).resolve();out=(ROOT/a.output).resolve()
    if not source.is_relative_to(ROOT/'.build') or not out.is_relative_to(ROOT) or out.exists():raise ValueError('local state cache and fresh output required')
    data=np.load(source,allow_pickle=False);meta=json.loads(str(data['metadata']))
    set_num_threads(2)
    with threadpool_limits(limits=2):
        system=ParentGH(meta['degree'],meta['epsilon'],gauge=meta['gauge'],backend=meta.get('backend','dense'))
        rows=[]
        for t,y in zip(data['times'],data['states']):
            J=system.jets(y)
            diag=wave_phase_diagnostics(J.g,J.dg,J.q,J.p,J.dq,J.d2q,J.dp,omega_ref=system.parent.reference_frequency)
            weights=system.basis.weights;N=system.N
            defined=diag['defined'][:N]
            full=diag['reference_gradient_squared'][:N]
            phase=diag['phase_gradient_contribution'][:N]
            amplitude=diag['amplitude_gradient_squared'][:N]
            integral=float(weights@full);usable=float(weights@(full*defined))
            power=float(weights@diag['amplitude_squared'][:N])
            rows.append(dict(coordinate_t=float(t),observer_proper_duration=float(y[system.dof+8]),
                phase_fraction_of_defined_quadrature_gradient=float(weights@phase/usable) if usable else None,
                amplitude_fraction_of_defined_quadrature_gradient=float(weights@amplitude/usable) if usable else None,
                undefined_gradient_weight_fraction=1-usable/integral if integral else 0.,
                weighted_reference_phase_wavenumber_over_r=(float(weights@phase/power)**.5/system.r) if power else None,
                masked_identity_RMS=float(np.sqrt(weights@(diag['identity_residual'][:N]**2)/sum(weights)))))
    r=dict(schema='nsc-atlas-gh-wave-phase-diagnosis-v1',degree=meta['degree'],epsilon=meta['epsilon'],
        scope='Diagnosis of supplied resolved numerical states; not validation of their Einstein constraints, conserved energy fractions, or baryon charge.',
        source_run_hashes=meta['source_sha256'],snapshots=rows)
    out.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(r,indent=2))


if __name__=='__main__':main()
