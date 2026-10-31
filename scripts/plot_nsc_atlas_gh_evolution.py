"""Plot the preserved finite-parent controls; no raw state arrays are needed."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    comparison=json.loads((ROOT/'results/summaries/nsc-atlas-gh-comparison-v1.json').read_text(encoding='utf-8'))
    campaign=json.loads((ROOT/'results/summaries/nsc-atlas-gh-campaign-v1.json').read_text(encoding='utf-8'))
    cases=[case for case in comparison['cases'] if case['temporal_tolerances']['rtol']==3e-8]
    cases.sort(key=lambda case:case['degree'])
    degree=[case['degree'] for case in cases]
    finest=cases[-1]
    candidates=[case['run'] for case in campaign['cases'] if case['run']['epsilon']==comparison['epsilon']
        and case['run']['degree']==finest['degree']
        and case['run']['temporal_tolerances']==finest['temporal_tolerances']]
    if len(candidates)!=1:raise ValueError('one matching finest run is required')
    run=candidates[0]
    points=run['snapshots'];initial=points[0]['excitation_normal_energy_over_initial_round_vacuum_energy']
    fig,axes=plt.subplots(2,2,figsize=(11.6,7.8),layout='constrained')
    a,b,c,d=axes.ravel()
    a.plot(degree,[case['final']['observer_phase_wrapped'] for case in cases],'o-',color='#225b91')
    a.set(xlabel='Maximum harmonic degree',ylabel='Direct endpoint phase [rad]',
        title='A  Same geodesic, same proper endpoint')
    b.plot(degree,[case['final']['geodesic_observer_total_density_over_initial_vacuum'] for case in cases],
        'o-',color='#a85b20')
    b.set(xlabel='Maximum harmonic degree',ylabel=r'Local $T(u,u)/V_*$',
        title='B  Local density converges separately')
    for name,label,marker in [('ADM_H_RMS_over_H_squared','Hamiltonian','o'),
            ('ADM_M_RMS_over_H_squared','Momentum','s'),('GH_C_RMS_over_H','GH gauge','^')]:
        c.semilogy(degree,[case['final_errors'][name] for case in cases],marker+'-',label=label)
    c.set(xlabel='Maximum harmonic degree',ylabel='Normalized endpoint RMS',
        title='C  Independent gravitational constraints')
    c.legend(fontsize=9)
    proper=np.array([point['observer_proper_duration'] for point in points])
    energy=np.array([point['excitation_normal_energy_over_initial_round_vacuum_energy'] for point in points])/initial
    predicted=1+np.array([point['integrated_excitation_work_normalized'] for point in points])/initial
    d.plot(proper,energy,'o-',color='#225b91',label='ADM excitation energy')
    d.plot(proper,predicted,'x--',color='#a85b20',label='Initial energy + integrated work')
    d.set(xlabel='Observer proper duration [model units]',ylabel=r'$E_{exc}/E_{exc}(0)$',
        title='D  Independent energy/work calculations')
    relative=abs(points[-1]['excitation_work_balance_residual']/initial)
    d.text(.98,.04,f'Endpoint relative mismatch: {relative:.2e}',transform=d.transAxes,fontsize=9,ha='right')
    d.legend(fontsize=8,loc='upper right')
    for axis in axes.ravel():
        axis.grid(alpha=.2);axis.ticklabel_format(axis='x',style='plain',useOffset=False)
    for axis in (a,b):axis.set_xticks(degree)
    b.ticklabel_format(axis='y',style='plain',useOffset=False)
    fig.suptitle('One nonspherical finite preparation: geometry, local signals and work\n'
        r'$\epsilon=0.25$ control; empirical convergence, not a selected cosmological state',fontsize=12)
    out=ROOT/'paper/figures/nsc-atlas-gh-evolution.png'
    fig.savefig(out,dpi=180);plt.close(fig)
    print(out.relative_to(ROOT))


if __name__=='__main__':main()
