"""Finite response separates density deformation from phase sensitivity."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,folders):
    h=json.loads((Path(folders['fine'])/'history.json').read_text(encoding='utf-8'))
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.3));blue='#176d8a';orange='#a56338'
    radial=h['radial'];tau=[x['tau'] for x in radial]
    axes[0,0].plot(tau,[x['field_norm_gain'] for x in radial],color=blue,label='Field-norm amplification')
    axes[0,0].plot(tau,[x['core_density_gain'] for x in radial],color=orange,label='Density response per initial field perturbation')
    axes[0,0].set(title='Radial field sensitivity is not density growth',xlabel='Slow time',ylabel='Declared-span response coefficient')
    axes[0,0].legend(frameon=False,fontsize=8)
    for name,label,color in (('quadrupole','Seeded at time 0',blue),('late_quadrupole','Seeded at time 480',orange)):
        rows=h[name];times=[x['tau'] for x in rows]
        axes[0,1].plot(times,[x['field_norm_gain'] for x in rows],color=color,label=label)
        axes[1,0].plot(times,[x['core_density_gain'] for x in rows],color=color,label=label)
    axes[0,1].axhline(1,color='#777777',ls=':',lw=1,label='Passive unitary reference')
    axes[0,1].set(title='Quadrupole field response stays modest',xlabel='Slow time',ylabel='Field-norm amplification')
    axes[0,1].legend(frameon=False,fontsize=8)
    axes[1,0].set(title='All quadrupole orientations share this bound',xlabel='Slow time',ylabel='Core density response coefficient')
    axes[1,0].legend(frameon=False,fontsize=8)
    checks=record['nonlinear_radial_derivative_checks'];epsilon=[r['epsilon'] for r in checks]
    for key,label,color in (('relative_field_derivative_error','Complex field',blue),
            ('relative_core_density_derivative_error','Core density',orange)):
        axes[1,1].loglog(epsilon,[r[key] for r in checks],'o-',color=color,label=label)
    axes[1,1].set(title='Independent nonlinear radial validation',xlabel='Initial fractional field perturbation',ylabel='Relative centred-difference error')
    axes[1,1].legend(frameon=False,fontsize=8)
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-hierarchy-perturbations.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
