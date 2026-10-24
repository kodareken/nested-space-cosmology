"""The core shrinks faster than its conserved supply region."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record):
    rows=[r for name,r in record['runs'].items() if name.startswith('frequency_')]
    w=np.array([r['frequency'] for r in rows]);rc=np.array([r['core']['cutoff_radius'] for r in rows])
    rd=np.array([r['reservoir_outer_radius'] for r in rows])
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(12,4.2))
    axes[0].loglog(w,rc,'o-',color='#176d8a',label='Core cutoff')
    axes[0].loglog(w,rd,'s-',color='#a56338',label='Solved supply radius')
    axes[0].loglog(w,rc[0]*np.ones_like(w)*6.62607707610/rc[0],':',color='#777777',label='Parent reference radius')
    axes[0].set(title='One parent, two local scales',xlabel='Inherited scalar frequency',ylabel='Radius in atlas units')
    axes[0].legend(frameon=False,fontsize=8,loc='lower left')
    axes[1].loglog(w,rd/rc,'o-',color='#176d8a')
    axes[1].axhline(2,color='#777777',ls=':',label='Fixed radius factor 2')
    axes[1].axvline(record['fixed_ratio_obstruction']['necessary_frequency_ceiling'],color='#a56338',ls='--',
        label='Ideal supply ceiling at factor 2')
    axes[1].set(title='A fixed proportion cannot supply the core',xlabel='Inherited scalar frequency',ylabel='Supply radius / core cutoff')
    axes[1].legend(frameon=False,fontsize=7,loc='upper left')
    ratio=[r['Lambda_reservoir_radius_cubed_over_core_mass'] for r in rows]
    axes[2].semilogx(w,ratio,'o-',color='#176d8a',label='Exact finite constraints')
    axes[2].axhline(rows[0]['asymptotic_Lambda_reservoir_radius_cubed_over_core_mass'],color='#a56338',ls='--',
        label='Derived limit 18.330588')
    axes[2].set(title='Vacuum, mass and supply radius stay linked',xlabel='Inherited scalar frequency',
        ylabel=r'$\Lambda_v R_d^3/\mathfrak{M}_c$')
    axes[2].legend(frameon=False,fontsize=8)
    for ax in axes:ax.grid(alpha=.15)
    fig.tight_layout();path=ROOT/'docs/figures/nsc-atlas-scale-hierarchy.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)
