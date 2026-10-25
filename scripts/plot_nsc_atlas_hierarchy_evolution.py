"""Localized mass, expansion, conservation and the full-equation comparison."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,folders):
    history={name:json.loads((Path(folders[name])/'history.json').read_text(encoding='utf-8'))
        for name in ('fine','no_gravity')}
    fine=history['fine'];free=history['no_gravity']
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.3))
    axes[0,0].semilogy([x['tau'] for x in fine],[x['core_mass_ratio'] for x in fine],color='#176d8a',label='Self-gravity retained')
    axes[0,0].semilogy([x['tau'] for x in free],[x['core_mass_ratio'] for x in free],color='#a56338',label='Same field, self-gravity removed')
    axes[0,0].set(title='The supplied core stays localized',xlabel='Slow time',ylabel='Mass inside fixed physical cutoff / initial mass')
    axes[0,0].legend(frameon=False,fontsize=8,loc='lower left')
    inset=axes[0,0].inset_axes([.56,.52,.4,.26])
    inset.plot([x['tau'] for x in fine],[x['core_mass_ratio'] for x in fine],color='#176d8a',lw=1)
    inset.axhline(1,color='#777777',ls=':',lw=.7);inset.set_title('Small inward transfer',fontsize=8)
    inset.tick_params(labelsize=7)
    initial=record['preparation']['core']['half_mass_radius']*1e6
    axes[0,1].plot([x['tau'] for x in fine],[x['physical_half_mass_radius'] for x in fine],color='#176d8a',label='Evolved core half-mass radius')
    axes[0,1].plot([x['tau'] for x in fine],[initial*x['scale_factor'] for x in fine],color='#a56338',ls='--',label='Initial radius multiplied by parent expansion')
    axes[0,1].set(title='Core size separates from parent expansion',xlabel='Slow time',ylabel=r'Physical half-mass radius $\omega\zeta R_{1/2}$')
    axes[0,1].legend(frameon=False,fontsize=8)
    labels=('coarse','medium','fine');x=np.array([1,2,4])
    for key,label,color in (('maximum_energy_work_gap_over_initial_energy','Energy plus expansion work','#a56338'),
            ('maximum_core_mass_flux_gaps_over_initial_core','Fixed-radius mass flux','#176d8a')):
        values=[max(record['runs'][name][key]) if isinstance(record['runs'][name][key],list) else record['runs'][name][key] for name in labels]
        axes[1,0].semilogy(x,values,'o-',color=color,label=label)
    axes[1,0].set(title='Two conservation checks improve',xlabel='Resolution multiplier',ylabel='Relative discrepancy',xticks=x)
    axes[1,0].legend(frameon=False,fontsize=8)
    comp=record['short_full_equation_comparison'];xx=[comp[k]['zeta']**2 for k in ('small','large')]
    yy=[abs(comp[k]['difference']) for k in ('small','large')]
    axes[1,1].loglog(xx,yy,'o-',color='#176d8a',label='Full EKG minus envelope')
    axes[1,1].axhline(abs(comp['small_full_refinement_gap']),color='#a56338',ls='--',label='Full-EKG numerical refinement gap')
    axes[1,1].set(title='Short full-equation check: slow time 5.12',xlabel=r'Weak-field parameter $\zeta^2$',ylabel='Core-mass-ratio difference')
    axes[1,1].legend(frameon=False,fontsize=8,loc='center left')
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-hierarchy-evolution.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
