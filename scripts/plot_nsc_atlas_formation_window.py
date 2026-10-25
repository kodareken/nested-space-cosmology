"""Formation budget instead of another preinstalled-core evolution."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def render(r):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(12,3.8))
    names=('dust_comparison','compensation_turnover','Jeans_scale','core_half_radius')
    labels=('Positive comparison','Supply turnover','Initial Jeans scale','Core half-radius')
    colors=('#777777','#176d8a','#a56338','#926299')
    for name,label,color in zip(names,labels,colors):
        c=r['curves'][name];u=np.array(c['reciprocal_scale'])
        ax[0].plot(1-u,c['density'],label=label,color=color,lw=1.4)
    ax[0].axhline(r['elementary_all_wavelength_bound']['density_gain_bound'],color='#555555',ls=':',lw=1)
    ax[0].set(xlabel=r'Future progress $1-1/A$',ylabel=r'$\delta/\delta_i$',title='A finite late-time growth budget')
    ax[0].legend(frameon=False,fontsize=7,loc='lower left')
    selected=('earlier_A_0p001','earlier_A_0p01','earlier_A_0p1','compensation_turnover')
    Ai=[r['runs'][n]['initial_scale'] for n in selected]
    gains=[r['runs'][n]['future_density_response'][0] for n in selected]
    ax[1].semilogx(Ai,gains,'o',color='#176d8a',label='Density seed, zero current')
    ax[1].axhline(0,color='#777777',lw=.7)
    ax[1].set(xlabel='Supplied starting scale A',ylabel='Signed final linear gain',title='Earlier preparation is different data')
    ax[1].legend(frameon=False,fontsize=7)
    e=r['compact_seed']['runs']['time_control'];c=e['curve'];u=np.array(c['reciprocal_scale']);valid=u>0
    ax[2].semilogx(1/u[valid],np.array(c['central_linear_contrast'])[valid],color='#176d8a')
    ax[2].axhline(.1,color='#777777',ls=':',label='Small-contrast checkpoint')
    ax[2].axhline(1,color='#a56338',ls='--',label='Formal nonlinear threshold')
    ax[2].axvline(e['scale_at_formal_linear_delta_1'],color='#a56338',ls=':',lw=.8)
    ax[2].set(xlabel='Scale factor A',ylabel='Central linear density contrast',title='A compact 5% earlier seed')
    ax[2].legend(frameon=False,fontsize=7,loc='upper left')
    for a in ax:a.grid(alpha=.15)
    fig.tight_layout();path=Path(__file__).resolve().parents[1]/'docs/figures/nsc-atlas-formation-window.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)
