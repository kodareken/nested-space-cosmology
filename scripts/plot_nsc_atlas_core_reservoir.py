"""Where the compact core's geometric mass comes from."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,columns):
    r=record['runs']['compensated'];R,M,f,Pi,beta=columns.T
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.4))
    cuts=r['annulus_boundaries'];colors=['#e8ecf4','#d9eeea','#ffffff','#f5e7d5']
    for ax in axes:
        for a,b,color in zip(cuts[:-1],cuts[1:],colors):ax.axvspan(a,b,color=color)
        ax.set_xlim(0,cuts[-1]);ax.set_xlabel('Areal radius (atlas units)');ax.grid(alpha=.15)
    axes[0].plot(R,M*1e4,color='#173c60',lw=2)
    axes[0].axhline(0,color='#555555',lw=.7)
    axes[0].set(title='A core supplied by the inherited parent',ylabel=r'Enclosed geometric excess $M\times10^4$')
    axes[0].text(.03,.9,'core',transform=axes[0].transAxes)
    axes[0].text(.29,.9,'reservoir',transform=axes[0].transAxes)
    axes[0].text(.79,.9,'pulse',transform=axes[0].transAxes)
    # Separate the small radiating shell so its nonzero flux remains legible.
    F=np.zeros_like(R)
    from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir
    model=CoreReservoir();_,F,_=model.fields(R,r['depletion_fraction'])
    current=-Pi*beta*F
    axes[1].plot(R,current*1e4,color='#9b5328',lw=1.4)
    axes[1].set(title='A real outward scalar current',ylabel=r'Outward energy current $j\times10^4$')
    axes[1].text(.04,.89,'Same exterior mass\nSame initial walls\nPositive absolute energy',transform=axes[1].transAxes,fontsize=10)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-core-reservoir.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
