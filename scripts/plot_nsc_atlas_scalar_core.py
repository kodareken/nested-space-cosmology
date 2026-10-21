"""Scalar cooling and its independently computed self-gravitating profile."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,plots,ground):
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.5))
    x=np.linspace(0,10,400)
    axes[0].plot(x,ground.fields(x)[0],color='#122033',lw=2.8,label='Independent nodeless solution')
    for name,color,label in (('gaussian_fine','#176d8a','Gaussian initial cloud'),
                              ('supergaussian_fine','#9a4764','Super-Gaussian initial cloud')):
        data=plots[name];run=record['cooling_runs'][name]
        axes[0].plot(data['radius'],data['profile'],color=color,ls='--',lw=2,
                     label=f"Late core from {label.lower()}")
        history=data['snapshots'];times=[s['tau'] for s in history]
        axes[1].plot(times,[s['inner_mass']/run['initial_mass'] for s in history],color=color,lw=2,label=label)
    axes[0].set(title='Self-gravity selects a common core shape',xlabel='Radius after mass normalization',
        ylabel='Time-averaged envelope amplitude, normalized',xlim=(0,10))
    axes[0].legend(frameon=False,fontsize=8)
    axes[1].set(title='The core settles while scalar radiation escapes',xlabel='Envelope time',
        ylabel='Scalar norm inside radius 20 / initial norm',xlim=(0,600),ylim=(.7,1.01))
    axes[1].legend(frameon=False,fontsize=9,loc='upper right')
    for ax in axes:ax.grid(alpha=.18)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-scalar-core.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)


if __name__=='__main__':
    from check_nsc_atlas_scalar_core import calculation
    render(*calculation())
