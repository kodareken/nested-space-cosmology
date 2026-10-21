"""Same-source density growth and mode dependence in the protected parent."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from recursive_horizons.nsc_atlas_parent_perturbations import ParentPacket,jeans_comparison

ROOT=Path(__file__).resolve().parents[1]


def main():
    packet=ParentPacket();evolution,_,initial=packet.evolve(128)
    theta=np.linspace(0,1.7,601)
    end=packet.reconstruct(evolution,initial,evolution.solution.t[-1],theta)
    first=packet.profiles(theta)
    data=evolution.sample([evolution.solution.t[-1]])
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.5))
    axes[0].plot(theta,first['Delta'],color='#888888',lw=2,label='Initial constrained packet')
    axes[0].plot(theta,end['Delta'],color='#176d8a',lw=2.5,label='After one e-fold')
    axes[0].axhline(0,color='#777777',lw=.7)
    axes[0].set(xlabel='Comoving angle from parent pole',ylabel='Density contrast / initial central contrast',
        title='The same inherited field develops structure',xlim=(0,1.7))
    axes[0].text(.65,.8,'Central contrast: +15.19%\nFixed-metric control: -1.61%',fontsize=10)
    axes[0].legend(frameon=False,fontsize=9,loc='upper right')
    axes[1].plot(evolution.modes,data['Delta'][:,0,0],color='#176d8a',lw=2,label='Unit initial density, zero field velocity')
    axes[1].plot(evolution.modes,data['Delta'][:,1,0],color='#9a4764',ls='--',lw=1.7,label='Unit initial scalar phase, zero density')
    axes[1].axhline(1,color='#777777',lw=.7,ls=':')
    axes[1].axhline(0,color='#777777',lw=.7)
    kj=jeans_comparison()['initial_physical_Jeans_wavenumber']*packet.model.initial['a']
    nj=-1+np.sqrt(1+kj*kj)
    axes[1].axvline(nj,color='#d28a2e',ls=':',label='Initial averaged Jeans estimate')
    axes[1].set(xlabel='Closed-space harmonic index n',ylabel='Final density response coefficient',
        title='Growth depends on scale and preparation',xscale='log',xlim=(2,128))
    axes[1].legend(frameon=False,fontsize=8.3,loc='lower left')
    for ax in axes:ax.grid(alpha=.18)
    fig.tight_layout()
    path=ROOT/'docs/figures/nsc-atlas-parent-perturbations.png'
    fig.savefig(path,dpi=170);plt.close(fig);print(path)


if __name__=='__main__':main()
