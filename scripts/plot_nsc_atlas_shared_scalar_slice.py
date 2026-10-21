"""Draw the compact spatial atlas and its causally protected parent cap."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice

ROOT=Path(__file__).resolve().parents[1]


def main():
    s=SharedScalarSlice();c=s.cap_geometry;l=s.lengths_and_volumes;d=s.parent_causal_domain()
    x1=np.linspace(0,c['inner_length'],400)
    a=np.sqrt(-c['inner_curvature']);r1=np.sinh(a*x1)/a
    r2=np.linspace(s.inner_radius,s.outer_radius,1001)
    x2=x1[-1]+cumulative_trapezoid(1/s.middle(r2)['beta'],r2,initial=0)
    x2[-1]=l['wall_positions'][1]
    rp=np.linspace(c['parent_length'],0,1001)
    x3=l['total_length']-rp
    r3=np.sin(np.sqrt(c['parent_curvature'])*rp)/np.sqrt(c['parent_curvature'])
    protected_start=l['total_length']-d['forever_protected_angle_lower_bound']*d['initial_scale_factor']
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.5),gridspec_kw={'width_ratios':[1.6,1]})
    colors=('#9a4764','#d28a2e','#176d8a')
    for ax in axes:
        for x,r,col,label in zip((x1,x2,x3),(r1,r2,r3),colors,('Grandchild','Shared middle','Parent')):
            ax.plot(x,r,color=col,lw=2.5,label=label)
        for x in l['wall_positions']:ax.axvline(x,color='#555555',lw=.8,ls=':')
        ax.set_xlabel('Proper distance along the initial slice')
        ax.set_ylabel('Areal radius R')
        ax.grid(alpha=.18)
    axes[0].axvspan(protected_start,l['total_length'],color='#176d8a',alpha=.10,label='Protected parent cap')
    axes[0].set_title('A compact atlas with two regular poles')
    axes[0].set_xlim(0,l['total_length']);axes[0].set_ylim(0,7.1)
    axes[0].text(11,2.3,'At least 90.7% of initial parent volume\nis outside both walls\' future influence',ha='center',fontsize=10)
    axes[0].legend(frameon=False,loc='upper right',fontsize=9)
    axes[1].set_xlim(0,2.3);axes[1].set_ylim(0,3.)
    axes[1].set_title('The two moving-wall joins, enlarged')
    axes[1].text(.38,.5,'q = 0.5',color=colors[0],fontsize=10)
    axes[1].text(1.02,1.0,'q = 1',color=colors[1],fontsize=10)
    axes[1].text(1.66,1.6,'q = 1.5',color=colors[2],fontsize=10)
    fig.tight_layout()
    output=ROOT/'docs/figures/nsc-atlas-shared-scalar-slice.png'
    fig.savefig(output,dpi=170);plt.close(fig)
    print(output)


if __name__=='__main__':main()
