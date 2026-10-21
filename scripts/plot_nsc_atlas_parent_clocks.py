"""Free-fall clock protocol and held-out relation in the protected parent."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from recursive_horizons.nsc_atlas_parent_clocks import ParentClocks

ROOT=Path(__file__).resolve().parents[1]


def main():
    m=ParentClocks();positions=(.04,.08,.12);rows=[m.response(s) for s in positions]
    R=np.array([r['log_frequency_response'] for r in rows]);q=np.array([1e-5,-4e-6]);y=R@q
    inferred=np.linalg.solve(R[:2],y[:2]);prediction=float(R[2]@inferred)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.7,4.5));colors=('#176d8a','#d28a2e','#9a4764')
    eta_o=m.background(m.to)['eta']
    for row,color in zip(rows,colors):
        t=np.linspace(row['background_emission_time'],m.to,180)
        theta=np.array([eta_o-m.background(z)['eta'] for z in t])
        if row is rows[-1]:axes[0].plot(theta,t,color='#7e8a97',lw=2,zorder=0)
        axes[0].scatter([row['initial_angle']],[row['background_emission_time']],color=color,s=25,
            label=f"Source at {row['initial_angle']:.2f} Lref initially")
        axes[0].axvline(row['initial_angle'],ls=':',color=color,lw=.8)
    axes[0].scatter([0],[m.to],marker='*',s=90,color='#122033')
    axes[0].set(xlabel='Unperturbed comoving radial angle',ylabel='Background proper time',
        title='Three signals, one receiver proper age',xlim=(-.004,.13))
    axes[0].legend(frameon=False,fontsize=8,loc='lower left')
    axes[1].scatter([1,2],y[:2]*1e7,s=70,color='#176d8a',label='Two synthetic calibration clocks')
    axes[1].scatter([3],[y[2]*1e7],s=100,marker='o',facecolors='none',edgecolors='#9a4764',lw=2,label='Withheld synthetic clock')
    axes[1].scatter([3],[prediction*1e7],s=65,marker='x',color='#122033',label='Predicted from the first two')
    axes[1].text(1.05,3.65,r'$y_3=-2.889905\,y_1+2.942327\,y_2$',fontsize=11)
    axes[1].set(xlabel='Emitter index',ylabel='Log frequency departure, in units of 1e-7',
        title='Eliminating the two unknown amplitudes',xticks=[1,2,3],xlim=(.7,3.3),ylim=(1.5,4.2))
    axes[1].legend(frameon=False,fontsize=8,loc='lower right')
    for ax in axes:ax.grid(alpha=.18)
    fig.tight_layout();path=ROOT/'docs/figures/nsc-atlas-parent-clocks.png'
    fig.savefig(path,dpi=170);plt.close(fig);print(path)


if __name__=='__main__':main()
