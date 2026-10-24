"""One receiving worldline: image pair, delay and the inherited vacuum scale."""
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record, raw, reduced):
    plt.rcParams.update({'font.size':11, 'axes.spines.top':False, 'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.6))
    colors={1:'#176d8a', -1:'#a56338'}
    images=record['central_images']; pair=record['central_pair']
    for im in images:
        branch=im['branch'];color=colors[branch]
        rows=[r for r in raw['fine']['rays'] if r['emission_proper_age']==2.8 and r['branch']==branch]
        axes[0,0].plot([math.degrees(r['sky_angle']) for r in rows],
            [1e3*(r['landing_angle']-math.pi) for r in rows],'o-',color=color,ms=4,
            label='Positive orientation' if branch==1 else 'Reflected orientation')
        axes[0,0].axhline(1e3*(im['target_unwrapped_angle']-math.pi),color=color,ls=':',lw=1)
        axes[0,0].plot(math.degrees(im['launch_angle']),1e3*(im['target_unwrapped_angle']-math.pi),
            'x',color=color,ms=9,mew=2)
        angle=math.degrees(im['received_sky_angle'])
        axes[0,1].scatter([angle],[0],s=230*im['fluence_weight']/images[0]['fluence_weight'],
            color=color,zorder=3)
        axes[0,1].annotate(f"{'Positive' if branch==1 else 'Negative'} parity\n{angle:+.6f} deg",
            (angle,0),xytext=(0,28),textcoords='offset points',ha='center',fontsize=10)
    axes[0,0].set(title='Two roots identify the same receiver',xlabel='Inward launch angle (degrees)',
        ylabel=r'Unwrapped landing angle minus $\pi$ (mrad)')
    axes[0,0].legend(frameon=False,fontsize=8)
    axes[0,1].axhline(0,color='#cccccc',lw=1)
    axes[0,1].axvline(0,color='#cccccc',lw=.8,ls=':')
    axes[0,1].set(title='Two resolved images of one pulse',xlabel='Signed angle in the receiving sky (degrees)',
        ylim=(-.5,.65),xlim=(-.48,.48),yticks=[])
    axes[0,1].text(.5,.13,f"Separation {math.degrees(pair['image_separation']):.6f} deg\n"
        f"Pulse fluence ratio {pair['fluence_ratio_negative_over_positive']:.5f}",
        transform=axes[0,1].transAxes,ha='center',fontsize=10)
    names=('coarse','medium','fine');x=np.array([1,2,4])
    measured=[1e6*record['independent_checks'][name]['delay_drift']['delay_derivative'] for name in names]
    predicted=[1e6*record['independent_checks'][name]['delay_drift']['predicted'] for name in names]
    axes[1,0].plot(x,measured,'o-',color='#176d8a',label='Arrival-time derivative')
    axes[1,0].plot(x,predicted,'s--',color='#a56338',label='Inverse-frequency difference')
    axes[1,0].set(title='One clock identity checks two calculations',xlabel='Resolution multiplier',
        ylabel=r'Delay change per source proper time ($10^{-6}$)',xticks=x)
    axes[1,0].legend(frameon=False,fontsize=8)
    delay=[1e6*record['refinement'][name]['delay_negative_minus_positive'] for name in names]
    axes[1,1].plot(x,delay,'o-',color='#176d8a')
    axes[1,1].set(title='Image delay and parent vacuum share a scale',xlabel='Resolution multiplier',
        ylabel=r'Received pulse delay ($10^{-6}$ atlas units)',xticks=x)
    conn=record['parent_vacuum_connection']
    axes[1,1].text(.04,.12,r'$\Lambda_{\rm vac}(c\Delta t)^2$'+f" = {conn['Lambda_c_delay_squared']:.6g}\n"
        'Conditional on the inherited parent and supplied cloud',transform=axes[1,1].transAxes,
        fontsize=9,bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
    for ax in axes.flat:ax.grid(alpha=.14)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-core-images.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
