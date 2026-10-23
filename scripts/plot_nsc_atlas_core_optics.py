"""Two epochs of the same evolving scalar lens, with independent checks."""
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,raw,homogeneous):
    reference={(r['emission_proper_age'],round(r['sky_angle'],12)):r for r in homogeneous}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,7.1))
    colors={2.8:'#176d8a',3.4:'#a56338'}
    for age,color in colors.items():
        rows=[r for r in raw['fine']['rays'] if r['emission_proper_age']==age]
        angle=[math.degrees(r['sky_angle']) for r in rows]
        refs=[reference[(age,round(r['sky_angle'],12))] for r in rows]
        deflection=[math.degrees(r['landing_angle']-h['landing_angle']) for r,h in zip(rows,refs)]
        area=[r['signed_beam_area']/h['signed_beam_area'] for r,h in zip(rows,refs)]
        axes[0,0].plot(angle,deflection,color=color,lw=2,label=f'Emission at screen age {age}')
        axes[0,1].plot(angle,area,color=color,lw=2,label=f'Age {age}')
    axes[0,0].set(title='The formed structure bends off-centre light',xlabel='Inward launch angle (degrees)',ylabel='Landing-angle departure (degrees)')
    axes[0,0].legend(frameon=False,fontsize=8)
    axes[0,1].axhline(0,color='#777777',lw=.7)
    axes[0,1].axhline(1,color='#999999',lw=.7,ls=':')
    axes[0,1].set(title='Beam focusing changes as the core breathes',xlabel='Inward launch angle (degrees)',ylabel='Signed beam area / homogeneous area')
    axes[0,1].legend(frameon=False,fontsize=9,loc='upper left')
    inset=axes[0,1].inset_axes([.57,.12,.39,.29])
    near=[r for r in raw['fine']['rays'] if r['emission_proper_age']==2.8 and r['sky_angle']<=.008]
    inset.plot([math.degrees(r['sky_angle']) for r in near],[r['jacobi_perpendicular'] for r in near],
               'o-',color=colors[2.8],lw=1.2,ms=3)
    inset.axhline(0,color='#555555',lw=.7)
    inset.axvspan(math.degrees(.006),math.degrees(.008),color=colors[2.8],alpha=.08)
    inset.set(title='Early caustic bracket',ylabel=r'$D_\perp$',xlabel='Launch angle (degrees)')
    inset.tick_params(labelsize=7);inset.title.set_fontsize(8)
    inset.xaxis.label.set_size(7);inset.yaxis.label.set_size(8)
    for key,label,color in (('rotation_jacobi_gap','Rotation Jacobi field','#176d8a'),('angular_jacobi_gap','In-plane ray derivative','#a56338')):
        values=[max(abs(r[key]) for r in record['runs'][name]['probes']) for name in ('coarse','medium','fine')]
        axes[1,0].semilogy([1,2,4],values,'o-',color=color,label=label)
    axes[1,0].set(title='Curvature focusing agrees with the ray map',xlabel='Resolution multiplier',ylabel='Maximum principal beam-scale difference',xticks=[1,2,4])
    axes[1,0].legend(frameon=False,fontsize=8)
    for j,(age,color) in enumerate(colors.items()):
        rows=[r for r in record['predictions'] if r['emission_proper_age']==age]
        axes[1,1].bar(np.arange(3)+(j-.5)*.32,[r['relative_beam_flux'] for r in rows],width=.31,color=color,
            label=f'Age {age}',yerr=[r['numerical_comparison_envelope']['relative_beam_flux'] for r in rows],capsize=4)
    axes[1,1].axhline(1,color='#777777',lw=.7)
    axes[1,1].set(title='Received beam flux from the same source law',ylabel='Flux / homogeneous control',
        xticks=[0,1,2],xticklabels=['0.02 rad','0.04 rad','0.08 rad'],xlabel='Declared inward launch angle')
    axes[1,1].legend(frameon=False,fontsize=9)
    for ax in axes.flat:ax.grid(alpha=.15)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-core-optics.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
