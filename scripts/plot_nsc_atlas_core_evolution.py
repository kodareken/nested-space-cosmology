"""The inherited diffuse cloud, its gravity and the received clock."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def render(record,plots):
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(11.8,7.1))
    h=plots['diffuse_fine'];ctrl=plots['fixed_geometry_control'];r=record['runs']
    ax[0,0].plot(h[:,0],h[:,-1],color='#166b85',lw=2,label='Coupled scalar and gravity')
    ax[0,0].axhline(r['diffuse_fine']['initial_half_mass_radius'],color='#999999',ls='--',label='Initial half-mass radius')
    ax[0,0].set(title='The diffuse cloud contracts and breathes',xlabel='Parent time (atlas units)',ylabel='Half excess-mass radius inside R = .12')
    ax[0,0].legend(frameon=False,fontsize=9)
    ax[0,1].semilogy(h[:,0],h[:,11],color='#166b85',lw=2,label='Self-consistent gravity')
    ax[0,1].semilogy(ctrl[:,0],ctrl[:,11],color='#b06b37',lw=1.7,label='Same field on homogeneous geometry')
    ax[0,1].set(title='Gravity drives the concentration',xlabel='Parent time (atlas units)',ylabel='Central scalar excitation density')
    ax[0,1].legend(frameon=False,fontsize=8)
    for family,color in (('core','#9c4360'),('diffuse','#166b85')):
        err=[max(r[family+'_'+level]['maximum_mass_flux_error']) for level in ('coarse','medium','fine')]
        ax[1,0].semilogy([1,2,4],err,'o-',color=color,label=family.capitalize())
    ax[1,0].set(title='Independent mass-flux discrepancy shrinks',xlabel='Resolution multiplier',ylabel='Maximum integrated geometric-mass gap',xticks=[1,2,4])
    ax[1,0].legend(frameon=False,fontsize=9)
    signals=record['observer_predictions'];values=[1000*signals[k]['frequency_fractional_departure'] for k in ('core','diffuse')]
    errors=[1000*signals[k]['fractional_departure_numerical_error'] for k in ('core','diffuse')]
    ax[1,1].bar([0,1],values,color=['#9c4360','#166b85'],width=.55,yerr=errors,capsize=6)
    ax[1,1].axhline(0,color='#999999',lw=.7)
    ax[1,1].set(title='A clock signal from the evolved geometry',ylabel='Frequency departure from control (per mille)',xticks=[0,1],
        xticklabels=['Core control\nemitted at central age 0','Diffuse cloud\nemitted at central age 4'])
    for a in ax.flat:a.grid(alpha=.15)
    fig.tight_layout();out=ROOT/'docs/figures/nsc-atlas-core-evolution.png'
    fig.savefig(out,dpi=170);plt.close(fig);print(out)
