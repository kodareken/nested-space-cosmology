"""One formed source, fixed observers and three connected readouts."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(record):
    rows=record['runs'];f=rows['fine'];p=record['conditional_physical_conversion']
    time_factor=p['proper_delay_days']/f['proper_delay_excess'];mas=180/np.pi*3600*1000
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(12,3.7));fan=f['fixed_receiver_fan']
    x=np.array([r['receiver_offset'] for r in fan])/2
    for a,key,factor,title,label in zip(ax,
        ['image_direction_change','proper_delay_excess','log_frequency_change'],
        [mas,time_factor,1e12],
        ['Apparent direction','Proper travel-time excess','Received frequency'],
        ['Image shift (milliarcseconds)','Delay (days)',r'$10^{12}\,\Delta\ln\nu$']):
        y=np.array([r[key] for r in fan])*factor
        a.plot(x,y,'o-',color='#176d8a',label='Evolving source')
        a.errorbar([.04],[f[key]*factor],yerr=[record['numerical_comparison_ranges'][key]*abs(factor)],color='#176d8a',capsize=4)
        a.plot([.04],[f['frozen_comoving_potential'][key]*factor],'x',color='#a56338',ms=8,label='Frozen potential control')
        a.set(title=title,xlabel='Half receiver offset (radians)',ylabel=label)
        a.grid(alpha=.15)
    ax[2].legend(frameon=False,fontsize=8)
    fig.tight_layout();path=Path(__file__).resolve().parents[1]/'docs/figures/nsc-atlas-formed-observer.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)
