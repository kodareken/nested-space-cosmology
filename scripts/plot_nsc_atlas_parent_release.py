"""Plot the source-bound passive release result; no dynamics or fitted curve."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    record=json.loads((ROOT/'results/summaries/nsc-atlas-parent-release-v2.json').read_bytes())
    row=next(r for r in record['runs'] if r['degree']==24)
    path=(ROOT/row['raw_path']).resolve()
    if not path.is_relative_to(ROOT/'.build') or hashlib.sha256(path.read_bytes()).hexdigest()!=row['raw_sha256']:
        raise ValueError('Exact retained passive raw series required')
    raw=json.loads(path.read_bytes()); z=raw['response']['series']; p=raw['point_series']; r=raw['region_series']
    t=np.asarray(z['tau']); s=np.asarray(z['s'])
    labels=z['force_labels']; work=np.asarray(z['cumulative_work'])
    mechanical=np.asarray(z['mechanical_radion_energy'])
    top=s[-1]+z['bare_top_distance'][-1]
    with plt.rc_context({'font.size':10,'axes.spines.top':False,'axes.spines.right':False}):
        fig,axes=plt.subplots(1,3,figsize=(12.8,4.0),layout='constrained')
        a,b,c=axes
        a.plot(t,s,color='#006b75',lw=2,label='Measured internal log size s')
        a.axhline(top,color='#a54b20',ls='--',label='Bare potential top')
        a.set(title='A   Internal size and bare barrier',ylabel='Internal log size')
        a.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,loc='center left',bbox_to_anchor=(.02,.76),fontsize=8)
        b.plot(t,work[:,labels.index('carrier_force')],color='#006b75',lw=2,label='Wave work')
        b.plot(t,work[:,labels.index('spatial_laplacian')],color='#a54b20',lw=2,label='Spatial-gradient work')
        b.plot(t,mechanical-mechanical[0],color='#3f4c55',lw=2,label='Mechanical energy gain')
        b.axhline(z['bare_mechanical_margin'][0],color='#88949b',ls='--',label='Initial barrier gap')
        b.set(title='B   Measured mechanical work',ylabel='Local mechanical energy / work')
        b.legend(frameon=True,facecolor='white',edgecolor='white',framealpha=.95,loc='upper left',fontsize=8)
        vol=np.array([x['volume'] for x in r]); energy=np.array([x['energy'] for x in r])
        c.plot(t,np.exp(12*(s-s[0])),color='#006b75',lw=2,label='Source-centre internal volume')
        c.plot(t,vol/vol[0],color='#496d9c',lw=2,label='Region Einstein-frame volume')
        c.plot(t,(energy/vol)/(energy[0]/vol[0]),color='#a54b20',lw=2,label='Region Einstein-frame density')
        c.axhline(1,color='#88949b',lw=.8)
        c.set(title='C   Internal volume and regional density',ylabel='Ratio to initial value')
        c.legend(frameon=False,loc='upper left',fontsize=8)
        for ax in axes:
            ax.set_xlabel('Source-anchor proper duration')
            ax.grid(alpha=.18)
            ax.set_xlim(0,.9)
        fig.suptitle('Unchanged parent: measured geometric driving, no barrier exit through 0.9',fontsize=13)
        output=ROOT/'docs/figures/nsc-atlas-parent-release.png'
        output.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(output,dpi=180)
        plt.close(fig)
        print(str(output))


if __name__=='__main__':main()
