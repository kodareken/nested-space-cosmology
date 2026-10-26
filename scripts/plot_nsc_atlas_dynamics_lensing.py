"""The precheck, dynamical calibration and genuinely withheld numerical check."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT=Path(__file__).resolve().parents[1]


def render():
    f=json.loads((ROOT/'results/summaries/nsc-p004-dynamics-lensing-forecast-v1.json').read_text())
    v=json.loads((ROOT/'results/summaries/nsc-p004-dynamics-lensing-validation-v1.json').read_text())
    mas=180/np.pi*3600*1000
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,3,figsize=(12,3.8))
    r=[x['radius'] for x in f['precheck']];a=[x['force_change_per_nominal_orbit'] for x in f['precheck']]
    ax[0].semilogx(r,a,'o-',color='#a56338');ax[0].axhline(.1,color='#777777',ls=':')
    ax[0].set(title='Check stationarity first',xlabel='Comoving radius / L',ylabel='Instantaneous change / nominal orbit')
    probes=f['calibration']['probes'];p=sorted(probes,key=lambda x:x['radius'])
    ax[1].loglog([x['radius'] for x in p],[x['inferred_relative_acceleration']*x['force_conversion_factor'] for x in p],
                '.',color='#176d8a',ms=4)
    ax[1].set(title='Calibrate from free-fall motion',xlabel='Comoving radius / L',ylabel='Reconstructed potential gradient')
    ax[1].set_xticks([.015,.03,.06,.12,.25],labels=['.015','.03','.06','.12','.25'])
    ax[1].xaxis.set_minor_formatter(NullFormatter())
    prediction=f['prediction'];actual=v['runs']['fine']['image_shift']*mas
    ax[2].errorbar([0],[prediction['value_milliarcseconds']],yerr=[prediction['comparison_tolerance_milliarcseconds']],
                   fmt='o',capsize=5,color='#176d8a',label='Frozen prediction')
    ax[2].errorbar([1],[actual],yerr=[v['withheld_numerical_envelope']*mas],
                   fmt='s',capsize=5,color='#a56338',label='Withheld null ray')
    ax[2].set(xticks=[0,1],xticklabels=['Motion prediction','Optical check'],ylabel='Image shift (milliarcseconds)',
              title='Predict before tracing the ray',xlim=(-.5,1.5))
    for a in ax:a.grid(alpha=.15)
    fig.tight_layout();path=ROOT/'docs/figures/nsc-atlas-dynamics-lensing.png'
    fig.savefig(path,dpi=180);plt.close(fig);print(path)


if __name__=='__main__':render()
