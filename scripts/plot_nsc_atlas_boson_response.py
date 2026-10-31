"""Plot the same-parent onset, finite response and unresolved-volume memory."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    r=json.loads((ROOT/'results/summaries/nsc-atlas-boson-response-v1.json').read_text(encoding='utf-8'))
    onset=json.loads((ROOT/'results/summaries/nsc-atlas-boson-onset-v1.json').read_text(encoding='utf-8'))
    shells=json.loads((ROOT/'results/summaries/nsc-atlas-boson-shells-v1.json').read_text(encoding='utf-8'))
    cases=sorted([v['run'] for v in r['cases'] if v['run']['status']=='endpoint_reached'
        and not v['run']['decoupled_control'] and v['run']['max_step']==.02],key=lambda v:v['probe_degree'])
    degrees=[v['probe_degree'] for v in cases]
    colors=['#225b91','#a85d20']
    fig,axes=plt.subplots(2,2,figsize=(11.6,8),layout='constrained')
    a,b,c,d=axes.ravel()
    for source,color in zip((1,2),colors):
        measured=[next(p for p in v['onset_cubic_response'] if p['source_harmonic']==source)
            ['local_log_g_KK_cubic_coefficient'] for v in cases]
        exact=next(p for p in onset['response'] if p['source_harmonic']==source)['local_log_g_KK_proper_cubic_coefficient']
        a.plot(degrees,measured,'o-',color=color,label=f'Source Y{source}')
        a.axhline(exact,color=color,ls='--',alpha=.7)
        ratio=[next(p for p in v['checkpoints'][-1]['response'] if p['source_harmonic']==source)
            ['volume_to_scalar_reference_response_ratio'] for v in cases]
        b.plot(degrees,ratio,'o-',color=color,label=f'Source Y{source}')
    a.set(xlabel='Probe harmonic cutoff',ylabel='Local cubic log-coupling coefficient',
        title='A  Direct field-equation onset (dashed)')
    b.set(xlabel='Probe harmonic cutoff',ylabel='Volume / scalar reference response norm',
        title='B  Endpoint reference response ratio')
    a.legend(fontsize=9);b.legend(fontsize=9)
    memory=r['schur_memory'];matrix=np.array(memory['projected_pchi_from_qchi_memory_over_mu_cubed'])*1e5
    scale=float(np.max(abs(matrix)))
    im=c.imshow(matrix,cmap='coolwarm',norm=TwoSlopeNorm(vmin=-scale,vcenter=0,vmax=scale))
    for i in range(2):
        for j in range(2):c.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',fontsize=12)
    c.set(xticks=[0,1],yticks=[0,1],xticklabels=['Y1','Y2'],yticklabels=['Y1','Y2'],
        xlabel='Source mode',ylabel='Return mode',
        title=f"C  Volume memory at (T,T/2), cutoff {memory['probe_degree']}")
    fig.colorbar(im,ax=c,shrink=.8,label=r'$10^5\mathcal{K}_{pq}/\mu_*^3$ (chart time)')
    spectral=sorted([v for v in shells['cases'] if not v['decoupled_control'] and v['max_step']==.02],key=lambda v:v['probe_degree'])
    for field,color in zip(('S','chi'),colors):
        x=[];y=[]
        for v in spectral:
            point=next(p for p in v['responses'] if p['source_harmonic']==2 and p['field']==field)
            value=point['last_two_degrees_power_fraction']
            if value is not None:x.append(v['probe_degree']);y.append(value)
        d.semilogy(x,y,'o-',color=color,label=field)
    d.set(xlabel='Probe harmonic cutoff',ylabel='Last two degrees: reference-norm power fraction',
        title='D  Source Y2: omitted-mode warning')
    d.legend(fontsize=9)
    for ax in (a,b,d):ax.grid(alpha=.2);ax.set_xticks(degrees)
    fig.suptitle('The same inherited interaction fixes response and memory\n'
        'Actual evolving parent; noise and entropy require the incoming covariance',fontsize=12)
    out=ROOT/'paper/figures/nsc-atlas-boson-response.png'
    fig.savefig(out,dpi=180);plt.close(fig)
    print(out.relative_to(ROOT))


if __name__=='__main__':main()
