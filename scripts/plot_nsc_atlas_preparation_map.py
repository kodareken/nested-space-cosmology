"""Render the physical dependency map and retained state-control evidence."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
import numpy as np
from recursive_horizons.nsc_atlas_volume_preparation import PreparedBackground

ROOT=Path(__file__).resolve().parents[1]

def main():
    r=json.loads((ROOT/'results/summaries/nsc-atlas-exit-preparations-v1.json').read_text(encoding='utf-8'))
    colors=['#176d8a','#b65e22','#765198']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(9.4,7.2))
    fig.suptitle('One preparation must determine curvature and the volume history',fontsize=13,y=.975)
    ax=fig.add_axes([.035,.555,.93,.36]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
    boxes=[(.16,.75,'Candidate action\nq, e, tension, V(s)\nSpecified','#dceaf4'),
           (.50,.75,'Exit to reference state\ns, v, fields at E = 3/8\nPreparation missing','#ffebcf'),
           (.84,.75,'Wall + bulk data\nR, M, normal gradients\nC/E0 still conditional','#ffebcf'),
           (.50,.23,'Evolution + light cone\nw_s(z), H(z), distances\nComputed conditionally','#d6ece5'),
           (.16,.23,'SN + DESI\nOne common ruler\nTested','#dceaf4')]
    for x,y,label,color in boxes:
        ax.add_patch(FancyBboxPatch((x-.137,y-.17),.274,.34,boxstyle='round,pad=.012',
                                  facecolor=color,edgecolor='#526778',lw=.8))
        ax.text(x,y,label,ha='center',va='center',fontsize=9.6,linespacing=1.5)
    for p,q in [((.31,.75),(.35,.75)),((.65,.75),(.69,.75)),((.84,.56),(.65,.27)),((.35,.23),(.31,.23))]:
        ax.add_patch(FancyArrowPatch(p,q,arrowstyle='-|>',mutation_scale=13,color='#526778',lw=1.2))
    ax.text(.82,.19,'E = 3/8 labels a\nreference event, not\na derived birth state.',ha='center',va='center',fontsize=9.5,color='#87511d')
    left=fig.add_axes([.09,.16,.40,.29]);right=fig.add_axes([.59,.16,.36,.29])
    labels=['Matched','HM + matter','Minimum + matter'];preps=['matched','HM_location','minimum']
    bas=r['matched_flat_baseline'];rows=[bas,r['models']['HM_location_flat'],r['models']['minimum_flat']]
    z=np.linspace(0,2,241);deltas=[];omegas=[]
    lcdm=r['flat_LCDM_control']['joint_best']['chi_squared_total']
    for prep,label,color,row in zip(preps,labels,colors,rows):
        bg=PreparedBackground(2.33,preparation=prep)
        No=row['joint_best']['theta'];w=bg.eos(No-np.log1p(z))['w_s']
        left.plot(z,w,color=color,lw=1.8,label=label)
        deltas.append(row['joint_best']['chi_squared_total']-lcdm)
        omegas.append(row['joint_best']['Omega_volume'])
    left.set(title='Pressure history of the volume component',xlabel='Redshift z',ylabel=r'$w_s=p_s/\rho_s$',ylim=(-1.12,1.12))
    left.axhline(-1,color='#777777',lw=.7,ls='--');left.grid(alpha=.15)
    handles,names=left.get_legend_handles_labels()
    fig.legend(handles,names,fontsize=9,frameon=False,ncol=3,loc='center',bbox_to_anchor=(.5,.505))
    bars=right.bar(range(3),deltas,color=colors,width=.58)
    for bar,value in zip(bars,deltas):right.text(bar.get_x()+bar.get_width()/2,value+3,f'{value:.2f}',ha='center',fontsize=9)
    right.set(title='Full fit with flatness imposed',ylabel=r'$\Delta\chi^2$ versus flat $\Lambda$CDM',ylim=(0,max(deltas)*1.18),
              xticks=range(3),xticklabels=[f'{label}\n'+rf'$\Omega_{{s,0}}={omega:.3f}$' for label,omega in zip(labels,omegas)])
    right.tick_params(axis='x',labelsize=8);right.grid(axis='y',alpha=.15);right.set_axisbelow(True)
    fig.text(.5,.055,'C = 0 controls; fixed E_i = 3/8; no initial displacement fitted.',ha='center',fontsize=9)
    fig.text(.5,.028,'A physical exit and wall-compatible curvature remain to be derived.',ha='center',fontsize=9,color='#76511e')
    path=ROOT/'docs/figures/nsc-atlas-preparation-map.png';fig.savefig(path,dpi=190);plt.close(fig);print(path)

if __name__=='__main__':main()
