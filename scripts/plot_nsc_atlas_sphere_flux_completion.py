"""Visualize the public v2 completion record, without raw fields or new solves."""
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/summaries/nsc-atlas-sphere-flux-completion-v2.json'
PRODUCER=ROOT/'scripts/plot_nsc_atlas_sphere_flux_completion.py'
OUTPUT=ROOT/'docs/figures/nsc-atlas-sphere-flux-completion.png'
MANIFEST=OUTPUT.with_suffix('.manifest.json')
PAIRS=tuple(combinations(range(4),2))


def identity(path):
    raw=Path(path).read_bytes()
    return dict(path=Path(path).relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def load_public():
    binding=identity(SOURCE);record=json.loads(SOURCE.read_text(encoding='utf-8'))
    if record.get('schema')!='nsc-atlas-sphere-flux-completion-v1' or record['provenance'].get('source_unchanged') is not True:
        raise ValueError('unchanged published sphere-flux completion record required')
    if 'overall_normalization' not in record['action']:
        raise ValueError('the v2 quantum-normalization disclosure is required')
    examples={e['name']:e for e in record['inheritance_embedding_comparison']['examples']}
    diagonal=examples['diagonal_six_cycle'];single=examples['single_pair_cycle']
    for example,period in ((diagonal,Fraction(1,96)),(single,Fraction(1,16))):
        if (Fraction(example['period_unit_squared'])!=period or
            Fraction(example['reduced_flux_energy'])!=Fraction(25,32) or
            Fraction(example['inherited_step_norm_squared'])!=Fraction(1,16)):
            raise ValueError('unexpected conditional flux-unit/energy comparison')
        if len(example['integer_periods'])!=6 or len(example['initial_shape_accelerations'])!=4:
            raise ValueError('six pair periods and four initial shape forces required')
    if diagonal['integer_periods']!=[5]*6 or single['integer_periods']!=[5,0,0,0,0,0]:
        raise ValueError('published conditional period vectors changed')
    slices=record['slices']
    if len(slices)!=29:raise ValueError('the29 recorded frozen slices are required')
    times=np.array([r['proper_time'] for r in slices])
    local=np.array([r['sampled_mass_range'][0] for r in slices])
    retained=np.array([r['first_four_eigenvalues'][0] for r in slices])
    if np.any(np.diff(times)<=0) or not np.all(np.isfinite([times,local,retained])):
        raise ValueError('finite ordered frozen-slice diagnostics required')
    accelerations=np.array([diagonal['initial_shape_accelerations'],single['initial_shape_accelerations']],float)
    if not np.all(np.isfinite(accelerations)):raise ValueError('finite initial shape forcing required')
    return record,binding,diagonal,single,times,local,retained,accelerations


def draw_embedding(ax,example,title,color):
    ax.set_axis_off();ax.set_title(title,loc='left',pad=12)
    graph=ax.inset_axes([.005,.005,.50,.94]);graph.set_aspect('equal');graph.set_axis_off()
    graph.set_xlim(-1.16,1.16);graph.set_ylim(-1.12,1.12)
    positions=np.array([[-.72,.72],[.72,.72],[.72,-.72],[-.72,-.72]])
    label_positions=[(0.,.80),(-.27,.27),(-.82,0.),(.82,0.),(.27,.27),(0.,-.80)]
    for (i,j),number,position in zip(PAIRS,example['integer_periods'],label_positions):
        if number==0:continue
        graph.plot(*positions[[i,j]].T,color=color,lw=2.,zorder=1)
        graph.text(*position,str(number),ha='center',va='center',color=color,fontsize=11,
            bbox=dict(boxstyle='round,pad=.10',fc='white',ec='none'),zorder=3)
    for i,(x,y) in enumerate(positions):
        graph.add_patch(Circle((x,y),.185,fc='#f1f5f8',ec=color,lw=1.6,zorder=4))
        graph.text(x,y,rf'$S^3_{{{i+1}}}$',ha='center',va='center',fontsize=11,zorder=5)
    vector=', '.join(map(str,example['integer_periods']))
    ax.text(.56,.81,f'N = ({vector})',transform=ax.transAxes,fontsize=10.5,va='top')
    ax.text(.56,.61,rf'Period unit$^2$: ${example["period_unit_squared"]}$',
        transform=ax.transAxes,fontsize=11)
    ax.text(.56,.43,r'$T=25/32$',transform=ax.transAxes,fontsize=11.5)
    ax.text(.56,.25,r'Effective unit norm: $1/4$',transform=ax.transAxes,fontsize=11)
    ax.text(.56,.06,'Edge labels: integer periods' if all(example['integer_periods']) else 'Unmarked pair cycles: zero flux',
        transform=ax.transAxes,color='#5d6772',fontsize=9)


def main():
    producer=identity(PRODUCER)
    record,binding,diagonal,single,times,local,retained,accelerations=load_public()
    blue='#245b80';orange='#b56326';gray='#8d99a5';muted='#586471'
    style={'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,
        'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,
        'axes.spines.top':False,'axes.spines.right':False}
    with plt.rc_context(style):
        fig,axes=plt.subplots(2,2,figsize=(12.8,8.0))
        fig.subplots_adjust(left=.083,right=.975,bottom=.165,top=.835,wspace=.29,hspace=.55)
        a,b,c,d=axes.ravel()
        draw_embedding(a,diagonal,'A  Conditional diagonal six-cycle embedding',blue)
        draw_embedding(b,single,'B  Conditional single-pair embedding',orange)
        c.plot(times,retained,color=blue,lw=2.,marker='o',ms=2.5,label=r'Lowest retained radial $\lambda_0$')
        c.plot(times,local,color=orange,lw=1.9,ls='--',label=r'Sampled local minimum $m^2(s)$')
        c.axhline(0.,color='#525d68',lw=.9)
        c.set(title='C  Diagonal history: local mass and retained mode',xlim=(0,2.8),
            xlabel='Synchronous proper time t (model units)',
            ylabel=r'$m^2$, $\lambda_0$ (model inverse-length$^2$)')
        c.legend(loc='upper right',frameon=False,fontsize=9.5)
        c.text(.47,.55,r'Local $m^2<0$ does not establish'+'\na global growing mode.',
            transform=c.transAxes,fontsize=9.5,color=muted,va='bottom')
        c.grid(axis='y',alpha=.18);c.set_axisbelow(True)
        x=np.arange(1,5);width=.34
        d.bar(x-width/2,accelerations[0],width=width,color=gray,label='Diagonal embedding',zorder=3)
        d.scatter(x-width/2,accelerations[0],s=20,color=gray,zorder=4)
        d.bar(x+width/2,accelerations[1],width=width,color=orange,label='Single-pair embedding',zorder=3)
        d.axhline(0.,color='#525d68',lw=.9)
        for coord,value in zip(x+width/2,accelerations[1]):
            d.text(coord,value+(.10 if value>=0 else -.10),f'{value:+.3f}',ha='center',
                va='bottom' if value>=0 else 'top',fontsize=10)
        d.set(title='D  INITIAL relative-radius forcing only',xticks=x,
            xticklabels=[rf'$S^3_{{{i}}}$' for i in x],xlabel='Internal sphere factor',
            ylabel=r'Initial $\ddot{\beta}_i$ (model time$^{-2}$)',ylim=(-2.5,3.0),yticks=[-2,-1,0,1,2])
        d.legend(loc='upper center',ncol=2,frameon=False,fontsize=9)
        d.grid(axis='y',alpha=.18);d.set_axisbelow(True)
        fig.suptitle('One scalar flux term, distinct internal embeddings',y=.966,fontsize=17,fontweight='bold')
        fig.text(.5,.920,'Conditional charge normalizations differ: these are not two states of one fixed microscopic lattice.',
            ha='center',fontsize=11,color=muted)
        fig.text(.5,.079,'Topology and wrapped cycles are supplied inputs; frozen diagnostics do not select an SI scale or demonstrate a child.',
            ha='center',fontsize=10.5,color='#37414c')
        fig.text(.5,.047,'Only the diagonal embedding supplies the plotted history. Single-pair bars are initial forcing, not an evolved history.',
            ha='center',fontsize=10.5,color='#37414c')
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(OUTPUT,dpi=125,facecolor='white',metadata={'Software':f'Matplotlib {matplotlib.__version__}'})
        plt.close(fig)
    if identity(SOURCE)!=binding or identity(PRODUCER)!=producer:raise RuntimeError('public input or plot producer changed during rendering')
    manifest=dict(schema='nsc-atlas-sphere-flux-completion-figure-v1',producer=producer,
        input=binding,historical_source_commit=record['provenance']['producer_commit'],
        output=identity(OUTPUT),render=dict(size_pixels=[1600,1000],dpi=125,matplotlib_version=matplotlib.__version__),
        plotted=dict(pair_order_1_based=[[i+1,j+1] for i,j in PAIRS],
            period_unit_squared=[diagonal['period_unit_squared'],single['period_unit_squared']],
            period_vectors=[diagonal['integer_periods'],single['integer_periods']],
            common_T='25/32',common_effective_unit_norm='1/4',initial_relative_radius_accelerations=accelerations.tolist(),
            frozen_slices=len(times),proper_domain=[float(times[0]),float(times[-1])],
            minimum_local_mass_squared=float(np.min(local)),minimum_retained_eigenvalue=float(np.min(retained))),
        scope='Two conditional different microscopic charge lattices with equal reduced scalar data. Only diagonal-background frozen slices are plotted; initial single-pair forcing is not its subsequent history. No full-continuum stability, SI selection or child claim.',
        data_reading_boundary='Public v2 summary and its recorded historical source identity; no raw cache or current physical decoder read.',
        raw_fields_read=False,new_source_calculation=False,PDE_evolution_performed=False,
        reproduction='python scripts/run_python.py scripts/plot_nsc_atlas_sphere_flux_completion.py')
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(figure=manifest['output']['path'],sha256=manifest['output']['sha256'],
        bytes=manifest['output']['bytes'],manifest=str(MANIFEST.relative_to(ROOT)),
        source_commit=manifest['historical_source_commit']),indent=2))


if __name__=='__main__':main()
