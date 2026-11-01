"""Plot the public directional connection summaries; no raw fields or evolution.

The synchronous minimum and continued constant-X observer have separate axes.
Pressure work is not reconstructed by subtracting the plotted energy curves.
"""
import ast
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PRODUCER=ROOT/'scripts/plot_nsc_atlas_directional_connection.py'
OUTPUT=ROOT/'docs/figures/nsc-atlas-directional-connection.png'
MANIFEST=OUTPUT.with_suffix('.manifest.json')
FILES={key:ROOT/'results/summaries'/name for key,name in (
    ('initial','nsc-atlas-valve-directional-v1.json'),
    ('synchronous_initial','nsc-atlas-valve-directional-evolution-v1.json'),
    ('synchronous_continued','nsc-atlas-valve-directional-evolution-v2.json'),
    ('conformal','nsc-atlas-valve-directional-conformal-v1.json'))}
MODEL=ROOT/'src/recursive_horizons/nsc_atlas_parent_scale_selector.py'


def identity(path):
    raw=Path(path).read_bytes()
    return dict(path=Path(path).relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def declared_dimension():
    # Read only the declared default; do not instantiate/evolve a model or
    # rebind the summaries to a current physical decoder.
    tree=ast.parse(MODEL.read_text(encoding='utf-8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='VolumeParent')
    node=next(n for n in cls.body if isinstance(n,ast.AnnAssign) and isinstance(n.target,ast.Name) and n.target.id=='d')
    d=ast.literal_eval(node.value)
    if d!=12:raise ValueError('the published volume records declare the d=12 candidate')
    return d


def load_public():
    bindings={key:identity(path) for key,path in FILES.items()}
    records={key:json.loads(path.read_text(encoding='utf-8')) for key,path in FILES.items()}
    initial=records['initial'];provenance=initial['provenance']
    if (initial['schema']!='nsc-atlas-valve-directional-initial-v1' or
        provenance.get('source_unchanged') is not True or provenance['source_sha256_before']!=provenance['source_sha256_after']):
        raise ValueError('unchanged published directional preparation required')
    if initial['preparation']['S']!=0.:raise ValueError('initial log-volume point requires the declared S=0 state')
    rows_by_time={};coarse=None;fine=None
    for key in ('synchronous_initial','synchronous_continued'):
        record=records[key]
        if record['schema']!='nsc-atlas-valve-directional-evolution-v1':raise ValueError('unexpected synchronous summary schema')
        for run in record['runs']:
            if run['authentication'].get('historical_Git_source_verified') is not True:
                raise ValueError('published run must record its historical source verification')
        order=record['comparison']['run_order']
        if coarse is None:coarse,fine=min(order),max(order)
        if sorted(order)!=[coarse,fine]:raise ValueError('synchronous matched grids changed across the handoff')
        table=record['comparison']['matched_times'];columns=table['columns']
        for values in table['rows']:
            row=dict(zip(columns,values));time=round(row['proper_time'],10)
            previous=rows_by_time.get(time)
            if previous is not None:
                for quantity in ('minimum_radius','fixed_cap_energy','net_material_received'):
                    name=f'{quantity}_N{fine}'
                    if not math.isclose(previous[name],row[name],rel_tol=2e-12,abs_tol=2e-10):
                        raise ValueError('published synchronous handoff is discontinuous: '+name)
            rows_by_time[time]=row
    rows=[rows_by_time[t] for t in sorted(rows_by_time)]
    time=np.array([r['proper_time'] for r in rows])
    radius=np.array([r[f'minimum_radius_N{fine}'] for r in rows])
    energy=np.array([r[f'fixed_cap_energy_N{fine}'] for r in rows])
    material=np.array([r[f'net_material_received_N{fine}'] for r in rows])
    baseline=float(initial['initial_cap_ledger']['plus_cap_proper_matter_energy'])
    if not math.isclose(energy[0],baseline,rel_tol=2e-12,abs_tol=2e-10):
        raise ValueError('published fixed-cap baseline differs from initial preparation')
    if not math.isclose(time[0],0.,abs_tol=1e-12) or not math.isclose(time[-1],2.8,abs_tol=1e-12):
        raise ValueError('expected matched admitted synchronous domain0..2.8')
    conformal=records['conformal']
    if conformal['schema']!='nsc-atlas-directional-conformal-comparison-v1':raise ValueError('unexpected conformal summary schema')
    aligned=conformal['aligned_center_proper_time'];times=np.array(aligned['times'],float)
    values=np.array(aligned['fine'],float);index=aligned['columns'].index('radius')
    if len(times)!=9 or values.shape[0]!=9 or not math.isclose(times[0],time[-1],abs_tol=1e-12):
        raise ValueError('nine aligned observer samples from the same handoff required')
    observer_radius=values[:,index]
    if not math.isclose(observer_radius[0],radius[-1],rel_tol=2e-12,abs_tol=2e-10):
        raise ValueError('conformal Cauchy radius does not match the synchronous handoff')
    continued=max(records['synchronous_continued']['runs'],key=lambda r:r['modes'])
    modulus=[continued['initial']['volume_modulus'],continued['endpoint']['volume_modulus']]
    d=declared_dimension();sstar=[]
    for m in modulus:
        sstar.extend([s-v/d for s,v in zip(m['sampled_s_range'],m['log_internal_volume_ratio_range'])])
    if max(sstar)-min(sstar)>3e-13:raise ValueError('published volume normalization is inconsistent')
    barrier_s=modulus[0]['bare_potential_barrier_s']
    if not math.isclose(barrier_s,modulus[1]['bare_potential_barrier_s'],abs_tol=2e-14):
        raise ValueError('published bare barrier changed between epochs')
    barrier=d*(barrier_s-sstar[0])
    log_volume=np.array([0.,modulus[0]['log_internal_volume_ratio_range'][1],modulus[1]['log_internal_volume_ratio_range'][1]])
    arrays=[time,radius,energy,material,times,observer_radius,log_volume]
    if any(not np.all(np.isfinite(a)) for a in arrays):raise ValueError('finite published plot values required')
    return dict(records=records,bindings=bindings,time=time,radius=radius,
        delta_energy=energy-baseline,material=material,baseline=baseline,
        observer_elapsed=times-times[0],observer_radius=observer_radius,
        log_volume=log_volume,barrier=barrier,barrier_s=barrier_s,sstar=sstar[0],d=d,coarse=coarse,fine=fine)


def main():
    producer=identity(PRODUCER);model=identity(MODEL);data=load_public()
    blue='#245b80';orange='#b56326';purple='#71508f';muted='#5d6772';red='#9c4345'
    style={'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12.5,
        'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,
        'axes.spines.top':False,'axes.spines.right':False,'axes.titlelocation':'left'}
    with plt.rc_context(style):
        fig,axes=plt.subplots(2,2,figsize=(12.8,8.0))
        fig.subplots_adjust(left=.082,right=.974,bottom=.15,top=.84,wspace=.28,hspace=.60)
        a,b,c,d=axes.ravel()
        a.plot(data['time'],data['delta_energy'],color=blue,lw=2.1,marker='o',ms=2.6,label='Fixed-cap energy change')
        a.plot(data['time'],data['material'],color=orange,lw=2.1,label='Cumulative material received')
        a.set(title='A  Fixed cap: material supply and total energy',
            xlabel='Synchronous proper time t (model units)',ylabel='Energy (model units)',xlim=(0,2.8))
        a.legend(loc='upper left',frameon=False,fontsize=9.5)
        a.text(.97,.07,'Energy change includes pressure work',transform=a.transAxes,
            ha='right',va='bottom',fontsize=9,color=muted)
        b.plot(data['time'],data['radius'],color=blue,lw=2.1,marker='o',ms=3)
        b.set(title='B  Migrating spatial minimum',xlabel='Synchronous proper time t (model units)',
            ylabel='Minimum radius R (model length units)',xlim=(0,2.8))
        b.text(.04,.10,f'Matched N={data["coarse"]} / N={data["fine"]}; plotted N={data["fine"]}',
            transform=b.transAxes,color=muted,fontsize=9.5)
        for ax in (a,b):ax.axvline(2.,color='#a7adb4',lw=.9,ls=':',zorder=0)
        c.plot(data['observer_elapsed'],data['observer_radius'],color=purple,lw=2.1,marker='o',ms=4)
        c.set(title='C  Continued constant-X observer',
            xlabel=r'New observer elapsed proper time $\tau-2.8$ (model units)',
            ylabel='Observer radius R (model length units)',xlim=(0,float(data['observer_elapsed'][-1])))
        c.text(.04,.10,'Same Cauchy handoff; this curve does not track\nthe migrating minimum.',
            transform=c.transAxes,color=muted,fontsize=9.5,va='bottom')
        epochs=np.array([0.,2.,2.8]);heights=data['log_volume']
        d.bar(epochs,heights,width=.34,color=['#8f9aa6','#8a91ab',purple],zorder=3)
        d.plot([0.],[0.],'o',color='#616c77',ms=4,zorder=4)
        d.axhline(data['barrier'],color=red,lw=1.4,ls='--',label=f'Bare-potential barrier ({data["barrier"]:.3f})')
        for x,value in zip(epochs,heights):d.text(x,value+.09,f'{value:.3f}',ha='center',fontsize=10.5)
        d.set(title='D  Auxiliary internal-volume response',xlabel='Synchronous proper time t (model units)',
            ylabel=r'Sampled $\max_x\,\ln(V_{\rm int}/V_{{\rm int},*})$',
            xticks=epochs,xticklabels=['0','2','2.8'],xlim=(-.4,3.2),ylim=(0,float(max(heights))+.45))
        d.legend(loc='upper left',frameon=False,fontsize=9.5)
        d.text(.05,.64,'d = 12',transform=d.transAxes,color=muted,fontsize=10)
        for ax in axes.ravel():
            ax.grid(axis='y',alpha=.19,lw=.7);ax.set_axisbelow(True)
            ax.tick_params(direction='out')
        fig.suptitle('Directional energy transfer and geometric continuation',y=.965,fontsize=17,fontweight='bold')
        fig.text(.5,.918,'Synchronous minimum tracking and conformal observer continuation use different clocks and paths.',
            ha='center',fontsize=11,color=muted)
        fig.text(.5,.065,'Finite records of one declared preparation; no formed child or full higher-dimensional solution is established.',
            ha='center',fontsize=10.5,color='#37414c')
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(OUTPUT,dpi=125,facecolor='white',metadata={'Software':f'Matplotlib {matplotlib.__version__}'})
        plt.close(fig)
    if identity(PRODUCER)!=producer or identity(MODEL)!=model or any(identity(FILES[k])!=v for k,v in data['bindings'].items()):
        raise RuntimeError('figure input or producer changed while rendering')
    manifest=dict(schema='nsc-atlas-directional-connection-figure-v1',producer=producer,inputs=data['bindings'],
        dimension_reference=dict(**model,symbol='VolumeParent.d',value=data['d'],method='Read declared default by AST; no model evolution.'),
        output=identity(OUTPUT),render=dict(size_pixels=[1600,1000],dpi=125,matplotlib_version=matplotlib.__version__),
        plotted=dict(synchronous_points=len(data['time']),synchronous_domain=[float(data['time'][0]),float(data['time'][-1])],
            fixed_cap_energy_baseline=data['baseline'],final_cap_energy_change=float(data['delta_energy'][-1]),
            final_cumulative_material_received=float(data['material'][-1]),
            observer_samples=len(data['observer_elapsed']),observer_elapsed_domain=[0.,float(data['observer_elapsed'][-1])],
            observer_radius_end=float(data['observer_radius'][-1]),
            log_volume_epochs=[0.,2.,2.8],sampled_max_log_volume_ratios=data['log_volume'].tolist(),
            bare_barrier_log_ratio=data['barrier'],bare_barrier_s=data['barrier_s'],record_inferred_sstar=data['sstar']),
        scope='Published summary values and their recorded historical provenance; no current decoder rebinding, raw field loading, PDE, independent work curve, formed-child or complete uplift claim.',
        pressure_work_curve_inferred=False,raw_fields_read=False,physical_sources_or_PDE_recomputed=False,
        reproduction='python scripts/run_python.py scripts/plot_nsc_atlas_directional_connection.py')
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(figure=str(OUTPUT.relative_to(ROOT)),manifest=str(MANIFEST.relative_to(ROOT)),
        sha256=manifest['output']['sha256'],bytes=manifest['output']['bytes'],plotted=manifest['plotted']),indent=2))


if __name__=='__main__':main()
