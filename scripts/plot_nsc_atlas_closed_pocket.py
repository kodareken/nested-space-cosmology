"""Plot authenticated saved closed-pocket diagnostics and signed ledgers.

No current physical decoder is imported; all five raw archives are verified
before any archive is opened. States themselves are not decoded or evolved.
"""
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/summaries/nsc-atlas-closed-pocket-v1.json'
REFERENCE=ROOT/'results/summaries/nsc-atlas-shape-transition-v1.json'
PRODUCER=ROOT/'scripts/plot_nsc_atlas_closed_pocket.py'
OUTPUT=ROOT/'docs/figures/nsc-atlas-closed-pocket.png'
MANIFEST=OUTPUT.with_suffix('.manifest.json')
LEDGERS=('north_carrier','north_radion','south_carrier','south_radion','north_pressure_work','south_pressure_work')


def identity(path):
    raw=Path(path).read_bytes()
    return dict(path=Path(path).resolve().relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def load_saved():
    summary_binding=identity(SOURCE);reference_binding=identity(REFERENCE)
    record=json.loads(SOURCE.read_text(encoding='utf-8'));reference=json.loads(REFERENCE.read_text(encoding='utf-8'))
    if record.get('schema')!='nsc-atlas-closed-pocket-v1' or record['provenance'].get('source_bound_to_commit') is not True:
        raise ValueError('historically source-bound closed-pocket public record required')
    recovery=record['provenance'].get('recovery',{})
    required=('all_raw_histories_complete_and_admitted','initial_states_recomputed_bitwise',
        'all_public_selected_events_recomputed_exactly','physical_decoder_symbols_match_original_commit')
    if not all(recovery.get(key) is True for key in required) or recovery.get('physics_replayed') is not False:
        raise ValueError('published raw-history recovery checks required; no telemetry flag is inferred')
    prep=record['preparation']
    if prep['carrier_scale']!=.25 or prep['compression']!=.01 or prep['reflection_projection'] is not False:
        raise ValueError('the plotted declared preparation changed')
    if tuple(record['ledger_convention']['order'])!=LEDGERS:raise ValueError('saved physical ledger order changed')
    # Authenticate every archive BEFORE any numerical array is read.
    bindings={};paths={}
    for run in record['runs']:
        path=(ROOT/run['raw_cache']).resolve()
        if not path.is_relative_to(ROOT/'.build') or not path.is_file():raise FileNotFoundError('saved closed-pocket archive required; no evolution started')
        binding=identity(path)
        if binding['sha256']!=run['raw_cache_sha256']:raise ValueError('raw archive SHA256 mismatch: '+run['name'])
        if run['status']!='endpoint_reached':raise ValueError('completed saved run required')
        paths[run['name']]=path;bindings[run['name']]=binding
    if len(bindings)!=5:raise ValueError('all five declared raw runs are required')
    runs={}
    for run in record['runs']:
        with np.load(paths[run['name']],allow_pickle=False) as raw:
            times=raw['times'].copy();ledgers=raw['ledgers'].copy()
            diagnostics=json.loads(str(raw['diagnostics_json']))
        if (times.ndim!=1 or len(times)!=len(diagnostics) or ledgers.shape!=(len(times),6) or
            not np.all(np.isfinite(times)) or not np.all(np.isfinite(ledgers)) or
            not math.isclose(times[0],0.,abs_tol=1e-12) or not math.isclose(times[-1],run['endpoint'],abs_tol=1e-12)):
            raise ValueError('saved time/ledger/diagnostic alignment failed')
        if not all(row['sampled_admitted'] for row in diagnostics):raise ValueError('unadmitted raw diagnostic found')
        if not np.allclose(times,[row['proper_time'] for row in diagnostics],rtol=0,atol=1e-12):
            raise ValueError('saved diagnostics use different proper times')
        index={round(row['proper_time'],9):row for row in diagnostics}
        for selected in run['selected_events']:
            actual=index[round(selected['proper_time'],9)]
            if actual['pole_order']!=['north','center','south'] or not np.array_equal(actual['pole_s'],selected['pole_s']):
                raise ValueError('public selected pole data differ from authenticated diagnostics')
            for field in ('theta_plus','theta_minus'):
                if actual['central_area_maximum'][field]!=selected['central_area_maximum'][field]:
                    raise ValueError('public selected cut differs from the saved diagnostic')
        runs[run['name']]=dict(times=times,ledgers=ledgers,diagnostics=diagnostics,
            selected=run['selected_events'],public=run)
    def stitch(modes):
        past=runs[f'N{modes}_past'];future=runs[f'N{modes}_future']
        if not np.array_equal(past['ledgers'][0],future['ledgers'][0]):raise ValueError('zero-time cumulative ledger mismatch')
        times=np.r_[past['times'][::-1],future['times'][1:]]
        diagnostics=past['diagnostics'][::-1]+future['diagnostics'][1:]
        signed=np.r_[past['ledgers'][::-1,2],future['ledgers'][1:,2]]
        received=signed-past['ledgers'][-1,2]
        expected=record['comparison']['physical_past_to_future'][str(modes)]['caps']['south']['carrier']
        if not math.isclose(received[-1],expected,rel_tol=2e-13,abs_tol=3e-15):raise ValueError('physical past-to-future ledger offset failed')
        return dict(times=times,diagnostics=diagnostics,carrier_received=received,
            offset=float(past['ledgers'][-1,2]))
    joined={modes:stitch(modes) for modes in (64,96)}
    fine=joined[96];times=fine['times'];pole_s=np.asarray([row['pole_s'] for row in fine['diagnostics']],float)
    volume=np.array([row['caps']['south']['volume'] for row in fine['diagnostics']])
    zero=int(np.flatnonzero(times==0.)[0]);s0=pole_s[zero,2]
    volume_ratio=volume/volume[zero];internal_ratio=np.exp(12*(pole_s[:,2]-s0))
    if np.any(volume_ratio<=0) or np.any(internal_ratio<=0):raise ValueError('positive normalized volumes required')
    sectors=reference['stationary_volume_sectors']
    reference_sector=next(row for row in sectors if row['q']==1.5)
    compact=next(row['s'] for row in reference_sector['extrema'] if row['kind']=='minimum')
    barrier=record['flat_horizon_comparison']['barrier_s']
    selected={round(row['proper_time'],9):row for name in ('N96_past','N96_future') for row in runs[name]['selected']}
    events=[selected[t] for t in (-25.,-20.,-10.,0.,10.,20.,25.)]
    if not all(row['central_area_maximum']['found'] and row['central_area_maximum']['spatial_maximum'] for row in events):
        raise ValueError('the seven plotted cuts must be the saved spatial maxima')
    control=runs['N64_future_time_control']
    control_received=control['ledgers'][:,2]-joined[64]['offset']
    return dict(record=record,summary_binding=summary_binding,reference_binding=reference_binding,paths=paths,raw_bindings=bindings,
        runs=runs,joined=joined,times=times,pole_s=pole_s,barrier=barrier,compact=compact,
        events=events,volume_ratio=volume_ratio,internal_ratio=internal_ratio,south_initial_s=float(s0),
        control_times=control['times'],control_received=control_received)


def main():
    producer=identity(PRODUCER);data=load_saved()
    blue='#245b80';orange='#b56326';purple='#71508f';green='#538078';muted='#596572'
    style={'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12,
        'axes.labelsize':10.3,'xtick.labelsize':10,'ytick.labelsize':10,
        'axes.spines.top':False,'axes.spines.right':False}
    with plt.rc_context(style):
        fig,axes=plt.subplots(2,2,figsize=(12.8,8.0))
        fig.subplots_adjust(left=.087,right=.974,bottom=.15,top=.84,wspace=.29,hspace=.56)
        a,b,c,d=axes.ravel();t=data['times'];s=data['pole_s']
        a.plot(t,s[:,0],color=orange,lw=1.7,ls='--',label='North pole')
        a.plot(t,s[:,1],color=blue,lw=2.,label='Center')
        a.plot(t,s[:,2],color=green,lw=1.4,ls=':',label='South pole')
        a.axhline(data['barrier'],color='#7d8893',lw=1.,ls='--',label='Reference barrier (chi=0)')
        a.axhline(data['compact'],color=purple,lw=1.,ls=':',label='Compact minimum')
        a.set(title='A  Saved radion field through the closed history',ylabel='Radion s (dimensionless)')
        a.legend(loc='upper right',frameon=False,fontsize=8.9)
        event_t=np.array([row['proper_time'] for row in data['events']])
        plus=np.array([row['central_area_maximum']['theta_plus'] for row in data['events']])
        minus=np.array([row['central_area_maximum']['theta_minus'] for row in data['events']])
        b.plot(event_t,plus,color=blue,lw=2.,marker='o',ms=4,label=r'$\theta_+$')
        b.plot(event_t,minus,color=orange,lw=1.4,ls='--',marker='s',ms=3,mfc='white',label=r'$\theta_-$ (coincident)')
        b.axhline(0.,color='#7b8590',lw=.8)
        b.set(title='B  Selected central areal maximum: null signs',ylabel=r'$\theta_\pm$ (model inverse-length units)')
        b.legend(loc='upper left',frameon=False,fontsize=10)
        b.text(.97,.07,'Cosmological trapping is not\na localized black-hole classification.',
            transform=b.transAxes,ha='right',fontsize=9.5,color=muted,va='bottom')
        c.plot(data['joined'][96]['times'],data['joined'][96]['carrier_received'],color=blue,lw=2.3,label='N=96')
        c.plot(data['joined'][64]['times'],data['joined'][64]['carrier_received'],color=orange,lw=1.5,ls='--',label='N=64')
        c.plot(data['control_times'],data['control_received'],color='#4a5057',lw=1.2,ls=':',
            label='N=64 future, half-step control')
        c.axhline(0.,color='#7b8590',lw=.8)
        c.set(title='C  South cap: carrier received from t=-25',ylabel='Cumulative carrier energy (model units)')
        c.legend(loc='upper left',frameon=False,fontsize=8.7)
        c.text(.97,.08,f'N96 net: {data["joined"][96]["carrier_received"][-1]:.8f}',
            transform=c.transAxes,ha='right',fontsize=9.5,color=muted)
        d.semilogy(t,data['volume_ratio'],color=blue,lw=2.,label='South-cap Einstein 3-volume')
        d.semilogy(t,data['internal_ratio'],color=orange,lw=1.8,ls='--',label='South-pole internal 12-volume')
        d.set(title='D  Distinct normalized volume measurements',ylabel='Ratio to its t=0 value (log scale)')
        d.legend(loc='upper left',frameon=False,fontsize=8.7)
        for ax in axes.ravel():
            ax.set(xlim=(-25,25),xlabel='Synchronous proper time t (model units)')
            ax.axvline(0.,color='#a3abb3',lw=.7,ls=':',zorder=0)
            ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
        fig.suptitle('A closed continuation with regional carrier transfer',y=.966,fontsize=17,fontweight='bold')
        fig.text(.5,.920,'Conditional preparation: supplied carrier scale 0.25 and compression 0.01; model units.',
            ha='center',fontsize=11,color=muted)
        fig.text(.5,.064,'Fixed coordinate caps are not causal child boundaries. This conditional closed continuation does not establish offspring.',
            ha='center',fontsize=10.5,color='#37414c')
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(OUTPUT,dpi=125,facecolor='white',metadata={'Software':f'Matplotlib {matplotlib.__version__}'})
        plt.close(fig)
    if (identity(SOURCE)!=data['summary_binding'] or identity(REFERENCE)!=data['reference_binding'] or
        identity(PRODUCER)!=producer or any(identity(data['paths'][name])!=binding for name,binding in data['raw_bindings'].items())):
        raise RuntimeError('summary, raw archive or plot producer changed while rendering')
    manifest=dict(schema='nsc-atlas-closed-pocket-figure-v1',producer=producer,
        inputs=dict(scientific_record=data['summary_binding'],compact_reference=data['reference_binding'],raw_archives=data['raw_bindings']),
        historical_source_commit=data['record']['provenance']['producer_commit'],
        original_source_unchanged_telemetry=data['record']['provenance']['source_unchanged'],
        recovery_physics_replayed=data['record']['provenance']['recovery']['physics_replayed'],output=identity(OUTPUT),
        render=dict(size_pixels=[1600,1000],dpi=125,matplotlib_version=matplotlib.__version__),
        plotted=dict(primary_modes=96,history_points=len(t),proper_domain=[-25.,25.],selected_null_events=event_t.tolist(),
            barrier_reference_chi=0.,barrier_s=data['barrier'],compact_minimum_s=data['compact'],
            carrier_ledger_index=2,carrier_ledger_stitch='L_raw2(t)-L_raw2(-25); future time-control usesN64past baseline.',
            physical_south_carrier_net_N64=float(data['joined'][64]['carrier_received'][-1]),
            physical_south_carrier_net_N96=float(data['joined'][96]['carrier_received'][-1]),
            south_initial_s=data['south_initial_s'],
            internal_volume_rule='exp(12*(s_south(t)-s_south(0)))',
            Einstein_cap_volume_rule='V3_south(t)/V3_south(0)'),
        raw_arrays_used=['times','ledgers','diagnostics_json'],states_decoded=False,current_physical_decoder_used=False,
        physical_source_or_PDE_recomputed=False,
        scope='Authenticated original diagnostics and signed ledgers with published recovery checks. One declared closed benchmark, no localized-black-hole/causal-offspring or selected-state claim.',
        reproduction='python scripts/run_python.py scripts/plot_nsc_atlas_closed_pocket.py')
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(figure=manifest['output']['path'],bytes=manifest['output']['bytes'],sha256=manifest['output']['sha256'],
        manifest=str(MANIFEST.relative_to(ROOT)),south_carrier=data['joined'][96]['carrier_received'][-1]),indent=2))


if __name__=='__main__':main()
