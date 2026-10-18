"""Record the frozen wall-radar causal chart, optionally drawing its geometry."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from recursive_horizons.nsc_atlas_dirac_radar import WallRadarChart, reference_radar_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_dirac_radar.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'tests/test_nsc_atlas_dirac_radar.py','scripts/measure_nsc_atlas_dirac_radar.py')


def draw_chart(chart,path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels = np.linspace(0.,chart.duration,129)
    U,V = np.meshgrid(labels,labels,indexing='ij')
    nested = chart.fields(U,V)
    cap = chart.fields(U,V,region='cap')
    x,T = (V-U)/2,(U+V)/2
    fig,axes = plt.subplots(1,2,figsize=(9.2,4.8),layout='constrained')
    panels = ((nested['lapse'],'Wall radar clock','N = d proper time / dT','viridis'),
        (nested['potential']-cap['potential'],'Change in Dirac potential',r'$\Delta W\,\ell$ vs smooth cap','RdBu_r'))
    for ax,(values,title,label,cmap) in zip(axes,panels):
        if cmap == 'RdBu_r':
            extent = max(float(np.max(np.abs(values))),1e-12)
            artist = ax.pcolormesh(x,T,values,shading='gouraud',cmap=cmap,vmin=-extent,vmax=extent,rasterized=True)
        else:
            artist = ax.pcolormesh(x,T,values,shading='gouraud',cmap=cmap,rasterized=True)
        ax.plot([0.,0.],[0.,chart.duration],color='black',lw=1.2)
        ax.plot([0.,chart.duration/2,0.,-chart.duration/2,0.],
                [0.,chart.duration/2,chart.duration,chart.duration/2,0.],color='black',lw=.8)
        ax.set(title=title,xlabel=r'signed wall-radar position $x/\ell$',ylabel=r'wall-radar time $T/\ell$')
        ax.set_aspect('equal')
        ax.text(-chart.duration*.20,chart.duration*.52,'cap',ha='center',fontsize=9)
        ax.text(chart.duration*.20,chart.duration*.52,'exterior',ha='center',fontsize=9)
        fig.colorbar(artist,ax=ax,location='bottom',shrink=.9,pad=.05,label=label)
    event = chart.from_wall_radii(.475,.525)
    u,v,xp,tp = [float(event[key]) for key in ('clock_u','clock_v','radar_position','time')]
    axes[0].plot([0.,xp,0.],[u,tp,v],color='white',lw=1.5,ls='--')
    axes[0].scatter([xp],[tp],color='white',edgecolor='black',s=30,zorder=4)
    axes[0].annotate('radar receiver',(xp,tp),xytext=(9,-18),textcoords='offset points',fontsize=8,
                     bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=2))
    fig.suptitle(r'Quiet nested geometry: $H\ell=1$, $k_{wall}\ell=2$, $M/\ell=0.4$, angular $\kappa=1$',fontsize=11)
    fig.savefig(path,dpi=160)
    plt.close(fig)
    return 129**2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--figure',type=Path)
    args = parser.parse_args()
    for path in (args.output,args.figure):
        if path is not None and (path.exists() or not path.parent.is_dir()):
            raise ValueError('output must have an existing parent and unused path: '+str(path))
    commit = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT) != raw:
            raise ValueError('freeze the radar method before its reference: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    started = time.process_time()
    chart = WallRadarChart()
    report = reference_radar_report(chart)
    figure = None
    if args.figure is not None:
        visual_points = draw_chart(chart,args.figure)
        raw = args.figure.read_bytes()
        figure = dict(path=str(args.figure).replace('\\','/'),sha256=hashlib.sha256(raw).hexdigest(),
                      bytes=len(raw),visual_metric_points=visual_points)
    elapsed = time.process_time()-started
    if elapsed > 30 or report['counts']['wall_ODE_rhs'] > 200000:
        raise RuntimeError('radar geometry reference exceeded its compute budget')
    result = dict(schema='nsc-atlas-dirac-radar-v1',status='finite_causal_radar_chart_checked',
        source_commit=commit,source_sha256=bindings,reference=report,figure=figure,
        runtime_cpu_seconds=elapsed,limits=dict(CPU_seconds=30,wall_ODE_rhs=200000,summary_bytes=10000))
    raw = json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode()) > 10000:
        raise RuntimeError('radar geometry summary exceeds its size budget')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=result['status'],duration=chart.duration,
        lapse_bound=chart.bounds['lapse_lower'],checks=report['checks'],
        radar_probe=report['radar_probe'],CPU_seconds=elapsed,figure=figure)))


if __name__ == '__main__':
    main()
