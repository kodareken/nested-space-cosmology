#!/usr/bin/env python3
"""Plot authenticated saved initial Dirac operator coefficients; never solve.

The squared matrix magnitudes and normalized column split are reference-
amplitude diagnostics. They are not probabilities, occupation or energy.
"""
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np


ROOT=Path(__file__).resolve().parents[1]
RECORD=ROOT/'results/summaries/nsc-atlas-valve-dirac-source-v1.json'
OUTPUT=ROOT/'paper/figures/nsc-atlas-dirac-source.png'
MANIFEST=OUTPUT.with_suffix('.manifest.json')
PRODUCER='scripts/plot_nsc_atlas_valve_dirac_source.py'
REGENERATE='python scripts/run_python.py scripts/check_nsc_atlas_valve_dirac_source.py'


def identity(path):
    raw=Path(path).read_bytes()
    return dict(path=Path(path).resolve().relative_to(ROOT).as_posix(),
                sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def load_saved():
    if not RECORD.is_file():raise FileNotFoundError(f'Missing source record. Regenerate with: {REGENERATE}')
    record_identity=identity(RECORD)
    record=json.loads(RECORD.read_text(encoding='utf-8'))
    if record.get('schema')!='nsc-atlas-valve-dirac-source-v1':raise ValueError('unexpected source-operator schema')
    provenance=record['provenance']
    if provenance.get('source_unchanged') is not True or provenance['source_sha256_before']!=provenance['source_sha256_after']:
        raise ValueError('source record does not bind an unchanged scientific producer')
    levels=[item for item in record['levels'] if item['k']==0]
    if not levels:raise ValueError('saved record has no k=0 source matrix')
    finest=max(levels,key=lambda item:item['degree'])
    binding=finest['raw_matrices']
    cache=(ROOT/str(binding['path']).replace('\\','/')).resolve()
    if not cache.is_relative_to(ROOT):raise ValueError('matrix cache escapes the repository')
    if not cache.is_file():
        raise FileNotFoundError(f'Missing authenticated matrix cache {cache}. Regenerate with: {REGENERATE}. No solver was started.')
    cache_identity=identity(cache)
    if cache_identity['sha256']!=binding['sha256']:raise ValueError('saved matrix cache SHA256 mismatch')
    # Authenticate BEFORE reading/plotting any matrix values. The archive's
    # source identities are compared to its saved record, not silently rebound
    # to whatever scientific source happens to be checked out while plotting.
    with np.load(cache,allow_pickle=False) as data:
        metadata=json.loads(str(data['metadata']))
        pair=np.array(data['positive_negative_pair_hessian'],copy=True)
        eigenvalues=np.array(data['spatial_eigenvalues'],dtype=float,copy=True)
        index=int(data['ground_index'])
    if (metadata['source_sha256']!=provenance['source_sha256_before'] or
            metadata['input_record_sha256']!=provenance['input_record_sha256'] or
            metadata['defined_metric_sha256']!=provenance['defined_metric_cache']['sha256']):
        raise ValueError('cache scientific input/source metadata differs from its record')
    if (metadata['degree']!=finest['degree'] or metadata['k']!=0 or
            not math.isclose(metadata['mass'],record['preparation']['mass'],rel_tol=2.e-14,abs_tol=2.e-14)):
        raise ValueError('cache preparation differs from the saved source record')
    n=finest['spatial_block_dimension']
    if (pair.shape!=(n,n) or eigenvalues.shape!=(n,) or not 0<=index<n or
            not np.all(np.isfinite(pair)) or not np.all(np.isfinite(eigenvalues)) or
            index!=finest['pair_matrix']['ground_column_index']):
        raise ValueError('matrix/label shape or ground index differs from its record')
    squared=np.abs(pair)**2
    column=squared[:,index];total=float(np.sum(column))
    if not math.isfinite(total) or total<=0:raise ValueError('ground-column squared norm must be positive')
    same=float(column[index]/total);other=1-same
    for actual,expected in ((total,finest['ground_column_squared_coefficient_norm']),
            (same,finest['same_mode_squared_coefficient_fraction']),
            (other,finest['other_mode_squared_coefficient_fraction'])):
        if not math.isclose(actual,expected,rel_tol=2.e-11,abs_tol=2.e-13):
            raise ValueError('matrix-derived plot quantity disagrees with the saved summary')
    return record,finest,cache,pair,eigenvalues,index,same,other,record_identity,cache_identity


def main():
    script_identity=identity(ROOT/PRODUCER)
    record,finest,cache,pair,eigenvalues,index,same,other,record_identity,cache_identity=load_saved()
    squared=np.abs(pair)**2
    floor=float(np.max(squared))*1.e-8
    style={'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':14,
           'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,
           'axes.spines.top':False,'axes.spines.right':False}
    with plt.rc_context(style):
        fig,axes=plt.subplots(1,2,figsize=(11.8,6.6),gridspec_kw={'width_ratios':[1.15,1.]})
        fig.subplots_adjust(left=.10,right=.97,bottom=.28,top=.79,wspace=.52)
        heat,bars=axes
        palette=plt.get_cmap('Blues').copy();palette.set_bad('#f4f6f8')
        image=heat.imshow(np.ma.masked_less_equal(squared,floor),origin='lower',
            cmap=palette,norm=LogNorm(vmin=floor,vmax=float(np.max(squared))),interpolation='nearest')
        positions=np.unique(np.r_[np.linspace(0,len(eigenvalues)-1,5,dtype=int),index])
        labels=[f'{eigenvalues[i]:+.2f}' for i in positions]
        heat.set(xticks=positions,yticks=positions,xticklabels=labels,yticklabels=labels,
                 xlabel='Negative initial-energy labels\nsigned spatial Dirac eigenvalue',
                 ylabel='Positive initial-energy labels\nsigned spatial Dirac eigenvalue',
                 title=f'A  Full k=0 matrix, cutoff L={finest["degree"]}')
        heat.add_patch(Rectangle((index-.5,-.5),1,len(eigenvalues),fill=False,
            edgecolor='#d18a17',linewidth=1.6))
        heat.plot(index,index,'o',ms=5,mfc='none',mec='#d18a17',mew=1.4)
        colorbar=fig.colorbar(image,ax=heat,fraction=.047,pad=.035)
        colorbar.ax.tick_params(labelsize=10)
        colorbar.set_label(r'$|M_{pn}|^2$  (log scale)',fontsize=11)
        values=100*np.array([same,other]);ys=np.array([1.,0.])
        bars.barh(ys,values,height=.42,color=['#245b8d','#739a94'])
        for y,value,label in zip(ys,values,('Same spatial mode','All other spatial modes')):
            bars.text(value+2,y,f'{value:.4f}%',ha='left',va='center',fontsize=14,fontweight='bold')
            bars.text(0,y+.30,label,ha='left',va='bottom',fontsize=12)
        bars.set(yticks=[],xlim=(0,103),
                 ylim=(-.55,2.),xticks=[0,25,50,75,100],
                 xlabel='Normalized squared coefficient (%)',title='B  Ground-column coefficient split')
        bars.grid(axis='x',alpha=.2);bars.set_axisbelow(True)
        for spine in ('left','top','right'):bars.spines[spine].set_visible(False)
        bars.tick_params(axis='y',length=0,pad=8)
        bars.text(0,1.92,f'{len(eigenvalues)} signed spatial labels\n'
            f'Column squared norm: {finest["ground_column_squared_coefficient_norm"]:.8f}',
            ha='left',va='top',fontsize=11,color='#4b5563')
        fig.suptitle('Initial Dirac coupling on the inherited parent',y=.965,fontsize=18,fontweight='bold')
        fig.text(.5,.90,'Squared reference-amplitude coefficients; no occupation is specified',
                 ha='center',va='center',fontsize=12,color='#4b5563')
        fig.text(.29,.125,'Gold outline: ground input column\nDisplay floor: 10⁻⁸ of the maximum',
                 ha='center',va='center',fontsize=10.5,color='#4b5563')
        fig.text(.5,.035,'Fixed initial energy basis, finite k=0 block. '
                 'The normalized split is not an energy fraction or probability.',
                 ha='center',va='center',fontsize=11)
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(OUTPUT,dpi=220,facecolor='white',metadata={'Software':f'Matplotlib {matplotlib.__version__}'})
        plt.close(fig)
    if identity(RECORD)!=record_identity or identity(cache)!=cache_identity or identity(ROOT/PRODUCER)!=script_identity:
        raise RuntimeError('record, matrix cache or plot producer changed during rendering')
    manifest=dict(schema='nsc-atlas-dirac-source-figure-v1',producer=script_identity,
        inputs={'scientific_record':record_identity,'raw_matrices':cache_identity},
        outputs={OUTPUT.relative_to(ROOT).as_posix():identity(OUTPUT)},
        render=dict(dpi=220,figsize_inches=[11.8,6.6],matplotlib_version=matplotlib.__version__,
                    quantity='Squared magnitudes; arbitrary eigenvector phases do not change the figure.',
                    heatmap_floor_relative_to_maximum=1.e-8),
        data=dict(degree=finest['degree'],k=0,spatial_labels=len(eigenvalues),ground_column_index=index,
            same_mode_squared_coefficient_fraction=same,other_mode_squared_coefficient_fraction=other,
            ground_column_squared_coefficient_norm=finest['ground_column_squared_coefficient_norm']),
        scientific_scope='Reference-amplitude coefficient norms of the saved initial Hamiltonian Hessian; not energy fractions, probabilities, occupation or a renormalized source.',
        scientific_sources_or_PDE_recomputed=False)
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(figure=str(OUTPUT.relative_to(ROOT)),manifest=str(MANIFEST.relative_to(ROOT)),
        same_percent=100*same,other_percent=100*other,figure_sha256=identity(OUTPUT)['sha256']),indent=2))


if __name__=='__main__':main()
