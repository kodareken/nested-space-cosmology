"""Plot the admitted fixed-charge transition record and authenticated PPoly.

No physical decoder, field solver, raw Lorentzian history or new state is used.
The two pole endpoints are the saved regular-boundary data, not extrapolation.
"""
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PPoly

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'results/summaries/nsc-atlas-shape-transition-v1.json'
PRODUCER=ROOT/'scripts/plot_nsc_atlas_shape_transition.py'
OUTPUT=ROOT/'docs/figures/nsc-atlas-shape-transition.png'
MANIFEST=OUTPUT.with_suffix('.manifest.json')


def identity(path):
    raw=Path(path).read_bytes()
    return dict(path=Path(path).resolve().relative_to(ROOT).as_posix(),sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))


def load_saved():
    binding=identity(SOURCE);record=json.loads(SOURCE.read_text(encoding='utf-8'))
    if (record.get('schema')!='nsc-atlas-shape-transition-v1' or
        record['provenance'].get('source_unchanged') is not True or
        record['provenance'].get('source_bound_to_commit') is not True):
        raise ValueError('unchanged historically source-bound transition record required')
    transition=record['source_free_transition'];fine=transition['fine'];chart=transition['connecting_chart']
    if transition['q']!=1.5 or transition['membrane_inserted'] or transition['charge_jump']!=0:
        raise ValueError('the plotted candidate must be source-free at fixed charge3/2')
    raw=(ROOT/record['provenance']['raw_cache']).resolve()
    if not raw.is_relative_to(ROOT/'.build') or not raw.is_file():
        raise FileNotFoundError('saved transition PPoly cache required; no solver was started')
    raw_identity=identity(raw)
    if raw_identity['sha256']!=record['provenance']['raw_cache_sha256']:
        raise ValueError('transition PPoly cache SHA256 mismatch')
    with np.load(raw,allow_pickle=False) as cache:
        knots=cache['fine_knots'].copy();coefficients=cache['fine_coefficients'].copy()
        length=float(cache['fine_length']);cutoff=float(cache['fine_cutoff'])
    if (knots.ndim!=1 or len(knots)<2 or np.any(np.diff(knots)<=0) or
        coefficients.shape!=(4,len(knots)-1,6) or
        not np.all(np.isfinite(knots)) or not np.all(np.isfinite(coefficients))):
        raise ValueError('saved six-field cubic PPoly shape is invalid')
    if (not math.isclose(length,fine['length'],rel_tol=2e-14) or cutoff!=fine['cutoff'] or
        not math.isclose(knots[0],cutoff,abs_tol=2e-14) or not math.isclose(knots[-1],1-cutoff,abs_tol=2e-14)):
        raise ValueError('PPoly normalized domain differs from the public record')
    profile=PPoly(coefficients,knots,extrapolate=False)
    x=np.linspace(knots[0],knots[-1],1601);values=profile(x)
    pole_fields=np.asarray(fine['pole_fields'],float)
    if pole_fields.shape!=(2,2):raise ValueError('north/south pair-field values required')
    if not np.allclose(profile([knots[0],knots[-1]])[:,:2],pole_fields,rtol=0,atol=2e-6):
        raise ValueError('PPoly pole-cutoff fields disagree with the saved regular pole values')
    equator_radius=float(chart['equator_radius'])
    if not math.isclose(float(profile(chart['equator_x'])[4]),equator_radius,rel_tol=3e-13,abs_tol=2e-12):
        raise ValueError('published equator radius differs from the authenticated profile')
    if chart['equator_spatial_second_derivative']>=0:
        raise ValueError('the plotted equator must be the recorded spatial maximum')
    # Append the saved regular pole limits explicitly; do not extrapolate
    # the cutoff PPoly outside its authenticated knot interval.
    normalized=np.r_[0.,x,1.]
    euclidean_radius=np.r_[0.,values[:,4],0.]
    fields=np.vstack((pole_fields[0],values[:,:2],pole_fields[1]))
    table=transition['interior_table'];columns=table['columns'];rows=np.asarray(table['rows'],float)
    interior_t=rows[:,columns.index('proper_time')];scale=rows[:,columns.index('scale')]
    inside=transition['interior']
    if (inside['endpoint_reason']!='first_scale_factor_turn' or
        not math.isclose(interior_t[-1],inside['proper_endpoint'],rel_tol=2e-14) or
        not math.isclose(scale[-1],inside['scale_at_endpoint'],rel_tol=2e-14)):
        raise ValueError('published interior table must end at its first scale-factor turn')
    interior_t=np.r_[0.,interior_t];scale=np.r_[0.,scale]
    T=np.linspace(-1.,1.,201);theta=2*np.tanh(T)/equator_radius
    for key,pair in chart['equator_null_expansions'].items():
        expected=2*math.tanh(float(key))/equator_radius
        if not all(math.isclose(v,expected,rel_tol=2e-14,abs_tol=2e-14) for v in pair):
            raise ValueError('saved null-screen normalization differs from the plot formula')
    if any(not np.all(np.isfinite(a)) for a in (normalized,euclidean_radius,fields,interior_t,scale,theta)):
        raise ValueError('all plotted profile/table/formula values must be finite')
    return dict(record=record,binding=binding,raw=raw,raw_identity=raw_identity,length=length,
        normalized=normalized,radius=euclidean_radius,fields=fields,pole_fields=pole_fields,
        interior_t=interior_t,scale=scale,inside=inside,T=T,theta=theta,equator_radius=equator_radius)


def main():
    producer=identity(PRODUCER);data=load_saved()
    blue='#245b80';orange='#b56326';purple='#71508f';muted='#586471'
    style={'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':12.3,
        'axes.labelsize':10.5,'xtick.labelsize':10,'ytick.labelsize':10,
        'axes.spines.top':False,'axes.spines.right':False}
    with plt.rc_context(style):
        fig,axes=plt.subplots(2,2,figsize=(12.8,8.0))
        fig.subplots_adjust(left=.085,right=.975,bottom=.165,top=.835,wspace=.29,hspace=.64)
        a,b,c,d=axes.ravel()
        x=data['normalized']
        a.plot(x,data['radius'],color=blue,lw=2.1)
        a.plot([0.,1.],[0.,0.],'o',color=blue,ms=4)
        a.set(title='A  Two-pole Euclidean profile',xlabel=r'Normalized coordinate $\xi/L$',
            ylabel=r'$R_E$ (model length units)',xlim=(0,1))
        a.text(.04,.89,rf'$L={data["length"]:.6f}$',transform=a.transAxes,fontsize=11)
        a.text(.045,.055,'N',transform=a.transAxes,fontsize=10,color=muted)
        a.text(.955,.055,'S',transform=a.transAxes,ha='right',fontsize=10,color=muted)
        a.text(.50,.075,'Two regular poles',transform=a.transAxes,ha='center',fontsize=9.5,color=muted)
        b.plot(x,data['fields'][:,0],color=blue,lw=2.,label=r'$u_{\rm single}$')
        b.plot(x,data['fields'][:,1],color=orange,lw=2.,ls='--',label=r'$u_{\rm other}=u_2=u_3=u_4$')
        b.set(title='B  Internal log-radius fields',xlabel=r'Normalized coordinate $\xi/L$',
            ylabel='Internal log radii u (dimensionless)',xlim=(0,1))
        b.legend(loc='upper right',frameon=False,fontsize=10)
        north,south=data['pole_fields']
        b.text(.5,-.29,f'Pole values (single, other): N=({north[0]:+.5f}, {north[1]:+.5f})\n'
            f'S=({south[0]:+.5f}, {south[1]:+.5f}); other represents three equal fields.',
            transform=b.transAxes,ha='center',va='top',fontsize=9.2,color=muted)
        c.plot(data['interior_t'],data['scale'],color=purple,lw=2.1)
        c.plot(data['interior_t'][-1],data['scale'][-1],'o',color=purple,ms=5)
        c.set(title='C  Open-region scale to its first turn',xlabel='Interior proper time t (model units)',
            ylabel='Scale factor a (model length units)',xlim=(0,float(data['interior_t'][-1])))
        c.text(.04,.94,f'Exit density = {data["inside"]["exit_density"]:.7f}\n'
            rf'Minimum $\Omega_k$ = {data["inside"]["minimum_curvature_fraction"]:.6f}',
            transform=c.transAxes,va='top',fontsize=10,color=muted)
        c.text(.97,.34,f'First turn: t={data["interior_t"][-1]:.5f}',
            transform=c.transAxes,ha='right',fontsize=10,color=muted)
        c.text(.97,.07,'This is not an observed cosmology.',transform=c.transAxes,ha='right',fontsize=9.5,color=muted)
        d.plot(data['T'],data['theta'],color=blue,lw=2.5,label=r'$\theta_+$')
        d.plot(data['T'],data['theta'],color=orange,lw=1.5,ls='--',label=r'$\theta_-$ (identical)')
        d.axhline(0.,color='#737e88',lw=.8);d.axvline(0.,color='#737e88',lw=.8,ls=':')
        d.plot(0.,0.,'o',color='#37414c',ms=4)
        d.set(title='D  Equator null-screen expansions',xlabel='Equator chart time T (dimensionless)',
            ylabel=r'$\theta_\pm$ (model inverse-length units)',xlim=(-1,1))
        d.legend(loc='upper left',frameon=False,fontsize=10)
        d.text(.04,.62,'Regular marginal cut at T=0',transform=d.transAxes,color=muted,fontsize=9.5)
        d.text(.97,.08,r'Spatial maximum of $R_E$;'+'\nnot a persistent neck.',transform=d.transAxes,
            ha='right',va='bottom',color=muted,fontsize=9.5)
        for ax in axes.ravel():ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)
        fig.suptitle('A regular shape transition and its Lorentzian continuations',y=.966,fontsize=16.5,fontweight='bold')
        fig.text(.5,.920,'Saved fixed-charge Euclidean saddle and derived Lorentzian charts (model units)',
            ha='center',fontsize=11,color=muted)
        fig.text(.5,.077,'Source-free at fixed q=3/2. Internal radii are fields; spatial dimension and SI scale are not selected.',
            ha='center',fontsize=10.5,color='#37414c')
        fig.text(.5,.045,'Full physical negative-mode index and incoming parent flow are unproved; no black-hole trigger or present-day density fraction is inferred.',
            ha='center',fontsize=10,color='#37414c')
        OUTPUT.parent.mkdir(parents=True,exist_ok=True)
        fig.savefig(OUTPUT,dpi=125,facecolor='white',metadata={'Software':f'Matplotlib {matplotlib.__version__}'})
        plt.close(fig)
    if identity(SOURCE)!=data['binding'] or identity(data['raw'])!=data['raw_identity'] or identity(PRODUCER)!=producer:
        raise RuntimeError('public record, authenticated PPoly or producer changed while rendering')
    manifest=dict(schema='nsc-atlas-shape-transition-figure-v1',producer=producer,
        inputs=dict(scientific_record=data['binding'],raw_profile=data['raw_identity']),
        historical_source_commit=data['record']['provenance']['producer_commit'],output=identity(OUTPUT),
        render=dict(size_pixels=[1600,1000],dpi=125,matplotlib_version=matplotlib.__version__),
        plotted=dict(Euclidean_length=data['length'],pole_fields=data['pole_fields'].tolist(),
            PPoly_keys=['fine_knots','fine_coefficients','fine_length','fine_cutoff'],
            profile_endpoint_rule='Saved regular pole limits explicitly appended; PPoly is not extrapolated outside its knots.',
            interior_points=len(data['interior_t']),interior_first_turn=float(data['interior_t'][-1]),
            exit_density=data['inside']['exit_density'],minimum_interior_Omega_k=data['inside']['minimum_curvature_fraction'],
            equator_radius=data['equator_radius'],equator_formula='theta_plus=theta_minus=2*tanh(T)/R_equator',
            equator_T_domain=[-1.,1.]),
        scope='One numerically admitted source-free fixed-charge saddle, its derived fields and chart continuations. No physical negative-index, incoming material flow, black-hole trigger, SI or observed-cosmology claim.',
        current_physical_decoder_used=False,physical_source_or_PDE_recomputed=False,
        reproduction='python scripts/run_python.py scripts/plot_nsc_atlas_shape_transition.py')
    MANIFEST.write_text(json.dumps(manifest,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(figure=manifest['output']['path'],bytes=manifest['output']['bytes'],sha256=manifest['output']['sha256'],
        manifest=str(MANIFEST.relative_to(ROOT))),indent=2))


if __name__=='__main__':main()
