"""Plot shell limit, regular horizon, phase charge and horizon response."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from recursive_horizons.nsc_atlas_parke_nucleation import thin_wall
from recursive_horizons.nsc_atlas_smooth_horizon import metric_f,horizons

ROOT=Path(__file__).resolve().parents[1]
record=json.loads((ROOT/'results/summaries/nsc-atlas-smooth-horizon-v1.json').read_text())
lam=record['admitted_vacuum']['Lambda'];rc=record['shell_limit']['reference_radius']
mass=record['radius_inversion_control']['required_mass']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,2,figsize=(9,8),constrained_layout=True)
ax=axes[0,0];s=np.geomspace(.0001,.2,120)
ax.loglog(s,[thin_wall(sigma=float(v))['radius'] for v in s],color='#2860A4',label='Parke shell, fixed contrast')
ax.axhline(record['admitted_vacuum']['horizon_m0'],color='#C27517',ls='--',label='Same-action de Sitter horizon')
ax.scatter([.2],[rc],color='#2860A4')
ax.set(xlabel='Surface tension (model units)',ylabel='Radius (model units)',title='A  Removing the shell')
ax.legend(fontsize=8,loc='upper left');ax.grid(alpha=.15)
ax=axes[0,1];r=np.linspace(.55,4.05,220)
ax.plot(r,[metric_f(float(v),mass,lam) for v in r],color='#2860A4',label='Kottler: mass inverted to put root at old Rc')
ax.plot(r,[metric_f(float(v),0,lam) for v in r],color='#C27517',ls='--',label='de Sitter: no central mass')
ax.axhline(0,color='#5B626C',lw=.8)
ax.text(.04,.08,'Metric determinant stays nonzero at f = 0.\nFor E = 1: dr/dtau = -1, dv/dtau = 1/2.',transform=ax.transAxes,fontsize=9)
ax.set(xlabel='Areal radius',ylabel='Static metric function f',title='B  A zero with a smooth continuation',ylim=(-1.25,1.1))
ax.legend(fontsize=7.7,loc='upper right');ax.grid(alpha=.15)
ax=axes[1,0];t=np.linspace(0,2*np.pi,450)
for row,color in zip(record['offset_bridge']['controls'],['#2860A4','#C27517']):
    delta=row['delta'];label='Action phase' if row['phase_source']=='retained_action_phase' else 'Golden angle'
    ax.plot(np.cos(t),np.cos(t-delta),color=color,label=f"{label}: Q/(E/omega) = {np.sin(delta):.4f}")
ax.plot([-1,1],[-1,1],color='#777777',ls=':',label='Zero phase: zero oriented charge')
ax.set(xlabel='Oscillator component x / A',ylabel='Oscillator component y / B',title='C  Phase offset as oriented charge',aspect='equal')
ax.legend(fontsize=7,loc='lower right');ax.grid(alpha=.15)
ax=axes[1,1];mu=np.linspace(.003,.332,140)
roots=[horizons(float(v),1.) for v in mu]
ax.plot(mu,[x['black_hole'] for x in roots],color='#2860A4',label='Black-hole root')
ax.plot(mu,[x['cosmological'] for x in roots],color='#C27517',label='Cosmological root')
ax.axhline(1,color='#777777',ls=':',lw=.8)
ax.annotate('decreasing mass',xy=(.09,1.05),xytext=(.25,1.05),arrowprops={'arrowstyle':'->'},ha='center',fontsize=9)
ax.set(xlabel='Mass parameter times sqrt(Lambda)',ylabel='Horizon radius times sqrt(Lambda)',title='D  One mass change, opposite responses')
ax.legend(fontsize=8,loc='lower right');ax.grid(alpha=.15)
fig.suptitle('From a material shell to a smooth horizon and a charge current\nExact conditional constructions; no selected drain rate or observed asymmetry',fontsize=12)
out=ROOT/'paper/figures/nsc-atlas-smooth-horizon.png'
fig.savefig(out,dpi=160,metadata={'Software':'NSC smooth horizon plot'});plt.close(fig)
print(out)
