"""Source-labelled route map and actual finite flow/work/flux readouts."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
ROOT=Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

fig,ax=plt.subplots(figsize=(10,6.7));ax.set(xlim=(0,10),ylim=(0,6.7));ax.axis('off')
boxes=[(0.2,4.55,'SPHERICAL BALANCE','Fixed volume, isotropic tension\nA cubed >= 36 pi V squared'),
       (3.55,4.55,'PREPARED CIRCULATION','Mean velocity can vanish\nCovariance and energy remain'),
       (6.9,4.55,'BOUNDARY GEOMETRY','Five-dimensional fluid/gravity\nComputed energy + curvature'),
       (0.2,1.85,'INHERITED FIELD','One real higher-dimensional field\nCharged harmonic pair'),
       (3.55,1.85,'HORIZON TRANSFER','Offset: frequency - m angular speed\nEnergy and angular-momentum flux'),
       (6.9,1.85,'EXPANDING CHILD','Regular interior, shared source\nTransition still to be constructed')]
for x,y,title,body in boxes:
    ax.add_patch(FancyBboxPatch((x,y),2.9,1.25,boxstyle='round,pad=0.09',facecolor='#EDF3F9',edgecolor='#456C96',lw=1.2))
    ax.text(x+1.45,y+.92,title,ha='center',va='center',fontsize=10,fontweight='bold',color='#173E70')
    ax.text(x+1.45,y+.4,body,ha='center',va='center',fontsize=8.7,linespacing=1.5)

def arrow(a,b,label,solid=True,offset=(0,0)):
    color='#245B9C' if solid else '#A76E19'
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,lw=1.7,color=color,linestyle='-' if solid else '--'))
    ax.text((a[0]+b[0])/2+offset[0],(a[1]+b[1])/2+offset[1],label,ha='center',va='center',fontsize=8.5,color=color)
arrow((3.15,5.15),(3.46,5.15),'prepare',False,(0,.48))
arrow((6.5,5.15),(6.81,5.15),'lift',True,(0,.48))
arrow((3.15,2.45),(3.46,2.45),'derive',True,(0,.48))
arrow((6.5,2.45),(6.81,2.45),'match',False,(0,.48))
arrow((8.35,4.43),(8.35,3.2),'Source, dimension and\nobserver matching',False,(-1.23,0))
ax.text(3.3,3.75,'Calculated averaged boundary work:\nb2 = 2(E0-E)/(3 V nu)',ha='center',va='center',fontsize=9,color='#245B9C')
ax.text(5,6.4,'The atlas connects equations across regions and resolutions',ha='center',fontsize=14,fontweight='bold')
ax.text(5,.95,'Blue: calculated connection in the stated domain.  Dashed: required physical attachment.',ha='center',fontsize=9)
ax.text(5,.45,'Existing parent-source chapters separately connect fields, finite evolution, clocks and observations.\nThis map does not identify every branch as one completed cosmic history.',ha='center',fontsize=9,color='#4D5968')
fig.savefig(ROOT/'paper/figures/nsc-atlas-connected-route.png',dpi=160,bbox_inches='tight',metadata={'Software':'NSC connected route'});plt.close(fig)

flow=json.loads((ROOT/'results/summaries/nsc-atlas-flow-evolution-v1.json').read_bytes())
work=json.loads((ROOT/'results/summaries/nsc-atlas-boundary-work-v1.json').read_bytes())
fine=flow['cases'][-1];history=fine['history'];w=work['cases'][-1]['history']
fig,axes=plt.subplots(2,2,figsize=(9,7.6),constrained_layout=True)
times=[r['t'] for r in history]
ax=axes[0,0]
for i,(label,color) in enumerate((('inner radius .75','#245B9C'),('outer radius 1.25','#C27517'))):
    ax.plot(times,[r['mean_theta'][i]-history[0]['mean_theta'][i] for r in history],'-o',ms=3,color=color,label=label)
ax.axhline(0,color='#888888',lw=.8);ax.set(title='A  Circulation develops from fine motion',xlabel='Slow fluid time',ylabel='Generated mean tangential velocity');ax.legend(fontsize=8);ax.grid(alpha=.15)
ax=axes[0,1]
ax.plot(times,[2*r['energy']/.001**1.5 for r in history],'-o',ms=3,color='#245B9C')
ax.set(title='B  Same state: gravitational boundary energy',xlabel='Slow fluid time',ylabel='Integrated coefficient of epsilon squared')
ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0));ax.grid(alpha=.15)
ax=axes[1,0]
ax.plot([r['t'] for r in w],[r['b2'] for r in w],'-o',ms=3,color='#245B9C')
ax.set(title='C  Energy loss fixes averaged boundary work',xlabel='Slow fluid time',ylabel='Boundary scale coefficient b2')
ax.text(.04,.82,'b = 1 + epsilon squared x b2\nAveraged Ward identity; bulk closure remains.',transform=ax.transAxes,fontsize=8.4)
ax.grid(alpha=.15)
ax=axes[1,1];x=np.linspace(.02,2,250)
ax.plot(x,x*(x-1),color='#245B9C',label='Energy into horizon')
ax.plot(x,(x-1)**2,color='#C27517',ls='--',label='Co-rotating energy (nonnegative)')
ax.axhline(0,color='#888888',lw=.8);ax.axvline(1,color='#888888',lw=.8,ls=':')
ax.set(title='D  Rotating horizon: a frequency threshold',xlabel='Mode frequency / horizon angular speed',ylabel='Flux / (amplitude norm x angular speed squared)')
ax.legend(fontsize=8);ax.grid(alpha=.15)
fig.suptitle('Computed connections: circulation, geometry, work and horizon transfer\nFinite periodic flow; hydrodynamic coefficients; separately specified horizon modes',fontsize=12)
fig.savefig(ROOT/'paper/figures/nsc-atlas-flow-gravity.png',dpi=160,metadata={'Software':'NSC flow gravity'});plt.close(fig)
print('Wrote connected route and finite flow/gravity figures.')
