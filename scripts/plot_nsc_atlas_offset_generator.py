"""Visualize the explicit one-seed map, finite coverage and metric measure."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from recursive_horizons.nsc_atlas_offset_generator import phase_history,directions
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-offset-generator-v1.json').read_bytes())
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(9,8),constrained_layout=True)
phases=phase_history(100,'action');v=directions(phases)
ax=fig.add_subplot(221)
for theta in 2*np.pi*phases[:32,0]:
    ax.plot([-np.cos(theta),np.cos(theta)],[-np.sin(theta),np.sin(theta)],color='#245B9C',alpha=.4,lw=.8)
t=np.linspace(0,2*np.pi,300);ax.plot(np.cos(t),np.sin(t),color='#C27517',lw=1.5)
ax.set(aspect='equal',title='A  First accumulation\nDisk completion; finite lines',xlabel='Derived planar coordinate',ylabel='Derived planar coordinate')
ax=fig.add_subplot(222,projection='3d')
for d in v[:65]:ax.plot([-d[0],d[0]],[-d[1],d[1]],[-d[2],d[2]],color='#245B9C',alpha=.22,lw=.55)
u=np.linspace(0,np.pi,15);p=np.linspace(0,2*np.pi,29)
U,P=np.meshgrid(u,p,indexing='ij')
ax.plot_wireframe(np.sin(U)*np.cos(P),np.sin(U)*np.sin(P),np.cos(U),color='#8094AB',alpha=.3,lw=.5)
ax.scatter(v[:,0],v[:,1],v[:,2],s=5,color='#C27517',alpha=.75)
ax.set_box_aspect((1,1,1),zoom=1.35);ax.view_init(22,35);ax.set_axis_off()
ax.set_xlim(-1.05,1.05);ax.set_ylim(-1.05,1.05);ax.set_zlim(-1.05,1.05)
ax.set_title('B  Second accumulation\nBall completion; finite spokes')
ax=fig.add_subplot(223)
for case,color in zip(r['cases'],['#245B9C','#C27517']):
    ax.loglog([x['passes'] for x in case['coverage']],[x['ball_Hausdorff_distance_over_radius'] for x in case['coverage']],'-o',color=color,label=case['offset']+' offset')
ax.set(title='C  Finite resolution improves',xlabel='Retained passes',ylabel='Maximum ball distance / radius')
ax.legend(fontsize=9);ax.grid(alpha=.2,which='both')
ax=fig.add_subplot(224);row=r['cases'][0]['coverage'][-1];x=np.arange(3)
ax.bar(x-.18,row['raw_direction_second_moment'],.34,color='#9EB8D7',label='Uniform pass count')
ax.bar(x+.18,row['area_weighted_direction_second_moment'],.34,color='#245B9C',label='Metric area weight')
ax.axhline(1/3,color='#C27517',ls='--',lw=1)
ax.set(xticks=x,xticklabels=['x direction','y direction','z direction'],ylim=(0,.56),title='D  Metric area versus pass count',ylabel='Direction second moment')
ax.legend(fontsize=8.7,loc='upper left');ax.grid(axis='y',alpha=.15)
fig.suptitle('One spatial seed + accumulated offset history\nMemory, distance kernel and completion are explicit construction rules',fontsize=13)
fig.savefig(ROOT/'paper/figures/nsc-atlas-offset-generator.png',dpi=170,metadata={'Software':'NSC offset generator'});plt.close(fig)
print('Wrote offset-generator figure.')
