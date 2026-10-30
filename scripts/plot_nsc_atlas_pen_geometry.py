"""Plot the constructed sphere, finite resolution, dipole and probe readout."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from recursive_horizons.nsc_atlas_pen_geometry import angular_gaps,phase_residual,orientations

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'results/summaries/nsc-atlas-pen-geometry-v1.json').read_text())
omega=r['construction']['omega'];delta=r['construction']['delta']
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(9,8.4),constrained_layout=True)
ax=fig.add_subplot(221,projection='3d')
theta=np.linspace(0,np.pi,90);phi=np.linspace(0,2*np.pi,181)
t,p=np.meshgrid(theta,phi,indexing='ij')
x,y,z=np.sin(t)*np.cos(p),np.sin(t)*np.sin(p),np.cos(t)
paint=np.zeros_like(x,dtype=bool)
for k,v in enumerate(orientations(18,delta)):
    paint|=np.abs(-np.sin(v)*x+np.cos(v)*y)<=np.sin(omega**(-k-1))
colors=np.empty(x.shape+(4,));colors[:]=matplotlib.colors.to_rgba('#D5E2F0')
colors[paint]=matplotlib.colors.to_rgba('#D88B2A')
ax.plot_surface(x,y,z,facecolors=colors,rstride=1,cstride=1,linewidth=0,shade=True)
u=np.linspace(0,2*np.pi,280)
for v in orientations(8,delta):
    ax.plot(1.004*np.sin(u)*np.cos(v),1.004*np.sin(u)*np.sin(v),1.004*np.cos(u),color='#173E70',lw=.55,alpha=.65)
ax.set_box_aspect((1,1,1));ax.view_init(23,40);ax.set_axis_off()
coverage=r['persistent_coverage'][0]
ax.set_title(f"A  Persistent drawing\nDrawn: {coverage['area_fraction']:.5f}; undrawn: {coverage['undrawn_fraction']:.5f}")

ax=fig.add_subplot(222)
n=np.arange(1,201)
actual=np.array([-np.log(max(angular_gaps(int(k),delta))/2) for k in n])
ax.plot(n,n*np.log(omega),label='Assigned log resolution: n ln(Omega)',color='#2860A4')
ax.plot(n,n*np.log(omega)/2,label='Assigned: n ln(Omega)/2',color='#7096C5',ls='--')
ax.plot(n,actual,label='Achieved uniform angular resolution',color='#C27517',lw=2)
ax.axhspan(136.84,138.26,color='#758276',alpha=.3,label='Requested band')
ax.set(xlabel='Number of circles drawn',ylabel='Log resolving ratio',ylim=(-5,285),title='B  Assigned refinement\nversus actual coverage')
ax.legend(fontsize=8,loc='upper left');ax.grid(alpha=.16)

ax=fig.add_subplot(223)
res=np.array([phase_residual(int(k),delta) for k in n])
ax.semilogy(n,res,color='#2860A4',lw=1)
ax.axhline(6.1e-10,color='#C27517',ls='--',label='Requested eta benchmark')
ax.scatter([100,200],[phase_residual(100,delta),phase_residual(200,delta)],color='#173E70')
ax.text(.04,.15,'Scalar signed charge stays exactly zero.\nCurve is a geometric dipole residual.',transform=ax.transAxes,fontsize=9)
ax.set(xlabel='Pass count',ylabel='Normalized dipole magnitude',title='C  Irrational offset:\nan endpoint dipole',ylim=(1e-10,2))
ax.legend(fontsize=8,loc='upper right');ax.grid(alpha=.16)

ax=fig.add_subplot(224)
probe=r['probe'];vals=[1.,probe['ratio_depth_2_over_1'],probe['depth_100_over_1']]
ax.bar(['1 pass','2 passes','100 passes'],vals,color=['#9EB8D7','#5D87BA','#245B9C'],width=.55)
ax.axhline(1,color='#646C75',lw=.8)
for i,val in enumerate(vals):ax.text(i,val+.004,f'{val:.5f}',ha='center',fontsize=10)
ax.set(ylim=(.97,1.13),ylabel='Resolved coupling / first-pass coupling',title='D  One quadratic probe\non the same drawn set')
ax.text(.04,.03,'Common bare coupling cancels.\nThis is a geometric probe, not a force identification.',transform=ax.transAxes,fontsize=9)
ax.grid(axis='y',alpha=.15)
fig.suptitle(f"Pen map from the late-vacuum mode clocks: delta = {delta:.9f}\nFixed source omega0 = 20; declared footprint; no target fitting",fontsize=12)
out=ROOT/'paper/figures/nsc-atlas-pen-geometry.png'
fig.savefig(out,dpi=170,metadata={'Software':'NSC pen geometry plot'});plt.close(fig)
print(out)
