"""Batched independent real-frequency modes of the existing PG Dirac owner.

The same matched collar, massive Jost ratio, current normalization and PG
frames are retained. Vectorization adds no coupling between frequencies.
"""
import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits

from .nsc_massive_jost_modes import outgoing_ratio
from .nsc_paired_horizon_preparation import source_covariance
from .nsc_unruh_state import horizon_frame
from .nsc_pg_massive_modes import MassivePGModeResolution
from .nsc_lorentzian import geometry


def _metric_A(delta,background):
    delta=np.asarray(delta,float);rho=background.horizon_rho+delta;A=np.empty_like(delta)
    near=delta<1e-4;far=rho>20;middle=~(near|far)
    a,b,c,d=background._near_coefficients()
    A[near]=a*delta[near]+b*delta[near]**2/2+c*delta[near]**3/6+d*delta[near]**4/24
    A[far]=1-sum(6*(-1)**n/((2*n+1)*(2*n+3)*rho[far]**(2*n+1)) for n in range(6))
    A[middle]=geometry(rho[middle])[2]
    return A


def real_mode_family(energies,mass,angular,preparation,*,rtol=3e-14,atol=3e-16):
    E=np.asarray(energies,float);n=len(E);p=preparation;b=p.background
    if E.ndim!=1 or np.any(E<=0) or np.any(E==mass):raise ValueError('positive real frequencies away from threshold required')
    open_channel=E>mass
    end=[];ratios=[];initial_current=[]
    for energy in E:
        k=np.sqrt(complex(energy*energy-mass*mass))
        radius=max(60.,30*max(1.,abs(angular),mass)/max(abs(k),.2))
        for _ in range(16):
            r,dr,_=outgoing_ratio(energy,mass,angular,radius)
            A=b.A_from_offset(radius-b.horizon_rho);v1=angular*np.sqrt(A)/np.sqrt(1+radius*radius);v2=mass*np.sqrt(A)
            residual=abs(A*dr-1j*(v1+1j*v2)+2j*energy*r-1j*(v1-1j*v2)*r*r)
            current=1-abs(r)**2
            if residual<3e-13 and (energy>mass or abs(current)<3e-13):break
            radius*=2
        else:raise ArithmeticError('batched Jost outer expansion unresolved')
        end.append(radius);ratios.append(r);initial_current.append(current)
    end=np.array(end);ratios=np.array(ratios);initial_current=np.array(initial_current)
    # Closed channels evolve a phase to preserve their zero current.
    state=np.stack((np.where(open_channel,ratios,np.angle(ratios)),np.zeros(n,complex)),axis=1)
    def jost_rhs(y,U,dy):
        delta=np.exp(y);rho=b.horizon_rho+delta;A=_metric_A(delta,b)
        W1=angular*np.sqrt(A)/np.sqrt(1+rho*rho);W2=mass*np.sqrt(A)
        ratio=np.where(open_channel,U[:,0],np.exp(1j*U[:,0].real))
        dr=delta/A*(1j*(W1+1j*W2)-2j*E*ratio+1j*(W1-1j*W2)*ratio*ratio)
        first=np.where(open_channel,dr,(dr/(1j*ratio)).real)
        second=1j*delta/A*(E-(W1-1j*W2)*ratio)
        return np.stack((first,second),axis=1)*np.asarray(dy)[:,None]
    y0=np.log(end-b.horizon_rho);y3=np.log(3.-b.horizon_rho);yh=np.log(p.collar_delta_q)
    with threadpool_limits(limits=1):
        def first_stage(t,u):return jost_rhs(y0+t*(y3-y0),u.reshape(n,2),y3-y0).ravel()
        run=solve_ivp(first_stage,(0.,1.),state.ravel(),method='DOP853',rtol=rtol,atol=atol,t_eval=[1.])
        if not run.success:raise ArithmeticError(run.message)
        at3=run.y[:,-1].reshape(n,2)
        def second_stage(y,u):return jost_rhs(np.full(n,y),u.reshape(n,2),np.ones(n)).ravel()
        run=solve_ivp(second_stage,(y3,yh),at3.ravel(),method='DOP853',rtol=rtol,atol=atol,t_eval=[yh])
        if not run.success:raise ArithmeticError(run.message)
        atH=run.y[:,-1].reshape(n,2)
    rH=np.where(open_channel,atH[:,0],np.exp(1j*atH[:,0].real))
    r3=np.where(open_channel,at3[:,0],np.exp(1j*at3[:,0].real))
    rh=p.horizon_radius;phase=np.arctan2(mass*rh,angular);orient=np.diag(np.exp(np.array([-1j,1j])*phase/2))
    frames=np.array([orient@horizon_frame(e,p.surface_gravity,np.hypot(angular,mass*rh),
               np.sqrt(2*p.collar_delta_q/p.surface_gravity)/rh,b.near_tortoise(p.collar_delta_q)) for e in E])
    coefficients=np.linalg.solve(frames,np.stack((np.ones(n),rH),axis=1)[...,None])[...,0]
    R=coefficients[:,1]/coefficients[:,0]
    T=np.where(open_channel,initial_current*np.exp(-2*atH[:,1].real)/abs(coefficients[:,0])**2,0.)
    current_error=float(np.max(abs(abs(R)**2+T-1)))
    if np.any(T<0) or current_error>3e-11:raise ArithmeticError(f'batched current error {current_error}')
    up3=np.exp(at3[:,1]-atH[:,1])/coefficients[:,0]
    up3=up3[:,None]*np.stack((np.ones(n),r3),axis=1)
    initial_in=frames[:,:,1]*np.sqrt(T)[:,None]
    initial_ks=np.array([p.working_frame(yh,e,mass,angular) for e in E])
    with threadpool_limits(limits=1):
        def incoming_rhs(y,u):
            delta=np.exp(y);rho=b.horizon_rho+delta;A=b.A_from_offset(delta)
            W1=angular*np.sqrt(A)/np.sqrt(1+rho*rho);W2=mass*np.sqrt(A);f=u.reshape(n,2)
            return (1j*delta/A*np.stack((E*f[:,0]-(W1-1j*W2)*f[:,1],(W1+1j*W2)*f[:,0]-E*f[:,1]),axis=1)).ravel()
        run=solve_ivp(incoming_rhs,(yh,y3),initial_in.ravel(),method='DOP853',rtol=rtol,atol=atol,t_eval=[y3])
        if not run.success:raise ArithmeticError(run.message)
        incoming3=run.y[:,-1].reshape(n,2)
        target=np.log(b.horizon_q-np.pi/2)
        def inner_rhs(y,u):
            delta=np.exp(y);W=b.interior_W(delta);radius=1/np.sin(b.horizon_q-delta)
            hx=-mass*radius*delta/np.sqrt(W);hy=angular*delta/np.sqrt(W);hz=-E*delta/W
            U=u.reshape(n,2,2)
            return (-1j*np.stack((hz[:,None]*U[:,0]+(hx-1j*hy)*U[:,1],(hx+1j*hy)*U[:,0]-hz[:,None]*U[:,1]),axis=1)).ravel()
        run=solve_ivp(inner_rhs,(yh,target),initial_ks.ravel(),method='DOP853',rtol=rtol,atol=atol,t_eval=[target])
        if not run.success:raise ArithmeticError(run.message)
        inside=run.y[:,-1].reshape(n,2,2)
    resolution=MassivePGModeResolution(p);F0=resolution.frame(0.)[0];F3=resolution.frame(3.)[0]
    sewing=np.zeros((n,2,3),complex);sewing[:,0,1]=1.;sewing[:,1,0]=R;sewing[:,1,2]=np.sqrt(T)
    mode0=F0@inside@sewing
    mode3=np.zeros((n,2,3),complex);mode3[:,:,0]=up3@F3.T;mode3[:,:,2]=incoming3@F3.T
    mode3*=np.exp(1j*E*resolution.clock(np.array(3.)))[:,None,None]
    C=np.array([source_covariance(e,p.surface_gravity,p.inheritance_ratio,mass) for e in E])
    P=np.zeros((n,3,3));P[:,0,0]=1;P[:,1,1]=1;P[:,2,2]=open_channel
    return {'mode_at_zero':mode0,'mode_at_three':mode3,'source':C,'projector':P,'reflection':R,'transmission':T,
            'current_residual':current_error,'inner_unitarity_residual':float(np.max(np.linalg.norm(inside@inside.swapaxes(-1,-2).conj()-np.eye(2),axis=(1,2))))}
