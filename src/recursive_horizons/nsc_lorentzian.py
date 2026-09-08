"""Causal radial Dirac propagation on the fixed black-universe PG background.

The characteristic basis diagonalizes sigma_2. This is the fermionic probe on
the given four-dimensional geometry, not a stationary solution of S_one.
Incoming data are imposed only where a characteristic enters the finite box.
The trapped left boundary is outflow for both fields, not a reflecting wall.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import block_diag, bmat, diags, lil_matrix


def geometry(rho):
    rho=np.asarray(rho,dtype=float)
    angle=np.pi/2-np.arctan(rho)
    beta_squared=3*((1+rho*rho)*angle-rho)
    beta=np.sqrt(beta_squared)
    derivative=-3*(1-rho*angle)/beta
    return beta,derivative,1-beta_squared


def horizon():
    return float(brentq(lambda x: float(geometry(x)[2]),0,8,xtol=1e-13))


def sbp_operator(points=401,left=-8.,right=8.,kappa=1.):
    if not isinstance(points,int) or points<9 or points%2!=1 or not left<right:
        raise ValueError('use an odd grid with at least nine ordered points')
    x=np.linspace(left,right,points);h=float(x[1]-x[0])
    weights=np.full(points,h);weights[[0,-1]]=h/2
    derivative=lil_matrix((points,points),dtype=float)
    for i in range(1,points-1):derivative[i,i-1]=-1/(2*h);derivative[i,i+1]=1/(2*h)
    derivative[0,0]=-1/h;derivative[0,1]=1/h
    derivative[-1,-2]=-1/h;derivative[-1,-1]=1/h
    derivative=derivative.tocsr()
    beta,beta_prime,A=geometry(x)
    speeds=np.array([1-beta,-1-beta])
    transports=[];sat=[];outgoing=[];incoming=[]
    for channel in range(2):
        coeff=diags(speeds[channel])
        transports.append(-.5*(coeff@derivative+derivative@coeff))
        penalty=np.zeros(points)
        for node,normal in ((0,-1),(-1,1)):
            normal_speed=float(normal*speeds[channel,node])
            entry=(channel,node,abs(normal_speed))
            if normal_speed<0:
                penalty[node]=normal_speed/weights[node];incoming.append(entry)
            elif normal_speed>0:outgoing.append(entry)
        sat.append(diags(penalty))
    w=kappa/np.sqrt(1+x*x)
    # U† sigma_1 U = sigma_2, hence -i w sigma_2 is a real rotation.
    mixing=bmat([[None,diags(-w)],[diags(w),None]],format='csr')
    volume=block_diag(transports,format='csr')+mixing
    evolution=volume+block_diag(sat,format='csr')
    norm_weights=np.tile(weights,2)
    norm=diags(norm_weights)
    expected=np.zeros(2*points)
    for channel,node,magnitude in outgoing+incoming:
        index=channel*points+(node%points);expected[index]=-magnitude
    residual=norm@evolution+evolution.T@norm-diags(expected)
    error=float(np.max(np.abs(residual.data))) if residual.nnz else 0.
    return {'x':x,'h':h,'weights':weights,'speeds':speeds,'L':evolution,
            'L_volume':volume,'incoming':incoming,'outgoing':outgoing,
            'boundary_form_residual':error,'beta_prime':beta_prime,'A':A}


def packet(x,center=.4,half_width=.25,momentum=4.):
    z=(np.asarray(x)-center)/half_width
    bump=np.zeros_like(z,dtype=float);inside=np.abs(z)<1
    bump[inside]=np.exp(-1/(1-z[inside]**2))
    return bump*np.exp(1j*momentum*(np.asarray(x)-center))


def characteristic(endpoint,time,channel=0):
    def field(t,state):
        beta,_,_=geometry(state[0]);return [(-1 if channel else 1)-float(beta)]
    result=solve_ivp(field,(0,time),[endpoint],method='DOP853',rtol=1e-11,atol=1e-12)
    if not result.success:raise RuntimeError(result.message)
    return float(result.y[0,-1])


def transport_control(x,time,channel=0,center=.4,half_width=.25,momentum=4.):
    """Independent exact-characteristic half-density solution when kappa=0."""
    x=np.asarray(x);n=len(x)
    def backward(t,state):
        positions=state[:n];beta,db,_=geometry(positions)
        speed=(-1 if channel else 1)-beta
        # Backwards flow Y_s=-a(Y), d log(dY/dx)/ds=-a'(Y)=beta'.
        return np.concatenate([-speed,db])
    result=solve_ivp(backward,(0,time),np.concatenate([x,np.zeros(n)]),method='DOP853',
                     rtol=2e-11,atol=2e-12)
    if not result.success:raise RuntimeError(result.message)
    return np.exp(result.y[n:,-1]/2)*packet(result.y[:n,-1],center,half_width,momentum)


def evolve(points=401,final_time=1.2,kappa=1.,center=.4,half_width=.25,momentum=4.,cfl=.45):
    op=sbp_operator(points=points,kappa=kappa)
    x=op['x'];weights=op['weights'];speeds=op['speeds'];L=op['L']
    state=np.zeros((2,points),dtype=complex)
    state[0]=packet(x,center,half_width,momentum)
    normalization=float(np.sqrt(np.sum(weights*np.abs(state[0])**2)))
    if normalization==0:raise ValueError('packet is unresolved or outside box')
    state/=normalization;state=state.ravel()
    initial=state.copy()
    count=int(np.ceil(final_time/(cfl*op['h']/np.max(np.abs(speeds)))))
    dt=final_time/count
    middle=points//2
    if abs(x[middle])>1e-12:raise ValueError('transmission diagnostic requires the symmetric box')
    def rates(y):
        field=y.reshape(2,points)
        physical=sum(magnitude*abs(field[ch,node])**2 for ch,node,magnitude in op['outgoing'])
        penalty=sum(magnitude*abs(field[ch,node])**2 for ch,node,magnitude in op['incoming'])
        left=sum(magnitude*abs(field[ch,node])**2 for ch,node,magnitude in op['outgoing'] if node==0)
        flux=0.
        for ch in range(2):
            a=speeds[ch];u=field[ch]
            flux+=((a[middle-1]+a[middle])*np.real(np.conjugate(u[middle-1])*u[middle])
                  +(a[middle]+a[middle+1])*np.real(np.conjugate(u[middle])*u[middle+1]))/4
        return np.array([physical,penalty,left,-flux])
    integrated=np.zeros(4)
    for _ in range(count):
        k1=L@state;one=state+dt*k1/2
        k2=L@one;two=state+dt*k2/2
        k3=L@two;three=state+dt*k3
        k4=L@three
        integrated+=dt*(rates(state)+2*rates(one)+2*rates(two)+rates(three))/6
        state=state+dt*(k1+2*k2+2*k3+k4)/6
    field=state.reshape(2,points)
    density=np.sum(np.abs(field)**2,axis=0)
    child_weights=weights.copy();child_weights[x>0]=0;child_weights[middle]/=2
    total=float(weights@density);child=float(child_weights@density)
    initial_child=float(child_weights@np.sum(np.abs(initial.reshape(2,points))**2,axis=0))
    lower=characteristic(center-half_width,final_time,1)
    upper=characteristic(center+half_width,final_time,0)
    outside=(x<lower-4*op['h'])|(x>upper+4*op['h'])
    return {'points':points,'time':final_time,'kappa':kappa,'x':x,'field':field,
            'initial':initial.reshape(2,points),'initial_normalization':normalization,'h':op['h'],
            'steps':count,'dt':dt,'probability':total,'child_probability':child,
            'physical_outflow':float(integrated[0]),'inflow_penalty_debit':float(integrated[1]),
            'left_child_outflow':float(integrated[2]),'throat_transmission':float(integrated[3]),
            'probability_balance_residual':float(total+integrated[0]+integrated[1]-1),
            'child_balance_residual':float(child+integrated[2]-integrated[3]-initial_child),
            'causal_interval':[lower,upper],'outside_characteristic_cone_probability':float(np.sum(weights[outside]*density[outside])),
            'outside_horizon_probability':float(np.sum(weights[x>horizon()]*density[x>horizon()])),
            'boundary_form_residual':op['boundary_form_residual'],
            'incoming_channels':[(ch,int(node)) for ch,node,_ in op['incoming']],
            'outgoing_channels':[(ch,int(node)) for ch,node,_ in op['outgoing']]}
