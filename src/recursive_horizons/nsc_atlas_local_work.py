"""Local fourth-order scalar constraint from a specified finite flow.

Double padding resolves the cubic source without aliasing for the retained
strict 2/3 fluid band. Its nonzero mean needs affine expansion or boundary
flux; only the zero-mean part has a periodic Poisson inverse.
"""
import numpy as np
from .nsc_atlas_flow_evolution import PeriodicFlow


def padded_spectrum(h,size):
    n=h.shape[-3]
    if size<2*n or size%2:raise ValueError('even padding at least twice the fluid grid required')
    out=np.zeros((*h.shape[:-3],size,size,size//2+1),dtype=complex)
    h0=n//2
    for sx,tx in ((slice(0,h0),slice(0,h0)),(slice(h0,n),slice(size-h0,size))):
        for sy,ty in ((slice(0,h0),slice(0,h0)),(slice(h0,n),slice(size-h0,size))):
            out[...,tx,ty,:h0+1]=h[...,sx,sy,:]*(size/n)**3
    return out


def constraint_fields(flow,h):
    big=PeriodicFlow(2*flow.n,length=flow.length,viscosity=flow.nu,workers=flow.workers)
    H=padded_spectrum(h,big.n);u=big.real(H);curl=big.real(big.curl(H))
    loss=flow.nu*flow.length**3*float(np.mean(np.sum(curl*curl,axis=0)))
    cross=np.array([u[1]*curl[2]-u[2]*curl[1],u[2]*curl[0]-u[0]*curl[2],u[0]*curl[1]-u[1]*curl[0]])
    N=big.fft(cross);del cross,curl
    p_hat=-1j*np.sum(big.k*N,axis=0)*big.inv-big.fft(np.sum(u*u,axis=0)/2);del N
    p_hat[0,0,0]=0
    lap=big.real(-big.k2*H)
    source=-np.sum(u*lap,axis=0);del lap
    pressure_term=np.zeros_like(source)
    for i in range(3):pressure_term+=u[i]*big.real(1j*big.k[i]*p_hat)/flow.nu
    source+=pressure_term
    strain_term=np.zeros_like(source)
    for i in range(3):
        strain_term+=2*big.real(1j*big.k[i]*H[i])**2
        for j in range(i+1,3):strain_term+=big.real(1j*(big.k[i]*H[j]+big.k[j]*H[i]))**2
    source+=strain_term;del strain_term,u,H
    mean=float(np.mean(source));hubble=mean/3
    centered=source-mean
    potential_hat=-big.fft(centered)*big.inv;potential_hat[0,0,0]=0
    divergence=big.real(-big.k2*potential_hat)
    residual=mean+divergence-source
    correction_norm=np.zeros_like(source)
    for i in range(3):correction_norm+=big.real(1j*big.k[i]*potential_hat)**2
    return dict(source=source,periodic_divergence=divergence,potential=big.real(potential_hat),
        mean=mean,affine_rate=hubble,loss=loss,pressure_mean=float(np.mean(pressure_term)),
        residual_rms=float(np.sqrt(np.mean(residual*residual))),
        residual_max=float(np.max(np.abs(residual))),
        periodic_correction_rms=float(np.sqrt(np.mean(correction_norm))),
        periodic_correction_max=float(np.sqrt(np.max(correction_norm))),
        padded_size=big.n)


def summary(flow,h):
    fields=constraint_fields(flow,h);source=fields.pop('source')
    fields.pop('periodic_divergence');fields.pop('potential')
    target=2*fields['loss']/(flow.length**3*flow.nu)
    fields.update(expected_mean=target,mean_identity_error=fields['mean']-target,
        total_outward_coefficient=flow.length**3*fields['mean'],
        source_rms=float(np.sqrt(np.mean(source*source))),
        source_min=float(np.min(source)),source_max=float(np.max(source)),
        local_residual_relative_rms=fields['residual_rms']/max(float(np.sqrt(np.mean(source*source))),1e-30),
        periodic_only_mean_residual=-fields['mean'],
        source_slice_z0=source[::max(1,fields['padded_size']//24),::max(1,fields['padded_size']//24),fields['padded_size']//2].tolist())
    return fields
