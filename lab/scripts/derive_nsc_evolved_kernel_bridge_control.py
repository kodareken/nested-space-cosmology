#!/usr/bin/env python3
"""Saved-state same-kernel E/k/eta control; no complete-spectrum promotion."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_common_subtracted_ks_source import (
    SpectrumDeclaration,physical_source_kernel,reference_symbol_kernel,
    raw_ks_vertex_coefficients,symmetric_pairing_density)
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_spacetime_geometry_jets import spacetime_geometry_jets
from recursive_horizons.nsc_scaled_reference_projector import symbol_metric_from_geometry
from recursive_horizons.nsc_spatial_reference_symbol import reference_projector
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction

OUTPUT=ROOT/'results/development/nsc-evolved-kernel-bridge-control.json'
FINE=ROOT/'results/development/nsc-ks-fine-trajectory.json'
OWNED=('scripts/derive_nsc_evolved_kernel_bridge_control.py','docs/nsc-evolved-kernel-bridge-control.md',
       'src/recursive_horizons/nsc_common_subtracted_ks_source.py','src/recursive_horizons/nsc_evolved_incoming_constraints.py',
       'src/recursive_horizons/nsc_scaled_reference_projector.py','src/recursive_horizons/nsc_ks_spacetime_geometry_jets.py',
       'src/recursive_horizons/nsc_spatial_reference_symbol.py','src/recursive_horizons/nsc_local_incoming_family.py')


def digest(path):
    p=Path(path);return sha256((p if p.is_absolute() else ROOT/p).read_bytes()).hexdigest()


def compute():
    record=json.loads(FINE.read_text());payload=record['payload']
    if digest(payload['path'])!=payload['sha256']:raise ValueError('saved incoming fields changed')
    with np.load(ROOT/payload['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
    meta=json.loads(a['metadata_json'].tobytes());parameters=meta['parameters']
    E=a['source/energies'].reshape(-1,3)
    weights=a['source/column_weights'].reshape(-1,3)
    if not np.all(E==E[:,:1]) or not np.all(weights==weights[:,:1]):raise ValueError('complete coherent three-source fibers required')
    E=E[:,0];weights=2*np.pi*weights[:,0]**2
    n=len(E);cov=a['source/covariance']
    blocks=np.array([cov[3*i:3*i+3,3*i:3*i+3] for i in range(n)])
    ix=np.arange(3*n)//3
    if np.any(cov[ix[:,None]!=ix[None,:]]!=0):raise ValueError('cross-energy source coherence requires a different kernel owner')
    z=a['target_z'];middle=len(z)//2;center=float(z[middle])
    envelope=LocalAxialFunction((1.,),center,inner=.02,outer=.04)
    if not meta['period_origin']<center-.04<center+.04<meta['period_origin']+meta['period_length']:
        raise ValueError('compact test variation must vanish before the cell ends')
    F=a['prepared/columns'].reshape(len(z),2,n,3).transpose(0,2,1,3)
    Fz=a['prepared/axial_columns'].reshape(len(z),2,n,3).transpose(0,2,1,3)
    Aref=a['prepared/reference_amplitudes']
    phase=np.exp(-1j*a['source/energies']*z[:,None])[:,None,:]
    ref=(phase*Aref[None]).reshape(len(z),2,n,3).transpose(0,2,1,3)
    refz=-1j*E[None,:,None,None]*ref
    # Independent k quadrature: neither its nodes nor weights equal the source table.
    x,w=leggauss(24);k=4*x;kw=4*w
    families=(P.family(P.ALPHA),P.family(0.))
    projectors=[]
    for family in families:
        geometry=spacetime_geometry_jets(family,1.,center,order=4)
        fields=symbol_metric_from_geometry(geometry)
        result=reference_projector(*fields,k,np.full(len(k),parameters['mass']),np.full(len(k),parameters['angular']))
        projectors.append(result['orders'])
    source_declaration=SpectrumDeclaration('control_only','four original signed source energies from the saved evolved history')
    reference_declaration=SpectrumDeclaration('control_only','independent24-node canonical k quadrature on[-4,4]')
    metric=np.array([1.,0.,parameters['axial_scale'],parameters['radius']])
    vertex=raw_ks_vertex_coefficients(metric,parameters['mass'],parameters['angular'],envelopes=np.ones(4))
    mu=parameters['multiplicity'];rows=[]
    for offset in (0,1,2):
        plus,minus=middle+offset,middle-offset;eta=float(z[plus]-z[minus])
        plus_vertex=raw_ks_vertex_coefficients(metric,parameters['mass'],parameters['angular'],
            envelopes=np.ones(4)*envelope(float(z[plus]),0))
        minus_vertex=raw_ks_vertex_coefficients(metric,parameters['mass'],parameters['angular'],
            envelopes=np.ones(4)*envelope(float(z[minus]),0))
        values=[]
        sources=[];refs=[]
        for field,derivative,projector in ((F,Fz,projectors[0]),(ref,refz,projectors[1])):
            source=physical_source_kernel(field[plus],field[minus],derivative[plus],derivative[minus],blocks,E,weights,
                separation=eta,declaration=source_declaration)
            reference=reference_symbol_kernel(projector,k,kw,separation=eta,declaration=reference_declaration)
            pair=symmetric_pairing_density(source,reference,plus_vertex,minus_vertex,compact_spatial_boundary=True)
            values.append(mu*pair.action_gradient_density[:2]);sources.append(source);refs.append(reference)
        rows.append({'eta':eta,'paired_delta_N_beta':(values[0]-values[1]).tolist(),
            'physical_delta_kernel_max':float(np.max(abs(sources[0].kernel-sources[1].kernel))),
            'reference_delta_kernel_max':float(np.max(abs(refs[0].kernel-refs[1].kernel))),
            'source_scope':'control_only','reference_scope':'control_only'})
    p=parameters
    def raw(field,derivative):
        f=field[middle].transpose(1,0,2).reshape(1,2,-1)*a['source/column_weights']
        fz=derivative[middle].transpose(1,0,2).reshape(1,2,-1)*a['source/column_weights']
        zero=np.zeros((0,*f.shape),complex)
        return source_column_matter(f,fz,cov,zero,zero,**p)['action_gradient'][0]
    raw_delta=raw(F,Fz)-raw(ref,refz)
    def raw_reference(Pj):
        H=vertex.multiplication[:2,None]+vertex.momentum[:2,None]*k[None,:,None,None]
        return -mu*np.einsum('bkij,kji,k->b',H,Pj.sum(axis=0),kw/(2*np.pi)).real
    reference_delta=raw_reference(projectors[0])-raw_reference(projectors[1])
    identity=float(np.max(abs(np.array(rows[0]['paired_delta_N_beta'])-(raw_delta-reference_delta))))
    return {'schema':'NSC-EVOLVED-KERNEL-BRIDGE-CONTROL-v1','accountable_author':'Douglas Ek',
        'status':('PASS' if identity<=3e-11 else 'OPEN')+': finite same-kernel identity; full ordered limits OPEN',
        'identity_residual':identity,'identity_tolerance':3e-11,'rows':rows,
        'raw_source_delta_N_beta':raw_delta.tolist(),'partial_reference_delta_N_beta':reference_delta.tolist(),
        'finite_eta_minus_zero':[(np.array(r['paired_delta_N_beta'])-np.array(rows[0]['paired_delta_N_beta'])).tolist() for r in rows[1:]],
        'reference_k_interval':[-4.,4.],'reference_k_count':24,'source_energies':E.tolist(),
        'compact_test_variation':{'center':center,'inner':.02,'outer':.04,'polynomial_coefficients':[1.]},
        'scope':{'physical_local_gate':'OPEN','source_tail_bound':None,'reference_tail_bound':None,
            'coincidence_remainder_bound':None,'finite_band_exchange_bound':None,
            'certified_local_reference_coefficients_replaced':False,'band_bulk_identity_reopened':False,
            'control_reference_evolution':'saved jointly evolved homogeneous reference',
            'compact_test_vertex':'unit plateau around the probe, vanishing before axial cell ends',
            'new_field_or_source_evolutions':0,'metric_timestep':False},
        'source_hashes':{p:digest(p) for p in OWNED},
        'input_hashes':{str(FINE.relative_to(ROOT)):digest(FINE),payload['path']:payload['sha256']},
        'reproducer':'python scripts/derive_nsc_evolved_kernel_bridge_control.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();result=compute()
    if args.record:
        if OUTPUT.exists():raise FileExistsError('kernel bridge control exists; use --check')
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif result!=json.loads(OUTPUT.read_text()):raise ValueError('kernel bridge replay differs')
    print(json.dumps({k:result[k] for k in ('status','identity_residual','finite_eta_minus_zero')},indent=2))
