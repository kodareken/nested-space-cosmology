"""Finite-operation replay of authenticated coherent matter contractions."""
import numpy as np
from .nsc_common_subtracted_ks_source import raw_ks_vertex_coefficients
from .nsc_weighted_receipt_replay import gamma, require_record


def matter_bounds(arrays, name, parameters):
    """Uncancelled bounds for F C F* and four tangent products.

    Eight real operations per complex product/sum and two matrix products
    give 2*gamma_(8n+16)+gamma^2. The extra16 charges column-weight products
    and momentum subtraction/division. gamma48 charges four-entry traces,
    addition, multiplicity and endpoint vertices (one division/fixed Pauli
    entries for N/beta). Both executions are charged by the comparator.
    This is a finite arithmetic guard, not a field or physical error bound.
    """
    weights=np.asarray(arrays['source/column_weights'])
    C=abs(np.asarray(arrays['source/covariance']))
    F,Fz,dF,dFz=(abs(np.asarray(arrays[name+'/'+key])*weights) for key in ('F','Fz','dF','dFz'))
    adj=lambda x:np.swapaxes(x,-1,-2)
    rho=F@C@adj(F)
    momentum=.5*(Fz@C@adj(F)+F@C@adj(Fz))
    drho=dF@C@adj(F)+F@C@adj(dF)
    dmomentum=.5*(dFz@C@adj(F)+Fz@C@adj(dF)+dF@C@adj(Fz)+F@C@adj(dFz))
    vertices=raw_ks_vertex_coefficients(
        np.array([1.,0.,parameters['axial_scale'],parameters['radius']]),
        parameters['mass'],parameters['angular'],envelopes=np.ones(4))
    M,V=abs(vertices.multiplication[:2]),abs(vertices.momentum[:2])
    trace=lambda D,J:np.einsum('bij,...ji->...b',M,D)+np.einsum('bij,...ji->...b',V,J)
    g=gamma(8*len(C)+16)
    coefficient=abs(parameters['multiplicity'])*(2*g+g*g+gamma(48))
    return coefficient*trace(rho,momentum),coefficient*trace(drho,dmomentum)


def require_envelope_record(saved,current,arrays,meta):
    parameters=meta['parameters'];bounds={};values={};errors={}
    comparison=saved['schema']=='NSC-KS-ENVELOPE-COMPARISON-v1'
    difference=saved['schema']=='NSC-KS-DIFFERENCE-CONTROL-v1'
    names=meta['completed'] if comparison else ['previous512' if difference else 'previous256',*meta['completed']]
    from .nsc_ks_source_envelope import comparison_control as C
    for name in set(names)|({'analytic_reference','pg_cf4_240'} if comparison else set()):
        values[name]=C.matter(arrays,name,parameters)
        errors[name]=matter_bounds(arrays,name,parameters)
    def gap(a,b,tangent=False):
        field='action_gradient_tangent' if tangent else 'action_gradient'
        index=1 if tangent else 0
        error=errors[a][index]+errors[b][index]+np.finfo(float).eps*(abs(values[a][field])+abs(values[b][field]))
        return np.max(error,axis=(0,1) if tangent else 0)
    if comparison:
        for name in names:
            bounds['measurements','per_case',name,'matter']=errors[name][0]
            if errors[name][1].size:
                bounds['measurements','per_case',name,'matter_tangent']=errors[name][1]
        checks=current['measurements']['checks']
        if 'reference_matter_drift' in checks:
            bounds['measurements','checks','reference_matter_drift']=gap('reference64','analytic_reference')
        for a,b in (('n64','n128'),('n128','n256')):
            if a+'_'+b in checks:
                bounds['measurements','checks',a+'_'+b,'matter_change']=gap(a,b)
        if 'tight' in names:
            for field,other in (('time_tolerance_change',meta['tight_base']),('PG_CF4_matter_difference','pg_cf4_240'),('local_matter_correction_max','analytic_reference')):
                bounds['measurements','checks',field]=gap('tight',other)
    elif saved['schema'] in ('NSC-KS-ENVELOPE-RESOLUTION-v1','NSC-KS-DIFFERENCE-CONTROL-v1'):
        pairs=(('n256','n512','space'),('n512','tighter512','time'),('previous512','tighter512','representation')) if difference else (('previous256','n512','space'),('n512','tighter512','time'))
        for a,b,label in pairs:
            if label in current['measurements']:
                bounds['measurements',label,'matter_change']=gap(a,b)
                bounds['measurements',label,'tangent_change']=gap(a,b,True)
    else:
        raise ValueError('explicit envelope receipt recipe required')
    return require_record(saved,current,bounds)
