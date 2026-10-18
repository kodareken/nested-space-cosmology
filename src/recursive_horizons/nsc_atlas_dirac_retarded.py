"""Retarded endpoint maps extracted from the existing Dirac field solver.

Incoming cardinal sample columns are numerical probes of a linear map,
not occupied quantum modes or a spectral cutoff. Each field grid is
discarded after its endpoint row and analytic spatial row are read.
The complete root covariance is paired separately in null coordinates.
"""
from concurrent.futures import ThreadPoolExecutor
import math

import numpy as np
from numpy.polynomial.legendre import leggauss,legvander

from .nsc_atlas_dirac_mean import MeanModeGrid
from .nsc_atlas_dirac_scattering import _inputs,_real_array
from .nsc_atlas_dirac_root_pairing import NullResponse


def volterra_endpoint_jets(W,Wu,Wv,square_integrals,mixed_integrals):
    """One-sided value, observation-x and source-s endpoint jets.

    All come from the same Volterra equation, including its moving lower
    limit. The source derivative is distinct from the observation derivative.
    Rows are incoming legs; columns are outgoing spinor components.
    """
    I=np.asarray(square_integrals,dtype=float);J=np.asarray(mixed_integrals,dtype=float)
    if (I.shape!=(2,) or J.shape!=(2,) or not np.isfinite(I).all() or np.any(I<0)
            or not np.isfinite(J).all() or not np.isfinite([W,Wu,Wv]).all()):
        raise ValueError('finite potential jets and two declared boundary integrals required')
    value=np.array([[-I[0]/4,-.5j*W],[-.5j*W,-I[1]/4]])
    spatial=np.array([[-W*W/4+J[0]/4-I[0]**2/32,-.5j*Wv-.125j*W*I[0]],
        [.5j*Wu+.125j*W*I[1],W*W/4-J[1]/4+I[1]**2/32]])
    source=np.array([[-J[0]/4-I[0]**2/32,-.5j*Wu-.125j*W*I[0]],
        [-.5j*Wv-.125j*W*I[1],-J[1]/4-I[1]**2/32]])
    return dict(value=value,spatial=spatial,source=source)


def _project_cells(row,length,width,order,endpoint):
    """Project a continuous linear reconstruction of the one-sided cell map.

    Regular source cells have width h; the final cell ends at the last
    cell CENTER and has width h/2. A coefficient row is its integrated
    ordinary kernel, after direct evaluation has been removed. Densities
    live at the source sample centers; their exact causal endpoint is
    supplied. A step reconstruction only gives first-order L2 accuracy.
    """
    cells=row.shape[1]
    edges=np.r_[np.arange(cells)*width,length]
    cell_widths=np.diff(edges)
    density=(row/cell_widths).T.copy()
    density[-1]=endpoint
    knots=np.r_[0.,(np.arange(cells)+.5)*width]
    values=np.vstack((1.5*density[0]-.5*density[1],density))
    nodes,weights=leggauss((order+2)//2)
    fraction=(nodes+1)/2
    points=knots[:-1,None]+np.diff(knots)[:,None]*fraction
    g=values[:-1,None,:]*(1-fraction[None,:,None])+values[1:,None,:]*fraction[None,:,None]
    qweights=np.diff(knots)[:,None]*weights/2
    coefficients=np.einsum('sq,sqn,sqi->ni',qweights,legvander(2*points/length-1,order-1),g)
    coefficients*=(2*np.arange(order)+1)[:,None]/length
    cell_norm=float(np.sum(qweights*np.sum(abs(g)**2,axis=-1)))
    polynomial_norm=float(np.sum(length/(2*np.arange(order)+1)*np.sum(abs(coefficients)**2,axis=1)))
    return coefficients,math.sqrt(max(0.,cell_norm-polynomial_norm))


def retarded_null_responses(potential,du,dv,*,potential_u,potential_v,
                           boundary_square_integrals,boundary_mixed_integrals,
                           order=12,workers=8):
    """Field and analytic x response at the final cell center.

    boundary_square_integrals supplies the independently evaluated pair
    (integral_0^V W(U,t)^2 dt, integral_0^U W(t,V)^2 dt). These fix the
    actual moving-endpoint terms, not a fitted or differentiated waveform.
    boundary_mixed_integrals supplies (integral W*W_u along V,
    integral W*W_v along U), fixing the spatial kernel's endpoint traces.
    The returned projection residual concerns the cell map only; it is
    NOT a continuum field, vacuum stress or angular-tail error bound.
    """
    w,du,dv,_=_inputs(potential,du,dv,None,2)
    wu,wv=(_real_array(v,'potential derivative',2) for v in (potential_u,potential_v))
    nu,nv=w.shape
    integrals=np.asarray(boundary_square_integrals,dtype=float)
    mixed=np.asarray(boundary_mixed_integrals,dtype=float)
    if (min(nu,nv)<4 or wu.shape!=w.shape or wv.shape!=w.shape
            or isinstance(order,bool) or not isinstance(order,int) or not 2<=order<=min(nu,nv)
            or isinstance(workers,bool) or not isinstance(workers,int) or not 1<=workers<=8
            or integrals.shape!=(2,) or not np.isfinite(integrals).all() or np.any(integrals<0)
            or mixed.shape!=(2,) or not np.isfinite(mixed).all()):
        raise ValueError('resolved field jets, two nonnegative boundary integrals, polynomial order and one to eight workers required')
    # Four complex spinor jets per worker, plus its outgoing boundary data.
    live_estimate=workers*(4*nu*nv*2*16+(nu+nv)*16)
    if live_estimate>256*1024**2:
        raise ValueError('streamed endpoint extraction exceeds the admitted 256 MiB field allowance')
    zeros_u=np.zeros(nu,complex);zeros_v=np.zeros(nv,complex)

    def column(index):
        plus=zeros_u.copy();minus=zeros_v.copy()
        if index<nu:
            plus[index]=1.
        else:
            minus[index-nu]=1.
        grid=MeanModeGrid(w,du,dv,potential_u=wu,potential_v=wv,
            incoming_plus=plus,incoming_minus=minus,plus_u=zeros_u,minus_v=zeros_v)
        return np.stack((grid.jets[0][-1,-1],grid.jets[2][-1,-1]-grid.jets[1][-1,-1]))

    first=column(0)  # Compile before starting independent GIL-free columns.
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows=np.stack([first,*pool.map(column,range(1,nu+nv))],axis=-1)
    direct=np.eye(2,dtype=complex)
    W=float(w[-1,-1])
    direct_x=np.array([[integrals[0]/4,-.5j*W],[.5j*W,-integrals[1]/4]])
    rows[0,:,nu-1]-=direct[:,0];rows[0,:,-1]-=direct[:,1]
    rows[1,:,nu-1]-=direct_x[:,0];rows[1,:,-1]-=direct_x[:,1]
    lengths=((nu-.5)*du,(nv-.5)*dv)
    jets=volterra_endpoint_jets(W,wu[-1,-1],wv[-1,-1],integrals,mixed)
    endpoint_values=np.stack((jets['value'],jets['spatial']))
    responses=[];residuals=[]
    for which in range(2):
        green=[];projection=[]
        for leg,(length,width,part) in enumerate(zip(lengths,(du,dv),(slice(0,nu),slice(nu,nu+nv)))):
            coefficients,residual=_project_cells(rows[which,:,part],length,width,order,endpoint_values[which,leg])
            green.append(coefficients);projection.append(residual)
        responses.append(NullResponse(lengths,tuple(green),direct if which==0 else direct_x,
            np.zeros((2,2)) if which==0 else np.diag([-1.,1.])))
        residuals.append(projection)
    return dict(cell_field=responses[0],cell_spatial=responses[1],endpoint_values=endpoint_values,
        source_endpoint_values=jets['source'],
        cell_projection_residual=np.asarray(residuals),
        live_field_bytes_estimate=live_estimate,workers=workers,
        scope='numerical_retarded_maps_from_existing_characteristic_solver')


def retarded_response(sample,kind,order):
    """Restrict a cell map while enforcing its known endpoint and mean.

    Among polynomial corrections of zero integral, this is the smallest
    L2 correction giving the exact Volterra endpoint trace. It enforces a
    field identity, not a target covariance. Reapply at each retained order.
    """
    if kind not in ('field','spatial'):
        raise ValueError('field or spatial retarded map required')
    response=sample['cell_'+kind]
    if isinstance(order,bool) or not isinstance(order,int) or not 2<=order<=min(map(len,response.green)):
        raise ValueError('admitted polynomial order required')
    endpoint=sample['endpoint_values'][0 if kind=='field' else 1]
    green=[]
    for leg in (0,1):
        coefficients=response.green[leg][:order].copy()
        correction=endpoint[leg]-np.sum(coefficients,axis=0)
        coefficients[1:]+=(2*np.arange(1,order)+1)[:,None]*correction/(order*order-1)
        green.append(coefficients)
    return NullResponse(response.lengths,tuple(green),response.direct,response.direct_derivative)
