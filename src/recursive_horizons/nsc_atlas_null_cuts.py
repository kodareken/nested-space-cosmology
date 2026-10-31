"""Parametric closed S2 cuts and complete initial optical data on S3 slices.

A cut is one common coordinate-time surface, sampled in the existing
Gauss-theta/periodic-phi ordering. Its four embedding components are
projected onto real S2 harmonics, differentiated analytically and then
normalized to S3 with the complete first/second derivative chain. No
star-shaped assumption, radial zero or independent per-ray event is used.

The physical regular frame obeys E_i=r*A_i*X and [E_i,E_j]=2*r*eps_ijk E_k.
K=-Lie_n(h)/2. The second form b_AB=h(D_eA s,eB) gives optical matrices
B_plus=b-K and B_minus=-b-K in an orthonormal screen. Their traces are
theta_plus and theta_minus for k=n+/-s. Initial null/screen components
are in (partial_t,E_i), matching the parent geodesic transport interface.

Immersion, coincidence and orientation checks below are sampled controls,
not a proof of global injectivity or absence of caustics between samples.
No result from this module admits a cut as a child universe.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.spatial import cKDTree

from .nsc_atlas_trapped_surfaces import (
    _angular_basis, _connection, _degree, _embedding, _geometry, _order, _unit_frame,
)


def _points(value,count):
    X=np.asarray(value)
    if X.dtype.kind not in 'iuf' or X.shape!=(count,4) or not np.all(np.isfinite(X)):
        raise ValueError('cut points require finite real shape (n_theta*n_phi,4)')
    X=np.asarray(X,dtype=float)
    if np.max(abs(np.sum(X*X,axis=1)-1))>2.e-8:
        raise ValueError('input cut points must lie on unit S3')
    return X


def _normalize_jets(raw,first,second):
    norm=np.linalg.norm(raw,axis=1)
    if np.min(norm)<=1.e-8:
        raise ValueError('fitted embedding approaches zero; normalization is singular')
    X=raw/norm[:,None]
    dn=np.einsum('ni,nai->na',X,first)
    ddn=(np.einsum('nai,nbi->nab',first,first)
         +np.einsum('ni,nabi->nab',raw,second)-dn[:,:,None]*dn[:,None,:])/norm[:,None,None]
    dX=(first-X[:,None,:]*dn[:,:,None])/norm[:,None,None]
    ddX=(second-dX[:,:,None,:]*dn[:,None,:,None]-dX[:,None,:,:]*dn[:,:,None,None]
          -X[:,None,None,:]*ddn[:,:,:,None])/norm[:,None,None,None]
    return X,dX,ddX,norm


def _provider_data(provider,X):
    raw=provider(X)
    h,dh,K=_geometry(lambda _:raw,X)
    data=dict(raw) if isinstance(raw,dict) else {}
    data.update(h=h,dh=dh,K=K)
    return data


def _lapse_shift(data,count):
    try:
        lapse=np.broadcast_to(np.asarray(data.get('lapse',1.),dtype=float),(count,))
        shift=np.broadcast_to(np.asarray(data.get('shift',0.),dtype=float),(count,3))
    except ValueError as exc:raise ValueError('lapse/shift must broadcast to (N,) and (N,3)') from exc
    if not np.all(np.isfinite(lapse)) or np.min(lapse)<=0 or not np.all(np.isfinite(shift)):
        raise ValueError('finite positive lapse and finite shift required')
    return lapse,shift


def _physical_surface(X,dX,ddX,angular,provider,r,*,orientation=1,outgoing_k=None,
                      immersion_tolerance=1.e-18):
    r=float(r);tolerance=float(immersion_tolerance)
    if not math.isfinite(r) or r<=0 or not math.isfinite(tolerance) or tolerance<0:
        raise ValueError('positive frame scale and nonnegative immersion tolerance required')
    if isinstance(orientation,(bool,np.bool_)) or orientation not in (-1,1):
        raise ValueError('orientation must be +1 or -1 for the continuous parameter sheet')
    refgram=np.einsum('nai,nbi->nab',dX,dX)
    reference_density_squared=np.linalg.det(refgram)/np.sin(angular['theta'])**2
    if np.min(reference_density_squared)<=tolerance:
        raise ValueError('sampled cut loses immersion or reaches its declared immersion threshold')
    # This detects sampled nonlocal coincidences, not all triangle crossings
    # or global injectivity. Near angular coordinate poles are local on S2.
    coincident=cKDTree(X).query_pairs(1.e-10,output_type='ndarray')
    if len(coincident):
        separation=np.linalg.norm(angular['directions'][coincident[:,0]]-angular['directions'][coincident[:,1]],axis=1)
        if np.any(separation>1.e-5):
            raise ValueError('sampled nonlocal self-coincidence; cut cannot be treated as embedded')
    data=_provider_data(provider,X);h,dh,K=(data[key] for key in ('h','dh','K'))
    hi,Gamma=_connection(h,dh,r)
    frame=_unit_frame(X)
    V=np.einsum('nik,nak->nai',frame,dX)/r
    dV=np.empty((len(X),2,2,3))
    for a in range(2):
        derivative_frame=_unit_frame(dX[:,a])
        for b in range(2):
            dV[:,a,b]=(np.einsum('nik,nk->ni',derivative_frame,dX[:,b])
                       +np.einsum('nik,nk->ni',frame,ddX[:,a,b]))/r
    induced=np.einsum('nai,nij,nbj->nab',V,h,V)
    gramdet=np.linalg.det(induced)
    if np.min(gramdet)<=0:raise ValueError('physical surface Gram determinant is not positive')
    covnormal=orientation*np.cross(V[:,0],V[:,1])
    normal=np.einsum('nij,nj->ni',hi,covnormal)
    scale=np.sqrt(np.einsum('ni,ni->n',covnormal,normal))
    normal/=scale[:,None];covnormal/=scale[:,None]
    normality=null_defect=None;orientation_alignment=None
    if outgoing_k is not None:
        k=np.asarray(outgoing_k,dtype=float)
        if k.shape!=(len(X),4) or not np.all(np.isfinite(k)) or np.min(k[:,0])<=0:
            raise ValueError('outgoing_k must be finite future-directed (N,4) frame components')
        lapse,shift=_lapse_shift(data,len(X))
        spatial=k[:,1:]+shift*k[:,0,None]
        dot=np.einsum('ni,nij,nj->n',normal,h,spatial)
        if np.all(dot<0):normal=-normal;covnormal=-covnormal;dot=-dot
        if np.min(dot)<=0:
            raise ValueError('mixed outgoing orientation across the cut; folded/caustic sheet is not admitted')
        energy=lapse*k[:,0]
        tangent=spatial-dot[:,None]*normal
        normality=float(np.max(np.sqrt(np.maximum(0,np.einsum('ni,nij,nj->n',tangent,h,tangent)))/energy))
        null_defect=float(np.max(abs(np.einsum('ni,nij,nj->n',spatial,h,spatial)-energy**2)/energy**2))
        orientation_alignment=float(np.min(dot/energy))
    acceleration=dV+np.einsum('nkij,nai,nbj->nabk',Gamma,V,V)
    b=-np.einsum('nk,nabk->nab',covnormal,acceleration)
    # Orthonormalize in the physical induced metric, without a fitted QR.
    T=np.zeros((len(X),2,2));T[:,0,0]=1/np.sqrt(induced[:,0,0])
    schur=induced[:,1,1]-induced[:,0,1]**2/induced[:,0,0]
    if np.min(schur)<=0:raise ValueError('screen orthonormalization lost immersion')
    T[:,1,1]=1/np.sqrt(schur)
    T[:,1,0]=-induced[:,0,1]/induced[:,0,0]*T[:,1,1]
    screen=np.einsum('nAa,nai->nAi',T,V)
    b_screen=np.einsum('nAa,nab,nBb->nAB',T,b,T)
    K_screen=np.einsum('nAi,nij,nBj->nAB',screen,K,screen)
    optical_plus=b_screen-K_screen;optical_minus=-b_screen-K_screen
    area_weights=angular['weights']*np.sqrt(gramdet)/np.sin(angular['theta'])
    return dict(data=data,normal=normal,screen=screen,b_screen=b_screen,K_screen=K_screen,
                optical_plus=optical_plus,optical_minus=optical_minus,
                theta_plus=np.trace(optical_plus,axis1=1,axis2=2),
                theta_minus=np.trace(optical_minus,axis1=1,axis2=2),area_weights=area_weights,
                immersion_min_gramdet=float(np.min(gramdet)),
                immersion_min_reference_density_squared=float(np.min(reference_density_squared)),
                outgoing_normality_max_relative=normality,outgoing_null_defect_max_relative=null_defect,
                outgoing_alignment_min=orientation_alignment)


@dataclass(frozen=True)
class ParametricCutEvaluation:
    points: np.ndarray
    theta_plus: np.ndarray
    theta_minus: np.ndarray
    normal: np.ndarray
    screen: np.ndarray
    second_form: np.ndarray
    K_screen: np.ndarray
    area_weights: np.ndarray
    angular_degree: int
    order: tuple
    diagnostics: dict
    expansion_tolerance: float = 1.e-7

    @property
    def area(self):return float(np.sum(self.area_weights))

    @property
    def plus_rms(self):return float(np.sqrt(np.dot(self.area_weights,self.theta_plus**2)/self.area))

    @property
    def minus_rms(self):return float(np.sqrt(np.dot(self.area_weights,self.theta_minus**2)/self.area))

    @property
    def future_trapped_sampled(self):
        return bool(np.max(self.theta_plus)<-self.expansion_tolerance and np.max(self.theta_minus)<-self.expansion_tolerance)

    @property
    def future_mots_sampled(self):
        return bool(np.max(abs(self.theta_plus))<=self.expansion_tolerance and np.max(self.theta_minus)<-self.expansion_tolerance)

    def summary(self):
        return dict(self.diagnostics,area=self.area,theta_plus_rms=self.plus_rms,theta_minus_rms=self.minus_rms,
            theta_plus_range=[float(np.min(self.theta_plus)),float(np.max(self.theta_plus))],
            theta_minus_range=[float(np.min(self.theta_minus)),float(np.max(self.theta_minus))],
            future_trapped_sampled=self.future_trapped_sampled,future_mots_sampled=self.future_mots_sampled,
            angular_degree=self.angular_degree,order=list(self.order),expansion_tolerance=self.expansion_tolerance,
            global_injectivity_certified=False,continuum_expansion_signs_certified=False,child_formation_proved=False,
            scope='One common coordinate-time cut, finite harmonic reconstruction and angular samples. No global injectivity, between-sample caustic exclusion or child admission.')


def evaluate_parametric_cut(points,provider,*,frame_scale,angular_degree,order,orientation=1,
                            outgoing_k=None,caustic_encountered=False,ray_area_determinant=None,
                            immersion_tolerance=1.e-18,expansion_tolerance=1.e-7):
    """Reconstruct and measure one transported cut without a radial-graph model.

    ``outgoing_k`` can orient the continuous sheet using the transported
    future null generators. Its tangential/normality and null defects are
    reported, not hidden. A known upstream caustic or nonpositive optical
    determinant rejects the cut. Their omission supplies no history proof.
    """
    degree=_degree(angular_degree);order=_order(degree,order)
    angular=_angular_basis(degree,order);X0=_points(points,len(angular['theta']))
    if caustic_encountered:raise ValueError('known upstream caustic: cut cannot be admitted as an unbroken sheet')
    determinant_min=None
    if ray_area_determinant is not None:
        determinant=np.asarray(ray_area_determinant,dtype=float)
        if determinant.shape!=(len(X0),) or not np.all(np.isfinite(determinant)) or np.min(determinant)<=0:
            raise ValueError('nonpositive/invalid optical determinant; caustic cut rejected')
        determinant_min=float(np.min(determinant))
    tolerance=float(expansion_tolerance)
    if not math.isfinite(tolerance) or tolerance<0:raise ValueError('finite nonnegative expansion tolerance required')
    coefficients=angular['values'].T@(angular['weights'][:,None]*X0)
    raw=angular['values']@coefficients
    first=np.einsum('nam,mi->nai',angular['first'],coefficients)
    second=np.einsum('nabm,mi->nabi',angular['second'],coefficients)
    X,dX,ddX,norm=_normalize_jets(raw,first,second)
    geometry=_physical_surface(X,dX,ddX,angular,provider,frame_scale,orientation=orientation,
                               outgoing_k=outgoing_k,immersion_tolerance=immersion_tolerance)
    error=np.linalg.norm(X-X0,axis=1)
    diagnostics={name:geometry[name] for name in ('immersion_min_gramdet','immersion_min_reference_density_squared',
        'outgoing_normality_max_relative','outgoing_null_defect_max_relative','outgoing_alignment_min')}
    diagnostics.update(embedding_interpolation_RMS=float(np.sqrt(np.dot(angular['weights'],error**2)/(4*math.pi))),
        embedding_interpolation_max=float(np.max(error)),fitted_embedding_norm_range=[float(np.min(norm)),float(np.max(norm))],
        sampled_nonlocal_coincidences=0,ray_area_determinant_min=determinant_min,
        immersion_threshold=float(immersion_tolerance),known_upstream_caustic=False)
    return ParametricCutEvaluation(X,geometry['theta_plus'],geometry['theta_minus'],geometry['normal'],geometry['screen'],
        geometry['b_screen'],geometry['K_screen'],geometry['area_weights'],degree,order,diagnostics,tolerance)


def initial_null_bundle(surface,provider,*,frame_scale,order=None,orientation=1,immersion_tolerance=1.e-18):
    """Exact graph jets supply both complete initial optical matrices.

    ``outgoing`` and ``ingoing`` are transport-ready mappings with X,k,
    spatial screen,D=I,P=B. The screen is orthonormal in the actual supplied
    spatial metric. Normals use n=(1/N,-beta/N), so k0 is future positive.
    """
    order=_order(surface.degree,order);angular=_angular_basis(surface.degree,order)
    rho,X,dX,ddX,radial=_embedding(surface,angular)
    geometry=_physical_surface(X,dX,ddX,angular,provider,frame_scale,orientation=orientation,
                               immersion_tolerance=immersion_tolerance)
    lapse,shift=_lapse_shift(geometry['data'],len(X))
    n=np.c_[1/lapse,-shift/lapse[:,None]]
    spatial=np.c_[np.zeros(len(X)),geometry['normal']]
    outgoing=dict(X=X.copy(),k=n+spatial,screen=geometry['screen'].copy(),
                  D=np.broadcast_to(np.eye(2),(len(X),2,2)).copy(),P=geometry['optical_plus'].copy())
    ingoing=dict(X=X.copy(),k=n-spatial,screen=geometry['screen'].copy(),
                 D=np.broadcast_to(np.eye(2),(len(X),2,2)).copy(),P=geometry['optical_minus'].copy())
    return dict(X=X,normal=geometry['normal'],screen=geometry['screen'],b_screen=geometry['b_screen'],
                K_screen=geometry['K_screen'],theta_plus=geometry['theta_plus'],theta_minus=geometry['theta_minus'],
                area_weights=geometry['area_weights'],area=float(np.sum(geometry['area_weights'])),
                outgoing=outgoing,ingoing=ingoing,order=order,orientation=orientation)
