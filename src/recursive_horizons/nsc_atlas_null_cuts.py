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


def _spacetime_surface(times,dt,ddt,X,dX,ddX,angular,provider_events,r,outgoing_k,
                       *,immersion_tolerance=1.e-18):
    """Event/jet core; Gamma[a,b,c] is (nabla_(E_b) E_c)^a.

    A separate core permits independent local-patch sign controls, while
    the public interface below always reconstructs one complete S2 map.
    """
    r=float(r);tolerance=float(immersion_tolerance)
    if not math.isfinite(r) or r<=0 or not math.isfinite(tolerance) or tolerance<0:
        raise ValueError('positive frame scale and nonnegative immersion tolerance required')
    data=provider_events(times,X)
    if not isinstance(data,dict):raise ValueError('event provider must return g and Gamma in a mapping')
    g=np.asarray(data['g'],dtype=float);Gamma=np.asarray(data['Gamma'],dtype=float)
    count=len(X)
    if g.shape!=(count,4,4) or Gamma.shape!=(count,4,4,4) or not np.all(np.isfinite(g)) or not np.all(np.isfinite(Gamma)):
        raise ValueError('finite event g(N,4,4) and physical Gamma(N,4,4,4) required')
    if not np.allclose(g,g.swapaxes(1,2),rtol=2.e-12,atol=2.e-14):
        raise ValueError('event metric must be symmetric')
    try:
        np.linalg.cholesky(g[:,1:,1:]);inverse=np.linalg.inv(g)
    except np.linalg.LinAlgError as exc:raise ValueError('event metric has invalid spatial slicing or inverse') from exc
    if np.any(inverse[:,0,0]>=0):raise ValueError('event slices must have a positive real lapse')
    frame=_unit_frame(X)
    T=np.empty((count,2,4));T[:,:,0]=dt
    T[:,:,1:]=np.einsum('nik,nak->nai',frame,dX)/r
    dT=np.empty((count,2,2,4));dT[:,:,:,0]=ddt
    for a in range(2):
        dframe=_unit_frame(dX[:,a])
        for b in range(2):
            dT[:,a,b,1:]=(np.einsum('nik,nk->ni',dframe,dX[:,b])
                          +np.einsum('nik,nk->ni',frame,ddX[:,a,b]))/r
    induced=np.einsum('nAi,nij,nBj->nAB',T,g,T)
    try:np.linalg.cholesky(induced)
    except np.linalg.LinAlgError as exc:raise ValueError('spacetime cut must be spacelike everywhere sampled') from exc
    determinant=np.linalg.det(induced)
    regular_determinant=r**4*determinant/np.sin(angular['theta'])**2
    if np.min(regular_determinant)<=tolerance:
        raise ValueError('spacetime cut reaches the sampled immersion threshold')
    qinverse=np.linalg.inv(induced)
    # Orthogonal projection subtracts the complete spacetime tangent plane.
    def normal_projection(vector):
        tangent_dot=np.einsum('nAi,nij,nj->nA',T,g,vector)
        coefficients=np.einsum('nAB,nB->nA',qinverse,tangent_dot)
        return vector-np.einsum('nA,nAi->ni',coefficients,T)
    lapse=1/np.sqrt(-inverse[:,0,0])
    Eulerian=-lapse[:,None]*inverse[:,:,0]
    projected=normal_projection(Eulerian)
    square=np.einsum('ni,nij,nj->n',projected,g,projected)
    if np.any(square>=0):raise ValueError('normal-plane projection did not retain a timelike future normal')
    u=projected/np.sqrt(-square)[:,None]
    if np.min(u[:,0])<=0:raise ValueError('normal-plane timelike normal is not future directed')
    k=np.asarray(outgoing_k,dtype=float)
    if k.shape!=(count,4) or not np.all(np.isfinite(k)) or np.min(k[:,0])<=0:
        raise ValueError('outgoing_k must be finite future-directed (N,4) frame components')
    kperp=normal_projection(k)
    energy=-np.einsum('ni,nij,nj->n',u,g,k)
    if np.min(energy)<=0:raise ValueError('outgoing generator has nonpositive normal-plane energy')
    vraw=kperp-energy[:,None]*u
    vsquare=np.einsum('ni,nij,nj->n',vraw,g,vraw)
    if np.any(vsquare<=1.e-20*energy**2):
        raise ValueError('outgoing generator does not select a spacelike normal direction')
    v=vraw/np.sqrt(vsquare)[:,None]
    plus=u+v;minus=u-v
    acceleration=dT+np.einsum('nijk,nAj,nBk->nABi',Gamma,T,T)
    mean_vector=np.einsum('nAB,nABi->ni',qinverse,acceleration)
    theta_plus=-np.einsum('ni,nij,nj->n',plus,g,mean_vector)
    theta_minus=-np.einsum('ni,nij,nj->n',minus,g,mean_vector)
    tangent_k=k-kperp
    normality=np.sqrt(np.maximum(0,np.einsum('ni,nij,nj->n',tangent_k,g,tangent_k)))/energy
    null_defect=abs(np.einsum('ni,nij,nj->n',k,g,k))/energy**2
    alignment=np.sqrt(vsquare)/energy
    screen_transform=np.zeros((count,2,2));screen_transform[:,0,0]=1/np.sqrt(induced[:,0,0])
    schur=induced[:,1,1]-induced[:,0,1]**2/induced[:,0,0]
    screen_transform[:,1,1]=1/np.sqrt(schur)
    screen_transform[:,1,0]=-induced[:,0,1]/induced[:,0,0]*screen_transform[:,1,1]
    screen=np.einsum('nAa,nai->nAi',screen_transform,T)
    area_weights=angular['weights']*np.sqrt(determinant)/np.sin(angular['theta'])
    return dict(theta_plus=theta_plus,theta_minus=theta_minus,normal_u=u,normal_v=v,
        null_plus=plus,null_minus=minus,screen=screen,area_weights=area_weights,normal_energy=energy,
        immersion_min_gramdet=float(np.min(determinant)),
        immersion_min_scaled_density_squared=float(np.min(regular_determinant)),
        outgoing_normality_max_relative=float(np.max(normality)),outgoing_null_defect_max_relative=float(np.max(null_defect)),
        outgoing_alignment_range=[float(np.min(alignment)),float(np.max(alignment))],
        lapse_range=[float(np.min(lapse)),float(np.max(lapse))])


@dataclass(frozen=True)
class SpacetimeCutEvaluation:
    times: np.ndarray
    points: np.ndarray
    theta_plus: np.ndarray
    theta_minus: np.ndarray
    normal_u: np.ndarray
    normal_v: np.ndarray
    null_plus: np.ndarray
    null_minus: np.ndarray
    screen: np.ndarray
    normal_energy: np.ndarray
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
            normal_energy_range=[float(np.min(self.normal_energy)),float(np.max(self.normal_energy))],
            future_trapped_sampled=self.future_trapped_sampled,future_mots_sampled=self.future_mots_sampled,
            angular_degree=self.angular_degree,order=list(self.order),expansion_tolerance=self.expansion_tolerance,
            global_injectivity_certified=False,continuum_expansion_signs_certified=False,child_formation_proved=False,
            scope='One closed spacetime cut, finite event/harmonic reconstruction and angular samples. Normals are projected-Eulerian u +/- v; no global injectivity or child admission.')


def evaluate_spacetime_cut(times,points,provider_events,*,frame_scale,angular_degree,order,outgoing_k,
                           caustic_encountered=False,ray_area_determinant=None,
                           immersion_tolerance=1.e-18,expansion_tolerance=1.e-7):
    """Measure an angularly time-tilted closed cut with full spacetime jets.

    provider_events(t_N,X_N) returns g and the actual physical frame
    connection Gamma[a,b,c]=(nabla_(E_b)E_c)^a. The reconstructed times are
    never clamped to history bounds. A common-time input remains exactly
    common-time, avoiding floating projection noise outside an endpoint.
    ``normal_energy`` is -u.k, with u the projected-Eulerian normal, and
    is the affine/unit-cut expansion factor when k is normal and null.
    """
    degree=_degree(angular_degree);order=_order(degree,order)
    angular=_angular_basis(degree,order);count=len(angular['theta']);X0=_points(points,count)
    t0=np.asarray(times,dtype=float)
    if t0.shape!=(count,) or not np.all(np.isfinite(t0)):
        raise ValueError('cut times require finite shape (n_theta*n_phi,)')
    if caustic_encountered:raise ValueError('known upstream caustic: spacetime cut is not an unbroken sheet')
    determinant_min=None
    if ray_area_determinant is not None:
        determinant=np.asarray(ray_area_determinant,dtype=float)
        if determinant.shape!=(count,) or not np.all(np.isfinite(determinant)) or np.min(determinant)<=0:
            raise ValueError('nonpositive/invalid optical determinant; spacetime caustic cut rejected')
        determinant_min=float(np.min(determinant))
    tolerance=float(expansion_tolerance)
    if not math.isfinite(tolerance) or tolerance<0:raise ValueError('finite nonnegative expansion tolerance required')
    weights=angular['weights'];Y=angular['values']
    if np.all(t0==t0[0]):
        t=t0.copy();dt=np.zeros((count,2));ddt=np.zeros((count,2,2))
    else:
        tc=Y.T@(weights*t0);t=Y@tc
        dt=np.einsum('nam,m->na',angular['first'],tc)
        ddt=np.einsum('nabm,m->nab',angular['second'],tc)
    coefficients=Y.T@(weights[:,None]*X0)
    raw=Y@coefficients
    first=np.einsum('nam,mi->nai',angular['first'],coefficients)
    second=np.einsum('nabm,mi->nabi',angular['second'],coefficients)
    X,dX,ddX,norm=_normalize_jets(raw,first,second)
    events=np.c_[float(frame_scale)*t,X]
    coincident=cKDTree(events).query_pairs(1.e-10,output_type='ndarray')
    if len(coincident):
        separation=np.linalg.norm(angular['directions'][coincident[:,0]]-angular['directions'][coincident[:,1]],axis=1)
        if np.any(separation>1.e-5):raise ValueError('sampled nonlocal event coincidence; self-intersecting cut rejected')
    geometry=_spacetime_surface(t,dt,ddt,X,dX,ddX,angular,provider_events,frame_scale,outgoing_k,
                                immersion_tolerance=immersion_tolerance)
    error=np.linalg.norm(X-X0,axis=1);terror=abs(t-t0)
    diagnostics={name:geometry[name] for name in ('immersion_min_gramdet','immersion_min_scaled_density_squared',
        'outgoing_normality_max_relative','outgoing_null_defect_max_relative','outgoing_alignment_range','lapse_range')}
    diagnostics.update(embedding_interpolation_RMS=float(np.sqrt(np.dot(weights,error**2)/(4*math.pi))),
        embedding_interpolation_max=float(np.max(error)),time_interpolation_RMS=float(np.sqrt(np.dot(weights,terror**2)/(4*math.pi))),
        time_interpolation_max=float(np.max(terror)),time_range=[float(np.min(t)),float(np.max(t))],
        fitted_embedding_norm_range=[float(np.min(norm)),float(np.max(norm))],ray_area_determinant_min=determinant_min,
        sampled_nonlocal_event_coincidences=0,known_upstream_caustic=False,immersion_threshold=float(immersion_tolerance))
    return SpacetimeCutEvaluation(t,X,geometry['theta_plus'],geometry['theta_minus'],geometry['normal_u'],geometry['normal_v'],
        geometry['null_plus'],geometry['null_minus'],geometry['screen'],geometry['normal_energy'],geometry['area_weights'],
        degree,order,diagnostics,tolerance)
