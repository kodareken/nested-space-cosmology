"""Closed, deformed S2 null expansions on a supplied regular-frame S3 slice.

The provider supplies ``h, dh, K`` at unit embedding points X, where
dh[:, k, i, j] = E_k h_ij and [E_i,E_j] = 2*r*epsilon_ijk E_k.
K = -Lie_n(h)/2, so theta_plus = H_surface + K_ss - tr(K) and
theta_minus = -H_surface + K_ss - tr(K). The future null normals are n+s
and n-s; their common positive normalization is fixed by that convention.

Surfaces are genuine closed star-shaped graphs about an arbitrary center,
not a collection of independent round-sphere zeros. A bounded logistic of
real S2 harmonics gives their angular radius and analytic first/second
derivatives. The finder solves all angular residuals and then evaluates
the full expansions on a separate denser grid. Finite sampling, finite
shape degree and the supplied slice remain explicit limitations: a found
MOTS is not a proof of the outermost apparent horizon or its formation.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, gammaln, lpmv


DEFAULT_RADIUS_BOUNDS = (1.e-4, math.pi - 1.e-4)


def _degree(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
        raise ValueError('degree must be a nonnegative integer')
    return int(value)


def _unit_frame(X):
    """The linear quaternion fields A_i X; this also accepts tangent vectors."""
    x, y, z, w = np.asarray(X).T
    return np.stack((np.c_[w,z,-y,-x], np.c_[-z,w,x,-y], np.c_[y,-x,w,-z]), axis=1)


def _center(value):
    c = np.asarray(value, dtype=float)
    if c.shape != (4,) or not np.all(np.isfinite(c)) or np.linalg.norm(c) == 0:
        raise ValueError('center must be a finite nonzero four-vector')
    return c / np.linalg.norm(c)


@dataclass(frozen=True)
class StarShapedSurface:
    """Angular-radius graph; coefficients expand its bounded radius LOGIT.

    Labels are (l, m, kind), in degree order, m=0 then cos/sin pairs.
    Radius = lower + (upper-lower)*expit(sum coefficients*Y_lm).
    Bounds describe the searched chart, not new physical parameters.
    """
    center: np.ndarray
    coefficients: np.ndarray
    degree: int
    radius_bounds: tuple = DEFAULT_RADIUS_BOUNDS

    def __post_init__(self):
        degree = _degree(self.degree)
        c = _center(self.center).copy()
        coeff = np.asarray(self.coefficients, dtype=float).copy()
        if coeff.shape != ((degree+1)**2,) or not np.all(np.isfinite(coeff)):
            raise ValueError('coefficients must be finite with shape ((degree+1)**2,)')
        lo, hi = map(float, self.radius_bounds)
        if not (0 < lo < hi < math.pi):
            raise ValueError('radius bounds must satisfy 0 < lower < upper < pi')
        c.setflags(write=False); coeff.setflags(write=False)
        object.__setattr__(self, 'center', c)
        object.__setattr__(self, 'coefficients', coeff)
        object.__setattr__(self, 'degree', degree)
        object.__setattr__(self, 'radius_bounds', (lo, hi))


def round_surface(center, radius, degree=0, *, radius_bounds=DEFAULT_RADIUS_BOUNDS):
    """Constant-radius initial graph, optionally embedded in a larger shape span."""
    degree = _degree(degree)
    lo, hi = map(float, radius_bounds)
    radius = float(radius)
    if not lo < radius < hi:
        raise ValueError('initial radius must lie strictly inside radius bounds')
    fraction = (radius-lo)/(hi-lo)
    coeff = np.zeros((degree+1)**2)
    coeff[0] = math.sqrt(4*math.pi) * math.log(fraction/(1-fraction))
    return StarShapedSurface(center, coeff, degree, radius_bounds)


@lru_cache(maxsize=32)
def _angular_basis(degree, order):
    """Exact analytic harmonic jets on pole-free Gauss/periodic angular nodes."""
    ntheta, nphi = order
    z, wz = np.polynomial.legendre.leggauss(ntheta)
    theta = np.repeat(np.arccos(z), nphi)
    phi = np.tile(2*math.pi*np.arange(nphi)/nphi, ntheta)
    weights = np.repeat(wz, nphi)*(2*math.pi/nphi)
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    directions = np.c_[st*cp, st*sp, ct]
    first = np.stack((np.c_[ct*cp,ct*sp,-st], np.c_[-st*sp,st*cp,np.zeros_like(st)]), axis=1)
    second = np.empty((len(theta),2,2,3))
    second[:,0,0] = -directions
    second[:,0,1] = second[:,1,0] = np.c_[-ct*sp,ct*cp,np.zeros_like(st)]
    second[:,1,1] = np.c_[-st*cp,-st*sp,np.zeros_like(st)]
    labels = []; values = []; derivative = []; hessian = []
    for ell in range(degree+1):
        for m in range(ell+1):
            P = lpmv(m,ell,ct)
            previous = lpmv(m,ell-1,ct) if m < ell else np.zeros_like(ct)
            dP = (ell*ct*P-(ell+m)*previous)/st
            ddP = -ct/st*dP-(ell*(ell+1)-m*m/(st*st))*P
            norm = math.exp(.5*(math.log((2*ell+1)/(4*math.pi))+gammaln(ell-m+1)-gammaln(ell+m+1)))
            if m: norm *= math.sqrt(2)
            for kind in (('real',) if m == 0 else ('cos','sin')):
                if kind == 'sin': F=np.sin(m*phi); dF=m*np.cos(m*phi)
                else: F=np.cos(m*phi); dF=-m*np.sin(m*phi)
                labels.append((ell,m,kind))
                values.append(norm*P*F)
                derivative.append(np.stack((norm*dP*F,norm*P*dF),axis=1))
                hessian.append(np.stack((np.stack((norm*ddP*F,norm*dP*dF),axis=1),
                                         np.stack((norm*dP*dF,-m*m*norm*P*F),axis=1)),axis=1))
    return dict(theta=theta, phi=phi, weights=weights, directions=directions,
                direction_first=first, direction_second=second, labels=tuple(labels),
                values=np.stack(values,axis=1), first=np.stack(derivative,axis=2),
                second=np.stack(hessian,axis=3))


def _order(degree, value, *, validation=False):
    if value is None:
        n = max(10,2*degree+4)
        return (2*n+1,4*n+3) if validation else (n,2*n+1)
    if isinstance(value,(int,np.integer)) and not isinstance(value,(bool,np.bool_)):
        n=int(value); order=(n,2*n+1)
    else:
        if len(value)!=2: raise ValueError('quadrature order must be an integer or (n_theta,n_phi)')
        order=tuple(value)
    if any(isinstance(n,(bool,np.bool_)) or not isinstance(n,(int,np.integer)) or n<2 for n in order):
        raise ValueError('angular quadrature orders must be integers >= 2')
    if order[0] < degree+1 or order[1] < 2*degree+1:
        raise ValueError('quadrature must resolve the requested surface harmonics')
    return tuple(map(int,order))


def _embedding(surface, angular):
    u = angular['values'] @ surface.coefficients
    du = np.einsum('nam,m->na',angular['first'],surface.coefficients)
    ddu = np.einsum('nabm,m->nab',angular['second'],surface.coefficients)
    lo, hi = surface.radius_bounds
    f = expit(u); fp=(hi-lo)*f*(1-f); fpp=fp*(1-2*f)
    rho=lo+(hi-lo)*f
    drho=fp[:,None]*du
    ddrho=fp[:,None,None]*ddu+fpp[:,None,None]*du[:,:,None]*du[:,None,:]
    B = _unit_frame(surface.center[None,:])[0]
    v = angular['directions']@B
    dv = angular['direction_first']@B
    ddv = angular['direction_second']@B
    s,c=np.sin(rho),np.cos(rho)
    X=c[:,None]*surface.center+s[:,None]*v
    radial=-s[:,None]*surface.center+c[:,None]*v
    dX=radial[:,None,:]*drho[:,:,None]+s[:,None,None]*dv
    ddX=(radial[:,None,None,:]*ddrho[:,:,:,None]
         -X[:,None,None,:]*drho[:,:,None,None]*drho[:,None,:,None]
         +c[:,None,None,None]*(drho[:,:,None,None]*dv[:,None,:,:]+drho[:,None,:,None]*dv[:,:,None,:])
         +s[:,None,None,None]*ddv)
    return rho,X,dX,ddX,radial


def _geometry(provider, points):
    data=provider(points)
    if isinstance(data,dict): h,dh,K=(data[name] for name in ('h','dh','K'))
    else: h,dh,K=data
    arrays=[]
    for value,name,shape in ((h,'h',(len(points),3,3)),(dh,'dh',(len(points),3,3,3)),(K,'K',(len(points),3,3))):
        a=np.asarray(value,dtype=float)
        if a.shape!=shape or not np.all(np.isfinite(a)):
            raise ValueError(name+' must be a finite array of shape '+str(shape))
        if not np.allclose(a,np.swapaxes(a,-1,-2),rtol=2e-12,atol=2e-14):
            raise ValueError(name+' must be symmetric in its last two indices')
        arrays.append(a)
    h,dh,K=arrays
    try: np.linalg.cholesky(h)
    except np.linalg.LinAlgError as exc: raise ValueError('spatial metric must be positive definite') from exc
    return h,dh,K


def _connection(h, dh, r):
    eps=np.zeros((3,3,3))
    for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
        eps[i,j,k]=1.;eps[j,i,k]=-1.
    C=2*r*eps
    lower=np.empty((len(h),3,3,3))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                value=dh[:,i,j,k]+dh[:,j,i,k]-dh[:,k,i,j]
                for ell in range(3):
                    value=value-h[:,i,ell]*C[j,k,ell]+h[:,j,ell]*C[k,i,ell]+h[:,k,ell]*C[i,j,ell]
                lower[:,k,i,j]=value/2
    inverse=np.linalg.inv(h)
    return inverse,np.einsum('nlk,nkij->nlij',inverse,lower)


@dataclass(frozen=True)
class SurfaceEvaluation:
    points: np.ndarray
    radius: np.ndarray
    theta_plus: np.ndarray
    theta_minus: np.ndarray
    mean_curvature: np.ndarray
    area_weights: np.ndarray
    quadrature_order: tuple
    expansion_tolerance: float
    normal: np.ndarray | None = None

    @property
    def outward_normal(self):
        """Physical unit spatial normal, contravariant in the supplied E frame."""
        return self.normal

    @property
    def area(self): return float(np.sum(self.area_weights))

    @property
    def plus_rms(self): return float(np.sqrt(np.dot(self.area_weights,self.theta_plus**2)/self.area))

    @property
    def minus_rms(self): return float(np.sqrt(np.dot(self.area_weights,self.theta_minus**2)/self.area))

    @property
    def plus_max_abs(self): return float(np.max(abs(self.theta_plus)))

    @property
    def future_trapped_sampled(self):
        return bool(np.max(self.theta_plus)<-self.expansion_tolerance and np.max(self.theta_minus)<-self.expansion_tolerance)

    @property
    def future_mots_sampled(self):
        return bool(self.plus_max_abs<=self.expansion_tolerance and np.max(self.theta_minus)<-self.expansion_tolerance)

    def summary(self):
        return dict(area=self.area,theta_plus_rms=self.plus_rms,theta_minus_rms=self.minus_rms,
                    theta_plus_max_abs=self.plus_max_abs,theta_plus_range=[float(np.min(self.theta_plus)),float(np.max(self.theta_plus))],
                    theta_minus_range=[float(np.min(self.theta_minus)),float(np.max(self.theta_minus))],
                    radius_range=[float(np.min(self.radius)),float(np.max(self.radius))],
                    future_trapped_sampled=self.future_trapped_sampled,future_mots_sampled=self.future_mots_sampled,
                    expansion_tolerance=self.expansion_tolerance,quadrature_order=list(self.quadrature_order),
                    scope='Finite angular samples of a supplied slice; no continuum sign bound, outermost-horizon or formation theorem.')


def evaluate_surface(surface, provider, *, frame_scale, quadrature_order=None, expansion_tolerance=1.e-7):
    """Physical expansions and area; raw frame derivatives, not covariant dh."""
    r=float(frame_scale); tolerance=float(expansion_tolerance)
    if not math.isfinite(r) or r<=0 or not math.isfinite(tolerance) or tolerance<0:
        raise ValueError('positive frame scale and finite nonnegative expansion tolerance required')
    order=_order(surface.degree,quadrature_order)
    angular=_angular_basis(surface.degree,order)
    rho,X,dX,ddX,radial=_embedding(surface,angular)
    h,dh,K=_geometry(provider,X)
    hi,Gamma=_connection(h,dh,r)
    frame=_unit_frame(X)
    V=np.einsum('nik,nak->nai',frame,dX)/r
    dV=np.empty((len(X),2,2,3))
    for a in range(2):
        dframe=_unit_frame(dX[:,a])
        for b in range(2):
            dV[:,a,b]=(np.einsum('nik,nk->ni',dframe,dX[:,b])+np.einsum('nik,nk->ni',frame,ddX[:,a,b]))/r
    induced=np.einsum('nai,nij,nbj->nab',V,h,V)
    area_density=np.sqrt(np.linalg.det(induced))
    covnormal=np.cross(V[:,0],V[:,1])
    radial_components=np.einsum('nik,nk->ni',frame,radial)/r
    covnormal*=np.where(np.einsum('ni,ni->n',covnormal,radial_components)>0,1.,-1.)[:,None]
    normal=np.einsum('nij,nj->ni',hi,covnormal)
    norm=np.sqrt(np.einsum('ni,ni->n',covnormal,normal))
    normal/=norm[:,None];covnormal/=norm[:,None]
    acceleration=dV+np.einsum('nkij,nai,nbj->nabk',Gamma,V,V)
    H=-np.einsum('nab,nk,nabk->n',np.linalg.inv(induced),covnormal,acceleration)
    transverse=np.einsum('nij,ni,nj->n',K,normal,normal)-np.einsum('nij,nji->n',hi,K)
    area_weights=angular['weights']*area_density/np.sin(angular['theta'])
    plus,minus=H+transverse,-H+transverse
    if not all(np.all(np.isfinite(a)) for a in (area_weights,plus,minus)) or np.min(area_weights)<=0:
        raise ValueError('surface is not a regular finite positive-area graph')
    return SurfaceEvaluation(X,rho,plus,minus,H,area_weights,order,tolerance,normal=normal)


@dataclass(frozen=True)
class MOTSResult:
    surface: StarShapedSurface
    evaluation: SurfaceEvaluation
    collocation_rms: float
    optimizer_success: bool
    nfev: int
    message: str

    @property
    def converged(self): return self.evaluation.plus_max_abs<=self.evaluation.expansion_tolerance

    @property
    def future_mots(self): return self.evaluation.future_mots_sampled

    def summary(self):
        return dict(self.evaluation.summary(),converged=self.converged,optimizer_success=self.optimizer_success,
                    collocation_rms=self.collocation_rms,nfev=self.nfev,message=self.message,
                    shape_degree=self.surface.degree,shape_coefficients=self.surface.coefficients.tolist(),
                    center=self.surface.center.tolist(),radius_bounds=list(self.surface.radius_bounds))


def find_mots(provider, *, center=None, initial_radius=None, frame_scale, degree=None,
              radius_bounds=None,quadrature_order=None,validation_order=None,
              max_nfev=100,expansion_tolerance=1.e-7,initial_surface=None):
    """Solve theta_plus=0 in a deformed shape span, then independently resample.

    ``converged`` refers to the denser sampled expansion residual. ``future_mots``
    additionally requires theta_minus < -tolerance everywhere sampled. A
    cosmological marginal sphere can solve theta_plus=0 while failing the latter.
    Failure to find a root excludes no surfaces outside this chart/span/start.
    ``initial_surface`` warm-starts its complete deformed graph on another
    metric. Omit initial_radius; explicitly supplied center/degree/bounds
    must agree. Without it, the old round seed defaults remain degree4 and
    DEFAULT_RADIUS_BOUNDS.
    """
    if initial_surface is None:
        if center is None or initial_radius is None:
            raise ValueError('center and initial_radius are required without initial_surface')
        initial=round_surface(center,initial_radius,4 if degree is None else degree,
                              radius_bounds=DEFAULT_RADIUS_BOUNDS if radius_bounds is None else radius_bounds)
    else:
        if not isinstance(initial_surface,StarShapedSurface):
            raise ValueError('initial_surface must be a StarShapedSurface')
        if initial_radius is not None:
            raise ValueError('initial_radius cannot be combined with initial_surface')
        if center is not None and not np.allclose(_center(center),initial_surface.center,rtol=0,atol=2.e-14):
            raise ValueError('center conflicts with initial_surface')
        if degree is not None and _degree(degree)!=initial_surface.degree:
            raise ValueError('degree conflicts with initial_surface')
        if radius_bounds is not None and tuple(map(float,radius_bounds))!=initial_surface.radius_bounds:
            raise ValueError('radius bounds conflict with initial_surface')
        initial=initial_surface
    order=_order(initial.degree,quadrature_order)
    vorder=_order(initial.degree,validation_order,validation=validation_order is None)
    if vorder[0]<=order[0] or vorder[1]<=order[1]:
        raise ValueError('validation quadrature must be strictly denser in both angles')
    angular=_angular_basis(initial.degree,order)
    weights=np.sqrt(angular['weights']/(4*math.pi))
    def residual(coeff):
        trial=StarShapedSurface(initial.center,coeff,initial.degree,initial.radius_bounds)
        evaluation=evaluate_surface(trial,provider,frame_scale=frame_scale,quadrature_order=order,
                                    expansion_tolerance=expansion_tolerance)
        return weights*evaluation.theta_plus
    optimum=least_squares(residual,initial.coefficients,method='trf',max_nfev=max_nfev,
                          ftol=1.e-11,xtol=1.e-11,gtol=1.e-11)
    surface=StarShapedSurface(initial.center,optimum.x,initial.degree,initial.radius_bounds)
    final=evaluate_surface(surface,provider,frame_scale=frame_scale,quadrature_order=vorder,
                           expansion_tolerance=expansion_tolerance)
    return MOTSResult(surface,final,float(np.linalg.norm(optimum.fun)),bool(optimum.success),
                      int(optimum.nfev),str(optimum.message))
