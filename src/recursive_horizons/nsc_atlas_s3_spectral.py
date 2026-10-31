"""Real scalar harmonics and full-volume quadrature on the unit round S3.

The retained symmetry sector has even j, cosine azimuthal modes and even
l+m. It excludes the l=1 conformal kernel. Orthonormality uses the actual
volume 2*pi**2, not probability-normalized volume. Labels and quadrature
orders are mathematical approximation data, not selected physical cutoffs.

Coordinate convention: X4=u, X1=sqrt(1-u**2)*v, with X2/X3 the cosine/sine
azimuthal components. Associated Legendre functions include SciPy's
Condon--Shortley phase. No filesystem cache or runtime setting is changed.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_gegenbauer,gammaln,lpmv,roots_jacobi


def _positive_order(value,name):
    if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<1:
        raise ValueError(name+' must be a positive integer')
    return int(value)


def _degree(value):
    if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<0:
        raise ValueError('degree must be a nonnegative integer')
    return int(value)


def _labels(degree):
    return tuple((l,j,m) for l in range(degree+1) for j in range(0,l+1,2)
                 for m in range(j+1) if (l+m)%2==0)


def _validated_labels(labels):
    result=[]
    for label in labels:
        if len(label)!=3 or any(isinstance(n,bool) or not isinstance(n,(int,np.integer)) for n in label):
            raise ValueError('harmonic labels must be integer triples (l,j,m)')
        l,j,m=map(int,label)
        if not 0<=m<=j<=l:
            raise ValueError('harmonic labels require 0 <= m <= j <= l')
        result.append((l,j,m))
    return tuple(result)


@lru_cache(maxsize=2048)
def _radial_normalization(l,j):
    n=l-j;alpha=j+1
    log_h=(math.log(math.pi)+(1-2*alpha)*math.log(2)+gammaln(n+2*alpha)
           -gammaln(n+1)-math.log(n+alpha)-2*gammaln(alpha))
    return math.exp(-float(log_h)/2)


@lru_cache(maxsize=1024)
def _angular_normalization(j,m):
    log_n=(math.log(2*j+1)-math.log(4*math.pi)+gammaln(j-m+1)-gammaln(j+m+1))/2
    if m:log_n+=math.log(2)/2
    return math.exp(float(log_n))


def _factors(u,v,phi,labels):
    radial={};angular={};azimuth={}
    one_minus_u2=np.maximum(0.,1-u*u)
    for l,j,m in labels:
        if (l,j) not in radial:
            radial[l,j]=(_radial_normalization(l,j)*one_minus_u2**(j/2)
                          *eval_gegenbauer(l-j,j+1,u))
        if (j,m) not in angular:
            angular[j,m]=_angular_normalization(j,m)*lpmv(m,j,v)
        if m not in azimuth:azimuth[m]=np.cos(m*phi)
    return radial,angular,azimuth


def _nodes(n_u,n_v,n_phi):
    n_u=_positive_order(n_u,'n_u');n_v=_positive_order(n_v,'n_v')
    n_phi=_positive_order(n_phi,'n_phi')
    u,wu=roots_jacobi(n_u,.5,.5)
    v,wv=leggauss(n_v)
    phi=2*math.pi*np.arange(n_phi)/n_phi
    return u,wu,v,wv,phi


def _quadrature_from_nodes(u,wu,v,wv,phi):
    nu,nv,nphi=len(u),len(v),len(phi)
    ru=np.sqrt(np.maximum(0.,1-u*u));rv=np.sqrt(np.maximum(0.,1-v*v))
    points=np.empty((nu,nv,nphi,4),dtype=float)
    points[...,0]=ru[:,None,None]*v[None,:,None]
    points[...,1]=ru[:,None,None]*rv[None,:,None]*np.cos(phi)[None,None,:]
    points[...,2]=ru[:,None,None]*rv[None,:,None]*np.sin(phi)[None,None,:]
    points[...,3]=u[:,None,None]
    weights=np.repeat((wu[:,None]*wv[None,:]).ravel(),nphi)*(2*math.pi/nphi)
    return points.reshape(-1,4),weights


def sphere_quadrature(n_u,n_v,n_phi):
    """Return points(N,4), weights(N) for unnormalized unit-S3 volume.

    Jacobi(1/2,1/2), Legendre and periodic azimuthal quadratures are used.
    Weights sum to 2*pi**2. This allocates no harmonic-value matrix.
    """
    u,wu,v,wv,phi=_nodes(n_u,n_v,n_phi)
    return _quadrature_from_nodes(u,wu,v,wv,phi)


def harmonic_values(points,labels):
    """Evaluate supplied real cosine harmonics at unit points, in one batch.

    The caller controls batch size (e.g. 8192 points) for high quadrature
    validation. Only the small label normalization constants are cached.
    All point-dependent arrays are local to this call.
    """
    points=np.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[1]!=4 or not np.all(np.isfinite(points)):
        raise ValueError('points must be a finite (N,4) array on the unit sphere')
    if np.any(np.abs(np.sum(points*points,axis=1)-1)>2e-10):
        raise ValueError('harmonic evaluation requires unit-sphere points')
    labels=_validated_labels(labels)
    result=np.empty((len(points),len(labels)),dtype=float)
    if not len(points) or not labels:return result
    u=np.clip(points[:,3],-1.,1.);radius=np.sqrt(np.maximum(0.,1-u*u))
    v=np.zeros_like(u)
    np.divide(points[:,0],radius,out=v,where=radius>0)
    v=np.clip(v,-1.,1.)
    phi=np.arctan2(points[:,2],points[:,1])
    radial,angular,azimuth=_factors(u,v,phi,labels)
    for column,(l,j,m) in enumerate(labels):
        result[:,column]=radial[l,j]*angular[j,m]*azimuth[m]
    return result


@dataclass(frozen=True,slots=True)
class SurfaceBasis:
    degree: int
    labels: tuple
    eigenvalues: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    values: np.ndarray

    def _coefficients(self,c):
        c=np.asarray(c,dtype=float)
        if c.ndim not in (1,2) or c.shape[0]!=len(self.labels):
            raise ValueError('coefficients require shape (M,) or (M,K)')
        return c

    def evaluate(self,c):
        """Evaluate one expansion or a column batch of expansions."""
        return self.values@self._coefficients(c)

    def project(self,f):
        """Volume-quadrature projection; no fitted Gram correction is used."""
        f=np.asarray(f,dtype=float)
        if f.ndim not in (1,2) or f.shape[0]!=len(self.weights):
            raise ValueError('field samples require shape (N,) or (N,K)')
        weighted=self.weights*f if f.ndim==1 else self.weights[:,None]*f
        return self.values.T@weighted

    def laplacian_values(self,c):
        """Evaluate the intrinsic scalar Laplacian of the expansion."""
        c=self._coefficients(c)
        scaled=-self.eigenvalues*c if c.ndim==1 else -self.eigenvalues[:,None]*c
        return self.values@scaled


def build_surface_basis(degree,*,n_u=None,n_v=None,n_phi=None):
    """Build the symmetry-restricted orthonormal basis and its quadrature.

    Defaults n_u=n_v=3L+3 and n_phi=6L+5 integrate polynomial products
    through degree 6L+4 (up to floating-point error). User-supplied lower
    orders do not inherit that exactness. Through L=8 there are45 modes;
    L=16 has249 and its default values array is about523MB (decimal).
    """
    degree=_degree(degree);labels=_labels(degree)
    u,wu,v,wv,phi=_nodes(3*degree+3 if n_u is None else n_u,
                        3*degree+3 if n_v is None else n_v,
                        6*degree+5 if n_phi is None else n_phi)
    points,weights=_quadrature_from_nodes(u,wu,v,wv,phi)
    values=np.empty((len(u),len(v),len(phi),len(labels)),dtype=float)
    radial,angular,azimuth=_factors(u,v,phi,labels)
    # One u block stays small; all v/phi points are evaluated together.
    # This avoids repeated special-function calls at the tensor grid points.
    for iu in range(len(u)):
        block=values[iu]
        for column,(l,j,m) in enumerate(labels):
            block[...,column]=radial[l,j][iu]*angular[j,m][:,None]*azimuth[m][None,:]
    values=values.reshape(len(points),len(labels))
    eigenvalues=np.array([l*(l+2) for l,_,_ in labels],dtype=float)
    for array in (eigenvalues,points,weights,values):array.setflags(write=False)
    return SurfaceBasis(degree,labels,eigenvalues,points,weights,values)
