"""FFT transforms for the unchanged real Wigner basis on unit S3.

The Hopf grid has X1=-sin(eta)sin(xi2), X2=-sin(eta)cos(xi2),
X3=cos(eta)sin(xi1), X4=cos(eta)cos(xi1). Consequently the existing
quaternion representation has b=sin(eta)exp(+i xi2), and its Fourier
indices are (r+s-l, s-r). This sign is tested against the regular-frame
Wigner evaluator. Labels, characters and derivative generators are unchanged.

No point-by-mode value matrix is built or stored. Evaluation accumulates
radial mode coefficients into Fourier bins; projection reverses the FFT
and contracts only the radial quadrature. No fitted transform, field filter,
filesystem cache or runtime setting is introduced.
"""
from dataclasses import dataclass,field
import math

import numpy as np
from numba import njit
from numpy.polynomial.legendre import leggauss
from scipy import sparse
from scipy.fft import fft2,ifft2

from .nsc_atlas_s3_frame import _level,VOLUME


def _integer(value,name,minimum=0):
    if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<minimum:
        raise ValueError(name+' must be an integer >= '+str(minimum))
    return int(value)


def _orders(degree,order):
    if order is None:return int(math.ceil(.75*degree))+2,3*degree+1
    if isinstance(order,(tuple,list)):
        if len(order)!=2:raise ValueError('quadrature tuple requires (n_eta,n_phi)')
        n=_integer(order[0],'n_eta',1);nphi=_integer(order[1],'n_phi',1)
    else:
        n=_integer(order,'quadrature_order',1);nphi=max(3*degree+1,2*n+1)
    if 2*n-1<degree or nphi<2*degree+1:
        raise ValueError('quadrature is below the two-product exactness requirement')
    return n,nphi


def _radial_modes(degree,u):
    """Normalized Sym^l([[cos eta,sin eta],[-sin eta,cos eta]]) entries."""
    cos=np.sqrt((1+u)/2);sin=np.sqrt((1-u)/2)
    M=(degree+1)*(degree+2)*(2*degree+3)//6
    radial=np.empty((len(u),M));k1=np.empty(M,dtype=np.int64);k2=np.empty_like(k1)
    previous=np.ones((len(u),1,1));offset=0
    for l in range(degree+1):
        if l:
            current=np.zeros((len(u),l+1,l+1))
            for r in range(l+1):
                for s in range(l+1):
                    value=current[:,r,s]
                    if r<l and s<l:value+=math.sqrt((l-r)*(l-s))/l*cos*previous[:,r,s]
                    if r<l and s:value+=math.sqrt((l-r)*s)/l*sin*previous[:,r,s-1]
                    if r and s<l:value-=math.sqrt(r*(l-s))/l*sin*previous[:,r-1,s]
                    if r and s:value+=math.sqrt(r*s)/l*cos*previous[:,r-1,s-1]
        else:current=previous
        normalization=math.sqrt((l+1)/VOLUME)
        for r in range(l+1):
            for s in range(l+1):
                radial[:,offset]=normalization*current[:,r,s]
                k1[offset]=r+s-l;k2[offset]=s-r
                offset+=1
        previous=current
    return radial,k1,k2


def _complex_conversion(labels):
    """A with real value row R=Y_complex A, hence C_complex=A c_real."""
    offsets={};count=0
    for l in range(max(label[0] for label in labels)+1):
        offsets[l]=count;count+=(l+1)**2
    rows=[];columns=[];data=[];invsqrt2=1/math.sqrt(2)
    for column,(l,r,s,kind) in enumerate(labels):
        i=offsets[l]+r*(l+1)+s;j=offsets[l]+(l-r)*(l+1)+(l-s)
        sign=1 if (r-s)%2==0 else -1
        if kind=='real':rows.append(i);columns.append(column);data.append(1.)
        elif kind=='cos':
            rows.extend([i,j]);columns.extend([column,column]);data.extend([invsqrt2,sign*invsqrt2])
        else:
            rows.extend([i,j]);columns.extend([column,column]);data.extend([-1j*invsqrt2,1j*sign*invsqrt2])
    A=sparse.csr_matrix((data,(rows,columns)),shape=(count,len(labels)),dtype=complex)
    A.sum_duplicates();A.sort_indices()
    return A


@njit(nogil=True)
def _accumulate(radial,coefficients,bin1,bin2,nphi):
    n_eta,M=radial.shape;fields=coefficients.shape[1]
    fourier=np.zeros((n_eta,nphi,nphi,fields),dtype=np.complex128)
    for a in range(n_eta):
        for mode in range(M):
            factor=radial[a,mode]
            for f in range(fields):
                value=coefficients[mode,f]
                if value!=0:
                    fourier[a,bin1[mode],bin2[mode],f]+=factor*value
    return fourier


@njit(nogil=True)
def _radial_projection(radial,radial_weights,fourier,bin1,bin2):
    n_eta,M=radial.shape;fields=fourier.shape[3]
    coefficients=np.zeros((M,fields),dtype=np.complex128)
    angular_volume=(2*math.pi)**2
    for mode in range(M):
        for a in range(n_eta):
            factor=angular_volume*radial_weights[a]*radial[a,mode]
            for f in range(fields):
                coefficients[mode,f]+=factor*fourier[a,bin1[mode],bin2[mode],f]
    return coefficients


@dataclass(frozen=True,slots=True)
class HopfFrameBasis:
    degree: int
    labels: tuple
    characters: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    generators: tuple
    _radial: np.ndarray=field(repr=False)
    _radial_weights: np.ndarray=field(repr=False)
    _conversion: object=field(repr=False)
    _bin1: np.ndarray=field(repr=False)
    _bin2: np.ndarray=field(repr=False)
    _n_phi: int=field(repr=False)

    @property
    def quadrature_shape(self):return (len(self._radial_weights),self._n_phi,self._n_phi)

    def evaluate(self,c):
        """Real coefficients (M,) or (M,F) -> grid samples (N,) or (N,F)."""
        c=np.asarray(c,dtype=float);vector=c.ndim==1
        if c.ndim not in (1,2) or c.shape[0]!=len(self.labels):
            raise ValueError('coefficients require shape (M,) or (M,F)')
        if vector:c=c[:,None]
        complex_coefficients=self._conversion@c
        fourier=_accumulate(self._radial,complex_coefficients,self._bin1,self._bin2,self._n_phi)
        samples=ifft2(fourier,axes=(1,2),overwrite_x=True)
        samples*=self._n_phi**2
        real=samples.real.copy().reshape(len(self.points),c.shape[1])
        return real[:,0] if vector else real

    def project(self,f):
        """Grid samples -> the same full-volume real-basis coefficients."""
        f=np.asarray(f,dtype=float);vector=f.ndim==1
        if f.ndim not in (1,2) or f.shape[0]!=len(self.points):
            raise ValueError('samples require shape (N,) or (N,F)')
        if vector:f=f[:,None]
        grid=f.reshape(self.quadrature_shape+(f.shape[1],))
        fourier=fft2(grid,axes=(1,2))
        fourier/=self._n_phi**2
        complex_coefficients=_radial_projection(self._radial,self._radial_weights,fourier,self._bin1,self._bin2)
        real=(self._conversion.conjugate().T@complex_coefficients).real.copy()
        return real[:,0] if vector else real


def build_hopf_frame_basis(degree,quadrature_order=None):
    """Build FFT/radial transforms with the regular-frame basis unchanged.

    Default: n_eta=ceil(3L/4)+2, n_phi=3L+1, for quadratic dealiasing.
    An integer n overrides n_eta and uses max(3L+1,2n+1) angular points.
    A tuple (n_eta,n_phi) controls both, subject to basis-product exactness.
    Nonpolynomial nonlinear products still require caller convergence checks.
    """
    degree=_integer(degree,'degree');n,nphi=_orders(degree,quadrature_order)
    levels=[_level(l) for l in range(degree+1)]
    labels=tuple(label for level in levels for label in level[0])
    characters=np.concatenate([level[1] for level in levels],axis=0)
    generators=tuple(sparse.block_diag([level[2][i] for level in levels],format='csr') for i in range(3))
    for D in generators:
        D.sum_duplicates();D.eliminate_zeros();D.sort_indices()
        for array in (D.data,D.indices,D.indptr):array.setflags(write=False)
    u,wu=leggauss(n);cos=np.sqrt((1+u)/2);sin=np.sqrt((1-u)/2)
    angle=2*math.pi*np.arange(nphi)/nphi
    points=np.empty((n,nphi,nphi,4))
    points[...,0]=-sin[:,None,None]*np.sin(angle)[None,None,:]
    points[...,1]=-sin[:,None,None]*np.cos(angle)[None,None,:]
    points[...,2]=cos[:,None,None]*np.sin(angle)[None,:,None]
    points[...,3]=cos[:,None,None]*np.cos(angle)[None,:,None]
    points=points.reshape(-1,4)
    radial_weights=wu/4
    weights=np.repeat(radial_weights,nphi*nphi)*(2*math.pi/nphi)**2
    radial,k1,k2=_radial_modes(degree,u)
    bin1=k1%nphi;bin2=k2%nphi;conversion=_complex_conversion(labels)
    for array in (characters,points,weights,radial,radial_weights,bin1,bin2):array.setflags(write=False)
    for array in (conversion.data,conversion.indices,conversion.indptr):array.setflags(write=False)
    return HopfFrameBasis(degree,labels,characters,points,weights,generators,
                          radial,radial_weights,conversion,bin1,bin2,nphi)
