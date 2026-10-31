"""Full real Wigner harmonics with sparse regular unit-S3 frame derivatives.

U(q)=w I-i(x sigma1+y sigma2+z sigma3) represents q=(x,y,z,w).
The frame e_i(q)=q i_i has bracket [e_i,e_j]=2 epsilon_ijk e_k.
For a point-by-mode value matrix Y, e_i Y=Y D_i; coefficients therefore
differentiate as c -> D_i c, without a transpose.

Only representation formulas and the existing sphere quadrature are used.
There is no fitted orthogonalization, eigendecomposition, filesystem cache,
or runtime setting. D_i differentiate scalar component functions; tensor
covariant derivatives still require the regular-frame connection.
"""
from dataclasses import dataclass
from functools import lru_cache
import math

import numpy as np
from numba import njit
from scipy import sparse

from .nsc_atlas_s3_spectral import sphere_quadrature

VOLUME=2*math.pi**2


def _degree(value):
    if isinstance(value,bool) or not isinstance(value,(int,np.integer)) or value<0:
        raise ValueError('degree must be a nonnegative integer')
    return int(value)


@lru_cache(maxsize=32)
def _level(l):
    """Real-pair labels, G+ characters, and exact-formula derivative blocks."""
    size=l+1;labels=[];characters=[];rows=[];columns=[];data=[]
    invsqrt2=1/math.sqrt(2)
    for r in range(size):
        for s in range(size):
            i=r*size+s;j=(l-r)*size+(l-s)
            if i>j:continue
            sign=1 if (r-s)%2==0 else -1
            if i==j:
                column=len(labels);labels.append((l,r,s,'real'));characters.append((1,1))
                rows.append(i);columns.append(column);data.append(1.)
            else:
                cosine=len(labels);sine=cosine+1
                labels.extend([(l,r,s,'cos'),(l,r,s,'sin')])
                characters.extend([(sign,(-1)**l*sign),(sign,-(-1)**l*sign)])
                rows.extend([i,j,i,j]);columns.extend([cosine,cosine,sine,sine])
                data.extend([invsqrt2,sign*invsqrt2,-1j*invsqrt2,1j*sign*invsqrt2])
    A=sparse.csr_matrix((data,(rows,columns)),shape=(size*size,size*size),dtype=complex)
    if l:
        r=np.arange(1,l+1)
        raising=sparse.diags(np.sqrt(r*(l-r+1)),1,shape=(size,size),format='csr',dtype=complex)
    else:raising=sparse.csr_matrix((1,1),dtype=complex)
    lowering=raising.T
    generators=(-1j*(raising+lowering),lowering-raising,
                -1j*sparse.diags(l-2*np.arange(size),format='csr',dtype=complex))
    result=[]
    for T in generators:
        D=sparse.kron(sparse.eye(size,format='csr'),T,format='csr')
        real=(A.conjugate().T@D@A).tocsr();real.sum_duplicates()
        if real.nnz and np.max(abs(real.data.imag))>64*np.finfo(float).eps*(l+1):
            raise ArithmeticError('real harmonic conversion acquired an imaginary generator')
        real=real.real.tocsr();real.eliminate_zeros();real.sort_indices()
        result.append(real)
    return tuple(labels),np.asarray(characters,dtype=np.int8),tuple(result)


def _validated_labels(labels):
    result=[]
    for label in labels:
        if len(label)!=4:raise ValueError('frame labels require (l,r,s,kind)')
        l,r,s,kind=label
        if any(isinstance(n,bool) or not isinstance(n,(int,np.integer)) for n in (l,r,s)):
            raise ValueError('frame degree and matrix indices must be integers')
        l,r,s=map(int,(l,r,s))
        if not 0<=r<=l or not 0<=s<=l or kind not in ('cos','sin','real'):
            raise ValueError('invalid real Wigner harmonic label')
        i=r*(l+1)+s;j=(l-r)*(l+1)+(l-s)
        if (kind=='real' and i!=j) or (kind!='real' and i>=j):
            raise ValueError('label must name its canonical conjugate-pair representative')
        result.append((l,r,s,kind))
    return tuple(result)


def arbitrary_points_values(points,labels):
    """Real normalized Wigner harmonics at unit points, in supplied label order.

    Sym^l(U) is built in the normalized symmetric-tensor basis by the
    maximal-spin Clebsch--Gordan recursion. Point batches and label subsets
    are supported; point-dependent arrays are not cached.
    """
    points=np.asarray(points,dtype=float)
    if points.ndim!=2 or points.shape[1]!=4 or not np.all(np.isfinite(points)):
        raise ValueError('points require a finite (N,4) array')
    if np.any(abs(np.sum(points*points,axis=1)-1)>2e-10):
        raise ValueError('Wigner harmonics require unit-sphere points')
    labels=_validated_labels(labels)
    values=np.empty((len(points),len(labels)),dtype=float)
    if not labels or not len(points):return values
    selected={}
    for column,label in enumerate(labels):selected.setdefault(label[0],[]).append((column,label))
    x,y,z,w=points.T
    a=w-1j*z;b=-y-1j*x;c=y-1j*x;d=w+1j*z
    previous=np.ones((len(points),1,1),dtype=complex)
    for l in range(max(selected)+1):
        if l:
            current=np.zeros((len(points),l+1,l+1),dtype=complex)
            for r in range(l+1):
                for s in range(l+1):
                    entry=current[:,r,s]
                    if r<l and s<l:entry+=math.sqrt((l-r)*(l-s))/l*a*previous[:,r,s]
                    if r<l and s:entry+=math.sqrt((l-r)*s)/l*b*previous[:,r,s-1]
                    if r and s<l:entry+=math.sqrt(r*(l-s))/l*c*previous[:,r-1,s]
                    if r and s:entry+=math.sqrt(r*s)/l*d*previous[:,r-1,s-1]
        else:current=previous
        normalization=math.sqrt((l+1)/VOLUME)
        for column,(_,r,s,kind) in selected.get(l,()):
            value=current[:,r,s]
            if kind=='cos':values[:,column]=math.sqrt(2)*normalization*value.real
            elif kind=='sin':values[:,column]=math.sqrt(2)*normalization*value.imag
            else:values[:,column]=normalization*value.real
        previous=current
    return values


@lru_cache(maxsize=16)
def _point_label_arrays(labels):
    labels=_validated_labels(labels)
    kind={'cos':0,'sin':1,'real':2}
    return tuple(np.asarray([label[i] if i<3 else kind[label[3]] for label in labels],dtype=np.int64)
                 for i in range(4))


@njit(nogil=True)
def _point_recurrence(point,ll,rr,ss,kind):
    """Point-local recurrence: no Python point loops and no disk JIT cache."""
    L=int(np.max(ll));values=np.empty(len(ll),dtype=np.float64)
    previous=np.zeros((L+1,L+1),dtype=np.complex128)
    current=np.zeros_like(previous);previous[0,0]=1.
    norm=math.sqrt(np.sum(point*point))
    x,y,z,w=point[0]/norm,point[1]/norm,point[2]/norm,point[3]/norm
    a=w-1j*z;b=-y-1j*x;c=y-1j*x;d=w+1j*z
    for l in range(L+1):
        if l:
            for r in range(l+1):
                for s in range(l+1):
                    value=0j
                    if r<l and s<l:value+=math.sqrt((l-r)*(l-s))/l*a*previous[r,s]
                    if r<l and s:value+=math.sqrt((l-r)*s)/l*b*previous[r,s-1]
                    if r and s<l:value+=math.sqrt(r*(l-s))/l*c*previous[r-1,s]
                    if r and s:value+=math.sqrt(r*s)/l*d*previous[r-1,s-1]
                    current[r,s]=value
        else:current[0,0]=1.
        normalization=math.sqrt((l+1)/VOLUME)
        for k in range(len(ll)):
            if ll[k]!=l:continue
            value=current[rr[k],ss[k]]
            if kind[k]==0:values[k]=math.sqrt(2)*normalization*value.real
            elif kind[k]==1:values[k]=math.sqrt(2)*normalization*value.imag
            else:values[k]=normalization*value.real
        previous,current=current,previous
    return values


def point_values(point,labels):
    """Fast one-observer row, normalizing the quaternion at each RK stage.

    This defines a smooth directional extension for intermediate RK stages;
    it does not renormalize or otherwise modify the driver's stored state.
    """
    point=np.asarray(point,dtype=float)
    if point.shape!=(4,) or not np.all(np.isfinite(point)) or np.dot(point,point)<=0:
        raise ValueError('a finite nonzero quaternion point is required')
    labels=tuple(labels)
    if not labels:return np.empty(0,dtype=float)
    return _point_recurrence(point,*_point_label_arrays(labels))


@dataclass(frozen=True,slots=True)
class FrameBasis:
    degree: int
    labels: tuple
    characters: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    values: np.ndarray
    generators: tuple

    def _coefficients(self,c):
        c=np.asarray(c,dtype=float)
        if c.ndim not in (1,2) or c.shape[0]!=len(self.labels):
            raise ValueError('coefficients require shape (M,) or (M,K)')
        return c

    def evaluate(self,c):
        return self.values@self._coefficients(c)

    def project(self,f):
        f=np.asarray(f,dtype=float)
        if f.ndim not in (1,2) or f.shape[0]!=len(self.weights):
            raise ValueError('samples require shape (N,) or (N,K)')
        weighted=self.weights*f if f.ndim==1 else self.weights[:,None]*f
        return self.values.T@weighted


def build_frame_basis(degree,quadrature_order=None):
    """Build full real modes and sparse coefficient derivative generators.

    Default quadrature uses n_u=n_v=L+2,n_phi=2L+3. An explicit integer n
    uses n_u=n_v=n,n_phi=2n+1 for caller-controlled dealiasing. G+ characters
    are columns (g1 flip X1/X2, g2 flip X1/X4). G+ alone retains X3 at l=1;
    it is not the earlier extra-reflection, kernel-free constraint subspace.
    """
    degree=_degree(degree)
    if quadrature_order is None:n=degree+2;nphi=2*degree+3
    else:
        n=_degree(quadrature_order)
        if n<1:raise ValueError('quadrature_order must be positive')
        nphi=2*n+1
    levels=[_level(l) for l in range(degree+1)]
    labels=tuple(label for level in levels for label in level[0])
    characters=np.concatenate([level[1] for level in levels],axis=0)
    generators=tuple(sparse.block_diag([level[2][i] for level in levels],format='csr') for i in range(3))
    for D in generators:
        D.sum_duplicates();D.eliminate_zeros();D.sort_indices()
        for array in (D.data,D.indices,D.indptr):array.setflags(write=False)
    points,weights=sphere_quadrature(n,n,nphi)
    values=arbitrary_points_values(points,labels)
    for array in (characters,points,weights,values):array.setflags(write=False)
    return FrameBasis(degree,labels,characters,points,weights,values,generators)
