"""Longitudinal York completion of the inherited scalar momentum kicks.

All tensors use the physical round-reference orthonormal frame E_i=r e_i.
For h=psi^4*hbar, K_ij=psi^-2*Abar_ij and Pi=lambda*f/psi^6, the conformal
factor cancels from the momentum equation. This module solves its finite
polynomial reference equation. It supplies no freely specifiable TT tensor.
"""
from dataclasses import dataclass
import math

import numpy as np

from .nsc_atlas_finite_parent_data import source_profile
from .nsc_atlas_parent_circulation import Parent
from .nsc_atlas_s3_frame import arbitrary_points_values,build_frame_basis


def _scale(r):
    if isinstance(r,(bool,np.bool_)):raise ValueError('positive finite reference scale required')
    r=float(r)
    if not math.isfinite(r) or r<=0 or not math.isfinite(r*r) or r*r==0:
        raise ValueError('positive finite reference scale required')
    return r


def longitudinal_operator(basis,*,r=1.):
    """Return div(L W), with vector coefficients ordered (W1,W2,W3).

    L_ij=E_i W_j+E_j W_i-2/3 delta_ij E_k W_k. The symmetric tensor's
    connection contraction vanishes in this regular frame, leaving
    div(L)_i=E_j L_ji. Derivative order Dj@Di is essential.
    """
    r=_scale(r);D=[value.toarray() for value in basis.generators]
    lap=sum(value@value for value in D);n=len(basis.labels)
    blocks=[[r*r*((lap if i==j else np.zeros((n,n)))+D[j]@D[i]-(2/3)*D[i]@D[j])
             for j in range(3)] for i in range(3)]
    result=np.block(blocks)
    scale=float(np.max(abs(result)))
    if np.max(abs(result-result.T))>5e-12*scale:
        raise ArithmeticError('York operator lost self-adjointness')
    return (result+result.T)/2


def longitudinal_tensor(basis,vector_coefficients,*,r=1.):
    """Return symmetric tracefree Abar coefficients with shape(3,3,modes)."""
    r=_scale(r);W=np.asarray(vector_coefficients,dtype=float);n=len(basis.labels)
    if W.shape!=(3,n) or not np.all(np.isfinite(W)):
        raise ValueError('finite vector coefficients with shape(3,modes) required')
    derivative=np.array([[(r*D)@W[j] for j in range(3)] for D in basis.generators])
    divergence=sum(derivative[i,i] for i in range(3))
    A=derivative+derivative.swapaxes(0,1)
    for i in range(3):A[i,i]-=(2/3)*divergence
    return A


def tensor_divergence(basis,tensor_coefficients,*,r=1.):
    r=_scale(r);A=np.asarray(tensor_coefficients,dtype=float);n=len(basis.labels)
    if A.shape!=(3,3,n) or not np.all(np.isfinite(A)):
        raise ValueError('finite tensor coefficients with shape(3,3,modes) required')
    return np.array([sum(r*(basis.generators[j]@A[j,i]) for j in range(3)) for i in range(3)])


def inherited_polynomial_inverse(basis,source_coefficients,*,r=1.):
    """Closed inverse on the inherited source's four spectral sectors.

    For Q=div(L)/r^2 the two polynomial sources obey
    (Q+5)(Q+16)(Q+20/3)(Q+28)j=0. Thus
    W=-(3Q^3+167Q^2+2984Q+20080)j/(44800*r^2).
    This is a source-specific identity, not a pseudoinverse for arbitrary j.
    """
    r=_scale(r);source=np.asarray(source_coefficients,dtype=float);n=len(basis.labels)
    if source.shape!=(3,n) or not np.all(np.isfinite(source)):
        raise ValueError('finite source coefficients with shape(3,modes) required')
    Q=longitudinal_operator(basis,r=1.)
    f=source.ravel();q1=Q@f;q2=Q@q1;q3=Q@q2
    annihilator=Q@q3+(167/3)*q3+(2984/3)*q2+(20080/3)*q1+(44800/3)*f
    norm=float(np.linalg.norm(f));relative=float(np.linalg.norm(annihilator))/norm if norm else 0.
    operator_norm=float(np.max(np.sum(abs(Q),axis=1)))
    polynomial_bound=operator_norm**4+(167/3)*operator_norm**3+(2984/3)*operator_norm**2+(20080/3)*operator_norm+44800/3
    # The fourth-degree certificate magnifies harmless projection roundoff.
    # Its scaled backward error, not its unscaled cancellation residual,
    # tests whether a supplied vector belongs to the stated sectors.
    if relative>2e-12*polynomial_bound:
        raise ValueError('source is outside the inherited polynomial inverse sectors')
    vector=-(3*q3+167*q2+2984*q1+20080*f)/(44800*r*r)
    return vector.reshape(3,n),relative


def solve_longitudinal(basis,source_coefficients,*,r=1.,compatibility_tolerance=2e-11):
    """Solve only after checking all conformal-Killing compatibility charges.

    The ten zero modes are gauge freedom of W, not removed physical tensor
    modes. A=L W is unique in the longitudinal range. Every additional
    round-reference divergence-free tracefree tensor increases its
    round-reference L2 norm orthogonally. This does not minimize the
    conformally weighted physical norm or assert absence of radiation.
    """
    r=_scale(r);source=np.asarray(source_coefficients,dtype=float);n=len(basis.labels)
    if basis.degree<2:raise ValueError('degree>=2 is required to contain all ten conformal Killing fields')
    if source.shape!=(3,n) or not np.all(np.isfinite(source)):
        raise ValueError('finite source coefficients with shape(3,modes) required')
    if not math.isfinite(compatibility_tolerance) or compatibility_tolerance<=0:
        raise ValueError('positive finite compatibility tolerance required')
    operator=longitudinal_operator(basis,r=r)
    eigenvalues,eigenvectors=np.linalg.eigh(operator)
    scale=float(np.max(abs(eigenvalues)));zero=np.abs(eigenvalues)<=2e-11*scale
    if int(np.sum(zero))!=10 or float(np.max(eigenvalues))>2e-11*scale:
        raise ArithmeticError('unexpected conformal-Killing kernel or positive York spectrum')
    f=source.ravel();charges=eigenvectors[:,zero].T@f
    norm=float(np.linalg.norm(f));charge_norm=float(np.linalg.norm(charges))
    if charge_norm>compatibility_tolerance*max(norm,np.finfo(float).tiny):
        raise ValueError('source violates conformal-Killing compatibility')
    transformed=eigenvectors.T@f;amplitudes=np.zeros_like(transformed)
    amplitudes[~zero]=transformed[~zero]/eigenvalues[~zero]
    vector=(eigenvectors@amplitudes).reshape(3,n)
    tensor=longitudinal_tensor(basis,vector,r=r)
    residual=tensor_divergence(basis,tensor,r=r)-source
    tensor_norm2=float(np.sum(tensor*tensor))
    identity=float(np.sum(vector*source)+tensor_norm2/2)
    return dict(vector_coefficients=vector,tensor_coefficients=tensor,
        source_coefficients=source,kernel_eigenvalues=eigenvalues[zero],
        diagnostics=dict(kernel_dimension=int(np.sum(zero)),
            source_coefficient_norm=norm,CKV_charge_norm=charge_norm,
            CKV_charge_over_source_norm=charge_norm/norm if norm else 0.,
            divergence_residual_coefficient_norm=float(np.linalg.norm(residual)),
            trace_coefficient_norm=float(np.linalg.norm(np.trace(tensor,axis1=0,axis2=1))),
            tensor_reference_L2_norm_squared=tensor_norm2/r**3,
            vector_reference_L2_norm_squared=float(np.sum(vector*vector))/r**3,
            integration_by_parts_residual=identity/r**3,
            minimum_nonzero_absolute_eigenvalue=float(np.min(abs(eigenvalues[~zero]))),
            TT_choice='Zero freely specified conformal TT datum; Abar has minimal round-reference L2 norm among tracefree tensors with this divergence. No physical radiation-state selection is claimed.'))


@dataclass
class YorkKick:
    parent: object
    epsilon: float
    source_harmonic: int
    basis: object
    vector_coefficients: np.ndarray
    tensor_coefficients: np.ndarray
    source_coefficients: np.ndarray
    diagnostics: dict

    def values(self,points):
        Y=arbitrary_points_values(points,self.basis.labels)
        return (Y@self.tensor_coefficients.reshape(9,-1).T).reshape(len(Y),3,3)

    def derivatives(self,points):
        """Raw physical E_k component derivatives, shape(N,k,i,j)."""
        Y=arbitrary_points_values(points,self.basis.labels)
        flat=self.tensor_coefficients.reshape(9,-1).T
        return np.stack([(Y@(self.parent.H*(D@flat))).reshape(len(Y),3,3)
                         for D in self.basis.generators],axis=1)

    def squared_norm(self,points):
        A=self.values(points)
        return np.sum(A*A,axis=(1,2))

    def momentum_profile(self,points):
        x=np.asarray(points,dtype=float)
        if x.ndim!=2 or x.shape[1]!=4 or not np.all(np.isfinite(x)):
            raise ValueError('finite unit-sphere points with shape(N,4) required')
        if np.any(abs(np.sum(x*x,axis=1)-1)>2e-10):raise ValueError('unit-sphere points required')
        profile=math.sqrt(2)/math.pi*x[:,0] if self.source_harmonic==1 else math.sqrt(12)/math.pi*x[:,1]*x[:,3]
        return self.parent.H**1.5*profile

    def source_values(self,points):
        return arbitrary_points_values(points,self.basis.labels)@self.source_coefficients.T


def solve_york_kick(*,epsilon=.25,source_harmonic=1,degree=4,parent=None):
    """Complete one fixed normalized scalar momentum kick, with no fitted data."""
    if isinstance(source_harmonic,(bool,np.bool_)) or not isinstance(source_harmonic,(int,np.integer)) or source_harmonic not in (1,2):
        raise ValueError('source_harmonic must be1 or2')
    if isinstance(degree,(bool,np.bool_)) or not isinstance(degree,(int,np.integer)) or not 4<=degree<=8:
        raise ValueError('this polynomial completion requires integer4<=degree<=8')
    if isinstance(epsilon,(bool,np.bool_)) or not math.isfinite(epsilon) or not 0<=epsilon<=1:
        raise ValueError('finite0<=epsilon<=1 required')
    parent=parent or Parent.retained();r=parent.H
    basis=build_frame_basis(int(degree),int(degree)+3)
    profile,_=source_profile(basis.points,parent)
    chi_coefficients=basis.project(float(epsilon)*profile)
    gradient=np.stack([basis.evaluate(r*(D@chi_coefficients)) for D in basis.generators],axis=1)
    points=basis.points
    f=(math.sqrt(2)/math.pi*points[:,0] if source_harmonic==1 else math.sqrt(12)/math.pi*points[:,1]*points[:,3])*r**1.5
    source=-parent.kappa*f[:,None]*gradient
    coefficients=basis.project(source).T
    result=solve_longitudinal(basis,coefficients,r=r)
    closed,annihilator=inherited_polynomial_inverse(basis,coefficients,r=r)
    projection=basis.evaluate(coefficients.T)-source
    result['diagnostics'].update(
        source_projection_L2_residual=float(np.sqrt(np.sum(basis.weights[:,None]*projection*projection))),
        closed_form='W=-(3Q^3+167Q^2+2984Q+20080)j/(44800*r^2), Q=div(L)/r^2',
        source_annihilator='(Q+5)(Q+16)(Q+20/3)(Q+28)j=0',
        source_annihilator_relative_roundoff=annihilator,
        closed_inverse_relative_difference=float(np.linalg.norm(closed-result['vector_coefficients']))/max(float(np.linalg.norm(result['vector_coefficients'])),np.finfo(float).tiny),
        source_polynomial_maximum_degree=4,
        conformal_identity='D_j K^ij=psi^-10 Dbar_j Abar^ij and j^i=-psi^-10 f Dbar^i chi; hence divbar(Abar)=-kappa f gradbar(chi).',
        degree=int(degree),source_harmonic=int(source_harmonic),epsilon=float(epsilon),
        physical_reference_scale=r,freely_specified_TT=False)
    return YorkKick(parent,float(epsilon),int(source_harmonic),basis,
        result['vector_coefficients'],result['tensor_coefficients'],coefficients,result['diagnostics'])
