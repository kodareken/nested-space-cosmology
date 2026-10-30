"""Exact closure and independent tensor identities for the cubic source maps."""
import numpy as np
import sympy as s

from recursive_horizons.nsc_atlas_s3_polynomials import (
    X,R2,laplacian,gradient_dot,reduce_on_sphere,inner_product)
from recursive_horizons.nsc_atlas_cubic_operators import (
    scalar_basis,cubic_basis,exact_operator_matrices,numeric_operator_matrices)

P=s.eye(4)-X*X.T


def grad(f):return P*s.Matrix([s.diff(f,x) for x in X])


def hess(f):
    return P*s.hessian(f,X)*P-sum(x*s.diff(f,x) for x in X)*P


def div(tensor):
    raw=s.Matrix([sum(P[a,b]*s.diff(tensor[b,j],X[a]) for a in range(4) for b in range(4)) for j in range(4)])
    return P*raw


def delta_laplacian(tensor,f):
    trace=s.trace(tensor)
    Hess=hess(f)
    return -sum(tensor[i,j]*Hess[i,j] for i in range(4) for j in range(4))-(div(tensor)-grad(trace)/2).dot(grad(f))


def polynomial(coefficients):return sum(c*f for c,f in zip(coefficients,cubic_basis()))


def test_laplacian_columns_close_exactly_and_have_the_required_spectrum():
    D=exact_operator_matrices()['laplacian']
    for j,F in enumerate(cubic_basis()):
        assert reduce_on_sphere(polynomial(D[:,j])-laplacian(F))==0
    assert D.eigenvals()=={-3:1,-8:1,-15:2,-24:2,-35:1,-48:1}
    gram=s.Matrix([[inner_product(a,b) for b in cubic_basis()] for a in cubic_basis()])
    assert gram*D==D.T*gram


def test_products_and_metric_traces_reconstruct_exactly():
    exact=exact_operator_matrices();S=scalar_basis();T=exact['trace']
    Y=(s.sqrt(2)*X[0]/s.pi,s.sqrt(12)*X[1]*X[3]/s.pi)
    for a,F in enumerate(S):
        for j,Yj in enumerate(Y):
            coefficients=[exact['multiplication'][k,a,j] for k in range(8)]
            assert reduce_on_sphere(polynomial(coefficients)-F*Yj)==0
        assert sum(T[k,a]*S[k] for k in range(5))==3*F
    for a,F in enumerate(S[1:],start=5):
        assert reduce_on_sphere(sum(T[k,a]*S[k] for k in range(5))-laplacian(F))==0
    assert T[:,9]==s.zeros(5,1) and T[:,10]==s.zeros(5,1)


def test_conformal_metric_map_matches_full_covariant_tensor_variation():
    exact=exact_operator_matrices();F=scalar_basis()[1]
    Y=s.sqrt(12)*X[1]*X[3]/s.pi
    full=delta_laplacian(P*F,Y)
    expected=-F*laplacian(Y)+gradient_dot(F,Y)/2
    mapped=polynomial([exact['metric_delta_laplacian'][k,1,1] for k in range(8)])
    assert s.simplify(reduce_on_sphere(full-expected))==0
    assert s.simplify(reduce_on_sphere(full-mapped))==0


def test_hessian_coordinate_pullback_map_matches_tensor_variation():
    exact=exact_operator_matrices();F=scalar_basis()[2]
    Y=s.sqrt(12)*X[1]*X[3]/s.pi
    full=delta_laplacian(hess(F),Y)
    mapped=polynomial([exact['metric_delta_laplacian'][k,6,1] for k in range(8)])
    assert s.simplify(reduce_on_sphere(full-mapped))==0


def test_vector_coordinate_pullback_and_tt_normalization():
    exact=exact_operator_matrices();z,x,_,y=X
    C1=s.sqrt(2)/s.pi;C2=s.sqrt(12)/s.pi
    W=2*s.sqrt(6)/(3*s.pi**2)*s.Matrix([-4*x*y,2*z*y,0,2*z*x])
    vector=P*(W.jacobian(X)+W.jacobian(X).T)*P
    mapped=polynomial([exact['metric_delta_laplacian'][k,9,0] for k in range(8)])
    assert s.simplify(reduce_on_sphere(delta_laplacian(vector,C1*z)-mapped))==0
    S2,S4=scalar_basis()[2],scalar_basis()[4]
    Y2=C2*x*y
    Z=grad(Y2)*grad(Y2).T-C2*C2*(P/9+P*S2/2+hess(S2)/8+4*P*S4/3+hess(S4)/3)
    assert s.simplify(reduce_on_sphere(s.trace(Z)))==0
    assert all(s.simplify(reduce_on_sphere(value))==0 for value in div(Z))
    contraction=-sum(Z[i,j]*hess(Y2)[i,j] for i in range(4) for j in range(4))
    expected=2*C2*C2*Y2/9
    assert s.simplify(reduce_on_sphere(contraction-expected))==0
    tt_mapped=polynomial([exact['metric_delta_laplacian'][k,10,1] for k in range(8)])
    assert s.simplify(reduce_on_sphere(tt_mapped-expected))==0
    assert all(exact['metric_delta_laplacian'][k,10,0]==0 for k in range(8))


def test_every_metric_image_is_in_the_declared_finite_polynomial_space():
    exact=exact_operator_matrices();S=scalar_basis();Y=(s.sqrt(2)*X[0]/s.pi,s.sqrt(12)*X[1]*X[3]/s.pi)
    for a,F in enumerate(S):
        for j,Yj in enumerate(Y):
            conformal=-F*laplacian(Yj)+gradient_dot(F,Yj)/2
            actual=polynomial([exact['metric_delta_laplacian'][k,a,j] for k in range(8)])
            assert reduce_on_sphere(conformal-actual)==0
    for a,F in enumerate(S[1:],start=5):
        for j,Yj in enumerate(Y):
            commutator=(gradient_dot(F,laplacian(Yj))-laplacian(gradient_dot(F,Yj)))/2
            actual=polynomial([exact['metric_delta_laplacian'][k,a,j] for k in range(8)])
            assert reduce_on_sphere(commutator-actual)==0


def test_float_interfaces_preserve_exact_normalization_and_index_order():
    exact=exact_operator_matrices();numeric=numeric_operator_matrices()
    for key,shape in (('laplacian',(8,8)),('metric_delta_laplacian',(8,11,2)),
                      ('multiplication',(8,5,2)),('trace',(5,11))):
        assert numeric[key].shape==shape
        assert numeric[key].dtype==np.float64
        assert not numeric[key].flags.writeable
        np.testing.assert_allclose(numeric[key],np.array(exact[key].tolist(),float),rtol=0,atol=0)
    assert numeric['metric_delta_laplacian'][0,0,0]==float(3*s.sqrt(2)/s.pi)
    assert numeric['metric_delta_laplacian'][1,0,1]==float(8*s.sqrt(12)/s.pi)
    assert len(numeric['labels']['metric'])==11
