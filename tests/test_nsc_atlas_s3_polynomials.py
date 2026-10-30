"""Exact quotient, harmonic and integral controls for the scalar S3 adapter."""
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_s3_polynomials import (
    X, R2, ambient_laplacian, laplacian, gradient_dot, reduce_on_sphere,
    harmonic_projection_homogeneous, spherical_components, mean, inner_product)


def test_known_harmonic_projections_through_degree_six():
    x1,x2,_,x4=X
    H2=x1*x1-R2/4
    H4=x2*x2*x4*x4-R2*(x2*x2+x4*x4)/8+R2*R2/48
    H5=x1*(x2*x2*x4*x4-R2*(x2*x2+x4*x4)/10+R2*R2/80)
    H6=x2*x4*(x2*x2*x4*x4-R2*(x2*x2+x4*x4)/4+3*R2*R2/40)
    for polynomial,degree,expected in ((x1*x1,2,H2),(x2*x2*x4*x4,4,H4),
                                       (x1*x2*x2*x4*x4,5,H5),(x2**3*x4**3,6,H6)):
        H=harmonic_projection_homogeneous(polynomial,degree)
        assert s.expand(H-expected)==0
        assert ambient_laplacian(H)==0
        assert s.expand(laplacian(H)+degree*(degree+2)*H)==0


def test_sphere_decomposition_reconstructs_degree_eight_with_symbolic_coefficients():
    a,b,c=s.symbols('a b c',real=True)
    polynomial=a*X[0]**2*X[1]**2*X[2]**2*X[3]**2+b*X[0]*X[1]**2*X[3]**2+c*R2*X[2]
    components=spherical_components(polynomial)
    assert list(components)==sorted(components)
    assert max(components)==8
    assert reduce_on_sphere(sum(components.values())-polynomial)==0
    for degree,H in components.items():
        assert ambient_laplacian(H)==0
        assert s.expand(laplacian(H)+degree*(degree+2)*H)==0


def test_radial_extensions_have_the_same_sphere_harmonics():
    f=3*X[0]*X[1]+X[0]**3*X[3]-2
    base=spherical_components(f)
    radial=spherical_components(R2*f)
    assert base.keys()==radial.keys()
    for degree in base:assert s.expand(base[degree]-radial[degree])==0
    assert spherical_components(0)=={}
    assert spherical_components(7)=={0:s.Integer(7)}


def test_quotient_reduction_preserves_intrinsic_differentiation():
    f=X[0]**3+2*X[1]*X[2]
    extension=f+(R2-1)*(X[0]**2*X[3]+X[1])
    assert reduce_on_sphere(extension)==reduce_on_sphere(f)
    assert reduce_on_sphere(laplacian(extension)-laplacian(f))==0
    assert reduce_on_sphere(gradient_dot(extension,X[3])-gradient_dot(f,X[3]))==0
    assert reduce_on_sphere(R2)==1


def test_intrinsic_gradient_pairing_obeys_product_rule():
    f=X[0]*X[1]**2
    g=X[2]*X[3]+X[0]**3
    assert s.expand(laplacian(f*g)-f*laplacian(g)-g*laplacian(f)-2*gradient_dot(f,g))==0
    assert gradient_dot(X[0],X[0])==1-X[0]**2
    assert gradient_dot(X[0],X[1])==-X[0]*X[1]
    assert reduce_on_sphere(gradient_dot(R2,g))==0


def test_normalized_moments_and_harmonic_orthogonality():
    assert mean(1)==1
    assert mean(X[0]**2)==s.Rational(1,4)
    assert mean(X[0]**4)==s.Rational(1,8)
    assert mean(X[0]**2*X[1]**2)==s.Rational(1,24)
    assert mean(X[0]**2*X[1]**2*X[2]**2*X[3]**2)==s.Rational(1,1920)
    assert mean(X[0]*X[1]**2)==0
    assert mean(R2**4)==1
    components=spherical_components(X[0]**4*X[1]**2)
    for degree,H in components.items():
        if degree:assert mean(H)==0
        for other,G in components.items():
            if degree!=other:assert inner_product(H,G)==0


def test_algebraic_and_symbolic_parameters_are_coefficients_not_coordinates():
    a,b,m2=s.symbols('a b m2',real=True)
    f=a*X[0]+s.sqrt(2)*b*X[1]*X[3]+m2
    assert mean(f)==m2
    assert inner_product(f,f)==a*a/4+b*b/12+m2*m2
    parts=spherical_components(f)
    assert parts=={0:m2,1:a*X[0],2:s.sqrt(2)*b*X[1]*X[3]}


def test_invalid_polynomials_and_inconsistent_homogeneous_degrees_are_rejected():
    for expression in (s.sin(X[0]),s.exp(X[1]),1/(1+X[2]**2)):
        with pytest.raises(ValueError,match='polynomial'):mean(expression)
    for degree in (-1,1.5,True):
        with pytest.raises(ValueError,match='degree'):harmonic_projection_homogeneous(X[0],degree)
    with pytest.raises(ValueError,match='homogeneous'):
        harmonic_projection_homogeneous(X[0]+X[1]**2,2)
    assert harmonic_projection_homogeneous(0,8)==0
    # A same-named Symbol with different assumptions must not silently become
    # a coefficient and disappear from the spatial derivative.
    with pytest.raises(ValueError,match='exported coordinate'):
        laplacian(s.Symbol('x1')**2)
