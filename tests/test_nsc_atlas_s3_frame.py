"""Low-degree representation, orientation and full-sphere frame controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_s3_frame import build_frame_basis,arbitrary_points_values,point_values


@pytest.fixture(scope='module')
def frame():return build_frame_basis(4)


def right_rotate(points,axis,t):
    x,y,z,w=points.T;c,s=math.cos(t),math.sin(t)
    if axis==0:return np.c_[c*x+s*w,c*y+s*z,c*z-s*y,c*w-s*x]
    if axis==1:return np.c_[c*x-s*z,c*y+s*w,c*z+s*x,c*w-s*y]
    return np.c_[c*x+s*y,c*y-s*x,c*z+s*w,c*w-s*z]


def test_full_l2_and_l4_real_harmonics_are_orthonormal():
    for L,count in ((2,14),(4,55)):
        b=build_frame_basis(L)
        assert len(b.labels)==count
        G=b.values.T@(b.weights[:,None]*b.values)
        np.testing.assert_allclose(G,np.eye(count),atol=1e-13,rtol=0)
        assert sum(b.weights)==pytest.approx(2*math.pi**2,rel=1e-14)


def test_first_degree_real_modes_are_the_expected_ambient_coordinates(frame):
    b=frame;indices=[i for i,label in enumerate(b.labels) if label[0]==1]
    x,y,z,w=b.points.T
    expected=math.sqrt(2)/math.pi*np.c_[w,-z,-y,-x]
    np.testing.assert_allclose(b.values[:,indices],expected,atol=3e-16,rtol=1e-14)
    np.testing.assert_allclose(b.values[:,0],1/math.sqrt(2*math.pi**2),atol=2e-16)
    # Directly differentiate [w,-z,-y,-x] along the stated quaternion fields.
    first_degree_generators=(
        [[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]],
        [[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]],
        [[0,-1,0,0],[1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    for D,known in zip(b.generators,first_degree_generators):
        np.testing.assert_allclose(D[np.ix_(indices,indices)].toarray(),known,atol=5e-16,rtol=0)


def test_sparse_generators_are_antisymmetric_with_positive_bracket_and_casimir(frame):
    b=frame;D=b.generators;M=len(b.labels)
    assert all(x.format=='csr' for x in D)
    for matrix in D:np.testing.assert_allclose((matrix+matrix.T).toarray(),0,atol=5e-15,rtol=0)
    for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
        np.testing.assert_allclose((D[i]@D[j]-D[j]@D[i]-2*D[k]).toarray(),0,atol=1.8e-14,rtol=0)
    casimir=sum((matrix@matrix for matrix in D))
    target=-np.array([l*(l+2) for l,_,_,_ in b.labels],dtype=float)
    np.testing.assert_allclose(casimir.toarray(),np.diag(target),atol=3e-14,rtol=0)
    assert all(x.nnz<5*M for x in D)


def test_derivative_matrices_act_on_coefficient_columns_without_transpose(frame):
    b=frame;points=b.points[::19];Y=arbitrary_points_values(points,b.labels)
    h=2e-4
    for i,D in enumerate(b.generators):
        plus=arbitrary_points_values(right_rotate(points,i,h),b.labels)
        minus=arbitrary_points_values(right_rotate(points,i,-h),b.labels)
        plus2=arbitrary_points_values(right_rotate(points,i,2*h),b.labels)
        minus2=arbitrary_points_values(right_rotate(points,i,-2*h),b.labels)
        derivative=(-plus2+8*plus-8*minus+minus2)/(12*h)
        np.testing.assert_allclose(derivative,Y@D.toarray(),atol=1.8e-10,rtol=0)


def test_global_flip_characters_are_diagonal_and_derivatives_change_sector(frame):
    b=frame;points=b.points[::17];Y=arbitrary_points_values(points,b.labels)
    for g,signs in enumerate(((-1,-1,1,1),(-1,1,1,-1))):
        moved=arbitrary_points_values(points*np.array(signs),b.labels)
        np.testing.assert_allclose(moved,Y*b.characters[None,:,g],atol=2e-15,rtol=1e-13)
    for D,signs in zip(b.generators,((-1,1),(-1,-1),(1,-1))):
        coo=D.tocoo();nonzero=abs(coo.data)>1e-13
        np.testing.assert_array_equal(b.characters[coo.row[nonzero]],b.characters[coo.col[nonzero]]*np.array(signs))
    # The proper group differs from the older extra-X3-reflection subspace.
    z=b.labels.index((1,0,0,'sin'))
    np.testing.assert_array_equal(b.characters[z],[1,1])


def test_addition_theorem_and_arbitrary_label_subsets(frame):
    b=frame
    subset=b.labels[::3]
    np.testing.assert_allclose(arbitrary_points_values(b.points[::13],subset),b.values[::13,::3],atol=2e-16,rtol=1e-14)
    for l in range(b.degree+1):
        select=np.array([label[0]==l for label in b.labels])
        expected=(l+1)**2/(2*math.pi**2)
        np.testing.assert_allclose(np.sum(b.values[:,select]**2,axis=1),expected,atol=1.3e-14,rtol=0)


def test_evaluation_projection_and_dealiased_quadrature_interface(frame):
    b=frame;c=np.arange(1,len(b.labels)+1,dtype=float)/100
    C=np.c_[c,-c/2]
    np.testing.assert_allclose(b.project(b.evaluate(c)),c,atol=1e-13,rtol=0)
    np.testing.assert_allclose(b.project(b.evaluate(C)),C,atol=1e-13,rtol=0)
    alias=build_frame_basis(2,quadrature_order=6)
    assert len(alias.points)==6*6*13
    assert not b.values.flags.writeable and not b.characters.flags.writeable
    with pytest.raises(ValueError,match='coefficients'):b.evaluate(np.zeros(2))
    with pytest.raises(ValueError,match='unit'):arbitrary_points_values(np.array([[0,0,0,2.]]),b.labels)


def test_fast_observer_row_matches_arbitrary_evaluation_and_rk_direction(frame):
    b=frame
    for point in b.points[::89]:
        expected=arbitrary_points_values(point[None,:],b.labels)[0]
        np.testing.assert_allclose(point_values(point,b.labels),expected,atol=1e-15,rtol=1e-14)
        np.testing.assert_allclose(point_values(1.00001*point,b.labels),expected,atol=2e-15,rtol=1e-14)
    subset=b.labels[::4];point=b.points[5]
    np.testing.assert_allclose(point_values(point,subset),arbitrary_points_values(point[None,:],subset)[0],atol=1e-15)
    with pytest.raises(ValueError,match='nonzero'):point_values(np.zeros(4),b.labels)
