"""Independent small-tree, closed-form and spectral tests of the theorem."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_infinite_inheritance import InheritedTree


def test_matrix_recursion_matches_noncommuting_full_tree_and_subtree_scale():
    H=np.array([[4.,.2j],[-.2j,3.]])
    B=np.array([[.3,.1j],[.07,.2]])
    tree=InheritedTree(H,B,3.,2)
    D,levels=tree.finite_matrix(2)
    for lam in (0.,.7,3.):
        direct=np.linalg.inv(D+lam*np.eye(len(D)))[:2,:2]
        np.testing.assert_allclose(tree.response(lam,2)['response'],direct,atol=2e-15)
    # Breadth-first indices: first child1 has children3,4.
    indices=[2*v+i for v in (1,3,4) for i in (0,1)]
    smaller,_=tree.finite_matrix(1)
    np.testing.assert_allclose(D[np.ix_(indices,indices)],tree.omega*smaller,rtol=2e-15,atol=2e-15)
    eigenvalues=np.linalg.eigvalsh(D)
    np.testing.assert_array_less(tree.lower_form_constant*np.sort(levels)-1e-12,eigenvalues)
    np.testing.assert_array_less(eigenvalues,tree.upper_form_constant*np.sort(levels)+1e-12)
    # Losing the branch/scale factor or reversing a noncommuting link
    # changes the response, rather than passing through a scalar control.
    child=tree.response(.7/tree.omega,1)['response']
    expected=tree.response(.7,2)['response']
    wrong_scale=np.linalg.inv(H+.7*np.eye(2)-B@child@B.conj().T/tree.omega)
    wrong_order=np.linalg.inv(H+.7*np.eye(2)-tree.branching/tree.omega*B.conj().T@child@B)
    assert np.linalg.norm(expected-wrong_scale)>1e-5
    assert np.linalg.norm(expected-wrong_order)>1e-5


def test_scalar_infinite_response_and_uniform_depth_bound():
    tree=InheritedTree(np.array([[3.]]),np.array([[1.]]),4.,3)
    exact=2/(3+math.sqrt(6))
    previous=0.
    for depth in (0,1,3,6):
        result=tree.response(0.,depth);value=result['response'][0,0].real
        assert previous<value<exact
        assert exact-value<=result['analytic_depth_bound']+2e-16
        previous=value
    for lam in (.1,1.,100.):
        a=tree.response(lam,3);b=tree.response(lam,10)
        assert np.linalg.norm(a['response']-b['response'])<=a['analytic_depth_bound']+2e-16
    assert tree.counting_exponent==pytest.approx(math.log(3)/math.log(4))
    assert 2*tree.counting_exponent==pytest.approx(math.log(3)/math.log(2))


def test_uncoupled_limit_and_dense_resource_boundary():
    tree=InheritedTree(np.diag([2.,3.]),np.zeros((2,2)),3.,2)
    for depth in (0,10,100):
        result=tree.response(.2,depth)
        np.testing.assert_allclose(result['response'],np.diag([1/2.2,1/3.2]),atol=0.)
        assert result['analytic_depth_bound']==0.
    with pytest.raises(ValueError,match='256'):
        tree.finite_matrix(8)


def test_domain_assumptions_are_enforced_without_changing_the_inputs():
    for H,B,omega,branch in (([[1]],[[1]],2.,2),([[2]],[[.1]],1.,2),
            ([[2]],[[.1]],2.,1),([[1j]],[[.1]],2.,2)):
        with pytest.raises(ValueError):
            InheritedTree(np.array(H),np.array(B),omega,branch)
