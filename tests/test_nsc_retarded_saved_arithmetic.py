"""Mutants cannot enter the finite saved-diagnostic portability policy."""
import copy

import numpy as np
import pytest

from recursive_horizons.nsc_retarded_saved_arithmetic import (
    require_diagnostic_receipt, require_linear_system, require_matrix_product,
    require_singular_invariants)


def test_matrix_equation_rejects_mutated_solution_and_right_hand_side():
    A=np.array([[1.,.2j],[-.2j,1.1]])
    B=np.array([[.3,.1j],[.2j,.4]])
    X=np.linalg.solve(A,B)
    require_linear_system(A,B,X,X,label='sealed example')
    bad=X.copy();bad[0,0]+=1e-8
    with pytest.raises(ValueError,match='residual failed'):
        require_linear_system(A,B,X,bad,label='sealed example')
    with pytest.raises(ValueError,match='residual failed'):
        require_linear_system(A,B+1e-8,X,X,label='sealed example')


def test_product_check_charges_uncancelled_terms_and_rejects_mutants():
    A=np.array([[1e5,1e5]],complex);B=np.array([[1.],[-1.]],complex)
    value=A@B
    # The zero result has a nonzero arithmetic scale from its actual operands.
    assert require_matrix_product(value,value,A,B,label='cancelled sum')>0
    with pytest.raises(ValueError,match='product replay failed'):
        require_matrix_product(value,value+1e-5,A,B,label='cancelled sum')


def test_spectral_invariants_reject_a_changed_singular_value():
    A=np.array([[1.,.2,.1],[.1,2.,.3],[.2,.1,3.]])
    s=np.linalg.svd(A,compute_uv=False)
    require_singular_invariants(A,s)
    bad=s.copy();bad[0]+=1e-5
    with pytest.raises(ValueError,match='invariant residual failed'):
        require_singular_invariants(A,bad)


def test_only_named_float_diagnostics_can_use_an_operation_allowance():
    saved={'diagnostic':1.,'physical_tolerance':3e-11,'status':'OPEN','count':3,'sha':'original'}
    current=copy.deepcopy(saved);current['diagnostic']+=float(np.finfo(float).eps)
    allowed={('diagnostic',):4*np.finfo(float).eps}
    assert len(require_diagnostic_receipt(current,saved,allowed))==1
    for key,value in [('diagnostic',1.1),('physical_tolerance',4e-11),('status','PASS'),('count',3.),('sha','mutant')]:
        bad=copy.deepcopy(saved);bad[key]=value
        with pytest.raises(ValueError,match='receipt'):
            require_diagnostic_receipt(bad,saved,allowed)
