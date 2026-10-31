"""Independent Clifford/sphere controls of the projected 16D contact vertex.

Expected tensors are evaluated from raw intrinsic Clifford products, not
from a second copy of the coefficient table in the production function.
All matrix comparisons use exact arithmetic.
"""
import itertools
import math

import pytest
import sympy as s
from sympy.physics.matrices import mgamma

from recursive_horizons.nsc_atlas_valve_source_gate import (
    complete_ground_contact_vertex, killing_ground_contact_block,
)


def _product(matrices,dimension):
    out=s.eye(dimension)
    for matrix in matrices:out=out*matrix
    return out


def _clifford(ambient_dimension):
    k=(ambient_dimension-1)//2
    X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1);I=s.eye(2)
    result=[]
    for j in range(k):
        result.extend(s.kronecker_product(*([Z]*j+[P]+[I]*(k-j-1))) for P in (X,Y))
    result.append(s.kronecker_product(*([Z]*k)))
    return result


def _is_zero(matrix):return all(entry==0 for entry in matrix)


@pytest.mark.parametrize('internal_dimension',[2,4])
def test_all_killing_sign_blocks_against_raw_sphere_projection(internal_dimension):
    N=internal_dimension+1;gamma=_clifford(N);dimension=gamma[0].rows
    identity=s.eye(dimension);tensor_zero=s.zeros(dimension**2)
    for i,j in itertools.product(range(N),repeat=2):
        assert _is_zero(gamma[i]*gamma[j]+gamma[j]*gamma[i]
                        -(2*identity if i==j else s.zeros(dimension)))
    # Cq contracts ordered indices. q! here follows directly by enumerating
    # one representative of each antisymmetric Clifford component.
    ambient={rank:math.factorial(rank)*sum((s.kronecker_product(G,G)
        for indices in itertools.combinations(range(N),rank)
        for G in [_product([gamma[i] for i in indices],dimension)]),tensor_zero)
        for rank in range(5)}
    sign_pairs=list(itertools.product((1,-1),repeat=2))
    raw={(r,a,b,c,d):tensor_zero.copy() for r in range(4)
         for a,b in sign_pairs for c,d in sign_pairs}
    # Signed coordinate axes integrate n_I n_J exactly. Antisymmetry
    # removes all higher moments from the contracted pointwise products.
    for axis in range(N):
        for orientation in (1,-1):
            Gn=orientation*gamma[axis];mode_map={1:identity,-1:Gn}
            tangent=[gamma[i] for i in range(N) if i!=axis]
            intrinsic=[s.I*Gn*G for G in tangent]
            for r in range(4):
                for indices in itertools.combinations(range(internal_dimension),r):
                    G=_product([intrinsic[i] for i in indices],dimension)
                    bilinear={(a,b):mode_map[a].H*G*mode_map[b] for a,b in sign_pairs}
                    for (a,b),(c,d) in itertools.product(sign_pairs,repeat=2):
                        raw[r,a,b,c,d]+=math.factorial(r)*s.kronecker_product(
                            bilinear[a,b],bilinear[c,d])/(2*N)
    for r in range(4):
        for signs in itertools.product((1,-1),repeat=4):
            block=killing_ground_contact_block(r,signs,internal_dimension=internal_dimension)
            coefficient=s.Rational(block['ordered_tensor_coefficient'])
            predicted=coefficient*ambient[block['ambient_Clifford_rank']]
            assert _is_zero(raw[(r,)+signs]-predicted),(internal_dimension,r,signs)
            assert block['coefficient']==pytest.approx(float(coefficient),abs=1e-15)


def test_external_three_form_weights_against_full_six_dimensional_clifford():
    g4=[mgamma(i) for i in range(4)];gamma5=mgamma(5);eta=[-1,1,1,1]
    # This controls the overall conversion from the general three-form
    # operator to the already fixed mostly-plus external axial contact.
    J3=sum((6*math.prod(eta[i] for i in indices)*s.kronecker_product(G,G)
        for indices in itertools.combinations(range(4),3)
        for G in [_product([g4[i] for i in indices],4)]),s.zeros(16))
    axial=sum((eta[i]*s.kronecker_product(g4[i]*gamma5,g4[i]*gamma5)
               for i in range(4)),s.zeros(16))
    assert _is_zero(J3-6*axial)
    ambient=_clifford(3);normal=ambient[2]
    intrinsic=[s.I*normal*ambient[i] for i in (0,1)]
    gamma6=[s.kronecker_product(G,s.eye(2)) for G in g4]
    gamma6.extend(s.I*s.kronecker_product(gamma5,G) for G in intrinsic)
    signs=eta+[1,1];direct=s.zeros(64)
    for indices in itertools.combinations(range(6),3):
        G=_product([gamma6[i] for i in indices],8)
        direct+=6*math.prod(signs[i] for i in indices)*s.kronecker_product(G,G)
    decomposed=s.zeros(64)
    for r in range(4):
        p=3-r
        weight=killing_ground_contact_block(r,(1,1,1,1),internal_dimension=2)['external_three_form_weight']
        for external_indices in itertools.combinations(range(4),p):
            for internal_indices in itertools.combinations(range(2),r):
                external=_product([g4[i] for i in external_indices],4)*gamma5**r
                internal=_product([intrinsic[i] for i in internal_indices],2)
                G=s.kronecker_product(external,internal)
                decomposed+=weight*math.factorial(p)*math.factorial(r)*math.prod(
                    eta[i] for i in external_indices)*s.kronecker_product(G,G)
    assert _is_zero(direct-decomposed)


@pytest.mark.parametrize('Killing_sign',[1,-1])
def test_positive_mass_rotation_for_every_internal_index_count(Killing_sign):
    gamma=[mgamma(i) for i in range(4)];gamma5=mgamma(5);identity=s.eye(4)
    U=(identity-s.I*Killing_sign*gamma5)/s.sqrt(2)
    # Rotate the actual pseudoscalar KK mass to a positive scalar mass.
    assert _is_zero(s.simplify(U*(-s.I*Killing_sign*gamma5)*U+identity))
    expected=[_product(gamma[:3],4),Killing_sign*_product(gamma[:2],4),
              -gamma[0],-Killing_sign*identity]
    for r in range(4):
        actual=s.I**r*U*_product(gamma[:3-r],4)*gamma5**r*U
        assert _is_zero(s.simplify(actual-expected[r]))


def test_complete_vertex_covers_all_ground_signs_without_truncation_claim():
    vertex=complete_ground_contact_vertex();blocks=vertex['blocks']
    keys={(block['r'],tuple(block['Killing_signs'])) for block in blocks}
    assert keys=={(r,signs) for r in range(4) for signs in itertools.product((1,-1),repeat=4)}
    assert len(blocks)==64 and vertex['all_internal_sums_ordered']
    assert vertex['ground_vertex_projection_derived']
    assert not vertex['independent_contact_parameter_added']
    assert not vertex['nonlinear_KK_truncation_proved']
    assert not vertex['occupied_state_selected']
    assert not vertex['radial_stress_computed']
    assert not vertex['regular_child_proved']


@pytest.mark.parametrize('r,signs,dimension',[
    (-1,(1,1,1,1),12),(4,(1,1,1,1),12),(True,(1,1,1,1),12),
    (1,(1,1,1,1),3),(1,(1,1,1,1),True),
    (1,(1,1,1),12),(1,(1,1,1,0),12),(1,(1,1,1,True),12),
])
def test_invalid_clifford_channel_or_killing_sign_rejected(r,signs,dimension):
    with pytest.raises(ValueError):killing_ground_contact_block(r,signs,internal_dimension=dimension)
