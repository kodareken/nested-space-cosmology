"""Independent small-tree, closed-form and spectral tests of the theorem."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_infinite_inheritance import (
    InheritedTree,geometric_inherited_tree,homogeneous_dilation_source,
    volume_root_observables,
)


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


def test_geometric_action_derives_the_operator_and_root_compliance():
    b,q,L=3,.5,2.5;tree=geometric_inherited_tree(b,q,L)
    D,scales=tree.finite_matrix(2);count=len(D);K=np.zeros_like(D)
    levels=[0]*count
    for child in range(1,count):
        parent=(child-1)//b;levels[child]=levels[parent]+1
        conductance=1/(L*q**levels[parent])
        K[parent,parent]+=conductance;K[child,child]+=conductance
        K[parent,child]-=conductance;K[child,parent]-=conductance
    # The fixed parent port and Dirichlet child ports are part of the
    # same physical edge action; they are not fitted diagonal entries.
    K[0,0]+=q/L
    for vertex,depth in enumerate(levels):
        if depth==2:
            K[vertex,vertex]+=b/(L*q**depth)
    mass0=.5*L*(1/q+b);masses=mass0*q**np.asarray(levels)
    direct=K/np.sqrt(masses[:,None]*masses[None,:])
    np.testing.assert_allclose(D,direct,rtol=2e-15,atol=2e-15)
    assert tree.gamma==pytest.approx(b/(L*mass0))
    assert tree.contraction==pytest.approx(q/b)
    assert 1/(tree.gamma*mass0)==pytest.approx(L/b)
    response=tree.response(0.,6)
    assert abs(L/b-response['response'][0,0]/mass0)<=response['analytic_depth_bound']/mass0


def test_spatial_volume_action_and_both_infinite_response_branches():
    # Both sides of b*q=1, with finite total three-dimensional edge volume.
    for b,q in ((2,.6),(3,.25)):
        L,A=2.5,1.7
        tree=geometric_inherited_tree(b,q,L,spatial_dimension=3,cross_section=A)
        D,_=tree.finite_matrix(2);count=len(D)
        K=np.zeros_like(D);mass=np.zeros(count);levels=[0]*count
        for child in range(1,count):
            parent=(child-1)//b;levels[child]=levels[parent]+1
            length=L*q**levels[parent];section=A*q**(2*levels[parent])
            conductance=section/length;half_volume=section*length/2
            mass[parent]+=half_volume;mass[child]+=half_volume
            K[parent,parent]+=conductance;K[child,child]+=conductance
            K[parent,child]-=conductance;K[child,parent]-=conductance
        K[0,0]+=(A/q**2)/(L/q);mass[0]+=(A/q**2)*(L/q)/2
        for vertex,depth in enumerate(levels):
            if depth==2:
                length=L*q**depth;section=A*q**(2*depth)
                K[vertex,vertex]+=b*section/length
                mass[vertex]+=b*section*length/2
        np.testing.assert_allclose(D,K/np.sqrt(mass[:,None]*mass[None,:]),atol=2e-14)
        expected=L/(A*max(b,1/q));m0=mass[0]
        assert 1/(tree.gamma*m0)==pytest.approx(expected)
        response=tree.response(0.,24)
        assert abs(response['response'][0,0]/m0-expected)<=response['analytic_depth_bound']/m0+2e-15
        dilated=geometric_inherited_tree(b,q,2*L,spatial_dimension=3,cross_section=4*A)
        D2,_=dilated.finite_matrix(2)
        np.testing.assert_allclose(D2,D/4,atol=2e-14)
        assert 1/(dilated.gamma*8*m0)==pytest.approx(expected/2)
    with pytest.raises(ValueError,match='critical'):
        geometric_inherited_tree(2,.5,spatial_dimension=3)


def test_pressure_from_independent_lapse_variation_and_fixed_canonical_work():
    import sympy as sp

    N,a,Vc,Tdot,G=sp.symbols('N a Vc Tdot G',positive=True)
    action=a**3*Tdot/N-N*a*G;volume=a**3*Vc
    rho=-sp.diff(action,N)/volume
    pressure=sp.diff(action,a)/(N*sp.diff(volume,a))
    substitutions={N:sp.Rational(4,3),a:sp.Rational(3,2),Vc:sp.Rational(7,5),
                   Tdot:sp.Rational(2,7),G:sp.Rational(3,8)}
    proper_kinetic=a**3*Tdot/N**2
    T0=float((proper_kinetic*a**3).subs(substitutions))
    G0=float(G.subs(substitutions));scale=float(substitutions[a]);vc=float(substitutions[Vc])
    result=homogeneous_dilation_source(T0,G0,scale,vc)
    assert result['rho_geom']==pytest.approx(float(rho.subs(substitutions)))
    assert result['p_geom']==pytest.approx(float(pressure.subs(substitutions)))
    step=1e-5
    plus=homogeneous_dilation_source(T0,G0,scale+step,vc)
    minus=homogeneous_dilation_source(T0,G0,scale-step,vc)
    derivative=(plus['energy_length']-minus['energy_length'])/(plus['physical_volume']-minus['physical_volume'])
    assert derivative==pytest.approx(-result['p_geom'],rel=2e-9)
    for kinetic,gradient,expected_w,expected_q in ((0.,2.,-1/3,0.),(2.,0.,1.,2.),(2.,2.,1/3,1.)):
        source=homogeneous_dilation_source(kinetic,gradient)
        assert source['p_geom']/source['rho_geom']==pytest.approx(expected_w)
        assert source['flat_deceleration']==pytest.approx(expected_q)
        assert source['vacuum_form_defect']>0
        assert source['rho_plus_3p']==pytest.approx(source['rho_geom']+3*source['p_geom'])
    zero=homogeneous_dilation_source(0.,0.)
    assert zero['flat_deceleration'] is None
    assert zero['vacuum_form_defect']==0.


def test_static_root_work_distinguishes_fixed_field_and_fixed_force():
    import sympy as sp

    a,xi,chi,f=sp.symbols('a xi chi f',positive=True)
    energy=a*xi**2/(2*chi);force=sp.diff(energy,xi)
    assert sp.simplify(sp.diff(energy,a)-energy/a)==0
    force_controlled=xi.subs(xi,chi*f/a)
    stored=energy.subs(xi,force_controlled)
    work=(energy/a).subs(xi,force_controlled)+f*sp.diff(force_controlled,a)
    assert sp.simplify(sp.diff(stored,a)-work)==0
    assert sp.simplify(sp.diff(stored,a)+stored/a)==0
    assert sp.simplify(force.subs(xi,force_controlled)-f)==0
    legendre=(energy-force*xi).subs(xi,force_controlled)
    assert sp.simplify(legendre+stored)==0


def test_root_inertia_matches_the_full_static_profile_and_schur_derivative():
    for b,q in ((3,.5),(3,.4),(3,.25)):
        L,A=2.5,1.7
        tree=geometric_inherited_tree(b,q,L,spatial_dimension=3,cross_section=A)
        readout=volume_root_observables(b,q,L,A)
        m0=.5*A*L*(q**-3+b)
        profile_ratio=1/(b*q) if b*q>1 else 1.
        # Independent physical profile: sum m_n*phi_n^2 over all vertices.
        mass_ratio=b*q**3*profile_ratio**2
        independent_inertia=m0/(1-mass_ratio)
        assert readout['effective_root_inertia']==pytest.approx(independent_inertia)
        R=0.;slope=0.;c2=b/tree.omega*float(abs(tree.link[0,0])**2)
        for _ in range(240):
            updated=1/(tree.local[0,0].real-c2*R)
            slope=-updated**2*(1-c2*slope/tree.omega)
            R=updated
        schur_inertia=-m0*slope/R**2
        assert readout['effective_root_inertia']==pytest.approx(schur_inertia,rel=2e-13)
        # A finite, calibrated Laplace readout lies in the proved interval.
        x=.1;response=tree.response((x/L)**2,128)
        actual=1/(tree.gamma*response['response'][0,0].real)
        upper=1+readout['dimensionless_inertial_scale_squared']*x**2
        lower=upper-readout['dimensionless_s4_remainder_coefficient']*x**4
        depth_error=response['analytic_depth_bound']/(tree.gamma*response['response'][0,0].real**2)
        assert lower-2e-14<=actual<=upper+depth_error+2e-14
    half=volume_root_observables(3,.5)
    other=volume_root_observables(3,.4)
    assert half['static_compliance']==pytest.approx(other['static_compliance'])
    assert half['dimensionless_inertial_scale_squared']==pytest.approx(2.2)
    assert other['dimensionless_inertial_scale_squared']==pytest.approx(18.625/5.2)
    half_upper=1+half['dimensionless_inertial_scale_squared']*.01
    other_lower=1+other['dimensionless_inertial_scale_squared']*.01-other['dimensionless_s4_remainder_coefficient']*.0001
    assert other_lower>half_upper
