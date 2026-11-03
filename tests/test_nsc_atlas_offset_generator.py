import math
import mpmath as mp
import numpy as np
import pytest
import sympy as s
from recursive_horizons.nsc_atlas_offset_generator import (
    offset_value,phase_history,directions,euclidean_kernel_distance,closed_parent_point,sphere_cover)


def test_one_offset_iteration_matches_the_accumulated_closed_form():
    with mp.workdps(80):
        a=offset_value('golden');theta=phi=mp.mpf(0)
        for n in range(400):
            assert abs(mp.sin(mp.pi*(theta-n*a)))<mp.mpf('1e-74')
            assert abs(mp.sin(mp.pi*(phi-n*(n-1)*a/2)))<mp.mpf('1e-73')
            oldtheta=theta;theta,phi=mp.frac(theta+a),mp.frac(phi+theta)
            recovered=mp.frac(theta-a)
            assert abs(mp.sin(mp.pi*(recovered-oldtheta)))<mp.mpf('1e-75')


def test_polynomial_character_and_its_difference_have_the_claimed_orders():
    n,h,a,k,l=s.symbols('n h a k l')
    P=a*(k*n+l*n*(n-1)/2)
    assert s.Poly(P,n).coeff_monomial(n**2)==a*l/2
    assert s.expand(P.subs(n,n+h)-P).coeff(n)==a*l*h


def test_distance_kernel_is_a_gram_metric_not_an_extra_axis_assumption():
    p=phase_history(12,'action');v=directions(p)
    np.testing.assert_allclose(np.sum(v*v,axis=1),1,atol=5e-16)
    G=v@v.T
    assert np.linalg.matrix_rank(G,tol=1e-11)==3
    for i,j in ((1,5),(4,8),(0,9)):
        d=euclidean_kernel_distance(.4,.8,*(2*np.pi*p[i]),*(2*np.pi*p[j]))
        assert d==pytest.approx(np.linalg.norm(.4*v[i]-.8*v[j]),abs=4e-16)


def test_closed_parent_metric_and_volume_from_the_embedding():
    r,R,t,p=s.symbols('r R t p',positive=True)
    x=s.Matrix([R*s.cos(r/R),R*s.sin(r/R)*s.sin(t)*s.cos(p),
                R*s.sin(r/R)*s.sin(t)*s.sin(p),R*s.sin(r/R)*s.cos(t)])
    J=x.jacobian([r,t,p]);g=(J.T*J).applyfunc(s.trigsimp)
    assert g==s.diag(1,R**2*s.sin(r/R)**2,R**2*s.sin(r/R)**2*s.sin(t)**2)
    assert s.trigsimp(x.dot(x))==R**2
    np.testing.assert_allclose(closed_parent_point(0,[0,0,1]),[1,0,0,0])


def test_finite_line_cover_controls_ball_distance():
    # Octahedral endpoints: maximum direction angle acos(1/sqrt3).
    v=np.r_[np.eye(3),-np.eye(3)];h=sphere_cover(v)
    assert h==pytest.approx(math.acos(1/math.sqrt(3)))
    assert math.sin(h)==pytest.approx(math.sqrt(2/3))


def test_phase_sampling_and_metric_area_are_different_measures():
    t=s.symbols('t',real=True)
    assert s.integrate(s.cos(t)**2,(t,0,2*s.pi))/(2*s.pi)==s.Rational(1,2)
    assert s.integrate(s.sin(t)*s.cos(t)**2,(t,0,s.pi))/2==s.Rational(1,3)
    # A one-pass line retains zero transverse rank; the completion is essential.
    assert np.linalg.matrix_rank(directions(phase_history(1)))==1


def test_rational_offset_returns_a_finite_history_instead_of_dense_completion():
    from fractions import Fraction
    a=Fraction(1,4);theta=phi=Fraction(0);states=[]
    for _ in range(8):
        states.append((theta,phi));theta,phi=(theta+a)%1,(phi+theta)%1
    assert (theta,phi)==states[0]
    assert len(set(states))==8


def test_closed_meridian_cover_has_the_claimed_radial_dependence():
    direction=np.ones(3)/math.sqrt(3);cover=math.acos(1/math.sqrt(3))
    for radial in (.2,math.pi/2,2.8):
        target=closed_parent_point(radial,direction)
        best=math.atan2(math.sin(radial)/math.sqrt(3),math.cos(radial))
        nearest=closed_parent_point(best,[1,0,0])
        distance=math.acos(np.clip(target@nearest,-1,1))
        assert distance<=cover+2e-15
        if radial==math.pi/2:assert distance==pytest.approx(cover)


def test_parent_laplacian_modes_follow_the_generated_metric():
    r,R,t,p=s.symbols('r R t p',positive=True)
    def lap(f):
        return s.diff(f,r,2)+2*s.cot(r/R)*s.diff(f,r)/R+(s.diff(f,t,2)+s.cot(t)*s.diff(f,t)+s.diff(f,p,2)/s.sin(t)**2)/(R**2*s.sin(r/R)**2)
    first=s.cos(r/R)
    second=s.cos(r/R)*s.sin(r/R)*s.cos(t)
    assert s.trigsimp(lap(first)+3*first/R**2)==0
    assert s.trigsimp(lap(second)+8*second/R**2)==0


def test_first_three_irrational_offset_directions_have_rank_three():
    a=s.symbols('a',real=True)
    v0=s.Matrix([0,0,1]);v1=s.Matrix([s.sin(a),0,s.cos(a)])
    v2=s.Matrix([s.sin(2*a)*s.cos(a),s.sin(2*a)*s.sin(a),s.cos(2*a)])
    assert s.simplify(s.Matrix.hstack(v0,v1,v2).det()-s.sin(a)**2*s.sin(2*a))==0


def test_opposite_loops_have_integer_degrees_and_positive_equal_variation():
    from recursive_horizons.nsc_atlas_offset_generator import opposite_winding_path
    x = np.linspace(0, 2, 8001)
    for turns in (1, 3, -2):
        z = opposite_winding_path(x, 0, turns, 1.7)
        phase = np.unwrap(np.angle(z))
        assert (phase[4000]-phase[0])/(2*np.pi) == pytest.approx(turns)
        assert (phase[-1]-phase[4000])/(2*np.pi) == pytest.approx(-turns)
        dx = x[1]-x[0]
        derivative = np.diff(z)/dx
        halves = [.5*dx*np.sum(np.abs(d)**2) for d in np.split(derivative, 2)]
        np.testing.assert_allclose(halves, [2*np.pi**2*1.7**2*turns**2]*2, rtol=2e-6)


def test_unwinding_stays_closed_and_nonzero_including_the_join():
    from recursive_horizons.nsc_atlas_offset_generator import opposite_winding_path
    x = np.linspace(0, 2, 1001)
    for h in (0, .17, .5, .93, 1):
        z = opposite_winding_path(x, h)
        np.testing.assert_allclose(np.abs(z), 1, atol=3e-16)
        assert z[0] == pytest.approx(z[-1])
        assert z[500] == pytest.approx(np.exp(2j*np.pi*(1-h)))
    np.testing.assert_allclose(opposite_winding_path(x, 1), np.ones(x.size))


def test_unwinding_variation_matches_independent_integration_and_sharp_bound():
    from recursive_horizons.nsc_atlas_offset_generator import winding_variation
    # Cauchy-Schwarz on each length-one half: integral f'^2 >= delta_f^2.
    h,m,R,x = s.symbols('h m R x',real=True)
    slope = 2*s.pi*m*(1-h)
    direct = s.integrate(R**2*slope**2/2,(x,0,1))*2
    assert s.simplify(direct-4*s.pi**2*R**2*m**2*(1-h)**2)==0
    assert winding_variation(.5,3,2)==pytest.approx(float(direct.subs({h:.5,m:3,R:2})))
    assert winding_variation(1)==0


def test_two_balanced_positive_outputs_force_equal_antipodal_geometry():
    p,q = s.symbols('p q',positive=True)
    # If p*u+q*v=0 and both directions have unit norm, p^2=q^2.
    assert s.factor(p**2-q**2)==(p-q)*(p+q)
    assert s.solve([p+q-1,p-q],[p,q])=={p:s.Rational(1,2),q:s.Rational(1,2)}
    from recursive_horizons.nsc_atlas_offset_generator import balanced_output_directions
    rays = balanced_output_directions(2)
    np.testing.assert_allclose(rays[0]+rays[1],0,atol=2e-16)
    assert rays[0]@rays[1]==pytest.approx(-1)


def test_balance_alone_also_allows_more_than_two_outputs():
    from recursive_horizons.nsc_atlas_offset_generator import balanced_output_directions
    for count in (3,4,5,17):
        rays=balanced_output_directions(count)
        np.testing.assert_allclose(np.linalg.norm(rays,axis=1),1,atol=2e-16)
        np.testing.assert_allclose(rays.mean(axis=0),0,atol=2e-16)


def test_sum_product_and_concatenation_are_distinct_contact_rules():
    t=s.symbols('t',real=True)
    forward=s.exp(s.I*t);reverse=s.exp(-s.I*t)
    assert s.simplify(forward*reverse)==1
    assert s.simplify(s.expand_complex(forward+reverse)-2*s.cos(t))==0
    assert (forward+reverse).subs(t,0)==2
    assert s.simplify((forward+reverse).subs(t,s.pi/2))==0


@pytest.mark.parametrize('args',[(0,0,0,1),(0,0,1,-1),(0,1.1,1,1),(3,0,1,1),(0,0,1,float('nan'))])
def test_winding_domain_rejects_unstated_or_invalid_geometry(args):
    from recursive_horizons.nsc_atlas_offset_generator import opposite_winding_path
    with pytest.raises(ValueError):opposite_winding_path(*args)


def test_conservative_winding_motion_solves_wave_equation_away_from_corners():
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport
    # A smooth initial phase checks the derivation independently of the
    # piecewise-linear reference solution (whose second derivatives vanish locally).
    x,t,k=s.symbols('x t k',real=True)
    phase=(s.cos(k*(x-t))+s.cos(k*(x+t)))/2
    assert s.simplify(s.diff(phase,t,2)-s.diff(phase,x,2))==0
    for closed in (False,True):
        v=winding_transport(np.array([.13,.73,1.33]),.19,closed=closed)
        dx=1e-6
        numerical=(winding_transport(np.array([.13,.73,1.33])+dx,.19,closed=closed)['phase']-
                   winding_transport(np.array([.13,.73,1.33])-dx,.19,closed=closed)['phase'])/(2*dx)
        np.testing.assert_allclose(numerical,v['drawing_derivative'],atol=1e-8)


def test_local_conservative_ledger_and_characteristic_densities_are_exact():
    x,t,R=s.symbols('x t R',real=True)
    f=s.Function('f')(x-t);g=s.Function('g')(x+t)
    phase=(f+g)/2
    velocity=s.diff(phase,t);spatial=s.diff(phase,x)
    density=R**2*(velocity**2+spatial**2)/2
    flux=-R**2*velocity*spatial
    assert s.simplify(s.diff(density,t)+s.diff(flux,x))==0
    assert s.simplify(density-R**2*(s.diff(f,x)**2+s.diff(g,x)**2)/4)==0


def test_open_and_closed_ledger_conserve_the_same_initial_budget():
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport_ledger
    for closed in (False,True):
        for offset in (0,.13,.25,.5,.83,1,1.7,2,3.19):
            row=winding_transport_ledger(offset,turns=2,radius=.7,closed=closed)
            assert row['density']==pytest.approx(4*np.pi**2*4*.7**2,rel=3e-15)
            assert row['right_density']==pytest.approx(row['initial_budget']/2,rel=3e-15)
            assert row['left_density']==pytest.approx(row['initial_budget']/2,rel=3e-15)
            assert row['kinetic']+row['variation']==pytest.approx(row['density'],rel=2e-15)


def test_open_packets_separate_and_transport_the_budget_in_opposite_directions():
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport
    t=1.3
    x=np.linspace(.001,1.999,2000)
    for sign in (-1,1):
        packet=winding_transport(x+sign*t,t)
        np.testing.assert_allclose(packet['flux'],sign*packet['density'],atol=1e-14)
        assert np.all(packet['density']>0)
    middle=winding_transport(1.,t)
    assert middle['density']==0
    assert middle['phase']==0


def test_closed_contact_flattens_then_returns_instead_of_losing_the_budget():
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport,winding_transport_ledger
    x=np.linspace(.001,1.999,137)
    np.testing.assert_allclose(winding_transport(x,.5,closed=True)['phase'],np.pi,atol=2e-15)
    mid=winding_transport_ledger(.5,closed=True)
    assert mid['variation']==pytest.approx(0,abs=1e-28)
    assert mid['kinetic']==pytest.approx(mid['initial_budget'])
    np.testing.assert_allclose(winding_transport(x,2.,closed=True)['phase'],
                               winding_transport(x,0.,closed=True)['phase'],atol=2e-15)
    np.testing.assert_allclose(winding_transport(x,2.,closed=True)['offset_derivative'],0)


def test_initial_velocity_and_both_winding_halves_are_retained():
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport,opposite_winding_path
    x=np.linspace(0,2,101)
    for closed in (False,True):
        row=winding_transport(x,0.,turns=-3,closed=closed)
        np.testing.assert_allclose(row['circle'],opposite_winding_path(x,0.,turns=-3),atol=1e-14)
        np.testing.assert_allclose(row['offset_derivative'],0)


def test_receiving_amounts_match_boundary_flux_and_the_remaining_initial_region():
    from scipy.integrate import quad
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport,winding_transport_ledger
    for offset in (.3,1.3,2.,3.):
        row=winding_transport_ledger(offset)
        cuts=[0.,*sorted({a for a in (1.,2.) if a<offset}),offset]
        for x,orientation,key in ((0.,-1.,'received_left_of_initial_interval'),
                                  (2.,1.,'received_right_of_initial_interval')):
            delivered=sum(quad(lambda t:orientation*float(winding_transport(x,t)['flux']),a,b)[0]
                          for a,b in zip(cuts[:-1],cuts[1:]))
            assert delivered==pytest.approx(row[key],abs=2e-12)
        spatial_cuts=sorted({0.,2.,*[a for a in
                            [k+sign*offset for k in (0.,1.,2.) for sign in (-1,1)] if 0<a<2]})
        retained=sum(quad(lambda x:float(winding_transport(x,offset)['density']),a,b)[0]
                     for a,b in zip(spatial_cuts[:-1],spatial_cuts[1:]))
        assert retained==pytest.approx(row['retained_in_initial_interval'],abs=2e-12)


@pytest.mark.parametrize('offset,closed',[(-.1,False),(float('nan'),False),(.2,'open')])
def test_transport_rejects_ambiguous_domain(offset,closed):
    from recursive_horizons.nsc_atlas_offset_generator import winding_transport
    with pytest.raises(ValueError):winding_transport(0.,offset,closed=closed)


def test_lifted_offset_flow_matches_every_discrete_update_and_composes():
    from recursive_horizons.nsc_atlas_offset_generator import lifted_offset_flow
    alpha=math.sqrt(2)/3
    p,q=np.array([0.,.2,-.3]),np.array([.1,-.4,.7])
    oldp,oldq=p.copy(),q.copy()
    for n in range(21):
        gotp,gotq=lifted_offset_flow(oldp,oldq,alpha,n)
        np.testing.assert_allclose(gotp,p,atol=3e-14)
        np.testing.assert_allclose(gotq,q,atol=3e-13)
        p,q=p+alpha,q+p
    once=lifted_offset_flow(oldp,oldq,alpha,.17+.63)
    twice=lifted_offset_flow(*lifted_offset_flow(oldp,oldq,alpha,.17),alpha,.63)
    np.testing.assert_allclose(once,twice,atol=3e-16)


def test_exact_offset_generator_and_phase_area_preservation_symbolically():
    p,q,a,t=s.symbols('p q a t',real=True)
    P=p+a*t;Q=q+t*p+a*t*(t-1)/2
    assert s.diff(P,t)==a
    assert s.simplify(s.diff(Q,t)-(P-a/2))==0
    assert s.Matrix([P,Q]).jacobian([p,q]).det()==1
    generator=(P-a/2)**2/2-a*Q
    assert s.simplify(s.diff(generator,t))==0


def test_different_retained_lifts_agree_at_steps_but_not_between_them():
    from recursive_horizons.nsc_atlas_offset_generator import lifted_offset_flow
    # theta and theta+1 are the same torus label; fractional interpolation
    # distinguishes their retained histories. At integer offsets it does not.
    for t in (1.,2.,3.):
        _,q=lifted_offset_flow(.2,.3,.17,t)
        _,q2=lifted_offset_flow(1.2,.3,.17,t)
        assert q2-q==pytest.approx(t)
    _,q=lifted_offset_flow(.2,.3,.17,.5)
    _,q2=lifted_offset_flow(1.2,.3,.17,.5)
    assert q2-q==pytest.approx(.5)


def test_uncoupled_update_changes_common_phase_not_the_winding_shape():
    from recursive_horizons.nsc_atlas_offset_generator import offset_driven_winding,opposite_winding_path
    x=np.linspace(0,2,500)
    for n in (0.,1.,7.,100.):
        row=offset_driven_winding(x,n,.17,0.)
        background=.17*n*(n-1)/2
        np.testing.assert_allclose(row['circle']*np.exp(-2j*np.pi*background),
                                   opposite_winding_path(x,0.),atol=2e-12)
        np.testing.assert_allclose(row['theta_lift'],.17*n)
        np.testing.assert_allclose(row['relative_offset_derivative'],0)


def test_coupled_extension_obeys_the_first_offset_equation():
    from recursive_horizons.nsc_atlas_offset_generator import offset_driven_winding
    x=np.array([.13,.73,1.33]);t=.19;dt=1e-6;alpha=.17
    for speed in (0.,.5,1.,2.):
        row=offset_driven_winding(x,t,alpha,speed)
        numerical=(offset_driven_winding(x,t+dt,alpha,speed)['phi_lift']-
                   offset_driven_winding(x,t-dt,alpha,speed)['phi_lift'])/(2*dt)
        np.testing.assert_allclose(numerical,row['theta_lift']-alpha/2,atol=1e-9)


def test_wave_solution_plus_background_gives_the_coupled_offset_equations():
    x,t,a,c=s.symbols('x t a c',real=True)
    relative=s.cos(x-c*t)/2+s.cos(x+c*t)/2
    phi=relative+a*t*(t-1)/2
    theta=s.diff(phi,t)+a/2
    assert s.simplify(s.diff(theta,t)-a-c**2*s.diff(phi,x,2))==0


def test_original_uniform_history_does_not_select_a_transfer_fraction():
    from recursive_horizons.nsc_atlas_offset_generator import offset_connection_calculation
    record=offset_connection_calculation();rows=record['same_history_different_transfer']
    assert [r['fraction_received_outside_initial_interval'] for r in rows]==pytest.approx([0.,.25,.5,1.])
    assert len({(r['uniform_theta'],r['uniform_phi']) for r in rows})==1
    assert record['original_pointwise_stationary_winding']['transported_fraction']==0


def test_parent_child_compression_composes_across_three_nested_regions():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    u=np.array([0.,.13,.5,1.])
    maps=[np.exp(-k*u)[:,None,None]*np.eye(3) for k in (1.,2.,4.)]
    a=inherited_compression_offset(maps[0],maps[1])['compression']
    b=inherited_compression_offset(maps[1],maps[2])['compression']
    c=inherited_compression_offset(maps[0],maps[2])['compression']
    np.testing.assert_allclose(a,u,atol=4e-16)
    np.testing.assert_allclose(b,2*u,atol=4e-16)
    np.testing.assert_allclose(c,a+b,atol=5e-16)


def test_common_frame_change_cancels_from_relative_compression():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    parent=np.array([np.eye(3),np.diag([.7,.9,.8])])
    child=np.array([np.eye(3),np.diag([.4,.6,.9])])
    frame=np.array([np.diag([2.,.3,1.7]),np.array([[1.,.4,0.],[0.,2.,.3],[0.,0.,.5]])])
    before=inherited_compression_offset(parent,child)['compression']
    after=inherited_compression_offset(frame@parent,frame@child)['compression']
    np.testing.assert_allclose(before,after,atol=2e-16)


def test_isotropic_compression_is_log_length_not_an_independent_phase_knob():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    parent=np.array([np.eye(3),np.eye(3)])
    child=np.array([np.eye(3),np.eye(3)/2])
    row=inherited_compression_offset(parent,child)
    assert row['compression'][-1]==pytest.approx(math.log(2))
    np.testing.assert_allclose(row['increments'],[math.log(2)])
    reverse=inherited_compression_offset(child,parent)
    np.testing.assert_allclose(reverse['compression'],-row['compression'])


def test_common_compression_adds_no_relative_offset_but_ancestor_offset_survives():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    root=np.array([np.eye(3),np.eye(3)])
    parent=np.array([np.eye(3),np.eye(3)/2])
    same_child=parent.copy()
    local=inherited_compression_offset(parent,same_child)['compression']
    inherited=inherited_compression_offset(root,parent)['compression']
    total=inherited_compression_offset(root,same_child)['compression']
    np.testing.assert_allclose(local,0)
    np.testing.assert_allclose(total,inherited)
    assert total[-1]>0


def test_log_jacobians_handle_strong_compression_without_determinant_underflow():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    parent=np.array([np.eye(3),np.eye(3)])
    child=np.array([np.eye(3),np.eye(3)*math.exp(-400)])
    assert inherited_compression_offset(parent,child)['compression'][-1]==pytest.approx(400)


def test_pure_rotation_does_not_masquerade_as_volume_compression():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    a=.7
    rotation=np.array([[math.cos(a),-math.sin(a),0],[math.sin(a),math.cos(a),0],[0,0,1]])
    parent=np.array([np.eye(3),np.eye(3)])
    child=np.array([np.eye(3),rotation])
    np.testing.assert_allclose(inherited_compression_offset(parent,child)['compression'],0,atol=1e-16)


def test_relative_volume_rate_has_the_logarithmic_derivative():
    u,a,b=s.symbols('u a b',real=True)
    p=s.exp(-3*a*u);c=s.exp(-3*b*u)
    rate=(s.diff(p,u)/p-s.diff(c,u)/c)/3
    assert s.simplify(rate-(b-a))==0


@pytest.mark.parametrize('child',[np.zeros((2,3,3)),np.array([np.eye(3),np.diag([-1.,1.,1.])]),np.ones((2,3,2))])
def test_inherited_offset_requires_admissible_corresponding_patches(child):
    from recursive_horizons.nsc_atlas_offset_generator import inherited_compression_offset
    with pytest.raises(ValueError):inherited_compression_offset(np.array([np.eye(3),np.eye(3)]),child)


def test_intrinsic_circle_phase_cost_is_invariant_under_reparameterization():
    from scipy.integrate import quad
    # theta=u+e sin(u), metric R^2 theta'^2, phase=m theta.
    for radius,turns,e in ((1.,1,.3),(.7,3,.6),(2.,-2,-.4)):
        value=quad(lambda u:.5*(turns*(1+e*math.cos(u)))**2 /
                   (radius**2*(1+e*math.cos(u))**2)*radius*(1+e*math.cos(u)),0,2*math.pi)[0]
        assert value==pytest.approx(math.pi*turns**2/radius,rel=2e-14)


def test_affine_intrinsic_cost_matches_direct_transformed_phase_gradient():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost
    A=np.array([[.7,.2,0.],[0.,1.3,.1],[.1,0.,.8]])
    # psi=x^2+y^2+z^2 on the unit cube centred at zero has M=I/3.
    value=intrinsic_phase_cost(A,np.eye(3)/math.sqrt(3))['cost']
    from numpy.polynomial.legendre import leggauss
    nodes,weights=leggauss(3);direct=0.
    for i,x in enumerate(nodes/2):
        for j,y in enumerate(nodes/2):
            for k,z in enumerate(nodes/2):
                gradient=np.linalg.solve(A.T,2*np.array([x,y,z]))
                direct+=weights[i]*weights[j]*weights[k]/8*np.linalg.det(A)*np.dot(gradient,gradient)/2
    assert value==pytest.approx(direct,rel=2e-15)


def test_intrinsic_cost_uniform_scaling_and_rigid_rotation_controls():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost
    for dimension in (1,2,3):
        for scale in (.25,.5,2.):
            value=intrinsic_phase_cost(scale*np.eye(dimension),np.eye(dimension))['cost']
            assert value==pytest.approx(dimension/2*scale**(dimension-2),rel=3e-15)
    R=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    A=np.diag([.7,1.3,2.]);C=np.diag([1.,2.,3.])
    assert intrinsic_phase_cost(R@A,C)['cost']==pytest.approx(intrinsic_phase_cost(A,C)['cost'])


def test_directional_response_matches_cost_derivatives_not_pressure_labels():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost
    A=np.diag([.7,1.3,2.]);C=np.diag([1.,2.,3.])
    row=intrinsic_phase_cost(A,C);eps=1e-6
    for i in range(3):
        positive=np.eye(3);negative=np.eye(3)
        positive[i,i]=math.exp(eps);negative[i,i]=math.exp(-eps)
        derivative=(intrinsic_phase_cost(positive@A,C)['cost']-
                    intrinsic_phase_cost(negative@A,C)['cost'])/(2*eps)
        assert derivative==pytest.approx(-row['cost']*row['logarithmic_force_per_cost'][i],abs=2e-9)


def test_isotropic_pattern_holds_equal_axes_at_fixed_volume_but_directional_pattern_does_not():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost
    for shape in (.1,.5,1.,2.):
        A=np.diag([math.exp(-shape/2),math.exp(-shape/2),math.exp(shape)])
        assert np.linalg.det(A)==pytest.approx(1)
        assert intrinsic_phase_cost(A,np.eye(3))['cost']>1.5
        assert intrinsic_phase_cost(A,[[0.,0.,1.]])['cost']==pytest.approx(.5*math.exp(-2*shape))
    shape=s.symbols('shape',real=True)
    isotropic=(2*s.exp(shape)+s.exp(-2*shape))/2
    assert s.diff(isotropic,shape).subs(shape,0)==0
    assert s.diff(isotropic,shape,2).subs(shape,0)==3


def test_conservative_neck_motion_satisfies_shape_equation_and_energy_budget():
    from recursive_horizons.nsc_atlas_offset_generator import conservative_neck_shape
    t=s.symbols('t',real=True);a=s.sqrt(s.Rational(2,3))
    shape=s.log(s.cosh(a*t))
    assert s.simplify(s.diff(shape,t,2)-s.Rational(2,3)*s.exp(-2*shape))==0
    for t0 in (0.,.3,1.,2.,4.,1000.):
        row=conservative_neck_shape(t0)
        assert row['total']==pytest.approx(.5,abs=3e-16)
        assert row['log_volume']==0
        assert row['log_axial_scale']>=0
        assert row['log_transverse_scale']<=0
    assert conservative_neck_shape(0)['shape_velocity']==0


def test_intrinsic_phase_zero_source_and_large_compression_remain_explicit():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost
    assert intrinsic_phase_cost(np.eye(3),np.zeros((1,3)))['cost']==0
    row=intrinsic_phase_cost(math.exp(-300)*np.eye(3),np.eye(3))
    assert row['log_cost']==pytest.approx(math.log(1.5)-300)


def test_inherited_compression_drives_intrinsic_cost_without_an_independent_alpha():
    from recursive_horizons.nsc_atlas_offset_generator import intrinsic_phase_cost,inherited_compression_offset
    for dimension in (1,2,3):
        parent=np.array([np.eye(dimension),math.exp(-.2)*np.eye(dimension)])
        child=np.array([np.eye(dimension),math.exp(-.7)*np.eye(dimension)])
        chi=inherited_compression_offset(parent,child)['compression'][-1]
        C=np.eye(dimension)
        qp=intrinsic_phase_cost(parent[-1],C)['cost']
        qc=intrinsic_phase_cost(child[-1],C)['cost']
        assert qc/qp==pytest.approx(math.exp(-(dimension-2)*chi),rel=2e-15)


def test_generator_least_stretch_cost_resists_compression_in_all_dimensions():
    from recursive_horizons.nsc_atlas_offset_generator import (
        generator_least_stretch_cost,inherited_compression_offset)
    for dimension in (1,2,3):
        parent=np.array([np.eye(dimension),math.exp(-.2)*np.eye(dimension)])
        child=np.array([np.eye(dimension),math.exp(-.9)*np.eye(dimension)])
        chi=inherited_compression_offset(parent,child)['compression'][-1]
        C=np.eye(dimension)
        qp=generator_least_stretch_cost(parent[-1],C)['stretch_cost']
        qc=generator_least_stretch_cost(child[-1],C)['stretch_cost']
        assert qc/qp==pytest.approx(math.exp(chi),rel=2e-15)
    R=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
    A=np.diag([.6,1.1,1.9]);C=np.diag([1.,2.,3.])
    assert generator_least_stretch_cost(R@A,C)['stretch_cost']==pytest.approx(
        generator_least_stretch_cost(A,C)['stretch_cost'],rel=2e-15)


def test_least_resistance_compression_minimum_matches_symbolic_stationarity_and_hessian():
    from recursive_horizons.nsc_atlas_offset_generator import (
        least_resistance_compression_minimum,generator_least_stretch_cost)
    l1,l2,l3,m1,m2,m3,P=s.symbols('l1 l2 l3 m1 m2 m3 P',positive=True)
    E=m1/l1+m2/l2+m3/l3+P*l1*l2*l3
    g=(m1*m2*m3*P)**s.Rational(1,4)
    sol={l1:m1/g,l2:m2/g,l3:m3/g}
    for var in (l1,l2,l3):
        assert s.simplify(s.diff(E,var).subs(sol))==0
    iso=least_resistance_compression_minimum([.5,.5,.5],1.)
    assert iso['outcome']=='uniform_contraction'
    assert iso['equilibrium_scales']==pytest.approx([2**(-.25)]*3,rel=2e-15)
    assert iso['virial_ratio_stretch_to_work']==pytest.approx(3.,rel=2e-15)
    raw=least_resistance_compression_minimum([.25,.25,.5],1.)
    assert raw['outcome']=='anisotropic_neck_equilibrium'
    assert raw['axial_to_transverse_aspect_ratio']==pytest.approx(2.,rel=2e-15)
    assert raw['log_axial_neck_parameter']==pytest.approx(math.log(2)/3,rel=2e-15)
    # Perturbing any principal scale around the minimum strictly raises total cost.
    opt=np.array(raw['equilibrium_scales'])
    C=np.diag(np.sqrt(2*np.array([.25,.25,.5])))
    base=generator_least_stretch_cost(np.diag(opt),C,enclosing_load=1.)['total_cost']
    assert base==pytest.approx(raw['total_cost'],rel=2e-15)
    for delta in (np.array([.05,-.03,.02]),np.array([-.04,-.04,-.04]),np.array([-.02,.02,.04])):
        perturbed=generator_least_stretch_cost(np.diag(opt*(1+delta)),C,enclosing_load=1.)['total_cost']
        assert perturbed>base
    seed=least_resistance_compression_minimum([0.,0.,.5],1.)
    assert seed['outcome']=='neck_pinch_opening_break'
    assert seed['equilibrium_scales'] is None


def test_neck_break_continuation_transforms_throat_pinch_into_child_friedmann_expansion():
    from recursive_horizons.nsc_atlas_offset_generator import continue_across_neck_break
    t=s.symbols('t',real=True);a=s.sqrt(s.Rational(2,3))
    Rc=s.cosh(a*t)
    Hc=s.diff(Rc,t)/Rc
    assert s.simplify(Hc**2+s.Rational(2,3)/Rc**2-s.Rational(2,3))==0
    at_six=continue_across_neck_break(6.)
    assert at_six['parent_frame']['velocity_fraction_of_asymptotic_wall']>0.9998
    assert at_six['parent_frame']['transverse_radius']<0.13
    assert at_six['child_frame']['scale_factor']==pytest.approx(
        1/at_six['parent_frame']['transverse_radius']**2,rel=2e-15)
    assert at_six['child_frame']['friedmann_residual']==pytest.approx(0.,abs=2e-15)
    assert at_six['child_frame']['effective_child_lambda']==pytest.approx(2.,rel=2e-15)
    assert at_six['child_frame']['transferred_budget_fraction']>0.9997


def test_spiral_decay_to_turn_comparison_does_not_impose_golden_ratio():
    from recursive_horizons.nsc_atlas_offset_generator import (
        spiral_decay_to_turn_comparison,inherited_offset_calculation,intrinsic_compression_calculation)
    comp=spiral_decay_to_turn_comparison()
    assert comp['static_compression_produces_rotation'] is False
    phi=(1+math.sqrt(5))/2
    assert comp['golden_spiral_decay_per_radian']==pytest.approx(2*math.log(phi)/math.pi,rel=2e-15)
    ga=comp['golden_angle_increment']
    assert abs(ga['difference_from_golden_spiral_ratio'])>0.1
    assert ga['lobe_asymmetry_eta']==pytest.approx(
        (1-ga['opposite_lobe_amplitude_ratio'])/(1+ga['opposite_lobe_amplitude_ratio']),rel=2e-15)
    inh=inherited_offset_calculation()
    assert inh['mechanism_evidence'] is False
    assert 'bookkeeping' in inh['bookkeeping_note']
    intr=intrinsic_compression_calculation()
    assert 'supplied' in intr['coupling_c_status']
    assert 'supplied' in intr['increment_alpha_status']


def test_helical_buckling_derives_nonzero_turning_above_critical_compression():
    from recursive_horizons.nsc_atlas_offset_generator import helical_buckling_under_compression
    v,lam,Q0,P=s.symbols('v lam Q0 P',positive=True)
    E=Q0/s.sqrt(lam**2+v**2)+P*v**2/2
    v_star=s.sqrt((Q0/P)**s.Rational(2,3)-lam**2)
    assert s.simplify(s.diff(E,v).subs(v,v_star))==0
    pre=helical_buckling_under_compression(-0.3,stretch_budget=1.,lateral_confinement=1.)
    assert pre['buckled'] is False
    assert pre['angular_rate_omega']==0.
    assert pre['straight_transverse_hessian']>0
    post=helical_buckling_under_compression(0.8,stretch_budget=1.,lateral_confinement=1.,radius=1.)
    assert post['buckled'] is True
    assert post['straight_transverse_hessian']<0
    assert post['angular_rate_omega']==pytest.approx(math.sqrt(1-math.exp(-1.6)),rel=2e-15)
    assert post['derived_alpha']==pytest.approx(post['angular_rate_omega']/(2*math.pi),rel=2e-15)
    assert post['minimum_cost']<post['straight_cost']
    assert post['cost_relief']>0


def test_damped_horizon_oscillation_conserves_budget_and_matches_lobe_integrals():
    from scipy.integrate import quad
    from recursive_horizons.nsc_atlas_offset_generator import (
        damped_horizon_oscillation,buckling_to_goldilocks_calculation)
    t,g,w,A=s.symbols('t g w A',positive=True)
    q=A*s.exp(-g*t)*s.sin(w*t)
    qt=s.diff(q,t);qtt=s.diff(q,t,2)
    w0_sq=g**2+w**2
    assert s.simplify(qtt+2*g*qt+w0_sq*q)==0
    ew=(qt**2+w0_sq*q**2)/2
    assert s.simplify(s.diff(ew,t)+2*g*qt**2)==0
    gamma=0.35;omega=2*math.pi;A0=1.3
    for tau in (0.4,1.1,2.7,5.0):
        row=damped_horizon_oscillation(tau,gamma,omega,amplitude=A0,enclosing_load=1.2)
        assert row['conservation_residual']==pytest.approx(0.,abs=2e-15)
        dissipated=quad(lambda u:2*gamma*damped_horizon_oscillation(u,gamma,omega,amplitude=A0)['velocity']**2,0,tau)[0]
        assert dissipated==pytest.approx(row['transferred_geometry_budget'],rel=1e-10)
    half_period=math.pi/omega
    plus_0=quad(lambda u:A0*math.exp(-gamma*u)*math.sin(omega*u),0,half_period)[0]
    minus_1=quad(lambda u:-A0*math.exp(-gamma*u)*math.sin(omega*u),half_period,2*half_period)[0]
    ref=damped_horizon_oscillation(4.,gamma,omega,amplitude=A0)
    assert plus_0==pytest.approx(ref['first_plus_lobe_area'],rel=1e-12)
    assert minus_1==pytest.approx(ref['first_minus_lobe_area'],rel=1e-12)
    assert (plus_0-minus_1)/(plus_0+minus_1)==pytest.approx(ref['lobe_asymmetry_eta'],rel=1e-12)
    def deposited_work(t):
        velocity = damped_horizon_oscillation(t, gamma, omega, amplitude=A0)['velocity']
        return 2 * gamma * velocity**2

    plus_work = quad(deposited_work, 0, half_period)[0]
    minus_work = quad(deposited_work, half_period, 2*half_period)[0]
    assert ref['first_plus_lobe_work'] == pytest.approx(plus_work)
    assert ref['first_minus_lobe_work'] == pytest.approx(minus_work)
    assert ref['work_asymmetry_eta'] == pytest.approx((plus_work-minus_work)/(plus_work+minus_work))
    assert (ref['cumulative_plus_area']-ref['cumulative_minus_area'])==pytest.approx(
        ref['net_asymmetric_remainder'],rel=2e-15)
    assert damped_horizon_oscillation(0.,0.,omega)['lobe_asymmetry_eta']==0.
    pipeline=buckling_to_goldilocks_calculation()
    assert pipeline['stage_3_goldilocks_equilibrium']['completion_at_offset_4']>0.98


def test_two_force_zero_crossing_smoothly_converts_plus_implosion_into_minus_expansion():
    from recursive_horizons.nsc_atlas_offset_generator import two_force_zero_crossing
    t,R0,Fin,Fout=s.symbols('t R0 Fin Fout',positive=True)
    Fnet=Fin-Fout
    R=R0-Fnet*t**2/2
    assert s.simplify(s.diff(R,t,2)+Fnet)==0
    plus=two_force_zero_crossing(0.5,3.,1.,initial_radius=1.)
    assert plus['region']=='+'
    assert plus['signed_radius']>0
    assert plus['physical_expansion_rate']<0
    assert plus['orientation_sign']==1
    assert plus['conservation_residual']==pytest.approx(0.,abs=2e-15)
    horizon=two_force_zero_crossing(1.,3.,1.,initial_radius=1.)
    assert horizon['region']=='0'
    assert horizon['signed_radius']==pytest.approx(0.,abs=2e-15)
    assert horizon['jacobian_determinant_ratio']==pytest.approx(0.,abs=2e-15)
    assert horizon['horizon_kinetic_budget']==pytest.approx(2.,rel=2e-15)
    minus=two_force_zero_crossing(2.,3.,1.,initial_radius=1.)
    assert minus['region']=='-'
    assert minus['signed_radius']==pytest.approx(-3.,rel=2e-15)
    assert minus['physical_radius']==pytest.approx(3.,rel=2e-15)
    assert minus['physical_expansion_rate']==pytest.approx(4.,rel=2e-15)
    assert minus['orientation_sign']==-1
    assert minus['jacobian_determinant_ratio']==pytest.approx(-27.,rel=2e-15)
    assert minus['conservation_residual']==pytest.approx(0.,abs=2e-15)


def test_spaghettification_reduces_infall_to_1d_and_isoperimetric_selects_3d_sphere():
    from recursive_horizons.nsc_atlas_offset_generator import spaghettification_to_spherical_winner
    at_zero=spaghettification_to_spherical_winner(0.)
    assert at_zero['infall_spaghettification']['spaghettification_aspect_ratio']==1.
    assert at_zero['isoperimetric_spherical_winner']['deformed_isoperimetric_quotient']==1.
    assert at_zero['cooker_volume_takeover']['child_minus_drop_fraction']==0.
    for chi in (.25,.5,1.,2.,4.):
        row=spaghettification_to_spherical_winner(chi)
        inf=row['infall_spaghettification']
        assert inf['transverse_stretch']**2*inf['radial_stretch']==pytest.approx(1.,rel=2e-15)
        assert inf['spaghettification_aspect_ratio']==pytest.approx(math.exp(3*chi),rel=2e-15)
        iso=row['isoperimetric_spherical_winner']
        assert iso['deformed_isoperimetric_quotient']>1.
        assert iso['sphere_area_saving_fraction']>0.
    deep=spaghettification_to_spherical_winner(4.)
    assert deep['infall_spaghettification']['normalized_deformation_diagonal'][0]<1e-5
    assert deep['cooker_volume_takeover']['child_minus_drop_fraction']>0.99999


def test_unified_single_cost_minimizes_all_five_variables_without_fluid_assumptions():
    from recursive_horizons.nsc_atlas_offset_generator import unified_least_resistance_cycle
    l3,tl1,tl2,tl3,th,X=s.symbols('l3 tl1 tl2 tl3 th X',positive=True)
    E=(2*(1-th)**2*X+(1-th)**2/l3+l3/X**2)+(th**2/tl1+th**2/tl2+th**2/tl3+tl1*tl2*tl3)
    vars5=[l3,tl1,tl2,tl3,th]
    grad=[s.diff(E,v) for v in vars5]
    hess=s.Matrix([[s.diff(E,u,v) for v in vars5] for u in vars5])
    # Exact algebraic check at chi = ln(2)/2 (X = sqrt(2)), where sqrt(16*X^2+17) = 7:
    sol_half={X:s.sqrt(2),th:s.Rational(1,2),l3:s.sqrt(2)/2,
              tl1:1/s.sqrt(2),tl2:1/s.sqrt(2),tl3:1/s.sqrt(2)}
    for g in grad:
        assert s.simplify(g.subs(sol_half))==0
    hess_num=np.array(hess.subs(sol_half).evalf(),dtype=float)
    assert np.all(np.linalg.eigvalsh(hess_num)>0)
    at_thresh=unified_least_resistance_cycle(.5*math.log(2))
    assert at_thresh['transferred_winding_fraction_theta']==pytest.approx(.5,rel=2e-15)
    assert at_thresh['plus_infall']['filament_aspect_ratio']==pytest.approx(1.,rel=2e-15)
    for chi in (.5,1.,2.,4.):
        row=unified_least_resistance_cycle(chi)
        assert row['fluid_incompressibility_assumed'] is False
        assert row['surface_tension_assumed'] is False
        assert row['volume_conservation_shortcut_assumed'] is False
        assert row['plus_infall']['parcel_volume_shrinks'] is True
        assert row['plus_infall']['is_elongated_filament'] is True
        assert row['minus_child_drop']['aspect_ratio']==1.
        assert row['cost_comparison']['unified_minimum_cost']<row['cost_comparison']['stay_all_in_plus_cost']
        assert row['cost_comparison']['unified_minimum_cost']<row['cost_comparison']['cross_all_to_minus_cost']
    deep=unified_least_resistance_cycle(4.)
    assert deep['plus_infall']['filament_aspect_ratio']>79.
    assert deep['derived_cooker_takeover']['minus_volume_fraction']>0.9994


def test_nested_zero_crossing_map_flips_orientation_and_approaches_virial_partition():
    from recursive_horizons.nsc_atlas_offset_generator import nested_zero_crossing_map
    res=nested_zero_crossing_map(4,3.)
    signs=[row['orientation_sign'] for row in res['nested_chain']]
    assert signs==[-1,1,-1,1]
    row=res['nested_chain'][0]
    assert row['dark_energy_like_fraction']+row['dark_matter_like_fraction']+row['visible_filament_like_fraction']==pytest.approx(1.,rel=2e-15)
    assert .70<row['dark_energy_like_fraction']<.72
    assert .23<row['dark_matter_like_fraction']<.25
    assert .945<row['dark_sector_combined_fraction']<.955
    assert .045<row['visible_filament_like_fraction']<.055
    lock=res['two_valve_steady_state_lock']
    assert lock['three_dimensional_e_fold_dark_horizon']==pytest.approx(1-math.exp(-3),rel=2e-15)
    assert lock['three_dimensional_e_fold_visible_floor']==pytest.approx(math.exp(-3),rel=2e-15)
    assert lock['virial_3_to_1_dark_energy_at_95_horizon']==pytest.approx(.75*(1-math.exp(-3)),rel=2e-15)
    assert lock['virial_3_to_1_dark_matter_at_95_horizon']==pytest.approx(.25*(1-math.exp(-3)),rel=2e-15)
    asym=nested_zero_crossing_map(1,8.)['nested_chain'][0]
    assert asym['dark_energy_like_fraction']==pytest.approx(.75,abs=5e-4)
    assert asym['dark_matter_like_fraction']==pytest.approx(.25,abs=5e-4)
    assert asym['visible_filament_like_fraction']==pytest.approx(0.,abs=5e-4)
    assert res['antimatter_and_baryon_excess_check']['matches_observed_eta_without_scale_separation'] is False


def test_spectral_band_shifts_with_compression_and_reproduces_baryon_eta_inside_band():
    from recursive_horizons.nsc_atlas_offset_generator import (
        spectral_band_and_baryon_harmonic,branching_cosmos_chain_evolution)
    parent=spectral_band_and_baryon_harmonic(0.)
    child=spectral_band_and_baryon_harmonic(2.)
    assert child['absolute_zero_max_wavelength']==pytest.approx(parent['absolute_zero_max_wavelength']*math.exp(-2.),rel=2e-15)
    assert child['planck_cutoff_min_wavelength']==pytest.approx(parent['planck_cutoff_min_wavelength']*math.exp(-2.),rel=2e-15)
    assert child['frequency_scale_factor']==pytest.approx(math.exp(2.),rel=2e-15)
    assert child['within_room_spectrum'] is True
    assert child['eta_at_baryon_harmonic']==pytest.approx(6.1e-10,rel=1e-12)
    assert child['eta_at_absolute_zero_horizon']==pytest.approx(math.tanh(math.pi/2),rel=2e-15)
    chain=branching_cosmos_chain_evolution(generations=3,branch_factor=3,compression_per_step=2.5)
    assert chain['history'][0]['pattern_signature']=='oOOOooooooooo'
    for step_row in chain['history']:
        assert step_row['conservation_residual']==pytest.approx(0.,abs=2e-15)
        assert step_row['active_levels'][0]['state']=='o'
        assert step_row['active_levels'][-1]['state']=='O'


def test_cyclic_vortex_diffusion_ring_diffuses_and_respins_into_hurricane():
    from recursive_horizons.nsc_atlas_offset_generator import cyclic_vortex_diffusion_ring
    start=cyclic_vortex_diffusion_ring(0.)
    half=cyclic_vortex_diffusion_ring(.5)
    full=cyclic_vortex_diffusion_ring(1.)
    assert start['node_0_is_giant_hurricane'] is True
    assert half['node_0_is_giant_hurricane'] is False
    assert full['node_0_is_giant_hurricane'] is True
    assert start['node_states'].count('O')==start['node_states'].count('o')==4
    for row in (start,half,full):
        assert row['conservation_residual']==pytest.approx(0.,abs=2e-15)
    assert full['pure_diffusion_without_offset_envelope']<full['active_envelope']


def test_nonlocal_to_local_atlas_readout_solves_fixed_point_and_k7_harmonic():
    from recursive_horizons.nsc_atlas_offset_generator import nonlocal_to_local_atlas_readout
    atlas=nonlocal_to_local_atlas_readout()
    fp=atlas['two_valve_fixed_point']
    assert fp['visible_filament_fraction']==pytest.approx(math.exp(-3),rel=1e-14)
    assert fp['combined_dark_sector_fraction']==pytest.approx(1-math.exp(-3),rel=1e-14)
    assert fp['dark_energy_fraction']+fp['dark_matter_fraction']+fp['visible_filament_fraction']==pytest.approx(1.,rel=2e-15)
    assert 0.711<fp['dark_energy_fraction']<0.712
    assert 0.238<fp['dark_matter_fraction']<0.239
    ladder=atlas['discrete_3d_harmonic_asymmetry_ladder']
    assert ladder['spaghettification_factor']==pytest.approx(0.5,rel=2e-15)
    assert ladder['stage_K_7_eta']==pytest.approx(math.tanh(math.pi/4*math.exp(-21)),rel=2e-15)
    assert ladder['within_bbn_window'] is True
    assert len(atlas['translation_dictionary'])==4


def test_signed_circle_pendulum_duality_proves_minus_closure_and_cusp_free_zero_crossing():
    from recursive_horizons.nsc_atlas_offset_generator import signed_circle_pendulum_duality
    res=signed_circle_pendulum_duality(1.,2.*math.pi)
    closure=res['negative_half_required_for_closure']
    assert closure['missing_arc_fraction_without_minus']==pytest.approx(0.5,rel=1e-15)
    assert closure['closes_without_minus'] is False
    assert closure['closes_with_minus'] is True
    flips=res['minus_times_minus_equals_plus']
    assert flips['single_crossing_3d_inversion_det']==pytest.approx(-1.,rel=1e-15)
    assert flips['two_crossings_3d_inversion_det']==pytest.approx(1.,rel=1e-15)
    assert flips['matrix_identity_residual']==pytest.approx(0.,abs=1e-15)
    assert flips['parent_to_child_to_grandchild_signs']==['+','-','+']
    diag=res['singularity_diagnosis_at_zero']
    assert diag['signed_R_velocity_jump']==pytest.approx(0.,abs=1e-15)
    assert diag['signed_R_acceleration_at_zero']==pytest.approx(0.,abs=1e-15)
    assert diag['unsigned_size_r_velocity_cusp_jump']==pytest.approx(4.*math.pi,rel=1e-15)
    assert diag['signed_crossing_is_smooth'] is True
    assert diag['unsigned_size_creates_artificial_cusp_singularity'] is True
    for pt in res['swing_points']:
        assert pt['energy_conservation_residual']==pytest.approx(0.,abs=2e-14)


def test_vacuum_bubble_crystallization_surrounding_flow_and_mass_ratios():
    from recursive_horizons.nsc_atlas_offset_generator import vacuum_bubble_crystallization_and_mass_spectrum
    res=vacuum_bubble_crystallization_and_mass_spectrum()
    eq=res['surrounding_medium_equilibrium']
    assert eq['surrounding_flow_energy_share']==pytest.approx(0.75,rel=1e-15)
    assert eq['cavity_pressure_work_share']==pytest.approx(0.25,rel=1e-15)
    assert eq['virial_ratio_surrounding_over_pressure']==pytest.approx(3.0,rel=1e-15)
    assert eq['energy_in_point_singularity']==0.0
    assert eq['radial_shell_energy_fractions']['R_to_2R']==pytest.approx(0.5,rel=1e-15)
    ann=res['annihilation_to_back_to_back_ripples']
    assert ann['combined_winding']==0
    assert ann['outgoing_ripple_energies_over_E1']==[1.,1.]
    assert ann['separation_angle_degrees']==180.
    res_int=res['shared_flow_resonance_interaction']
    assert res_int['off_resonance_time_average']==0.
    assert res_int['in_resonance_aligned_energy_at_d_1']==-2.
    burst=res['discrete_burst_spectrum']
    assert burst['emissions_are_sharp_lines_not_smear'] is True
    assert burst['transition_2_to_1_energy_over_E1']==pytest.approx(2.**1.5-1.,rel=1e-15)
    masses=res['untuned_mass_ratio_tests']
    cand_b=masses['candidate_B_dimensional_accumulation_1D_to_3D']
    assert cand_b['bare_6_pi_5']==pytest.approx(6.*math.pi**5,rel=1e-15)
    assert cand_b['bare_ppm_error']<20.0
    assert cand_b['dressed_ppm_error']<1.0
    np_gap=masses['neutron_proton_gap_comparisons']
    assert abs(np_gap['threshold_capture_dressed_relative_error'])<0.002


def test_thermodynamic_gravity_and_cmb_temperature_bridges_planck_and_horizon():
    from recursive_horizons.nsc_atlas_offset_generator import thermodynamic_gravity_and_cmb_temperature
    res=thermodynamic_gravity_and_cmb_temperature(67.4)
    assert res['room_dynamic_range_R_H_over_l_P']>8e60
    assert 19.0<res['bare_geometric_mean_sqrt_TP_TH_K']<20.0
    assert 31.0<res['critical_blackbody_temperature_K']<32.0
    assert res['horizon_2d_efold_relative_error']<0.04
    assert res['two_valve_fixed_point_relative_error']<0.002


def test_derived_three_accumulator_bubble_ladder_selects_n3_for_3d_volume():
    from recursive_horizons.nsc_atlas_offset_generator import derived_three_accumulator_bubble_ladder, nonlocal_to_local_atlas_readout
    res=derived_three_accumulator_bubble_ladder(6)
    assert res['n1_is_1d_loop'] is True
    assert res['n2_is_2d_sheet'] is True
    assert res['minimum_n_for_3d_volume']==3
    assert res['untuned_n3_to_n1_energy_ratio']==pytest.approx(6.0,rel=1e-15)
    assert res['hits_3d_threshold_at_n_equals_3'] is True
    assert res['hits_1836_without_extra_factor'] is False
    assert res['ladder'][2]['volumetric_factor_sqrt_W1_W2_W3']==pytest.approx(3.0,rel=1e-15)
    assert res['ladder'][3]['volumetric_factor_sqrt_W1_W2_W3']==pytest.approx(math.sqrt(96.0),rel=1e-15)
    atlas=nonlocal_to_local_atlas_readout()
    assert 'derived_theorems_from_axioms' in atlas
    assert 'numerical_coincidences_not_derived' in atlas


def test_exponential_wall_bypasses_and_triad_winner_proves_111_ground_and_hard_core():
    from recursive_horizons.nsc_atlas_offset_generator import exponential_wall_bypasses_and_triad_winner
    res=exponential_wall_bypasses_and_triad_winner()
    for row in res['compression_wall_bypasses']:
        if row['chi']>0:
            assert row['bypass_1d_to_helix']<row['wall_1d_straight']
            assert row['bypass_3d_cross_zero_to_minus']<row['wall_3d_stay_in_plus']
    assert res['unique_3d_ground_state']==[1,1,1]
    assert res['ground_state_energy']==pytest.approx(4.0,rel=1e-15)
    assert res['ground_state_cost_per_winding']==pytest.approx(4.0/3.0,rel=1e-15)
    assert res['first_excited_energy']==pytest.approx(4.0*math.sqrt(2.0),rel=1e-15)
    assert res['short_range_hard_core_repulsion_barrier']==pytest.approx(8.0*(math.sqrt(2.0)-1.0),rel=1e-15)
    assert res['short_range_hard_core_repulsion_barrier']>3.31
    assert res['nuclear_potential_shape_derived'] is True


@pytest.mark.parametrize('winding,load',[(1.,1.),(-2.,3.),(.7,4.)])
def test_receiving_shape_spectrum_matches_independent_cost_second_variation(winding,load):
    from recursive_horizons.nsc_atlas_offset_generator import receiving_shape_spectrum
    row=receiving_shape_spectrum(winding,enclosing_load=load)
    y0=np.full(3,math.log(row['equilibrium_radius']))

    def cost(y):
        return winding**2*np.exp(-y).sum()+load*np.exp(y.sum())

    # Differentiate the original cost, rather than rebuilding its reported Hessian.
    h=2e-4
    e0=cost(y0)
    eye=np.eye(3)
    hessian=np.empty((3,3))
    for i in range(3):
        assert (cost(y0+h*eye[i])-cost(y0-h*eye[i]))/(2*h)==pytest.approx(0.,abs=e0*1e-8)
        hessian[i,i]=(cost(y0+h*eye[i])-2*e0+cost(y0-h*eye[i]))/h**2
        for j in range(i):
            hessian[i,j]=hessian[j,i]=(cost(y0+h*(eye[i]+eye[j]))
                -cost(y0+h*(eye[i]-eye[j]))-cost(y0+h*(-eye[i]+eye[j]))
                +cost(y0-h*(eye[i]+eye[j])))/(4*h*h)
    np.testing.assert_allclose(row['stiffness_matrix'],hessian,rtol=8e-8,atol=e0*1e-8)
    np.testing.assert_allclose(row['frequency_squared'],np.linalg.eigvalsh(hessian),rtol=1e-7)
    # Two independent trace-free distortions are degenerate; uniform dilation is separate.
    shear1=np.array([1.,-1.,0.])/math.sqrt(2)
    shear2=np.array([1.,1.,-2.])/math.sqrt(6)
    breathing=np.ones(3)/math.sqrt(3)
    shear_rate=shear1@hessian@shear1
    assert shear2@hessian@shear2==pytest.approx(shear_rate,rel=1e-7)
    assert breathing@hessian@breathing==pytest.approx(4*shear_rate,rel=1e-7)
    assert np.all(np.array(row['frequency_squared'])>0)


def test_receiving_mode_ratio_survives_winding_reversal_and_scale_changes():
    from recursive_horizons.nsc_atlas_offset_generator import receiving_shape_spectrum
    base=receiving_shape_spectrum(1.,1.)
    reversed_winding=receiving_shape_spectrum(-1.,1.)
    np.testing.assert_allclose(reversed_winding['stiffness_matrix'],base['stiffness_matrix'],rtol=1e-14)
    for winding,load,frequency_squared_factor in ((2.,1.,2**1.5),(1.,16.,2.)):
        row=receiving_shape_spectrum(winding,load)
        np.testing.assert_allclose(row['frequency_squared'],np.array(base['frequency_squared'])*frequency_squared_factor,rtol=2e-14)
        assert row['breathing_to_shear_ratio']==pytest.approx(2.,rel=2e-14)


def test_generalized_five_variable_minimum_is_stationary_in_original_cost():
    from recursive_horizons.nsc_atlas_offset_generator import unified_least_resistance_cycle
    for chi,moment,load in ((0.,4.,1.),(1.,4.,1.),(1.,1.,4.),(2.,1.,4.)):
        row=unified_least_resistance_cycle(chi,winding_budget=moment,enclosing_load=load)
        y=np.array([math.log(row['plus_infall']['optimal_radial_scale']),
                    *np.log(row['minus_child_drop']['optimal_principal_scales']),
                    row['transferred_winding_fraction_theta']])

        def cost(z):
            theta=z[-1]
            return (moment*(1-theta)**2*(2*math.exp(chi)+math.exp(-z[0]))
                    +load*math.exp(z[0]-2*chi)
                    +moment*theta**2*np.exp(-z[1:4]).sum()+load*math.exp(z[1:4].sum()))

        original_cost=cost(y)
        assert row['cost_comparison']['unified_minimum_cost']==pytest.approx(original_cost,rel=2e-13)
        for direction in np.eye(5):
            derivative=(cost(y+1e-5*direction)-cost(y-1e-5*direction))/(2e-5)
            assert derivative==pytest.approx(0.,abs=original_cost*2e-8)


def test_free_compression_minimum_is_not_the_imposed_visible_fraction_root():
    from recursive_horizons.nsc_atlas_offset_generator import unified_least_resistance_cycle
    theta=s.symbols('theta',positive=True)
    # Free scale minimization in BOTH regions leaves this one strictly convex function.
    reduced=4*((1-theta)**s.Rational(3,2)+theta**s.Rational(3,2))
    assert s.simplify(s.diff(reduced,theta).subs(theta,s.Rational(1,2)))==0
    assert s.simplify(s.diff(reduced,theta,2))==3/s.sqrt(theta)+3/s.sqrt(1-theta)
    chi=.5*math.log(2)
    row=unified_least_resistance_cycle(chi)
    assert row['transferred_winding_fraction_theta']==pytest.approx(.5,rel=2e-15)
    assert row['cost_comparison']['unified_minimum_cost']==pytest.approx(2*math.sqrt(2),rel=2e-15)
    # The envelope's outward derivative remains positive beyond this minimum.
    for later_chi in (1.,2.,3.):
        later=unified_least_resistance_cycle(later_chi)
        retained=1-later['transferred_winding_fraction_theta']
        derivative=2*retained**2*math.exp(later_chi)-2*retained*math.exp(-later_chi)
        assert derivative>0
        assert later['cost_comparison']['unified_minimum_cost']>row['cost_comparison']['unified_minimum_cost']


def test_free_parent_causal_length_matches_independent_optical_integral():
    from scipy.integrate import quad
    from recursive_horizons.nsc_atlas_offset_generator import one_room_geometric_continuation
    a=math.sqrt(4*math.pi**2/6)
    for offset in (0.,.25,.5):
        row=one_room_geometric_continuation(offset)
        # Integrate the actual metric's null speed, without using its tail formula.
        optical_tail=quad(lambda u:1/math.cosh(a*u)**2,offset,10.,epsabs=1e-13)[0]
        assert row['observer_causal_horizon_proper_length']==pytest.approx(
            row['stretch']*optical_tail,rel=2e-12)
        assert row['phase_energy']+row['geometric_kinetic_energy']==pytest.approx(4*math.pi**2,rel=3e-15)
    large=one_room_geometric_continuation(1000.)
    assert large['stretch'] is None  # use the retained logarithm beyond float exponent range
    assert math.isfinite(large['log_stretch'])
    assert large['observer_causal_horizon_proper_length']==pytest.approx(1/(2*a),rel=2e-15)


def test_characteristic_fork_initial_transmission_matches_contact_constraints():
    from recursive_horizons.nsc_atlas_offset_generator import parent_receiver_characteristics
    offset=1e-6
    row=parent_receiver_characteristics(parents=(-1,0,0),duration=offset,optical_step=.01,sample_count=2)['rows'][-1]
    # Initial root incoming power is pi^2. This probes the implemented junction,
    # independently of a separately constructed scattering matrix.
    received=np.array(row['gross_received_through_contact'][1:])/(math.pi**2*offset)
    np.testing.assert_allclose(received,[4/9,4/9],rtol=2e-10)
    assert 1-received.sum()==pytest.approx(1/9,abs=2e-10)


@pytest.mark.parametrize('parents',[(-1,0,1),(-1,0,0)])
def test_parent_receiver_forces_and_regional_power_follow_one_hamiltonian(parents):
    from recursive_horizons.nsc_atlas_offset_generator import _parent_receiver_system
    system=_parent_receiver_system(parents,6)
    state=system['initial'].copy()
    nf=system['free_phase_count'];nr=system['room_count']
    # Exercise junction inertia and moving, unequal regions away from the supplied start.
    state[:nf]+=.07*np.sin(np.arange(nf))
    state[nf:2*nf]=np.linspace(-.12,.17,nf)
    state[2*nf:2*nf+nr]+=[.01,-.005,.003]
    state[2*nf+nr:2*nf+2*nr]=[.07,-.03,.02]
    rhs,values=system['evaluate'](state)
    h=1e-6
    for indices,expected in ((range(nf),-rhs[nf:2*nf]),
                             (range(nf,2*nf),rhs[:nf]),
                             (range(2*nf,2*nf+nr),-values['geometry_force'])):
        derivatives=[]
        for index in indices:
            delta=np.zeros_like(state);delta[index]=h
            derivatives.append((system['evaluate'](state+delta)[1]['total_energy']
                               -system['evaluate'](state-delta)[1]['total_energy'])/(2*h))
        np.testing.assert_allclose(derivatives,expected,rtol=2e-7,atol=2e-8)
    energy_rate=(system['evaluate'](state+h*rhs)[1]['phase_energy']
                -system['evaluate'](state-h*rhs)[1]['phase_energy'])/(2*h)
    np.testing.assert_allclose(energy_rate,values['transfer']+values['work'],rtol=2e-7,atol=2e-8)
    assert values['transfer'].sum()==pytest.approx(0.,abs=2e-13)
    geometric_power=np.dot(state[2*nf+nr:2*nf+2*nr]/system['inertia'],values['geometry_force'])
    assert geometric_power+values['work'].sum()==pytest.approx(0.,abs=2e-13)


def test_parent_receiver_zero_phase_has_no_hidden_pressure_or_motion():
    from recursive_horizons.nsc_atlas_offset_generator import _parent_receiver_system
    system=_parent_receiver_system((-1,0,1),6)
    state=system['initial'].copy()
    state[:system['free_phase_count']]=0
    rhs,values=system['evaluate'](state)
    np.testing.assert_array_equal(rhs,np.zeros_like(rhs))
    assert values['total_energy']==0
    np.testing.assert_array_equal(values['pressure'],np.zeros(3))
    np.testing.assert_array_equal(values['transfer'],np.zeros(3))


def test_parent_receiver_characteristics_recovers_exact_free_one_room_solution():
    from recursive_horizons.nsc_atlas_offset_generator import parent_receiver_characteristics
    q0=4*math.pi**2
    a=math.sqrt(q0/6)
    results=[parent_receiver_characteristics(parents=(-1,),duration=.5,optical_step=h,sample_count=3)
             for h in (.02,.01)]
    for result in results:
        assert result['reached_offset']==pytest.approx(.5,abs=1e-15)
        assert result['initial_total_budget']==pytest.approx(q0,rel=2e-15)
        assert result['positive_volume_guard_reached'] is False
    fine=results[1]
    for row in fine['rows']:
        t=row['offset'];lam=math.cosh(a*t)**2
        # The returned intermediate rows use linear history interpolation.
        assert row['outer_radius'][0]==pytest.approx(lam,rel=5e-6)
        assert row['outer_log_radius_velocity'][0]==pytest.approx(2*a*math.tanh(a*t),rel=5e-6,abs=1e-14)
        assert row['continuum_optical_distance'][0]==pytest.approx(math.tanh(a*t)/a,rel=5e-6,abs=1e-14)
        assert row['regional_phase_energy'][0]==pytest.approx(q0/lam,rel=5e-6)
    endpoint=fine['rows'][-1]
    assert endpoint['outer_radius'][0]==pytest.approx(math.cosh(a*.5)**2,rel=1e-9)
    assert endpoint['total_energy']==pytest.approx(q0,rel=1e-9)
    assert fine['max_accepted_and_sampled_total_budget_residual']<results[0]['max_accepted_and_sampled_total_budget_residual']
    assert fine['max_accepted_and_sampled_regional_ledger_residual']<results[0]['max_accepted_and_sampled_regional_ledger_residual']


def test_characteristic_fork_preserves_equal_receivers_and_contact_balance():
    from recursive_horizons.nsc_atlas_offset_generator import parent_receiver_characteristics
    result=parent_receiver_characteristics(parents=(-1,0,0),duration=.3,optical_step=.02,sample_count=3)
    assert result['reached_offset']==pytest.approx(.3,abs=1e-15)
    assert result['max_accepted_and_sampled_contact_cancellation']<2e-12
    assert result['max_accepted_and_sampled_interface_phase_gap']<2e-12
    for row in result['rows']:
        for key in ('regional_volume','regional_phase_energy','outer_radius',
                    'net_received_through_contact','field_geometry_work'):
            assert row[key][1]==pytest.approx(row[key][2],rel=2e-13,abs=2e-13)
        assert row['net_received_through_contact'][0]+sum(row['net_received_through_contact'][1:])==pytest.approx(0.,abs=2e-12)
    assert result['rows'][-1]['regional_phase_energy'][1]>0


def test_inherited_four_center_geometry_preserves_first_split_and_opens_transversely():
    from recursive_horizons.nsc_atlas_offset_generator import (
        _four_center_volume,_neck_volume_factor,inherited_daughter_calculation)

    result=inherited_daughter_calculation()
    separation=result['parent_half_separation']
    assert _four_center_volume(separation,0.,64)==pytest.approx(_neck_volume_factor(separation),abs=2e-12)
    assert result['inherited_shape'][-1]['axial_volume_derivative']>0
    # Differentiate the four-source volume using whole upper-lobe rays,
    # independently of the coarea derivative and Voronoi quadrature.
    mu,weight=np.polynomial.legendre.leggauss(48)
    phi=2*np.pi*(np.arange(96)+.5)/96
    nx=np.sqrt(1-mu[:,None]**2)*np.cos(phi)[None,:]
    def whole_lobe_volume(half_split):
        low=np.full_like(nx,.1);high=np.full_like(nx,separation)
        for _ in range(48):
            radius=(low+high)/2
            value=np.zeros_like(radius)
            for dx,dz in ((half_split,0.),(-half_split,0.),
                          (half_split,2*separation),(-half_split,2*separation)):
                value+=1/np.sqrt(radius**2+dx*dx+dz*dz-2*dx*radius*nx+2*dz*radius*mu[:,None])
            low=np.where(value>4,radius,low)
            high=np.where(value>4,high,radius)
        return np.dot(weight,np.mean(((low+high)/2)**3,axis=1))
    derivative=(whole_lobe_volume(.005)-whole_lobe_volume(0.))/.005**2
    expected=result['inherited_shape'][-1]['transverse_volume_derivative']
    assert derivative==pytest.approx(expected,rel=.002)
    b=result['second_neck_half_separation'];z=result['second_neck_axial_position']
    potential=sum(1/math.sqrt(b*b+(z-sign*separation)**2) for sign in (-1,1) for _ in range(2))
    assert potential==pytest.approx(4.,abs=1e-11)
    force=sum((z-sign*separation)/(b*b+(z-sign*separation)**2)**1.5 for sign in (-1,1))
    assert force==pytest.approx(0.,abs=1e-11)
    state=result['inherited_from_first_split']
    assert state['shape_velocity']>0
    assert any(any(abs(x)>0 for x in mode['real']) for mode in state['outgoing_memory'])
    c=math.sqrt(state['shape_coordinate'])
    speed=state['shape_velocity']/(2*c)
    point=np.array([5.,1.,2.]);a=math.sqrt(2)
    def field_jet(point, center, velocity, charge):
        phase=0.;rate=0.
        for sign in (-1,1):
            delta=point-np.array([0.,0.,sign*center])
            distance=np.linalg.norm(delta)
            phase+=charge/distance
            rate+=charge*sign*velocity*delta[2]/distance**3
        return phase,rate
    original=field_jet(point/a,c,speed,1.)
    scaled=field_jet(point,a*c,speed,2.)
    assert scaled[0]==pytest.approx(a*original[0],rel=2e-14)
    assert scaled[1]==pytest.approx(original[1],rel=2e-14)


def test_whole_exterior_neck_onset_obeys_wave_equation_and_shape_boundary():
    from recursive_horizons.nsc_atlas_offset_generator import split_onset_field,split_onset_calculation

    radius,rate=s.symbols('radius rate',positive=True)
    radial=s.exp(-rate*(radius-1))*(rate**2*radius**2+3*rate*radius+3)/radius**3
    assert s.simplify(rate**2*radial-s.diff(radial,radius,2)
                      -2*s.diff(radial,radius)/radius+6*radial/radius**2)==0
    at_boundary=s.simplify((s.diff(radial,radius)+4*radial).subs(radius,1))
    assert s.expand(at_boundary)==-rate**3+3*rate+3
    result=split_onset_calculation()
    assert result['growth_rate']**3-3*result['growth_rate']-3==pytest.approx(0.,abs=5e-14)
    assert max(abs(row['energy_residual']) for row in result['rows'])<1e-15
    quiet=split_onset_field(.4,1.5)
    assert quiet['time_rate']==0.
    arrived=split_onset_field(.8,1.5)
    assert arrived['wavefront_reached']
    assert arrived['time_rate']>0
    endpoint=split_onset_field(1.,1.)
    assert endpoint['radial_rate']+4*endpoint['phase']==pytest.approx(0.,abs=1e-15)
    far=split_onset_field(10.,10.)
    assert np.isfinite(far['phase'])


def test_neck_motion_generates_causal_wave_with_independent_energy_balance():
    from recursive_horizons.nsc_atlas_offset_generator import split_waveform_evolution

    result = split_waveform_evolution()
    radius = result['matching_radius']
    assert result['zero_shape_inertia'] == pytest.approx(8/15*(1-radius**-3), rel=2e-12)
    assert result['pinch']['time'] < result['endpoint']['time']
    assert result['pinch']['coordinate'] > 0
    assert result['endpoint']['separation'] == pytest.approx(1.2, abs=1e-8)
    assert result['maximum_budget_residual'] < 1e-6
    # Reconstructed fields integrated in space independently check boundary work.
    assert abs(result['independent_exterior_energy_gap']) < 2e-8
    assert result['detector']['integrated_outward_energy'] > 0
    assert result['detector']['polar_phase_departure'][-1] > 0
    assert result['detector']['equatorial_phase_departure'][-1] < 0
    assert result['detector']['arrival_offsets'][0] == pytest.approx(10-radius)


def test_hybrid_echo_cutoffs_and_late_redshift():
    from recursive_horizons.nsc_atlas_offset_generator import hybrid_receiving_echo

    result = hybrid_receiving_echo(samples=65, tolerance=2e-9, endpoint=16.)
    assert result['max_budget_defect'] < 1e-7
    assert .11 < result['detector_cutoff'] < .12
    assert .22 < result['mouth_cutoff'] < .23
    assert result['mouth_cutoff'] < result['horizon_enters_receiver_at']
    for region, key, arrival in ((0, 'mouth_echoes', 'mouth_return_offset'),
                                 (1, 'detector_echoes', 'detection_offset')):
        before, after = result[key][-2:]
        rate = result['asymptotic_rate_bracket'][region][0]
        assert after[arrival] - before[arrival] == pytest.approx(math.log(10)/rate, rel=.002)
        assert after['energy_ratio']/after['distance_to_cutoff'] == pytest.approx(rate, rel=.001)


def test_characteristic_chain_has_exactly_zero_field_before_grandchild_front():
    from recursive_horizons.nsc_atlas_offset_generator import parent_receiver_characteristics
    result=parent_receiver_characteristics(parents=(-1,0,1),duration=.3,optical_step=.02,sample_count=3)
    assert result['reached_offset']==pytest.approx(.3,abs=1e-15)
    assert result['rows'][-1]['continuum_optical_distance'][1]<2
    for row in result['rows']:
        assert row['regional_phase_energy'][2]==0
        assert row['net_received_through_contact'][2]==0
        assert row['field_geometry_work'][2]==0
    assert result['rows'][-1]['regional_phase_energy'][1]>0
    assert result['new_child_formation_proved'] is False
    assert result['one_way_horizon_established'] is False


def test_three_axis_orthogonal_cascade_and_joint_kinetic_handoff():
    from recursive_horizons.nsc_atlas_offset_generator import inherited_daughter_calculation

    result = inherited_daughter_calculation()
    cascade = result['three_axis_orthogonal_cascade']
    assert cascade['first_to_second_split_ratio'] == pytest.approx(-.5, abs=1e-15)
    assert cascade['second_to_third_symmetric_ratio'] == pytest.approx(-2., rel=1e-12)
    for row in cascade['lobe_directional_scan']:
        assert abs(row['harmonic_trace_residual']) < 1e-14
        assert row['y_is_unique_global_minimum']
        assert row['volume_derivatives']['y'] < 0.
        if row['daughter_half_separation'] >= 1.:
            assert row['both_in_plane_axes_uphill']
    sym = cascade['lobe_directional_scan'][3]
    assert sym['daughter_half_separation'] == pytest.approx(sym['parent_half_separation'])
    assert sym['volume_derivatives']['x'] == pytest.approx(sym['volume_derivatives']['z'], rel=1e-12)
    assert sym['volume_derivatives']['y'] == pytest.approx(-2. * sym['volume_derivatives']['x'], rel=1e-12)

    pinch1 = cascade['third_neck_pinch_at_b_1']
    assert pinch1['third_neck_half_separation_y'] == pytest.approx(.360449324968, rel=1e-9)
    assert abs(pinch1['potential_residual']) < 1e-11
    pinch_sym = cascade['third_neck_pinch_at_b_equal_s']
    assert pinch_sym['third_neck_half_separation_y'] == pytest.approx(.347227979949, rel=1e-9)
    assert abs(pinch_sym['potential_residual']) < 1e-11

    handoff = result['joint_two_coordinate_kinetic_handoff']
    matrices = handoff['quasistatic_full_matrices_by_matching_radius']
    ref_mat = np.array(matrices[0]['matrix'])
    for entry in matrices[1:]:
        assert np.allclose(entry['matrix'], ref_mat, rtol=1e-9, atol=1e-9)
    for mode in handoff['multipole_mode_ratios']:
        l = mode['degree']
        assert mode['ratio'] == pytest.approx(-.5 * (l - 1), abs=1e-15)
    assert handoff['projected_determinant'] > .1
    assert handoff['initial_transverse_conjugate_momentum_p2'] > .05
    assert handoff['transverse_potential_slope_dU_dq2'] < 0.
    assert handoff['transverse_wave_memory_reaction_R2'] < 0.
    assert handoff['net_transverse_potential_plus_memory_drive'] > .18


def test_exact_budget_curve_invariant_and_zero_parameter_inverse():
    from recursive_horizons.nsc_atlas_offset_generator import exact_budget_curve_invariant

    result = exact_budget_curve_invariant()
    assert result['maximum_cubic_residual'] < 1e-15
    assert result['maximum_inverse_chi_error'] < 1e-8
    for row in result['sampled_curve']:
        if row['chi'] <= 4.:
            assert row['chi_recovery_error'] < 1e-12
    assert result['sampled_curve'][-1]['B_over_V_squared'] == pytest.approx(16. / 27., rel=1e-4)
    for pred in result['zero_parameter_predictions_from_V'].values():
        assert abs(pred['cubic_identity_residual']) < 1e-15
        assert pred['derived_D'] + pred['derived_M'] + pred['input_visible_share_V'] == pytest.approx(1., abs=1e-15)
        assert pred['derived_chi'] >= 3. - 1e-12
        assert pred['tangent_dD_dV'] + pred['tangent_dM_dV'] == pytest.approx(-1., abs=1e-15)
        assert -0.82 < pred['tangent_dD_dV'] < -0.75


def test_covariant_action_to_a0_and_friedmann_eos():
    from recursive_horizons.nsc_atlas_offset_generator import covariant_action_to_a0_and_friedmann_eos

    res = covariant_action_to_a0_and_friedmann_eos()
    assert res['covariant_to_a0_reduction']['reduction_is_exact_on_harmonic_curve_bundles'] is True
    eos = res['friedmann_equations_of_state']
    assert eos['closed_valve_partial_dE_dV'] == pytest.approx(-1., abs=1e-15)
    assert eos['open_valve_turn_injection_dE_dV'] == pytest.approx(4., abs=1e-15)
    assert eos['open_valve_total_dE_dV'] == pytest.approx(3., abs=1e-15)
    assert eos['w_DE'] == pytest.approx(-1., abs=1e-15)
    assert eos['w_DE_equals_minus_one'] is True
    assert eos['w_DM'] == pytest.approx(0., abs=1e-15)
    assert eos['w_DM_equals_zero'] is True
    assert eos['baryon_triad_w'] == pytest.approx(0., abs=1e-15)
    assert eos['baryon_triad_w_equals_zero'] is True
    reg = res['curvature_regularization_3plus1d']
    assert reg['kinematic_cusp_removed_by_signed_R'] is True
    assert reg['standing_wave_hard_core_radius_R_core'] == pytest.approx(math.sqrt(2. * math.pi), rel=1e-15)
    assert reg['kretschmann_static_at_R_zero_12_over_R_core_4'] == pytest.approx(3. / (math.pi ** 2), rel=1e-15)
    assert reg['kretschmann_dynamic_at_R_zero'] == pytest.approx(12. / (math.pi ** 2), rel=1e-15)
    assert reg['proper_birth_distance_finite'] is True
    vir = res['gauss_divergence_exterior_virial']
    assert vir['bulk_to_throat_virial_ratio'] == pytest.approx(3., abs=1e-15)
    assert vir['bulk_share_of_exterior_virial'] == pytest.approx(.75, abs=1e-15)
    assert vir['throat_share_of_exterior_virial'] == pytest.approx(.25, abs=1e-15)
    card = res['cardano_golden_ratio_neck_onset']
    assert abs(card['cubic_cardano_residual'] if 'cubic_cardano_residual' in card else card['cubic_residual']) < 1e-14
    assert card['p_plus'] == pytest.approx(2.1038034027355366, rel=1e-14)
    freeze = res['three_axis_cascade_freezeout']
    assert freeze['cumulative_freezeout_chi'] == pytest.approx(3., abs=1e-15)



