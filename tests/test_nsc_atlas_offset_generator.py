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
