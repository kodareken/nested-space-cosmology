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
