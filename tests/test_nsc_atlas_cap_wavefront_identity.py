"""Exact internal wavefront laws and finite-interval geometric evolution.

Right-minus-left jumps; no numerical time steps or retained field samples.
The two-front solution is a manufactured alpha=0 control, not a finite-alpha
solution or a new physical preparation for the existing cap-wall candidate.
"""
import pytest
import sympy as s

from control_budget import ControlBudget


CONTROLS=ControlBudget(maximum_evaluations=50_000,cpu_seconds=20.)


@pytest.fixture(scope='module',autouse=True)
def budget():
    with CONTROLS.active() as value:
        value.check(18_000)
        yield
    print('internal wavefront identities:',CONTROLS.budget.evaluations,
          'charged algebra values; zero field samples;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def test_nonlinear_mass_interface_term_pairs_the_two_scalar_jump_residuals():
    r,A,B,v,tx,P=s.symbols('r A B v tx P',positive=True)
    pl,pr,hl,hr=s.symbols('pl pr hl hr',real=True)
    rho=lambda p,h:r*r*(P/2+A*(p*p+h*h)/4)
    flux=lambda p,h:r*r*A*B*p*h/2
    ep=v*(pr-pl)+tx*B*(hr-hl)
    eh=v*(hr-hl)+tx*B*(pr-pl)
    direct=v*(rho(pr,hr)-rho(pl,hl))+tx*(flux(pr,hr)-flux(pl,hl))
    pairing=A*r*r*((pl+pr)*ep+(hl+hr)*eh)/4
    assert s.expand(direct-pairing)==0


def test_nonzero_jumps_require_the_actual_characteristic_speed():
    B,v,tx,j=s.symbols('B v tx j',nonzero=True,real=True)
    matrix=s.Matrix([[v,tx*B],[tx*B,v]])
    assert s.expand(matrix.det())==v*v-tx*tx*B*B
    for sign in (-1,1):
        assert matrix.subs(v,sign*tx*B)*s.Matrix([-sign*j,j])==s.zeros(2,1)
    # A background front moving at +/-A0 is not a finite-B characteristic
    # just because its polarization matches that background front.
    A0=s.symbols('A0',positive=True)
    defect=matrix.subs(v,tx*A0)*s.Matrix([-j,j])
    assert s.expand(defect[0]-tx*j*(B-A0))==0
    assert s.expand(defect[1]-tx*j*(A0-B))==0


def test_zero_mass_interface_work_does_not_imply_scalar_compatibility():
    u,B,r,A=s.symbols('u B r A',positive=True)
    pl,pr,hl,hr=-u,u,0,0
    # A stationary discontinuity of this type conserves the quadratic
    # density across the face but violates the Phi equation.
    rho_jump=A*r*r*(pr*pr+hr*hr-pl*pl-hl*hl)/4
    flux_jump=A*B*r*r*(pr*hr-pl*hl)/2
    assert rho_jump==flux_jump==0
    assert B*(pr-pl)==2*B*u


def test_ramp_evolves_as_two_compatible_fronts_and_preserves_full_cap_mass():
    z,t,ell,Z=s.symbols('z t ell Z',positive=True)
    q,beta=s.symbols('q beta',real=True)
    # 0<t<min(ell,Z-ell): no axis or wall interaction. These are the
    # three open pieces of 3q[(z-ell-t)_+ + (z-ell+t)_+].
    fronts=(ell-t,ell+t)
    chi=(s.Integer(0),3*q*(z-ell+t),6*q*(z-ell))
    k=tuple(s.diff(c,t) for c in chi)
    slope=tuple(s.diff(c,z) for c in chi)
    density=tuple((p*p+h*h)/4 for p,h in zip(k,slope))
    flux=tuple(p*h/2 for p,h in zip(k,slope))
    for c,p,h in zip(chi,k,slope):
        assert s.diff(c,t)-p==0
        assert s.diff(p,t)-s.diff(h,z)==0
        assert s.diff(h,t)-s.diff(p,z)==0
    for i,front in enumerate(fronts):
        velocity=s.diff(front,t)
        assert s.expand((chi[i+1]-chi[i]).subs(z,front))==0
        assert velocity*(k[i+1]-k[i])+slope[i+1]-slope[i]==0
        assert velocity*(slope[i+1]-slope[i])+k[i+1]-k[i]==0
        assert s.expand(velocity*(density[i+1]-density[i])+flux[i+1]-flux[i])==0
    edges=(0,*fronts,Z)
    wall=chi[-1].subs(z,Z)
    mass=sum(s.integrate(density[i],(z,edges[i],edges[i+1])) for i in range(3))-beta*wall**2/4
    expected=9*q*q*(Z-ell)-9*beta*q*q*(Z-ell)**2
    assert s.expand(mass-expected)==0
    assert s.diff(mass,t)==0
    assert k[-1]*(slope[-1]-beta*wall)/2==0
    # The t->0 data are precisely 6q(z-ell)_+ with zero velocity.
    assert chi[0]==k[0]==k[-1]==0
    assert chi[-1]==6*q*(z-ell)


def test_same_solution_obeys_the_moving_wall_shape_and_incoming_trace():
    L,q,beta,A0,Zx=s.symbols('L q beta A0 Zx',real=True)
    mass=9*q*q*L-9*beta*q*q*L*L
    # L=Z-ell and beta_Z=-2A0-beta^2 for R=tanh(Z).
    shape=s.diff(mass,L)-(2*A0+beta*beta)*s.diff(mass,beta)
    wall_chi,wall_k,wall_slope=6*q*L,0,6*q
    density_z=A0*wall_chi**2/2+(wall_k**2+(wall_slope-beta*wall_chi)**2)/4
    assert s.expand(shape-density_z)==0
    assert s.expand(Zx*shape-Zx*density_z)==0
    incoming_drive=wall_k+wall_slope-beta*wall_chi
    assert incoming_drive==6*q-6*beta*q*L
    assert wall_k*(incoming_drive-wall_k)/2==0
