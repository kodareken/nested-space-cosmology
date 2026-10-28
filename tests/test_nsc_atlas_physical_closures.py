"""Exact dimensional freedom and surface-state ownership in the used actions."""
import sympy as s


def test_scalar_flux_parameter_family_preserves_shape_but_not_absolute_scale():
    m,mu,q,phi,scale=s.symbols('m mu q phi scale',real=True)
    U=(m*m*phi*phi+(q+mu*phi)**2)/2
    transformed=U.subs({m:scale*m,mu:scale*mu,q:scale*q},simultaneous=True)
    assert s.expand(transformed-scale**2*U)==0
    assert s.simplify((m*m/(m*m+mu*mu)).subs({m:scale*m,mu:scale*mu},simultaneous=True)-m*m/(m*m+mu*mu))==0
    rv,E,C,A=s.symbols('rho_v E0 C A',positive=True)
    H2=rv/3+E/(3*A**3)-C/A**2
    changed=H2.subs({rv:scale**2*rv,E:scale**2*E,C:scale**2*C},simultaneous=True)
    assert s.expand(changed-scale**2*H2)==0


def test_surface_scalar_conservation_needs_its_own_inherited_charge():
    R,velocity,charge,Rdot=s.symbols('R velocity Q Rdot',positive=True)
    acceleration=-2*Rdot*velocity/R
    sigma=velocity**2/2
    sigma_dot=s.diff(sigma,velocity)*acceleration
    assert s.simplify(sigma_dot*R**4+4*sigma*R**3*Rdot)==0
    assert s.simplify(sigma.subs(velocity,charge/(4*s.pi*R**2))-charge**2/(32*s.pi**2*R**4))==0
