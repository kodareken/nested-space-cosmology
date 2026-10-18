"""Hilbert/time, covariant Ward, canonical variation and transmission controls."""
import math

import numpy as np
import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_evolution import ALPHA
from recursive_horizons.nsc_atlas_dirac_stress import (
    SIGMA1, canonical_metric_vertices, stress_matrices,
)
from recursive_horizons.nsc_atlas_collar_field import QuietCollar


U = np.array([[1+.3j, -.2+.4j], [.7-.1j, .4+.6j]])
UN = np.array([[.2+.5j, .6-.4j], [-.1+.8j, .9+.2j]])


def test_normal_time_Hilbert_energy_flux_and_angular_multiplet():
    # A point of the actual exterior normal chart, not a flat operator.
    g = QuietCollar(.5)
    r, rn, N, Nn, w, wn = g.at([.01])[:, 0]
    Ut = -ALPHA @ (N*UN+Nn*U/2)-1j*SIGMA1 @ U*N/r
    T = stress_matrices(U, UN, radial_length=1., radius=r)
    weight = 2/(4*math.pi)
    rho = weight*.5j*(U.conj().T @ Ut-Ut.conj().T @ U)/(N*r*r)
    flux = weight*.25j*(-U.conj().T @ UN+UN.conj().T @ U
        +(U.conj().T @ ALPHA @ Ut-Ut.conj().T @ ALPHA @ U)/N)/(r*r)
    assert T['energy_density'] == pytest.approx(rho, abs=2e-15)
    assert T['outward_energy_flux'] == pytest.approx(flux, abs=2e-15)
    assert -T['energy_density']+T['radial_pressure']+2*T['angular_pressure'] == pytest.approx(np.zeros((2,2)), abs=1e-15)
    # Explicit lowest angular multiplet has pointwise summed norm 2/(4*pi).
    theta, phi = .7, .4
    omega = np.array([[-math.cos(theta), -math.sin(theta)*np.exp(-1j*phi)],
                      [-math.sin(theta)*np.exp(1j*phi), math.cos(theta)]])/math.sqrt(4*math.pi)
    assert np.trace(omega.conj().T @ omega).real == pytest.approx(weight)
    assert T['charge_density'] == pytest.approx(weight*(U.conj().T @ U)/r**2)


def test_four_dimensional_covariant_conservation_from_on_shell_jets():
    # Exact local data; the field time jets follow the Dirac equation.
    # Contract all four coordinate indices against independently built Gamma.
    t, n = s.symbols('t n', real=True)
    I = s.I
    alpha, sigma = s.Matrix([[0,-I],[I,0]]), s.Matrix([[0,1],[1,0]])
    u = s.Matrix([[1+I/3, -s.Rational(1,5)+I/7], [s.Rational(2,3)-I/4, 1+I/2]])
    un = s.Matrix([[s.Rational(2,7)+I/5, 1-I/3], [-I/2, s.Rational(3,5)+I]])
    unn = s.Matrix([[I/7, s.Rational(1,3)], [1-I/2, I/4]])
    N0, Nn, Nnn, Nt = s.Rational(5,4), -s.Rational(2,3), s.Rational(7,8), s.Rational(9,10)
    r0, rn, rt = s.Rational(3,5), -s.Rational(1,4), s.Rational(2,7)
    W, Wn = N0/r0, Nn/r0-N0*rn/r0**2
    ut = -alpha*(N0*un+Nn*u/2)-I*sigma*W*u
    unt = -alpha*(N0*unn+3*Nn*un/2+Nnn*u/2)-I*sigma*(Wn*u+W*un)
    field = u+t*ut+n*un+t*n*unt+n*n*unn/2
    N, r = N0+Nt*t+Nn*n+Nnn*n*n/2, r0+rt*t+rn*n
    dn = field.diff(n)
    qa = -I*(field.H*alpha*dn-dn.H*alpha*field)/2
    qi = -I*(field.H*dn-dn.H*field)/2
    angular = field.H*sigma*field
    rho, pr, pt, j = (qa+angular/r)/r**2, qa/r**2, angular/(2*r**3), qi/r**2
    zero = s.zeros(2)
    tensor = [[zero for _ in range(4)] for _ in range(4)]
    tensor[0][0], tensor[0][1], tensor[1][0], tensor[1][1] = rho/N**2, j/N, j/N, pr
    tensor[2][2] = tensor[3][3] = pt/r**2
    metric = s.diag(-N*N, 1, r*r, r*r)  # theta=pi/2
    at0 = {t:0, n:0}
    g0 = metric.subs(at0)
    inverse = g0.inv()
    dg = [metric.diff(x).subs(at0) for x in (t,n)]+[s.zeros(4),s.zeros(4)]
    Gamma = [[[sum(inverse[i,k]*(dg[j][k,m]+dg[m][k,j]-dg[k][j,m])/2
                   for k in range(4)) for m in range(4)] for j in range(4)] for i in range(4)]
    T0 = [[entry.subs(at0) for entry in row] for row in tensor]
    for nu in range(4):
        divergence = tensor[0][nu].diff(t).subs(at0)+tensor[1][nu].diff(n).subs(at0)
        for mu in range(4):
            for lam in range(4):
                divergence += Gamma[mu][mu][lam]*T0[lam][nu]+Gamma[nu][mu][lam]*T0[mu][lam]
        assert divergence.applyfunc(s.expand) == s.zeros(2)


def test_canonical_metric_variation_and_nonzero_local_contacts():
    N, a, r, beta, qa, qi, angular = s.symbols('N a r beta qa qi angular', real=True)
    coordinates = (N,a,r,beta)
    h = 2*(N*(qa/a+angular/r)-beta*qi)
    g = dict(lapse=1.4, radial_length=1.2, radius=.6, shift=-.13)
    vertex = canonical_metric_vertices(U,UN,**g)
    T = stress_matrices(U,UN,radial_length=g['radial_length'],radius=g['radius'])
    assert vertex['first']['lapse'] == pytest.approx(4*math.pi*g['radial_length']*g['radius']**2*T['energy_density'])
    assert vertex['first']['radial_length'] == pytest.approx(-4*math.pi*g['lapse']*g['radius']**2*T['radial_pressure'])
    assert vertex['first']['radius'] == pytest.approx(-8*math.pi*g['lapse']*g['radial_length']*g['radius']*T['angular_pressure'])
    assert vertex['first']['shift'] == pytest.approx(-4*math.pi*g['radial_length']**2*g['radius']**2*T['outward_energy_flux'])
    # Independent polynomial differentiation and a four-parameter physical strain.
    qas = -.5j*(U.conj().T @ ALPHA @ UN-UN.conj().T @ ALPHA @ U)
    qis = -.5j*(U.conj().T @ UN-UN.conj().T @ U)
    ss = U.conj().T @ SIGMA1 @ U
    keys = ('lapse','radial_length','radius','shift')
    for i,x in enumerate(coordinates):
        for j,y in enumerate(coordinates):
            function = s.lambdify((N,a,r,beta,qa,qi,angular),s.diff(h,x,y),'numpy')
            expected = function(*g.values(),qas,qis,ss)*np.ones((2,2))
            assert vertex['second'][keys[i]][keys[j]] == pytest.approx(
                expected,abs=3e-15,rel=16*np.finfo(float).eps)
    direction = dict(lapse=.3,radial_length=-.2,radius=.4,shift=.15)
    step = 1e-6
    before = canonical_metric_vertices(U,UN,**{k:g[k]-step*direction[k] for k in g})
    after = canonical_metric_vertices(U,UN,**{k:g[k]+step*direction[k] for k in g})
    derivative = sum(direction[key]*vertex['first'][key] for key in keys)
    assert (after['hamiltonian_density']-before['hamiltonian_density'])/(2*step) == pytest.approx(derivative,abs=2e-9)


def test_transmission_orientation_and_physical_units():
    radius, a = .5, 1.3
    minus = stress_matrices(U,UN,radial_length=1.,radius=radius)
    plus = stress_matrices(U,UN+U,radial_length=1.,radius=radius)  # [xi_n]=k_wall*xi/2, k_wall=2
    for key in minus:
        assert minus[key] == pytest.approx(plus[key],abs=1e-15)
    reflected = stress_matrices(SIGMA1 @ U,-SIGMA1 @ UN,radial_length=1.,radius=radius)
    for key in minus:
        sign = -1 if key.startswith('outward_') else 1
        assert reflected[key] == pytest.approx(sign*minus[key],abs=1e-15)
    base = stress_matrices(U,UN,radial_length=a,radius=radius,angular_k=-2)
    L = 7.
    scaled = stress_matrices(U/math.sqrt(L),UN/(L*math.sqrt(L)),
                             radial_length=a,radius=radius*L,angular_k=-2)
    for key in base:
        power = 3 if 'charge' in key else 4
        assert scaled[key] == pytest.approx(base[key]/L**power,abs=2e-18)


def test_matrix_shapes_sample_axes_and_admitted_inputs():
    values = np.stack((U, U/2))
    rows = stress_matrices(values,np.stack((UN, UN/2)),radial_length=[1.,1.2],radius=[.5,.6])
    assert rows['energy_density'].shape == (2,2,2)
    for row in rows.values():
        assert row == pytest.approx(row.conj().swapaxes(-1,-2),abs=1e-15)
    for change in (dict(radius=0.),dict(radial_length=True),dict(angular_k=0)):
        with pytest.raises(ValueError):
            stress_matrices(U,UN,**(dict(radius=.5,radial_length=1.,angular_k=1)|change))
    with pytest.raises(ValueError): stress_matrices(U,UN[:,0],radius=.5,radial_length=1.)
