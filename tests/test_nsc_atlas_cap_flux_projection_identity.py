"""Exact mass adjoint, reaction cancellation and flux projection controls.

No saved fields, numerical radial solves, trajectories or grid scans are used.
These check the reference mathematics, not a repaired native cap evolution.
"""
import pytest
import sympy as s

from control_budget import ControlBudget

CONTROLS = ControlBudget(maximum_evaluations=1, cpu_seconds=20.)


@pytest.fixture(autouse=True)
def exact_budget():
    with CONTROLS.active():
        yield


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("cap flux projection algebra:", CONTROLS.budget.evaluations,
          "field samples;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def test_linearized_radial_constraint_and_integrating_factor():
    r, alpha, n, f, p, q = s.symbols("r alpha n f p q", nonzero=True)
    dn, df, dp, dq, U = s.symbols("dn df dp dq U")
    P = s.Function("P")(f)
    A = 1-r*r-2*alpha*n/r
    a = alpha*(p*p+(q-f)**2)/(2*r)
    rho = r*r*P/2+A*(p*p+(q-f)**2)/4
    varied = sum(s.diff(rho, v)*dv for v, dv in
                 [(n, dn), (f, df), (p, dp), (q, dq)])
    varied = varied.subs(s.diff(P, f), U)
    source = r*r*U*df/2+A*(p*dp+(q-f)*(dq-df))/2
    assert s.expand(varied+a*dn-source) == 0
    W, wr, dnr = s.symbols("W W_r dn_r")
    product = wr*dn+W*dnr
    assert s.expand(product.subs({wr: a*W, dnr: source-a*dn})-W*source) == 0
    assert s.expand(varied-source) != 0  # backreaction cannot be dropped


def test_nonlocal_primitive_adjoint_against_exact_polynomial_integrals():
    r, t = s.symbols("r t", positive=True)
    R = s.Rational(3, 5)
    f, p, q = 1+3*r*r, r+2*r**3, 1+9*r*r
    A, W, U = 1-r*r, (2+r)/3, 2*f+4*f**3
    dp, dq = r-r**3, 1+2*r*r
    df = s.integrate(dq, (r, 0, r))/r
    hp = W*A*p/2
    hq = W*A*(q-f)/2
    hf = W*(r*r*U-A*(q-f))/2
    tail = s.integrate((hf/r).subs(r, t), (t, r, R))
    direct = s.integrate(s.expand(hp*dp+hq*dq+hf*df), (r, 0, R))
    adjoint = s.integrate(s.expand(hp*dp+(hq+tail)*dq), (r, 0, R))
    assert direct == adjoint
    assert s.integrate(s.expand(tail*dq), (r, 0, R)) != 0


def test_moment_gradient_includes_cross_dependence_of_both_reconstructions():
    r, ip, iq = s.symbols("r ip iq", real=True)
    R = s.Rational(1, 2)
    # Both integral moments are exact, but each shape also depends on the
    # other moment, as can happen through characteristic reconstruction.
    p = 2*ip*r/R**2+iq*r*(r*r-R*R/2)
    q = iq/R+ip*(r*r-R*R/3)
    f = s.integrate(q, (r, 0, r))/r
    assert s.integrate(p, (r, 0, R)) == ip
    assert s.integrate(q, (r, 0, R)) == iq
    A = 1-r*r
    rho = r*r*f*f/2+A*(p*p+(q-f)**2)/4
    mass = s.integrate(s.expand(rho), (r, 0, R))
    gp, gq = A*p/2, A*q/2-(1-R*R)*f.subs(r, R)/2
    for variable in [ip, iq]:
        pair = gp*s.diff(p, variable)+gq*s.diff(q, variable)
        weight = s.integrate(s.expand(pair), (r, 0, R))
        assert s.expand(s.diff(mass, variable)-weight) == 0
    omitted_cross = s.integrate(s.expand(gq*s.diff(q, ip)), (r, 0, R))
    assert omitted_cross != 0


def test_reaction_tail_follows_from_mass_lapse_and_primitive_constraints():
    r, alpha = s.symbols("r alpha", positive=True)
    n, f, p, q, W = (s.Function(k)(r) for k in ["n", "f", "p", "q", "W"])
    P, U = s.symbols("P U")
    A = 1-r*r-2*alpha*n/r
    a = alpha*(p*p+(q-f)**2)/(2*r)
    rho = r*r*P/2+A*(p*p+(q-f)**2)/4
    hf = W*(r*r*U-A*(q-f))/2
    # Differentiate gq-W*A*q/2, using gq=hq+integral_r^R hf/s.
    difference_rate = s.diff(-W*A*f/2, r)-hf/r
    rules = {s.diff(W, r): a*W, s.diff(n, r): rho,
             s.diff(f, r): (q-f)/r}
    reaction_without_lapse = r*((U-2*f)+alpha*(2*n/r**3-P)*f)
    assert s.expand(difference_rate.subs(rules)+W*reaction_without_lapse/2) == 0
    # W'=a W and delta'=a imply W/exp(delta) is constant in r.
    delta = s.Function("delta")(r)
    assert s.simplify(s.diff(W*s.exp(-delta), r).subs(
        {s.diff(W, r): a*W, s.diff(delta, r): a})) == 0


def test_reaction_tail_cancels_in_the_full_continuum_mass_flux():
    r = s.symbols("r")
    c, C = s.symbols("c C")
    u, v, H, S = (s.Function(k)(r) for k in ["u", "v", "H", "S"])
    pairing = c*u*(s.diff(v, r)-S)/2+(c*(v+H)/2-C)*s.diff(u, r)
    primitive = c*u*(v+H)/2-C*u
    assert s.expand((pairing-s.diff(primitive, r)).subs(s.diff(H, r), -S)) == 0
    missing_reaction = pairing+c*u*S/2
    assert s.expand((missing_reaction-s.diff(primitive, r)).subs(
        s.diff(H, r), -S)) == c*u*S/2
    A, lapse, p, q, f = s.symbols("A lapse p q f", positive=True)
    wall = primitive.subs({u: A*lapse*p, v: A*lapse*q, H: 0,
                          c: 1/lapse, C: A*f/2})
    assert s.expand(wall-A*A*lapse*p*(q-f)/2) == 0


def test_mass_gradient_generates_both_original_field_equations():
    r = s.symbols("r")
    lapse_wall, C = s.symbols("lapse_wall C", nonzero=True)
    u, v, H, S = (s.Function(k)(r) for k in ["u", "v", "H", "S"])
    gp, gq = u/(2*lapse_wall), (v+H)/(2*lapse_wall)-C
    p_rate = (2*lapse_wall*s.diff(gq, r)).subs(s.diff(H, r), -S)
    q_rate = 2*lapse_wall*s.diff(gp, r)
    assert s.expand(p_rate-s.diff(v, r)+S) == 0
    assert s.expand(q_rate-s.diff(u, r)) == 0
    assert s.expand(gp*p_rate+gq*q_rate-s.diff(2*lapse_wall*gp*gq, r).subs(
        s.diff(H, r), -S)) == 0


def test_fixed_mesh_gradient_flux_telescopes_without_an_added_mass_projection():
    lapse = s.symbols("lapse")
    wp = s.symbols("wp0:4")
    wq = s.symbols("wq0:4")
    fp0, fq0, fpn, fqn = s.symbols("fp0 fq0 fpn fqn")
    fp = [fp0]+[lapse*(wq[j-1]+wq[j]) for j in range(1, 4)]+[fpn]
    fq = [fq0]+[lapse*(wp[j-1]+wp[j]) for j in range(1, 4)]+[fqn]
    direct = sum(wp[j]*(fp[j+1]-fp[j])+wq[j]*(fq[j+1]-fq[j]) for j in range(4))
    endpoints = (wp[-1]*fpn+wq[-1]*fqn-2*lapse*wp[-1]*wq[-1]
                 -wp[0]*fp0-wq[0]*fq0+2*lapse*wp[0]*wq[0])
    assert s.expand(direct-endpoints) == 0
    # Arbitrary physical boundary traces do not automatically give the
    # desired wall flux. This identity exposes, rather than erases, that gap.
    assert s.expand(endpoints) != 0


def test_interface_conjugate_jump_and_exact_upwind_projection_split():
    A, W, f, tail = s.symbols("A W f tail")
    pl, pr, ql, qr = s.symbols("pL pR qL qR")
    wp, wq, sp, sm = s.symbols("delta_wp delta_wq abs_sp abs_sm")
    jp, jm = pr+qr-pl-ql, pr-qr-pl+ql
    gpl, gpr = W*A*pl/2, W*A*pr/2
    gql, gqr = W*A*(ql-f)/2+tail, W*A*(qr-f)/2+tail
    for sign, jump in [(1, jp), (-1, jm)]:
        assert s.expand(gpr-gpl+sign*(gqr-gql)-W*A*jump/2) == 0
    up = (sp*jp+sm*jm)/4
    uq = (sp*jp-sm*jm)/4
    direct = wp*up+wq*uq
    eta_plus, eta_minus = wp+wq+W*A*jp/2, wp-wq+W*A*jm/2
    trace_dissipation = -W*A*(sp*jp**2+sm*jm**2)/8
    projection = (sp*eta_plus*jp+sm*eta_minus*jm)/4
    assert s.expand(direct-trace_dissipation-projection) == 0


def test_positive_upwind_work_for_an_actual_geometric_mass_functional():
    r = s.symbols("r")
    faces = [s.Rational(1, 4), s.Rational(1, 2), s.Rational(3, 4)]
    p = [s.S.One, s.Rational(9, 8)]
    # An interior two-cell, piecewise-constant representation, q=f=0.
    # This is not a claim about the native limited-linear reconstruction.
    weights = [s.integrate((1-r*r)*value/2, (r, lo, hi))/(hi-lo)
               for value, lo, hi in zip(p, faces[:-1], faces[1:])]
    A = 1-faces[1]**2
    jump = p[1]-p[0]
    work = (weights[0]-weights[1])*A*jump/2
    trace = -A*A*jump**2/4  # both characteristic speeds have magnitude A
    assert work == s.Rational(67, 16384) > 0
    assert trace == -s.Rational(9, 4096) < 0
    assert work-trace == s.Rational(103, 16384)


def test_zero_backreaction_energy_is_mass_plus_an_explicit_wall_term():
    r = s.symbols("r", positive=True)
    f = s.Function("f")(r)
    p = s.symbols("p")
    q = f+r*s.diff(f, r)
    A = 1-r*r
    rho = r*r*f*f/2+A*(p*p+(q-f)**2)/4
    energy = A*(p*p+q*q)/4
    boundary = A*r*f*f/4
    assert s.expand(rho-energy+s.diff(boundary, r)) == 0
    assert s.expand(rho-energy) != 0


def test_geometric_coordinate_preserves_moments_and_transforms_the_wave_system():
    z, T = s.symbols("z T", real=True)
    r, A = s.tanh(z), 1-s.tanh(z)**2
    u, v = s.Function("u")(z, T), s.Function("v")(z, T)
    p, q = u/A, v/A
    assert s.simplify(s.diff(r, z)-A) == 0
    assert s.simplify(p*s.diff(r, z)-u) == 0
    assert s.simplify(q*s.diff(r, z)-v) == 0
    assert s.simplify(A*(p*p+q*q)*s.diff(r, z)/4-(u*u+v*v)/4) == 0
    # d_r=(1/A)d_z; A is T independent in the alpha=0 limit.
    assert s.simplify(A*(s.diff(p, T)-s.diff(A*q, z)/A)
                      -(s.diff(u, T)-s.diff(v, z))) == 0
    assert s.simplify(A*(s.diff(q, T)-s.diff(A*p, z)/A)
                      -(s.diff(v, T)-s.diff(u, z))) == 0


def test_moving_wall_energy_and_scalar_boundary_term_recover_mass_flux():
    R, chi, u, v, rx, tx = s.symbols("R chi u v R_x T_x", nonzero=True)
    A, f = 1-R*R, chi/R
    K = A*chi*chi/(4*R)
    # chi_T=u, chi_z=v, z_x=R_x/A; retain both wall motions.
    k_rate = s.diff(K, R)*rx+s.diff(K, chi)*(tx*u+rx*v/A)
    energy_rate = tx*u*v/2+rx*(u*u+v*v)/(4*A)
    rho = R*R*f*f/2+A*((u/A)**2+(v/A-f)**2)/4
    flux = u*(v-A*f)/2
    assert s.cancel(energy_rate-k_rate-rx*rho-tx*flux) == 0
    assert s.cancel(energy_rate-rx*rho-tx*flux) != 0
    assert s.cancel(energy_rate-k_rate-tx*flux) != 0


def test_moving_geometric_cells_have_a_signed_energy_flux_without_projection():
    ul, ur, vl, vr, tx, speed, sp, sm = s.symbols("uL uR vL vR Tx speed abs_sp abs_sm")
    eleft, eright = (ul*ul+vl*vl)/4, (ur*ur+vr*vr)/4
    fu = tx*(vl+vr)/2+speed*(ul+ur)/2
    fv = tx*(ul+ur)/2+speed*(vl+vr)/2
    # For E_i=(Iu_i^2+Iv_i^2)/(4*h_i), geometry contributes
    # -h_i' e_i; at this interior face the pair is speed*[e].
    central = (ul-ur)*fu/2+(vl-vr)*fv/2+speed*(eright-eleft)
    assert s.expand(central+tx*(ur*vr-ul*vl)/2) == 0
    jp, jm = ur+vr-ul-vl, ur-vr-ul+vl
    up, uq = (sp*jp+sm*jm)/4, (sp*jp-sm*jm)/4
    dissipation = (ul-ur)*up/2+(vl-vr)*uq/2
    assert s.expand(dissipation+(sp*jp**2+sm*jm**2)/8) == 0
    # Omitting the moving-cell volume derivative fails even centrally.
    assert s.expand(central-speed*(eright-eleft)+tx*(ur*vr-ul*vl)/2) != 0


def test_incoming_wall_and_reflecting_centre_have_explicit_energy_penalties():
    u, v, incoming, tx, speed = s.symbols("u v incoming Tx speed")
    uh, vh = (incoming+u-v)/2, (incoming-u+v)/2
    fu, fv = tx*vh+speed*uh, tx*uh+speed*vh
    boundary = u*fu/2+v*fv/2-tx*u*v/2-speed*(u*u+v*v)/4
    physical = tx*uh*vh/2+speed*(uh*uh+vh*vh)/4
    penalty = (tx+speed)*(u+v-incoming)**2/8
    assert s.expand(boundary-physical+penalty) == 0
    # At the left endpoint, preserve the outgoing u+v and impose u_hat=0.
    centre_fu, centre_fv = tx*(u+v), s.S.Zero
    centre = -u*centre_fu/2-v*centre_fv/2+tx*u*v/2
    assert s.expand(centre+tx*u*u/2) == 0


def test_full_moving_geometric_cell_mass_identity_including_wall_primitive():
    R, chi, tx, incoming = s.symbols("R chi Tx incoming", nonzero=True)
    u, v = s.symbols("u0:3"), s.symbols("v0:3")
    speeds = [s.S.Zero, *s.symbols("speed1:4")]
    # Every face is timelike: Tx>0 and |speed_j|<Tx. This permits
    # |s_plus|=Tx+speed and |s_minus|=Tx-speed in the upwind flux.
    fu, fv = [tx*(u[0]+v[0])], [s.S.Zero]
    losses = tx*u[0]**2/2
    for j in [1, 2]:
        jp = u[j]+v[j]-u[j-1]-v[j-1]
        jm = u[j]-v[j]-u[j-1]+v[j-1]
        sp, sm = tx+speeds[j], tx-speeds[j]
        fu.append(tx*(v[j-1]+v[j])/2+speeds[j]*(u[j-1]+u[j])/2+(sp*jp+sm*jm)/4)
        fv.append(tx*(u[j-1]+u[j])/2+speeds[j]*(v[j-1]+v[j])/2+(sp*jp-sm*jm)/4)
        losses += (sp*jp*jp+sm*jm*jm)/8
    uh, vh = (incoming+u[-1]-v[-1])/2, (incoming-u[-1]+v[-1])/2
    fu.append(tx*vh+speeds[-1]*uh)
    fv.append(tx*uh+speeds[-1]*vh)
    losses += (tx+speeds[-1])*(u[-1]+v[-1]-incoming)**2/8
    energy_rate = sum(u[j]*(fu[j+1]-fu[j])/2+v[j]*(fv[j+1]-fv[j])/2
                      -(speeds[j+1]-speeds[j])*(u[j]**2+v[j]**2)/4 for j in range(3))
    A, f = 1-R*R, chi/R
    rx = A*speeds[-1]
    K = A*chi*chi/(4*R)
    # chi_W=sum I_v; its rate is the SAME conservative boundary flux.
    k_rate = s.diff(K, R)*rx+s.diff(K, chi)*(fv[-1]-fv[0])
    rho = R*R*f*f/2+A*((uh/A)**2+(vh/A-f)**2)/4
    flux = uh*(vh-A*f)/2
    defect = energy_rate-k_rate-rx*rho-tx*flux
    assert s.cancel(s.expand(defect+losses)) == 0
    assert s.expand(losses) != 0  # finite-cell mass flux is dissipative
