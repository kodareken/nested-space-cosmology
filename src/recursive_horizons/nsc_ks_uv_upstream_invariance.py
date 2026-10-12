"""Shared upstream scalar normalization drops out of the Sigma C4 difference.

This is a formal occupied-vacuum coefficient identity, not a bound on the
physical source, the actual C4 value, or the higher-order remainder C_M.
It uses the existing Dirac operator, continuity, fixed axial scale and matched
incoming radius. No physical source column or covariance is rescaled.
"""
from functools import lru_cache
import sympy as sp

from .nsc_ks_finite_history_uv_coefficients import (
    _apply_L0, _minor_A1, _sym_pauli, _sym_potential, _unit,
    vacuum_upstream_conditions, confirm_raw_ks_vertices,
)


def _zero(expr):
    real_derivatives = [d for d in expr.atoms(sp.Derivative)
                        if d.expr.is_real is True and all(v.is_real is True for v in d.variables)]
    expr = expr.xreplace({d: sp.Dummy(real=True) for d in real_derivatives})
    value = sp.simplify(sp.expand_complex(expr))
    if value != 0:
        raise ArithmeticError('upstream invariance identity failed: '+str(value))
    return '0'


@lru_cache(maxsize=1)
def upstream_normalization_identities():
    rho, z = sp.symbols('rho z', real=True)
    a = sp.Function('a', positive=True)(rho)
    r = sp.Function('r', positive=True)(rho, z)
    q = sp.Function('q', real=True)(rho, z)
    mass, ell, moment = sp.symbols('m ell k', real=True)
    real_a2, imag_a2 = sp.symbols('R h', real=True)
    _, _, s3 = _sym_pauli()
    potential = _sym_potential(mass, ell, r)
    Q = mass**2+ell**2/r**2
    identities = {}
    for s in (-1, 1):
        tag = f's{s:+d}'
        major, minor = (0, 1) if s == 1 else (1, 0)
        a0 = _unit(s)
        a1 = sp.zeros(2, 1)
        a1[major], a1[minor] = sp.I*q, _minor_A1(a, mass, ell, r, s)
        a2 = sp.zeros(2, 1)
        a2[major] = real_a2+sp.I*imag_a2
        l2 = sp.simplify(sp.expand_complex(sp.conjugate(a1[minor])*a1[minor]))
        identities[tag+'/minor_A1_norm'] = _zero(l2-a**2*Q/4)
        # Continuity at e^-2 gives T_s n2 = -2s/a^2 d_z |l|^2.
        # Since a_z=0, this is s*T_s q_z, with T_s q=-Q/2.
        identities[tag+'/transport_n2_minus_s_qz'] = _zero(
            -2*s/a**2*sp.diff(l2,z) - s*sp.diff(-Q/2,z))
        n2 = (a0.H*a2+a2.H*a0+a1.H*a1)[0]
        n2 = sp.simplify(sp.expand_complex(n2))
        identities[tag+'/density_n2'] = _zero(n2-(2*real_a2+q**2+l2))
        # k is transported common upstream data; it is NOT assigned zero.
        r_sub = (s*q.diff(z)+moment-q**2-l2)/2
        j1 = sp.im((a0.H*a1.diff(z))[0])-s*n2
        k1 = sp.im((a0.H*s3*a1.diff(z))[0])-(n2-2*l2)
        v1 = (a0.H*potential*a1+a1.H*potential*a0)[0]
        identities[tag+'/identity_current_e_minus1'] = _zero(
            j1.subs(real_a2,r_sub)+s*moment)
        identities[tag+'/potential_e_minus1'] = _zero(v1+a*Q)
        identities[tag+'/N_bracket_e_minus1'] = _zero(
            (v1+k1/a).subs(real_a2,r_sub)+a*Q/2+moment/a)
        identities[tag+'/identity_current_leading'] = _zero(-s*(a0.H*a0)[0]+s)
        identities[tag+'/N_bracket_leading'] = _zero(-s*(a0.H*s3*a0)[0]/a+1/a)
        identities[tag+'/order_zero'] = _zero(
            (a0.H*a1+a1.H*a0)[0])

    # A scalar formal series independent of rho,z commutes with the ACTUAL L0.
    c = sp.Symbol('c', complex=True)
    psi = sp.Matrix([sp.Function('u')(rho,z), sp.Function('v')(rho,z)])
    commutator = (_apply_L0(c*psi,a,potential,s3,rho,z)
                  -c*_apply_L0(psi,a,potential,s3,rho,z))
    for i in range(2):
        identities[f'constant_scalar_commutes_with_L0/{i}'] = _zero(commutator[i])
    variable_commutator = (_apply_L0(z*psi,a,potential,s3,rho,z)
                           -z*_apply_L0(psi,a,potential,s3,rho,z))
    for i in range(2):
        identities[f'spatial_scalar_commutator_control/{i}'] = _zero(
            variable_commutator[i]+(s3*psi)[i]/a**2)

    # Explicit normalization fixes formal A2 in the declared real-major gauge.
    eps = sp.Symbol('eps', real=True)
    x,y,u,v = sp.symbols('x y u v', real=True)
    ratio = eps*(x+sp.I*y)+eps**2*(u+sp.I*v)
    ratio_norm = sp.expand(ratio*sp.conjugate(ratio))
    major = sp.series((1+ratio_norm)**sp.Rational(-1,2),eps,0,4).removeO().expand()
    identities['normalized_real_major_A2'] = _zero(major.coeff(eps,2)+(x*x+y*y)/2)
    identities['normalized_imag_major_A2'] = _zero(sp.im(major.coeff(eps,2)))
    phase = sp.Symbol('phase', real=True)
    phased = sp.series(sp.exp(sp.I*phase*eps**2)*major,eps,0,3).removeO().expand()
    identities['normalization_alone_does_not_fix_phase'] = _zero(
        sp.im(phased.coeff(eps,2))-phase)

    # General Hermitian bilinear with one carrier derivative: lowest power eps^-1.
    c2r,c2i,c3r,c3i,c4r,c4i = sp.symbols('c2r c2i c3r c3i c4r c4i',real=True)
    f = 1+(c2r+sp.I*c2i)*eps**2+(c3r+sp.I*c3i)*eps**3+(c4r+sp.I*c4i)*eps**4
    jm,j0,j1,j2,j3 = sp.symbols('jm j0 j1 j2 j3',real=True)
    current = jm/eps+j0+j1*eps+j2*eps**2+j3*eps**3
    shifted = sp.expand(f*sp.conjugate(f)*current).coeff(eps,3)
    shift = 2*c2r*j1+2*c3r*j0+(c2r**2+c2i**2+2*c4r)*jm
    identities['cubic_current_under_common_scalar'] = _zero(shifted-j3-shift)
    identities['cubic_difference_when_lower_differences_vanish'] = _zero(
        shift.subs({jm:0,j0:0,j1:0}))

    # Direct massless check: an upstream real A2 shift also changes transported
    # h3 and n4. Holding those fixed would create a spurious uncancelled term.
    c, qz, rr2, rr2z, b, bz, B, Bz, h3z, n4 = sp.symbols(
        'c qz R2 R2z b bz B Bz h3z n4', real=True)
    for s in (-1, 1):
        current3 = h3z+rr2*qz-q*rr2z+s*(b*Bz-B*bz)-s*n4
        shifted3 = current3.xreplace({rr2:rr2+c, h3z:h3z+c*qz,
                                     n4:n4+2*c*s*qz})
        identities[f's{s:+d}/transported_h3_n4_cancel_Rup_in_J3'] = _zero(shifted3-current3)
        identities[f's{s:+d}/n4_shift_matches_continuity'] = _zero(
            -s/a**2*sp.diff(c*a**2*ell**2/r**2,z)
            -2*c*s*ell**2*r.diff(z)/r**3)

    rg,rr,aa,kk = sp.symbols('r_g r_ref a_common k_common',positive=True)
    delta_n1 = -aa*ell**2/2*(1/rg**2-1/rr**2)
    identities['Sigma_matching_radius_eliminates_delta_N1'] = _zero(delta_n1.subs(rg,rr))
    normal, chi, w, U = sp.symbols('normal chi w U', real=True)
    identities['declared_history_matches_radius_on_Sigma'] = _zero(
        (chi*(normal*w+normal**3*U/6)).subs(normal,0))
    if delta_n1.subs({aa:1,ell:1,rg:2,rr:1}) != sp.Rational(3,8):
        raise ArithmeticError('radius-mismatch negative control failed')
    return {
        'identities': identities,
        'conditional_normalized_gauge_major_A2': '-a_up^2*(m^2+ell^2/r_up^2)/8',
        'conditional_normalized_gauge_imag_major_A2': 0,
        'normalized_gauge_premise': 'entire major component real positive (1+|v|^2)^(-1/2), not only A0',
        'invariance_uses_n2_up_equals_zero': False,
        'invariance_uses_normalized_gauge_example': False,
        'declared_upstream_conditions': dict(vacuum_upstream_conditions()),
        'continuity_constant': 'k=n2-s*q_z; common transported data, not assumed zero',
        'identity_current_lower_coefficients': ['-s', '0', '-s*k'],
        'N_bracket_lower_coefficients': ['-1/a', '0', '-a*(m^2+ell^2/r^2)/2-k/a'],
        'raw_action_convention': 'N=-mu/(2*pi)*(V density+S3 current/a); beta=+mu/(2*pi)*I current',
        'raw_vertex_owner': dict(confirm_raw_ks_vertices()),
        'scalar_multiplier': '1+c2/e^2+c3/e^3+c4/e^4+O(e^-5), common and independent of rho,z',
        'cubic_coefficient_change': '2 Re(c2) J1 + 2 Re(c3) J0 + (|c2|^2+2 Re(c4)) Jminus1',
        'validity': 'formal occupied-vacuum history-minus-reference coefficient on Sigma, fixed a and matched r',
        'C4_difference_independent_of_common_scalar_upstream_constants': True,
        'all_upstream_state_dependence_eliminated': False,
        'whole_slab_invariance_claimed': False,
        'spatial_scalar_negative_control': '[L0,z]psi=-S3 psi/a^2, nonzero',
        'radius_mismatch_control_delta_N1': '3/8',
        'normalization_only_phase_control': 'Im(A2_up)=phase; unit norm alone does not fix it',
        'physical_source_matching_verified': False,
        'source_columns_changed': False,
        'numerical_C4': None, 'numerical_C_M': None,
        'physical_upstream_error': None, 'physical_local_gate': 'OPEN',
    }
