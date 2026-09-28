"""Independent countercontrols for a surface-scoped UV cancellation."""
import numpy as np
import sympy as sp

from recursive_horizons.nsc_ks_uv_upstream_invariance import upstream_normalization_identities
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter


def test_actual_operator_and_continuity_identities():
    result = upstream_normalization_identities()
    assert set(result['identities'].values()) == {'0'}
    for s in (-1, 1):
        for key in ('transport_n2_minus_s_qz', 'identity_current_e_minus1',
                    'N_bracket_e_minus1', 'transported_h3_n4_cancel_Rup_in_J3'):
            assert result['identities'][f's{s:+d}/'+key] == '0'
    assert result['physical_source_matching_verified'] is False
    assert result['physical_upstream_error'] is None
    assert result['numerical_C4'] is None and result['numerical_C_M'] is None


def test_pointwise_partial_R_derivative_does_not_represent_boundary_change():
    # Holding h3,n4 fixed leaves +c*qz. Their actual linked shifts remove it.
    c,qz,R,Rz,q,h3z,n4 = sp.symbols('c qz R Rz q h3z n4', real=True)
    for s in (-1,1):
        J = h3z+R*qz-q*Rz-s*n4
        assert sp.expand(J.subs(R,R+c)-J) == c*qz
        transformed = J.xreplace({R:R+c,h3z:h3z+c*qz,n4:n4+2*s*c*qz})
        assert sp.expand(transformed-J) == 0


def test_scalar_series_countercontrols_require_matched_lower_coefficients():
    eps = sp.Symbol('eps', real=True)
    f = 1+eps**2
    # A nonzero first-order history difference changes the cubic coefficient.
    delta = sp.Rational(3,8)*eps
    assert sp.expand((f*f-1)*delta).coeff(eps,3) == sp.Rational(3,4)
    # When the first difference starts at order three, that coefficient stays.
    delta = sp.Rational(7,11)*eps**3+sp.Rational(5,13)*eps**4
    assert sp.expand(f*f*delta).coeff(eps,3) == sp.Rational(7,11)
    # A scalar change at order e^-1 is excluded by the fixed major-A1 datum.
    assert sp.expand((1+eps)**2*eps**2).coeff(eps,3) == 2


def test_real_major_gauge_is_stronger_than_unit_norm():
    eps, phase = sp.symbols('eps phase', real=True)
    major = (1+eps**2)**sp.Rational(-1,2)
    phased = sp.series(sp.exp(sp.I*phase*eps**2)*major,eps,0,3).removeO()
    assert sp.im(sp.expand(phased).coeff(eps,2)) == phase
    assert upstream_normalization_identities()['normalization_only_phase_control'].startswith('Im(A2_up)=phase')


def test_action_signs_match_real_source_column_owner():
    for s in (-1,1):
        F = np.zeros((1,2,1),complex)
        F[0,0 if s==1 else 1,0] = 1
        Fz = -2j*s*F
        empty = np.zeros((0,*F.shape),complex)
        result = source_column_matter(F,Fz,np.eye(1),empty,empty,
            mass=0,angular=0,axial_scale=1,radius=1,multiplicity=1)
        # Columns are already weighted; no second 1/(2*pi) is inserted.
        assert np.array_equal(result['action_gradient'][0], [2.,-2.*s])
