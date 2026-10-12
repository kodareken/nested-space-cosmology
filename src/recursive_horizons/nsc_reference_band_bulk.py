"""Exact bulk closure of the existing paired, unweighted reference band trace.

This evaluates only small exact Pauli coefficient identities and the order
bounds of the owned projector recursion. It neither samples large momenta
nor computes a physical source. Finite-k and physical endpoint terms are
not removed by the bulk result.
"""
from fractions import Fraction
from functools import lru_cache
from math import comb, isfinite

import sympy as sp


FIELDS = ('N', 'beta', 'a', 'r')
I = sp.eye(2)
S1 = sp.Matrix([[0, 1], [1, 0]])
S2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
S3 = sp.diag(1, -1)


def _residual(value):
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    reduced = [sp.simplify(v) for v in entries]
    if any(v != 0 for v in reduced):
        raise ArithmeticError('nonzero exact band coefficient: '+str(reduced))
    return [str(v) for v in reduced]


def _polar_coefficients(P, Pi, epsilon):
    """Differentiate the actual owned S/sqrt(tr(Pi P)) at inverse |k|=0."""
    S = Pi*P+(I-Pi)*(I-P)
    U = S/sp.sqrt(sp.trace(Pi*P))
    return (sp.simplify(U.subs(epsilon, 0)),
            sp.simplify(U.diff(epsilon).subs(epsilon, 0)))


@lru_cache(maxsize=1)
def exact_leading_coefficients():
    eps = sp.Symbol('epsilon', real=True)
    u, v, du, dv, uz, vz, dvz = sp.symbols('u v du dv u_z v_z delta_v_z', real=True)
    c, beta, dc, db = sp.symbols('c beta delta_c delta_beta', real=True)
    Pi = (I+S1)/2
    result = {}
    for sign in (-1, 1):
        P = (I+(u*eps*S1-v*eps*S2-sign*S3)/sp.sqrt(1+(u*u+v*v)*eps*eps))/2
        U0, U1 = _polar_coefficients(P, Pi, eps)
        expected0 = (I+sp.I*sign*S2)/sp.sqrt(2)
        expected1 = (u*(I-sp.I*sign*S2)/2-sp.I*v*S3)/sp.sqrt(2)
        dU1 = U1.diff(u)*du+U1.diff(v)*dv
        Uz1 = U1.diff(u)*uz+U1.diff(v)*vz
        D1 = sp.simplify(dU1*U0.H)
        Hlead = sign*(c*S3-beta*I)
        deltaHlead = sign*(dc*S3-db*I)
        hlead = sp.simplify(U0*Hlead*U0.H)
        connection = sp.simplify(sp.trace(Pi*D1))
        push_trace = sp.simplify(sp.trace(Pi*Uz1*deltaHlead*U0.H))
        band_trace = sp.simplify((-c-sign*beta)*sp.diff(connection, dv)*dvz)
        wanted_push = -sp.I*(sign*dc+db)*vz/2
        wanted_band = -sp.I*(sign*c+beta)*dvz/2
        # The first-Moyal k primitive is i times these traces, modulo the
        # separately retained spatial boundary. These do not vanish singly.
        push_primitive, band_primitive = sp.simplify(sp.I*push_trace), sp.simplify(sp.I*band_trace)
        odd = {v: -v, dv: -dv, vz: -vz, dvz: -dvz}
        # The extra half-commutator [delta H,P] in the SYMMETRIC insertion
        # has no constant k-boundary term. Its O(k) term is a total z
        # derivative under the compact test variation.
        qx, qy = sp.symbols('delta_Hx_z delta_Hy_z', real=True)
        Pinf = P.subs(eps, 0); Pone = P.diff(eps).subs(eps, 0)
        projector_constant = sp.trace(deltaHlead*Pone+(qx*S1+qy*S2)*Pinf)
        residuals = {
            'polar_frame_leading': _residual(U0-expected0),
            'polar_frame_first': _residual(U1-expected1),
            'unitarity_leading': _residual(U0*U0.H-I),
            'unitarity_first': _residual(U1*U0.H+U0*U1.H),
            'intertwining_leading': _residual(U0*Pinf*U0.H-Pi),
            'intertwining_first': _residual(U1*Pinf*U0.H+U0*Pone*U0.H+U0*Pinf*U1.H),
            'projected_connection': _residual(connection-sp.I*sign*dv/2),
            'push_primitive_trace': _residual(push_trace-wanted_push),
            'band_primitive_trace': _residual(band_trace-wanted_band),
            'paired_push_primitive': _residual(push_primitive+push_primitive.subs(odd, simultaneous=True)),
            'paired_band_primitive': _residual(band_primitive+band_primitive.subs(odd, simultaneous=True)),
            'symmetric_projector_constant': _residual(projector_constant),
        }
        result[str(sign)] = {
            'frame_leading': str(U0), 'frame_inverse_K_coefficient': str(U1),
            'connection_inverse_K_coefficient': str(connection),
            'push_k_primitive_constant': str(push_primitive),
            'band_k_primitive_constant': str(band_primitive),
            'unpaired_primitive_generically_nonzero': bool(push_primitive != 0 and band_primitive != 0),
            'residuals': residuals,
        }

    planar = {}
    w, dw, wz = sp.symbols('w delta_w w_z', real=True)
    for angular_sign in (-1, 1):
        target = (I-angular_sign*S2)/2
        for sign in (-1, 1):
            P = (I+(-angular_sign*w*eps*S2-sign*S3)/sp.sqrt(1+w*w*eps*eps))/2
            U0, U1 = _polar_coefficients(P, target, eps)
            dU1 = U1.diff(w)*dw
            connection = sp.trace(target*dU1*U0.H)
            push = sp.trace(target*U1.diff(w)*wz*sign*(dc*S3-db*I)*U0.H)
            planar[f'{angular_sign}/{sign}'] = {
                'projected_connection': _residual(connection),
                'push_primitive_constant': _residual(push),
                'band_primitive_constant': _residual((-c-sign*beta)*connection),
            }
    chiral = {}
    for sign in (-1, 1):
        target = (I-sign*S3)/2
        U0, U1 = _polar_coefficients(target, target, eps)
        chiral[str(sign)] = {'frame': _residual(U0-I), 'frame_variation': _residual(U1),
                            'normal_band_commutator': _residual(target*sign*(c*S3-beta*I)-sign*(c*S3-beta*I)*target)}
    return {
        'definitions': {'epsilon': '1/abs(k)', 's': 'sign(k)', 'u': 'm*a',
                        'v': 'lambda*a/r', 'c': 'N/a', 'D': 'delta U star U^dagger'},
        'massive_chart': result, 'massless_nonzero_angular_chart': planar,
        'LLL': {'U': 'I in each fixed chiral k-sign chart', 'delta_U': '0',
                'band_exchange_all_orders': '0', 'gap_zero_floor_introduced': False,
                'band_exchange_scope': 'each fixed nonzero-k chiral chart of the owned canonical reference',
                'distributional_gap_zero_bridge_proved': False,
                'residuals': chiral},
        'chart_overlap': {
            'm_positive': '(1+m*a/sqrt(k^2+(m*a)^2+(lambda*a/r)^2))/2 >= 1/2',
            'massless_nonzero_angular': '(1+abs(lambda)*a/r/sqrt(k^2+(lambda*a/r)^2))/2 >= 1/2',
        },
    }


def _p_bound(order, spatial=0, momentum=0):
    if order == 0:
        return -1-momentum if spatial or momentum else 0
    return -order-1-momentum


def _h_bound(momentum):
    # H and delta H are exactly at most linear in canonical k, including beta*k.
    return 1-momentum if momentum <= 1 else None


def _star_power(left, right, order):
    powers = []
    for i in range(order+1):
        for j in range(order-i+1):
            ell = order-i-j
            for q in range(ell+1):
                a, b = left(i, ell-q, q), right(j, q, ell-q)
                if a is not None and b is not None:
                    powers.append(a+b)
    return max(powers)


def _series_bounds(powers, *, constant_leading=True):
    def bound(order, coordinate, momentum):
        if order >= len(powers):
            return None
        if order == 0 and constant_leading and (coordinate or momentum):
            return -1-momentum
        return powers[order]-momentum
    return bound


def recursion_power_certificate():
    """Integer valuations of every existing recursion term through order four."""
    rows = []
    for n in range(1, 5):
        F = [_p_bound(n-1, spatial=1)]
        for j in range(n):
            ell = n-j
            for h_momentum in range(ell+1):
                h = _h_bound(h_momentum)
                if h is not None:
                    F.append(h+_p_bound(j, spatial=h_momentum, momentum=ell-h_momentum))
        G = []
        for i in range(n):
            for j in range(n):
                ell = n-i-j
                if ell >= 0:
                    for momentum_left in range(ell+1):
                        G.append(_p_bound(i, spatial=ell-momentum_left, momentum=momentum_left)
                                 +_p_bound(j, spatial=momentum_left, momentum=ell-momentum_left))
        f, g = max(F), max(G)
        offdiagonal = 1+f-2  # spin commutator divided by the quadratic normal gap
        p = max(g, offdiagonal)
        if (f, g, p) != (-n, -n-2, -n-1):
            raise ArithmeticError('projector-recursion order bound failed')
        primitive = []
        for ell in range(1, n+1, 2):
            for i in range(n-ell+1):
                j = n-ell-i
                a = -i-1  # a z derivative removes U0's constant limit; D0 is O(1/K)
                b = 1 if j == 0 else -j  # deltaH star U^dagger or effective h
                primitive.append(a+b-(ell-1))
        bound = max(primitive)
        rows.append({'formal_order': n, 'F_power': f, 'G_power': g,
                     'P_power': p, 'U_and_deltaU_power': -n-1,
                     'D_power': -n-1, 'effective_h_power': -n,
                     'k_boundary_power_before_pairing': bound,
                     'higher_boundary_decays': n > 1 and bound < 0,
                     'recursion_residuals': [f+n, g+n+2, p+n+1]})
    # Propagate the existing polar star-normalization, rather than assigning
    # its higher-order rates. Dividing by its gapped leading root is O(1).
    S = [0]+[row['P_power'] for row in rows]
    A = [_star_power(_series_bounds(S), _series_bounds(S), n) for n in range(5)]
    Z = [0]
    for n in range(1, 5):
        ZA = [_star_power(_series_bounds(Z), _series_bounds(A), j) for j in range(n+1)]
        Z.append(_star_power(_series_bounds(ZA), _series_bounds(Z), n))
    U = [_star_power(_series_bounds(S), _series_bounds(Z), n) for n in range(5)]
    deltaU = [-1]+U[1:]
    D = [_star_power(_series_bounds(deltaU, constant_leading=False), _series_bounds(U), n) for n in range(5)]
    H = lambda n, z, k: _h_bound(k) if n == 0 else None
    UH = [_star_power(_series_bounds(U), H, n) for n in range(5)]
    # UH's leading coefficient is O(k), not a metric-independent constant.
    uh = lambda n, z, k: (1-k if n == 0 else UH[n]-k) if n < len(UH) else None
    for row in rows:
        n = row['formal_order']
        h = max(_star_power(uh, _series_bounds(U), n), _series_bounds(U)(n-1, 1, 0))
        residuals = [A[n]+n+1, Z[n]+n+1, U[n]+n+1, D[n]+n+1, h+n]
        if any(residuals):
            raise ArithmeticError('polar-frame/connection order propagation failed')
        row.update(polar_norm_power=A[n], inverse_star_root_power=Z[n],
                   U_and_deltaU_power=U[n], D_power=D[n], effective_h_power=h,
                   frame_recursion_residuals=residuals)
    return {
        'power_convention': 'O(abs(k)^power), at either fixed k sign',
        'leading_P_U': 'metric-independent constant plus O(1/abs(k))',
        'momentum_derivative': 'improves a decaying power by one; H has no k derivatives >=2',
        'coordinate_and_metric_variation': 'preserve decay powers for smooth finite jets',
        'scalar_shift': '-beta*k is retained; it cancels only in the ordinary spin commutator',
        'rows': rows,
    }


def exact_star_trace_coefficients():
    """Leibniz/IBP coefficients of the existing spatial Weyl trace identity.

    integral_z tr[A,B]_star = 2 sum_odd ell (i/2)^ell/ell!
        partial_k^ell integral_z tr[(partial_z^ell A) B].
    Spatial endpoint terms are kept outside this compact-test bulk identity.
    """
    rows = []
    for ell in range(1, 5):
        prefactor = (sp.I/2)**ell/sp.factorial(ell)
        residuals = []
        for j in range(ell+1):
            # Moving B's j spatial derivatives onto A cancels (-1)^j in
            # star_order. Reversing the commutator leaves 1-(-1)^ell.
            from_star = prefactor*comb(ell, j)*(1-(-1)**ell)
            from_total_derivative = (2 if ell % 2 else 0)*prefactor*comb(ell, j)
            residuals += _residual(from_star-from_total_derivative)
        rows.append({'Moyal_order': ell, 'k_primitive_coefficient': str((2 if ell % 2 else 0)*prefactor),
                     'primitive_k_derivative_order': ell-1, 'coefficient_residuals': residuals})
    return rows


def angular_pairing_certificate(channels):
    """Use the retained degeneracy contract, including both actual angular signs.

    nsc_incoming_source_assembly selects (-1,+1) when lambda!=0.
    nsc_incoming_state_moments.incoming_group_factor assigns deg/n_signs,
    with copy count once. This exact certificate concerns that same trace.
    """
    rows = []
    seen = set()
    for channel in channels:
        index = channel['index']
        if index in seen:
            raise ValueError('distinct retained group indices required')
        seen.add(index)
        copies, degeneracy = channel['copy_count'], channel['degeneracy']
        if any(isinstance(x, bool) or not isinstance(x, int) or x <= 0 for x in (copies, degeneracy)):
            raise ValueError('positive integer retained copy and angular multiplicities required')
        angular = channel['angular_eigenvalue']; mass = channel['compact_mass']
        if (not isinstance(angular, (float, int)) or not isinstance(mass, (float, int))
                or not isfinite(angular) or not isfinite(mass) or angular < 0 or mass < 0):
            raise ValueError('retained nonnegative angular magnitudes and compact masses required')
        signs = [-1, 1] if angular else [1]
        weight = Fraction(copies*degeneracy, len(signs))
        rows.append({'group': index, 'angular_signs': signs,
                     'weight_per_sign': {'numerator': weight.numerator, 'denominator': weight.denominator},
                     'chart': ('positive_S1' if mass else 'opposite_angular_S2' if angular else 'exact_chiral_LLL'),
                     'odd_angular_weight_sum': str(sum(weight*s for s in signs)) if angular else 'not required: lambda=0'})
    if not rows:
        raise ValueError('the retained nonempty channel inventory is required')
    return {'groups': rows, 'group_count': len(rows),
            'signed_family_count': sum(len(r['angular_signs']) for r in rows),
            'paired_group_count': sum(len(r['angular_signs']) == 2 for r in rows),
            'sign_owner': 'nsc_incoming_source_assembly.assemble',
            'weight_owner': 'nsc_incoming_state_moments.incoming_group_factor'}


def paired_unweighted_bulk_closure(channels, *, full_momentum_line, unweighted,
                                   compact_test_variations, smooth_positive_metric_jets):
    """Return the proved BULK Euler exchange, subject to its explicit domain.

    The result is not a finite-cutoff cancellation, an EndpointBranchJets
    completion, a source-convergence proof, or a statement about Gamma_rest.
    """
    if any(flag is not True for flag in (full_momentum_line, unweighted,
                                         compact_test_variations, smooth_positive_metric_jets)):
        raise ValueError('full unweighted k trace and compact tests on smooth positive N,a,r jets required')
    leading = exact_leading_coefficients()
    powers = recursion_power_certificate()
    star = exact_star_trace_coefficients()
    inventory = angular_pairing_certificate(channels)
    return {
        'bulk_action_gradient': {field: 0 for field in FIELDS},
        'formal_orders': [0, 1, 2, 3, 4],
        'leading_exact_certificate': leading, 'recursion_power_certificate': powers,
        'star_trace_certificate': star, 'angular_pairing_certificate': inventory,
        'remainder_identity': 'i*epsilon*dT tr(Pi D) + tr[Pi U, deltaH star Udag]_star + tr[Pi D,h]_star + tr[deltaH,P]_star/2',
        'boundary_accounting': {'finite_k_endpoints': 'required separately',
                                'physical_time_endpoints': 'required separately',
                                'physical_spatial_endpoints': 'required separately'},
        'scope': {'paired_unweighted_full_k_bulk_exchange_closed': True,
                  'finite_momentum_window_exchange_zero': False,
                  'extra_momentum_weight_allowed': False,
                  'EndpointBranchJets_completed': False,
                  'physical_source_spectral_convergence_proved': False,
                  'stationarity_or_constraint_solution_claimed': False,
                  'Gamma_rest_assigned': False,
                  'raw_heat_equivalence_required': False,
                  'light_geometric_allocation_changed': False,
                  'LLL_trace_convention': 'sum of the two owned fixed nonzero-k chiral charts; retain the separate c=4 geometric allocation',
                  'LLL_distributional_gap_zero_bridge_proved': False,
                  'normal_band': 'lower normal energy; never sign of shifted coordinate Hamiltonian',
                  'metric_domain': 'smooth finite coordinate and variation jets; N,a,r positive and nondegenerate',
                  'infinite_inventory_limit': 'closure is per paired group; spectral summation convergence is separate'},
    }
