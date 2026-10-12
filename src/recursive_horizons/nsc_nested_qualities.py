"""Exact finite-window controls of the nested-qualities proposition."""
from __future__ import annotations

import sympy as sp


def finite_window(local, link, omega, first: int, depth: int):
    h, b, scale = sp.Matrix(local), sp.Matrix(link), sp.sympify(omega)
    if h.rows == 0 or h.rows != h.cols or h != h.H or b.shape != h.shape:
        raise ValueError('use a Hermitian local matrix and a matching link')
    if scale.is_positive is not True or scale.is_finite is not True or isinstance(first, bool) or not isinstance(first, int):
        raise ValueError('use a positive scale and integer starting label')
    if isinstance(depth, bool) or not isinstance(depth, int) or depth < 1:
        raise ValueError('depth must be a positive integer')
    d = h.rows
    result = sp.zeros(d * depth)
    for j in range(depth):
        a = slice(j*d, (j+1)*d)
        result[a, a] = scale**(first+j)*h
        if j+1 < depth:
            nxt = slice((j+1)*d, (j+2)*d)
            result[a, nxt] = scale**(first+j)*b
            result[nxt, a] = scale**(first+j)*b.H
    return result


def exact_controls():
    h = sp.Matrix([[1, sp.I/5], [-sp.I/5, 2]])
    b = sp.Matrix([[sp.Rational(1,4), sp.I/7], [sp.Rational(1,9), sp.Rational(1,6)]])
    omega, z, depth = sp.Rational(3,2), 1+2*sp.I, 3
    j = finite_window(h,b,omega,0,depth)
    shifted = finite_window(h,b,omega,1,depth)
    simplify = lambda m: m.applyfunc(sp.simplify)
    k = z*sp.eye(j.rows)-j
    f = k[4:6,4:6]
    middle = k[2:4,2:4]-k[2:4,4:6]*f.inv()*k[4:6,2:4]
    schur = k[:2,:2]-k[:2,2:4]*middle.inv()*k[2:4,:2]
    gamma = z/omega**2*sp.eye(2)-h
    wrong = gamma
    for n in (1,0):
        gamma = z/omega**n*sp.eye(2)-h-b*gamma.inv()*b.H/omega
        wrong = z/omega**n*sp.eye(2)-h-b*wrong.inv()*b.H
    c = sp.diag(*[sp.Rational(n,5) for n in range(6)])
    cstat = sp.eye(6)/2
    noise = b*c[2:4,2:4]*b.H
    checks = {
        'hermitian_coupled_window': j == j.H and b != sp.zeros(2),
        'normalized_shift_exact': simplify(shifted/omega-j) == sp.zeros(6),
        'unscaled_shift_rejected': simplify(shifted-j) != sp.zeros(6),
        'joined_equals_nested_resolvent': simplify(k.inv()[:2,:2]-schur.inv()) == sp.zeros(2),
        'normalized_recurrence_exact': simplify(gamma-schur) == sp.zeros(2),
        'omitted_omega_factor_rejected': simplify(wrong-schur) != sp.zeros(2),
        'reversed_coupling_rejected': simplify(b*middle.inv()*b.H-b.H*middle.inv()*b) != sp.zeros(2),
        'distinct_admissible_occupations': all(0<=x<=1 for x in c.diagonal()) and c[:2,:2] != c[2:4,2:4],
        'same_law_stationary_state': j*cstat-cstat*j == sp.zeros(6),
        'same_law_nonstationary_state': j*c-c*j != sp.zeros(6),
        'exterior_state_changes_noise': noise != sp.zeros(2),
        'total_energy_derivative_zero': sp.simplify(sp.trace((-sp.I*(j*c-c*j))*j)) == 0,
    }
    if not all(checks.values()):
        raise AssertionError({k:v for k,v in checks.items() if not v})
    return {'dimension_per_region':2,'depth':depth,'omega':str(omega),
            'spectral_parameter':str(z),'checks':{k:bool(v) for k,v in checks.items()}}
