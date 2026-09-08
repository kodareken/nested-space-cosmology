#!/usr/bin/env python3
"""Check the spin-frame/domain map for the imported free KK reduction.

This is a convention check, not a new spectrum scan or a selected NSC state.
See docs/nsc-compact-mass-map.md for the metric and domain assumptions.
"""
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_spinor_bridge import weyl_matrices


def check_map():
    frame = weyl_matrices()
    beta, gamma5 = frame["gamma"][0], frame["gamma5"]
    left, right = frame["left"], frame["right"]
    alpha_y = sp.I * beta * gamma5
    y = sp.Symbol("Y", real=True)
    length = sp.Symbol("L_star", positive=True)
    n = sp.Symbol("n", integer=True, positive=True)
    mass = n * sp.pi / (2 * length)
    phase = mass * (y + length)
    f_left = sp.cos(phase) / sp.sqrt(length)
    f_right = sp.sin(phase) / sp.sqrt(length)
    embedding = f_left * left + f_right * right
    sigma = -sp.Rational(18, 1015) * (y / length) ** 2
    momenta = sp.symbols("p_x p_y p_z", real=True)
    shift = sp.Symbol("shift", real=True)
    kinetic = sum((a * p for a, p in zip(frame["alpha"], momenta)), sp.zeros(4))
    kinetic -= shift * momenta[0] * sp.eye(4)

    identities = {
        "compact_H_intertwines_exact_mass": (
            -sp.I * alpha_y * embedding.diff(y) - embedding * beta * mass
        ),
        "4D_massless_kinetic_preserves_embedding": kinetic * embedding - embedding * kinetic,
        "chiral_endpoint_current_form_vanishes": left * alpha_y * left,
        "warp_spin_connection_cancels": (
            (sp.exp(-2 * sigma) * embedding).diff(y)
            + 2 * sigma.diff(y) * sp.exp(-2 * sigma) * embedding
            - sp.exp(-2 * sigma) * embedding.diff(y)
        ),
        "odd_profile_at_both_ends": sp.Matrix([
            f_right.subs(y, -length), f_right.subs(y, length)
        ]),
        "even_profile_derivative_at_both_ends": sp.Matrix([
            f_left.diff(y).subs(y, -length), f_left.diff(y).subs(y, length)
        ]),
    }
    for name, value in identities.items():
        residual = sp.simplify(value)
        if residual != sp.zeros(*value.shape):
            raise AssertionError((name, residual))
    norms = [sp.simplify(sp.integrate(f ** 2, (y, -length, length)))
             for f in (f_left, f_right)]
    if norms != [1, 1]:
        raise AssertionError(("massive_mode_normalization", norms))
    zero_gram = sp.integrate(1 / (2 * length), (y, -length, length)) * left
    if zero_gram != left:
        raise AssertionError(("zero_mode_normalization", zero_gram))
    return {
        "scope": "published free KK reduction applied to the declared length-2 interval",
        "exact_residuals": {name: "0" for name in identities},
        "massive_mode_chiral_norms": [str(v) for v in norms],
        "zero_mode_gram": "P_L",
        "mass_in_inverse_length": "n*pi/(2*L_star)",
        "mass_over_cutoff": "n*pi/(2*sqrt(zeta))",
        "sigma_prime_at_endpoints": [str(sp.simplify(sigma.diff(y).subs(y, end)))
                                     for end in (-length, length)],
        "physical_domain_or_occupations_selected": False,
        "quantum_determinant_or_backreaction_computed": False,
    }


if __name__ == "__main__":
    print(json.dumps(check_map(), indent=2))
