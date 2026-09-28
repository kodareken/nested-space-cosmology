#!/usr/bin/env python3
"""Exact Weyl-sequence control for the existing spatial throat Dirac operator."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from check_nsc_scale_closure import compare

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/nsc-3-radial-spectrum.json"


def calculate():
    rho, kappa, s = sp.symbols('rho kappa s', real=True)
    w = kappa / sp.sqrt(1+rho**2)
    primitive = kappa*sp.asinh(rho)
    plus = sp.simplify(w**2+sp.diff(w,rho))
    minus = sp.simplify(w**2-sp.diff(w,rho))
    bump = sp.sqrt(sp.Rational(15,16))*(1-s**2)
    norm = sp.integrate(bump**2,(s,-1,1))
    derivative_norm_squared = sp.integrate(sp.diff(bump,s)**2,(s,-1,1))
    assert norm == 1 and derivative_norm_squared == sp.Rational(5,2)
    assert sp.simplify(sp.diff(primitive,rho)-w) == 0
    n = sp.symbols('n',positive=True)
    bound = sp.sqrt(sp.Rational(5,2))/n + abs(kappa)/sp.sqrt(1+(n**2-n)**2)
    assert sp.limit(bound,n,sp.oo) == 0
    sigma2 = sp.Matrix([[0,-sp.I],[sp.I,0]])
    spinor = sp.Matrix([1,sp.I])/sp.sqrt(2)
    assert sigma2*spinor == spinor
    return {
        "schema": "nsc-radial-spectrum-v1", "artifact_id": "NSC-3-RADIAL-SPECTRUM",
        "classification": "gapless_essential_spectrum_of_the_frozen_spatial_throat_not_a_physical_particle_prediction",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "operator": {"D": "-i sigma_2 partial_rho + sigma_1 kappa/sqrt(1+rho²)",
                     "Hilbert_space": "L²(R,d_rho) tensor C²", "domain": "H¹(R) tensor C²",
                     "self_adjointness": "free self-adjoint Dirac plus bounded Hermitian multiplication",
                     "potential": str(w), "squared_partner_plus": str(plus), "squared_partner_minus": str(minus)},
        "Weyl_sequence": {
            "energy": "any real E", "spinor": "(1,i)/sqrt(2), sigma_2 eigenvalue +1",
            "bump": "f(s)=sqrt(15/16)(1-s²) for |s|<=1; zero otherwise; f in H¹",
            "bump_norm_squared": str(norm), "bump_derivative_norm_squared": str(derivative_norm_squared),
            "state": "psi_n(rho)=n^(-1/2)f((rho-n²)/n)exp(i E rho)v_plus",
            "support": "[n²-n,n²+n] escapes to +infinity; normalized states converge weakly to zero",
            "residual_norm_upper_bound": str(bound), "limit": "0",
            "bounds_at_kappa_1": [{"n": value,"upper_bound": float(bound.subs({n:value,kappa:1}))}
                                  for value in (4,8,16,32,64,128)],
            "conclusion": "every real E belongs to essential spectrum; spec_ess(D)=R",
        },
        "zero_energy": {"primitive_of_w": str(primitive),
                        "solutions": ["u=C exp(-kappa asinh(rho))", "v=C exp(+kappa asinh(rho))"],
                        "L2_global_zero_mode": False,
                        "reason": "for nonzero kappa each solution grows at one end; kappa=0 gives nonintegrable constants"},
        "consequence": "smooth radial throat and current-conserving gluing alone do not generate a nonzero spectral mass gap; finite box levels are not physical masses",
        "nonclaims": {"full_recursive_or_interacting_spectrum_gapless": False,
                      "compact_warped_or_Lorentzian_spectrum_proved_here": False,
                      "all_possible_embedded_eigenvalues_excluded_here": False,
                      "dynamical_Phi_or_stationary_scale_determined": False,
                      "nested_architecture_refuted": False,
                      "new_to_world_mathematical_priority_established": False},
        "terminal": True,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--output',type=Path)
    args=parser.parse_args(); record=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),record)
        print('Radial essential-spectrum proof reproduced; every exact and numerical field checked.')
    elif args.output:
        with args.output.open('x') as stream:
            json.dump(record,stream,sort_keys=True,indent=2,allow_nan=False); stream.write('\n')
        print(args.output)
    else: print(json.dumps(record,sort_keys=True,indent=2))


if __name__=='__main__': main()
