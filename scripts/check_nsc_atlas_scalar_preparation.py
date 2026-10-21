"""Coupled scalar-gravity preparation, physical DtN reduction and child state."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation

ROOT=Path(__file__).resolve().parents[1]


def identities():
    a,A,z,v,w,k,rv=sp.symbols("a A z v w k rv",nonzero=True,real=True)
    U=rv+w*w*z*z/2
    C=A*A-1-k*a*a*(v*v/2-U)/3
    derivative=lambda f: (sp.diff(f,a)*A-sp.diff(f,A)*k*a*(v*v+U)/3
                           +sp.diff(f,z)*v+sp.diff(f,v)*(w*w*z-3*A*v/a))
    si,so,d=sp.symbols("si so d",positive=True)
    zi,zo=-d*so/(si+so),d*si/(si+so)
    R,sigma,bi,bo=sp.symbols("R sigma bi bo",positive=True)
    jump=k*sigma*R/2
    deltaU=3*(bi*bi-bo*bo)/(k*R*R)
    force=sigma*3*(bi-bo)/(2*R)-deltaU
    results={
        "Euclidean_constraint_propagation":derivative(C),
        "bulk_to_boundary_scalar_action":derivative(a**3*z*v)-a**3*(v*v+w*w*z*z),
        "Schur_value_matching":zi-zo+d,
        "Schur_normal_matching":si*zi+so*zo,
        "Schur_quadratic":si*zi*zi+so*zo*zo-d*d*si*so/(si+so),
        "normal_force_from_constraints":force.subs(bo,jump-bi),
    }
    out={name:str(sp.factor(value)) for name,value in results.items()}
    if any(value!="0" for value in out.values()):
        raise RuntimeError(out)
    return out


def calculation():
    start=time.process_time()
    model=ScalarPreparation()
    coarse,mc=model.solve(tolerance=1e-7)
    fine,mf=model.solve(tolerance=1e-9)
    d=model.diagnostics(fine)
    tc=model.pole_transfer(coarse)
    tf=model.pole_transfer(fine)
    tp=model.pole_transfer(fine,side=1)
    reduced=model.scalar_reduction(fine)
    child=model.child_response(tf["pole_displacement"])
    radius_gap=float(fine.y[0,-1]-coarse.y[0,-1])
    pole_relative=abs(tf["pole_displacement"]/tc["pole_displacement"]-1)
    if (max(d["Hamilton_max"])>1e-9 or abs(d["force_residual"])>1e-9
            or abs(radius_gap)>1e-8 or pole_relative>1e-6
            or abs(reduced["bulk_boundary_action_gap"])>1e-8):
        raise RuntimeError("coupled preparation check failed")
    files=("src/recursive_horizons/nsc_atlas_scalar_preparation.py",
           "src/recursive_horizons/nsc_atlas_flux_conversion.py",
           "scripts/check_nsc_atlas_scalar_preparation.py",
           "tests/test_nsc_atlas_scalar_preparation.py",
           "scripts/plot_nsc_atlas_scalar_preparation.py")
    record=dict(
        schema="nsc-atlas-scalar-preparation-v1",
        action="Einstein + canonical scalar U=rho_b+m^2 phi^2/2+(q+mu phi)^2/2 + charged constant-tension membrane",
        inputs=dict(epsilon=model.scalar.epsilon,omega=model.scalar.omega,
                    parent_flux=model.parent_flux,child_flux=model.child_flux,
                    tension=model.tension,bare_density=model.scalar.bare_density,
                    normalization="Same normalized classical variables as flux conversion; no physical SI ruler selected."),
        preparation="Regular compact O(4) scalar/gravity saddle with two poles and full charged-wall conditions. This selects one candidate preparation, not a nucleation probability.",
        solved=d,scalar_reduction=reduced,child_pole=tf,parent_pole=tp,
        lorentzian_wall=model.lorentzian_wall(d["radius"],1.),
        child_response=child,
        comparison=dict(previous_instant_transfer_energy=model.scalar.prepare()["oscillation_energy"],
                        prepared_to_instant_transfer_energy=tf["scalar_energy_at_pole"]/model.scalar.prepare()["oscillation_energy"],
                        uniform_child_vacuum_density=model.child_vacuum,
                        maximum_prepared_excitation_to_vacuum=tf["scalar_energy_at_pole"]/model.child_vacuum,
                        interpretation="The selected regular branch puts scalar excitation mainly in its spatial layer, not in a large homogeneous child condensate. It does not realize the previous sudden-transfer preparation."),
        checks=dict(symbolic_zero_identities=identities(),coarse=mc,fine=mf,
                    radius_tolerance_difference=radius_gap,pole_log_transfer_relative_difference=pole_relative,
                    direct_collocation_pole_relative_difference=float(fine.y[2,0]/tf["pole_displacement"]-1),
                    sampled_min_positive_cap_scale_factor=float(min(np.min(fine.sol(np.linspace(.001,1,257))[0]),np.min(fine.sol(np.linspace(.001,1,257))[4])))),
        validity="Regular cap branch, O(3,1) wall continuation |eta|<=1, child proper-time response0<tau<=1. The analytic field/constraint/junction identities are exact; the nonlinear saddle is a converged numerical construction, not interval-certified existence or a uniqueness proof.",
        limits="No negative-mode count, determinant, tunnelling probability, preferred quantum state, black-hole collapse history, complete matter production, measured cosmological scale or infinite nesting theorem. The finite P002 dust parent is not imported into this different prepared geometry.",
        next="Determine a physically specified excited-parent or nonequilibrium preparation if a substantial child matter component is claimed; do not infer it from the previous instantaneous quench.",
        numerical_error="Residuals, two collocation tolerances, direct bulk/boundary action and independent Riccati transfer; not observational or model error bars.",
        process_cpu_seconds=time.process_time()-start,
        source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_scalar_preparation.py")
    if record["process_cpu_seconds"]>120:
        raise RuntimeError("finite preparation calculation exceeded its CPU budget")
    return record


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    encoded=json.dumps(calculation(),indent=2,allow_nan=False)+"\n"
    if len(encoded.encode())>18000:
        raise RuntimeError("compact evidence budget exceeded")
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8",newline="\n")
    else:
        print(encoded,end="")
