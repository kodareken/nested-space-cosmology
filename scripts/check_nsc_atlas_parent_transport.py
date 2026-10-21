"""Bounded causal and energy transfer calculation on the prepared geometry."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import sympy as sp
from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
from recursive_horizons.nsc_atlas_parent_transport import (
    RadialCapTransport, independent_child_arrival, de_sitter_embedding_control)

ROOT = Path(__file__).resolve().parents[1]


def identities():
    z, c = sp.symbols("z c", positive=True)
    phiL, phiR, qp, eL, eR, mu, w = sp.symbols("phiL phiR qp eL eR mu w", real=True, nonzero=True)
    minimum = lambda q: -mu*q/w**2
    phi_overlap = minimum(qp-eL)+minimum(qp-eR)-minimum(qp)
    K, b, a = sp.symbols("K b a", positive=True)
    raw = {
        "regularized_timelike_primitive": sp.diff(-sp.log(z+sp.sqrt(z*z+c)),z)+1/sp.sqrt(z*z+c),
        "wall_chart_mass_shell": (K*K/a**2-1)-K*K/a**2+1,
        "child_chart_mass_shell": -(1+K*K/b**2)+K*K/b**2+1,
        "additive_flux_passage_lands_on_minimum": phi_overlap-minimum(qp-eL-eR),
    }
    out = {name:str(sp.simplify(value)) for name,value in raw.items()}
    if any(v != "0" for v in out.values()):
        raise RuntimeError(out)
    return out


def calculation():
    start = time.process_time()
    model = ScalarPreparation()
    coarse, _ = model.solve(tolerance=1e-7)
    fine, meta = model.solve(tolerance=1e-9)
    tc = RadialCapTransport.from_preparation(model, coarse)
    t = RadialCapTransport.from_preparation(model, fine)
    rows = [t.massive(.3,g) for g in (1.1,1.5,2.,5.)]
    primary = rows[2]
    source_angle = math.asin(t.child_h*primary["source"]["a_source"])
    control = de_sitter_embedding_control(t.child_h,source_angle,2.)
    control["energy_ratio"] = control["gamma_received"]/2
    old = tc.massive(.3,2.)
    threshold, old_threshold = t.central_threshold(.3), tc.central_threshold(.3)
    ode_gaps = [independent_child_arrival(t.child_h,r["radial_charge"],r["log_label"])-r["arrival_time"]
                for r in rows if r["central_arrival"]]
    spreads = {key:abs(primary[key]-old[key])
               for key in ("arrival_time","energy_ratio","gamma_received")}
    spreads["threshold_gamma"] = abs(threshold-old_threshold)
    if max(spreads.values()) > 2e-8 or max(map(abs,ode_gaps)) > 2e-8:
        raise RuntimeError("transport checks exceed declared numerical comparison tolerance")
    files = ("src/recursive_horizons/nsc_atlas_parent_transport.py",
             "src/recursive_horizons/nsc_atlas_scalar_preparation.py",
             "src/recursive_horizons/nsc_atlas_flux_conversion.py",
             "scripts/check_nsc_atlas_parent_transport.py",
             "tests/test_nsc_atlas_parent_transport.py",
             "scripts/plot_nsc_atlas_parent_transport.py")
    record = dict(
        schema="nsc-atlas-parent-transport-v1",
        question="Which pre-existing parent sources can reach the child center, when, and with what locally measured energy?",
        background="The existing scalar/four-form/gravity saddle and constant-tension membrane; unchanged couplings and solved spatial layer.",
        action_parameters=dict(epsilon=model.scalar.epsilon,omega=model.scalar.omega,
            kappa=model.scalar.kappa,parent_flux=model.parent_flux,charge=model.charge,
            tension=model.tension,bare_density=model.scalar.bare_density),
        preparation="Specified inward radial test bodies at parent proper depth .3 L_ref on eta=0; gamma is source input. No abundance or added material species is selected.",
        units="c=1, lengths and proper times in L_ref; physical time is tau_hat L_ref/c. The SI ruler remains unselected.",
        geometry=dict(wall_radius=float(fine.y[0,-1]),child_length=t.child_length,
                      parent_length=t.parent_length,wall_x=t.wall_x,child_h=t.child_h,
                      vacuum_Lambda_Lref_squared=model.scalar.kappa*model.child_vacuum),
        central_threshold_gamma=threshold,reference_trajectories=rows,
        null_reference=t.null(.3),parent_null_collar=t.null_collar(),
        matched_unwalled_control=dict(
            geometry="Single open de Sitter space with the same child h, no membrane or parent vacuum change.",
            matching="Same initial areal radius a_source, inward gamma=2 and time-symmetric initial slice. Source is on the hemisphere nearest the receiving pole.",
            source_angle=source_angle,**control,
            nested_minus_control_arrival=primary["arrival_time"]-control["arrival_time"],
            nested_minus_control_energy_ratio=primary["energy_ratio"]-control["energy_ratio"],
            relative_received_energy_change=primary["energy_ratio"]/control["energy_ratio"]-1,
            limit="Distinguishes this supplied source/geometry pair, not all alternatives or arbitrary unknown preparations."),
        checks=dict(symbolic_zero_identities=identities(),fine_background=meta,
                    two_background_tolerance_spreads=spreads,independent_child_ODE_time_gaps=ode_gaps,
                    conservative_numerical_comparison_tolerance=2e-8),
        child_approximation=dict(initial_scalar_energy=model.pole_transfer(fine)["scalar_energy_at_pole"],
            bracket="b_hminus(tau)<=b_actual(tau)<=b_hplus(tau); hplus^2-hminus^2=kappa E0/3. At fixed transmitted K and label, tau_hminus<=tau_actual<=tau_hplus when both arrive. The zero-source de Sitter comparison is used numerically.",
            meaning="Analytic background comparison, not a bound on finite-density source backreaction or scalar-wave scattering."),
        collision_control=dict(parent_flux=1.5,equal_charge=.5,overlap_flux=.5,
            overlap_vacuum=.03125,homogeneous_excitation_at_minima=0.,
            domain="Linear scalar equation on one prescribed metric with additive swept flux and continuing prescribed membranes. This is not a solved self-gravitating collision or an annihilation rule."),
        validity="Radial freely falling test bodies with no nongravitational membrane interaction. Primary reception tau<6; analytic central-access thresholds include the asymptotic vacuum comparison. O(3,1) background continuation is reused, without a new evolution grid.",
        limits="No population normalization, coupled scalar-wave transmission, finite-density backreaction, isotropization, cold-matter abundance, clustering spectrum, new vacuum selection, or measured cosmological correlation. Entering the child does not imply reaching its center.",
        next="Use a constraint-compatible parent scalar packet in the accessible collar and derive its full scalar/metric/membrane perturbation matching before computing deposited stress and clustering.",
        process_cpu_seconds=time.process_time()-start,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_parent_transport.py")
    if record["process_cpu_seconds"] > 120:
        raise RuntimeError("120 CPU-second calculation budget exceeded")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    encoded = json.dumps(calculation(),indent=2,allow_nan=False)+"\n"
    if len(encoded.encode()) > 18000:
        raise RuntimeError("compact record budget exceeded")
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding="utf-8",newline="\n")
    else:
        print(encoded,end="")
