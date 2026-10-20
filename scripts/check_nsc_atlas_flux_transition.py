"""A charged parent-child source, two shared walls and a cosmological join."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import mpmath as mp
import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_flux_transition import reference_geometry, evolve_reference

ROOT = Path(__file__).resolve().parents[1]


def symbolic_checks():
    r, m, dl, li, k = sp.symbols("r m dl li k", positive=True)
    dm = sp.symbols("dm", real=True)
    bm = dm/(k*r*r)+(dl/(6*k)+k/2)*r
    bp = bm-k*r
    fm = 1-2*m/r-li*r*r/3
    fp = 1-2*(m+dm)/r-(li+dl)*r*r/3
    V = fm-bm*bm
    a = -sp.diff(V, r)/2
    km = (a+sp.diff(fm, r)/2)/bm+2*bm/r
    kp = (a+sp.diff(fp, r)/2)/bp+2*bp/r
    n, f, e, kap, bare = sp.symbols("n f e kap bare", real=True)
    L = lambda j: bare+kap*(f-j*e)**2/2
    checks = {
        "unsquared_constraint": bm-bp-k*r,
        "same_velocity_on_both_sides": bm*bm-fm-bp*bp+fp,
        "full_temporal_Israel": kp-km+3*k,
        "normal_force": k*(km+kp)-dl,
        "tidal_trace_from_metric": -sp.diff(fm,r,2)/2-sp.diff(fm,r)/r-li,
        "inheritance_second_difference": L(n+2)-2*L(n+1)+L(n)-kap*e*e,
    }
    result = {name: str(sp.factor(expr)) for name, expr in checks.items()}
    if any(value != "0" for value in result.values()):
        raise RuntimeError(result)
    return result


def independent_quadratures(end=.25):
    """40-digit separated first integral, not the evolution RHS.

    R=R0+sign*u^2 removes the turning-point square root exactly. The
    polynomial divided difference avoids subtraction of nearly equal V.
    """
    with mp.workdps(40):
        k = mp.mpf(1)/10
        lambdas = [mp.mpf(9)/8, mp.mpf(1)/2, mp.mpf(1)/8]
        masses = [mp.mpf(0)]*3
        for j, r in [(2, mp.mpf(2)/5), (1, mp.mpf(4)/5)]:
            beta = mp.sqrt(1-2*masses[j]/r-lambdas[j]*r*r/3)
            masses[j-1] = (masses[j]+(lambdas[j]-lambdas[j-1])*r**3/6
                           +k*r*r*beta-k*k*r**3/2)
        answers = []
        for j, r0, sign, seed in ((1, mp.mpf(4)/5, 1, .16),
                                   (2, mp.mpf(2)/5, -1, .29)):
            A = (masses[j-1]-masses[j])/k
            B = (lambdas[j-1]-lambdas[j])/(6*k)+k/2
            D, C = masses[j]+A*B, lambdas[j]/3+B*B
            def integrands(u):
                r = r0+sign*u*u
                power = lambda n: sum(r**(n-1-i)*r0**i for i in range(n))
                Q = power(4)-2*D*power(3)-C*power(6)
                Z = -sign*Q
                dtau = 2*r*r/mp.sqrt(Z)
                v2 = u*u*Z/r**4
                fmid = 1-2*masses[1]/r-lambdas[1]*r*r/3
                other = 0 if j == 1 else 2
                fother = 1-2*masses[other]/r-lambdas[other]*r*r/3
                return dtau*mp.sqrt(v2+fmid)/fmid, dtau, dtau*mp.sqrt(v2+fother)/fother
            u = mp.findroot(lambda z: mp.quad(lambda x: integrands(x)[0], [0,z])-end,
                            (mp.mpf(seed)*.99, mp.mpf(seed)*1.01))
            answers.append(dict(radius=float(r0+sign*u*u),
                                proper_time=float(mp.quad(lambda x: integrands(x)[1], [0,u])),
                                other_static_time=float(mp.quad(lambda x: integrands(x)[2], [0,u]))))
        return dict(masses=[float(x) for x in masses], walls=answers, decimal_digits=40)


def calculation():
    started = time.process_time()
    law, outer, inner, parent = reference_geometry()
    exact = symbolic_checks()
    a, meta = evolve_reference()
    b, meta_fine = evolve_reference(rtol=2e-13, atol=2e-15)
    nodes = np.linspace(0,.25,65)
    difference = float(np.max(np.abs(a.sol(nodes)-b.sol(nodes))))
    y = b.y[:,-1]
    independent = independent_quadratures()
    qi = independent["walls"]
    gaps = [y[0]-qi[0]["radius"],y[2]-qi[0]["proper_time"],y[3]-qi[0]["other_static_time"],
            y[4]-qi[1]["radius"],y[6]-qi[1]["proper_time"],y[7]-qi[1]["other_static_time"]]
    residuals = {key: 0. for key in ("angular_residual", "temporal_residual", "force_residual", "energy_residual")}
    min_f, clearance = 1., 1.
    for column in b.sol(nodes).T:
        for wall, offset in ((outer,0),(inner,4)):
            g = wall.geometry(column[offset], column[offset+1])
            for key in residuals:
                residuals[key] = max(residuals[key],abs(g[key]))
            min_f = min(min_f,g["f_inside"],g["f_outside"])
        clearance = min(clearance,column[0]-column[4])
    state = parent.state(y[8])
    clock_gap = parent.static_time(y[8])-y[3]
    if max(abs(x) for x in gaps) > 1e-10 or max(residuals.values()) > 1e-9 or difference > 1e-9:
        raise RuntimeError("independent evolution / constraint tolerance failed")
    if abs(clock_gap) > 1e-10 or min_f <= 0 or clearance <= 0:
        raise RuntimeError("common cosmological chart check failed")
    files = ("src/recursive_horizons/nsc_atlas_flux_transition.py",
             "scripts/check_nsc_atlas_flux_transition.py", "tests/test_nsc_atlas_flux_transition.py",
             "scripts/plot_nsc_atlas_flux_transition.py")
    record = dict(
        schema="nsc-atlas-flux-transition-v1", candidate="Einstein/canonical four-form/charged pure-tension membranes/dust",
        borrowed="Brown-Teitelboim flux transition; Bousso-Polchinski action/boundary ensemble; Hirano 1804.09985 Eqs2.1-2.4 and2.28-2.31; Kottler and Einstein-Straus matching",
        connection="One action determines regional vacuum jumps, both wall forces, a shared intermediate mass, cosmological matching and a mass-independent tidal trace.",
        inputs=dict(kappa=law.kappa, bare_Lambda_Lref_squared=0.,
                    sqrt_kappa_parent_flux_Lref=1.5, sqrt_kappa_charge_Lref=.5,
                    k_Lref=law.k, initial_radii_Lref=[.8,.4], initial_proper_speeds=[0.,0.],
                    central_geometric_mass=0., initial_cosmological_join_Lref=1.2,
                    physical_scale="Lref is an independent reference unit; no SI value, charge or parent flux is inferred from observations."),
        inherited=[law.inherit(1.5,j) for j in range(3)],
        vacuum_changes_Lref_squared=[-.625,-.375], curvature_second_difference_Lref_squared=.25,
        regions_outside_to_inside=[dict(mass=reg.mass,Lambda=reg.Lambda) for reg in (outer.outside,outer.inside,inner.inside)],
        common_clock=dict(name="middle Kottler T/Lref", interval=[0.,.25],
                          outer=outer.geometry(y[0],y[1]), inner=inner.geometry(y[4],y[5]),
                          outer_proper_time=float(y[2]),inner_proper_time=float(y[6]),
                          parent_static_time=float(y[3]),core_static_time=float(y[7])),
        cosmological_parent=state,
        observer=dict(kind="Geodesic-deviation trace in static orthonormal laboratories, excluding walls and dust",
                      radii_Lref=[1.1,.6,.1],
                      local=[reg.tidal(r) for reg,r in zip((outer.outside,outer.inside,inner.inside),(1.1,.6,.1))],
                      relation="A_radial+2 A_transverse=Lambda; successive trace differences have second difference kappa*e^2",
                      SI="Relative acceleration per separation is c^2*A. Divide its trace by c^2 before comparing with Lambda.",
                      uncharged_control="With the same parent flux, gravity, tension and radii but e=0, all regional tidal traces equal the parent Lambda despite altered masses."),
        checks=dict(symbolic_zero_identities=exact, max_sampled_junction_residuals=residuals,
                    sampled_min_wall_lapse=min_f, sampled_min_wall_separation=clearance,
                    tolerance_comparison_max_state_difference=difference,
                    independent_quadratures=independent, max_quadrature_endpoint_gap=max(abs(x) for x in gaps),
                    common_parent_clock_gap=clock_gap, ode_evaluations=[meta["evaluations"],meta_fine["evaluations"]]),
        validity="Ordinary spherical timelike thin walls on 0<=T_middle/Lref<=.25, regular center, nonintersecting shells and timelike cosmological join. Existence is supported analytically by speed bounds and positive static collars; sampled residuals alone are not a proof.",
        uncertainty="Numerical differences check this benchmark, not model uncertainty, observational error bars or interval-certified endpoints.",
        remaining="Parent flux, charge, bare vacuum, tension, initial radii and physical ruler remain inputs. No nucleation probability, selection of observed Lambda, reheating, cosmological matter fractions or nonspherical stability is derived. This action does not identify the earlier inherited Schur operator.",
        field_grids=0, process_cpu_seconds=time.process_time()-started,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_flux_transition.py",
    )
    if record["process_cpu_seconds"] > 120:
        raise RuntimeError("finite calculation budget exceeded")
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(calculation(), indent=2, allow_nan=False)+"\n"
    if len(encoded.encode()) > 18000:
        raise RuntimeError("compact result budget exceeded")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8", newline="\n")
    else:
        print(encoded, end="")
