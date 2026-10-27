#!/usr/bin/env python3
"""Bounded finite nonlinear transport calculation and source record."""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import time

from scipy.integrate import quad

from recursive_horizons.nsc_atlas_vortex_transfer import (
    angular_advection_average, bump, direct_integrals, dissipation, energy,
    exact_theta_rate, mean_theta_rate, moments, radial_fluxes,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/"results/summaries/nsc-atlas-vortex-transfer-v1.json"
MODES = (4, 8, 16)
ORDERS = (64, 128, 256)
SOURCE_PATHS = ["src/recursive_horizons/nsc_atlas_vortex_transfer.py",
                "scripts/check_nsc_atlas_vortex_transfer.py",
                "tests/test_nsc_atlas_vortex_transfer.py"]


def task(job):
    mode, order = job
    cpu = time.process_time()
    direct = direct_integrals(mode, order)
    return dict(mode=mode, order=order, **direct,
                separated_energy=energy(mode, order), separated_loss=dissipation(mode, order=order),
                cpu_seconds=time.process_time()-cpu)


def main():
    started = time.perf_counter()
    jobs = [(m, q) for m in MODES for q in ORDERS]
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(task, jobs))
    total_points = sum(r["quadrature_points"] for r in rows)+sum(2*q for q in ORDERS)
    cpu = sum(r["cpu_seconds"] for r in rows)
    if cpu > 60 or total_points > 500_000:
        raise RuntimeError("finite transport budget exhausted before evidence completion")
    coefficients = [dict(order=q, **moments(q)) for q in ORDERS]
    I_adaptive, I_adaptive_error = quad(lambda z: float(bump(z)[0])**2, -1., 1.,
                                       epsabs=1e-13, epsrel=1e-13)
    result = []
    for mode in MODES:
        levels = [r for r in rows if r["mode"] == mode]
        c, m, f = levels
        indicators = {}
        for name in ("energy", "curl_loss"):
            first, second = abs(c[name]-m[name]), abs(m[name]-f[name])
            if not second < first:
                raise RuntimeError("quadrature did not contract; diagnose before another calculation")
            cross = abs(f[name]-f["separated_energy" if name == "energy" else "separated_loss"])
            if name == "curl_loss":
                cross = max(cross, abs(f["curl_loss"]-f["gradient_loss"]))
            indicators[name] = 2*max(second, cross, 1e-13)
        rates = []
        for r in (.75, 1.25):
            actual = float(mean_theta_rate(r, mode))
            direct = angular_advection_average(r, 0., mode)
            exact = float(exact_theta_rate(r))
            if max(abs(actual-exact), abs(direct-exact)) > 5e-13:
                raise RuntimeError("independent nonlinear angular force check failed")
            rates.append(dict(radius=r, rate=actual, direct_advection_rate=direct,
                              exact_rate=exact, comparison_tolerance=5e-13))
        result.append(dict(mode=mode, levels=levels, indicators=indicators,
                           initial_rates=rates, energy_over_initial_loss=f["energy"]/f["curl_loss"]))

    flux = {k: float(v) for k, v in radial_fluxes(1., order=256).items()}
    source = dict(
        author="OpenAI", title="Finite time blowup for Navier–Stokes", release_date="2026-09-08",
        announcement="https://openai.com/index/navier-stokes-solution/",
        paper="https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf",
        public_pdf_pages=166, formalization_repository="https://github.com/openai/NavierStokesAndEuler",
        formalization_commit="f9e8bc5b38b6e212696e8a30e3e91517af887bbd",
        supplied_text_bytes=523447,
        supplied_text_sha256="46c37997b150309499f6c4383ad0df2d7b0ac91fdd9ba9b52ceb4a3296af36fd",
        reading_boundary="Theorem1.1, Sections2.1-2.2 and selected3.2-3.3; independent Grok selected7.5/AppendixC review; not full proof audit",
        theorem_scope="3D incompressible NS, nu>0, zero initial velocity, smooth compact forcing, bounded energy and finite-time velocity blowup",
        lean_rebuilt=False, entire_theorem_independently_verified=False,
    )
    L, U, rho, nu = 1., .1, 1000., 1e-4
    record = dict(
        schema="NSC-ATLAS-VORTEX-TRANSFER-v1", date="2026-10-27", source=source,
        scope="Exact initial rotational-mean acceleration of locally classical unforced Navier-Stokes flow from specified nonzero compact data",
        observational_validation=False, blowup_demonstrated=False, sustained_cascade_demonstrated=False,
        preparation=dict(radial_support=[.5, 1.5], height_support=[-1., 1.],
                         radial_bump="exp(1-1/(1-4(r-1)^2)) inside support", phase="m theta-2m log r",
                         potential="(0,-a(r)chi(z)/(2m),a(r)chi(z)/m)*cos(phase)",
                         viscosity=.001, modes=MODES, quadrature_orders=ORDERS,
                         initial_rotational_mean=0, full_initial_velocity="nonzero prepared curl field"),
        exact_rates=dict(inner="-1024 exp(-2/3)/81", outer="1024 exp(-2/3)/225"),
        coefficients=coefficients, adaptive_Ichi=I_adaptive, adaptive_Ichi_error=I_adaptive_error,
        mode_results=result, flux_at_unit_radius=flux,
        source_core_scalings=dict(velocity_exponent="-1/2-h", volume_exponent="3/2-h",
                                 energy_exponent="1/2-3h", core_loss_exponent="-1/2-3h",
                                 critical_L3_exponent="-4h/3", h_domain="0<h<1/100",
                                 note="power counting checked; not the global blowup construction"),
        physical_example=dict(length_m=L, speed_m_s=U, density_kg_m3=rho, viscosity_m2_s=nu,
                              reference_time_s=L/U, reference_energy_J=rho*U*U*L**3,
                              inner_mean_acceleration_m_s2=float(exact_theta_rate(.75))*U*U/L,
                              outer_mean_acceleration_m_s2=float(exact_theta_rate(1.25))*U*U/L,
                              angular_momentum_flux_Nm=flux["angular"]*rho*U*U*L**3,
                              axial_momentum_flux_N=flux["axial"]*rho*U*U*L**2,
                              note="supplied finite fluid benchmark, not an astronomical identification"),
        compute=dict(points=total_points, point_limit=500_000, solver_cpu_seconds=cpu,
                     cpu_limit_seconds=60, elapsed_seconds=time.perf_counter()-started, workers=4,
                     note="CPU excludes process imports/startup; tensor-product quadrature, no NS time grid"),
        source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(record, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(dict(output=str(OUTPUT.relative_to(ROOT)), points=total_points,
                          cpu_seconds=cpu, rates=record["physical_example"], coefficients=coefficients[-1])))


if __name__ == "__main__":
    main()
