# Dirac vacuum layer at the neck, and the relation between a layer's pressure/energy ratio and its energy growth

7 Oct 2026. Script `casimir_layer_check.py`, data `casimir_layer_check.json`, figure `casimir_layer_check.png`. Continues `torsion-layer-note.md`. Needs `israel_matching.py`, `junction_neck.py` and the restored `child_*` files.
Setup as before: 10 solar-mass parent, N_need child (radiation), thread path from `israel-matching-note.md`. G = c = 1 in formulas.

## 0. Correction to the torsion-layer note (own error)

The lifetime ratios (`3e33` at `R = 1 m` up to `3.7e40`) were computed on the thread path only. That path is the slowest (the `T = const` geodesic), so it has the longest shell proper time and the smallest energy. They are not a property of torsion matter.
Along any neck, `|E| tau ~ R^2`: `|E| tau = 0.667 R^2` on the thread path (`R/r_g` from `1e-12` to `0.3`) and `8.0 R^2` on the run in section 4. For a layer of fixed surface energy, `tau ~ 1/(6 pi |sigma|)` independent of `R`. A torsion layer at `2 n_b` that is one Planck length thick has `|sigma| = 30/l_P` and `tau ~ 1e-46 s`, which is 0.02 of the torsion dynamical time. So lifetime does not exclude torsion matter on a faster path. The pressure-sign result of the torsion-layer note (its section 1) stands and is the real obstruction (section 2 below).

## 1. Dirac vacuum layer

Massless Dirac field between plates, bag boundary conditions, zeta-regularised half-integer spectrum, four fermionic oscillators: `E/A = -(7 pi^2/2880) hbar c/d^3` (coefficient `0.023989`; the scalar Dirichlet check gives `-pi^2/1440`, matched to 16 digits). The energy density is negative and the tangential pressure is positive.
Boost invariance along the plates gives `p_parallel = -rho`, so a planar vacuum layer has surface pressure `P = -sigma`: ratio `w = -1`, for any field.

## 2. Exact relation along an outgoing neck

From the Codazzi identity (verified to 25 digits in `israel-matching-note.md`), with `q = d ln|E|/d ln R`, `E = 4 pi R^2 sigma`:
`sigma (q + 2w) Rdot/R = -T_un`, where `T_un = (eps + p) a chi' gamma^2` is the child fluid's energy flux through the neck, `<= 0` when `a chi' <= 0` (outgoing neck, child fluid at rest). With `sigma < 0` this gives `w <= -q/2`.
Checked on the thread path: `q = 0.5`, `w = -0.5`, identity to 4.6e-13 (median), `T_un < 0` everywhere.
Consequences: torsion matter (`w >= 1`) needs `q <= -2`, energy falling faster than `R^-2`, while every path has `|E| >= sqrt(2R|M_p - m_-|)`, which rises with `R`. A vacuum layer (`w = -1`) is allowed for `1/2 <= q < 2`.

## 3. Vacuum layer on the thread path (one Dirac field)

Thickness that carries the thread-path energy: `d = 7.7e-25 sqrt(R) m` (7.7e-36 m at `1e-22 m`, 7.7e-25 m at 1 m, 2.5e-23 m at 1 km, 1.3e-22 m at `r_g`). It is thin (`d/R <= 8e-14` for `R >= 1e-22 m`), at least one Planck length for `R >= 4.4e-22 m`, and adiabatic (`(d/c) d ln d/d tau < 1`, parameter `6.7e-23 m/R`) for `R >= 7.3e-23 m`. Several species scale `d` as `N^(1/3)`.

## 4. A `w = -1` layer launched on the thread

Integrating `sigma_dot = -T_un` for the radiation child from `R0 = 1e-25 m` with `beta_+ = 0` (the thread). The layer leaves the thread at once:
- `b/B -> 1`, `q = 1.75`, `Rdot = 2 pi |sigma| R` (ratio 1.000), a neck expanding at near light speed.
- At `R = r_g`: `|E| = 1.2e37 M_p` (`1.5e-3` of the child's bare mass `7.9e39 M_p`), shell proper time `1.3e-40 s` (`3e4` torsion dynamical times), child time `r_g/2` (49 microseconds, `R = 2t`).
- `|sigma|` falls from `4.3e38` to `1.7e31 m^-1`. The vacuum-layer thickness for that `sigma` is `0.015 l_P` at the start, one Planck length at `R = 5.1e-4 m`, and `4.5 l_P` at `r_g`.
- `|E|` at `r_g` scales as `R0^(-5/4)`; the value above is for `R0 = 1e-25 m`.

## Reading

The vacuum layer has the right sign pattern and, on the thread path, a comfortable thickness. Its ratio is `-1`, twice the `-1/2` the thread needs, and that factor does not stay small: the layer cannot hug the thread, and the neck opens at near light speed (child time 49 microseconds) with a negative energy far above the thread's minimum.

## Not accounted for

- What bounds the layer (the plates), and curvature corrections of order `d/R`.
- The first `5e-4 m` of the neck, where the layer is thinner than `l_P`, and the Planck-scale start (`R0` is a free choice here).
- A path between the thread and the near-light-speed run, with a layer of ratio between `-1` and `-1/2` (needs anisotropic stress or a flux that is not child fluid at rest).
- Whether `1.2e37 M_p` of negative energy at `r_g` is an acceptable fraction of the child's bare mass (it is `1.5e-3` of it).
