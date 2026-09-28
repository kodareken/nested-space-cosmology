# Controlled covariant transition model

## v0.7 retained baseline: closed de Sitter in Einstein-scalar gravity

This is the project’s first controlled covariant core sector. It turns a
limited part of the motivation into equations without claiming more than those
equations show.

**Status:** an exact toy-sector realization: global de Sitter spacetime in
closed FLRW slicing. It is **not** a black-hole exterior/interior embedding,
not a solution for stellar collapse, not a derivation of a child universe, not
a controlled Planck-curvature calculation, and not a derivation of H-SAT or
`Lambda = 3/r_s,parent^2`.

Version 0.4 keeps this exact core unchanged and adds the separate
[DSF-1/DST-1 synthesis](domain-settings.md). In particular, the unit choices
`c=hbar=M_Pl=L=1` used by machine benchmarks are conventions, not a result that
local `c` or any dimensional constant varies. The matter-like oscillatory
scalar and Kottler limits in DST-1 are controlled adjacent calculations, not
matter or a compact parent already present in this solved de Sitter core.

Version 0.7 adds the separate **GMF-1** matching/flux audit. It is an exact
restriction on one proposed direct bridge, not an embedding of this core:
same-`Lambda`, finite-`M` Kottler cannot Darmois-match directly to pure CCT-1
because their spherical Misner–Sharp masses differ by `M`. The canonical
record is [global-matching-flux-closure.md](global-matching-flux-closure.md)
and its deterministic reproduction command is `python3 scripts/reproduce_gmf.py`.

## Conventions and action

Use signature `(-,+,+,+)`, `c=hbar=1`, and reduced Planck mass
`M_Pl^2=(8 pi G)^-1`. Let `Phi` be a real scalar order parameter of mass
dimension one. Its potential has mass dimension four:

```text
S = integral d^4x sqrt(-g) [M_Pl^2 R/2 - (nabla Phi)^2/2 - V(Phi)].
```

The equations are second order:

```text
G_mn = M_Pl^-2 [nabla_m Phi nabla_n Phi
                 - g_mn (nabla Phi)^2/2 - g_mn V],
box Phi - dV/dPhi = 0.
```

There is no higher-derivative/Ostrogradsky mode. At the stationary state,
require `V'(Phi_0)=0`, `V''(Phi_0)>=0`, and `V_0=V(Phi_0)>0`; the scalar
kinetic coefficient and propagation speed are positive and unity. The tensor
kinetic coefficient is `M_Pl^2>0`. This is a positive-potential stationary
point that is linearly non-tachyonic in the potential: `V''>0` gives a strict
quadratic minimum, while `V''=0` requires a separate higher-order stability
analysis.

## Exact closed-FLRW sector

```text
ds^2 = -dt^2 + a(t)^2 [dchi^2 + sin(chi)^2 dOmega_2^2],
H^2 + 1/a^2 = [dot(Phi)^2/2 + V]/(3 M_Pl^2),
dot(H) - 1/a^2 = -dot(Phi)^2/(2 M_Pl^2).
```

At a bounce, `H=0` and `dot(H)>0` requires `V_b > dot(Phi)_b^2`. The exact
solution implemented in [`transition.py`](../src/recursive_horizons/transition.py)
is

```text
Phi=Phi_0,  V_0=3 M_Pl^2/L^2,
a(t)=L cosh(t/L),  H(t)=tanh(t/L)/L.
```

This exact solution contains no ordinary matter. If a separate `S_matter` is
added to the action, its stress tensor must be set to zero in this solved core,
with any vacuum-energy contribution absorbed into `V_0`; radiation, dust, or a
generic fluid would change the Friedmann solution.

It contracts for `t<0`, has `a(0)=L`, and expands for `t>0`. All implemented
invariants remain finite:

```text
R=12/L^2,
R_mn R^mn=36/L^4,
R_mnrs R^mnrs=24/L^4.
```

At scales where `Phi` is stationary or decouples, the gravitational sector is
exactly GR plus a cosmological constant. Ordinary matter may be coupled in an
extended solution, but is not present in this exact core. The EFT is controlled
only when curvature and mass scales are below its cutoff; calling `L` Planckian
does not by itself justify this low-derivative action.

## Null expansions and regions

The areal radius is `R_A=a sin(chi)`. With future null normals

```text
k_+ = partial_t + a^-1 partial_chi,
k_- = partial_t - a^-1 partial_chi,
k_+ . k_- = -2,
```

the expansions are

```text
theta_+ = 2[H + cot(chi)/a],
theta_- = 2[H - cot(chi)/a].
```

Future-trapped spheres have both expansions negative; future anti-trapped
spheres have both positive; normal spheres have opposite signs; and marginal
spheres obey

```text
theta_+ theta_- = 0  <=>  |a H| = |cot(chi)|.
```

At the bounce the equator is marginal and other non-polar spheres are normal.
After it, an anti-trapped region develops around the equator. An expanding
FLRW domain is not anti-trapped everywhere.

These CCT expansions are not the null expansions of a finite-mass parent. In
spherical symmetry their product is a local quasi-local identity,
`theta_+ theta_-=-4 nabla_a R nabla^a R/R^2`; it locates a marginal round
sphere in the declared geometry. It does not by itself provide a common null
congruence, an event horizon, a flux, or a causal path from a parent trapped
surface to a child anti-trapped surface. A horizon-regular global solution
must calculate both sides with one orientation and matching prescription.

## GMF-1: the direct smooth-match obstruction

For a non-null Darmois junction, the induced metric and extrinsic curvature
must agree. In spherical symmetry this includes the Misner–Sharp mass on a
smooth interface. For a Kottler region and the pure CCT-1 core,

```text
M_MS,Kottler = M + Lambda R^3/(6 G),
M_MS,CCT-1  = Lambda R^3/(6 G).
```

Thus, at the same `Lambda`, a finite `M>0` leaves a mass jump `M`; a direct
smooth Kottler-to-pure-CCT-1 match fails. This does **not** say that every
spherical collapse-to-expansion model fails. It says that a successful model
must supply an inhomogeneous transition region, a specified shell/layer, or
modified field equations and then solve their constraints.

Oppenheimer–Snyder is the relevant control: it smoothly matches a finite
closed dust-FLRW interior to Schwarzschild, but its classical collapse is
singular. A prescribed Israel shell can represent a distributional local
interface, including the required surface stress; it is not a finite-invariant
global completion, a topology change, a derived `Q` flux, or a stability
result. A null interface requires the distinct Barrabes–Israel formalism.
The GMF-1 audit does not identify dark energy as external, and it has no
implication for local `c`.

## Perturbative and observational scope

For the exact stationary scalar, the physical fluctuation has a positive
kinetic term and `c_s^2=1`, so there is no ghost or gradient instability. That
does not by itself guarantee bounded evolution of every infrared mode: for
light fields the effective frequency squared displayed in the perturbation
note can become negative over part of the closed core. Any bounded-mode claim
must specify the mass, harmonic, and finite interval being tested. Because
`dot(Phi)=0`, the usual single-clock comoving-curvature variable is not defined
at this exact solution. This baseline therefore does **not** derive `n_s`,
`A_s`, or the observed adiabatic spectrum. A specified rolling/reheating sector
must be added and propagated through the regular background first.

If the observable child is additionally postulated to remain globally closed
with `K=+1`, that strict branch has the sign-level dimensionless discriminator

```text
Omega_K = -1/(a H/c)^2 < 0.
```

It is undefined at the exact bounce because `H=0`. CCT-1 by itself does not
derive that global post-core topology. A robust positive-`Omega_K` inference
that survives the stated nonflat cosmological model and systematic checks
rules out this additional strict branch. An observational result statistically
consistent with zero does not confirm it because later expansion can make its
magnitude arbitrarily small. The discrete spectrum becomes quantitative only
after the post-bounce history and initial state are independently fixed.

## What remains required

This sector must not be identified with the Recursive Horizons bridge until a
horizon-regular spherical solution connects an actual trapped interior to it,
constraints remain finite under refinement, shear and anisotropy are
controlled, a gravitational information map is defined, and perturbations
supply an independently testable spectrum or waveform. In v0.7 that solution
must also pass GMF-1B: derive the action/layer, global causal structure,
quasi-local flux/conservation law, energy conditions, and radial/bulk stability
rather than treating the GMF-1A obstruction or a prescribed shell as success.
