# Which swaps of the throat keep the Dirac action, and what they do to matter and helicity

9 Oct 2026. Script `swap_symmetry.py`, data `swap_symmetry.json`, figure `swap_symmetry.png`. Continues `child-and-neck-status-note.md` (fork 4).

## Question

If the two sides of a throat are the same space seen from either side, a swap exchanges them. The Dirac field carries a charge density `n_V` (matter minus antimatter) and an axial density `n_5` (right-handed minus left-handed, the helicity-like quantity that the torsion vector-axial term `V.A` couples to). For each swap: does the Dirac equation survive with the same mass, does the `V.A` term survive, and what happens to `n_V`, `n_5` and the particle/antiparticle labels?

## Setup

Kruskal form `ds^2 = f(s)^2 (dT^2 - dX^2) - R(s)^2 dOmega^2`, `s = X^2 - T^2`, with `f`, `R` arbitrary even functions of `s`. Regions: I and III exteriors, II black-hole interior, IV white-hole interior. Three swaps of the `(T, X)` plane:

| Swap | Map | Exchanges |
|---|---|---|
| `Sigma_X` | `(T, X) -> (T, -X)` | I with III |
| `Sigma_T` | `(T, X) -> (-T, X)` | II with IV |
| `Sigma_XT` | `(T, X) -> (-T, -X)` | I with III and II with IV |

For each swap every spinor lift `psi'(x) = S psi(Rx)` (linear) or `B conj(psi(Rx))` (antilinear), with `S`, `B` one of the 16 products of Dirac matrices, is tested numerically for `D psi' = kappa S (D psi)(Rx)` at random points, with the Dirac operator built from the tetrad and spin connection of the metric above. `kappa = +1` keeps the mass term; `kappa = -1` is the massless-only chiral-rotation partner. Charge and axial densities of the mapped field are computed from the 4x4 algebra. The sign of the frequency flips for exactly one of {time reversal, antilinear}, and then the mapped solution is read as the antiparticle (extra `-1` from anticommutation).

## Checks

- Spin connection: local rotation of the tetrad (angle depending on `T, X, theta, phi`) gives `D'(S psi) = S D psi` to `3e-16`; the other sign of the spin generator gives `1`.
- Flat space (Cartesian tetrad): the same search finds the standard lifts of `P`, `T`, `C` and `CPT`, and the computed signs of `psibar psi`, `i psibar gamma5 psi`, `V^mu`, `A^mu` equal the standard C, P, T, CPT table (the table is quoted from memory, not re-fetched; the script asserts equality).
- Negative control: multiplying `f` by `(1 + T/4)` and `R` by `(1 + X/5)` leaves no lift for any swap.
- A second pair `f = 1/(1+s^2/3)`, `R = exp(s/4)+1` with a different test spinor returns the same lifts.

## Result (`kappa = +1`)

| Swap | Lift | Frequency flips | `n_V` | `n_5` | `V.A` | Chirality kept |
|---|---|---|---|---|---|---|
| `Sigma_T` | `T`, antilinear, `gamma^1 gamma^3` | no | `+` | `+` | kept | yes |
| `Sigma_T` | `C T`, linear, `gamma^1 gamma^2 gamma^3` | yes | `-` | `+` | flips | no |
| `Sigma_X` | `C P`, antilinear, `gamma^0 gamma^3` | yes | `-` | `-` | kept | yes |
| `Sigma_X` | `P`, linear, `gamma^0 gamma^2 gamma^3` | no | `+` | `-` | flips | no |
| `Sigma_XT` | `C R T`, linear, `gamma^0 gamma^1` | yes | `-` | `-` | kept | yes |
| `Sigma_XT` | `P T`, antilinear, `gamma^0 gamma^1 gamma^2` | no | `+` | `-` | flips | no |

Three swaps keep the Dirac equation, its mass and the `V.A` term together: `T` on `Sigma_T`, `CP` on `Sigma_X`, `CRT` on `Sigma_XT`.
- `II <-> IV` by `T`: matter density and axial density keep their sign. The two interiors carry the same sign of matter excess and of helicity.
- `I <-> III` by `CP` or by `CRT`: matter density and axial density both change sign. The other exterior carries antimatter where this one carries matter, and the opposite helicity.
- The same swaps without the charge conjugation (`P`, `PT`, `CT`) flip the sign of the `V.A` term. With a fixed coefficient `c` the splitting is `delta E ~ q c n_5`, so flipping `n_5` at fixed `c` and flipping `c` at fixed `n_5` are the same statement: either the two sides carry opposite `n_5`, or the coefficient has opposite sign on the two sides.

## What this gives for the hierarchy numbers

- The thread is a path with `R` growing from 0 to `r_g` inside the horizon, which is the time-reverse of an infalling geodesic. A `Sigma_T`-type relation therefore links the child side of the neck to the parent's infalling matter.
- If the child carries the `T` image of the parent's net baryon number `N0 = 1.187e58`, then `eta = N0 / N_need = 6.53e-21` per quantum, with the same sign as the parent's. The value needed for `eta = 6e-10` is `9.2e10` times larger. This is of the same order as the earlier bounds (`eta <= 8.9e-21`, `chi_h <= 8.9e-21`).
- If the two sides are related by `CP` or `CRT`, the other side carries `-N0`, and the sum over both sides is zero.
- With `eta = 6.5e-21` the net baryon number swept up by the neck (`0.65` to `1.8` times `N0` quanta) is about `1e-20 N0`.

## Not accounted for

- A source of net baryon number or helicity in the child about `1e11` times the inherited amount (a Sakharov-type process at `~160 MeV` and above is not modelled).
- Which swap relates the two sides of the NSC throat. The test says what each swap implies; it does not select one.
- The child is `5e19` times the parent's energy, so it is not the image of the parent's collapsing region. The statement above covers the inherited part only.
- Free fields on the Kruskal vacuum: no interactions, no rotation (a rotating parent adds `J` and the earlier `p_max`), no layer.
