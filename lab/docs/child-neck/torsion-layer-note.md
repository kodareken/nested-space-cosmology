# Can torsion matter above the bounce density supply the neck's negative-energy layer?

**Scope correction 7 Oct 2026 (own error): the lifetime ratios in section 3 (`3e33` to `3.7e40`) hold for the thread path only, the slowest and least energetic path; they are not a property of torsion matter. Along any neck `|E| tau ~ R^2`, so a more energetic path is shorter-lived (`casimir-layer-note.md`, section 0). The pressure-sign result of section 1 stands, and a path-independent form of it is in `casimir-layer-note.md`, section 2.**

7 Oct 2026. Script `torsion_layer_eos.py`, data `torsion_layer_eos.json`, figure `torsion_layer_eos.png`. Continues `israel-matching-note.md`, section 4. Needs `israel_matching.py`, `junction_neck.py` and the restored `child_*` files.
Fluid: degenerate neutron gas with the contact term, `eps_eff = eps_F(n) - alpha n^2`, `p_eff = p_F(n) - alpha n^2` (pressure formula `p = n d eps/dn - eps` checked to 6e-10). The bounce density is `n_b = 1.20e105 m^-3`, `eps_b = 20.25 eps_P`, and `eps_eff = 0` at `n_b`.
A perfect-fluid layer has surface energy `sigma = eps_eff d` and surface pressure `P = p_eff d`, so `w = P/sigma = p_eff/eps_eff`.

## What the junction asks for

Negative surface energy with `w = -1/2` on the thread path (for `R << r_g`), and `w` falling from `+4.3` to `-0.5` along the mass-matched path (`israel-matching-note.md`, section 3). So the layer needs negative energy and **positive** pressure.

## 1. The fluid cannot give that sign pattern

`eps_eff < 0` means `alpha n^2 > eps_F`. Since `p_F <= eps_F/3`, `p_eff < eps_F/3 - alpha n^2 < -(2/3) eps_F < 0`. So whenever the energy is negative the pressure is negative too, and `w > 0`.
Numbers: the largest `p_eff` with `eps_eff < 0` is `-0.667 eps_b` (at `n -> n_b` from above). `w` is `101` at `1.01 n_b`, `2.13` at `2 n_b`, `1.35` at `5 n_b`, `1.03` at `100 n_b`, and tends to 1. The null energy condition `eps_eff + p_eff < 0` holds for every `n >= n_b` (`-0.67 eps_b` at `n_b`, `-4.6 eps_b` at `2 n_b`).
`w = -1/2` is reached only at `n = 0.414 n_b`, where `eps_eff > 0`: the wrong sign of the energy.

## 2. Where the mass-matched path allows a torsion layer

`w >= 1` on the mass-matched path holds for `a` between `0.1925` and `0.2175 a_exit` (`R` from `0.111` to `0.131 r_g`), a window of `0.51 microseconds`. There a torsion layer would have `n = 1.35 n_b` (`w = 4`) to `2.15 n_b` (`w = 2`). Beyond it `w < 1` and no torsion layer fits.

## 3. Thickness and lifetime

- Thickness of a thin layer carrying the thread-path energy at `2 n_b`: `3.8e-36 m` (0.23 Planck lengths) at `R = 1e-23 m`, `1.2e-70 m` at `R = 1 m`, `2.2e-77 m` at `R = r_g`. It is at least one Planck length only for `R < 3.8e-24 m`, and thinner than `R` only for `R > 1.1e-28 m`. (At `1.5 n_b`: `7.5e-24 m` and `1.6e-28 m`; at `5 n_b`: `7.6e-25 m` and `4.1e-29 m`.)
- Filling the whole neck volume instead needs a density excess `delta = n/n_b - 1 = 7.9e-70` at `R = 1 m` (`7.9e-20` at `1e-20 m`, `5.1e-81` at `r_g`; `eps_eff = -0.667 delta eps_b`).
- Dynamical time of matter at `eps_b`: `4.1e-45 s`. The thread path needs the layer for the parent proper time up to the radius: `3.1` dynamical times at `R = 1e-22 m`, `3.1e6` at `1e-18 m`, `3.1e33` at `1 m`, `1.0e38` at `1 km`, `3.7e40` at `R = r_g` (155 microseconds in total). In one dynamical time the neck reaches `4.7e-23 m`. The mass-matched window lasts `1.2e38` dynamical times.

## Reading

Torsion matter above `n_b` has the right sign of energy and the right order of size, and violates the null energy condition (the property a flare-out neck needs). It lacks the positive pressure the junction needs, it would have to last `1e33` to `1e40` of its own dynamical times, and as a thin layer it is sub-Planckian for `R > 4e-24 m`.

## Not accounted for

- A source with negative energy and positive tangential pressure: an anisotropic stress, or a layer re-supplied continuously.
- Candidate with that sign pattern, not computed here: Casimir-type vacuum stress of the Dirac field across the neck (parallel plates have negative energy density and positive tangential pressure, `w = -1`).
- A treatment of the neck at the Planck scale (thickness and radius below `l_P` for the first `1e-23 m`).
