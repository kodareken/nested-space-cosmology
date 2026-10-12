# Composite neck layer: Dirac vacuum layer plus torsion matter, matched to the thread path

**Partly superseded 7 Oct 2026 (`balanced-layers-note.md`): the 125-field floor and the cancellation of `1e7` to `1e41` below came from leaving the normal stress unbalanced. With a vacuum-like fill that balances it, one field suffices for `R >= 1.1e-22 m` and the layers do not cancel; the tuning moves to the torsion layer's density.**

7 Oct 2026. Script `composite_layer.py`, data `composite_layer.json`, figure `composite_layer.png`. Continues `casimir-layer-note.md`. Needs `torsion_layer_eos.py`, `israel_matching.py`, `junction_neck.py` and the restored `child_*` files.
Geometric units. Target (thread path): surface energy `-s`, `s = |E|/(4 pi R^2)`, `|E| = sqrt(2R(M_p - m_-))`, surface pressure `+s/2` (`w = -1/2`).

## Construction

Vacuum layer (`m` Dirac fields), thickness `d_v`: `sigma_v = -x`, `P_v = +x`, `x = m (7 pi^2/2880) l_P^2/d_v^3`.
Torsion matter, thickness `d_m`, density `n`: `sigma_m = eps_eff(n) d_m`, `P_m = p_eff(n) d_m`.
Totals: `sigma_m = x - s`, `P_m = s/2 - x`, so `w_m = (x - s/2)/(s - x)`.
- **Branch A** (`x > s`): `sigma_m > 0`, `n < n_b`, `w_m < -1`. Torsion matter below `n_b` has positive energy and negative pressure and cancels the excess vacuum energy.
- **Branch B** (`3s/4 < x < s`): `sigma_m < 0`, `n > n_b`, `w_m > 1`.
Thickness `d_m = |sigma_m|/|eps_eff(n)|` with `G eps/c^4 = (eps/eps_b) 20.25/l_P^2`. The composite satisfies `w <= -q/2` of `casimir-layer-note.md` as the thread path does (`-1/2` against `-1/4`).

## Results (10 solar-mass parent, N_need radiation child)

- **Both layers can be thick enough.** Branch A with `n` near `0.55 n_b` gives a common thickness (the thinner of the two, maximised over `x`) that does not depend on `R` for `R >= 1e-22 m`: `d* = (m (7 pi^2/2880)/(eps_eff,max 20.25))^(1/4) l_P = 0.299 m^(1/4) l_P`, with `eps_eff,max = 0.148 eps_b` for `n < n_b`. That is `0.30 l_P` for one field, `0.65 l_P` for 22.5 fields (the paper's count), `1.0 l_P` for 123 fields and `1.7 l_P` for 1000.
- **Species floor.** Reaching one Planck length for both layers needs `m >= 125` Dirac fields for `R >= 1e-22 m`. Below `1e-22 m` the number rises as `R^(-3/2)`: `7.7e3` at `1e-24 m` and `7.0e6` at `1e-26 m`.
- **Cancellation depth.** At about 125 fields the vacuum and torsion layers nearly cancel: `x/s` is `16` at `R = 1e-22 m`, `1.4e7` at `1e-18 m`, `1.5e34` at `1 m`, `4.4e38` at `1 km`, `7.9e40` at `r_g`.
- **Why it works.** At `n = 0.55 n_b` the torsion fluid has `p_eff = -eps_eff` (a cosmological-constant-like fluid). The composite is two layers of that type with opposite signs of energy. The residual ratio `-1/2` comes from `w_m` differing from `-1` by `s/(2x)`, set by tuning `n` to that level (`1e-34` at `R = 1 m`).

- **Normal stress is not balanced.** The vacuum layer has normal pressure `3 rho_v` (plates), integrated `-3x`; the torsion layer is isotropic, integrated `P_m = s/2 - x`. The sum, `-4x + s/2`, is `-2.5 s` at `R = 1e-26 m`, `-63 s` at `1e-22 m`, `-6e34 s` at `1 m`, `-1.8e39 s` at `1 km` and `-3e41 s` at `r_g`. A thin-shell junction allows a normal stress only of the shell's own scale (`~ s/R` per unit thickness), so the composite needs a confining wall that carries `+4x` of integrated normal stress.

## Not accounted for

- The confining wall that balances the normal stress above, and what it adds to the surface energy and tangential pressure. (A positive vacuum-energy region with `p = -eps`, like the bag constant of the MIT bag model, is a candidate sign pattern; not computed.)
- The planar vacuum formula is used with curvature corrections of order `d/R` and without any interaction between the two layers; the coefficient depends on boundary conditions and on geometry (spherical shells differ from plates).
- The species floor of 125 is 5.6 times the paper's 22.5 fields.
- The tuning of `n` to `1e-34` to `1e-41` along the path, and whether the layers' internal dynamics can follow `s(R)`, were not checked. For comparison, the bounce itself needs `eps - alpha n^2` to vanish to about 38 digits.
- Only the thread path is matched. The Planck-scale start of the neck (`R < 1e-22 m`) is not described.
