# Normal-stress balance of the neck layer, and the layer whose own ratio is -1/2

7 Oct 2026. Script `balanced_layers.py`, data `balanced_layers.json`, figure `balanced_layers.png`. Continues `composite-layer-note.md`. Needs `casimir_layer_check.py`, `torsion_layer_eos.py`, `israel_matching.py`, `junction_neck.py` and the restored `child_*` files.
10 solar-mass parent, N_need radiation child, thread path; `s = |E|/(4 pi R^2)` is the surface energy the thread needs, with surface pressure `+s/2`. Geometric units.

## 1. Balance condition

A static planar layer in an exterior vacuum has `T^zz = 0` at every point (conservation of the stress tensor). A constituent with normal pressure `p_n` therefore sits in a vacuum-like fill (`rho_b`, `p_b = -rho_b`, isotropic) with `p_n + p_b = 0`.
- Dirac vacuum slab: `rho_f = -x/d`, `p_perp = 3 rho_f`, `p_par = -rho_f`. The fill is `rho_b = -3x/d`. Surface energy `-4x`, tangential pressure `+4x`, ratio `-1`.
- Isotropic torsion fluid: the fill is `rho_b = p_eff`. Surface energy `(eps_eff + p_eff) d_m`, tangential pressure `0`: a dust-like layer, with negative surface energy where `eps_eff + p_eff < 0`, i.e. `n > 0.5443 n_b` (zero at `(2/3)^(3/2) n_b`; the slope `d h/d f = -0.726` with `h = (eps_eff + p_eff)/eps_b`; `h(n_b) = -2/3`, `h(2 n_b) = -4.64`).
A balanced mixture has `w = -(fraction of the surface energy in the vacuum-like layer)`. The thread's `-1/2` is half and half.

## 2. The balanced composite replaces the earlier one

Thread path: `4x = s/2`, so `x = s/8`, and `(eps_eff + p_eff) d_m = -s/2`.
- The earlier composite (`composite-layer-note.md`: 125 fields, cancellation `1e7` to `1e41`) left the normal stress unbalanced. Those two features were consequences of that and are superseded.
- Vacuum layer, one field: thickness `d_v = (8 C l_P^2/s)^(1/3)` (`C = 7 pi^2/2880`). At least one Planck length for `R >= 1.1e-22 m`; `95 l_P` at `1e-18 m`, `9.5e10 l_P` (`1.5e-24 m`) at `R = 1 m`, `1.7e13 l_P` at `r_g`. 22.5 fields: threshold `1.4e-23 m`. No cancellation between layers.
- Torsion (dust-like) layer: `d_m = (s/2)/(|h| eps_b)` in geometric units. It reaches one Planck length only if `n` is within a fraction `delta f = s l_P/(2 eps_b-geometric |h'|)` of `0.5443 n_b`: `7.6e-3` at `R = 1e-22 m`, `7.5e-9` at `1e-18 m`, `7.6e-36` at `1 m`, `2.4e-40` at `1 km`, `1.4e-42` at `r_g`. Elsewhere it is thinner: `1.2e-36 l_P` at `R = 1 m` for `n = 2 n_b`.
Any negative-energy, dust-like layer of surface energy `-s/2` would serve; the torsion fluid is one example, and its tiny required `s/2` against the Planck density is the reason for the tuning.
- Adiabaticity: the relaxation time `l_P/c` of the normal stress against the neck proper time `(2/3) R^2/|E|` is `4e-6` at `1e-18 m`, `4e-33` at `1 m`; balance applies for `R >= 2.6e-22 m`.
- Energy exchange: with both layers scaling as `R^(-3/2)`, the vacuum layer must gain `0.75 s Rdot/R` per area and the dust-like layer lose net `0.25 s Rdot/R` (total `T_un = -0.5 s Rdot/R`). The mechanism is not checked.

## 3. A layer whose own ratio is -1/2

Integrating `sigma_dot + (sigma + P) Adot/A = -T_un` with `P = -sigma/2` from thread initial data (`R0 = 1e-25 m`): the layer stays on the thread to `E/E_thread - 1 = 1.0e-4` at `R = 2.3e-5 m`, `1.3e-3` at `4.3e-3 m`, `1.7e-2` at `0.79 m`, `0.22` at `145 m`, after which it leaves it. At `R = 0.9 r_g`: `|E| = 6.1 M_p`, `b/B = 0.92`, child time `1.16 G M_p/c^3` (57 microseconds), against `3.3` for the thread. It departs because the thread's own ratio drifts from `-1/2` (to `-0.37` at `0.39 r_g`).
A self-balanced layer with `w = -1/2` and no fill: an isotropic planar network of negative-tension lines (`T^zz = 0`, tangential stress `-lambda mu/2` per direction). With `g = G mu/c^4` the line length per area is `lambda = s/g`. The minimum tension for line spacing of one Planck length is `g = s l_P`: `0.22` at `R = 1e-22 m`, `2.2e-7` at `1e-18 m`, `2.2e-34` at `1 m`, `4e-41` at `r_g`. For Planck tension (`g = 1`) the spacing is `7 cm` at `R = 1 m`, with `172 m` of line in total.

## Not accounted for

- What the negative-energy vacuum-like fill is (energy density `-3x/d` for the slab, `p_eff` for the torsion layer, up to Planck density) and what holds it.
- The energy exchange between the two layers, and the Planck-scale start of the neck (`R < 3e-22 m`).
- For the line-network layer: what the lines are (a candidate stress type only, not a derivation).
- Self-gravity of Planck-thick layers (`G sigma d` is of order one at one Planck length), which the planar slab ignores.
