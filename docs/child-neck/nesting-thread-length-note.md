# What sets the child's thread length, and why it is not the lever for nesting

7 Oct 2026. Script `nesting_thread_length.py`, data `nesting_thread_length.json`, figure `nesting_thread_length.png`.
Needs `child_energy_floor.py`, `child_matter_fraction.py` (turnaround formula) and the saved JSON files. Follows `second-quantised-splitting-note.md`.
Thread length `L` is in units of `2GM/c^3` (the infall span of the parent; `L = 1` is `0.1 ms` for 10 solar masses).

## Correction first

`child-expansion-and-matter-fraction-note.md` (version 2) said the generation map has an unstable fixed point at `1.8e11 kg` with lighter children below it. That is wrong:
`M' = N m_n` holds only where the matter term dominates the turnaround. For small `N` the radiation energy released at the bounce sets it, and with the full turnaround formula the child outweighs its parent
at every mass above the Planck scale (child = parent only at `0.65` Planck masses, `N_f < 1`; check at `1e-20` solar masses: child `1.1e-14`). The comparison of a child's net matter with the parent's baryon number, `eta <= M*/M`, is a number comparison and stands.

## Three routes to `L`, all order one

| Route | `u_c = 0.5` | `1` | `2` |
|---|---|---|---|
| infall span for the paper's `N_need` | `0.94` | `1.37` | `2.44` |
| cap from thermal entropy, `S_th = 4.20 N_f <= S_BH` | `1.29` | `1.88` | `3.35` |
| cap from entanglement entropy of the pairs `<= S_BH` | `2.6` | `3.8` | `6.8` |

The paper's `N_need = sqrt(8/27) R^2` child has `S_th / S_BH = 0.73` for any `u_c`. Closure of a closed radiation child (isotropic collapse until `c tau* = x a_max`, then Kasner):
`L_c = (pi/3) [x(2-x)]^(3/2) / x`, which is `0.29` for `x = 0.01`, `0.87` for `0.1`, `1.36` for `0.5`, `1.05` for `1`. So `L` is `0.3 .. 2` by three independent routes.

## Generation map (u_c = 1)

`N_f / N_need` is the same every generation: `0.22` (`L = 0.3`), `0.73` (`L = 1`), `1.0` (`L = 1.37`). The paper's count relation reproduces itself at every level.
The scale ratio between generations does not: `r_g(M_{n+1}) / r_g(M_n)`.

| | generation 0 to 1 | generation 1 to 2 |
|---|---|---|
| all matter, `L = 1` | `1.1e20` | `1.2e40` |
| radiation only, `L = 1` | `1.1e13` | `2.5e17` |

The repository's `Omega = 3.97` is a constant ratio; this map gives a ratio that grows with the mass. The thread length cannot repair it: to keep generation 2 inside the observable radius (`4.4e26 m`) `L` would have to be `5e-6` (radiation) or `2e-13` (matter).

## Hubble mass (numerical observation, three free parameters)

For every fermion surviving as matter, a parent of `63` solar masses (`52 .. 84` over `u_c`, `L = 1`; `96 .. 154` for `L = 0.3`) gives a child at turnaround of the Hubble mass `M_H = c^3 / (2 G H0) = 4.5e22` solar masses (`r_g = c/H0`, `H0 = 70 km/s/Mpc`), lifetime `pi/H0 = 44 Gyr`.
With net matter `eta = 6e-10` plus hidden-antimatter dust the parent is `0.9 .. 2.6e6` solar masses; radiation only, `2 .. 7e7`. This is an order-of-magnitude coincidence of the model's parameters, not a constraint.

## Not accounted for

What caps the mass ratio between a parent and its child. The ADM mass seen from outside the throat is the parent's `M`; the child's mass is measured by its own interior. Which of the two enters the child's own collapse (and so its `N_f`) is not settled here.
This is the quantity that would turn the growing ratio into something constant.

## Limits

Thread length from a closed radiation child is crude (the Kasner onset is a free parameter `x`). Thermal entropy is an upper bound for the produced state. One Kasner epoch, `u_c` as before.
