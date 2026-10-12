# After the bounce: how far the child expands, what it weighs, and what the antimatter fraction does

7 Oct 2026. Scripts `child_expansion_ladder.py` (full EoS expansion), `child_matter_fraction.py` (matter fraction, generation map, figure).
Data `child_expansion_ladder.json`, `child_matter_fraction.json`. Follows `child-energy-floor-and-thread-length-note.md`.

## Model

Closed child (`k = +1`), `N` fermions, torsion contact term `eps_eff = eps - alpha n^2`, adiabatic expansion at fixed `N`:
`adot^2 = (kappa/3)(eps - alpha n^2) a^2 - 1`. Bounce at `eps = alpha n^2` (`n_b = 5.08 l_P^-3`, `3.99 E_P` per fermion).
Degenerate neutron-gas equation of state with `m_n`; exact relativistic form above `x = 0.05`, series below. Turnaround where `(kappa/3) eps a^2 = 1`.
There `eps V = a c^4 / (2G)`: **the child's size equals its own Schwarzschild radius**, the same condition that fixes the parent's turnaround `a0 = r_g`.

Check: the baseline (`N = M_BH / m_n`, 10 solar masses) turns around at `a_max = 1.027 r_g`, matching the turnaround solved separately in the previous note.

## Results (`N` fermions; paper's `N_need = 1.8e78`; shear count with thread `2M` is `1.3e78`)

| Quantity | `N = 1.3e78` | `N = N_need` |
|---|---|---|
| `a_b` | `6.4e-10 m` | `7.1e-10 m` |
| size where each fermion carries `1.67 J` (paper's recorded `M_b`) | `3.0 m` (`a/a_b = 4.7e9`) | `3.3 m` |
| size where each fermion carries `2 m_n c^2` (non-relativistic) | `1.9e10 m` (`3.0e19 a_b`) | `2.2e10 m` |
| turnaround size, massive fermions | `3.3e24 m` | `4.5e24 m` |
| turnaround mass, massive fermions | `1.12e21` solar masses `= 1.00 N m_n` | `1.53e21 = 1.00 N m_n` |
| turnaround size / mass, massless (radiation) | `3.3e17 m`, `1.1e14` solar masses | `4.1e17 m` (analytic `4.08e17`), `1.4e14` |
| expansion available, `a_max / a_b` | `5e33` | `6e33` |

- After the bounce the child accelerates only briefly (`a / a_b < sqrt 2`); then it is radiation-like. No inflation.
- The `~5e9` expansion that would turn `3.99 E_P` per fermion into the paper's `1.7 J` is far below what is available, so the paper's recorded `M_b` is one point on a
  `1/a` curve (size about 3 m), not a special epoch I can identify. The mass `N x (energy per fermion)` is constant only after the fermions are non-relativistic.

## What the antimatter fraction does

Shear production in a free field is symmetric between particles and antiparticles. If pairs annihilate, only the net matter `eta N m_n c^2` survives as dust and the rest is radiation (falls as `1/a`).
Turnaround from `a^2 - r_m a - r_rad^2 = 0` with `r_m = 2G eta N m_n / c^2`, `r_rad^2 = (2G/c^4) N 3.99 E_P a_b`:

| net matter fraction `eta` | 0 | 6e-10 | 1e-8 | 1e-7 | 1e-6 | 1e-3 | 1 |
|---|---|---|---|---|---|---|---|
| turnaround mass, `N = N_need` (solar masses) | `1.38e14` | `1.39e14` | `1.46e14` | `2.35e14` | `1.5e15` | `1.5e18` | `1.5e21` |

Matter and radiation contribute equally at `eta = 9.0e-8`. Below that the child is radiation dominated and its turnaround mass is `1.4e14` solar masses, nearly independent of `eta`.
The baryon-to-photon ratio of the observed universe, `6e-10`, is below this crossing (`eta` here is net matter per fermion in the child, not per photon, so the comparison is of order only). The calculation has no source of asymmetry; it is where an antimatter mechanism would enter, and it moves the child's
late mass over seven orders of magnitude.

## A first generation map

With `N_f = c1 (r_g / l_P)^2`, `c1 = 0.3975` (shear production, `u_c = 1`, thread `2M`) and the turnaround mass above:

| generation | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| all matter (`eta = 1`), solar masses | `10` | `1.1e21` | `1.4e61` | `2e141` |
| radiation (`eta = 0`), solar masses | `10` | `1.1e14` | `2.8e31` | `4.5e54` |

Correction (same day): the sentence that stood here claimed `M_{n+1} = K M_n^2` has an unstable fixed point at `1.8e11 kg` with lighter children below it. That was wrong. `M_{n+1} = N m_n` holds only where the matter term dominates the turnaround; for small `N` the radiation energy released at the bounce sets the turnaround mass, so with the full turnaround formula the child outweighs its parent at every mass down to the Planck scale (the map crosses child = parent only at `0.65` Planck masses, where `N_f < 1`). The number comparison `eta <= M*/M` (child's net matter against the parent's baryon count) is unaffected. See `nesting-thread-length-note.md`.
From the second generation the Schwarzschild radius exceeds the observable universe (`4e64 m` for matter, `8e34 m` for radiation). So applying the same rule at every level cannot be what nests the spaces: something else (thread length set by a different infall, a cut-off, BKL oscillation) must change it.

## Limits

Degenerate-gas equation of state (all fermions at the bounce momentum scale; real distributions are broader by order one). Thread length `2M`, `u_c = 1`, one Kasner epoch as before.
Surviving matter is taken as neutrons; other Standard Model fermions have other masses and lifetimes. Annihilation is assumed complete; entropy and thermalisation are not modelled.
Paper's recorded `M_b` is from memory of its Eq. 9 scaling.
