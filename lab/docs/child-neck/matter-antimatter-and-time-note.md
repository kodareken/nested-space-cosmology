# Matter and antimatter in the stretched thread, the parent's books, and what the imbalance does for time

7 Oct 2026. Script `child_lifetime_entropy.py`, data `child_lifetime_entropy.json`, figure `child_lifetime_entropy.png`.
Working assumption (the user's): the imbalance is part of what makes the child exist and last. This note derives what that implies.

## 1. Is matter versus antimatter part of what drives the stretch?

- The stretch itself is the parent's geometry: Kasner shear `sigma^2 = 1 / (3 tau^2)`, fixed by the parent's mass. The pairs' energy enters through the constraint
  `kappa eps = sigma^2`, which ends the stretch at `x_br = 0.26` (range `0.11 .. 0.70`). Matter and antimatter contribute equally to `eps`, so the imbalance does not change this scale at leading order.
- The torsion term that makes the bounce is axial-axial, `alpha n_f^2`, with `n_f` counting fermions and antifermions. The axial current is even under charge conjugation (argued, not computed here), so this term also treats them alike.
- Standard-literature handle on a difference (found by search, abstracts only, not yet in `docs/arxiv-studies.md`):
  Popławski, Phys. Rev. D 83, 084033 (arXiv:1101.4012): the Hehl-Datta term is cubic in the spinor and changes sign relative to the mass term under charge conjugation, so fermions have
  higher energy levels than antifermions in torsion; a decay asymmetry sends baryon number to normal matter and the opposite to hidden antimatter that forms dark matter, with total baryon number zero.
  Non-minimal coupling (Freidel-Minic-Takeuchi, Holst/Immirzi terms) can add vector-axial contact terms, odd under C.
- Size of the effect at the bounce (my estimate): torsion energy per fermion `2 alpha n_b = 8.0 E_P`, kinetic energy per fermion `4.0 E_P`. A level splitting would be of order the polarisation times these.
- Shear production leaves no net polarisation (argued): a right-handed Weyl pair with particle at `k` has spin `+k_hat`, at `-k` has `-k_hat`, left-handed ones the opposite; the distribution is symmetric in `k`. So a net spin or chirality has to come from outside this calculation.

## 2. Parent's books

Pair production adds equal numbers of fermions and antifermions, so the parent's baryon number is unchanged. The exterior sees only `M` (the child's energy is balanced by negative gravitational energy). If baryon number is conserved, a child with net matter `eta N_f` needs
net `-eta N_f` somewhere: sibling children in the same thread, the other side of the throat, or a hidden sector. If instead the child inherits net matter from the parent:

`eta <= M* / M`, `M* = m_P^2 / (4 c1 m_n) = 1.8e11 kg`, `c1 = 0.3975`.

| parent mass | 3 solar | 10 solar | 30 solar |
|---|---|---|---|
| maximum `eta` | `3.0e-20` | `8.9e-21` | `3.0e-21` |
| net baryons available | `3.6e57` | `1.2e58` | `3.6e58` |

`eta = 6e-10` needs a parent of `3.0e20 kg` (`r_g = 0.44 micron`); `eta = 9e-8` (matter equals radiation) needs `2.0e18 kg` (`r_g = 3 nm`).

## 3. Time

Entanglement entropy of the pairs when one partner is traced out (for example partners on opposite sides of a throat): `h(F) = -F ln F - (1-F) ln(1-F)` per mode block, same density of states as the count.
`S / N_pairs = 4.11` nats per pair (`4.09 .. 4.13` for cuts at `lambda <= 2 .. 10`). At `x_br`, thread `2M`, `S = 0.26 S_BH` (`0.15 .. 0.38` over `u_c`).
`S` grows linearly with the infall span, equalling `S_BH` after `3.8` times `2GM/c^3 = 0.38 ms` for 10 solar masses (`2.6 .. 6.8`, i.e. `0.26 .. 0.67 ms`).
`S` grows monotonically toward the crunch (`~ 1/x^2`) with `eta = 0`: the arrow of time in this model does not need an imbalance.

Child lifetime, bounce to recollapse (proper time), `N = N_need = 1.8e78`:

| `eta` | 0 | 6e-10 | 9e-8 | 1e-6 | 1e-3 | 1 |
|---|---|---|---|---|---|---|
| lifetime | `86 yr` | `87 yr` | `174 yr` | `1.5e3 yr` | `1.5e6 yr` | `1.5e9 yr` |

Checks: radiation `2 a_max / c = 86.3 yr`, dust `pi a_max / c = 1.50e9 yr`, both reproduced by the quadrature.

So with the imbalance assumed necessary, its quantitative content here is the child's duration and its ability to hold structure:
`86 yr` of radiation-only time without it, up to `1.5 Gyr` with all fermions surviving.

## Limits

Spin-polarisation and C-odd splitting are not computed. The entropy assumes a partner is traced out; for the whole thread in a pure state the net entropy is zero. Closed-child lifetime uses the degenerate-gas bounce and complete annihilation. The baryon bound assumes baryon number is conserved.
