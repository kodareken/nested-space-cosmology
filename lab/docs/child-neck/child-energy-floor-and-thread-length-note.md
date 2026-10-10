# Child energy, the no-production baseline, and what sets the thread length

7 Oct 2026. Script `child_energy_floor.py`, data `child_energy_floor.json`, figure `child_energy_floor.png`.
Follows `kasner-shear-production-note.md` (v3 carries a pointer here).

## 1. Fermions, not pairs

The paper's `N_need = 1.8e78` counts fermions. The shear count in the earlier note was in pairs. Each pair is a fermion plus an antifermion
and the torsion term counts both, so the comparable number is `2 N_pair`. This is an error in my own comparison, now corrected (10 solar-mass hole, thread length `2M`):

| u_c | 0.5 | 1 | 2 |
|---|---|---|---|
| fermions at `x_br` | `1.9e78` | `1.3e78` | `7.5e77` |
| fermions / `N_need` | `1.07` | `0.73` | `0.41` |
| thread length needed for `N = N_need` (units of `2M`) | `0.94` | `1.37` | `2.44` |
| exterior advanced-time span needed | `0.09 ms` | `0.14 ms` | `0.24 ms` |
| `x` where the fermion count reaches `N_need` | `0.11` | `0.22` | `0.45` |

The count mismatch is now a factor of order one (`0.4 .. 1.1`), and the two scales `x_br` (`0.11, 0.26, 0.70`) and `x` of the count are within a factor of about 1.5.

## 2. What sets the thread length

Inside a Schwarzschild hole the surface `r = const` is a homogeneous cylinder `R x S^2` along the Killing coordinate `t`.
In ingoing coordinates `v = t + r*`, and `r* -> 0` as `r -> 0`, so at the singularity `t = v`: the thread length is the span of advanced time
over which the hole has been fed. One end is cut by the collapsing matter itself. The thread length, and with it `N` and `E`, are linear in that span.
`2GM/c^3 = 9.85e-5 s` for 10 solar masses. The needed span, `0.09 .. 0.24 ms`, is the order of the collapse time of a stellar core
(not computed; an order-of-magnitude remark).

Rate at `u_c = 1`: `1.35e82` fermions per second of infall for a 10 solar-mass hole. If each child holds `N_need` fermions
(an added assumption), the number of children is `K = 0.73 x (advanced-time span) / (2GM/c^3)`, with a coefficient that does not depend on `M`.

## 3. No-production baseline (computed)

`N = M_BH / m_n` neutrons, free degenerate neutron gas (full relativistic EoS), spin-torsion `eps_eff = eps - alpha n^2`, `alpha = kappa (hbar c)^2 / 32`.
Turnaround radius solved from `(kappa/3) eps(n0) a0^2 = 1` (gives `a0 / r_g = 1.03` for 10 solar masses). The collapse stays real (`adot^2 > 0`) down to the bounce.

| Quantity | Value |
|---|---|
| bounce density `n_b` | `1.20e105 /m^3 = 5.08 l_P^-3` |
| bounce energy density | `20.2 eps_P` (the paper's stage 2: `15.4`) |
| energy per fermion `alpha n_b` | `3.99 E_P = 7.8e9 J` |
| child energy / parent rest energy | `5.18e19 = 3.99 m_P / m_n`, independent of `M_BH` |
| child energy, `M_BH = 3, 10, 30` solar masses | `1.6e20, 5.2e20, 1.6e21` solar masses (linear in `M_BH`) |
| bounce radius, 10 solar masses | `1.3e-16 m` |

The energy of the child comes from compression work on the degenerate gas. The paper's `M_b ~ M_BH^2` needs production; the baseline gives `M_BH^1`.

## 4. Floor: energy of a child is set by its fermion number

At `eps = alpha n^2` with degenerate fermions, each carries `3.99 E_P`, so `E_child >= 3.99 E_P x N` whatever made the fermions.

| Child | `N` | energy (solar masses) | per fermion |
|---|---|---|---|
| no production, 10 solar masses | `1.2e58` | `5.2e20` | `3.99 E_P` |
| shear-produced, stage 2, `L = 1` (`u_c = 1`) | `1.3e78` | `2.2e41` | `14.8 E_P` |
| floor at `N_need` | `1.8e78` | `7.9e40` | `3.99 E_P` |
| paper, `N_need m_n c^2` | `1.8e78` | `1.5e21` | `1.5e-10 J` |
| paper, recorded `M_b` (from memory of Eq. 9 scaling) | `1.8e78` | `1.7e31` | `1.7 J` |

The recorded `M_b` is `4.7e9` below the floor; `N_need m_n c^2` is `5e19` below. So the `~1e10` energy mismatch in the shear note is not specific to shear:
it is the gap between the paper's `N_need` and its mass. For radiation-like content the energy per fermion falls as `1/a` after the bounce, so a child that has
expanded by a factor `~5e9` from its bounce would have the paper's `1.7 J` per fermion. If so, the paper's `M_b` is read at a later epoch than the bounce. This is a candidate, not a result.

## Limits

Degenerate Fermi gas with `g = 2`; other species or statistics change `3.99` by order one, not by 10 orders. Paper's `M_b` is a recorded value whose prefactor I did not re-derive.
Thread length, `u_c` counting and the single Kasner epoch as in the earlier note. The child count `K` assumes a fixed fermion number per child.
