# The Hehl-Datta splitting at the bounce: how much spin the child needs, and what the parent can give

7 Oct 2026. Script `torsion_splitting_asymmetry.py`, data `torsion_splitting_asymmetry.json`.
Source read in full: Popławski, "Matter-antimatter asymmetry and dark matter from torsion", Phys. Rev. D 83, 084033, arXiv:1101.4012v2 (5 pages).
Follows `matter-antimatter-and-time-note.md`. Working assumption: the paper's mechanism holds at the bounce.

## What the paper does

- Hehl-Datta equation `i gamma^a psi_:a = m psi - (3 kappa / 8)(psibar gamma5 gamma^a psi) gamma5 gamma_a psi` (`hbar = c = 1`, `kappa = m_Pl^-2`, reduced Planck mass).
- For the charge-conjugate field the cubic term has the opposite sign relative to the mass term (his Eqs. 12-13). In a constant background of spin-aligned fermions of number density `n`, fermion levels are `omega = m +- (3/8) kappa n`, antifermion levels `m -+ (3/8) kappa n` (his Eqs. 14-15), so the corrections are C antisymmetric.
- Scenario: equal numbers of heavy fermions `X` and antifermions are produced; `X` is effectively heavier, decays faster; strong-mode decays give nucleons, weak-mode decays give dark-matter antiparticles; total baryon number stays zero; the hidden antimatter is dark matter. Observed ratio of visible to dark matter is reproduced for dark particles of 4-5 GeV.
- Asymmetry from Fermi-Dirac occupation with `E_X = |p| + alpha kappa N`, `E_Xbar = |p| - alpha kappa N`: `Delta n / s ~ kappa T_f^2`, giving `T_f ~ 1e13 GeV` and `m_X = 14.4 TeV` for `Delta n / s = 0.9e-10`.
- The paper itself states: with field operators the Dirac Lagrangian is C symmetric, and nevertheless the classical Hehl-Datta equation is used to get energy levels. That step carries the whole splitting.

## Putting it at the torsion bounce

At the bounce `eps = alpha n_f^2` with `alpha = kappa / 32`, so `alpha n_f` equals the energy per fermion. The level shift of the paper is `delta = (3/8) kappa n p = 12 alpha n_f p`, `p` the net polarisation fraction.
For a thermalised massless Fermi gas `<E> = 3.151 T`, so `delta / T = 12 x 3.151 p = 37.8 p`, the same at every stage.
Equilibrium asymmetry (Fermi-Dirac, levels `p +- delta`), `eta = (n_Xbar - n_X) / (n_X + n_Xbar)`:

| `delta / T` | 0.01 | 0.05 | 0.2 | 0.5 | 0.9 |
|---|---|---|---|---|---|
| `eta` | `0.0091` | `0.046` | `0.18` | `0.43` | `0.67` |

Linear slope `0.912`, so `eta = 34.5 p` for `p << 0.03`.

| target `eta` | 6e-10 | 1.6e-8 (crossover with hidden antimatter) | 9e-8 (crossover, visible matter only) | 1e-3 |
|---|---|---|---|---|
| polarisation needed | `1.7e-11` | `4.5e-10` | `2.6e-9` | `2.9e-5` |

## What the parent can supply

Child fermions `N_f = c1 L R^2` (`R = r_g / l_P`, `c1 = 0.40`, range `0.22 .. 0.58`, `L` = thread length in units of `2M`). Parent angular momentum at most `J_ext = G M^2 / c`, i.e. `J / hbar = (M / m_P)^2 = R^2 / 4` (`8.3e77` for 10 solar masses).
Polarisation the parent's rotation (spin parameter `a*`) can carry: `p_max = 2 J a* / (N_f hbar) = a* / (2 c1 L)`, which is `1.26 a* / L` (range `0.86 .. 2.2`; numerically equal to the mean pair energy times `tau`). It does not depend on the parent's mass.
The spin angular momentum needed for `eta = 6e-10` is `1.4e-11` of the extremal value. So the parent's rotation exceeds the needed polarisation by about ten orders of magnitude for any `a*` above `~1e-10`.

## Hidden antimatter as dust

If the compensating antimatter is dark (4.5 GeV, the paper's 4-5 GeV), each net unit adds `(m_n + m_DM) / m_n = 5.8` units of dust. Then, `N = 1.8e78`:

| `eta` | 6e-10 | 1.6e-8 | 1e-6 | 1e-3 |
|---|---|---|---|---|
| turnaround mass (solar masses) | `1.4e14` | `2.2e14` | `8.9e15` | `8.9e18` |
| lifetime | `89 yr` | `174 yr` | `8.7e3 yr` | `8.7e6 yr` |

The matter-equals-radiation crossover moves from `9e-8` to `1.6e-8`.

## Mismatches, sized

1. The C-odd splitting needs net polarisation; shear production gives none (argued). The parent's angular momentum is the available source, with room of ten orders of magnitude.
2. The classical-to-quantum step of the paper. In a second-quantised mean field with a C-even axial current, the spin-dependent shift is the same for a particle and an antiparticle of the same physical spin (my argument, not computed). A C-odd contact term does exist in the same family: the vector-axial interaction that non-minimal (Holst or Freidel-Minic-Takeuchi) coupling generates (found by search, not read). A splitting proportional to net matter number times axial density follows from it, and it would amplify an initial imbalance.
3. `eta` is net matter per fermion in the child; the observed `6e-10` is per photon. Order-of-magnitude comparison only.
4. The paper's `kappa` is the reduced-Planck-mass `kappa = m_Pl^-2` with `hbar = c = 1`; this note uses `kappa = 8 pi G / c^4`, same quantity.

## Draft entry for docs/arxiv-studies.md (not written into the repository)

`arXiv:1101.4012` (Popławski, PRD 83, 084033). Helps: C-asymmetric Hehl-Datta levels `m +- (3/8) kappa n`; hidden antimatter equals dark matter with total baryon number zero. Clue: net matter fraction at the torsion bounce is `eta = 34.5 p`. Constraint: the classical-versus-operator C step; polarisation must come from the parent's spin.

## Limits

Thermal equilibrium at the bounce is assumed. The splitting formula is the paper's, applied outside its density range. Polarisation is a free parameter here.
Spin-aligned background of one species per the paper; the mixed medium of this note is not derived.
