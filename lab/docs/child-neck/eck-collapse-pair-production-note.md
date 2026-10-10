# Pair production in the Einstein-Cartan black-hole universe: a self-consistency check

7 Oct 2026. Script: `eck_collapse_consistency_check.py`; numbers: `eck_collapse_consistency_check.json`.
Source model: Poplawski, arXiv:1110.5019 (Eqs. 1-4, closed FLRW with spin-torsion term).

## What was tested

The paper assumes the fermion number grows as `N ~ a^-3` during the collapse, through
gravitational pair production, so that stiff matter obeys `eps ~ a^-6`. This note takes that
assumption as given and asks whether the dynamics it implies would allow the production it needs.
It uses the paper's own Eq. (4), by step-by-step substitution. No parameter is fitted.

## Results (M_BH = 3, 10, 30 solar masses; identical where the entry has no mass dependence)

| Quantity | Value |
|---|---|
| Peak expansion rate `H_max` | `7.77e23 /s`, independent of `M_BH`, reached at `a = 1.12 a*` |
| Neutron Compton frequency `omega_n = m_n c^2 / hbar` | `1.43e24 /s` |
| `H_max / omega_n` | `0.544` |
| Suppression `exp(-pi omega_n / H_max)` at that single instant | `3.1e-3` (rough) |
| Fermions needed at the bounce, `M_BH = 10` | `1.8e78` (starting from `1.2e58`) |
| Upper bound on excitable pairs, `(2/3) adot_max^3` | `2.8e39` |
| Gap | `38.8` orders of magnitude, same for every `M_BH` |
| Number density at the bounce | `7.7e84 /m^3` |
| Energy per fermion implied by the scalings | `5.0e-11 J` (`= m_n c^2 / 3`) |
| Fermi energy per fermion at that density, `hbar c k_F` | `1.9e3 J` |
| Minimum Fermi energy density at that density vs `eps(bounce)` | `1.1e88` vs `3.9e74 J/m^3`, `13.5` orders |

## What does not follow, on its own terms

1. **The production law.** Within Eq. (4), the expansion rate never reaches the neutron
   Compton frequency. For mass `m` the standard adiabatic argument suppresses production when
   `H << m c^2 / hbar`. The mode-count bound caps the pairs at about `1e39`, against the `1e78`
   that `N ~ a^-3` needs. So `N ~ a^-3` is an extra assumption, not a consequence of the
   mechanism the paper names.
2. **The bounce state.** The scalings fix the energy per fermion at `m_n c^2 / 3`, while the
   number density `7.7e84 /m^3` forces a Fermi energy per fermion of order `1e3 J`. The
   stage-1 bounce state violates the Fermi bound by 13.5 orders. The nonrelativistic-neutron
   remark does not hold at that density (`hbar k_F / m_n c ~ 1e13`).
3. **Light versus stiff.** Fields light enough to be produced (`m c^2 / hbar <~ H`) have
   `p <= eps/3`, so work done on them gives `E ~ a^-1` at most, not `a^-3`. The stiff
   `p = eps` needed for `m ~ a^-3` comes from interacting nuclear matter, which is heavy.

## What was not tested

- Production by shear in an anisotropic (Bianchi-type) collapse. The mode-count bound covers it
  as an upper bound, but a rate calculation is a separate step.
- The paper's second stage (annihilation into radiation, bounce at `15.4 eps_Pl`).
- The paper's earlier work on the matter-antimatter and dark-energy readings.

## Conventions

Eqs. 5-6 agree with direct substitution only if `m_Pl` in them is the non-reduced Planck mass.
The constant in Eq. 8 does not reproduce. Scalings (`M_b ~ M^2`, `a_b ~ M^(2/3)`) are
trusted; prefactors are not.

## Where a new assumption can enter

- A different mass scale for the produced fermions in the daughter room (an inheritance relation).
- A stiffening interaction that does not require heavy particles.
- A no-production baseline: with `w < 1` matter and conserved particle number (`n ~ a^-3`),
  `alpha n^2 ~ a^-6` eventually exceeds `eps ~ a^-3` to `a^-4`, so a bounce still occurs while the
  mass gain `M_b ~ M^2` does not. This is argued from the scalings here, not yet computed.
