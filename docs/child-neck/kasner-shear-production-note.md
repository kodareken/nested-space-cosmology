# Pair production by shear in a stretching (Kasner) collapse, and where the stretch stops

7 Oct 2026, version 2. Scripts: `kasner_shear_production.py` (F, C, energy), `kasner_stopping_scales.py` (three scales, stage 2),
`kasner_shear_production_figure.py` (figure). Numbers: `kasner_shear_production.json`, `kasner_stopping_scales.json`.

## Update (version 3)

Counts of pairs below were compared with the paper's `N_need`, which counts fermions. Fermions are `2 x pairs`, so every ratio to `N_need`
in this note should be doubled (`N_br`: `0.73` of `N_need`, range `0.41 .. 1.07`; thread length needed `1.37`, range `0.94 .. 2.44`).
See `child-energy-floor-and-thread-length-note.md`, which also computes the no-production baseline and the energy floor.

## Correction to version 1

The pair-density coefficient in version 1 was a factor 2 too low. The Dirac field gives `2 F` pairs per momentum
(the 4x4 trace in the script confirms it: `0.2489 = 2 x 0.1245`), and the factor 2 was not carried into `C`. With both signs of `k_x`
and the checked density of states `d^3k = 4 pi lambda^2 r^-3 dlambda dr / tau0^3`:

`C = n (c tau)^3 = Lambda_2 u_c^2 / (2 pi^2)`, `Lambda_2 = 0.01424`. Version 1 had `4 pi^2`.

| u_c | 0.5 | 1 | 2 |
|---|---|---|---|
| `C`, version 1 | `9.0e-5` | `3.6e-4` | `1.4e-3` |
| `C`, corrected | `1.8e-4` | `7.2e-4` | `2.9e-3` |

Everything linear in `C` doubles. Version 1 also assumed each pair carries energy `2 hbar / tau`. The exact weighting from the same
two-level solution (`|q| tau = lambda u^(1/3) sqrt(1+u^2)`) gives `eps tau^4 = (2/pi^2) Lambda_3 J(u_c)` with `Lambda_3 = 0.00846`
(cut at `lambda <= 5`; the `F` tail beyond is finite-time dressing, 56% of the uncut `Lambda_3`) and `J(u_c) = integral_0^u_c u^(4/3) sqrt(1+u^2) du`:

| u_c | 0.5 | 1 | 2 |
|---|---|---|---|
| `eps tau^4` (units `hbar c`) | `1.6e-4` | `9.1e-4` | `6.5e-3` |
| mean pair energy `x tau` | `0.86` | `1.26` | `2.24` |

## Counts for a 10 solar-mass hole, one Dirac field (corrected)

`N(T) = 6 pi C (r_g / cT)^2`, thread length `2M` at the reference time, `L = 1`.

| Quantity | `u_c = 1` (range over `u_c`) |
|---|---|
| Pairs at `cT = l_P` | `4.5e76` (`1.1e76 .. 1.8e77`) |
| as fraction of the paper's needed `1.8e78` | `2.5%` (`0.6% .. 10%`) |
| as fraction of `S_BH = 1.05e79` | `0.43%` (`0.1% .. 1.7%`) |
| `cT` where `N = 1.8e78` | `0.158 l_P` (`0.079 .. 0.32`) |
| Pairs at the paper's stage-1 bounce density `3.9e74 J/m^3` (`cT = 1.65e-26 m`) | `4.4e58` (`2.6e58 .. 6.5e58`); with all Standard Model fermions (22.5 Dirac, from memory) `2.1e59`: 19 orders short |
| With all fermions at `cT = l_P` | `1.0e78` = `0.56` of needed (`0.14 .. 2.3`) |
| Energy of the pairs at `cT = l_P` | `6.3e38` solar masses (`1.1e38 .. 4.5e39`), scales as `M^2` |

## Three scales (all independent of `M_BH`)

`x = c tau / l_P`. Pair gas: `n = C / tau^3`, `eps = e / tau^4`.

1. **Shear energy used up.** Vacuum Kasner has `theta^2 / 3 = sigma^2 = 1 / (3 tau^2)`. The pairs back-react when `kappa eps = sigma^2`, i.e. `x_br^2 = 24 pi e`.
2. **Count reached.** `N(x_N) = N_need`.
3. **Torsion contact equals the pairs' energy.** `alpha n_f^2 = eps` with `alpha = kappa (hbar c)^2 / 32` and `n_f = 2 n` (fermions plus antifermions): `x_t^2 = pi C^2 / e`.

| u_c | 0.5 | 1 | 2 |
|---|---|---|---|
| `x_br` | `0.108` | `0.262` | `0.698` |
| `x_N` | `0.079` | `0.158` | `0.316` |
| `x_t` | `0.026` | `0.042` | `0.064` |

Order along the collapse: shear energy used up first, then (if production simply continued) the count would reach the paper's number,
then the torsion contact term would catch up. The three are within a factor of about 10 of each other.

## What the shear energy can supply

At `x_br`: `N_br = R^2 L / (4 <E tau>)` with `R = r_g / l_P` and `<E tau> = e / C`. It does not depend on the species count or on `C` separately.
The paper's own requirement is `N_need = sqrt(8/27) R^2 = 0.544 R^2` (algebra from its Eqs. 1-4; `m_n` cancels; checked to 13 digits against the saved check).
Both counts are `(r_g / l_P)^2` times a coefficient. The exponent agrees by construction; the coefficients differ:

| u_c | 0.5 | 1 | 2 |
|---|---|---|---|
| `N_br` | `9.7e77` | `6.6e77` | `3.7e77` |
| `N_br / N_need` | `0.53` | `0.37` | `0.21` |
| `N_br / S_BH` | `0.093` | `0.063` | `0.036` |
| thread length needed for `N_br = N_need`, in units of `2M` | `1.9` | `2.7` | `4.9` |

## Stage 2: isotropise at `x_br`, then torsion bounce

After `x_br` take isotropic radiation with `N` conserved (`eps ~ a^-4`, `n ~ a^-3`); the bounce is at `alpha n_f^2 = eps`:
`a_b / a_br = x_t / x_br`, `eps_b = e / x_t^4 = e^3 / (pi^2 C^4)`, `E_b = E_br x_br / x_t`.

| One Dirac field, `u_c` | 0.5 | 1 | 2 |
|---|---|---|---|
| `eps_b` (units `eps_P`) | `358` | `280` | `394` |
| `E_br` (solar masses) | `8.5e40` | `3.5e40` | `1.3e40` |
| `E_b` | `3.6e41` | `2.2e41` | `1.4e41` |
| `E_b` / paper's `M_b` (`1.7e31`, from memory of Eq. 9 scaling) | `2.1e10` | `1.3e10` | `8.5e9` |
| `E_b` / paper's stage-1 energy `N_need m_n c^2 = 1.5e21` | `2.3e20` | `1.4e20` | `9.4e19` |

With 22.5 Dirac species: `x_br = 0.51, 1.24, 3.3` (the backreaction starts at or above the Planck length for `u_c >= 1`),
`eps_b = 15.9, 12.4, 17.5 eps_P`, `E_b = 7.5e40, 4.5e40, 3.0e40` solar masses. `N_br` does not change.
The paper's stage-2 bounce density is `15.4 eps_P`. `eps_b` scales as `1 / species`, so this match is a coincidence at the level of the species count.

## Mismatches, sized, with what each would need

1. Count: `N_br` is `0.2 .. 0.5` of `N_need`. A thread length of `1.9 .. 4.9` times `2M`, or a coefficient from the missing pieces, closes it.
2. Energy: `E_b` is `~1e10` times the paper's `M_b` and `~1e20` times its stage-1 energy. At the bounce the effective energy
   `eps - alpha n_f^2` is zero by construction, so the matter energy `eps V` is balanced by the torsion contact term. The paper's `M_b` may be defined differently from `eps V`; not re-derived here.
3. Geometry: the thread at `cT = l_P` has transverse radius `2.6e-22 m` and length `3.2e17 m`. Its volume is `6 S_BH` Planck volumes,
   and the paper's child is a ball of radius `~4e-12 m`. The number of causal patches (size `c tau`) in the thread is of order `S_BH`.
4. Stage 1: the paper's neutron-like bounce (`n ~ a^-3`, stiff) is not produced; the produced quanta are massless radiation.

## Limits

One Kasner epoch (BKL oscillation ignored). Thread length set by hand. Counting cutoff `u_c` (factor 4 in `C`). Planck-length stop imposed by hand.
`eps = e / tau^4` comes from the sharp `u_c` cut and the adiabatic energy `2|q|`. Locally flat sphere, massless field, no torsion term inside the field evolution.
Single epoch backreaction criterion `kappa eps = sigma^2` has an O(1) ambiguity. Species count from memory.
