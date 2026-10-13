# Child and neck chain: what is computed, what is superseded, what is not yet accounted for

7 Oct 2026. One page of record for the chain that began with the shear-production test and ended with the neck-layer notes. Every number is taken from the note named in the last column.
Setting: parent of 10 solar masses (`r_g = 29.5 km`, `G M_p/c^3 = 49 microseconds`), one Dirac field, `u_c = 1` unless a range is given. The working assumption throughout is the user's picture: a stretched collapse that turns inward, with the child's space joined to the parent's through a neck. Nothing here grades that picture; mismatches are listed as pieces not yet accounted for.

## 1. Recorded exploratory results (read with the corrections below)

| Step | Result | Depends on | Note |
|---|---|---|---|
| Pair production by Kasner shear | `n (c tau)^3 = C = 7.2e-4` (`1.8e-4 .. 2.9e-3`); exact pair energy `x tau = 1.26` (`0.86 .. 2.24`) | single Kasner epoch, `u_c` cut, thread `2M` | `kasner-shear-production-note.md` (v3) |
| Where the stretch stops (units of `l_P`) | shear energy used up `0.26` (`0.11 .. 0.70`); fermion count reaches `N_need` `0.22` (`0.11 .. 0.45`); torsion contact equals pair energy `0.042` (`0.026 .. 0.064`) | same | same |
| Count against the paper's `N_need = 1.8e78` fermions | `1.3e78 = 0.73 N_need` (`0.41 .. 1.07`) | thread length `2M` | `child-energy-floor-and-thread-length-note.md` |
| No-production baseline | `n_b = 1.20e105 m^-3`, `eps_b = 20.2 eps_P`, `3.99 E_P` per fermion; child/parent energy `5.18e19`; floor `7.9e40` solar masses at `N_need` | degenerate gas with torsion contact | same |
| Thread length `L` (units of `2GM/c^3`) | order one: `0.3 .. 2` by three routes (infall span, thermal entropy cap, entanglement cap) | | `nesting-thread-length-note.md` |
| After the bounce | turnaround where the child's size equals its own Schwarzschild radius; mass `1.00 N m_n` (`1.5e21` solar masses) for matter, `1.4e14` for radiation; matter equals radiation at `eta = 9.0e-8`; no inflation | adiabatic, fixed `N` | `child-expansion-and-matter-fraction-note.md` (v2) |
| Parent's books | pair production leaves the parent's baryon number unchanged; inherited net matter `eta <= M*/M = 8.9e-21` (`M* = 1.8e11 kg`) | baryon number conserved | `matter-antimatter-and-time-note.md` |
| Time | entanglement entropy `4.11` nats per pair, `0.26 S_BH` at the stopping scale; arrow of time needs no imbalance; child lifetime `86 yr` (radiation) to `1.5e9 yr` (all matter) | one partner traced out | same |
| Level splitting | classical Hehl-Datta: `delta/T = 37.8 p`, `eta = 34.5 p` (`p = 1.7e-11` for `eta = 6e-10`). Second quantised: `A.A` splits by exactly 0; the vector-axial term splits by `4 kappa' (alpha/gamma) h`; `beta ~ 3.9e9` with inherited `h` | needs `alpha != 0` and a net helicity density | `torsion-splitting-and-inherited-spin-note.md`, `second-quantised-splitting-note.md` |
| Nesting | `N_f/N_need` repeats each generation; the scale ratio does not (all matter `1.1e20` then `1.2e40`); `Omega = 3.97` not reproduced. Observation with three free parameters: a `63` solar-mass parent (`52 .. 84`) gives a Hubble-mass child | | `nesting-thread-length-note.md` |
| Neck, topology | a ball child holds at most `M_p`; a large closed child with a small neck (open) needs a layer of negative energy `|E| >= sqrt(2R|M_p - m_-|)` | standard thin-shell junction, effective-fluid torsion | `israel-matching-note.md` |
| Neck, paths | mass-matched surface is timelike only from `t_exit/27` (`R = r_g/9`); exit at `t = G M_p/c^3`; thread path (slowest) ends at the throat after `pi G M_p/c^3` | | `israel-matching-note.md`, `junction-neck-note.md` (v2) |
| Swaps of the throat that keep the Dirac action | `T` on II<->IV: same-sign `n_V`, `n_5`; `CP` on I<->III and `CRT` on both pairs: opposite sign; child inherits `eta = 6.5e-21` under the `T` image | free Dirac field on Kruskal form | `swap-symmetry-note.md` |
| Quanta swept up by the neck | `0.65 .. 1.8 N0` quanta and `1.3 .. 3.2 M_p` of energy by `R = r_g`; net baryon number `~1e-9 N0` for `eta = 6e-10` | exterior vacuum, one Dirac field | `neck-flux-note.md` |
| Neck layer | torsion matter alone: pressure/energy `>= 1`, wrong sign for a thread layer; vacuum layer alone: ratio `-1`, runs the neck out at near light speed (child time `r_g/2`); balanced composite (half vacuum, half dust-like) or a layer with ratio `-1/2` fits | normal stress balanced by a vacuum-like fill | `torsion-layer-note.md` (v2), `casimir-layer-note.md`, `balanced-layers-note.md` |

## 2. Superseded or corrected (my own errors)

1. Shear pair density `C` was a factor 2 low (v1); corrected in v2.
2. The shear count was in pairs, compared with the paper's count of fermions; fermions are `2 x` pairs (`0.37` became `0.73` of `N_need`).
3. The generation map's "unstable fixed point at `1.8e11 kg`" was wrong; the child outweighs its parent at every mass above the Planck scale.
4. The neck at the bounce, `3.57e-23 m` independent of `N`, used the bare density; the geometric mass includes the negative torsion piece.
5. The lifetime ratios of the torsion layer (`3e33 .. 3.7e40`) held for the thread path only.
6. The composite layer's 125-field floor and `1e7 .. 1e41` cancellation came from leaving the normal stress unbalanced.

7. **10 October reset:** the Israel note conflated the high-density cancellation
   scale `eps=alpha*n^2` with the exact closed-FLRW bounce. The full equation
   has `adot=0` and `1+adot^2=1` at a real turning point, not zero. The claimed
   zero geometric mass there is withdrawn. Tiny near-cancellation neck-onset
   values require a consistent curvature treatment; the saved JSON is unchanged.
   See [the consolidated junction correction](#closed-child-junction).

## 3. Named missing ingredients (sized where possible)

- **Source of a net helicity density `h`** (a P-odd background) in the child; inherited helicity is number-bounded (`chi_h <= 8.9e-21`). A rotating parent supplies a pseudoscalar: `p_max = 1.26 a*/L` (`torsion-splitting-and-inherited-spin-note.md`).
- **What the negative-energy layer and its vacuum-like fill are made of.** Needed surface energy `s = sqrt(R r_g)/(4 pi R^2)`; the torsion part must sit within `7.6e-36` of `0.5443 n_b` at `R = 1 m` to be Planck-thick.
- **The neck's Planck-scale start** (`R < 3e-22 m`), the self-gravity of Planck-thick layers (`G sigma d ~ 1`), and energy exchange between the two layers.
- **Spin delta-function terms** at the neck (torsion is taken as an effective fluid), the neck's start time relative to the bounce, and the child's recollapse and a second neck.
- **What fixes the parent-to-child mass ratio** beyond topology; `Omega = 3.97` stays unexplained.
- **Single Kasner epoch** (no BKL oscillation) and the species count (22.5 Dirac fields is from memory).

## 4. Interface for the matter-antimatter question

Radiation relation `rho c^2 = 3 c^2/(32 pi G t^2)` at the child's time of neck opening `t = G M_p/c^3` (child clock; `49 .. 164 microseconds` over the neck paths, 10 solar masses):

| parent | `t_open` | `T` (one Dirac field, `g = 3.5`) | `T` (`g = 106.75`) |
|---|---|---|---|
| 3 solar masses | `15 microseconds` | `296 MeV` | `126 MeV` |
| 10 solar masses | `49 microseconds` | `162 MeV` (`89 .. 162` over paths) | `69 MeV` |
| 30 solar masses | `0.15 ms` | `94 MeV` | `40 MeV` |
| `1e6` solar masses | `4.9 s` | `0.51 MeV` | `0.22 MeV` |

`T_open` scales as `M_p^(-1/2)`; with `g = 106.75`, `T_open = 100 GeV` needs `M_p = 4.8e-6` solar masses and `T_open = 1 MeV` needs `4.8e4`. Energy inside the neck sphere at opening is `M_p c^2` for the mass-matched neck; at `T = 162 MeV` (mean energy `3.15 T`) that is `2.2e58` quanta, `1.8` times the parent's baryon number `1.2e58`. These are sizings, not results: the child's field content at opening is not known.

## 5. Open forks, in recommended order

1. Matter-antimatter with the neck: first pass done (`neck-flux-note.md`). By the time the neck reaches `r_g` it has swept `0.65 .. 1.8` times the parent's baryon number of quanta and `1.3 .. 3.2 M_p` of energy out of the child into the layer; net baryon number `eta x` that, so `~1e-9 N0` for `eta = 6e-10` and order `N0` only for `eta ~ 1`. Remaining: what the layer does with them, and the child's field content at the neck.
2. Net helicity density: weak interactions out of equilibrium in the child, the parent's `J` through the neck, or inheritance.
3. One outside kernel `Pi(k)` for both the late push (`~5.8e-10 J/m^3`) and the catch (`~3.9e74 J/m^3`).
4. Whether the two sides of a throat carry opposite `h` (swap symmetry): first pass done (`swap-symmetry-note.md`). Three swaps keep the Dirac equation, its mass and the `V.A` term: `T` exchanging the black-hole and white-hole interiors (same sign of matter and helicity), `CP` exchanging the two exteriors and `CRT` exchanging both pairs (opposite sign of both). Remaining: which swap the NSC throat realises; with the `T` image the child inherits `eta = N0/N_need = 6.5e-21`, `9.2e10` below `6e-10`.
5. The items in section 3 for the neck layer.

## 6. Consolidated derivations and assumptions

### Collapse and production

The cold, free-neutron-gas comparison gave a TOV maximum of `0.710` solar
masses at central density `4.0e18 kg/m^3`, radius `9.26 km`, versus the
Newtonian limiting mass `5.73` solar masses. These were session calculations,
not an independently bound repository result. Pressure contributes to GR
gravity. With the illustrative `H0 = 70 km/s/Mpc`, `Omega_Lambda = 0.7`, the
Sun's pull at 1 AU exceeds `Omega_Lambda H0^2 r` by `1.1e22`; that comparison
does not derive an NSC cosmological constant.

The Popławski closed-FLRW production check assumed `N ~ a^-3` and substituted
the paper's Eq. (4). It found `H_max = 7.77e23 /s`,
`H_max / omega_n = 0.544`, and the mode-count bound `(2/3) adot_max^3 = 2.8e39`
pairs, `38.8` orders below the required population. Its density
`7.7e84 /m^3` requires Fermi energy density `1.1e88 J/m^3`, against the
assumed `3.9e74 J/m^3`. The production law and bounce state therefore remain
unsupported within that model; shear production is a separate calculation.
The paper's non-reduced Planck-mass convention is required for Eqs. 5–6;
the Eq. 8 prefactor was not reproduced.

For one massless Dirac field in one Kasner epoch, the corrected pair density is
`C = n(c tau)^3 = Lambda_2 u_c^2/(2 pi^2)`, `Lambda_2 = 0.01424`.
At `u_c = 1`, `C = 7.2e-4`. Fermions are twice the pair count: the recorded
`0.37 N_need` comparison becomes `0.73 N_need`. The stop estimates use
`x_br^2 = 24 pi e` and `x_t^2 = pi C^2/e`, with `x = c tau/l_P`,
`eps = e/tau^4`, `n_f = 2C/tau^3`. Their ranges in section 1 depend on the
sharp `u_c` cut, a chosen thread length and omission of BKL epochs and torsion
from the field evolution. They do not supply a controlled Planck-scale history.

### Particle number, expansion and nesting

The free degenerate-gas contact model uses
`eps_eff = eps_F - alpha n^2`, `alpha = kappa (hbar c)^2/32`. At the
density-cancellation scale its no-production energy is `3.99 E_P` per
fermion, so the quoted floor at `N_need = 1.8e78` is `7.9e40` solar masses.
That is a matter-energy sizing under this equation of state, not the
geometric mass of an exact closed-FLRW bounce. The recorded paper mass
`1.7e31` solar masses was recalled from its Eq. 9 scaling without rederiving
the prefactor. Reading it at a later expansion epoch is an open interpretation.

With fixed `N`, adiabatic expansion and complete pair annihilation, the
turnaround equation is `a^2 - r_m a - r_rad^2 = 0`,
`r_m = 2G eta N m_n/c^2`, `r_rad^2 = (2G/c^4) N 3.99 E_P a_b`.
It gives the radiation/matter crossing `eta = 9.0e-8` and the masses in
section 1; surviving matter is assumed to be neutrons. Acceleration lasts
only to `a/a_b < sqrt(2)`, so this construction supplies no inflation.
The full generation map crosses child = parent only at `0.65` Planck masses,
where `N_f < 1`; the former fixed point `1.8e11 kg` was a false extension of
the matter-dominated formula. The thermal and entanglement caps are bounds
on a chosen state. The closed-radiation estimate
`L_c = (pi/3) [x(2-x)]^(3/2)/x` uses a free Kasner-onset parameter `x`.
None selects the constant recursive scale `Omega = 3.97`.

### Conserved numbers and operator splitting

Pair production leaves the parent's net baryon number unchanged. With
baryon conservation, inheritance alone gives `eta <= M*/M`,
`M* = m_P^2/(4 c1 m_n) = 1.8e11 kg`, `c1 = 0.3975`.
The `4.11` nats per pair assumes one partner is traced out; the whole pure
thread has zero entropy. Its arrow-of-time example already works at `eta=0`.

The classical, thermal Hehl–Datta estimate is `delta/T = 37.8 p`,
`eta = 34.5 p` for small polarisation; rotation could supply
`p_max = a*/(2 c1 L) = 1.26 a*/L`. The operator mean-field calculation
supersedes this as a particle/antiparticle splitting mechanism: minimal
`A.A` gives exactly zero splitting between equal-helicity partners.
The nonminimal FMT interaction gives `4 kappa' (alpha/gamma) h` from `A.V`
and `-4 kappa' alpha^2 n_V` from `V.V`,
`kappa' = (3 kappa/16) gamma^2/(gamma^2+1)`.
Here `alpha` is the nonminimal coupling, distinct from the contact coefficient
above. With `beta = 2 gamma alpha/(gamma^2+1)`,
`eta = 34.5 (beta/2) chi_h`; inherited `chi_h <= 8.9e-21` for the example
requires `beta ~ 3.9e9` to reach `eta = 6e-10`. Thermal equilibrium, a held
helicity background and one massless species are assumptions. No source of
that helicity or asymmetry was derived.

The free-Dirac Kruskal swap test preserves both mass and the `V.A` term for
`T` on II↔IV, `CP` on I↔III and `CRT` on both pairs. `T` preserves the signs
of `n_V,n_5`; `CP,CRT` reverse both. The inherited `T` image gives
`eta = N0/N_need = 6.53e-21`; the calculation does not choose the physical
swap and does not make the whole child's energy an inherited parent image.

<a id="closed-child-junction"></a>
### Closed-child junction and neck stress

In geometric units, `R = a sin(chi)` and
`2m = a sin(chi)^3(1 + adot^2)`. At a real turning point `adot=0`, so
`m = a sin(chi)^3/2`. Exact density cancellation gives `adot^2=-1` in the
closed equation and cannot establish zero mass at a real bounce. The tiny
near-cancellation onset numbers remain unverified approximation outputs.

For a timelike shell, `beta^2 = 1-2m/R+Rdot^2`,
`sigma = -[beta]/(4 pi R)` and
`P = ([K^tau_tau]+[K^theta_theta])/(8 pi)`, with the normal directed
from child to exterior. The radiation-era mass match uses
`R_j^3 = r_gp a^2/(1+adot^2)`; its timelike segment starts at
`t_exit/27`, `R=r_g/9`. A shell-free vacuum match requires a comoving
boundary and zero child pressure. A large closed child with
`beta_-<0,beta_+>0` needs negative surface energy, with
`|E| >= sqrt(2R|M_p-m_-|)`. These thin-shell statements treat torsion as an
effective fluid and omit spin delta functions, recollapse and a physical
layer model. Unread FMM papers do not independently support the example.

Codazzi gives `sigma(q+2w) Rdot/R = -T_un`,
`q = d ln|E|/d ln R`, `T_un = (eps+p) a chi' gamma^2`.
For the outgoing, child-rest-frame segment with `sigma<0,T_un<=0`,
`w <= -q/2`. Torsion with `eps_eff<0` has `w>=1`, incompatible with a
rising minimum `|E|`; its lifetime estimates apply only to the thread path.
The bag-boundary Dirac slab has
`E/A = -(7 pi^2/2880) hbar c/d^3`, `w=-1`. Normal-stress balance turns
its filled surface into a vacuum-like layer and the filled torsion constituent
into dust-like stress. Half of each gives the thread's `w=-1/2`, with the
`7.6e-36` density tuning at `R=1 m` in section 3. The fill, energy exchange,
Planck-scale onset and self-gravity remain missing.

With exterior vacuum, the sweep ledger is
`dE = A (4/3) eps v gamma dt`, `dN = A n_f v dt`, `dB = eta dN`,
`v=|a chi'|`. Quanta enter the layer, not the exterior. The `w=-1/2`
trajectory's end energy varies from `-6.1` to `-6.8 M_p` with integration
tolerance; the `1.3–3.2 M_p` sweep range excludes the `w=-1` run's much
larger energy. Layer transport and field content at the neck remain open.
