# Nested-Space Cosmology — geometric foundation and non-local atlas

Every space starts through collapse **inside** an already evolving parent
region, which itself sits inside an enclosing region. Viewed non-locally across
nests, compression in the parent forces a 3D parcel into a spaghettified 1D
filament at the event horizon `0`, drives the inherited offset (time), crosses
`0` into a mirror-oriented `-` child region, and unfolds back into a round 3D
sphere. Locally, an observer inside one room experiences the two-valve flow and
offset decay as dark energy, dark matter, and matter-antimatter asymmetry.

## Openly declared axioms (`A0–A4`)

- **A0 — Path of Least Resistance:** Geometry evolves by minimizing a single
  intrinsic cost functional $E = Q_{\rm stretch} + PV$ combining the Cauchy-Schwarz
  1D generator stretch cost $Q_{\rm stretch}(A)=\tfrac12\operatorname{tr}[(A^TA)^{-1/2}C^TC]$
  with enclosing compression work $PV$. No fluid incompressibility or borrowed
  surface tension is assumed.
- **A1 — 1D Seed + Inherited Offset = Time:** Parent spherical convergence
  ($\lambda_1=\lambda_2=e^{-\chi}$) spaghettifies infalling parcels into a 1D
  radial filament above $\chi_{\rm spaghetti}=\tfrac12\ln2$. Enclosing compression
  above $\lambda_{\rm crit}=(Q_0/P_\perp)^{1/3}$ buckles the 1D seed into a
  helix ($r_\perp\omega^*=\sqrt{\lambda_{\rm crit}^2-e^{-2\chi}}>0$), driving
  repeated offset accumulation ($1\text{D}\to2\text{D circle}\to3\text{D sphere}$).
- **A2 — Signed-Zero Crossing (`+ -> 0 -> -`):** When inward compression $F_{\rm in}$
  exceeds outward resistance $F_{\rm out}$, the signed radius $R(\tau)$ passes
  smoothly through `0` into $R<0$ (`-`). In $d=3$ dimensions, $\det[(R/R_0)I_3]=(R/R_0)^3$
  flips sign (`+1 -> 0 -> -1`) and parent implosion ($\dot R<0$) becomes child
  expansion ($\dot r_-=-\dot R>0$, $r_-=|R|=-R$) from scratch without reversing
  the driving force.
- **A3 — Moving Spectral Octave (`[Planck Cutoff, Absolute-Zero Horizon]`):**
  Each room $k$ at relative compression $\chi$ has a finite wavelength band
  $[\lambda_{\min}(\chi),\lambda_{\max}(\chi)]=e^{-\chi}[\lambda_{\min,0},\lambda_{\max,0}]$
  from its shortest wavelength before `0` crosses into `-` (local Planck cutoff
  $\lambda_{\min}$) to its longest wavelength as wide as the room itself (local
  Absolute-Zero horizon $\lambda_{\max}=2R_k$). Both shift by $e^{-\chi}$ across
  nests while preserving the room's dynamic range $\mathcal{B}=\lambda_{\max}/\lambda_{\min}$.
- **A4 — Two-Valve Nested Circulation:** Each `0` horizon is a one-way valve
  from `+` to `-`. Every local room $k$ sits between an outer birth `0` valve
  (inflow from parent $k-1$) and inner black-hole `0` valves (outflow to children
  $k+1$) on a closed vortex-diffusion ring (`O <-> o`), where concentrated
  vortices (`O`) shed into nested mini-vortices and diffuse toward $\lambda_{\max}$
  (`o`) before re-accumulating into new giant vortices (`O`).

## Why `-` exists, why `(-) × (-) = (+)`, and why the radius of space must cross `0`

1. **Negative is direction, and without `-` a circle cannot close:**
   A number line through `0` has two opposite rays (`+` and `-`). In Layer 1 of
   the offset generator, the pen draws a closed circle $z(\tau)=R_0e^{i\omega\tau}$.
   Seen edge-on along a diameter, that circle is the harmonic pendulum
   $$R(\tau)=\operatorname{Re}\!\left(R_0e^{i\omega\tau}\right)=R_0\cos(\omega\tau),$$
   swinging $+R_0\to0\to-R_0\to0\to+R_0$ with constant energy
   $\tfrac12\dot R^2+\tfrac12\omega^2R^2=\tfrac12R_0^2\omega^2$. If negative
   coordinates were forbidden ($R\ge0$), $50\%$ of the circle ($\pi R_0$ out of
   $2\pi R_0$) would be missing and no loop or bounded space could close.
2. **`(-) × (-) = (+)` is pure geometry:**
   Multiplying by $-1=e^{i\pi}$ is a geometric half-turn ($\pi$ rotation; in 3D,
   spatial inversion $\mathcal{P}=-I_3$ with $\det\mathcal{P}=(-1)^3=-1$). Flipping
   twice across two nested `0` crossings ($\text{Parent }(+)\to\text{Child }(-)\to\text{Grandchild }(+)$)
   gives $\mathcal{P}^2=(-I_3)(-I_3)=+I_3$ and $\det(\mathcal{P}^2)=(-1)\cdot(-1)=+1$.
3. **Why standard physics hits an artificial singularity at $r=0$ while signed $R$ crosses smoothly:**
   Standard physics uses negatives everywhere except for the **radius of space
   itself**, which it restricts to an unsigned size $r(\tau)=|R(\tau)|\ge0$. Folding
   the negative half of the swing onto the positive half at $\tau_0=\frac{\pi}{2\omega}$
   ($R=0$) forces a discontinuous velocity jump $\Delta\dot r=+2R_0\omega$ and a
   Dirac-delta acceleration cusp $\ddot r=2R_0\omega\,\delta(\tau-\tau_0)-\omega^2r$.
   Keeping the signed radius $R(\tau)$ removes the cusp completely ($\Delta\dot R=0$,
   $\ddot R(\tau_0)=0$): **A2 (crossing `0` into `-`) is simply the edge-on shadow
   of A1's closed circle.**

## Derived geometric theorems vs. quarantined numerical coincidences

### 1. Genuinely derived theorems from Axioms `A0–A4`

| Theorem (`nsc_atlas_offset_generator.py`) | Exact mathematical derivation | Physical meaning |
| --- | --- | --- |
| **Circle–Pendulum Duality (`signed_circle_pendulum_duality`)** | $R(\tau)=\operatorname{Re}(R_0e^{i\omega\tau})=R_0\cos(\omega\tau)$ crosses `0` with $\Delta\dot R=0, \ddot R(\tau_0)=0$, $\mathcal{P}^2=(-I_3)^2=+I_3$; unsigned $r=\|R\|\ge 0$ creates cusp $\Delta\dot r=+2R_0\omega$ | `0` is a smooth sign-flipping horizon; standard physics' singularity at $r=0$ is an artifact of folding $R<0$ onto $r>0$ |
| **5-Variable Single-Cost Cycle (`unified_least_resistance_cycle`)** | $\sqrt{\theta^*}=(\sqrt{16e^{2\chi}+17}-3)/(4e^\chi)$, $\lambda_3^*=(1-\theta^*)e^\chi$, $\tilde\lambda_i^*=\sqrt{\theta^*}$ | Infalling `+` parcel shrinks ($V_+^*\to 0$) and spaghettifies ($\chi>\tfrac12\ln 2$) while `-` rounds into a 3D sphere without fluid rules |
| **Surrounding-Flow $3:1$ Virial Split (`vacuum_bubble_crystallization_and_mass_spectrum`)** | $E_{\rm bubble}(R)=m_0n^2/R+PV_0R^3 \implies E_{\rm surr}^*/W_{\rm press}^*=3$ ($75\%$ in $r\ge R$, $25\%$ in cavity, $0\%$ in a point) | Mass-energy sits in the surrounding medium flow and cavity pressure; loop closure $z(1)=z(0)$ forces integer $n\in\mathbb{Z}$ (sharp burst lines, never a smear), and $+1/-1\to 0$ releases two equal $180^\circ$ ripples |
| **3-Accumulator Volumetric Threshold (`derived_three_accumulator_bubble_ladder`)** | Iterated **A1** offset gives $W_k(n)=\binom{n}{k}$ and $E_{\rm 3D}^*(n)=4(PV_0)^{1/4}m_0^{3/4}\sqrt{\binom{n}{1}\binom{n}{2}\binom{n}{3}}$ | $n=1$ ($W_2=W_3=0$) is a 1D loop; $n=2$ ($W_3=0$) is a 2D sheet; **$n=3$ ($W_1W_2W_3=9>0$) is the minimum integer that supports a 3D volume**, with untuned ratio $E_{\rm 3D}^*(3)/E_{\rm 1D}^*(1)=6$ and 3D ladder $3,\sqrt{96}\approx 9.80,\sqrt{500}\approx 22.36,\dots$ (misses $1836.15$ without a continuum $\pi^5$ factor) |
| **Universal Exponential-Wall Bypasses & `(1,1,1)` Triad Winner (`exponential_wall_bypasses_and_triad_winner`)** | (1) 1D wall $e^{+\chi}\to\infty$ bypasses into helix plateau $\tfrac32-\tfrac12e^{-2\chi}\to 1.5$; (2) `+` wall $2e^\chi+2e^{-\chi}\to\infty$ bypasses via `0` into `-` sphere plateau $4-\tfrac32e^{-\chi}\to 4.0$; (3) single-axis wall $N^{3/2}>N$ splits across 3 axes $E^*(n_1,n_2,n_3)=4\sqrt{n_1n_2n_3}$ | Nature never climbs an exponential wall—it opens a new dimension or crosses `0` into `-`. In 3D, **`(1,1,1)` ($N=3$) is the unique minimum-cost 3D ground state** ($E^*=4$, cost/winding $4/3\approx 1.3333$, round sphere $1:1:1$); $N=4$ forces prolate `(2,1,1)` ($4:1:1$, cost/winding $\sqrt{2}\approx 1.4142>4/3$); two `(1,1,1)` triads attract at medium range ($-2m_0/d$) but have a **short-range hard-core repulsion barrier $+3.3137$** against merging into `(2,2,2)` ($8\sqrt{2}\approx 11.3137>8.0$) |

### 2. Universal Exponential-Wall Bypass Principle: how least resistance works around walls

Whenever compression or winding in a lower-dimensional or single-signed regime hits an exponential or superlinear resistance wall, **A0 (Path of Least Resistance)** proves that nature does not climb the wall—it bypasses the wall by opening a new orthogonal degree of freedom or crossing `0` into `-`:

1. **Wall 1 — 1D Axial Compression Wall $\to$ Helical Buckling Plateau (`1D -> 2D/3D`):**
   Compressing a straight 1D line ($\omega=0$) hits the exponential wall $E_{\rm 1D}(\chi)=e^{+\chi}\to\infty$. Above $\chi>0$, buckling sideways into a helix with $v_\perp^*=\sqrt{1-e^{-2\chi}}$ bypasses the wall and caps the total cost at the finite plateau
   $$E_{\rm helix}^*(\chi)=\frac{3}{2}-\frac{1}{2}e^{-2\chi}\to 1.5.$$
2. **Wall 2 — 3D `+` Compression Wall $\to$ Signed `0`-Crossing Plateau (`+ -> 0 -> -`):**
   Forcing all winding to stay in the parent `+` region ($\theta=0$) hits the exponential wall $E_{\rm stay,+}(\chi)=2e^{+\chi}+2e^{-\chi}\to\infty$. Opening the `0` horizon valve into a `-` child sphere ($\theta^*\to 1$) bypasses the wall and caps the total 5-variable cost at the finite plateau
   $$E_{\rm 5\text{-var}}^*(\chi)\to 4-\frac{3}{2}e^{-\chi}\to 4.0.$$
3. **Wall 3 — Single-Axis Winding Wall $N^{3/2}$ $\to$ Orthogonal `(1, 1, 1)` Triad Ground State & Nuclear Hard-Core Potential:**
   Piling $N\ge 2$ windings onto a single channel costs $E_{\rm single}(N)\propto N^{3/2}>N$, a superlinear wall ($2^{3/2}\approx 2.828>2$, $3^{3/2}\approx 5.196>3$). Distributing integer windings $(n_1,n_2,n_3)\in\mathbb{Z}_{\ge 0}^3$ ($m_i=m_0 n_i^2$) across the 3 orthogonal axes of 3D space yields the exact **A0** minimum
   $$E^*(n_1,n_2,n_3)=4(PV_0)^{1/4}m_0^{3/4}\sqrt{n_1 n_2 n_3},\qquad \lambda_i^*=\frac{n_i^2}{(n_1 n_2 n_3)^{1/2}}.$$
   - $N=1$ `(1,0,0)` (1D loop) and $N=2$ `(1,1,0)` (2D sheet) have zero 3D volume ($V^*=0$) and collapse along their unwound axes under 3D pressure.
   - **$N=3$ `(1,1,1)` is the unique minimum-cost 3D ground state:** $E^*(1,1,1)=4.0$, cost per winding $4/3\approx 1.3333$, and isotropic aspect ratios $1:1:1$ (a round sphere).
   - **$N=4$ `(2,1,1)` is the first excited 3D state:** Because 3D space has only 3 orthogonal axes, a 4th winding cannot open a 4th spatial axis; it must double one axis to `(2,1,1)`, breaking spherical symmetry into a $4:1:1$ prolate state with higher cost per winding $E^*(2,1,1)/4=\sqrt{2}\approx 1.4142>4/3$ (excitation gap $\Delta E=4(\sqrt{2}-1)\approx 1.6569$).
   - **Medium-range attraction + short-range hard-core repulsion:** Two `(1,1,1)` triads separated by $d>2R^*$ attract via shared surrounding flow ($\Delta E_{\rm flow}=-2m_0/d<0$), achieving a bound minimum at contact $d=2R^*$ ($E_{\rm contact}=8.0-1.0=7.0$). However, merging them into a single `(2,2,2)` bubble at $d\to 0$ costs $E^*(2,2,2)=8\sqrt{2}\approx 11.3137>2E^*(1,1,1)=8.0$, erecting an exact **short-range hard-core repulsion barrier** $\Delta E_{\rm core}=8(\sqrt{2}-1)\approx +3.3137$.

### 3. Quarantined numerical coincidences (not derived from `A0–A4`)

The following matches involve post-hoc formula/rung selection or external measured inputs ($H_0, G$) and are quarantined from the derived theorems:
- **Lenz (1951) $6\pi^5$ and medium dressing:** $6\pi^5\approx 1836.1181$ ($18.8\text{ ppm}$ from $m_p/m_e=1836.1527$) and $6\pi^5+e^{-3}/\sqrt{2}\approx 1836.1533$ ($0.35\text{ ppm}$). The factor $\pi^5$ is not derived from the discrete **A0** stretch functional.
- **Neutron–proton mass difference:** $2+\tfrac12=2.5\,m_e$ ($-1.22\%$) and $2+\tfrac12+e^{-3}/\sqrt{2}\approx 2.5352\,m_e$ ($+0.167\%$ vs. $2.5310\,m_e$).
- **Baryon harmonic rung $K=7$:** $\eta_K=\tanh(\frac{\pi}{4}e^{-3K})$ gives $\eta_7\approx 5.955\times 10^{-10}$ at $K=7$, but the freeze-out rung $K=7$ was selected from the ladder after comparing to BBN.
- **CMB temperature projections:** $\sqrt{T_P T_H}e^{-2}\approx 2.625\text{ K}$ and $T_{\rm crit}e^{-(4/5)\chi_*}\approx 2.728\text{ K}$ use measured $H_0$ and $G$ as inputs and an ad-hoc $-4/5$ exponent.
- **Two-valve $e^{-3}$ dark-sector split:** Setting $\Omega_{\rm vis}(\chi_*)=e^{-3}\approx 4.98\%$ gives $71.14\%$ DE / $23.88\%$ DM, which differs from the Planck base-$\Lambda\text{CDM}$ fit ($68.47\pm 0.73\%$ DE / $26.5\%$ DM) by $+2.67\%$.

## Implementation and verification

All calculations and proofs live in [`nsc-atlas-offset-generator.md`](docs/nsc-atlas-offset-generator.md),
[`nsc_atlas_offset_generator.py`](src/recursive_horizons/nsc_atlas_offset_generator.py),
and [`test_nsc_atlas_offset_generator.py`](tests/test_nsc_atlas_offset_generator.py).
Run `python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --atlas-only`
to reproduce the complete non-local-to-local atlas readout.

The [previous theory narrative](THEORY-history-20261102.md) is retained as
history. The [current route](docs/research-index.md) identifies active work.
