# Nested-Space Cosmology — geometric foundation and non-local atlas

**Foundation Edition (`v0.27.0`)** · [Zenodo Record (`10.5281/zenodo.23890258`)](https://zenodo.org/records/23890258) · [Main Manuscript](paper/nested-space-cosmology.md) · [Illustrated PDF](paper/nested-space-cosmology.pdf)

Every space starts through collapse **inside** an already evolving parent
region, which itself sits inside an enclosing region. Viewed non-locally across
nests, compression in the parent forces a 3D parcel into a spaghettified 1D
filament at the event horizon `0`, drives the inherited offset (time), crosses
`0` into a mirror-oriented `-` child region, and unfolds back into a round 3D
sphere. Locally, an observer inside one room experiences the two-valve flow,
finite spectral resolution band, and offset decay as dark energy, dark matter,
baryonic triads, stellar black-hole collapse thresholds, and matter-antimatter
asymmetry.

## Openly declared axioms (`A0–A4`) as the fundamental generating law

In Nested-Space Cosmology, **Axioms A0–A4 are the primary geometric law that generates space, time, and boundary exchange**, whereas 3+1D continuum field equations and general relativity are coarse-grained local dictionary readouts (Propositions 1.1–1.2b) of what an interior observer measures before the single-horizon approximation breaks down (Part I, `NSC-P006`).

- **A0 — Path of Least Resistance:** Geometry evolves by minimizing a single
  intrinsic cost functional $E = Q_{\mathrm{stretch}} + PV$ combining the Cauchy-Schwarz
  1D generator stretch cost $Q_{\mathrm{stretch}}(A)=\tfrac{1}{2}\mathrm{tr}[(A^TA)^{-1/2}C^TC]$
  with enclosing compression work $PV$. No fluid incompressibility or borrowed
  surface tension is assumed.
- **A1 — 1D Seed + Inherited Offset = Time:** Parent spherical convergence
  ($\lambda_1=\lambda_2=e^{-\chi}$) spaghettifies infalling parcels into a 1D
  radial filament above $\chi_{\mathrm{spaghetti}}=\tfrac{1}{2}\ln 2$. Enclosing compression
  above $\lambda_{\mathrm{crit}}=(Q_0/P_\perp)^{1/3}$ buckles the 1D seed into a
  helix ($r_\perp\omega^\star=\sqrt{\lambda_{\mathrm{crit}}^2-e^{-2\chi}}>0$), driving
  repeated offset accumulation ($1\text{D}\to 2\text{D circle}\to 3\text{D sphere}$)
  and generating the effective 3+1D metric tensor $g_{\mu\nu}$ (Proposition 1.1b).
- **A2 — Signed-Zero Crossing (`+ -> 0 -> -`):** When inward compression $F_{\mathrm{in}}$
  exceeds outward resistance $F_{\mathrm{out}}$, the signed radius $R(\tau)$ passes
  smoothly through `0` into $R<0$ (`-`). In $d=3$ dimensions, $\det[(R/R_0)I_3]=(R/R_0)^3$
  flips sign (`+1 -> 0 -> -1`) and parent implosion ($\dot R<0$) becomes child
  expansion ($\dot r_-=-\dot R>0$, $r_-=|R|=-R$) from scratch without reversing
  the driving force.
- **A3 — Moving Spectral Octave (`[Planck Cutoff, Absolute-Zero Horizon]`):**
  Each room $k$ at relative compression $\chi$ has a finite wavelength band
  $[\lambda_{\min}(\chi),\lambda_{\max}(\chi)]=e^{-\chi}[\lambda_{\min,0},\lambda_{\max,0}]$
  from its shortest wavelength before `0` crosses into `-` (local Planck ceiling
  $\lambda_{\min}=\ell_P$) to its longest wavelength as wide as the room itself (local
  Absolute-Zero floor $\lambda_{\max}=2R_k$). Both shift by $e^{-\chi}$ across
  nests while preserving the room's dynamic range $\mathcal{B}=\lambda_{\max}/\lambda_{\min}$.
- **A4 — Two-Valve Nested Circulation:** Each `0` horizon is a one-way valve
  from `+` to `-`. Every local room $k$ sits between an outer birth `0` valve
  (inflow from parent $k-1$) and inner black-hole `0` valves (outflow to children
  $k+1$) on a closed vortex-diffusion ring (`O <-> o`), where concentrated
  vortices (`O`) shed into nested mini-vortices and diffuse toward $\lambda_{\max}$
  (`o`) before re-accumulating into new giant vortices (`O`).

## Why `-` exists, why `(-) × (-) = (+)`, and why a sphere has finite facets

1. **Negative is direction, and without `-` a circle cannot close:**
   A number line through `0` has two opposite rays (`+` and `-`). In Layer 1 of
   the offset generator, the pen draws a closed circle $z(\tau)=R_0e^{i\omega\tau}$.
   Seen edge-on along a diameter, that circle is the harmonic pendulum
   $$R(\tau)=\mathrm{Re}\!\left(R_0e^{i\omega\tau}\right)=R_0\cos(\omega\tau),$$
   swinging $+R_0\to 0\to -R_0\to 0\to +R_0$ with constant energy
   $\tfrac{1}{2}\dot R^2+\tfrac{1}{2}\omega^2R^2=\tfrac{1}{2}R_0^2\omega^2$. If negative
   coordinates were forbidden ($R\ge 0$), $50$% of the circle ($\pi R_0$ out of
   $2\pi R_0$) would be missing and no loop or bounded space could close.
2. **`(-) × (-) = (+)` is pure geometry:**
   Multiplying by $-1=e^{i\pi}$ is a geometric half-turn ($\pi$ rotation; in 3D,
   spatial inversion $\mathcal{P}=-I_3$ with $\det\mathcal{P}=(-1)^3=-1$). Flipping
   twice across two nested `0` crossings ($\text{Parent }(+)\to\text{Child }(-)\to\text{Grandchild }(+)$)
   gives $\mathcal{P}^2=(-I_3)(-I_3)=+I_3$ and $\det(\mathcal{P}^2)=(-1)\cdot(-1)=+1$.
3. **Why standard physics hits an artificial singularity at $r=0$ while signed $R$ crosses smoothly:**
   Standard physics uses negatives everywhere except for the **radius of space
   itself**, which it restricts to an unsigned size $r(\tau)=|R(\tau)|\ge 0$. Folding
   the negative half of the swing onto the positive half at $\tau_0=\frac{\pi}{2\omega}$
   ($R=0$) forces a discontinuous velocity jump $\Delta\dot r=+2R_0\omega$ and a
   Dirac-delta acceleration cusp $\ddot r=2R_0\omega\delta(\tau-\tau_0)-\omega^2r$.
   Keeping the signed radius $R(\tau)$ removes the cusp completely ($\Delta\dot R=0$,
   $\ddot R(\tau_0)=0$, finite Ricci scalar $R_{\mathrm{Ricci}}(0)=6\omega^2$ and Kretschmann scalar $K(0)=12\omega^4$): **A2 (crossing `0` into `-`) is simply the edge-on shadow of A1's closed circle.**
4. **Why a sphere does not have infinite sides (`Proposition 9.2b`):**
   Textbooks assume a sphere is an infinitely smooth continuum ($F\to\infty$). Under **Axiom A3**, every room is bounded between the Absolute-Zero floor $\lambda_{\max}=2R$ and the Planck ceiling $\ell_P=\sqrt{\hbar G/c^3}$. A polyhedral sphere with $F<\infty$ facets has excess isoperimetric ratio $A_F^3/(36\pi V_F^2)=1+\pi\sqrt{3}/F+O(F^{-2})>1$, decreasing monotonically toward $1$ as $F$ grows. A star distributes enclosing pressure across the room's frequency band until its surface facets hit the Planck ceiling ($\mathcal N_{1\mathrm{D}}=\lambda_n/\ell_P=M_P/m_n=1.29942\times 10^{19}$, $F_\star=(\lambda_n/\ell_P)^2=1.688\times 10^{38}$ facets, $N_\star=(\lambda_n/\ell_P)^3=2.194\times 10^{57}$ baryons, $\sqrt{R_\star\ell_P}=\lambda_n$). At that ceiling, the enclosing room forces spaghettification ($\chi_{\mathrm{spaghetti}}=\tfrac{1}{2}\ln 2$, neutron-star radius $R_{\mathrm{NS}}=4\lambda_n^2/\ell_P=10.916\ \mathrm{km}$, compactness $R_S/R_{\mathrm{NS}}=1/2=\theta^\star(\chi_{\mathrm{spaghetti}})$) and collapses the star inward across `0` into `-` above the Theorem 2(3) `(1,1,1)` hard-core barrier $M_{\mathrm{TOV}}\in[2.0414,2.1454]\ M_\odot$, inflating an internal child space. Because the Planck length $\ell_P$ is the smallest resolution step of the room ($\mathcal R(0)=R_{\mathrm{core}}=\ell_P>0$, Proposition 1.2), **a zero-volume point singularity does not exist**: the Planck ceiling is the exact pen tip where parent compression ends and the 1D pen (**Axiom A1**) starts drawing out the new child room (`-`), deflating in frequency (stretching in wavelength) from the Planck ceiling back down toward the new room's Absolute-Zero floor $\lambda_{\max,k+1}=2R_{k+1}$ (*setting the room again* with invariant octave capacity $\mathcal B_{k+1}=\mathcal B_k$).

## Derived geometric theorems vs. quarantined numerical coincidences

### 1. Genuinely derived theorems from Axioms `A0–A4`

| Theorem (`nsc_atlas_offset_generator.py`) | Exact mathematical derivation | Physical meaning |
| --- | --- | --- |
| **Circle–Pendulum Duality & Bounce (`signed_circle_pendulum_duality`, Prop. 1.2–1.2b)** | $R(\tau)=\mathrm{Re}(R_0e^{i\omega\tau})=R_0\cos(\omega\tau)$ crosses `0` with $\Delta\dot R=0, \ddot R(\tau_0)=0$, $R_{\mathrm{Ricci}}(0)=6\omega^2$, $K(0)=12\omega^4$; unsigned $r=\|R\|\ge 0$ creates cusp $\Delta\dot r=+2R_0\omega$ | `0` is a smooth sign-flipping horizon; standard physics' singularity at $r=0$ is an artifact of folding $R<0$ onto $r>0$ |
| **5-Variable Single-Cost Cycle & Open-Valve EoS (`unified_least_resistance_cycle`, Theorem 1)** | $\sqrt{\theta^\star}=(\sqrt{16e^{2\chi}+17}-3)/(4e^\chi)$, $9DVB=4(V+B/2)^3$, open-valve $dE_{\mathrm{DE}}^\star=(-\tfrac{1}{3}+\tfrac{4}{3})\rho_{\mathrm{DE}} dV_-=+\rho_{\mathrm{DE}} dV_-\Rightarrow w_{\mathrm{DE}}=-1$ | Infalling `+` parcel shrinks ($V_+^\star\to 0$) and spaghettifies ($\chi>\tfrac{1}{2}\ln 2$) while `-` rounds into a 3D sphere; open `0` valve maintains constant stretch density without a closed-box repulsive fluid |
| **Universal Exponential-Wall Bypasses & `(1,1,1)` Triad Winner (`exponential_wall_bypasses_and_triad_winner`, Theorem 2)** | (1) 1D wall $e^{+\chi}\to\infty$ bypasses into helix plateau $\tfrac{3}{2}-\tfrac{1}{2}e^{-2\chi}\to 1.5$; (2) `+` wall $2e^\chi+2e^{-\chi}\to\infty$ bypasses via `0` into `-` sphere plateau $4-\tfrac{9}{2}e^{-\chi}+O(e^{-2\chi})\to 4.0$; (3) single-axis wall $N^{3/2}>N$ splits across 3 axes $E^\star(n_1,n_2,n_3)=4\sqrt{n_1n_2n_3}$ | Nature never climbs an exponential wall—it opens a new dimension or crosses `0` into `-`. In 3D, **`(1,1,1)` ($N=3$) is the unique minimum-cost 3D ground state** ($E^\star=4$, cost/winding $4/3\approx 1.3333$, round sphere $1:1:1$); two `(1,1,1)` triads attract at medium range ($-2m_0/d$) and have a **short-range hard-core repulsion barrier $\Delta E_{\mathrm{core}}=8(\sqrt{2}-1)\approx +3.3137$** |
| **Finite-Facet Sphere & Stellar Collapse (`stellar_resolution_and_neutron_star_collapse`, Prop. 9.2b)** | $A_F^3/(36\pi V_F^2)=1+\pi\sqrt{3}/F+O(F^{-2})$, $\sqrt{R_\star\ell_P}=\lambda_n$, $M_{\star,0}=M_P^3/m_n^2=1.8481\ M_\odot$, $R_{\mathrm{NS}}=4\lambda_n^2/\ell_P=10.916\ \mathrm{km}$, $M_{\mathrm{TOV}}\in[2.0414,2.1454]\ M_\odot$ | A sphere has $F_\star=(\lambda_n/\ell_P)^2=1.688\times 10^{38}$ Planck facets; hitting the Planck ceiling forces spaghettification ($R_S/R_{\mathrm{NS}}=1/2$) and black-hole collapse into a child space |
| **Surrounding-Flow $3:1$ Virial Split (`vacuum_bubble_crystallization_and_mass_spectrum`, Theorem 3)** | $E_{\mathrm{bubble}}(R)=m_0n^2/R+PV_0R^3 \Rightarrow E_{\mathrm{surr}}^\star/W_{\mathrm{press}}^\star=3$ ($75$% in $r\ge R$, $25$% in cavity), $\Delta E_{\mathrm{bind}}^\star=-m_0/d$, $M_{\mathrm{eff}}=\tfrac{3}{2}M_{\mathrm{cav}}$ | Mass-energy sits in the surrounding medium flow and cavity pressure ($3:1$ exterior-to-cavity split by Gauss divergence); potential-flow added mass clusters with baryonic cores |
| **Surface-Area Conservation & Quadrupole Cubic (`Theorem 4 & Theorem 5`)** | $\sum_i 4\pi(R_{n_i}^\star)^2=4\pi(R_N^\star)^2$, $p^3-3p-3=0\Rightarrow p_+=\phi^{2/3}+\phi^{-2/3}=2.1038034027$ ($\phi=\tfrac{3+\sqrt{5}}{2}$) | Same-sign bubble fission conserves horizon area while releasing kinetic budget; exact causal exterior neck-instability growth rate |
| **Harmonicity Cascade & $\chi_{\mathrm{freeze}}=3$ Dark-Sector Closure (`Theorem 6, Theorem 7, Prop. 9.7b`)** | $\partial_{\zeta_x}v+\partial_{\zeta_y}v+\partial_{\zeta_z}v=0\Rightarrow \chi_{\mathrm{freeze}}=3\chi_1=3$, $\Omega_{\mathrm{DE}}=70.88\text{--}71.14$%, $\Omega_{\mathrm{DM}}=23.81\text{--}23.88$%, $\Omega_b=e^{-3}=4.98$%, $\Delta\Omega_{\mathrm{WHIM}}=0.33$%, $H_0=67.03$, $S_8=0.776$ | 3-axis harmonic cascade freezes at $\chi_{\mathrm{freeze}}=3$ ($\chi_{\mathrm{freeze,dyn}}=3.0754$), closing the cosmic budget, late WHIM excess, Hubble rate, and $S_8$ growth suppression with zero free parameters |
| **Thermodynamic Room Lifecycle & Optical Clock (`Proposition 9.7c`)** | $[\ell_P(\tau),2R(\tau)]=e^{-\chi(\tau)}[\ell_{P,0},2R_0]$, $ds=d\tau/\lambda(\tau)$, $\mathcal P^2=(-I_3)^2=+I_3$ | Two crossings close the cycle without an infinite tower or sterile branches; primary $1\to 2\to 4\to 8$ cascade nodes seed early SMBHs ($z>10$); BH `0`-valves ($E_{\mathrm{DM}}^\star=23.88$%) grip the horizon ledge to stabilize $\lambda_{\mathrm{Goldilocks}}$ so `(1,1,1)` matter, stars, and planets crystallize like ice in a pressure chamber at $\chi_{\mathrm{freeze}}=3$; optical clock accelerates during birth ($\ell_P$) and deflation (`0`) |

### 2. Universal Exponential-Wall Bypass Principle: how least resistance works around walls

Whenever compression or winding in a lower-dimensional or single-signed regime hits an exponential or superlinear resistance wall, **A0 (Path of Least Resistance)** proves that nature does not climb the wall—it bypasses the wall by opening a new orthogonal degree of freedom or crossing `0` into `-`:

1. **Wall 1 — 1D Axial Compression Wall $\to$ Helical Buckling Plateau (`1D -> 2D/3D`):**
   Compressing a straight 1D line ($\omega=0$) hits the exponential wall $E_{\mathrm{1D}}(\chi)=e^{+\chi}\to\infty$. Above $\chi>0$, buckling sideways into a helix with $v_\perp^\star=\sqrt{1-e^{-2\chi}}$ bypasses the wall and caps the total cost at the finite plateau
   $$E_{\mathrm{helix}}^\star(\chi)=\frac{3}{2}-\frac{1}{2}e^{-2\chi}\to 1.5.$$
2. **Wall 2 — 3D `+` Compression Wall $\to$ Signed `0`-Crossing Plateau (`+ -> 0 -> -`):**
   Forcing all winding to stay in the parent `+` region ($\theta=0$) hits the exponential wall $E_{\mathrm{stay},+}(\chi)=2e^{+\chi}+2e^{-\chi}\to\infty$. Opening the `0` horizon valve into a `-` child sphere ($\theta^\star\to 1$) bypasses the wall and caps the total 5-variable cost at the finite plateau
   $$E_{\mathrm{5\text{-var}}}^\star(\chi)\to 4-\frac{9}{2}e^{-\chi}+O(e^{-2\chi})\to 4.0.$$
3. **Wall 3 — Single-Axis Winding Wall $N^{3/2}$ $\to$ Orthogonal `(1, 1, 1)` Triad Ground State & Nuclear Hard-Core Potential:**
   Piling $N\ge 2$ windings onto a single channel costs $E_{\mathrm{single}}(N)\propto N^{3/2}>N$, a superlinear wall ($2^{3/2}\approx 2.828>2$, $3^{3/2}\approx 5.196>3$). Distributing integer windings $(n_1,n_2,n_3)\in\mathbb Z_{\ge 0}^3$ ($m_i=m_0 n_i^2$) across the 3 orthogonal axes of 3D space yields the exact **A0** minimum
   $$E^\star(n_1,n_2,n_3)=4(PV_0)^{1/4}m_0^{3/4}\sqrt{n_1 n_2 n_3},\qquad \lambda_i^\star=\frac{n_i^2}{(n_1 n_2 n_3)^{1/2}}.$$
   - $N=1$ `(1,0,0)` (1D loop) and $N=2$ `(1,1,0)` (2D sheet) have zero 3D volume ($V^\star=0$) and collapse along their unwound axes under 3D pressure.
   - **$N=3$ `(1,1,1)` is the unique minimum-cost 3D ground state:** $E^\star(1,1,1)=4.0$, cost per winding $4/3\approx 1.3333$, and isotropic aspect ratios $1:1:1$ (a round sphere).
   - **$N=4$ `(2,1,1)` is the first excited 3D state:** Because 3D space has only 3 orthogonal axes, a 4th winding cannot open a 4th spatial axis; it must double one axis to `(2,1,1)`, breaking spherical symmetry into a $4:1:1$ prolate state with higher cost per winding $E^\star(2,1,1)/4=\sqrt{2}\approx 1.4142>4/3$ (excitation gap $\Delta E=4(\sqrt{2}-1)\approx 1.6569$).
   - **Medium-range attraction + short-range hard-core repulsion:** Two `(1,1,1)` triads separated by $d>2R^\star$ attract via shared surrounding flow ($\Delta E_{\mathrm{flow}}=-2m_0/d<0$), achieving a bound minimum at contact $d=2R^\star$ ($E_{\mathrm{contact}}=8.0-1.0=7.0$). However, merging them into a single `(2,2,2)` bubble at $d\to 0$ costs $E^\star(2,2,2)=8\sqrt{2}\approx 11.3137>2E^\star(1,1,1)=8.0$, erecting an exact **short-range hard-core repulsion barrier** $\Delta E_{\mathrm{core}}=8(\sqrt{2}-1)\approx +3.3137$.

### 3. Quarantined numerical coincidences (not derived from `A0–A4`)

The following matches involve post-hoc formula/rung selection or external measured inputs ($H_0, G$) and are quarantined from the derived theorems:
- **Lenz (1951) $6\pi^5$ and medium dressing:** $6\pi^5\approx 1836.1181$ ($18.8\text{ ppm}$ from $m_p/m_e=1836.1527$) and $6\pi^5+e^{-3}/\sqrt{2}\approx 1836.1533$ ($0.35\text{ ppm}$). The factor $\pi^5$ is not derived from the discrete **A0** stretch functional.
- **Neutron–proton mass difference:** $2+\tfrac{1}{2}=2.5\ m_e$ ($-1.22$%) and $2+\tfrac{1}{2}+e^{-3}/\sqrt{2}\approx 2.5352\ m_e$ ($+0.167$% vs. $2.5310\ m_e$).
- **Baryon harmonic rung $K=7$:** $\eta_K=\tanh(\frac{\pi}{4}e^{-3K})$ gives $\eta_7\approx 5.955\times 10^{-10}$ at $K=7$, but the freeze-out rung $K=7$ was selected from the ladder after comparing to BBN.
- **CMB temperature projections:** $\sqrt{T_P T_H}e^{-2}\approx 2.625\text{ K}$ and $T_{\mathrm{crit}}e^{-(4/5)\chi_\star}\approx 2.728\text{ K}$ use measured $H_0$ and $G$ as inputs and an ad-hoc $-4/5$ exponent.

## Dynamical continuation of the same geometric picture

The [finite-parent construction](docs/nsc-atlas-offset-generator.md#coupled-finite-parent-receivers-and-characteristic-continuation)
now derives pressure, transmitted phase energy and geometric work from one
candidate action. Exact optical characteristics continue beyond the
exponentially stretched grid. The free-parent limit expands as
$\lambda=\cosh^2(a\tau)$, transfers phase energy into geometry, and has a finite
future causal reach. A symmetric branch transmits $4/9$ of incident energy to
each quiet receiver, with $1/9$ reflected. With radial-displacement inertia,
continuum wave-driven opening from zero radius obeys
$r(\tau)=\pi^2\tau^2-\frac{4}{3}\pi^4\tau^3+O(\tau^4)$, and the hybrid radial
inertia $M(r)=1/(1+r^2/3)$ satisfies both the finite-distance birth integral
and the finite-future optical-reach integral, producing a sharp last-echo
launch cutoff with logarithmic delay $-\ln(\delta)/\kappa$ and linear energy
vanishing $\propto\kappa\delta$.

Across the 5-variable cost family, eliminating $\chi$ yields the exact
parameter-free algebraic budget identity $9DVB=4(V+B/2)^3$ ($B=M-D/3$,
$B/V^2\to 16/27$ at deep compression), so any single share $V$ determines
$(D,M,\chi)$ with zero free parameters. For bubble fragmentation, same-sign
splits preserve equilibrium surface area $\sum_i 4\pi(R_{n_i}^\star)^2=4\pi(R_N^\star)^2$
while releasing budget $1-\sum_i(n_i/N)^{3/2}>0$. Solving the whole exterior
quadrupole wave at split onset gives the exact characteristic cubic
$p^3-3p-3=0$ ($p_+=\phi^{2/3}+\phi^{-2/3}\approx 2.1038034027$). At the first split's separated
endpoint $s=1.2$, the sibling field induces the daughters' multipole
preparation without a second supplied seed, and harmonicity ($\Delta\Phi=0$,
$\partial_{\zeta_x}v+\partial_{\zeta_y}v+\partial_{\zeta_z}v=0$) forces the
cascade through all three mutually orthogonal axes
$\mathbf e_z\to\mathbf e_x\to\mathbf e_y$ ($1\to 2\to 4\to 8$ lobes, with
exact directional derivative ratios $-1/2$ for $1\text{D}\to 2\text{D}$ and
$-2$ at $b=s$ for $2\text{D}\to 3\text{D}$), while the $B$-independent
$2\times 2$ kinetic matrix and outgoing-memory reaction drive the transverse
handoff dynamically, freezing out at $\chi_{\mathrm{freeze}}=3\chi_1=3$ ($\chi_{\mathrm{freeze,dyn}}=3\ln(1.2+2^{2/3})=3.075371$).

## Implementation and verification

All calculations and proofs live in [`nsc-atlas-offset-generator.md`](docs/nsc-atlas-offset-generator.md),
[`nsc_atlas_offset_generator.py`](src/recursive_horizons/nsc_atlas_offset_generator.py),
and [`test_nsc_atlas_offset_generator.py`](tests/test_nsc_atlas_offset_generator.py).
Run `python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --atlas-only`
to reproduce the complete non-local-to-local atlas readout.

The [previous theory narrative](THEORY-history-20261102.md) is retained as
history. The [current route](docs/research-index.md) identifies active work.
