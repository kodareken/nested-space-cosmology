# Current result — foundation reset, 2 November 2026

## Active: enclosing deformation drives the proposed child offset

The [source relation](nsc-atlas-offset-generator.md#active-start-parent-compression-not-an-isolated-offset)
now starts inside an evolving parent, itself enclosed. For corresponding
local deformation histories it computes
chi_pc=(log(Jp/Jp0)-log(Jc/Jc0))/d. This is dimensionless relative
log-length compression; its identification with the generating offset is
the proposed physical connection, not an already derived clock or phase.

Proved: chi_gc=chi_gp+chi_pc, common frame changes cancel, and common
compression adds no relative offset while preserving ancestor inheritance.
This additivity is an algebraic log-ratio telescoping identity (bookkeeping,
not dynamical evidence for the mechanism). The increment is obtained from the
deformation history rather than an independently supplied alpha. Analytic
three-level controls give u,2u,3u with zero sampled composition residual; they
are not formation histories.

The intrinsic-phase calculation compares the 3D volume-weighted functional
Q_vol(A)=det(A) trace(A^-T M A^-1)/2 against the A0 1D generator least-stretch
cost Q_stretch(A)=trace[(A^T A)^(-1/2) C^T C]/2 derived by Cauchy-Schwarz on
the generated 1D phase paths. Because Q_vol dilutes 1D phase stretch by 3D
volume lambda^3, Q_vol falls as exp(-chi) in d=3 and has no positive isotropic
holding radius under enclosing compression P>0 (a lambda->0 break, continued by
scale inversion lambda_child=lambda^-1=exp(chi)). By contrast, Q_stretch scales
as exp(+chi) in every dimension d>=1 and resists compression.

Minimizing Q_stretch + P V under enclosing load P V0 = 1 yields exact outcomes:
- Isotropic 3D winding m=(1/2,1/2,1/2): uniform contraction to
  lambda_* = 2^(-1/4) ≈ 0.8408964153 with Virial ratio Q_stretch^* / W_ext^* = 3.
- Raw two-accumulator moment m=(1/4,1/4,1/2): localized prolate neck equilibrium
  with aspect ratio lambda_3^* / lambda_perp^* = 2 (s_* = (1/3) log 2 ≈ 0.2310490602).
- 1D axial seed m=(0,0,1/2): transverse neck-pinch break (lambda_perp^* -> 0).

Applying the break rule across the axial neck motion s(tau)=log cosh(sqrt(2/3) tau)
(where parent transverse radius a=exp(-s/2)->0 and normalized velocity hits
0.9998888764 at tau=6) via throat inversion R_child = a^(-2) = exp(s) = cosh(sqrt(2/3) tau)
continues the trajectory into an expanding child frame obeying
(R_child_dot/R_child)^2 + 2/(3 R_child^2) = 2/3, with H_child -> sqrt(2/3) ≈ 0.8164965809
and effective Lambda_child = 2.

Both Q_vol and Q_stretch are SO(3)-invariant and produce zero rotational torque
from pure compression. When paired with an angular turn rate omega=2*pi*alpha per
unit chi, decay-per-radian is 1/(2*pi*alpha) (0.4166730505 for golden-angle alpha
with opposite-lobe asymmetry eta=tanh(1/(4*alpha))≈0.5746972340, versus
2*ln(phi)/pi ≈ 0.3063489625 for a golden spiral).

The connected pipeline (`--buckling-only`) now derives the entire cycle from a
single 5-variable least-resistance cost functional (`unified_least_resistance_cycle`)
with openly declared axioms (A0 least resistance, A1 1D+offset quadratic winding,
A4 turn conservation across `0`) and zero fluid-physics assumptions (no parcel
incompressibility, no surface tension, no assumed V_minus = V0(1-exp(-3*chi))):
- **Single-cost 5-variable global minimum:** Minimizing E(lambda_3, tilde_lambda_1..3, theta; chi) = E_+ + E_-
  gives sqrt(theta^*) = (sqrt(16*exp(2*chi)+17)-3)/(4*exp(chi)),
  lambda_3^* = (1-theta^*)*exp(chi), and tilde_lambda_1^* = tilde_lambda_2^* = tilde_lambda_3^* = sqrt(theta^*).
  The `+` parcel shrinks in volume (V_+^* ~ (3/2)*exp(-2*chi) -> 0) while stretching
  into a 1D filament above chi = ln(2)/2 ≈ 0.3465735903 (where theta^* = 1/2
  and aspect ratio = 1, growing as (3/2)*exp(+chi) -> infinity); the `-` region
  rounds into a 3D sphere (aspect ratio 1) by AM-GM on Q_stretch; and the `-`
  drop's derived volume share V_-^*/(V_+^*+V_-^*) reaches 84.63% at chi=1,
  97.73% at chi=2, and 99.95% at chi=4.
- **Nested signed-zero crossing map & ~95% exponential horizon (`nested_zero_crossing_map`):**
  Iterating the signed-zero crossing axiom (sign(det F_k) = (-1)^k) across
  nested levels places local space k between an outer parent `0` valve (inflow)
  and inner black-hole `0` valves (outflow), preventing multiple compressions
  from emptying local visible matter. At one 3D e-fold (chi = 3), the 3D
  exponential horizon is 1 - exp(-3) ≈ 95.0213% dark sector (split 3:1 by the
  3D Virial lock into 71.27% DE-like and 23.76% DM-like) with visible floor
  exp(-3) ≈ 4.9787%; the 5-variable single-cost functional gives
  (68.21% DE, 23.24% DM, 8.54% Vis) at chi=2.5 and (70.72% DE, 23.79% DM, 5.49% Vis,
  94.51% combined dark sector) at chi=3.0. Bulk orientation inside each child is
  100% single-signed, whereas single-mode damping eta = tanh(pi*gamma/(2*omega))
  is O(0.086-0.575) and requires omega/gamma ≈ 1.04e9 (or ~21.22 e-folds) to
  reach 6.1e-10.
- **Moving spectral band [Planck cutoff, Absolute-Zero horizon] and cyclic vortex-diffusion ring (`spectral_band_and_baryon_harmonic`, `branching_cosmos_chain_evolution`, `cyclic_vortex_diffusion_ring`):**
  Both the longest wavelength lambda_max(chi) = lambda_max_0 * exp(-chi)
  (Absolute Zero, as wide as the room) and the shortest wavelength before `0`
  crosses to `-` lambda_min(chi) = lambda_min_0 * exp(-chi) (Planck cutoff)
  shift by exp(-chi) across nests, while the spectral asymmetry curve
  eta(lambda) = tanh((pi/2)*(lambda/lambda_max)) spans from tanh(pi/2) ≈ 0.917152
  at the horizon mode down to 6.1e-10 at lambda_B/lambda_max ≈ 3.88338e-10
  (~31.26 octaves below the horizon). Across generations, child black-hole
  expansion drains the parent ('O' -> 'o') while children balloon ('o' -> 'O',
  'Oooo -> oOOOooooooooo'), and on a closed ring ('oo / OO oO / oO Oo / oo')
  diffused still air ('o') re-accumulates and spins back up into a giant
  vortex ('O') with zero conservation residual.
- **Non-local to local Atlas readout (`nonlocal_to_local_atlas_readout`, `--atlas-only`):**
  Solving the zero-parameter two-valve fixed point where the 5-variable
  single-cost visible filament share equals the 3D volumetric floor exp(-3) ≈ 4.9787%
  gives chi_* ≈ 3.06696835 (71.1441% Dark Energy-like, 23.8772% Dark Matter-like,
  95.0213% combined dark sector, 4.9787% visible filament; at integer chi=3,
  70.8827% DE, 23.8134% DM, 5.3039% Vis). Along the discrete 3D volumetric
  e-fold ladder eta_K = tanh((pi/4)*exp(-3*K)) anchored at the derived
  spaghettification onset exp(-2*chi_spaghetti)=1/2, K=7 (21 e-folds) yields
  eta_7 ≈ 5.955329e-10 inside the observed BBN window [5.8e-10, 6.5e-10].
- **Circle-pendulum duality and why the radius of space must cross `0` (`signed_circle_pendulum_duality`):**
  Seen edge-on, the Layer-1 closed circle z(tau)=R0*exp(i*omega*tau) is the
  pendulum R(tau)=R0*cos(omega*tau), swinging +R0 -> 0 -> -R0 -> 0 -> +R0.
  Forbidding R < 0 clips 50% of the circle (pi*R0 of 2*pi*R0) so no loop can
  close; multiplying by -1 = exp(i*pi) is a half-turn (in 3D, P = -I_3 with
  det(P) = -1), so two crossings give P^2 = +I_3 and (-1)*(-1) = +1. Treating
  radius as unsigned size r = |R| >= 0 folds R < 0 onto r > 0 and creates an
  artificial velocity cusp Delta(r_dot) = +2*R0*omega at r=0, whereas signed
  R(tau) crosses 0 into - smoothly with Delta(R_dot) = 0 and R_tt(tau_0) = 0.
- **Vacuum bubble crystallization, surrounding-medium flow, 180-deg annihilation, resonance attraction, and discrete bursts (`vacuum_bubble_crystallization_and_mass_spectrum`):**
  Minimizing E(R) = m0*n^2/R (surrounding flow r>=R) + P*V0*R^3 (cavity work)
  crystallizes equilibrium bubbles R_n^* ~ |n|^(1/2), E_n^* ~ |n|^(3/2) with 75%
  of mass-energy in the surrounding medium flow and 0% in a point singularity.
  Opposite bubbles (+1, -1) cancel to n_net=0 and release two equal 180-deg
  ripples (E_1^* each); shared-flow cross-energy averages to 0 off-resonance
  and attracts as -1/d^2 when omega1=omega2; integer loop closure z(1)=z(0)
  forces sharp line bursts rather than a smear.
- **Derived 3-accumulator volumetric bubble ladder vs. quarantined numerical coincidences (`derived_three_accumulator_bubble_ladder`):**
  Iterated A1 offset accumulation yields binomial windings W_k(n) = C(n, k) and
  3D equilibrium energy E_3D^*(n) = 4*(P*V0)^(1/4)*m0^(3/4)*sqrt(C(n,1)*C(n,2)*C(n,3)).
  This proves that n=1 (W_2=W_3=0) is a 1D loop, n=2 (W_3=0) is a 2D sheet, and
  n=3 (W_1*W_2*W_3=9>0) is the minimum integer supporting a 3D volume, with
  untuned ratio E_3D^*(3)/E_1D^*(1) = 6 and 3D ladder 1 (1D) : 6 : 19.60 : 44.72 : 84.85
  (missing 1836.15 by ~306.03x without a continuum/nesting scale factor).
  Post-hoc numerical coincidences (Lenz 1951 6*pi^5 + exp(-3)/sqrt(2), n-p gap
  2+1/2+exp(-3)/sqrt(2), rung selection K=7 on eta_K, CMB projections using
  measured H0/G and -(4/5)*chi_* in `thermodynamic_gravity_and_cmb_temperature`,
  and the 71.14%/23.88%/4.98% split vs LambdaCDM 68.5%) are strictly quarantined
  in `numerical_coincidences_not_derived`.
- **Universal exponential-wall bypasses, `(1, 1, 1)` orthogonal triad winner, and nuclear hard-core repulsion (`exponential_wall_bypasses_and_triad_winner`):**
  Whenever compression or winding hits an exponential or superlinear resistance
  wall, A0 bypasses the wall by opening a new spatial dimension or crossing `0`
  into `-`: (1) the 1D compression wall E_1D(chi)=exp(+chi)->infinity bypasses
  into the helical buckling plateau 3/2 - (1/2)*exp(-2*chi) -> 1.5; (2) the 3D
  `+` compression wall E_stay,+(chi)=2*exp(+chi)+2*exp(-chi)->infinity bypasses
  via `0` into the `-` child sphere plateau 4 - (3/2)*exp(-chi) -> 4.0; (3) the
  single-axis winding wall N^(3/2) > N splits integer windings across the 3
  orthogonal axes with E^*(n1,n2,n3) = 4*sqrt(n1*n2*n3), proving that `(1,0,0)`
  and `(1,1,0)` have zero 3D volume, **`(1,1,1)` (N=3) is the unique minimum-cost
  3D ground state** (E^*=4.0, cost/winding 4/3≈1.3333, round sphere 1:1:1), N=4
  forces prolate `(2,1,1)` (4:1:1, cost/winding sqrt(2)≈1.4142 > 4/3), and two
  `(1,1,1)` triads attract at medium range (-2*m0/d, bound contact E=7.0 at
  d=2*R^*) while erecting a **short-range hard-core repulsion barrier +3.3137**
  against merging into `(2,2,2)` (8*sqrt(2)≈11.3137 > 8.0).
0. **Two-force zero-crossing:** When inward compression F_in exceeds outward
   resistance F_out (F_net = F_in - F_out > 0), R_tt = -F_net < 0 drives R(tau)
   smoothly through R=0 into R<0. Busting back into `+` (R>0) is uphill against
   F_net>0, and stopping at R=0 requires infinite deceleration; crossing `0`
   into `-` (R<0) is the unique finite-cost path, flipping det(F)=(R/R0)^3 from
   `+` to `-` and turning implosion R_dot < 0 directly into expansion
   r_minus_dot = -R_dot > 0 of a new sphere r_minus = -R > 0 from scratch.
1. **Helical buckling:** A 1D seed compressed axially to lambda=exp(-chi) under
   lateral confinement P_perp has arc-length factor sqrt(lambda^2+(r_perp*omega)^2).
   Once lambda < lambda_crit=(Q0/P_perp)^(1/3), the straight seed (omega=0)
   loses transverse stability (d^2E/dv_perp^2 < 0) and buckles into a helix with
   derived turning speed r_perp*omega^* = sqrt(lambda_crit^2 - lambda^2) > 0.
2. **Damped `+`/`0`/`-` horizon oscillation:** Coupling q(tau)=A0*exp(-gamma*tau)*sin(omega*tau)
   to child geometry growth via E_geom_dot = 2*gamma*q_dot^2 conserves
   E_wave(tau) + E_geom(tau) = A0^2*omega^2/2 identically and yields the exact
   invariant `+`/`-` lobe asymmetry eta = (I_+ - I_-)/(I_+ + I_-) = tanh(pi*gamma/(2*omega)) > 0
   with net remainder Delta_I = A0*omega/(gamma^2+omega^2) > 0.
3. **Goldilocks equilibrium:** As E_geom(tau) populates the child's 3D stretch
   budget against enclosing load P V0, the child holding scale
   lambda_child(tau) = (E_geom(tau)/(3 P V0))^(1/4) expands from 0 and stabilizes
   at lambda_Goldilocks = (E_total/(3 P V0))^(1/4) (98.40% complete by tau=4 for
   omega=2*pi, gamma=ln(2)/2).

## Retained geometry and transport controls

The retained offset history has an explicit three-dimensional metric
completion, conditional on its distance kernel, retained history and
completion rule. It does not select those rules or a physical clock.

Opposite integer windings now have an explicit closed-loop contraction:
Q(h)=4 pi^2 R^2 m^2 (1-h)^2, the sharp least squared variation for the
declared midpoint and endpoint constraints. Two balanced positive outputs
have ratio 1 and angle 180 degrees. Balance also admits three or more
outputs; balance alone does not select the channels.

The next calculation now adds the circle's geometric inertia explicitly.
Its resulting phase equation has an exact conservative two-characteristic
solution. On an open drawing coordinate the two packets separate at
dimensionless offset 1, each carrying half the initial budget. By offset 2,
the original interval has transferred its entire budget to the two outside
intervals; the receiving amounts agree with integrated boundary flux. On the
periodic coordinate the spatial variation vanishes at offset 1/2 while
the full budget remains in motion; the initial state returns at offset 2.
This distinguishes transfer into an open domain from transient flattening.

The inertial rule, zero initial velocity and open receiving geometry are
declared candidate assumptions, not selected by the recurrence. Their next
connection must come from the enclosing deformation. No physical energy scale,
electron, horizon or receiving child is identified. Both geometry-only
CLI routes run without the conditional parent and write no result files.

The [previous calculation record](current-result-history-20261102.md)
preserves all numbers, domains and unsuccessful routes. In particular,
the old parent calculation produced finite internal-volume growth, not a
formed child or a derived present-day vacuum fraction. No new PDF is built.

No progress percentage is used. Progress means a stated connection derived
or a missing assumption isolated so the next calculation can address it.
