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
