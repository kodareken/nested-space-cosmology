# Current result — coupled geometric evolution, 2 November 2026

The existing damped oscillator already gives the signed lobe-area residue
and the K=7 example. Its same energy ledger gives a separate work residue
η_work=2η_area/(1+η_area²). At K=7 the area and work values are
5.955329e-10 and 1.191066e-9 respectively. The next connection to the
bubble split is a transition that produces the two lobe histories and
specifies how the receiver records them. See the coupled oscillation
section of the [offset owner](nsc-atlas-offset-generator.md).

The five-variable geometry now gives an exact relation across its whole
budget curve. With D the stretch share, M the load share, V the filament
share and B=M−D/3, stationarity gives 9DVB=4(V+B/2)³. The shares recover
the generating compression, and B/V² tends to 16/27 at deep compression.
See the [budget proof](nsc-atlas-offset-generator.md#exact-relation-between-the-three-budget-shares).

The coupled hybrid birth has a last-echo launch cutoff near .226942387
at the mouth and .116502935 at a parent detector .01 material units away.
Arrival delay diverges logarithmically and the returned energy ratio
vanishes linearly. A null boundary enters the receiving arc near .635108
and tends to proper distance 1.108043 from the mouth. See the
[proof and reproduction](nsc-atlas-offset-generator.md#last-echo-from-a-coupled-hybrid-opening).

The continuum mouth is timelike at positive offset. An incident packet
at .02 returns to the parent detector at .039385342912 and retains
.502902410597 of its energy at the mouth. The conserved budget gives
separate integral criteria for finite-offset birth and finite future
optical reach. A hyperbolic radial inertia satisfies both and gives a
free-room optical reach .261158980063, conditional on that inertia.
See [mouth geometry and inertia criteria](nsc-atlas-offset-generator.md#shared-mouth-geometry-and-the-returned-signal).

The continuum opening with radial-displacement inertia gives
r=π²τ²−(4/3)π⁴τ³+O(τ⁴). At τ=.05 the receiving arc holds .299289539
energy and has radius .016736471416. Curvature diverges at birth while
physical-distance coordinates admit a continuous curve metric extension.
See the [continuum derivation](nsc-atlas-offset-generator.md#continuum-wave-at-the-opening).

The zero-radius calculation finds infinite kinetic distance under logarithmic
inertia. A separate radial-displacement inertia admits an affine continuation
of the same triangular winding with r=2.9718260669 τ²+O(τ⁴). At τ=0.1 the
new arc carries energy 0.306711998348, with budget defect below 2e-13 on the
fine seed comparison. The arc is admitted by the reduction and its spatial
and causal completion remains open. See the zero-radius derivation in the
[offset owner](nsc-atlas-offset-generator.md#zero-radius-formation-with-finite-radial-inertia).

The active connection is **one inherited phase -> finite enclosing pressure ->
compression and receiving energy -> geometric response**. The
[coupled derivation](nsc-atlas-offset-generator.md#coupled-finite-parent-receivers-and-characteristic-continuation)
extends the existing offset module; no new particle fields or fitted cosmic
fractions enter it.

## Mathematical connections now established

- One action contains the A0 phase stretch cost, its phase inertia and an
  explicit logarithmic-boundary inertia. Finite parent pressure is
  P_i=E_i/(3V_i); the enclosing radius is free.
- Parent pressure drives inner compression. The receiving field supplies
  opposing pressure. Regional transfer, geometric work and their total
  conserved budget follow from that action. No P(t), chi(t) or damping
  coefficient is independently prescribed.
- Exact characteristic coordinates continue the calculation beyond the
  exponentially stretched material-grid representation. Quiet regions remain
  exactly quiet until their characteristic front arrives.
- A symmetric parent/two-receiver contact reflects energy 1/9 and transmits
  4/9 to each receiver. The contact scattering matrix is orthogonal.
- The same action has an exact free-parent limit lambda=cosh^2(a*tau),
  a=sqrt(Q0/6), Q0=4*pi^2. Phase energy decreases as Q0/lambda while geometric
  motion receives the complementary budget. Total future optical reach is
  0.3898484006, and proper causal-boundary length tends to 0.1949242003.
  Extending this scale to the flat 3D completion gives late curvature
  Lambda_curvature=8*pi^2 in the chosen units, without an added vacuum
  potential. This is a metric result, not an SI cosmological prediction.
- The receiving-sphere cost with the declared log-strain inertia predicts
  breathing/shear frequency ratio 2, conditional on fixed inherited winding.

## Finite calculations and numerical boundaries

The [small reference record](../results/summaries/nsc-finite-parent-receiver-v1.json)
and existing CLI --evolution-only reproduce the chain through offset 2.5 and
the pre-reflection symmetric fork through offset .9, at three optical step
sizes. Eight independent jobs use eight worker processes. No raw histories
are committed; the CLI regenerates them in memory.

In the fine chain the first receiver expands after compression; its available
volume reaches about 18.90255 and phase energy .452975. The next region
receives a genuine wave after its front arrives, retains phase energy .039102,
and also returns energy. Net transfer and geometric work are reported
separately. A negative net transfer does not mean that no energy arrived.

The fine chain's sampled total-budget defect is about 3.36e-4 against initial
39.4784. Adjacent phase-energy and volume differences contract. These are
empirical numerical controls, not certified continuum error bars. The fine
pre-reflection fork has phase energy about 1.054502 in each receiver, exact
symmetry and a total-budget defect about 7.29e-8. Later reflected fork
histories have noncontracting jump-interpolation errors and are not admitted
numerical claims. That diagnosis is retained; it does not reactivate a grid
campaign.

## What is still a premise or open connection

The regional tree, receiving domains and unit reference volumes are admitted
initial geometry. The volume-to-path stretch law and geometric inertia are
candidate dynamical premises; the offset recurrence does not uniquely select
them. A dynamically born child, a spatial embedding of several spherical
children, and a one-way null interface are not proved by this construction.
The regional +/- orientation convention remains exterior/horizon/interior,
with offset continuing on both sides. Local characteristic signs alone do
not determine the global event horizon.

The phase part has work-conjugate pressure/density ratio 1/3. Its budget
shares have not become cosmological matter/vacuum fractions. Charge,
baryon/photon conversion, an SI ruler and a selected observation epoch remain
open. No e^-3 floor, K=7 rung or measured abundance selects this evolution.

Next: derive the causal/contact geometry and admissible branch transition
from the generated paths, reusing the conserved work and characteristic
transfer just established. Do not resume GH, horizon-precision, Dirac, or
cosmic-number campaigns. The existing manuscript source includes the new
connection; the PDF is not rebuilt in this task.

The previous equilibrium, spectral and action-based constructions remain
in the [offset owner](nsc-atlas-offset-generator.md), the
[retained library](research-index-history-20261102.md) and Git history.
No progress percentage is used.
