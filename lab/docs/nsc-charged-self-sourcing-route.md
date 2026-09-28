# Reusing a charged-fermion self-sourced throat

A closer existing realization of the proposed matter–geometry feedback is
available. [Maldacena, Milekhin and Popov](https://arxiv.org/html/1807.04726v3)
construct a semiclassical Einstein–Maxwell throat supported by charged
massless-fermion vacuum energy. Their sections2 and5 derive the light
magnetic channels, the anomaly-corrected Casimir stress and its Einstein
matching. These results should be reused, not recomputed as an NSC discovery.

In their unit-charge example, magnetic flux q supplies q light two-dimensional
Dirac channels. A dynamical gauge Wilson line selects the energy-minimizing
antiperiodic phase. For the short-exterior-length regime, the imported
relations include

\[
r_e^2=\frac{\pi q^2G_N}{g^2},\qquad
E_{\rm quantum}=-\frac{q}{8\ell},\qquad
\ell=\frac{16r_e^3}{qG_N}.
\]

The energy includes the conformal anomaly; the same length follows from
Einstein matching, not only a variational estimate. The approximation uses
large flux, weak gauge coupling and specified exterior geometry. Longer
separations and mouth stabilization require their sections5.5–6. This is a
traversable-throat construction, not a demonstrated black-hole-to-child
cosmology or generic formation process. It is a prior theoretical example,
not an observed wormhole.

## The operator modification needed for this project

The current numerical throat operator is neutral. To use that physical
source, a local realization must contain an actual gauge connection,

\[
\mathscr D_{5,A}=i\Gamma^M(\nabla_M-iA_M),
\]

with its gauge action and coefficients obtained from the same functional.
This is a candidate U(1) extension of the existing carrier. It does not
derive the internal gauge algebra or identify the field as observed
electromagnetism. Its addition is motivated by the physical source being
matched and the existing charge/gauge targets; no scalar stress source is
introduced.

The old neutral angular spectrum has no physical kappa=0 spherical mode.
Magnetic flux changes the angular operator and its domain, allowing the
charged lowest-Landau channels used in the imported construction. The
previous neutral kappa=0 control cannot simply be relabeled as that result.
The full charged operator and flux state have to be represented explicitly.

## Minimal charged zero-mode compatibility

The [compact mass map](nsc-compact-mass-map.md) gives one Weyl zero mode
for a single bulk Dirac field with one chosen endpoint parity. By itself a
single charged Weyl field is not the anomaly-free four-dimensional Dirac
sector needed for the proposed gauge control.

Consider the already used two-copy spinor structure on a common development
background, now assigning the same unit gauge charge to both copies. This
field-content choice is explicit. Take a room-diagonal compact parity
Q=diag(eta_p,eta_c), with eta_p,eta_c=±1. Its internal reflection matrix
and compact even projector are

\[
R_Y=-Q\otimes\gamma^5,\qquad
P_{\rm even}=\frac{I_8-Q\otimes\gamma^5}{2}.
\]

The fifth-direction current matrix is I2 tensor i beta gamma5. Each
candidate endpoint projector is isotropic for this current. The additional
requirements are informative:

- A neutral scalar sheet link Phi tau1 tensor beta with an even canonical
  compact profile preserves the reflection exactly when {Q,tau1}=0.
- The four-dimensional cubic gauge and mixed gravitational gauge anomalies
  of the zero modes both cancel exactly when eta_p+eta_c=0, for these two
  equally charged copies and no additional light fields.

Thus the minimal charged candidate and its even neutral link both require
opposite compact parities. The two allowed zero-sector projectors are

\[
P_{\rm even}=\frac{I_8\mp\tau_3\otimes\gamma^5}{2}.
\]

They are the two invariant scalar-link sectors already identified in the
[spinor bridge](nsc-observable-bridge.md). The new check is their compact
parity and gauge-anomaly compatibility. The absolute orientation remains
unselected, and a different field content or link parity would define a
different candidate. This is not a derivation of the actual throat's scalar
coupling.

A right-handed charge-one field is counted as a left-handed charge-minus-one
field when computing anomalies. That bookkeeping does not assign opposite
physical charges to the two sheets or identify sheet exchange with charge
conjugation.

Before a link-induced mass is added, this pair has **one** four-dimensional
Dirac zero field. It has two Dirac fields at each nonzero compact level and
eight spinor components in the five-dimensional ultraviolet trace. These
different multiplicities must not be interchanged. An even link can lift
the zero pair; its actual profile and mass must be checked before applying
a massless-source formula. No value of Phi is chosen here.

## Same-operator gauge normalization

Use the existing proper-time window and the universal Dirac heat coefficient
from [Vassilevich, equation4.28](https://arxiv.org/html/hep-th/0306138v3).
For one complex unit-charge Dirac field, a4 contains (2/3)F_MN F^MN,
excluding the heat prefactor. The two bulk copies therefore supply

\[
\Gamma_{\Lambda,\nu}\supset
\frac{4(\Lambda-\nu)}{3(4\pi)^{5/2}}
\int d^5x\sqrt g\,F_{MN}F^{MN}.
\]

This is a checked leading contribution, not the complete renormalized gauge
coupling. Gauge-field variation and the reduction of its metric stress are
part of the same functional; an independent Maxwell or gravitational weight
is not fitted to hold open the throat.

## What must match before claiming an NSC self-sourced solution

The source construction supplies the reusable Einstein/Dirac feedback.
The following NSC identifications remain necessary:

1. The gauge representation, complete charged light sector and its five-
   dimensional determinant phase/boundary consistency. The zero-mode anomaly
   check is necessary and does not replace the full anomaly analysis.
2. The actual gauge-flux geometry and complete fermion return path through
   the outside region. A neutral throat or an arbitrarily closed axial circle
   is not that configuration.
3. The renormalized gravitational, gauge and vacuum coefficients from the
   common action. In particular the [UV map](nsc-torsion-uv-map.md) retains a
   volume term. It is not set to zero to obtain the imported asymptotically
   flat exterior.
4. The light-sector mass and physical state. The neutral P/AP controls and
   the massive compact tower do not already provide that source.
5. The stationary geometry and its admissible perturbations in the matched
   realization. Child expansion and cosmological predictions remain further
   questions for the same action and parameters.

This route makes the gauge sector part of the self-sourcing calculation.
The fundamental quantum-measure and fixed-point questions remain in the
programme, but the existing charged semiclassical construction should be
used when testing the source mechanism. There is no need to rediscover its
Einstein or Casimir calculation.

The [charged-sector record](../results/development/charged-sector.json) and
[focused runner](../scripts/check_nsc_charged_sector.py) verify only the new
parity, anomaly and charge/multiplicity normalization map. They do not
reproduce or claim to complete the imported wormhole solution. Every field
and source hash is checked by
`python3 scripts/check_nsc_charged_sector.py --check`.
