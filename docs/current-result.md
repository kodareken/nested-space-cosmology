# What the calculations show

**Retained-state integration inputs: PASS.** Two further input collections
are now complete for the full covariance construction:

- [Compact subgap panels](nsc-pg-retained-subgap.md): all **38 signed families**
  in the 20 compact groups; largest summed quadrature estimate
  $3.14\times10^{-11}$.
- [High-energy tail inputs](nsc-pg-retained-tail.md): all **61 signed families**
  for the remaining 31 groups. The largest covariance representation bound is
  $9.73\times10^{-10}$ and the Gram bound is $1.95\times10^{-9}$.

Five groups require the numerical transition at 320 instead of 160. This
changes the numerical integration split, with the same physical operator and
parameters. The [precision audit](nsc-pg-tail-moment-audit.json) independently
recomputed all 29 new energy tables with 50-digit process-local precision;
their stored numerical values were identical.

**Full retained C1b remains OPEN at covariance assembly and its error checks.**
The completed group-13 result below remains the first full-energy massive
covariance. PDF v0.25.0 and the metric-evolution gate are unchanged.


**Latest development — first full-energy massive PG covariance: PASS.**
The [new calculation](nsc-pg-group13-covariance.md) completes the retained
$m=\pi/2$, zero-angular-eigenvalue group on the common PG Cauchy surface.
The [record](../results/development/nsc-pg-group13-covariance.json) contains
the covariance, independently integrated CAR matrix, and lesser/greater/Keldysh
blocks, with the complete energy account and its numerical error checks.

| Quantity | Result |
|---|---:|
| CAR residual | $3.27\times10^{-11}$ |
| Covariance eigenvalue range | $[0.03570181,0.96429819]$ |
| Original-probe / bulk-probe correlation norm | $0.31816523$ |
| Finite-energy covariance refinement | $2.81\times10^{-15}$ |
| Tail mode-order change | $5.33\times10^{-15}$ |
| Tail spatial-grid change | $1.26\times10^{-12}$ |
| Fourier/normalization covariance-tail bound | $5.70\times10^{-10}$ |

The primary covariance uses the **calculated** CAR Gram; no identity correction
is used to fill a truncated frequency table. The eight probes inspect the full
mode-defined field and its unresolved-bulk correlations. They do not form a
closed eight-mode evolution model.

**Full retained C1b remains OPEN:** 31 massive groups remain. The locked action,
scales, seeds and all 209 earlier scientific files are unchanged. The existing
fourth-order reference/local allocation is preserved; no stress or metric
timestep is assigned here. **PDF v0.25.0 remains unchanged.**

Earlier entries below retain their checkpoint scope.


**Previous checkpoint — massive PG mode/Green construction: PASS.**
The [new derivation](nsc-pg-massive-mode-resolution.md) and
[record](../results/development/nsc-pg-massive-mode-resolution.json) evaluate
all **32 massive groups** on the complete existing node tables: **3,584
positive-energy/angular-sign combinations**, with signed-frequency partners
fixed by the existing operator map.

The input is the affine horizon covariance and inherited incoming occupation.
Complex Jost amplitudes give the exterior modes; the matched inner preparation
gives their transmitted components and independent interior horizon partner.
An open infinity channel remains present even when its occupation vanishes.

| Check | Maximum | Tolerance |
|---|---:|---:|
| Complex scattering isometry | $1.12\times10^{-12}$ | $3\times10^{-11}$ |
| Physical source law | $0$ | $3\times10^{-11}$ |
| Source CAR violation | $2.66\times10^{-17}$ | $3\times10^{-11}$ |
| Whole-line Green jump, relative | $8.29\times10^{-10}$ | $3\times10^{-9}$ |
| PG fields under outer-radius doubling | $7.79\times10^{-10}$ | $3\times10^{-9}$ |

Seven focused tests pass, including rejection of an omitted open infinity
channel. The old seed generations and all 201 previous scientific files are
unchanged. **Full C1b remains OPEN at the packet energy integral and its
contact/tail account.** A finite mode table is not used to fill the continuum
CAR identity. PDF v0.25.0 and the metric-evolution gate remain unchanged.

The earlier entries below retain their individual checkpoint scope.


**Previous checkpoint — full compact restart: H1 PASS.** The
[versioned restart](nsc-compact-matched-restart.md) now regenerates every
retained compact frequency from the original horizon/incoming data: 20 groups,
1,280 X blocks, 1,152 angular-partner Y blocks and their signed-frequency
partners. The historical seed is used only after construction, as a control.

| Start convention | Explicit control | Result |
|---|---|---|
| `historical_radial_r_h` | Immutable original compact seed | Maximum X residual `4.11e-14` |
| `matched_delta_q` | New `compact-seed-v2-matched-delta-q` checkpoint | Cold regeneration reproduces the stored artifact byte for byte |

The new matched control uses the corrected interior Frobenius argument and
shares the same inherited propagator, source and scattering with the historical
comparison. The two conventions are not made numerically equal. Their largest
X difference on the full grid is `2.73e-4`; it remains a diagnostic of the finite
start convention.

The [record](../results/development/nsc-compact-matched-restart.json) verifies
unitarity `5.97e-14`, matched modal normalization `5.99e-14`, T feedback
`4.37e-16`, and the prior independent matched modal control `9.20e-12`
(tolerance `3e-11`). The loader requires an explicit convention and hash.
The original seed, source file and previous result records are unchanged.

**H2 and full C1b remain OPEN.** The complete compact seed-frequency table is
now available under a consistent matched convention. It still needs the massive
global PG mode transform and current/Plancherel completeness before spatial
covariance assembly. No massive spatial C, transmitting endpoint jets, new
stress or metric evolution is inferred from this checkpoint. **PDF v0.25.0
is unchanged.**

**Preceding development — paired horizon construction and a compact normalization correction.**
The [new map](nsc-horizon-paired-pg-map.md) normalizes the magnetic angular
spinors, their gauge transition and their signed pairing. Under the already
selected affine-horizon preparation, orthogonality and angular-sign conservation
fix $Z(E)=0$. The horizon-partner correlations inside each angular sector remain.
The other marginal $Y(E)$ is computed by its own massive mode propagation,
without using stored seed covariance as input.

The [record](../results/development/nsc-horizon-paired-pg-map.json) evaluates
one existing frequency in each of the 32 massive groups. Magnetic angular
residuals are below `1.16e-13`; matched frame and mode residuals are below
`6.44e-12` (tolerance `3e-11`). These are finite modal maps on the trapped
probe chart, not a completed spatial PG covariance or all-frequency state.

**Finite compact collar mismatch identified.** The stored interior uses
$\delta_q=q_{h,\mathrm{geom}}-q_{\mathrm{geom}}$, with
$q_{\mathrm{geom}}=\pi/2+\arctan\rho$; this coordinate is distinct from the
fixed magnetic flux. Since $|\rho-\rho_h|=r_h^2\delta_q+O(\delta_q^2)$,
the interior Frobenius distance is $\sqrt{2\delta_q/\kappa_h}$.
The historical compact initializer additionally divides by $r_h$.
That factor is correct for an exterior radial offset, but not for this interior
coordinate. The new owner uses the matched normalization and leaves the old
files untouched.

All twelve angular-only controls recover their stored X. For all twenty
compact controls, the legacy frame recovers X while the matched frame changes
it by norms between `3.88e-5` and `1.47e-4`. This isolates the discrepancy to
the finite initial frame. Both conventions have the same zero-collar limit;
no new physical interaction or fitted coefficient is introduced.

**C1b-H remains OPEN.** The compact retained-state checkpoint must be reconciled
under the matched collar convention, and the global signed PG spectral/current
resolution must be completed before spatial covariance assembly. No corrected
modal counterpart is silently substituted into the locked seed. C1b-M2 and
C2–C4 remain unexecuted. **PDF v0.25.0 is unchanged.**

**Preceding development — massive signed-spin map and data-sufficiency gate.**
The [new operator map](nsc-massive-signed-preparation.md) connects the imported
paired angular representation to T's current basis. Its coefficient identities
are exact; the maximum numerical residual is `6.06e-15` (tolerance `3e-11`).
T normal-current and signed-frequency checks have zero residual.

Fourteen groups admit simple radial signed-frequency maps: the twelve
angular-only groups and two compact-only groups. The eighteen groups with
both compact mass and angular potential need the paired representation for
this local signed-frequency symmetry. Their pointwise two-component axis flip
carries a required basis-connection term; dropping it would change the operator.

**C1b-M OPEN, with a verified information certificate.** Across all 1,152 stored
mixed-group seed blocks, two algebraic completions preserve the observed
marginal exactly while differing in the missing angular/signed partner by
Frobenius norm `sqrt(2)`. CAR violations from numerical roundoff stay below
`4.18e-14`. These are information controls, not physical state choices.
The [record](../results/development/nsc-massive-signed-preparation.json) names
the required opposite-angular marginal and cross block, to be generated by
the horizon/infinity-to-paired-PG modes in the current/gauge basis. The
operator involution does not by itself prove physical state symmetry.
No massive PG covariance or transmitting branch jet is assigned. Full C1b
and C2–C4 remain OPEN; **PDF v0.25.0 is unchanged**.

**Preceding development — physical LLL PG Cauchy preparation.** The
[LLL preparation](nsc-pg-lll-preparation.md) uses the original horizon pair and
inherited incoming occupation to construct a covariance on the same common
PG-time slice as the retarded memory operator. It retains the exterior outgoing
partner and packet/bulk correlations; stored seed blocks are comparison data,
not copied spatial inputs.

The [record](../results/development/nsc-pg-lll-preparation.json) gives seed
reconstruction `5.56e-17`, spatial CAR-Gram residual `7.57e-15`, and quadrature
refinement `4.34e-12` (tolerance `3e-11`). The independent characteristic
response matches the locked LLL resolvent to `1.53e-14`. Its seven-probe
covariance has eigenvalues between `0.11534` and `0.88466`; the packet-to-sampled-
bulk correlation norm is `0.50269`. These are quantum correlation matrices,
not stress or energy fractions.

**LLL C1b PASS; full C1b and C2–C4 OPEN.** The 32 massive groups still require
the signed-frequency spin/angular mode map in T's current basis, or the full
paired four-spinor PG description and its state correlations. Degeneracy and
copy count do not supply that map. Full transmitting `EndpointBranchJets`,
remaining boundary variation and physical Weyl mismatch stay unevaluated.
This development record follows the unchanged **PDF v0.25.0**; no new paper
version or metric evolution is claimed.

**Preceding result — causal memory and temporal metric transfer.** The
[new CTP connection](nsc-transmitting-ctp-resolvent.md) places the computed
packet response in the convention $G_J^R=-\mathcal R$, $D_J^R=-\mathcal K$.
It preserves the energy-dependent memory and calculates the four metric
vertices between distinct frequencies, using the same transmitted Dirac
operator and archived incoming fields.

All 33 channels pass: driven-versus-weak-form residual `5.57e-13`, adjoint
frequency-exchange residual `6.19e-13`, and inverse-memory Dyson residual
`5.14e-16` (tolerance `3e-11`). The zero-transfer limit reconstructs the
locked static derivatives; reverse retarded response and its metric transfer
remain zero in the trapped probe cell. The harmonic source is a linear-response
probe, not a selected metric history or physical frequency prediction.

**Retarded-memory C1 component PASS; full C1 and C2–C4 OPEN.** The
[record](../results/development/nsc-transmitting-ctp-resolvent.json) retains the
exact missing preparation: the physical covariance on a common PG-time slice,
including its packet/complement cross correlations. That state, full history
reconstruction and the KS endpoint pullback are required before the memory
vertices become `EndpointBranchJets`. The existing consumers leave the
remaining boundary derivative and physical Weyl mismatch unevaluated.

**Preceding result — domain-correct transmitting response: Z2a PASS.** The
[cross-resolvent calculation](nsc-transmitting-cross-resolvent.md) evaluates
all 33 retained channels on the inherited whole-line Dirac domain. Continuous
seam data cancel the individual sharp-projector delta terms. Independent
parent/child spatial probes then give an energy-dependent inverse response
$\mathcal K(z)=[J^\dagger(H_D-z)^{-1}J]^{-1}$.

In the trapped probe cell, the retarded parent-to-child response is nonzero;
the reverse retarded response is zero by causal direction. This is a response
statement, not a choice $B=0$. The adjoint residual is `3.53e-14`, the true
resolvent-identity residual `1.18e-13`, and seam mismatch zero (tolerance
`3e-11`). The four static metric kernels also pass: independent weak-form
agreement is below `6.56e-13`; the finite-difference check is `1.01e-8`
(tolerance `3e-8`). The [record and payload](../results/development/nsc-transmitting-cross-resolvent.json)
store the response, its derivatives and uncompressed quadrature fields.

**Static Z2b PASS; time-history CTP bridge and Z4 OPEN.** The computed
resolvent derivatives need time/spectral reconstruction, physical state/memory
data and the endpoint-coordinate pullback before they become
`EndpointBranchJets`. The remaining boundary action and physical two-sided
Weyl mismatch are unevaluated. No stationary history or stress is assigned.
Moving-neck jets remain out of scope for the fixed coordinate seam.

**Preceding definition — common-time bulk split: Z1 PASS.** The
[new bulk domain](nsc-common-time-bulk-split.md) uses the inherited PG time
and spatial supports $\rho>0$ and $\rho<0$. Its finite witness has parent and
child ranks 1,904 each, with zero projector/CAR residual and rank defect zero.
The current matching to T has residual `1.78e-15` (tolerance `3e-11`).
The [duplicated-trace FAIL](nsc-hamiltonian-trace-representation.md) is preserved.
The ordinary sharp-block Z2 route remains OPEN because projecting nonzero
transmitting trace data leaves the Dirac operator domain. No delta coefficient is promoted
to a Hamiltonian link or source.

**Preceding definition — transmitting Dirac seed-seam domain: PASS in its finite scope.**
The [new domain](nsc-transmitting-dirac-domain.md) identifies the 33 retained
channels with two oriented traces of one field on the stored spacelike
surface. All 1,904 blocks are mapped using the actual coframe, quadrature
and spin-current basis. Normalization residual is `6.66e-16`, covariance
recovery `5.55e-16`, and the oriented boundary-form residual `1.08e-16`
(tolerance `3e-11`). The same-surface domain is now explicit; an instantaneous
Hamiltonian link, selected Cauchy history and remaining boundary action are
still separate definitions. No physical propagator or stress is assigned.

**Preceding calculation — known bulk endpoint jets.** The
[Hamiltonian-to-CTP chain](nsc-dirac-endpoint-jets.md) is now executable for
the already owned bulk vertices, using a fixed representative block in each
of the 33 channels on the frozen control. The derivative residual is
`1.0276e-10`; its contraction against the existing CTP action differs by
`4.9394e-11` (finite-difference tolerance `3e-8`). These are diagnostic bulk
jets, not complete physical transmitting jets. The
[J→M gate](nsc-physical-jet-extended-gate.md) remains **OPEN** for the link,
embedding and remaining boundary derivatives.

**Preceding calculation — transmitting CTP first differential.** The
[new action differential](nsc-transmitting-ctp-variation.md) agrees with the
existing determinant in all 33 retained channels, with maximum error
`5.551121183780314e-12` (finite-difference tolerance `3e-8`). Common-branch
variation is exactly zero; relative variation supplies the conditional
Gaussian source differential. Its [endpoint chain rule](nsc-transmitting-metric-pullback.md)
preserves the old covector and leaves the missing physical branch jets and
boundary remainder explicit. The [F→H decision](nsc-extended-action-completion-gate.md)
is **OPEN**; the old certificates are unchanged.

**Preceding binding result — full extended gate: OPEN.** The
[transmitting-action audit](nsc-full-extended-history-gate.md) finds
norm-preserving directions that change the stored covariance in all 33 channels,
with maximum tangent norm `0.679080121964124` and isometry residual zero
(tolerance `3e-11`). The [Weyl binding](nsc-weyl-endpoint-match.md) extracts
all eight diagnostic end-node coefficients exactly. The physical boundary
kernel selector, endpoint derivative and its pullback remain unevaluated.
The value `93.54264532195464` is the local diagnostic gradient, not an
evaluated two-sided physical mismatch.

The project starts from one simple claim: a wave, a particle, an unseen
boundary response, and a changing geometry can be different measurements of
one inherited spectrum. This page shows the equations and numbers behind that
claim.

| Plain-language idea | Calculated connection |
|---|---|
| A particle is a stable note of a field | One scalar link gives the exact massive Dirac dispersion and propagator self-energy |
| Another room can be unseen but still act here | Direct propagation and Schur elimination give the same curved boundary response |
| A boundary stores and transfers energy | The real and imaginary parts of one coherence give link energy and occupation transfer |
| Geometry can turn into particles | Work from a changing radius equals the energy of the created Dirac pairs |
| An inherited state can push geometry apart | Both computed radial null contractions are negative at the smooth neck |
| The completed stress has a child-time energy ledger | Proper volume and the child clock fix its first conserved density derivative |
| Backreaction begins with a constraint gate | The locked tensor and neck fail the Hamiltonian and homogeneous momentum constraints at $T=0$ |
| The tensor can select valid replacement geometry | Its Landau frame and density determine a constraint-complete Kantowski--Sachs neck |
| Future pressure needs the full Gaussian state | Four integrated stress moments do not determine the mode commutators under the changed metric |
| The seed Gaussian state is now inspectable | All 1,904 retained covariance blocks are serialized and reconstruct the completed tensor |
| Cauchy propagation is history dependent | Conditional maps are unitary, while equal endpoint data produce different propagated states |
| The present joint BVP is non-unique | The seed cannot start and the endpoint problem lacks an executable same-action history selector |
| The homogeneous history class is decided | The complete shift ledger excludes every smooth frequency-diagonal no-interface history from the locked seed |
| The larger tilted class now has real component owners | Reference, local bulk and weighted interface gates run separately before stationarity composition |
| Empty uniform space defines the zero | The recursive $a_0$ projector fixes $V_{\rm full}=0$ and preserves all gradient terms |
| Charged radius and spectral coefficients meet at one scale | The first cutoff-resolved checked branch gives $\Omega=3.9730743688$ and $\zeta=15.7853199397$ |
| The direct MMP-LLL substitution is decided | Its recursive power is fixed, but both neck null components have the wrong sign |

The construction connects the locally resolved field, the response across its
boundary, and the stress acting on its geometry through the same operator.
General research status is stated in the [README](../README.md#research-status).

The [September 10 development snapshot](development-update-2026-09-10.md),
through laboratory commit `3ba9e01`, now contains forty-six post-preprint records.
Its latest gates bind the charged field, parity, AP state and normalization to
the imported MMP throat; implement canonical CTP as the causal owner; and
apply the recursive zero-tadpole law that fixes the homogeneous unlinked
$V_{\rm full}=0$. The charged radius equation now supplies a fixed-$q$ scale
candidate. Applying its simplest inherited transparent LLL state to the actual
black-universe neck gives positive rather than required negative null sources,
so the missing owner is the full charged angular/compact CTP tensor on that
domain rather than another scale scan.

The subsequent [charged CTP neck calculation](nsc-charged-ctp-neck-source.md)
returns the retained child-frame tensor
$(\rho,T_{01},p_\parallel,p_\perp)=(0.0468138,0.00122387,-0.106079,0.0302511)$
and negative radial null components $-0.0568175$ and $-0.0617129$ at the
locked scale. This passes the retained sign test, while the hard full gate
remains open for the first positive compact levels' nonlocal CTP covariance
and four ADM variations; no cosmological projection is accepted yet.

That channel is now supplied by the
[positive compact CTP completion](nsc-compact-ctp-completion.md). It reads the
previous angular tensor without rerunning it and adds the two positive compact
levels below the locked cutoff with massive exterior scattering and
fourth-order superadiabatic subtraction. The completed free-Gaussian tensor is

$$
(\rho,T_{01},p_\parallel,p_\perp)
=(0.10074836289,0.00122387016,-0.34978809527,0.12958962232),
$$

with $T_{++}=-0.24659199207$ and $T_{--}=-0.25148747269$. The hard neck sign
gate therefore passes in this realization. Its homogeneous neck projection is
recorded.

The [proper-volume/clock projection](nsc-background-projection.md) now places
that tensor in the stored child geometry. At the neck,
$V=24.21233094747$, $\dot V=12\pi$, and conservation gives

$$
\boxed{\dot\rho_0=-\frac{\dot V_0}{V_0}\rho_0
-H_\parallel p_\parallel=0.38776013531.}
$$

The volume term is $-0.15686733379$ and pressure work is
$+0.54462746910$. The recorded parent Killing power satisfies
$P_{\rm parent}=-a_\parallel VT_{01}$ with residual
$4.4\times10^{-16}$: inside the child it is the opposite of a conserved
spatial-momentum charge, not an additional energy-deposition rate. The free
compact blocks give regular bulk $Q=0$. This fixes a local Cauchy jet; the
pressure derivatives and global background require backreaction and evolution
of the same CTP state before any $H(z)$ calculation.

The [child backreaction gate](nsc-child-metric-backreaction.md) applies the
existing ADM constraints before starting that evolution. With the locked
$A=0.04501936182826115$, the stored neck requires
$\rho=2A=0.0900387236565223$ and $T_{01}=0$. The completed tensor instead gives

$$
\boxed{\mathcal C_H=-0.118944813912720,
\qquad \mathcal C_M=-0.0135927088490713.}
$$

Both null components are still negative. The evolution stops at $T=0$
because stepping constraint-violating data would not solve the coupled
Einstein--CTP system. The required same-action completion is
$\Delta\rho=-0.0107096392302639$ and
$\Delta T_{01}=-0.00122387015580509$, or a source-derived initial geometry.
Neither repair is inserted in this result.

The [constraint-complete assembly](nsc-constraint-complete-neck.md) then takes
the permitted geometry-from-source route. The tensor's unique subluminal
Landau boost,

$$
v=0.00491447570673408,
$$

sets the normal-frame momentum to zero without changing occupations. At a
temporal minimum, the Hamiltonian constraint fixes

$$
\boxed{r_\star=\sqrt{\frac{2A}{\rho_L}}
=0.945328394434129.}
$$

The maximum constraint and frame-invariance residual is
$2.2\times10^{-16}$, with
$T_{++,L}=T_{--,L}=-0.249027703022043$. This PASS geometry remains in the
existing Kantowski--Sachs class but replaces the stored unit-radius Bronnikov
profile. Joint evolution of the same covariance and this selected metric is
therefore the next executable owner.

The [coupled CTP/metric evolution gate](nsc-coupled-ctp-metric-evolution.md)
checks the first required stress re-evaluation. Because

$$
\frac{1}{r_\star}-1=0.0578334533139646,
\qquad
\left.\frac{\partial\rho}{\partial\log r}\right|_C
=-0.460687999769215,
$$

the unit-radius tensor cannot be reused unchanged on the selected geometry.
The published source records retain integrated stresses and covariance bounds,
but not the complex covariance for each compact/angular/frequency mode.
Consequently the exact break occurs at $T=0$, during same-covariance source
evaluation before the first metric step. The fixed-tensor constraints and null
signs still pass; future pressure becomes determined only after the recorded
`ModeResolvedCauchyState` contract is supplied.

The [mode-resolved state artifact](nsc-mode-resolved-cauchy-state.md) now
supplies the seed side of that contract. A deterministic 407,306-byte NPZ
contains 33 channels and 1,904 physical covariance blocks. Reload gives

$$
\lambda_{\min}(C)=-2.67\times10^{-16},
\qquad
\lambda_{\max}(C)=1+4.15\times10^{-14},
$$

and reconstructs all completed unit-radius tensor components with maximum
residual $9.11\times10^{-13}$. The physical Landau Cauchy isometry and the
general-KS fourth-order/local source history are not determined by endpoint
$(v,\eta,r_\star)$, so finite selected-surface stress and metric evolution
remain gated on those two explicit owners.

The [Landau Cauchy-isometry owner](nsc-landau-cauchy-isometry.md) now computes
the blockwise time-ordered Dirac exponential for a supplied ADM history. Two
smooth histories with identical seed and Landau endpoint data give

$$
\max\lvert U_{0.5}-U_{1.0}\rvert=1.99966699255,
$$

while each satisfies $\lVert U^\dagger U-I\rVert<1.2\times10^{-14}$ and
preserves the CAR spectrum. This makes the selection issue finite and exact:
the endpoint tuple does not determine the physical $U_{L0}$. Neither control
is promoted to a metric history, so a unique finite $r_\star$ stress and its
updated constraints remain undefined.

The [joint history/state BVP gate](nsc-joint-history-state-bvp.md) now decides
the two allowed solve routes. Dynamic route A stops on

$$
(\mathcal C_H,\mathcal C_M)
=(-0.118944813913,-0.0135927088491)
$$

at the seed. Two-point route B has at least two endpoint-compatible solutions,
with

$$
\max|C_A-C_B|=0.192089815976105.
$$

The chosen result is route C: non-uniqueness bound and stop. The missing
selector is the executable general-KS history variation of the existing CTP
action, including its fourth-order, local and boundary terms. No finite stress
or metric history is assigned until that owner exists.

The [general-KS same-action owner](nsc-general-ks-same-action-history.md) now
executes the independent shift variation in the minimal homogeneous class.
Every already owned local, fourth-order-reference and magnetic contribution
has zero homogeneous $T_{01}$, whereas the serialized state yields

$$
\boxed{\mathcal E_\beta(0)
=-\frac{T_{01}}{2A}
=-0.0135927088490713.}
$$

Smooth unitary evolution makes this nonzero source continuous, while the
homogeneous momentum constraint requires zero on every slice. Hence no smooth
stationary history exists in the frequency-diagonal no-interface KS
truncation. A tilted transmitting interface or spatially inhomogeneous
spherical history lies outside this certificate.

The extended class now has separate executable
[fourth-order reference](nsc-general-ks-reference.md),
[node-wise local history](nsc-general-ks-local-history.md), and
[tilted interface](nsc-tilted-landau-interface.md) owners. Their composition
gives

$$
R_{\rm ref}\le 7.06\times10^{-12},\qquad
R_{\rm local,node}=2.08\times10^{-7},\qquad
R_{\rm interface}\le2.35\times10^{-14}.
$$

The [extended gate](nsc-extended-tilted-history-gate.md) nevertheless remains
open: the diagnostic Weyl end-node gradient is $93.54264532195464$, and
the weighted-isometry equation leaves $447488$ real kernel directions.
The new [endpoint binding](nsc-weyl-endpoint-match.md) preserves that gradient
and leaves the physical mismatch unevaluated, with target tolerance
$3\times10^{-11}$. The same transmitting boundary action
must select one $V_c$ and complete that endpoint variation before
$\delta\Gamma$ is single-valued. No optimizer, finite stress or existence
claim is launched before then.

## 1. Shared mass and visible response

In the invariant scalar-link sector,

$$
[H_8,\Pi_-]=0,\qquad
E^2=\lvert\mathbf p\rvert^2+\Phi^2,\qquad
\Sigma_p(E)=\Phi^2(E+\boldsymbol\alpha\cdot\mathbf p)^{-1}.
$$

The same $\Phi$ determines the mass shell and visible self-energy.
[The Dirac reduction](nsc-observable-bridge.md) specifies the projector, spin
frame and charge conventions. The complementary sector is spectrally
identical. The evaluated massless throat kernel has vector and axial
Clifford content; the scalar coupling and sector selection are action/domain
requirements.

## 2. Evaluated curved boundary maps

Parent and child maps are obtained from the smooth Dirac geometry with
oriented normal data. Direct inversion, Schur elimination and transfer agree
below $2\times10^{-15}$ relative error in the recorded finite system.
Independent radial integration differs by 0.00947%–0.03645% at the stated
resolution. The [full spinor calculation](nsc-chiral-boundary.md) retains both
angular sectors and the common spin frame.

At common dimensional energy, the recursive response is

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

The explicit $\Omega^{-1}$ is the first-order normalization factor.
[Finite recursion](nsc-regulated-recursion.md) and
[endpoint criteria](nsc-tail-limit.md) supply the recorded controls.

## 3. The state completes the energy account

For the finite quadratic system,

$$
z_B=\mathrm{Tr}_p(B C_{cp}),\qquad
E_{\rm link}=2\Re z_B,\qquad \dot N_p=2\Im z_B,
$$

$$
\dot E_{\rm tot}=\mathrm{Tr}(C\dot H).
$$

[The common-source derivation](nsc-common-source-derivation.md) includes
energy stored in the link. [The boundary-state implementation](nsc-boundary-state.md)
retains occupied and empty kernels and mixed initial correlations.
Its finite covariance reconstruction agrees with full evolution to
$1.1\times10^{-15}$.

A prescribed changing radius creates pairs with excitation energy equal to
the supplied work. In the representative $\kappa=1$ channel,
$\epsilon=0.02$, $\tau=0.5$ gives approximately
$1.40181\times10^{-5}$ pairs and energy $3.46640\times10^{-5}$
in throat units. [The vacuum-work record](nsc-vacuum-work.md) and
[normalized influence functional](nsc-influence.md) connect this response
to its geometric fluctuation weight.

## 4. Parent data produce a definite child stress

Affine-horizon data and the actual complex exterior reflection determine
the canonical massless state on the expanding child. In the child
orthonormal frame, with $\hbar=c=L_{\rm throat}=1$, $\mu=1$:

| Quantity | Recorded neck value |
|---|---:|
| Density $\rho$ | $-0.0033068962$ |
| Longitudinal pressure $p_\parallel$ | $-0.0494945966$ |
| Angular pressure $p_\perp$ | $-0.0429957713$ |
| Radial null $+$ | $-0.0528075900$ |
| Radial null $-$ | $-0.0527953956$ |
| Outward parent Killing power | $0.0001422206795$ |

Stress has dimensions $L^{-4}$, and power $L^{-2}$. Angular refinement
from 16 to 32 after the recorded tail correction changes tensor components
by less than $4\times10^{-6}$ across the sampled radii. This measures
numerical sensitivity under that asymptotic-tail approximation.
[Source and state](nsc-unruh-state.md).

The positive-Einstein static neck requires

$$
(\rho,p_\parallel,p_\perp,T_{\hat T\hat z})
=(2a,2a(1-3\pi),2a(1-3\pi),0),\quad a>0.
$$

The computed null signs are negative; density, pressure anisotropy and flux
leave explicit componentwise residuals. The [ADM frame map](nsc-adm-neck-source-map.md)
places those quantities in the metric equations. Under the stated asymptotic
Einstein normalization, the outgoing channel contributes
$dM_B/du=-P_0$, selecting a backreaction calculation.
[Flux and mass balance](nsc-parent-backreaction-gate.md).

## 5. One normalization and one metric variation

The [canonical–spectral bridge](nsc-canonical-spectral-bridge.md) partitions
the compact determinant into its physical tower, cutoff conversion and warp
contribution. The [Gaussian measure](nsc-gaussian-cutoff-measure.md) gives an
exact finite representation of the adopted proper-time action.
[Compact anomaly matching](nsc-compact-anomaly-bridge.md) fixes the domain and
the scope of the imported anomaly result.

The [ADM source equations](nsc-adm-source-constraints.md) retain lapse, shift,
radial metric and sphere radius through variation. Their Ward identities
bind stress, energy, flux and geometric work. The
[local warp source](nsc-warp-local-neck-source.md) evaluates an identified
summand on the actual neck; its nonlocal remainder has no small-error bound.

## 6. The causal and vacuum source now have explicit owners

Canonical Lorentzian evolution supplies the causal state-dependent response:

$$
\Gamma_{\rm one}^{\rm CTP}
=S_{\rm induced}[+]-S_{\rm induced}[-]
-i\log\det(I-C_0+C_0U_-^\dagger U_+).
$$

The [causal common functional](nsc-causal-common-functional.md) returns the
four metric forces, stress projections, boundary power, retarded response,
noise, and discrete Ward residuals from one covariance. Equal histories,
unitarity, covariance trace, and energy/work balance agree at approximately
$10^{-16}$ in its finite control.

The homogeneous gravitational zero is fixed separately by the same-action
relational law

$$
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},
\qquad
(V,A,C)\mapsto(0,A,C).
$$

This gives $V_{\rm full}=\lambda_4=\Xi=0$ while preserving Einstein, gauge,
curvature, Casimir, and link response. The
[normalization derivation](nsc-relational-vacuum-normalization.md) proves
idempotence, additivity, CTP normalization, and common room-scale covariance.

The physical geometry is the common stationary point:

$$
\left.\frac{\delta\Gamma_{\rm one}^{\rm CTP}}
{\delta g_\Delta^{\mu\nu}}\right|_{g_\Delta=0}=0.
$$

The same solution determines the following physical links:

| Target | Quantity to extract from the common solution |
|---|---|
| Constants and resolution | Characteristic cones, gravitational coupling, physical thermal state |
| Finite continuation | Evolved geometry, constraints, stress and boundary flux |
| Dark response | Background pressure and density, anisotropic stress, growth and lensing |
| Particles and measurement | Poles, residues, interactions, detector statistics |
| Mass and antimatter | Chiral interaction, sector selection, gauge and conjugate representations |
| Recursive infinity | State marginals, clock transfer and limiting observables |

## Inspect the record

The [preprint manifest](../results/manifest.json) contains 100 records;
[the later snapshot](../results/development-snapshot.json) contains forty-six.
[The development index](development-update-2026-09-10.md) links their
derivation notes. The canonical paper retains its technical appendices,
attribution and historical scale corrections.

The default demonstration displays stored evidence. Use
[the reproduction guide](reproducing.md) for an explicitly requested replay.
