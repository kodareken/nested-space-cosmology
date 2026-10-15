# Nested-Space Cosmology

## The Geometry of a Universe

**Douglas Ek**  
Working research manuscript, 15 October 2026

## Abstract

Nested-Space Cosmology asks how a region's geometry, fields and surrounding
structure determine the world available to an observer inside it. We develop
this question as a finite atlas: a parent with several black-hole boundaries,
possible unequal children and further generations. Local spherical geometry
provides tractable boundaries and observer charts; it does not make a
multi-hole parent globally spherical or establish the topology of our universe.
Published spacetime solutions, junction laws, branching mathematics and
regional reduction supply the foundations. The work then concentrates on
missing connections between them.

One connection is explicit. A fixed de Sitter/Schwarzschild wall joins a
regular expanding cap to a vacuum bridge and an asymptotic end on a checked
spacelike slice. A vacuum-end gluing theorem supplies a conditional
multi-body initial-data family preserving that body. In the isolated spherical
reference, the center, moving wall and right exterior have different causal
access: selected early signals connect them, while new emissions after the
chosen initial slice do not. A prescribed early conformal scalar pulse crosses
the wall and produces a calculable transient at the center. Its waveform alone
does not establish a unique exterior: the incident state is an independent
input. We now identify a finite characteristic preparation in the retained
white-hole sector and a bounded exterior-to-wall transfer operator. Its
existence is certified analytically. Evaluating it predicts a central pulse
followed by a small positive scattering tail within the declared finite
window. An independent positive-kernel argument fixes the tail's sign;
spatial and clock refinements test its numerical magnitude.
For this fixed geometry, a complete compatible signed central history also
determines a unique preparation. Positivity and a prescribed source cutoff
exclude some histories; geometry comparison requires those input restrictions
to be justified independently.
At quadratic order the canonically normalized probe adds no singular surface stress,
tangential flux jump or normal-stress jump at the fixed wall.
With fixed corner geometry and null normalization,
the same preparation requires a calculable correction to the initial
areal radius and its affine slope. Transverse constraint transport then
determines a geometric mass response on that initial null surface, including
a negative quadratic coefficient after the compact input. The existing wall
field traces also close an analytic exterior-side mass functional through
both null fluxes and the full conformal transformation. Keeping the cap vacuum at the actual radius
now supplies a compatible cap mass functional and a regular, quiet-anchored
wall displacement equation. Two-sided mass transport also preserves the
temporal and normal junction conditions. The same prepared field now drives
a computed quadratic exterior radius and lapse response on a finite
characteristic rectangle, checked against independent mass transport and
null constraints. Extending the same field evaluator to the retained wedge
now supplies a second mass path, the computed wall displacement and its
exterior null-clock correction, with independent numerical refinements.
The cap metric, central clock, next scalar response and finite-amplitude
remainder remain open.
A separate finite spherical field model demonstrates regional transfer,
geometric feedback and reduced memory, without identifying its mode-space
complement with a parent universe.

These results connect known mathematics in stated domains. They establish
neither regular black-hole formation of children nor a complete cosmology.
Physical inheritance of constants, quantum-state selection, persistent
gravitational effects, antimatter asymmetry and an arrow of time remain open.
The atlas makes those gaps explicit and identifies calculations that can
resolve them without requiring every microscopic question to be solved first.

## I. The question and the spherical idea

### 1. A universe viewed from within

An observer measures clocks, fields, curvature and correlations within a
limited domain. Such measurements need not uniquely identify the larger
structure containing that domain. Our question is whether this ambiguity
can become a constructive research programme: assemble possible surroundings,
determine which attachments satisfy their equations, and ask what each permits
an interior observer to experience.

The organizing sequence is **borrow established mathematics, connect compatible
regions, investigate branching and inheritance, calculate local consequences,
then compare with observations**. A universal action tying together all fields
and regions remains an ambitious candidate. It is not a prerequisite for
learning whether a particular boundary can be joined or whether a signal can
reach a particular observer.

A mature atlas would let a reader follow a possible universe from its
surroundings, through its boundaries, into an observer's world and onward
through compatible transitions. It would compare families of histories with
explicit geometry, dynamics and observables. Completion means dependable
coverage of declared model families, including exclusions and observational
degeneracies. Establishing that a family describes our universe is a further
empirical achievement. The present finite connections are steps toward that
destination, not a claim to have reached it.

The first motif contains two spatially separated black holes in a common
parent, possible unequal children, and one further generation. A separate
motif lets one hole have multiple offspring. These pictures suggest questions;
their arrows do not supply an equation of motion, a formation mechanism or
a communication channel. Figure 1 distinguishes the broad motif from the
specific attachment whose mathematics has been checked.

![Nested regions and the checked geometric connection](figures/geometry-overview.png)

*Figure 1. Possible descent (upper panel) remains hypothetical; the lower
panel shows the checked cap–wall–bridge body and observer/field reference.
Colors distinguish attachment, descent, influence and correlation. Circles
and spacing do not establish an embedding, relative scales or global
spherical topology. The [figure record](figures/geometry-overview.json)
identifies each connection's evidence and status.*

### 2. What spherical geometry means here

There is a useful borrowed energy argument with a precise domain: in Euclidean
space, a ball minimizes boundary area at fixed enclosed volume. A positive
constant isotropic surface tension therefore favors a spherical boundary when
surface energy is the quantity being minimized. Direction-dependent tension
instead selects a Wulff shape; see [Neumayer, 1503.06705v2, Section 1.1](https://arxiv.org/html/1503.06705v2).
This motivates testing rotational symmetry. It does not prove a cosmological
equilibrium, dynamical relaxation or a probability ranking over all geometries.
Gravity, rotation, matter stresses and the actual boundary law must be specified
before importing that variational conclusion into a spacetime model.

A spherical boundary is a two-dimensional surface with area $4\pi R^2$.
A spatial interior may be a three-ball; a closed spatial universe may instead
have three-sphere topology. These are different assertions. Drawing one circle
inside another establishes none of them.

In a locally spherically symmetric spacetime, a useful form is

$$
ds^2=g_{ab}(x)\,dx^a dx^b+R(x)^2d\Omega^2.
$$

Here $g_{ab}$ is a Lorentzian metric on the two-dimensional orbit space and
$R$ is areal radius. We use signature $(-+++)$ for the atlas. Radial
proper distance, areal radius and a diagram's vertical coordinate need not
coincide. In particular, a radial profile is not automatically a surface
that can be embedded in Euclidean three-space.

The physical spherical hypothesis is that useful cosmological regions or
interfaces may admit this structure, exactly or as a controlled approximation.
It is not a theorem inferred from statistical isotropy of observations.
Several off-center holes generally break global spherical symmetry.
We therefore use local spherical charts where justified and an explicitly
nonspherical parent where required. Deviations from symmetry are a question
for each model, rather than a reason to discard the organizing idea.

### 3. One thread through three bodies of work

The continuing question is: **how do a region's fields, boundary and
surroundings shape its state, geometry and measurements?**

| Body of work | What it contributes | Connection still missing |
| --- | --- | --- |
| [Preserved foundation](archive/foundation-v0.26.0/paper/nested-space-cosmology.md) | Spectral inheritance proposal, boundary operators, source and response constructions | A physical inheritance law and common cosmological source completion |
| [Finite spherical companion](finite-regeneration/main.tex) | A finite realization of regional transfer, geometric feedback and memory with declared action and state | Identification with an actual assembled parent–child spacetime |
| [Atlas evidence](../docs/nsc-cosmological-atlas.md) and this living manuscript | Geometric compatibility, a sourced slice, observer access, a calibration response and an evaluated finite exterior preparation | Remote preparation, later evolution, other attachments and observational discrimination |

Their mathematics is complementary. A finite operator's eliminated modes
are not automatically the exterior spacetime beyond a black hole. The atlas
must first identify a field equation, state, clock and matching law before
an operator reduction can acquire that geometric meaning.

The foundation's useful historical calculations remain in their source-bound
edition and owning notes. This manuscript changes the organizing story:
geometry and observer access lead; stronger unifying claims follow only when
their missing connections are supplied. No priority claim is made for junction
conditions, Schur complements, branching processes or the published wall
solutions used below.

## II. Regions, inheritance and the meaning of a connection

### 4. Four arrows that must not be confused

Each region in the atlas carries a geometry, source or state, clock
convention, dimensionless parameters and validity domain. Each connection
states the borrowed mathematics and the additional condition still needed.

| Connection | Mathematical question | What it does not establish by itself |
| --- | --- | --- |
| Descent | Where do a child's initial data originate? | A traversable path |
| Geometric attachment | Which boundaries share admissible metric and source data? | Formation or stability |
| Causal influence | Can an intervention reach an observer within the stated domain? | A nonzero field response without state and coupling |
| Correlation | Which information is shared through state or ancestry? | Communication between siblings |

Two regions may inherit correlated data without ever exchanging a later
signal. An attached interior can lie beyond the causal reach of a specified
parent observer. Conversely, communicating mouths require their own causal
and material assumptions; they cannot be inserted into a disconnected
offspring model as an unlabeled extra arrow.

The multiple-daughter geometry of [Barrabès and Frolov
(hep-th/9511136v1)](https://arxiv.org/html/hep-th/9511136v1), the charged
semiclassical multi-mouth route of [Emparan et al.
(2012.07821v2)](https://arxiv.org/html/2012.07821v2), and quantum-state
coupling distributions discussed by [Giddings and Turiaci
(2004.02900v2)](https://arxiv.org/html/2004.02900v2) remain separate library
entries. A state's coupling distribution is not itself a signal exchanged
between its branches. Their source cards and unresolved connections are
recorded in the [mathematical library](../docs/arxiv-studies.md).

### 5. A precise finite inheritance construction

<!-- nsc-claim:finite-normalized-inheritance -->

Evidence: [exact finite record](../results/development/nsc-nested-qualities-v1.json),
[frozen snapshot](nested-quality-evidence/snapshot.json), and
[proof and domain](../docs/nsc-nested-qualities.md).

The foundation proposes an inherited normalized law,
$\mathcal T_\Theta^*\mathbb D_\Theta=\mathbb D_\Theta$.
One controlled meaning is available without postulating an infinite cosmos.
Fix a Hermitian $d$-dimensional matrix $H$, a nonzero link $B$, a
positive scale ratio $\Omega$, and a finite depth $N$. For regional labels
$j=0,\ldots,N-1$, define a block Hermitian matrix by the following entries
where their indices exist:

$$
(J_{m,N})_{jj}=\Omega^{m+j}H,\qquad
(J_{m,N})_{j,j+1}=\Omega^{m+j}B.
$$

The lower off-diagonal blocks are the adjoints. Let $U_m$ consistently
relabel the equal-sized regional frames. Then

$$
J_{m+1,N}=\Omega U_mJ_{m,N}U_m^\dagger,
$$

$$
\Omega^{-1}U_m^\dagger J_{m+1,N}U_m=J_{m,N}.
$$

Each block acquires the same extra factor, so the identity follows directly.
For fixed depth,

$$
\mathrm{spec}(\Omega^{-m}J_{m,N})=\mathrm{spec}(J_{0,N}).
$$

This is a finite realization of normalized inheritance under chosen scaling
data. It does not derive $\Omega$, establish a physical transition through
a horizon, or show that the universe selects this law. Changing $N$
changes the boundary problem and can change the dressed spectrum.
The finite identity is not a proof of cosmic natural selection.

The shared law also does not fix identical states. A fermionic covariance
$0\leq C\leq I$ evolves as

$$
C(t)=e^{-iJt}C(0)e^{iJt}.
$$

Regional occupations and cross correlations can differ under the same $J$.
Thus “same laws,” “same state,” and “same observed constants” require separate
arguments. A physical transition would have to determine the relevant state
and dimensionless data, rather than obtaining them from the normalization
identity alone.

### 6. What an unresolved surrounding region can leave locally

Split a specified operator into retained and complementary blocks $A,D$
with link $B$. For an admissible complex spectral parameter $z$, elimination
gives the usual self-energy

$$
\Sigma(z)=B(zI-D)^{-1}B^\dagger.
$$

The retained resolvent is $[zI-A-\Sigma(z)]^{-1}$.
Nested elimination is associative when the required inverses exist.
For the finite inherited chain its normalized recurrence is

$$
\Gamma_j(x)=xI-H-\Omega^{-1}B
\Gamma_{j+1}(x/\Omega)^{-1}B^\dagger.
$$

The factor $1/\Omega$ is essential. This is borrowed finite linear
algebra applied to the declared family, not an additional force law.
To interpret $\Sigma$ as a cosmological exterior response requires a
physical operator, field state and boundary map. Sections 14–15 give a separate
finite Hilbert-space realization of this reduction, without identifying its
complementary modes with a cosmological exterior.

### 6.1. The local effective Lambda target

<!-- nsc-claim:target-local-effective-lambda -->

The intended cosmological target is a local effective $\Lambda$ arising
from inherited Schur memory: surrounding structure that a region does not
resolve contributes through $\Sigma(z)=B(zI-D)^{-1}B^\dagger$ to its
local response. A vacuum-form local gravitational part would be written
$T^{\mathrm{eff}}_{\mu\nu}=-(\Lambda_{\mathrm{eff}}/\kappa_4)g_{\mu\nu}$,
with $\kappa_4=8\pi G$. Identifying the physical exterior operator, its
state and boundary map, and obtaining that local part from the causal
response remain the missing connections. The finite Schur identity fixes
the reduction; it does not yet provide a measured cosmological value.
The [Schur-memory owner](../docs/schur-memory.md) keeps this target and
its retained spectral and geometric dictionary.

### 6.2. The retained local incoming gate

<!-- nsc-claim:local-source-insertion-majorant -->

Evidence: [majorant v7 record](../results/development/nsc-ks-source-operator-majorant-v7.json)
and [gate budget v5 record](../results/development/nsc-ks-gate-budget-v5.json).

The [local incoming source-insertion majorant v7](../docs/nsc-ks-source-operator-majorant-v7.md)
is a separate controlled application. For the fixed switched state law
$C_\Sigma[g]=U_gC_{\mathrm{up}}U_g^\dagger$ on $I=S(1)+[.12,.18]$,
family $(14,1)$ covers **1168/1168 signed rows**, with original weights
counted once and no remaining row keys. Its continuous insertion bounds are

$$
\epsilon_N\leq
3646674054457366029024452496769420151764188105426118876609\,2^{-229}.
$$

$$
\epsilon_\beta\leq
5793343250846686183515013928067141454707331929528967339969\,2^{-230}.
$$

These are approximately $4.226958072607968\times10^{-12}$ and
$3.357610065481873\times10^{-12}$ for this family's source-error
insertion. The physical local gate remains **OPEN**:
`physical_upstream_budget_component` is null. Source quadrature,
numerical field error, value-assembly arithmetic and every other source
family are excluded from these bounds. Complete row coverage in this
family does not supply those missing errors.

The [v5 gate budget](../docs/nsc-ks-gate-budget-v5.md) names nine distinct
components, with total error allocation $2\times10^{-11}$ and residual
reserve $10^{-11}$ against the physical tolerance $3\times10^{-11}$.
Its applicable component bounds remain null; 60 positive and 120 signed
families are declared requirements rather than authenticated coverage.
The gate is neither an existence nor a non-existence certificate and
does not determine $\Lambda_{\mathrm{eff}}$. Its paused campaign does not
invalidate the finite inheritance theorem or the scalar-chain results.
Retired supporting records are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

## III. A finite map of possible geometries

### 7. Start with two holes, without manufacturing two children

The explicit parent benchmark uses time-symmetric Brill–Lindquist vacuum
initial data,

$$
h_{ij}=\psi^4\delta_{ij},\qquad K_{ij}=0,\qquad
\psi=1+\frac{m_A}{2r_A}+\frac{m_B}{2r_B}.
$$

The recorded choice is $m_A=1$, $m_B=1.5$, separation $d=12$.
It has outer ADM mass $2.5$ and two inner-end masses $1.0625$ and
$1.5625$. These are three asymptotically Euclidean ends of an initial
slice. Extra ends do not by themselves become cosmological children.
The companion hole supplies nonspherical structure, so an exact isolated
Schwarzschild collar cannot simply be assumed at either puncture.
The [parent record](../results/development/nsc-cosmological-atlas-v2/parent.json)
and [owning analysis](../docs/nsc-cosmological-atlas.md) retain the
construction and checks.

![Two-hole parent slice and its distinct asymptotic ends](../results/development/nsc-cosmological-atlas-v2/parent-slice.png)

*Figure 2. The explicit two-hole initial-data benchmark. This is a
nonspherical common parent with unequal ends, not a solved evolution that
creates children. Local spherical models must earn an attachment to it or
belong to a separate parent family.*

The alternative parent route imports an existence result instead of
re-solving its gluing equations. [Chruściel, Corvino and Isenberg
(0909.1101v1)](https://arxiv.org/html/0909.1101v1) provide sufficiently
separated scalar-flat bodies with retained interior regions and a common
Schwarzschild outer end. Their more general
[vacuum-end theorem (1004.1355v1)](https://arxiv.org/html/1004.1355v1)
allows appropriate non-time-symmetric ends. These establish separate
families, not a modification that leaves the displayed Brill–Lindquist
data unchanged. No explicit sufficient separation is certified here.

### 8. A failed attachment can sharpen the map

<!-- nsc-claim:atlas-junction-screen -->

Evidence: [radiation junction](../results/development/nsc-cosmological-atlas-v2/junction.json),
[moving-boundary conservation/control](../results/development/nsc-cosmological-atlas-v3/moving-boundary.json),
and [source assumptions](../docs/nsc-atlas-moving-boundaries.md).

The first attachment screen pairs a closed radiation FLRW region with a
vacuum exterior and a pure constant-tension wall. For a comoving boundary,
$a=A\sin\eta$, $t=A(1-\cos\eta)$, and $R=a\sin\chi$.
The junction condition would require, over an interval of changing $R$,

$$
k^2R^4-2\cos\chi\,kR^3+2MR-A^2\sin^4\chi=0.
$$

For nonzero $A\sin\chi$, this cannot be a polynomial identity on that
interval. The recorded $\chi=\pi/3$, $M=0.4A$ test exposes a source
incompatibility, not a failure of all interiors.

Could a moving wall repair the same model? Intrinsic conservation answers
before a new simulation is needed. Fixed pure tension has zero intrinsic
stress divergence, so the Codazzi relation requires zero net normal energy
flux. In the fluid rest frame,

$$
T(U,n)=(\rho+p)\gamma^2v.
$$

A vacuum exterior contributes no flux and radiation has $\rho+p>0$;
therefore $v=0$. Motion alone cannot rescue these assumptions.
This is a compatibility application of an established conservation law.
Allowing a different surface source, exterior state or flux would define
a different candidate requiring its own test.

### 9. A distinct fixed-source wall that does connect

For the next branch, take de Sitter inside and Schwarzschild outside,

$$
f_-=1-H^2R^2,\qquad f_+=1-\frac{2M}{R}.
$$

In $G=c=1$ units let $\kappa=4\pi\sigma>0$, with surface pressure
$P=-\sigma$. The signed angular extrinsic curvatures obey

$$
\beta_--\beta_+=\kappa R,\qquad
\beta_\pm=\epsilon_\pm\sqrt{\dot R^2+f_\pm}.
$$

Here the dot is wall proper time and $\epsilon_\pm$ fixes the normal
orientation. Squaring and then forgetting those signs would lose physical
branches. The useful explicit expressions are

$$
\beta_-=\frac{M}{\kappa R^2}+\frac{\kappa^2-H^2}{2\kappa}R,
$$

$$
\beta_+=\frac{M}{\kappa R^2}-\frac{\kappa^2+H^2}{2\kappa}R,\qquad
\dot R^2=\beta_-^2-f_-.
$$

We hold $H\ell=1$, $\kappa\ell=2$, $M/\ell=0.4$ fixed.
On the tested expanding segment $1.2\leq R/\ell\leq1.8$,
$\beta_->0$ and $\beta_+<0$. The wall proper duration is
$0.3718706103\ell$. Direct unsquared and pressure checks, including a
separate Christoffel-based calculation, agree within the recorded numerical
residuals. The wall's areal-radius rate is not a local three-speed.
This published wall family is classified more fully in
Section 11; the source has not been retuned to make the segment work.

### 10. A sourced initial slice and a conditional common parent

<!-- nsc-claim:atlas-sourced-slice -->

Evidence: [slice constraints](../results/development/nsc-cosmological-atlas-v4/connection-slice.json),
[parent assembly](../results/development/nsc-cosmological-atlas-v4/parent-assembly.json),
[slice derivation](../docs/nsc-atlas-connection-slice.md), and
[end-only gluing argument](../docs/nsc-atlas-connection-gluing.md).

At $R_w=1.5\ell$, the interior is a flat expanding de Sitter cap: a
regular spatial three-ball with $h=dl^2+R^2d\Omega^2$, $R'=1$
and $K^i{}_j=-H\delta^i{}_j$. We use
$K=-\mathcal L_T h/2$, where $T$ is the future unit slice normal.
Define $V=T(R)$. The cap has $V_-=HR$.
Its slice determines the wall boost,

$$
\gamma=\frac{\beta_--HR\dot R}{1-H^2R^2},\qquad
\gamma v=\frac{\dot R-HR\beta_-}{1-H^2R^2}.
$$

The exterior inherits

$$
R'_+=\gamma\beta_+-\gamma v\dot R,\qquad
V_+=\gamma\dot R-\gamma v\beta_+.
$$

At this event $\gamma=1.00926203$, $v=0.13516599$,
$R'_+=-2.02778609$, and $V_+=1.90925372$.
The exterior slice is therefore not time symmetric at the material wall.

A fixed smooth slicing function tapers $V(R)$ to zero at $R_0=\ell$,
before the Schwarzschild throat $2M=0.8\ell$. In the transition,

$$
R'=-\sqrt{f_++V^2},\qquad
K^\theta{}_\theta=-\frac{V}{R},\qquad K^l{}_l=-\frac{dV}{dR}.
$$

These data describe a slice of the same vacuum Schwarzschild spacetime.
They do not interpolate the matter source. Its mass identity is

$$
M=\frac{R}{2}(1+V^2-R'^2).
$$

Where $V=K=0$, the slice continues smoothly through a bridge using
$R=2M\cosh^2u$, with the cap on the negative-$u$ side and the
selected parent end on the positive-$u$ side.

The distributional shell constraints are equally necessary. The surface
energy and radial momentum on this slice are
$E_\Sigma=\sigma\gamma$, $J_\Sigma=\sigma\gamma v$, giving

$$
[R']=-4\pi R E_\Sigma,\qquad [V]=4\pi R J_\Sigma.
$$

Their recorded residuals are below $2\times10^{-15}$ in normalized
units; independent finite-difference bulk residuals are below
$10^{-10}$. These are implementation checks, not a global stability
estimate.

![A regular cap, sourced wall, vacuum transition and bridge](../results/development/nsc-cosmological-atlas-v4/connection-slice.png)

*Figure 3. The verified initial connection. Spatial distances are not ages.
Because $|R'|>1$ on part of the chosen slice, the radius profile should
not be interpreted as a Euclidean surface of revolution. The protected
right-exterior observer interval is a separate causal-domain statement.*

The selected right end of this body $A$, of mass $0.4\ell$, can be
used with the selected right end of an ordinary complete Schwarzschild
bridge $B$, of mass $0.6\ell$, in the imported vacuum-end gluing
theorem. It acts on smooth vacuum ends and preserves an open overlap.
Reattaching unchanged interiors on identical full $(h,K)$ data leaves
the original shell delta source intact and creates no seam source.
We do not apply a smooth global-data corollary directly to a thin shell.

The resulting conditional initial-data family has one designated common
parent end **and the untouched other end of $B$**: two asymptotic ends
in total. That second end has not become a cosmological child.
Sufficiently large separation is required but no numerical threshold or
explicit multi-body development has been calculated. The achievement is a
compatible sourced initial connection, not a regular formation history,
a general Einstein–shell evolution theorem or a stable cosmological future.

## IV. What the observers can reach

### 11. Complete one reference before assigning its causal arrows

<!-- nsc-claim:atlas-causal-access -->

Evidence: [isolated causal completion](../results/development/nsc-cosmological-atlas-v5/causal-map.json),
[observer access](../results/development/nsc-cosmological-atlas-v5/observer-access.json),
and [protected parent interval](../results/development/nsc-cosmological-atlas-v4/causal-window.json).

The wall's isolated spherical completion is the flat-exterior limit of
solution 13 in [Aguirre and Johnson, gr-qc/0508093v2,
Figures 7 and 10](https://arxiv.org/html/gr-qc/0508093v2).
Writing $x=R/\ell$, its expanding branch has

$$
\dot R^2=\frac{25x^2}{16}-1+\frac{3}{10x}+\frac{1}{25x^4}>0.
$$

It starts at a finite proper-past $R=0$ boundary in the Schwarzschild
white-hole sector, crosses $R=2M$ into the left exterior, and crosses the
de Sitter horizon into an expanding interior chart. It reaches arbitrarily
large radius at infinite wall proper time. Its Schwarzschild curvature
$48M^2/R^6$ diverges at the past boundary. The regular cap on the
chosen slice therefore does not imply a regular past formation.
The retained reference also contains the usual future Schwarzschild
black-hole singular boundary.

Choose three clocks. $C$ is the cap center, with de Sitter flat proper
time $t_C$. $W$ carries wall proper time $\tau_W$, set to zero at
$R=1.5\ell$. $P$ is static in the right exterior at
$R_P=1.6\ell$, with proper time
$\tau_P=\sqrt{1-2M/R_P}\,t_{\mathrm{Schw}}$. Their origins coincide
with a selected slice convention; this is not an equality of universe ages.

For the flat de Sitter chart,

$$
\eta=-\ell e^{-t_C/\ell},\qquad R=a r,\qquad a=e^{t_C/\ell}.
$$

Radial null lines have constant $\eta-r$ or $\eta+r$.
At $W_0$, $\eta=-\ell$ and $r=1.5\ell$.
Its inward ray would reach $\eta=0.5\ell$, outside the chart's
future $\eta<0$; hence it never reaches $C$.
A central signal emitted at $t_C=0$ likewise cannot catch the future wall.

![Interior and exterior causal charts with named observers](../results/development/nsc-cosmological-atlas-v5/spacetime.png)

*Figure 4. The isolated reference in joined causal charts. Horizons are
clipped to retained domains; the central cosmological horizon continues
into the Schwarzschild chart along its appropriate constant null coordinate.
The faint past $R=0$ line also shows uncut chart context, part of which is
excised. Open wall endpoints at $R/\ell=0.02$ and $6$ mark finite
plot samples, not physical beginnings or endings.*

| Emission and destination | Result in this reference | Required distinction |
| --- | --- | --- |
| New $C\to W$ or $W\to C$ emission at or after the anchor slice | No future reception | Earlier signals may arrive later |
| Right parent $P\to W$ or $P\to C$ | No future causal route | Left and right exterior regions are different |
| Early $W\to C$, while $R<\ell$ | Allowed | Reception time depends on the emission event |
| Early $C\to W\to P$, crossing while $R<2M$ | Allowed | The limiting horizon ray has no finite parent reception |
| These same arrows in a glued multi-body future | Unresolved | Initial-data gluing does not supply that global future |

For a wall event with interior time $t_W$, the central emission reaching
it and the reception of its inward ray satisfy

$$
t_e=t_W-\ell\ln(1+R/\ell),\qquad
t_r=t_W-\ell\ln(1-R/\ell).
$$

The second relation requires $R<\ell$.
At the early event $R=0.5\ell$, the recorded central emission is
$-1.5242168241\ell$ and inward reception is
$-0.4256045355\ell$. An outward early route through that event reaches
$P$ at its own proper time $+2.5311656510\ell$.
Thus exclusion of **new emissions** after a slice must never be rewritten
as exclusion of all later receptions. Emissions approaching the central
horizon from below can arrive arbitrarily late.

The critical central emission for escape to the right parent is
$t_C=-1.2312944525\ell$. Equality follows the parent horizon and has
no finite parent reception. That horizon extends into the cap along
$\eta-r=\eta_W(2M)-r_W(2M)$; it is not the surface $R=2M$
everywhere. The central observer's de Sitter horizon $\eta+r=0$
is a different boundary. Detailed clock domains and finite checks are
retained in the [observer note](../docs/nsc-atlas-causal-observers.md).

### 12. Different exteriors can be locally indistinguishable

There is also a result that does not require the glued parent's full future.
Retain the exact right Schwarzschild $K=0$ slice for
$2M<R<R_{\mathrm{cut}}$, with $R_{\mathrm{cut}}=10M$.
For a static observer at $R_{\mathrm{obs}}=4M$, the exterior cutoff
cannot influence it within

$$
0<t_{\mathrm{Schw}}<
r_*(R_{\mathrm{cut}})-r_*(R_{\mathrm{obs}}),\qquad
r_*(R)=R+2M\ln(R/(2M)-1).
$$

The observer's proper duration is the right side multiplied by
$\sqrt{1-2M/R_{\mathrm{obs}}}$.
For $M=0.4\ell$, the bound is
$t_{\mathrm{Schw}}<3.5090354889\ell$, or
$\tau_P<2.4812627896\ell$.
The inner optical boundary is at negative infinity, so it introduces no
finite inner cutoff.

Shared smooth vacuum data give the same vacuum development in this domain.
The [owning argument](../docs/nsc-atlas-causal-observers.md) uses standard
Cauchy-development uniqueness, with [Sbierski
(1309.7591v2)](https://arxiv.org/html/1309.7591v2) as a source.
The bound ends a guaranteed interval; it does not predict that a signal
actually arrives then. Matching additional fields and their states requires
additional data.

This supplies one precise atlas of indistinguishable exteriors. Distinct
larger assemblies can share a protected local geometry. The next useful
question is which later measurement or extra state datum separates them.
Failure to distinguish exteriors can be a result, rather than a reason to
add unconstrained parameters until a preferred story appears.

## V. Fields, transfer and finite memory

### 13. A first field that reaches the interior

<!-- nsc-claim:atlas-field-response -->

Evidence: [field-transfer record](../results/development/nsc-cosmological-atlas-v6/field-transfer.json)
and [field-transfer owner](../docs/nsc-atlas-field-transfer.md).

On the fixed de Sitter/Schwarzschild wall, retain $H\ell=1$, $k\ell=2$,
$M/\ell=0.4$ and a real classical massless scalar in its spherical mode,
with conformal coupling $\xi=1/6$ and no additional scalar wall interaction.
For $\chi=R\phi$, the distributional scalar equation gives

$$
[\phi]=0,\qquad [n\phi]=6\xi k\phi,\qquad
[\chi]=0,\qquad [n\chi]=(6\xi-1)k\chi=0.
$$

The conformal interface barrier cancels; the smooth exterior potential
$f_+2M/R^3$ remains. A regular center has
$\chi=F(\eta+r)-F(\eta-r)$. For incoming wall characteristic
$I=(d/d\tau_W+n)\chi$, its exact conditional response is

$$
\phi_C[T(\tau_W)]=\frac{I(\tau_W)}{T'(\tau_W)},\qquad
T=t_W-\ell\log(1-R/\ell),\qquad
T'=\frac{\gamma(1+v_{\rm pec})}{1-R/\ell}>0.
$$

The prescribed wall-clock calibration has half-width $0.03\ell$, duration
stretch $2.63230219$ and normalized extrema $+23.84144800$ and
$-24.46533478$. Its central field vanishes after the received support
through $t_C=0$, conditional on zero incoming data afterward and the
proved no-return interval. These are field amplitudes, not energy fractions.
The calibration does not generate its exterior input.

### 13.1. Where the early signal can originate

<!-- nsc-claim:atlas-field-preparation -->

Evidence: [preparation record](../results/development/nsc-cosmological-atlas-v7/preparation-access.json)
and [preparation owner](../docs/nsc-atlas-field-preparation.md).

Independent data can instead be placed in retained white-hole sector IV.
The ordinary right-parent observer remains causally excluded. Let
$x=R/\ell$, $F(x)=(1-x/0.8)e^{x/0.8}=UV$, $V_a=V_W(0.45)$,
$V_b=V_W(0.55)$ and $U_0=F(0.4)/V_a$. Define $D$ by
$V_a\leq V\leq V_b$ and $U_0\leq U\leq U_W(V)$;
there $0.4\leq R/\ell\leq0.55$. Prescribe zero bottom and cap data and
$g(V)=\varepsilon B[4(V-V_m)/(V_b-V_a)]$ on $U=U_0$, where
$B(z)=(1-z^2)^4$ for $|z|<1$ and zero otherwise. This is a compact
$C^3$ input, independent of any desired central waveform.

With $L(u)=V_a$ below $U_a$ and the inverse wall graph above it,

$$
\chi=g-\int_{U_0}^{U}du\int_{L(u)}^Vdv\,K\chi,\qquad
K=\frac{16M^4}{R^4}e^{-R/(2M)},\qquad \chi_U|_W=0.
$$

The integral operator has analytic uniform-norm bound

$$
q_*=16e^{-1/2}[\sqrt{F(0.4)}-\sqrt{F(0.55)}]^2
=0.13882099234357653<1.
$$

It therefore determines a unique finite transfer. The central arrival
window is $-0.581758951595344\leq T/\ell\leq-0.250230700560700$;
the no-return null-label gap is $2.634591912177392\ell$. Remote physical
generation and data beyond that window remain unspecified.

### 13.2. A prepared pulse leaves a finite scattering tail

<!-- nsc-claim:atlas-prepared-response -->

Evidence: [prepared-response record](../results/development/nsc-cosmological-atlas-v8/prepared-response.json)
and [response owner](../docs/nsc-atlas-prepared-response.md).

For this unchanged input, let $A$ denote the positive integral operator
above. After the compact input ends, its Neumann series gives
$-Ag\leq\chi\leq-(1-q_*)Ag<0$. Hence the determined wall derivative

$$
J=g'(V)-\int_{U_0}^{U_W(V)}K\chi\,dU>0,
\qquad (1-q_*)\mathcal B(V)\leq J(V)\leq\mathcal B(V),
\quad \mathcal B=\int KAg\,dU,
$$

produces a positive central tail on
$-0.340772594886934<T/\ell\leq-0.250230700560700$.
The sampled normalized extrema are $47.29833694$ and $-45.69595941$;
the final field is $0.00129049806$. At the final point,
$J/\varepsilon=0.0008020390045$ lies in the independently evaluated
bracket $[0.0006965878198,0.0008088769160]$.
The 257-to-513 spatial and 129-to-257 geometry comparisons change the
central field by at most $2.79286258\times10^{-7}$ and
$3.40622250\times10^{-8}$. The calculation uses 5,597,264 kernel values
within its ten-million/60-CPU-second limits. These numerical controls
support the finite waveform; the sign follows from the analytic argument.
No permanent memory, continuum enclosure or later decay law is established.

### 13.3. What the signal determines when preparation is unknown

<!-- nsc-claim:atlas-observer-inverse -->

Evidence: [inverse record](../results/development/nsc-atlas-observer-inference-v1/inference.json)
and [inverse owner](../docs/nsc-atlas-observer-inference.md).

In this fixed geometry, put $f=\chi/\varepsilon$,
$h=f(U_W(V),V)$, $\theta=T/\ell$ and $Y=\ell\phi_C/\varepsilon$.
The zero cap fixes the integration constant:

$$
h(V)=\frac12\int_{\theta_a}^{\theta(V)}Y(s)ds,\qquad
f=(I-\mathsf{B})^{-1}Eh,\qquad g(V)=f(U_0,V),
$$
$$
(\mathsf{B}f)(U,V)=\int_U^{U_W(V)}dp\int_{L(p)}^Vdq\,K(p,q)f(p,q).
$$

The operator-norm bound $q_*<1$ gives a bijection between initially
quiet $C^2$ source traces and compatible $C^1$ signed central histories.
For nonzero nonnegative compact sources, the post-input restrictions are
$h<0$, $Y>0$ and $\int_{\theta_a}^{\theta(V)}Yds=2h<0$.
An arbitrary fitted trace need not satisfy them. This analytic source
inverse neither identifies geometry nor reconstructs a source numerically
from sampled or intensity-only data.

### 13.4. The transmitted probe and the wall's stress

<!-- nsc-claim:atlas-scalar-wall-consistency -->

Evidence: [stress record](../results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json)
and [wall-stress owner](../docs/nsc-atlas-scalar-wall-consistency.md).

Fix canonical kinetic normalization, $\kappa_4=8\pi G$ and
$f(\phi)=1/(2\kappa_4)-\phi^2/12>0$ in the weak-field interpretation.
The full stress includes the nonminimal Einstein term once:

$$
T^\phi_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
-\frac{1}{2}g_{\mu\nu}(\nabla\phi)^2
+\xi[G_{\mu\nu}\phi^2+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2].
$$

For the original sufficiently regular traces, $d=[n\phi]=k\phi$ and
$C_{ab}=[K_{ab}]-\gamma_{ab}[K]=2k\gamma_{ab}$ give
$S^\phi_{ab}=\xi(2\phi d\gamma_{ab}-\phi^2C_{ab})=0$ and
$[T^\phi_{na}]=(1-6\xi)k\phi D_a\phi=0$.
Thus no additional scalar surface source or tangential flux jump is
required at quadratic order. The individual fluxes and smooth bulk
gravitational source need not vanish.

### 13.5. The prepared pulse requires gravitational initial data

<!-- nsc-claim:atlas-characteristic-initialization -->

Evidence: [initialization record](../results/development/nsc-atlas-characteristic-initialization-v1/initialization.json)
and [initial-data owner](../docs/nsc-atlas-characteristic-initialization.md).

On the source leg, keep the reference lapse, radius $0.4\ell$, both corner
expansions and affine normalization. Use the affine background radius
$x=R_0(U_0,V)/\ell$, $V(x)=F(x)/U_0$ and
$b(x)=B(4[V(x)-V_m]/(V_b-V_a))$; primes here and in section 13.6 mean
$d/dx$. With $p=b(x)/x$,
$\phi=(\varepsilon/\ell)p+O(\varepsilon^3)$ and
$\alpha=\kappa_4\varepsilon^2/\ell^2$, write
$R/\ell=x+\alpha h+O(\alpha^2)$. The intrinsic constraint forces

$$
h''=-\frac{x}{3}p'^2+\frac{x}{6}pp'',\qquad h(0.4)=h'(0.4)=0.
$$

After the compact trace, $h'<0$ and $h<0$ by the integrated expressions;
$h'$ becomes constant. At $x_b=0.4815042082576364$,
$h=-6.12283354819419$ and $h'=-162.18055362491776$.
Direct integration, integration by parts and the transformed affine
Einstein-frame equation agree. These are initial radius coefficients,
not a mass shift or an evolved wall; no finite amplitude or perturbative
remainder bound is selected.

### 13.6. An invariant mass response on the same initial surface

<!-- nsc-claim:atlas-source-mass -->

Evidence: [source-mass record](../results/development/nsc-atlas-source-mass-v1/mass.json)
and [mass owner](../docs/nsc-atlas-source-mass.md).

With $m=M/\ell=0.4$, the fixed transverse corner datum and scalar equation supply
$d=(a/x-1)/4$, $a=0.8$, and $s'=-(s+dp')/x$, $s(0.4)=0$.
For the geometric mass length, use

$$
\mathfrak{m}=\frac{R}{2}[1-g^{ab}\partial_aR\partial_bR],\qquad
\frac{\mathfrak{m}_J}{\ell}=m+\alpha\mu_J+O(\alpha^2).
$$

The borrowed flux and full conformal frame conversion yield

$$
\mu_E=-\frac14\int_{0.4}^{x}t(a-t)p'(t)^2dt,\qquad
\mu_J=\mu_E+\frac{mp^2}{12}+\frac{x^2p}{3}(s+dp').
$$

After the input, the frame masses agree at the constant negative coefficient
$\mu_J=\mu_E=-28.88319705406449$. Their largest sampled difference
during input is $4.265459834842238$. Independent mixed-radius and improved
stress laws agree. This is a finite source-leg response relative to the
reference mass, not a negative total mass, vacuum future or asymptotic energy.

### 13.7. The scalar adds no normal-force jump

<!-- nsc-claim:atlas-scalar-normal-force -->

Evidence: [normal-force record](../results/development/nsc-atlas-scalar-normal-force-v1/normal-force.json)
and [normal-force owner](../docs/nsc-atlas-scalar-normal-force.md).

The Gauss jump $[G_{nn}]=-2k\bar K$ and complete scalar stress give

$$
[T^\phi_{nn}]=2\xi(6\xi-1)k\bar K\phi^2=0,\qquad
S^{ab}\bar K_{ab}=[T_{nn}].
$$

The cap contribution $G^-_{nn}=-3H^2$ is essential. Fixed vacuum and
material tension therefore require $\delta\bar K=0$ at order
$\varepsilon^2$. This mean-curvature variation includes the bulk metric,
induced metric and wall displacement. It does not require zero displacement
or an unchanged trajectory.

### 13.8. Exterior mass closes from the known wall field

<!-- nsc-claim:atlas-wall-mass-closure -->

Evidence: [analytic closure record](../results/development/nsc-atlas-wall-mass-closure-v1/wall-mass.json)
and [wall-mass owner](../docs/nsc-atlas-wall-mass-closure.md).

On the smooth exterior wall limit in the same $0.45\leq R/\ell\leq0.55$
window, put $s=\tau/\ell$, $x=R/\ell$, $v=dx/ds$, $\beta=\beta_+$,
$w=\chi_W/\varepsilon$, $j=dw/ds$ and $p=w/x$. The two known derivatives
are $A=(j-vp)/x$ and $B=(j-\beta p)/x$. Both null mass fluxes and the
complete conformal endpoint term give

$$
L=\frac{x^2}{2}\left[\beta AB-\frac{v}{2}(A^2+B^2)\right],\qquad
Q=\frac{mp^2}{12}+\frac{x^2p}{6}(vA-\beta B),\qquad
\mu_J(s)=\int_{s_a}^sL(s')ds'+Q(s),\quad Q(s_a)=0.
$$

The original-frame improved stress independently gives
$d\mu_J/ds=L+dQ/ds$. Neither expression divides by $\beta$.
This is an analytic conditional functional, not an evaluated mass curve
or a proof that a coupled continuation exists. The wall tail supplies no
late mass sign or plateau. Coupled bulk and wall evolution, the corrected
observer response, physical preparation selection, stability and a
finite-amplitude validity bound remain open. Detailed scalar-chain
derivations belong to their owners; historical supporting calculations
are cited at [Zenodo](https://doi.org/10.5281/zenodo.23090015).

### 13.9. Cap mass and a regular wall displacement functional

<!-- nsc-claim:atlas-cap-wall-response -->

Evidence: [cap-wall record](../results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json),
[model specification snapshot](../results/development/nsc-atlas-cap-wall-response-v1/model-specification.json)
and [cap-wall owner](../docs/nsc-atlas-cap-wall-response.md).
The common [model interface](../src/recursive_horizons/nsc_cosmological_atlas.py)
separates the fixed action, prescribed preparation and unit scale from
derived responses. Its parameters $h=H\ell=1$, $b=k\ell=2$ and
$m=M/\ell=2/5$ remain declared inputs. We retain the same canonical real
massless scalar, conformal coupling, scalar-independent vacuum terms and
constant material tension in the Jordan frame.

Let $r=x+\alpha\rho+O(\alpha^2)$ be the actual wall radius in units of
$\ell$, with $\alpha=\kappa_4\varepsilon^2/\ell^2$ and proper wall time
$s=\tau/\ell$. The two geometric masses have expansions

$$
\frac{\mathfrak{m}_+}{\ell}=m+\alpha\mu_++O(\alpha^2),\qquad
\frac{\mathfrak{m}_-}{\ell}=\frac{h^2r^3}{2}+\alpha\nu_-+O(\alpha^2).
$$

Here $\mu_+$ is the exterior coefficient $\mu_J$ of Section 13.8.
The cap vacuum is evaluated at the actual radius: its total quadratic
mass coefficient on the reference wall is
$\nu_-+3h^2x^2\rho/2$. Holding that vacuum mass constant would omit a
term required by the angular junction. The shared quiet corner fixes
$\mu_+(s_a)=\nu_-(s_a)=\rho(s_a)=0$.

Using the signed normals $\beta_-=\beta_++bx$, define
$B_-=(j-\beta_-p)/x$ with the same $p=w/x$, $j=w_s$ and $A=p_s$
as above. The cap's vacuum-subtracted mass closes as

$$
L_-=\frac{x^2}{2}\left[\beta_-AB_--\frac{v}{2}(A^2+B_-^2)\right]
+\frac{h^2x^2vp^2}{2},
$$
$$
Q_-=\frac{x^2p}{6}(vA-\beta_-B_-)-\frac{h^2x^3p^2}{12},\qquad
\nu_-(s)=\int_{s_a}^sL_-(s')\,ds'+Q_-(s).
$$

Direct differentiation reproduces the Jordan-frame improved-stress law
$\nu_{-,s}=x^2(\beta_-q_--vP_-)/2$. Here $q_-=t_{uN,-}$ and
$P_-=t_{NN,-}$ are the flux and normal pressure of
$t_{ab}=\ell^4T^\phi_{ab}/\varepsilon^2$. The extra vacuum terms
in $L_-$ and $Q_-$ are essential; the exterior's zero-vacuum shortcut
cannot be carried into the cap. The evaluation uses first scalar
derivatives, while covariant Hessians enter the independent transport
check. No derivative of sampled waveform data is used.

For this conformal scalar, $f+3f'^2=1/(2\kappa_4)$ reduces the
scalar-tensor pure-tension junction to the same signed angular law.
With $D=\mu_+-\nu_-$, its linearization gives

$$
\delta\beta_\pm=\beta_{\pm,x}\rho+\frac{D}{bx^2},\qquad
\rho_s=v_x\rho+\frac{\mu_+/x+\beta_+D/(bx^2)}{v}.
$$

The cap expression replaces the numerator by
$\nu_-/x+\beta_-D/(bx^2)$ and agrees exactly. Because $v=x_s>0$
on $x\in[0.45,0.55]$, the equation stays regular through
$\beta_+=0$. Its homogeneous solution is $C v(s)$, a relative
proper-time origin shift removed by the shared quiet event.

The existing scalar junction identities give equal one-sided fluxes and
pressures, $q_+=q_-=q$ and $P_+=P_-=P$. Subtracting the two mass laws
yields $D_s=-bx^3q/2$. Together with the displacement equation, this gives
equal corrections to the dimensionless normal accelerations $a_\pm$ and
zero change in the dimensionless mean extrinsic-curvature trace:

$$
\delta a_+=\delta a_-=\frac{6m\rho}{bx^4}-\frac{2D}{bx^3},\qquad
\delta\overline{K}=0.
$$

Thus the temporal jump and normal-force condition close with the angular
condition and mass transport, without setting $\rho=0$ or adding a freely
chosen cap mass or wall force. These are formal quadratic identities,
checked by independent algebra, not an evaluated displacement curve or
a completed sourced spacetime.

The next connection is the smooth bulk metric on both sides; the exterior
characteristic calculation below supplies its first evolved domain. Geometric
mass alone does not determine the null lapse, the wall's null-coordinate
embedding or the central observer's clock. Their reconstruction must use
the shared wall data and characteristic initial data before calculating
the corrected arrival map and subsequent scalar response. A nonzero
allowed amplitude, a controlled remainder, formation and stability remain
separate questions.

### 13.10. The prepared field drives a smooth exterior metric response

<!-- nsc-claim:atlas-bulk-response -->

Evidence: [metric refinement](../results/development/nsc-atlas-bulk-response-v1/metric/response.json),
[field refinement](../results/development/nsc-atlas-bulk-response-v1/field/response.json),
[derivative-quadrature comparison](../results/development/nsc-atlas-bulk-response-v1/derivative/response.json)
and [bulk-response owner](../docs/nsc-atlas-bulk-response.md), with separate
[geometry](../results/development/nsc-atlas-bulk-response-v1/geometry/response.json),
[source quadrature](../results/development/nsc-atlas-bulk-response-v1/source/response.json)
and [time tolerance](../results/development/nsc-atlas-bulk-response-v1/time/response.json)
comparisons.

The gravitational data of Sections 13.5--13.6 determine a computed smooth
exterior response at quadratic order on $[U_0,U_a]\times[V_a,V_b]$, whose two past null
sides are specified. This rectangle lies inside the retained wedge. Its
final transverse leg $U=U_a$ supplies data for a later extension toward
the moving wall; no exterior metric is continued through the excised cap.

Keep the original null coordinates and fixed Jordan source lapse. With
$r_0=R_0/\ell$, write $r_E=r_0+\alpha H$ and
$\log C_E=\log C_0+\alpha c$ to the retained order. The leading scalar
is $\phi=(\varepsilon/\ell)p$, with both first derivatives obtained
from its characteristic integral equations. The existing smooth exterior
Einstein equations simplify under $y=r_0H$, $z=c+H/r_0$ to

$$
y_{UV}=\frac{C_0}{4r_0^2}y-\frac{C_0}{4}z,\qquad
z_{UV}=-\frac{3C_0}{4r_0^4}y+\frac{C_0}{4r_0^2}z-p_Up_V.
$$

Integrating their right sides in $U$ gives a linear evolution in $V$.
The quiet bottom has $y=z=0$. At $U=U_0$, let $x=r_0(U_0,V)$ and
let $h_{\mathrm{init}}$ denote the original Jordan radius coefficient from
Section 13.5. The source values in the unchanged coordinate gauge are

$$
y_0=xh_{\mathrm{init}}-\frac{x^2p^2}{12},\qquad
z_0=\frac{h_{\mathrm{init}}}{x}-\frac{p^2}{4}.
$$

The frame and affine-parameter corrections must both be retained;
the Einstein-affine radius coefficient alone is not the source value at
fixed original coordinates. Returning to the Jordan metric gives
$h_J=H+r_0p^2/12$ and $c_J=c+p^2/6$, where
$r_J=r_0+\alpha h_J$ and $\log C_J=\log C_0+\alpha c_J$.

At $(U_a,V_b)$, the 65-node metric calculation with 257 field nodes
reports the following dimensionless quadratic coefficients:

| Quantity | Recorded coefficient |
| --- | ---: |
| Radius, $h_J$ | $-5.7131513827$ |
| Logarithmic null lapse, $c_J$ | $-1.5461215575$ |
| Geometric mass, $\mu_J$ | $-28.8833375991$ |

Changing the metric grid from 33 to 65 nodes reduces the largest sampled
mass-transport residual from $1.65\times10^{-6}$ to
$4.31\times10^{-7}$. On the common sampled grid, the changes in
$h_J,c_J,\mu_J$ are respectively at most $1.15\times10^{-7}$,
$3.82\times10^{-6}$ and $1.23\times10^{-6}$. A separate 129/257
field-grid comparison changes these coefficients by at most
$2.14\times10^{-8}$, $9.55\times10^{-8}$ and $5.78\times10^{-7}$;
the 32/64 derivative-quadrature comparison changes them by less than
$7.2\times10^{-15}$. Both null constraints are checked in integral
form, and geometric mass is compared with transport along both null
directions. The radial lapse reconstruction checks the same $U$ constraint
in different variables.

Changing geometry interpolation and source quadrature changes these three
coefficients by less than $7.2\times10^{-15}$ and $1.5\times10^{-13}$
respectively; tightening the time tolerance tenfold changes them by less
than $1.21\times10^{-11}$. Each complete two-solve comparison took
between 17.16 and 21.91 process-CPU seconds within its 120-second budget.

These are finite numerical comparisons, not rigorous error enclosures or
a selected nonzero physical amplitude. They establish a sourced exterior
metric response on the rectangle and provide its transverse radius/lapse
data. Section 13.11 continues to the moving wall using the same source
and material junction. The cap metric and central clock, full sourced
configuration, corrected observer signal and finite-amplitude remainder
have not yet been established.

### 13.11. The same field displaces the wall and changes its exterior clock

<!-- nsc-claim:atlas-wall-evolution -->

Evidence: [field refinement](../results/development/nsc-atlas-wall-evolution-v1/field/response.json),
[radial quadrature](../results/development/nsc-atlas-wall-evolution-v1/radial/response.json),
[clock integration](../results/development/nsc-atlas-wall-evolution-v1/time/response.json)
and the [wall-evolution owner](../docs/nsc-atlas-wall-evolution.md).

The leading field now supplies a computed response on the moving boundary,
with the action, source and quiet corner of Sections 13.5--13.10 unchanged.
Extending its first-derivative evaluator to the retained wedge requires the
correct lower endpoint: a column with $U>U_a$ begins where it meets the wall.
At that event $\chi_U=0$, while $p_U=-p r_{0,U}/x$ is generally nonzero.
Both derivatives are obtained from the integral equations of the original
field solution.

Define the areal-null lapse $Z=C/(2r_U)$. The Einstein-frame mass and
logarithmic lapse coefficients follow radial constraints at fixed $V$:

$$
\mu_E=\mu_{E,0}-\int_{U_0}^{U}\frac{r_0^2r_{0,V}}{C_0}p_U^2\,dq,
\qquad
\zeta_E=\zeta_{E,0}+\int_{U_0}^{U}\frac{r_0}{2r_{0,U}}p_U^2\,dq.
$$

Their source constants are already fixed by the characteristic data.
At the wall, frame conversion gives $\mu_+^{\mathrm{radial}}=\mu_{E,W}+Q_+$
and $\zeta_{J,W}=\zeta_{E,W}+p^2/4$. The independent timelike mass
functionals of Section 13.9 determine $\mu_+$ and $\nu_-$, and hence the
quiet-anchored displacement $\rho$. The radial and timelike mass values are
compared, rather than equated by the numerical implementation.

The exterior embedding also follows from these data. Write
$V_W=V_0(1+\alpha\eta)$, with $s$ the original Jordan wall proper time,
$d=v-\beta_+>0$, and $Z_0=-2a/V_0$, $a=2m$. Expanding
$V_{W,s}=1/[Z_J(r_{W,s}-\beta_+)]$ gives

$$
\eta_s=\frac{1}{2ad}\left[
\zeta_{J,W}+\frac{\rho_s-\delta\beta_+}{d}\right],
\qquad \eta(s_a)=0.
$$

The background lapse is evaluated at the moved $V_W$ in this expansion.
The equation stays regular through the zero of the signed $\beta_+$.

At the final reference radius $x=0.55$, the 257-node field calculation gives:

| Dimensionless quadratic coefficient | Recorded value |
| --- | ---: |
| Exterior mass, $\mu_+$ | $-28.8834223256$ |
| Cap excess mass, $\nu_-$ | $17.7590784205$ |
| Wall displacement, $\rho$ | $-5.2868390417$ |
| Exterior clock ratio, $\eta=V_1/V_0$ | $-0.8178125677$ |

For example, the radius expansion is $r_W=0.55+\alpha\rho+O(\alpha^2)$;
the tabulated number is its coefficient, with no finite $\alpha$ selected.
The cap excess subtracts its vacuum mass at the actual wall radius.

Changing the leading field from 129 to 257 nodes reduces the largest sampled
mass-path discrepancy from $6.49\times10^{-7}$ to $1.63\times10^{-7}$.
The maximum changes at nine common time cuts in $\mu_+,\nu_-,\rho,\eta$
are respectively $3.38\times10^{-7}$, $2.09\times10^{-7}$,
$3.64\times10^{-8}$ and $3.57\times10^{-9}$. Doubling the clock panels
from 8 to 16 changes $\eta$ by at most $2.41\times10^{-8}$.
The radial 2/4-point comparison changes it by less than $1.35\times10^{-13}$.
The angular jump, two forms of the displacement equation and clock
normalization are checked on the same trajectories.

The [32/64 field-derivative comparison](../results/development/nsc-atlas-wall-evolution-v1/derivative/response.json)
leaves the eight reported comparison quantities identical in this run.
Separate [geometry](../results/development/nsc-atlas-wall-evolution-v1/geometry/response.json),
[source](../results/development/nsc-atlas-wall-evolution-v1/source/response.json)
and [ODE tolerance](../results/development/nsc-atlas-wall-evolution-v1/ode/response.json)
comparisons change them by less than $5.64\times10^{-11}$,
$2.85\times10^{-14}$ and $8.25\times10^{-13}$ respectively. All seven
comparisons completed within their individual 120 process-CPU-second budgets.
An [independent rectangle overlap](../results/development/nsc-atlas-wall-overlap-v1/response.json)
finds mass and lapse differences of at most $4.08\times10^{-7}$ and
$1.30\times10^{-6}$ against the mixed metric solver. Different radial
quadratures expose error that a shared-rule constraint residual can cancel;
that residual alone is not a total metric error estimate.

These computations connect the sourced exterior and the two-sided mass
functionals to an actual wall displacement and exterior null embedding at
quadratic order. Their numerical differences are refinement diagnostics.
The smooth cap metric, central arrival correction, next scalar waveform and
finite-amplitude remainder are the remaining connections to the interior
observer.

### 14. A separate finite model with geometric feedback

<!-- nsc-claim:finite-transfer-memory -->

Evidence: [focused companion snapshot](finite-regeneration/evidence/snapshot.json),
[finite episode](../results/development/nsc-spherical-conformal-episode-v2.json),
[local response](../results/development/nsc-conformal-local-response-v2.json),
and [memory controls](../results/development/nsc-conformal-memory-control-v1.json).

The finite companion supplies another part of the programme: fields
transfer between selected regions while the same declared action evolves
the geometry. It uses a periodic radial chart, positive chart variables
and antiperiodic spinors. In the present signature convention its metric is

$$
ds^2=r^2[-L^2dt^2+Q^2(dx+\beta\,dt)^2+d\Omega^2].
$$

The physical lapse is $rL$ and radial scale is $rQ$.
Conformal gauge has $L=Q$, $\beta=0$; the reduced field Hamiltonian
is $H=\sigma_2P+\kappa Q\sigma_1$.
Here $\kappa=1$ labels the retained angular block and is independent of
the wall-tension parameter used in the atlas.
Six prepared Gaussian fermion columns determine

$$
C=\Phi\,\mathrm{diag}(0.75,0.75,0.5,0.5,0.25,0.25)\Phi^\dagger,
\qquad \Phi^\dagger\Phi=I_6.
$$

This is a finite admissible state with $0\leq C\leq I$ and
$\mathrm{Tr}\,C=3$. Unoccupied retained bands are part of the
preparation, not a derived physical vacuum. The action coefficients,
magnetic sector, clock and cutoff are recorded inputs.

Four matched discretizations, using 256/512 field modes and time steps
0.001/0.0005, reach coordinate $T=0.3$ without resetting the state.
The reported transfer and curvature comparisons have sub-percent refinement
indicators in their stated observables. The finite evolution supplies a
same-action example, not a continuum error theorem, an eternal engine or
the creation of new universes. Collars of this periodic geometry are
regions of one evolving model.

Its local reduction keeps the original first two mode columns and eliminates
their Hilbert-space complement. That complement includes modes with support
inside and outside the spatial observation window. Consequently it is not
the parent end of Section 10, and its Hamiltonian cannot be transferred to
the wall model without a field-theoretic identification.

### 15. Memory is a consequence of elimination, not an added reservoir

For a time-dependent Hermitian block Hamiltonian with retained block $A(t)$,
coupling $B(t)$ and complementary block $D(t)$, let $U_E(t,s)$
propagate $D$. Eliminating the complementary amplitude gives

$$
y(t)=U_E(t,t_*)y_*-
i\int_{t_*}^{t}U_E(t,s)B^\dagger(s)x(s)\,ds,
$$

$$
\dot x(t)=-iA(t)x(t)-iB(t)U_E(t,t_*)y_*-
\int_{t_*}^{t}K(t,s)x(s)\,ds,
$$

$$
K(t,s)=B(t)U_E(t,s)B^\dagger(s).
$$

The memory kernel, actual initial complementary drive and initial cross
correlations all matter. If $R_A,R_E$ denote the corresponding solution
maps, the retained covariance is

$$
C_A=R_A C_{AA}R_A^\dagger+R_E C_{EE}R_E^\dagger+
R_A C_{AE}R_E^\dagger+R_E C_{EA}R_A^\dagger.
$$

The companion uses the actual evolved state at $t_*=0.2$, then compares
the reduced and full histories through $t=0.3$. The selected occupations
change from approximately $(0.10557448,0.10557464)$ to
$(0.03214496,0.03226900)$. The maximum reduced/full occupation
discrepancy is $6.35\times10^{-6}$, compared with a physical-model
change of about $0.07343$.

| Conditional omission control | Maximum occupation effect | Its separate numerical reference discrepancy |
| --- | --- | --- |
| Omit the memory term | 0.0168185 | $2.20\times10^{-6}$ |
| Omit the actual initial complementary drive | 0.0485412 | $1.66\times10^{-5}$ |
| Omit initial cross correlations | 0.0920614 | $1.23\times10^{-5}$ |

These controls hold the generated geometry fixed. They are not separately
coupled physical evolutions; in particular, deleting the memory term is a
diagnostic alteration, not a new Hermitian theory. They demonstrate why an
interior description must retain the appropriate history and state
information. The next bridge to the atlas would derive such a regional
Hamiltonian and state from a particular matched spacetime and field, rather
than labeling these finite modes as universes.

## VI. What the map can ask next

### 16. Constants, spectral scales and physical inheritance

<!-- nsc-claim:target-constants -->

**Open physical target.** A shared normalized operator can preserve
dimensionless spectral relationships while regional clocks and absolute
frequency units differ. To turn this into a statement about constants,
the attachment must specify clocks, rods, fields and the transition of their
dimensionless parameters. Useful comparisons include
$\alpha$, mass ratios, and $Gm^2/(\hbar c)$.
A numerical change caused solely by unit rescaling is not a changing law.

Successive inheritance maps also need compatible output and input spaces:
the first transition must supply data on which the next law is defined.
Drawing another nested sphere does not establish this composition.

The atlas wall parameters $H\ell$, $M/\ell$ and $\kappa\ell$
are declared inputs, not values inherited from a parent calculation.
Physical inheritance might preserve a law while changing states, select
among admissible states, or alter an effective environment. These
possibilities need different maps and should remain distinguishable.

The foundation also contains a scalar-link construction whose invariant
sector has the ordinary Dirac relation
$E^2=|\mathbf p|^2+\Phi^2$. It is a useful benchmark relating a
coupling to a gap. It does not derive the observed particle species,
nuclear masses or their interactions. Particles are not assumed to be
little universes. A particle interpretation requires the appropriate field
representation, observables and dynamics, not resemblance between a spectral
scale and a cosmological radius.

### 17. Wavelength, preparation and the role of light speed

The derived scales

$$
\ell_P=\sqrt{\frac{\hbar G}{c^3}},\qquad t_P=\frac{\ell_P}{c}
$$

are defined combinations of constants, as tabulated by
[NIST CODATA](https://physics.nist.gov/cuu/Constants/).
Dimensional construction alone does not establish a smallest physical
wavelength or an indivisible spatial cell. A claimed minimum would need
a model of probes, observables and gravitational backreaction.
Minimal-length scenarios themselves distinguish limits on resolution,
distance products and other quantities; see
[Hossenfelder, 1203.6191v2, Sections 5.3–6](https://arxiv.org/html/1203.6191v2).

Absolute zero refers to thermodynamic temperature and the limiting thermal
state, not to an endpoint of allowed wavelengths; see
[NIST's temperature definition](https://www.nist.gov/si-redefinition/kelvin/kelvin-thermodynamic-temperature).
A relation between a regional temperature and a spectral scale requires a
specified state and Hamiltonian. Changing a temperature unit supplies neither.

A horizon is not automatically a reflecting cavity wall or a prohibition
on longer-wavelength modes. A particle-like exterior description likewise
requires an additional quantum state and measurement model; exterior mass,
spin and charge alone do not establish that description. Quantum observables
can have continuous as well as discrete spectra. A detector event or a wave
crest is not a general definition of a particle.

Separately, $T=0$ in an evolution record names a prepared state and a clock
origin. It is not the short-wavelength endpoint of a spectrum.
An ultraviolet cutoff, an initial time and a spatial boundary have different
mathematical roles.

A local freely falling observer's null-cone speed, coordinate propagation
rate, group velocity, phase velocity and a reduced response delay are also
different quantities. A varying coordinate speed near a horizon is not
evidence that the local causal constant changes. Claims about different
regional $c$'s require an operational, dimensionless comparison and a
consistent theory. [Ellis and Uzan, gr-qc/0305099v2](https://arxiv.org/html/gr-qc/0305099v2)
explain why the distinct roles often assigned the symbol $c$ should be
separated. The present atlas uses conventional local null cones; its clock
conversion and delayed response do not change that local law.

### 18. Dark gravitational effects require a stress and an observable

<!-- nsc-claim:target-dark-sector -->

**Open physical target.** The unexplained gravitational and expansion
phenomena motivating the project must be investigated after an exterior
model supplies consequences. A prescribed fraction, including “95%,” is
not inserted as an exterior abundance or a target that chooses the map.

Clustering and lensing require a stress or effective gravitational response
with the right spatial, temporal and scale dependence. Accelerated
expansion requires a separate background and perturbation analysis.
Antimatter and clock arrows require still other data.
The calibration scalar pulse has no calculated energy flux or backreaction and
leaves no persistent central field in its certified interval. It supplies
none of these dark components.

The foundation's Schur response and the library's shadow-matter and
brane examples remain conditional starting points. A successful next
connection would carry an actual source/state through the admitted geometry
into a gauge-invariant gravitational observable. It would then either
distinguish that exterior from alternatives or demonstrate a stated
degeneracy. Fitting a freely chosen incident waveform cannot supply that
discrimination.

### 19. Matter, antimatter and orientation

<!-- nsc-claim:target-antimatter -->

**Open physical target.** “Inside” and “outside” are geometric labels.
Charge conjugation transforms fields and their gauge representations.
Changing radial orientation or the sign of a chart coordinate does not,
by itself, convert matter into antimatter. A proposed boundary identification
must say which field components and charges it matches and whether the
source, state and conserved currents permit it.

A real wave's sign change is not charge conjugation. CPT symmetry is also
distinct from time reversal alone; neither supplies an entropy arrow.

The library preserves distinct CPT, black-mirror and torsion-production
proposals. They cannot be combined by assigning one side of every sphere
to antimatter. A useful next test would choose one actual attachment,
write its field transformation and determine whether it selects a charge
asymmetry from the stated initial data. The current real scalar has no
particle–antiparticle charge degree of freedom and cannot settle this.

### 20. Time, recycling and the next generation

<!-- nsc-claim:target-continuation -->

**Open physical target.** A continuing child needs a domain beyond a
single initial slice and a source whose evolution remains meaningful there.
The present wall reference has a singular past; the gluing theorem supplies
initial data only. Neither is a proof that an astrophysical collapsing
black hole forms the mapped cap, or that the cap is stable under disturbances.
Regular background curvature and perturbative stability are separate tests.

Proper times can be compared along specified worldlines with a matching
convention. Cross-universe ages cannot be read from unrelated coordinate
times or unrecorded clock origins. Expansion, entropy increase, causal
orientation and an observer's records may all participate in an arrow of
time, but a clock conversion alone establishes none of those thermodynamic
claims.

<!-- nsc-claim:target-recursion -->

**Open physical target.** Finite-depth inheritance and periodic regional
transfer do not prove cosmic recycling or an eternal source. A new
generation requires an identified transition, source/state transfer and
entropy accounting. Multiple holes motivate multiple possible offspring;
they do not require equal children or communication between siblings.

Known branching mathematics should be borrowed. In the finite-type
illustration used by [Altenberg, 1302.1293v2](https://arxiv.org/html/1302.1293v2),
an offspring law determines a mean matrix $N_{ij}$; conditional expected
populations evolve by $N^n$. Its long-time conclusions require the stated
branching and spectral hypotheses. The unresolved physics is the transition
law and its admissible parameter/state inheritance. Choosing that law to
favor an intended cosmic population would assume the result.
No infinite-generation computation is needed before one such transition
has been established.

Conformal cyclic cosmology remains a separate possible continuation, with
its own conformal matching and matter assumptions. It is not identified
with black-hole offspring or communicating mouths by a scale analogy.

### 21. Measurements that could decide between maps

<!-- nsc-claim:target-measurement -->

**Open physical target.** A discriminating prediction must survive changes
of units, chart labels and unspecified initial data. The protected vacuum
interval proves a finite geometric degeneracy. The scalar inverse in
Section 13.3 determines the preparation for a fixed geometry and supplies
preparation-class exclusions. It does not establish equivalence between
two physical exterior models. Both results identify what information the
next observer comparison must supply.

The next admissible calculations therefore have concrete decisions:

| Missing connection | Finite decision | What it would unlock |
| --- | --- | --- |
| Exterior rectangle to the shared wall and cap | Can the computed transverse metric data be extended with the two mass functionals and quiet-anchored displacement, then matched to the smooth cap? | The wall embedding and central clock needed for a corrected observer response |
| Preparation-restricted observer inference | Does the source reconstructed from a measured trace satisfy an independently justified preparation law in each candidate? | A discriminating response-range test or an explicit equivalence class |
| Alternative coupling or source | Does a fixed minimal-coupling or sourced model leave a distinguishable tail? | A comparison of interior response mechanisms |
| Physical regional Hamiltonian | Does a chosen matched field yield the required reduction and initial covariance? | Memory with an identified geometric exterior |
| Stability or formation of an attachment | Do fixed-source disturbances or collapse data remain in the claimed domain? | A stronger cosmological interpretation |
| Dimensionless inheritance | Which transition quantities are preserved or changed? | A meaningful comparison of constants across regions |

These are candidates, not a compulsory queue. Each proposed calculation
must state what is borrowed, what is missing, what it decides, what it
unlocks, its finite domain and budget, and the outcome that stops or
redirects it. The [existing admission owner](../docs/prior-art-and-open-claim.md)
keeps that rule. A failed attachment remains visible, while another
geometric route can proceed under its own assumptions.

![Research prerequisites and the questions they unlock](../results/development/nsc-cosmological-atlas-v7/research-dependencies.png)

*Figure 10. Research dependencies, distinct from spacetime causal arrows.
The [v8 programme record](../results/development/nsc-cosmological-atlas-v8/programme.json)
names prerequisites, decisions, budgets and stop conditions for eight themes.
Each calculation retains its own domain; arrows do not declare questions solved.*

The current result is a continuous mathematical story at a deliberately
finite scale: an exterior can be attached on a sourced slice; an isolated
completion determines observer access; a specified field state gives an
interior waveform; independent finite exterior data now predict a pulse
and a positive finite scattering tail; and finite regional reduction explains how unresolved
degrees of freedom can retain memory in another fully stated model.
The new inverse shows how a complete finite signal selects its preparation
within a fixed geometry, and why physical source restrictions must accompany
any comparison of maps.
The quadratic stress check now connects that preparation to compatible local
wall conditions. Its characteristic constraint requires a nonzero initial
radius correction under the declared corner and null normalization;
transverse constraint transport then supplies an invariant geometric mass coefficient
on that same surface. The vanishing normal-stress jump adds a necessary
condition on their later wall matching. Known field traces also determine
an analytic exterior-side wall-mass functional without a new field solve.
Subtracting the cap vacuum at the actual radius now supplies the other
mass functional and the quiet-anchored wall displacement equation.
The temporal and normal conditions agree with their common mass transport.
The prepared field now drives an evolved exterior radius and lapse on a
finite characteristic rectangle, with mass-transport and null-constraint
checks. The retained-wedge evaluation now compares independent radial and
timelike mass paths, computes the two-sided wall displacement and fixes its
exterior clock correction. Reconstruction of the cap and its central clock
is the next step toward the corrected arrival map and scalar signal.
The remaining work is to connect those statements physically where their
domains presently differ.

## Appendix A. Evidence, attribution and preserved mathematics

This living manuscript draws from source-bound atlas records without changing
their historical bytes. Record residuals and refinement comparisons are
reported at their declared scope; floating-point quadrature agreement is
not a rigorous error bound. The overview producer reads existing records
and renders their relationships without launching a new physics calculation.

| Result family | Owning derivation and preserved evidence |
| --- | --- |
| Finite inheritance and nested reduction | [Proof](../docs/nsc-nested-qualities.md); [snapshot](nested-quality-evidence/snapshot.json) |
| Two-hole parent and radiation attachment screen | [Atlas v2](../docs/nsc-cosmological-atlas.md); [manifest](../results/development/nsc-cosmological-atlas-v2/manifest.json) |
| Signed moving wall and gluing literature | [Atlas v3](../docs/nsc-cosmological-atlas.md); [manifest](../results/development/nsc-cosmological-atlas-v3/manifest.json) |
| Sourced slice and vacuum observer domain | [Atlas v4](../docs/nsc-atlas-connection-slice.md); [manifest](../results/development/nsc-cosmological-atlas-v4/manifest.json) |
| Isolated causal completion and clocks | [Atlas v5](../docs/nsc-atlas-causal-classification.md); [manifest](../results/development/nsc-cosmological-atlas-v5/manifest.json) |
| Classical incoming state and scalar transfer | [Atlas v6](../docs/nsc-atlas-field-transfer.md); [manifest](../results/development/nsc-cosmological-atlas-v6/manifest.json) |
| Finite exterior preparation and research dependencies | [Preparation](../docs/nsc-atlas-field-preparation.md); [programme](../docs/nsc-cosmological-atlas.md); [v7 manifest](../results/development/nsc-cosmological-atlas-v7/manifest.json) |
| Evaluated exterior preparation and finite tail | [Response](../docs/nsc-atlas-prepared-response.md); [v8 manifest](../results/development/nsc-cosmological-atlas-v8/manifest.json) |
| Fixed-model source inverse and preparation exclusions | [Proof](../docs/nsc-atlas-observer-inference.md); [record](../results/development/nsc-atlas-observer-inference-v1/inference.json) |
| Conformal probe surface stress and tangential flux | [Derivation](../docs/nsc-atlas-scalar-wall-consistency.md); [record](../results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json) |
| Characteristic gravitational initialization | [Constraint and frame check](../docs/nsc-atlas-characteristic-initialization.md); [record](../results/development/nsc-atlas-characteristic-initialization-v1/initialization.json) |
| Geometric mass on the prepared null surface | [Mixed geometry and stress checks](../docs/nsc-atlas-source-mass.md); [record](../results/development/nsc-atlas-source-mass-v1/mass.json) |
| Conformal scalar normal-force compatibility | [Projection and force sign](../docs/nsc-atlas-scalar-normal-force.md); [record](../results/development/nsc-atlas-scalar-normal-force-v1/normal-force.json) |
| Exterior wall-mass functional | [Both fluxes and full frame conversion](../docs/nsc-atlas-wall-mass-closure.md); [record](../results/development/nsc-atlas-wall-mass-closure-v1/wall-mass.json) |
| Cap mass and wall displacement functional | [Vacuum subtraction and junction closure](../docs/nsc-atlas-cap-wall-response.md); [record](../results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json) |
| Smooth exterior metric response | [Characteristic evolution and radial checks](../docs/nsc-atlas-bulk-response.md); [metric](../results/development/nsc-atlas-bulk-response-v1/metric/response.json), [field](../results/development/nsc-atlas-bulk-response-v1/field/response.json) and [derivative](../results/development/nsc-atlas-bulk-response-v1/derivative/response.json) comparisons |
| Moving wall and exterior clock | [Mass paths, displacement and clock](../docs/nsc-atlas-wall-evolution.md); [field](../results/development/nsc-atlas-wall-evolution-v1/field/response.json), [radial](../results/development/nsc-atlas-wall-evolution-v1/radial/response.json) and [clock-time](../results/development/nsc-atlas-wall-evolution-v1/time/response.json) comparisons |
| Finite spherical transfer, feedback and memory | [Companion](finite-regeneration/main.tex); [frozen build/evidence](finite-regeneration-manifest.json) |
| Foundation boundary, spectral and source constructions | [Archived v0.26.0](archive/foundation-v0.26.0/paper/nested-space-cosmology.md); [research inventory](../docs/research-inventory.md) |

The archived foundation retains, among other strands, its full-spinor
boundary maps, scalar-link Dirac benchmark, source and work ledgers,
curved spectral constructions and explicitly limited black-universe
benchmarks. Their theory, state and domain remain with their owning notes.
An RN pressure diagnostic is not promoted to mechanical balance, a reached
inner horizon is not promoted to central regularity, and stopped numerical
runs without completed summaries supply no new physical conclusion.
These distinctions preserve useful work without making historical routes
mandatory.

The [arXiv library](../docs/arxiv-studies.md) records paper versions,
reading depth, reusable equations, assumptions and unresolved connections.
It includes the sources for the supplementary multiple-daughter,
multi-mouth, CPT, torsion and regular-interior branches. A bounded search
cannot certify that no one has already closed an apparent gap.
The [paper guide](../docs/papers.md) locates the unchanged atlas working
drafts and the focused companion. The author's responsibility includes the
interpretations and synthesis here; borrowed results retain their source
attribution and are not presented as new laws.
