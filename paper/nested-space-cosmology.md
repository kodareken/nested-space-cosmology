# Nested-Space Cosmology

## The Geometry of a Universe

**Douglas Ek**  
Working research manuscript, 12 October 2026

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
the wall and produces a calculable transient at the center. Its waveform does
not uniquely identify the exterior because the incident state is an independent
input. A separate finite spherical field model demonstrates regional transfer,
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

The first motif contains two spatially separated black holes in a common
parent, possible unequal children, and one further generation. A separate
motif lets one hole have multiple offspring. These pictures suggest questions;
their arrows do not supply an equation of motion, a formation mechanism or
a communication channel. Figure 1 distinguishes the broad motif from the
specific attachment whose mathematics has been checked.

![Nested regions and the checked geometric connection](figures/geometry-overview.png)

*Figure 1. A finite research map and its current mathematical foothold.
The left panel records possible descent through two holes and a further
generation; those child identifications remain hypotheses. The right panel
follows the checked cap–wall–bridge body and its isolated observer/field
reference. Colors distinguish attachment, descent, causal influence and
correlation. Layout and circles are diagrammatic: they are not a spatial
embedding, scale comparison or proof of global spherical topology.
The [figure record](figures/geometry-overview.json) identifies the source
records and the status of each displayed connection.*

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
| [Atlas v1–v6 history](../docs/papers.md) and this living manuscript | Geometric compatibility, a sourced slice, observer access and a conditional scalar response | Exterior preparation, other attachment laws and observational discrimination |

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
and [owning analysis](../docs/nsc-cosmological-atlas-v2.md) retain the
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
residuals. The large value of $\dot R$ is an areal-radius rate, not a
local three-speed. This published wall family is classified more fully in
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
The [owning argument](../docs/nsc-atlas-connection-observer.md) uses standard
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

Evidence: [field transfer and checks](../results/development/nsc-cosmological-atlas-v6/field-transfer.json),
[distributional matching](../docs/nsc-atlas-field-matching.md), and
[finite initial-state construction](../docs/nsc-atlas-field-transfer.md).

A causal path permits a response; a field law and state determine one.
On the fixed wall geometry, take a classical real massless scalar in four
dimensions, its spherical $l=0$ sector, and conformal coupling
$\xi=1/6$. Its action and equation are

$$
S_\phi=-\frac12\int\sqrt{-g}\,[(\nabla\phi)^2+\xi\mathcal R\phi^2]\,d^4x,
\qquad (\Box-\xi\mathcal R)\phi=0.
$$

There is no independent scalar interaction added to the wall.
With the normal directed from cap to exterior, the pure-tension jump gives
$[K]=-3\kappa$ and singular curvature $6\kappa\delta$.
Distributional matching therefore requires

$$
[\phi]=0,\qquad [n\phi]=6\xi\kappa\phi.
$$

For $\chi=R\phi$, using $[nR]=-\kappa R$,

$$
[\chi]=0,\qquad [n\chi]=(6\xi-1)\kappa\chi.
$$

Conformal coupling cancels this radial delta barrier. It does not remove
the smooth Schwarzschild scattering potential $f_+2M/R^3$.
The jump law is borrowed from the distributional field equation; see
[Santiago et al., 1512.02120v2, equations 31–33](https://arxiv.org/html/1512.02120v2).
Their static charged-shell setup is not the moving solution used here.

Inside the conformally flat cap, regularity at the center gives

$$
\chi=F(\eta+r)-F(\eta-r),\qquad
\phi_C=\frac{2F'(\eta)}{a}.
$$

The outgoing term enforces central regularity; it is not a second
independently adjustable central source. The conformal reduction follows
the published field transformation of [Huguet, Queva and Renaud
(gr-qc/0603031v1)](https://arxiv.org/html/gr-qc/0603031v1).

Set the initial interior field and its normal derivative to zero on the
cap slice where the wall has $R=0.05\ell$. Prescribe one exterior-side
incoming characteristic at the wall,

$$
I=(d/d\tau_W+n)\chi.
$$

The plus derivative extracts the inward wave; its propagation direction
is $U_W-n$. No global exterior state or quantum vacuum is specified.
For the calibration pulse let $\delta=0.03\ell$,
$\tau_*=-1.1138558046\ell$, and $z=(\tau_W-\tau_*)/\delta$.
Define $B(z)=(1-z^2)^4$ for $|z|<1$, with $B=0$ outside, and set

$$
I=\frac{\varepsilon}{\delta}B'(z).
$$

The extended $B$ is $C^3$, and $I$ is $C^2$, not an infinitely
smooth bump. The amplitude $\varepsilon$ is a formal test-field
normalization. The incoming trace is stipulated zero outside the pulse
through the retained response interval.

Let $T(\tau_W)$ be central arrival proper time. On this early branch,

$$
T=t_W-\ell\ln(1-R/\ell),\qquad
T'=\frac{\gamma(1+v_{\mathrm{pec}})}{1-R/\ell}>0.
$$

The exact conditional central transfer is

$$
\phi_C[T(\tau_W)]=\frac{I(\tau_W)}{T'(\tau_W)}.
$$

This relation includes the wall boost and reception-clock Jacobian.
It is not a geometric-optics approximation.

![Incoming wall pulse and its central proper-time response](../results/development/nsc-cosmological-atlas-v6/field-response.png)

*Figure 5. A prescribed classical pulse, propagated on the fixed geometry.
The source record normalizes input as $\ell I/\varepsilon$ and response
as $\ell\phi_C/\varepsilon$. The incoming history is a calibration
state. Its remote preparation, energy flux, quantum state and backreaction
have not been calculated.*

The wall support lies at $R/\ell=0.47499005$ to $0.52285042$.
The central response lies at $t_C/\ell=-0.50593397$ to
$-0.34799584$, giving a duration stretch of $2.63230219$.
Its normalized peaks are $+23.84144800$ and $-24.46533478$;
the magnitude ratio is $1.02616816$. These are waveform amplitudes,
not particle abundances or energy-transfer efficiencies.

The center field is exactly zero after this received support through
$t_C=0$, conditional on the stipulated incoming history.
From the first central arrival until that endpoint, the wall's comoving
radius is at least $0.8\ell$, while the earliest reflected ray reaches
at most $0.65853382\ell$. The $0.14146618\ell$ gap excludes a
cap round trip in this interval. It does not exclude an uncomputed exterior
tail under different incoming data, erase an observer's record, or imply
that the whole field vanishes everywhere.

There is an explicit degeneracy. On any admissible monotone arrival interval,
a desired central trace $Y(T)$ can be produced by
$I(\tau)=T'(\tau)Y[T(\tau)]$, subject to the same regularity and state
restrictions. Within the zero-area pulse family, the corresponding trace
must also have zero area. Thus this freely specified pulse does not infer
a unique exterior. The right parent observer cannot prepare it through
a causal channel in this reference. Physical exterior preparation is the
next missing connection.

Minimal coupling remains a distinct alternative:
$[n\chi]=-\kappa\chi$ and the cap radial equation contains
$-2\chi/\eta^2$. Its known bulk de Sitter tail does not solve the
moving-boundary problem. The conformal transient therefore neither proves
nor excludes memory in other field theories.

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
The scalar pulse here has no calculated energy flux or backreaction and
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

### 21. Measurements that could decide between maps

<!-- nsc-claim:target-measurement -->

**Open physical target.** A discriminating prediction must survive changes
of units, chart labels and unspecified initial data. The protected vacuum
interval proves a finite geometric degeneracy. The scalar transfer proves
a second, input-versus-geometry degeneracy. Both identify information an
observer would need beyond the currently selected measurement.

The next admissible calculations therefore have concrete decisions:

| Missing connection | Finite decision | What it would unlock |
| --- | --- | --- |
| Exterior preparation of the incoming scalar | Can a specified admissible exterior state supply the wall trace after scattering? | A response predicted from exterior data |
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

The current result is a continuous mathematical story at a deliberately
finite scale: an exterior can be attached on a sourced slice; an isolated
completion determines observer access; a specified field state gives an
interior waveform; and finite regional reduction explains how unresolved
degrees of freedom can retain memory in another fully stated model.
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
| Two-hole parent and radiation attachment screen | [Atlas v2](../docs/nsc-cosmological-atlas-v2.md); [manifest](../results/development/nsc-cosmological-atlas-v2/manifest.json) |
| Signed moving wall and gluing literature | [Atlas v3](../docs/nsc-cosmological-atlas-v3.md); [manifest](../results/development/nsc-cosmological-atlas-v3/manifest.json) |
| Sourced slice and vacuum observer domain | [Atlas v4](../docs/nsc-cosmological-atlas-v4.md); [manifest](../results/development/nsc-cosmological-atlas-v4/manifest.json) |
| Isolated causal completion and clocks | [Atlas v5](../docs/nsc-cosmological-atlas-v5.md); [manifest](../results/development/nsc-cosmological-atlas-v5/manifest.json) |
| Classical incoming state and scalar transfer | [Atlas v6](../docs/nsc-cosmological-atlas-v6.md); [manifest](../results/development/nsc-cosmological-atlas-v6/manifest.json) |
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
