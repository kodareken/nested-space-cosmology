# A Finite Atlas of Nested Cosmologies: Branching, Junctions and Observer Response

**Working manuscript draft, atlas v3 - 10 October 2026.** This successor develops the
exterior-first investigation alongside the existing manuscripts. It is not a
submission-ready edition, a new release, or a replacement for their frozen
evidence. The [paper guide](../docs/papers.md) identifies those editions. The
[atlas v2 draft](atlas-working-draft.md) remains bound to its original evidence.

## Abstract

We construct a finite atlas for investigating possible relations between
black-hole exteriors, cosmological interiors and further generations of
regions. The atlas distinguishes descent, geometric attachment, causal
influence and correlation. Published geometries supply its starting points;
an unexplained line in the diagram remains an open mathematical connection.
Two-hole Brill–Lindquist initial data provide a shared spatial parent without
imposing spherical symmetry on the whole configuration. Candidate spherical
interiors are tested separately using junction conditions with a specified
surface source held fixed over a finite interval. An exact compatibility
identity excludes a constant-tension wall between a comoving closed-radiation
cap and a fixed-mass Schwarzschild region on any varying-radius interval.
Surface conservation also shows why allowing this wall to move through the
radiation cannot rescue the same source and bulk assumptions: its fixed pure
tension permits no net energy flux, while motion relative to a radiation
fluid produces a nonzero flux. A published de Sitter/Schwarzschild moving
wall supplies a separate candidate whose vacuum equation of state removes
that obstruction. These are scoped applications of known junction laws,
not new formation theorems. A published constraint-gluing theorem opens
an alternative multi-body parent with exactly retained spherical collars at
sufficiently large separation. Joining the moving wall to suitable body
initial data remains unresolved. Alternative multiple-daughter, communicating-mouth
and observationally equivalent exterior constructions remain independent.
The present work assembles imported mathematics and records which
compatibility questions it decides. It does not establish black-hole child
formation, a stable complete spacetime, inherited particle constants or a
dark-sector explanation.

## 1. Begin with the possible exterior

An interior observer measures fields, curvature and correlations within a
limited region. The larger structure that supplies those data may be harder
to determine. Nested-Space Cosmology asks whether a systematic map of possible
exteriors can make that ambiguity useful: which regions can be joined, which
can influence one another, and which leave the same local response?

Our first motif places two candidate black-hole interfaces in a common
parent. Their proposed children may differ, and one child may contain a
further interface. This is a finite question about two generations. No
infinite hierarchy is needed to ask whether a single attachment is admissible.

![The finite atlas and its distinct connections](../lab/results/development/nsc-cosmological-atlas-v3/map.png)

*Figure 1. A map of regions and the status of their proposed relations. Its
layout is diagrammatic; separation on the page is neither proper distance nor
an embedding of the complete spacetime. Candidate descendants remain visible
when an attempted source or attachment is rejected.*

The research sequence is to borrow a geometry, identify the physical data
required at its boundary, decide a finite compatibility question, and then
derive consequences for an observer. The proposed explanatory reach includes
gravitational response, expansion, antimatter and time. These are distinct
questions, each requiring a mechanism and an observable. A desired cosmic
abundance does not select the geometry or its parameters.

## 2. Four relations carried by the atlas

Each region records its geometry, matter or quantum state, clock convention,
dimensionless parameters and validity domain. A spatial initial-data set
does not claim the information content of an evolved four-geometry. Each
connection records its mathematical source, assumptions and missing conditions.

| Relation | Mathematical question | What the relation alone does not establish |
|---|---|---|
| Descent | Which transition supplies the putative child's initial data? | A joined metric or a continuing signal channel. |
| Geometric attachment | Which boundary data match, with which normal orientations and source? | Traversability or regular formation from earlier data. |
| Causal influence | Can an intervention reach the named observer along an allowed route? | That every correlated quantity can be controlled or transmitted. |
| Correlation | Which state or ancestry relates measurements in different regions? | A communication channel between them. |

A round two-surface has topology $S^2$. A candidate closed spatial slice
may have topology $S^3$. These statements have different dimensions and
neither follows from drawing one circle inside another. Similarly, multiple
ends of a spatial slice do not by themselves specify a parent-to-child
history. Causal interpretation requires a spacetime extension; the review by
[Ansoldi, Merali and Guendelman](https://arxiv.org/html/1801.04539v1) provides
the relevant distinction between a spatial bridge and a traversable route.

## 3. A shared parent with two punctures

We import the uncharged, time-symmetric Brill–Lindquist construction. In
geometrized units $G=c=1$, its initial data are

$$
\gamma_{ij}=\psi^4\delta_{ij},\qquad K_{ij}=0,\qquad
\psi=1+\frac{m_A}{2|\mathbf x-\mathbf x_A|}
       +\frac{m_B}{2|\mathbf x-\mathbf x_B|}.
$$

The domain is $\mathbb R^3\setminus\{\mathbf x_A,\mathbf x_B\}$.
Away from the punctures, the harmonic conformal factor satisfies the vacuum
Hamiltonian constraint. For coordinate separation
$d=|\mathbf x_A-\mathbf x_B|$, the three asymptotic-end ADM masses are

$$
M_\infty=m_A+m_B,\qquad
M_A=m_A\left(1+\frac{m_B}{2d}\right),\qquad
M_B=m_B\left(1+\frac{m_A}{2d}\right).
$$

These are borrowed results, not new solutions. The original source is
[Brill and Lindquist (1963)](https://doi.org/10.1103/PhysRev.131.471);
the uncharged conventions above appear in
[Pook-Kolb et al., v1, Section III, equations 17–21](https://arxiv.org/html/2308.10906v1).
The parameters $m_A,m_B$ must not be silently identified with all three
end masses or with horizon masses.

The finite example uses $m_A=1$, $m_B=1.5$ and $d=12$ in one shared length
unit. The corresponding end masses are $M_\infty=2.5$, $M_A=1.0625$ and
$M_B=1.5625$. Their inequality is recorded input and its geometric consequence;
it is not a predicted distribution of children or inherited constants.
The [parent adapter and controls](../lab/docs/nsc-atlas-parent-initial-data.md)
also compare two finite coordinate spheres with their local monopole
approximations. Those probes are not computed horizons.

![A section through the two-hole parent initial data](../lab/results/development/nsc-cosmological-atlas-v2/parent-slice.svg)

*Figure 2. A coordinate section of a borrowed initial-data geometry. The
conformal factor describes intrinsic spatial geometry; a colored contour is
not automatically a marginally trapped surface. Locating horizons and evolving
the slice are separate problems.*

The parent is axisymmetric for punctures on one axis, but it is not globally
spherically symmetric. Near one puncture, the companion contributes both a
leading offset and direction-dependent corrections to the conformal factor.
A local Schwarzschild description therefore requires a stated approximation
and an error criterion. Our spherical junction test does not supply that
criterion or remove the companion's tidal field.

On a coordinate sphere of radius $r$ about one puncture the induced metric is
$q_{AB}=r^2\psi(\theta)^4\Omega_{AB}$. The angular dependence is part of the
boundary data, even though the coordinate surface is called a sphere. The
adapter records its intrinsic Gaussian curvature as well as the conformal
factor's angular profile. These spatial quantities can test a proposed round
boundary; they are not four-dimensional observer tides.

Brill–Lindquist data carry three asymptotically flat ends. We do not relabel
the two puncture ends as expanding children. Doing so would replace a spatial
fact with a dynamical hypothesis. The atlas retains a separate, unresolved
transition from each candidate interface to its proposed cosmological region.

### 3.1 An alternative parent that preserves spherical collars

The nonround coordinate spheres of this particular parent do not exclude
exact spherical neighborhoods in every multi-body initial geometry.
[Chrusciel, Corvino and Isenberg, v1, Theorem 1.1 and its following discussion](https://arxiv.org/html/0909.1101v1)
provide a different construction. For positive-mass, scalar-flat,
asymptotically Euclidean body data, their gluing preserves selected compact
regions exactly while creating a common outer Schwarzschild end. Applied
to two Schwarzschild bodies, it can retain compact neighborhoods of their
throats. Unequal masses are allowed, and the body separation must be
sufficiently large. The proof adjusts mass and center parameters to remove
the finite-dimensional constraint obstruction. It supplies no usable numeric
separation threshold for our selected collars.

This imported existence result opens a separate parent family. It does not
make the unmodified Brill-Lindquist metric exactly spherical near either
puncture, nor certify its chosen $d=12$ example. We do not rerun the published
constraint-solving proof or label a schematic picture as its computed metric.

The remaining question is more specific: can a proposed child junction be
included in admissible body data while retaining the required vacuum end?
A moving timelike wall and a time-symmetric spatial collar are different
objects. Equality of a mass parameter does not match their induced metric,
extrinsic curvature, source and Cauchy data. The theorem supplies the parent
assembly route after those data are compatible; it does not supply that
compatibility itself. Extra ends retained from the individual bodies also
retain their spatial meaning, rather than becoming children by relabeling.

Time symmetry need not become a permanent restriction. The same authors'
[non-time-symmetric extension, v1, Theorem 1 and Corollary 1.1](https://arxiv.org/html/1004.1355v1)
retains specified body data $(g,K)$ with vacuum asymptotically Euclidean ends
and timelike ADM energy-momentum. Matter may remain away from the gluing
ends, and the common far exterior is a Kerr slice. Sufficient separation and
compatible initial data are still required. The stated smooth-data theorem
does not by itself establish applicability to an ideal Israel wall with a
distributional surface source. A source-supported smooth body construction,
or a justified extension to that regularity, is an explicit remaining
connection. The [gluing assessment](../lab/docs/nsc-atlas-parent-gluing.md)
records these routes and their reading boundaries.

## 4. Decide an attachment with its source held fixed

The first compatibility question is deliberately limited:

> Can the specified surface source support the proposed spherical
> parent–child boundary over the declared finite interval without changing
> the source law or its coefficients along the trajectory?

Let a timelike shell have induced metric

$$
ds_\Sigma^2=-d\tau^2+R(\tau)^2d\Omega^2,
\qquad S^a{}_b=\operatorname{diag}(-\sigma,p,p).
$$

For a fixed orientation convention, Israel matching requires the jump in
extrinsic curvature, including its temporal component, to reproduce this
surface stress tensor:

$$
[K^a{}_b]-\delta^a{}_b[K]=-8\pi S^a{}_b.
$$

The surface conservation equation also relates the shell's evolving energy
and pressure to flux from the bulk. When that flux vanishes it reduces to

$$
\dot\sigma+2\frac{\dot R}{R}(\sigma+p)=0.
$$

Thus an angular matching identity alone is insufficient. The density
inferred from that identity, the pressure required by the temporal equation,
and the chosen material law must agree on the same trajectory. In particular,
a pure domain wall has $p=-\sigma$; with zero flux, its density is constant.
Replacing that constant by a separately fitted value at every radius defines
a different problem.

The source of these identities is
[Israel, *Singular hypersurfaces and thin shells in general relativity*](https://doi.org/10.1007/BF02710419).
The classic de Sitter/Schwarzschild domain-wall problem was already analyzed
by [Blau, Guendelman and Guth](https://doi.org/10.1103/PhysRevD.35.1747).
Replaying a known admissible wall is an integration control. The present
investigation tests the stated atlas candidate and records a rejection just
as explicitly as a successful attachment.

![The finite junction-source comparison](../lab/results/development/nsc-cosmological-atlas-v2/junction.svg)

*Figure 3. The source comparison is evaluated along one specified trajectory
and over one finite domain. A failure of that candidate closes its attempted
connection, while leaving other sources, trajectories and geometries open.*

### 4.1 A comoving radiation cap and a constant-tension wall

The chosen interior is the exact closed, radiation-dominated FLRW solution

$$
ds_-^2=-dt^2+a(t)^2(d\chi^2+\sin^2\chi\,d\Omega^2),\qquad
a=A\sin\eta,\qquad t=A(1-\cos\eta).
$$

A boundary at fixed $\chi$ has $R=a\sin\chi$ and $\tau=t$. We compare it
with a Schwarzschild region of fixed mass $M$, restricting the finite
numerical example to an exterior static patch and a real timelike shell.
This $M$ belongs to the separate spherical test; it is not identified with
one of Section 3's asymptotic-end masses without an additional matching map.
The radiation fluid is comoving with the boundary, so its energy flux
through the boundary is zero; the exterior is vacuum. The trial wall has
$p=-\sigma$ with fixed $\sigma$.

For the ordinary cap orientation let $b=\cos\chi$. The angular extrinsic
curvatures can be expressed through

$$
\beta_-=b,\qquad
\beta_+^2=b^2+\frac{A^2\sin^4\chi}{R^2}-\frac{2M}{R}.
$$

The required angular jump is $\beta_- -\beta_+=4\pi\sigma R$.
Set $k=4\pi\sigma$. Squaring this necessary equation and multiplying by
$R^2$ gives

$$
k^2R^4-2bkR^3+2MR-A^2\sin^4\chi=0.
$$

For $A>0$ and $0<\chi<\pi$, the constant coefficient is nonzero. Therefore
this polynomial cannot vanish throughout an open interval of radii for any
fixed $k$ and $M$. This excludes a constant-tension wall on the stated
varying-radius, comoving radiation trajectory. Reversing the cap orientation
changes $b$'s sign but not that obstruction. Squaring is used only for a
necessary condition; no extra solutions from squaring are admitted.

The conclusion is stronger than observing a large numerical residual at a
few sample times, but it remains narrower than a no-child-universe theorem.
The polynomial alone tests the comoving boundary. Section 4.3 separately
checks whether noncomoving motion is available with the same bulk and wall
assumptions. Another material law, compatible bulk flux, or different bulk
geometries require their own matching problem. The radiation solution used
here has no regular bounce; its finite expanding segment is a matching test
rather than a formation history.

### 4.2 Finite witness and independent checks

The recorded test uses $M/A=0.4$, $\chi=\pi/3$ and
$\eta\in[1.25,1.45]$. Its specified positive wall density is
$\sigma=1/(64\pi A)$, chosen before the comparison. The parameters are
dimensionless model inputs, not inferred cosmological constants. The interval
has $\tau/A\in[0.684678,0.879497]$ and
$R/A\in[0.821845,0.859715]$, outside the Schwarzschild horizon $2M/A=0.8$
and before the FLRW turnaround.

The finite implementation evaluates the surface density and pressure
required by both components of the junction condition, checks conservation,
and compares them with the same fixed wall. A second normal orientation is
recorded separately. The exact polynomial obstruction prevents a different
constant wall density from rescuing the candidate on the whole interval.

| Ordinary-cap quantity | At $\eta=1.25$ | At $\eta=1.45$ |
|---|---:|---:|
| Required $A\sigma$ | 0.01638975 | 0.02001736 |
| Required $Ap$ | −0.05885743 | −0.05825351 |
| Fixed wall $A\sigma$ | 0.00497359 | 0.00497359 |
| Fixed wall $Ap$ | −0.00497359 | −0.00497359 |
| Required $A(p+\sigma)$ | −0.04246768 | −0.03823615 |

The required density varies, and the required pressure does not satisfy the
domain-wall equation of state. This is a rejection of the chosen attachment.
The required tensor is a geometric diagnostic, not an independently derived
new material source.

The [source test and its independent controls](../lab/docs/nsc-atlas-junction-sources.md)
compare the temporal extrinsic curvature with a direct Christoffel-symbol
calculation, check the surface conservation law using a five-point numerical
derivative, and check the normal Einstein jump relation. Their maximum
absolute discrepancies in the recorded dimensionless convention are
$1.21\times10^{-10}$, $9.16\times10^{-13}$ and $6.94\times10^{-18}$,
respectively. An imported Oppenheimer–Snyder dust boundary supplies a positive
control: across nine samples its required surface density and pressure are
zero to absolute discrepancies below $1.9\times10^{-17}$ and
$8.9\times10^{-18}$. These are numerical verification discrepancies, not
rigorous error bounds or physical perturbative-stability results. The
polynomial argument, rather than sampling density or a residual tolerance,
establishes the incompatibility over the interval.

Even a fully matched finite spherical shell would not establish its embedding
in the nonspherical parent of Section 3. It would also leave regular formation
and perturbative stability open. The formation obstruction of
[Farhi and Guth](https://doi.org/10.1016/0370-2693(87)90429-1) is conditional on
its classical energy and global assumptions. The nonspherical instabilities
studied by [Aguirre and Johnson, v2](https://arxiv.org/html/gr-qc/0508093v2)
show why a spherical background trajectory cannot substitute for a stability
analysis. These results guide separate checks; neither is applied outside its
stated model to exclude every atlas branch.

### 4.3 Is motion through the radiation an available freedom?

Changing the boundary trajectory is useful only if the chosen source permits
the associated energy transfer. The relevant restriction comes before a new
trajectory integration. For a pure wall with constant tension,

$$
S^a{}_b=-\sigma\delta^a{}_b,\qquad D_aS^a{}_b=0.
$$

The contracted Codazzi identity therefore requires
$[T_{\mu\nu}U^\mu n^\nu]=0$ across the wall. Here $U$ is the wall's unit
velocity and $n$ its consistently oriented unit normal. We import the
conservation identity from
[Lobo and Crawford, v2, equation 8](https://arxiv.org/html/gr-qc/0507063v2),
rather than assuming a transparent boundary independently of its matter.

Let $u$ be the interior perfect fluid's velocity and $e$ a radial unit vector
in its rest frame. A spherical wall moving at speed $v$ relative to that
fluid has

$$
U=\gamma(u+ve),\qquad n=\gamma(vu+e),\qquad
\gamma=(1-v^2)^{-1/2},\qquad |v|<1.
$$

Contracting $T_{\mu\nu}=(\rho+p)u_\mu u_\nu+p g_{\mu\nu}$ gives

$$
T_{\mu\nu}U^\mu n^\nu=(\rho+p)\gamma^2v.
$$

Reversing the normal changes the flux sign, not its zero. The Schwarzschild
side supplies no matter flux. Positive-density radiation has
$\rho+p=4\rho/3>0$, so conservation forces $v=0$. In the FLRW coordinates
of Section 4.1, $v=a\,d\chi/dt$, hence the wall must be comoving.
The polynomial obstruction then applies on every open interval of changing
radius in its stated domain. Motion alone cannot rescue this radiation,
vacuum exterior and fixed pure-tension combination.

This is a compatibility deduction from an established conservation law.
It does not rule out a wall that absorbs energy and changes its surface
state, a different prescribed material law, compensating exterior flux, or
a different interior. An imposed function $\sigma(\tau)$ chosen to cancel
the mismatch would require a physical exchange model before it could become
another admitted candidate. The atlas consequently stops the present route
without commissioning a noncomoving trajectory search under the same
assumptions.

### 4.4 A distinct moving-wall branch

The de Sitter case provides a useful comparison. Its vacuum stress satisfies
$\rho+p=0$, so the preceding fluid flux vanishes for any timelike wall
velocity. A moving de Sitter/Schwarzschild wall is an imported construction;
its availability does not repair the radiation interior. It identifies the
bulk assumption that changes the compatibility question.

In geometrized units, write

$$
f_-(R)=1-H^2R^2,\qquad f_+(R)=1-\frac{2M}{R},\qquad
\kappa=4\pi\sigma>0.
$$

With the normal directed from the minus to the plus region, the angular
junction equation is

$$
\beta_- -\beta_+=\kappa R,\qquad
\beta_\pm=\epsilon_\pm\sqrt{\dot R^2+f_\pm(R)},\qquad
\epsilon_\pm\in\{-1,+1\}.
$$

The signs specify which side of each spherical geometry is retained. They
cannot be recovered by plotting a squared effective potential alone.
Keeping those signs gives

$$
\beta_-=\frac{f_--f_++\kappa^2R^2}{2\kappa R},\qquad
\beta_+=\frac{f_--f_+-\kappa^2R^2}{2\kappa R},\qquad
\dot R^2=\beta_-^2-f_-.
$$

These are borrowed thin-wall dynamics, attributed to
[Blau, Guendelman and Guth](https://doi.org/10.1103/PhysRevD.35.1747)
and read through the sign-sensitive account of
[Ansoldi, Guendelman and Shilon, v1](https://arxiv.org/html/0711.2198v1).
The finite adapter tests an already specified branch; it does not claim to
rediscover the bubble solution. The temporal junction equation, the
unsquared angular equation and the allowed metric domain must agree with
the same tension throughout the recorded segment.

The [moving-boundary adapter](../lab/docs/nsc-atlas-moving-boundaries.md)
uses $H\ell=1$, $\kappa\ell=2$ and $M/\ell=0.4$, with
$R/\ell\in[1.2,1.8]$ on the expanding branch. Here $\ell$ is the shared
length unit. Throughout this interval $\epsilon_-=+1$ and
$\epsilon_+=-1$: the normal on the Schwarzschild side points toward smaller
areal radius. This supplies the signed local orientation, without deciding
which complete spacetime history creates it.

The de Sitter segment lies beyond its static-patch horizon, where the static
coordinate $t_-$ is spacelike. The decreasing signed Schwarzschild time
coordinate on the selected side is also not a negative observer age.
The wall's proper time is calculated from $d\tau=dR/\sqrt{\dot R^2}$;
the indefinite metrics and timelike normalization must hold on both sides.
No cross-universe clock comparison or entropy arrow is inferred from these
coordinate choices.

| Vacuum-wall quantity | At $R/\ell=1.2$ | At $R/\ell=1.8$ |
|---|---:|---:|
| $\dot R$ | 1.232594874 | 2.057420001 |
| $\beta_-$ | 1.038888889 | 1.411728395 |
| $\beta_+$ | -1.361111111 | -2.188271605 |
| Fixed $\ell\sigma$ | 0.159154943 | 0.159154943 |
| Fixed $\ell p$ | -0.159154943 | -0.159154943 |

The elapsed wall proper time is $\Delta\tau/\ell=0.3718706103$.
The quantity $\dot R$ is the rate of change of areal radius with proper time;
it is not a locally measured three-speed. Its values greater than one do
not violate the checked timelike normalization. Across 33 radii, the maximum
unsquared angular residual is $4.45\times10^{-16}$. The direct pressure
comparison differs by at most $2.78\times10^{-17}$ in the stated scaled
units; an independent Christoffel-symbol and five-point-derivative check
differs by at most $2.86\times10^{-12}$. The largest absolute tangent
normalization discrepancy is $3.56\times10^{-15}$.

Quadrature orders 32 and 64 agree on the duration within
$5.56\times10^{-17}\ell$. These are numerical verification discrepancies,
not rigorous error bounds or perturbative stability results. Replacing the
signed normal components by their absolute values produces an angular
mismatch of at least $2.72$ on this interval. The explicit orientation is
therefore consequential to the admitted attachment, rather than decorative
notation.

![Moving-boundary compatibility and its independent vacuum branch](../lab/results/development/nsc-cosmological-atlas-v3/moving-boundary.png)

*Figure 4. The moving-boundary comparison separates a conservation constraint
on the radiation candidate from a sourced vacuum-wall control. A compatible
finite segment supplies local junction data. It does not establish formation,
stability, sibling communication or attachment to the two-puncture parent.*

## 5. Alternative exteriors remain available

![Independent candidate families](../lab/results/development/nsc-cosmological-atlas-v2/candidates.svg)

*Figure 5. Alternatives have separate sources and validity domains. A failed
radiation boundary does not invalidate a limiting-curvature construction,
and a conditional traversable wormhole does not supply a neutral child's
formation mechanism.*

[Barrabès and Frolov, v1](https://arxiv.org/html/hep-th/9511136v1) supply a
one-hole multiple-daughter construction. Its null shells and collision
conditions differ from two distinct holes in a shared parent. In particular,
the cited disconnection criterion uses a near-singularity approximation,
equal nucleation radii and a scale hierarchy; it is not a universal rule for
sibling separation. Formation microphysics remains an additional question.

[Emparan et al., v2](https://arxiv.org/html/2012.07821v2) provide a
communicating-mouth candidate within charged semiclassical assumptions. This
branch offers a way to investigate an actual causal channel, provided its
matter, state and geometric conditions are retained. It does not make
communication a generic consequence of shared ancestry.

A third possibility is that several exteriors are indistinguishable to the
selected observer. The finite atlas includes a control in which changing an
exactly uncoupled exterior mode changes the full spectrum but leaves the
local resolvent unchanged. This is established Schur-complement algebra. Its
role is to expose a limit on reconstruction, then ask what extra observable
or coupling would distinguish the alternatives.

## 6. From a connection to an observer's response

For a finite Hermitian operator divided into visible and exterior sectors,
the projected resolvent has the familiar form

$$
P(z-H)^{-1}P
=\left[z-H_{\rm v}
-V(z-H_{\rm e})^{-1}V^\dagger\right]^{-1},
\qquad \operatorname{Im}z>0.
$$

The exterior enters through the self-energy
$\Sigma(z)=V(z-H_{\rm e})^{-1}V^\dagger$. This identity tells us how
to reduce a specified operator. It does not choose the physical $H$, the
quantum state, or the inter-region clock map. The
[existing finite-reduction framework](../lab/docs/nsc-cosmological-atlas.md)
also accepts a supplied time-dependent Hamiltonian. Its software fixtures
exercise that interface without identifying their coefficients with
cosmological couplings.

A causal prediction requires a compatible field problem on the admitted
geometry, a joint state, boundary conditions and an observer. Proper time is
attached to a worldline and metric. A coordinate interval, a lapse sum lacking
its trajectory, or two independently normalized clocks cannot establish a
cross-universe age. Likewise, an entropy arrow requires a state and
coarse-graining, and an antimatter interpretation requires an explicit field
transformation and interaction model.

The future observational questions are therefore conditional: which
admitted exterior produces which local response, and which competing
exterior can reproduce it? A lensing response, clustering response, expansion
history and charge asymmetry will have separate evidence and comparison
domains. No result in this draft sets their observed amplitudes.

## 7. Reproducibility and the next missing connection

The [versioned atlas bundle](../lab/results/development/nsc-cosmological-atlas-v3/manifest.json)
binds the recorded inputs, sources, tests and displayed artifacts. Numerical
replay verifies the stated calculation; byte hashes authenticate its
provenance. Neither is a substitute for a missing source model or a proof of
physical formation. The [reuse library](../docs/arxiv-studies.md) records
paper versions and reading boundaries, and the existing
[admission rule](../docs/prior-art-and-open-claim.md) requires a finite domain,
budget, independent check and branch-stopping outcome before further compute.

The retained [parent record](../lab/results/development/nsc-cosmological-atlas-v2/parent.json)
and [radiation-junction record](../lab/results/development/nsc-cosmological-atlas-v2/junction.json)
keep the unrounded predecessor values and their domains. The
[successor connection map](../lab/results/development/nsc-cosmological-atlas-v3/atlas.json)
adds the conservation decision and the separate vacuum-wall and gluing routes.
The [moving-boundary record](../lab/results/development/nsc-cosmological-atlas-v3/moving-boundary.json)
stores the source parameters, signed branch, proper-time convention and
independent checks behind Section 4.4.
The Brill-Lindquist adapter
uses the initial-data section of Pook-Kolb et al. v1; it does not reproduce
that paper's merger evolution or quasi-local-mass analysis. Historical
originals are credited without implying their full derivations were rechecked.
Versioned candidate sources retain the reading boundaries recorded in the
library. This bounded prior-art comparison makes no claim to exhaustive
coverage or priority.

The atlas gains scientific value when a connection changes status for a
stated reason. The radiation result identifies a freedom that was only
apparent: a wall cannot be set in motion through that fluid while keeping
its pure constant tension and vacuum exterior. The de Sitter case supplies
a different, already solved local dynamics. The gluing theorem supplies a
different way to assemble a common parent while preserving local geometry.
Keeping these implications separate makes the next question smaller and
more useful.

The next missing connection is a body initial-data set that contains the
selected wall and child segment and meets the asymptotic hypotheses of a
multi-body construction. A finite wall trajectory alone does not provide
that slice; a theorem retaining a vacuum collar does not supply the wall.
The task is to identify compatible slicing, surface stress and constraint
data, borrowing a published construction where one exists. Only after that
step would an evolution or a field-and-clock problem on the joined geometry
have the required physical inputs. Further material laws and flux models
remain separate admitted candidates. No large simulation is justified by
the diagram alone.
