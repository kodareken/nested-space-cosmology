# Shared parents with retained spherical collars

**Status, 10 October 2026:** imported existence route; no gluing PDE solved and
no moving child attached. The [v2 parent](nsc-atlas-parent-initial-data.md)
remains the unchanged Brill-Lindquist (BL) slice. A separate parent family can
retain exact Schwarzschild neighbourhoods while its common exterior is
nonspherical. This removes a need to invent multi-hole gluing mathematics;
it does not certify the existing BL parameters or select a cosmology.

## What can be borrowed

| Source and reading depth | Reusable result and conditions | Atlas use and remaining boundary |
| --- | --- | --- |
| [Chrusciel, Corvino and Isenberg, 0909.1101v1](https://arxiv.org/html/0909.1101v1), sections 1-3, Theorems 1.1 and 2.1, proof's obstruction equations; section 4 propositions inspected | Smooth scalar-flat asymptotically Euclidean (AE) ends with positive ADM masses admit a common scalar-flat exterior. Prescribed compact regions are retained isometrically; masses and placements need not be equal or reflection symmetric. Separation must be sufficiently large. The proof balances the four scalar-curvature moments against mass/centre parameters. | Apply to two Schwarzschild bodies while retaining their collars: **existence-backed alternative parent**, not BL. The threshold is not supplied numerically. The exterior ADM mass approaches the sum in the large-separation construction; it is not an independently fixed input. No certificate for separation 12 follows. |
| [Chrusciel, Corvino and Isenberg, 1004.1355v1](https://arxiv.org/html/1004.1355v1), section 1 Theorem 1/Corollary 1.1, section 3.1 charge-balance argument | Removes time symmetry. Smooth initial data with vacuum AE ends and timelike ADM energy-momenta can retain specified interior $(g,K)$ regions at sufficiently large separation; the shared outer end is Kerr. Matter may be supported away from the vacuum gluing ends. | Possible later route for already constructed moving-source bodies. The theorem does not manufacture their source, Cauchy slice or evolution. Its smooth-data hypothesis does not automatically include an ideal Israel delta-layer. |
| [Pook-Kolb and Giulini, 1809.00935v1](https://arxiv.org/html/1809.00935v1), section 2 and 2.1, equations 1-10; section 3.2 equations 11-13 | Corvino-type gluing replaces an AE end outside a large radius while preserving the inner data. An annular constraint correction is essential; a cutoff interpolation alone is generally not scalar flat. Their implemented example uses equal-mass BL data. The Dirichlet/Neumann compatibility and exterior mass must be solved together. | Supplies an existing numerical method if a future named uncertainty needs actual transition data. Gluing our outer end does not make the individual BL throats round. Applying the end construction after puncture inversion acts far down that AE end, not automatically at the existing throat. |
| [Chrusciel and Delay, gr-qc/0301073v2](https://arxiv.org/pdf/gr-qc/0301073v2), Theorem 8.6, section 8.7 Theorem 8.15, and section 8.9 construction | Theorem 8.6 preserves an interior and installs a time-symmetric reference end under decay/parity and nonzero-mass assumptions, allowing nearby mass/centre parameters. The localized IMP construction requires trivial Killing-initial-data (KID) kernel on the chosen annuli. The many-Kerr construction explicitly treats parity and finite charge obstructions. | Different sufficient hypotheses must remain separate. We have not certified their weighted estimates for a chosen BL annulus. The scalar, time-symmetric deformation problem and the full constraint/KID problem are distinct. The full 87-page analysis has not been audited. |
| [Chrusciel, Isenberg and Pollack, gr-qc/0403066v2](https://arxiv.org/html/gr-qc/0403066v2), section 1 equations 1.1-1.6 and Theorem 1.1 | Smooth vacuum data with no KIDs in each of two open neighbourhoods admit a localized neck construction; data outside those neighbourhoods remain exact. The theorem imposes no global completeness or asymptotic condition. | A useful topology-changing route for qualifying data. It neither identifies a cosmological child nor establishes traversal. The no-KID condition must be checked, not inferred from the word “generic.” |
| [Hintz, 2210.13960v1](https://arxiv.org/pdf/2210.13960v1), introduction, Theorems 1.1-1.2 and their qualification preceding section 1.1 | Inserts a small rescaled AE body near a point in a background with no local KIDs, with controlled limiting behaviour. Exact background agreement outside the chosen neighbourhood differs from convergence to the inserted geometry as the scale tends to zero. The Kerr exception described there permits constraint violation deep inside the black hole. | Successor route for scale-separated bodies, not an automatic exact-collar replacement in BL. We do not import the exceptional construction as an everywhere-vacuum certificate. Detailed microlocal proof and section 6.2 remain unread. |

The original Corvino paper, *Scalar curvature deformation and a gluing
construction for the Einstein constraint equations*, Commun. Math. Phys.
214 (2000), 137-189, [DOI 10.1007/PL00005533](https://doi.org/10.1007/PL00005533),
is attribution through the inspected primary implementations/theorems above;
its complete original proof has not been separately read here.

## A correction that changes which construction is usable

[Chrusciel and Delay, gr-qc/0203053v3](https://arxiv.org/html/gr-qc/0203053v3),
sections 2 and 4, constructs parity-symmetric many-Schwarzschild data with
sufficiently small masses relative to fixed geometric scales. The paired
placements and masses constrain that particular construction. Its section 5
claim that a generic locally parity-symmetric point can acquire a
positive-mass Einstein-Rosen bridge was withdrawn in the
[one-page erratum](https://homepage.univie.ac.at/piotr.chrusciel/papers/cs/erratum_proof.pdf),
read in full: the added mass could have the wrong sign because an error term
was underestimated. We use the later N-body theorem, not that withdrawn
bridge argument. The erratum does not withdraw section 4's separate
many-Schwarzschild construction.

There is also a direct hypothesis check for the unchanged BL parent. Its
axial rotation $Y=-y\partial_x+x\partial_y$ obeys
$\mathcal L_Yg=0$, and $K=0$. Consequently $(N,Y)=(0,Y)$ solves the full
KID equations quoted in gr-qc/0403066v2. Every three-dimensional open
neighbourhood contains points where this field is nonzero. Thus the exact
axisymmetric BL data fail that theorem's trivial-KID hypothesis. Exact
Schwarzschild data likewise have KIDs. This is an analytic check of a
sufficient hypothesis, **not a proof that gluing is impossible**. It does
not rule out the separate scalar-curvature or charge-balanced constructions.

## Connection meanings in the map

| Candidate | Status now | Evidence needed for a stronger statement |
| --- | --- | --- |
| Unchanged unequal BL parent | Explicit vacuum initial slice; its sampled coordinate spheres are not exactly round | A specific deformed boundary and compatible junction, or a declared change of parent data |
| Shared parent assembled from two Schwarzschild bodies | Imported existence family with retained spherical collars; no selected numerical separation or transition metric | Particular retained sets and theorem hypotheses; quantitative data only if needed by the next question |
| Moving de Sitter wall behind one collar | Local spherical candidate; shared-parent attachment unresolved | Source-supported Cauchy data, the selected Schwarzschild side and orientations, constraint matching, and an evolution domain |
| Parent with exterior Corvino modification | Exact imported spatial construction pattern; implementation deferred | Actual constraint-corrected annulus and its mass/centre if an observable depends on them |

The retained collar and its parent describe one initial geometry. An
asymptotic end or Einstein-Rosen bridge is not a formed child, a traversable
channel, or a causal message path. The theorem's placement scale is not a
cross-universe age. No local dark-sector fraction follows from these data.

## The moving-wall connection still to be supplied

Equal Schwarzschild mass parameters are necessary bookkeeping, not an
identification of four-geometries. In particular, a wall branch with
$\beta_+<0$ selects decreasing areal radius along its chosen normal. It must
be placed on the correct side of the retained Schwarzschild bridge, with a
complete spacelike slice through the proposed source and cap. A spatial
collar does not specify that continuation.

The next record needs the induced three-metric, its extrinsic curvature
$K_{ij}$, wall embedding and velocity, matter initial data, and both sides'
constraint and junction equations. A globally time-symmetric vacuum slice
cannot simply be relabelled as one containing a moving material wall.
The non-time-symmetric N-body result may become applicable after such a
body is independently constructed. For a delta-function thin shell, either
a suitable theorem for that regularity or a justified resolved source is
still required; arbitrary smoothing would introduce a new model.

If complete compatible initial data agree on a retained open set, local
uniqueness determines agreement inside its domain of dependence. Therefore
the useful finite question is whether the selected wall trajectory and
observer remain in that domain before signals from the modified exterior
can arrive. This conclusion requires the appropriate evolution problem,
including the source; vacuum uniqueness alone does not prove matter-wall
evolution. A Schwarzschild radial null integral can supply a conditional
collar crossing time after the slice, collar edges and clock are specified.
It cannot supply a global sibling delay or formation history.

## Proposed finite admission, before any new simulation

- **Borrow:** the relevant N-body existence theorem and an independently
  sourced local wall solution; import their equations and fixed parameters.
- **Missing connection:** one complete source-supported initial slice whose
  retained vacuum collar matches the selected wall branch, followed by a
  finite domain-of-dependence certificate.
- **Decision and unlock:** determine whether a local spherical trajectory
  can be insulated from a nonspherical parent for a stated observation
  window; otherwise identify the first incompatible data or incoming signal.
- **Domain:** one collar, one wall, one observer, one finite proper-time
  interval. The second collar remains a separate region, without an inferred
  communication link. Choose collar radii only after the Cauchy slice exists.
- **Budget and stop:** a source-and-constraint audit first; no gluing PDE or
  global evolution admitted by this note. A subsequent imported null-ray
  evaluation may use at most 100 sample points and 60 CPU seconds. Stop if
  the initial slice/source is unavailable, its regularity misses the theorem,
  the trajectory leaves the known chart, or a required delay crosses the
  domain boundary. Do not enlarge the collar after seeing an undesired result.
- **Check:** compare null propagation with an independent regular-coordinate
  calculation and verify the proper-time convention. Treat success as a
  scoped connection of imported results, pending focused prior-art comparison.

## Search coverage and unread boundary

The 10 October 2026 search covered arXiv and primary author/institution
copies for Corvino gluing, BL numerical gluing, many-black-hole data, local
KID obstructions, N-body successors, and errata for 0909.1101 and 1004.1355.
The 2002 erratum was found and incorporated; the bounded successor search
did not find an erratum to those two N-body papers. That is not a guarantee
of exhaustive coverage. Their references led to the later Hintz result.
Original IMP proofs, full weighted-space estimates, complete global horizon
analysis, characteristic gluing, and matter-layer regularity theorems remain
outside this reading. No priority claim or proof audit is made.

The existing [reuse rule](../../docs/prior-art-and-open-claim.md) remains the
calculation owner. This note changes the available parent candidate, not the
meaning of the frozen [v2 atlas](nsc-cosmological-atlas-v2.md) or its evidence.
