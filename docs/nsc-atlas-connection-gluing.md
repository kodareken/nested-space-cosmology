# Attaching verified shell bodies through their vacuum ends

**Status, 10 October 2026:** a conditional initial-data assembly using an
imported existence theorem. The [connection slice](nsc-atlas-connection-slice.md)
supplies the body to check: regular de Sitter cap, moving thin wall,
Schwarzschild bridge and an exact right asymptotically Euclidean (AE) end.
The construction below can retain that body while replacing its distant
exterior. It establishes neither a formation history nor a global evolution
theorem for ideal matter shells.

This resolves a narrower question than the earlier
[gluing review](nsc-cosmological-atlas.md). A direct application of a smooth
global-body theorem would miss the wall's regularity. Instead, apply the
theorem whose inputs are smooth **ends**, and use locality to retain an
already verified, piecewise-smooth body. The sealed earlier records retain
their original scope.

## The imported end theorem and its actual inputs

[Chrusciel, Corvino and Isenberg, 1004.1355v1](https://arxiv.org/html/1004.1355v1),
equation (1) and Theorem 1, take smooth vacuum AE ends $(E_i,h_i,K_i)$ with
timelike ADM four-momenta. Their required decay is

$$
\partial^\alpha(h-\delta)=O(r^{-1-|\alpha|}),\qquad
\partial^\beta K=O(r^{-2-|\beta|}),\qquad
|\alpha|\leq\ell+1,\quad |\beta|\leq\ell,\quad\ell\geq2.
$$

Prescribed precompact open neighbourhoods of the end boundaries are retained
isometrically. For sufficiently large separation the output is vacuum, has
one AE exterior, and agrees with Kerr near infinity. Its ADM four-momentum
approaches the sum of the input four-momenta. The separation threshold is an
existence statement, with no numerical value supplied here. Section 3's
construction acts outside the retained neighbourhoods. Neither a trivial
KID kernel for the exact Schwarzschild end nor input Regge-Teitelboim parity
is an extra hypothesis of Theorem 1.

For this application choose $E_i$ wholly inside the right Schwarzschild
vacuum end, beyond the shell and bridge. In asymptotically Cartesian
coordinates its time-symmetric data have the stated decay and
$(m_i,\mathbf p_i)=(M_i,0)$ with $M_i>0$. Its inner boundary and a surrounding
open overlap are ordinary exterior spheres, away from the horizon. This
specifies smooth inputs without changing the physical wall.

## Why the shell does not need to be smoothed for this assembly

The following is a locality argument for the constraint equations, not a
new gluing PDE result. Let each original body have a positive-definite,
continuous, locally Lipschitz three-metric $h$, smooth on each side of its
closed shell, and bounded piecewise-smooth extrinsic curvature $K$. Require
the bulk constraints and the shell jump constraints before assembly. In
particular, $K$ is not itself a delta distribution.

Use signature $(-,+,+,+)$, future slice normal $N$, and
$K=-\mathcal L_Nh/2$. Define $E=T(N,N)$ and $J_i=-T(N,e_i)$. With $G=c=1$,

$$
\mathcal H=R[h]+(\operatorname{tr}K)^2-K_{ij}K^{ij}-16\pi E=0,
\qquad
\mathcal M_i=D_j(K^j{}_i-\delta^j{}_i\operatorname{tr}K)-8\pi J_i=0.
$$

Near a shell use Gaussian proper distance $s$, increasing from the cap side
to the exterior side, so $h=ds^2+q_{AB}dy^A dy^B$. Write
$A_{AB}=\partial_s q_{AB}/2$ and $H_S=q^{AB}A_{AB}$, and let brackets mean
exterior minus cap. The singular coefficients are

$$
R[h]_{\rm sing}=-2[H_S]\delta_S,\qquad
\bigl(D_j\pi^j{}_i\bigr)_{\rm sing}=[\pi^s{}_i]\delta_S,
\qquad \pi^j{}_i=K^j{}_i-\delta^j{}_i\operatorname{tr}K.
$$

The first identity follows from the normal derivative of mean curvature in
[Miao, math-ph/0212025v1](https://arxiv.org/html/math-ph/0212025v1), section 2,
equation (3). The second follows by differentiating a step in $\pi$; the
connection coefficients are bounded and add no delta. Products of bounded
$K$ add no delta either. Thus the required shell conditions are

$$
[H_S]=-8\pi E_S,\qquad [\pi^s{}_i]=8\pi J_{S i}.
$$

This is consistent with the continuous-metric distributional junction
framework of [Mars and Senovilla, gr-qc/0201054v1](https://arxiv.org/html/gr-qc/0201054v1),
sections 4-5, equations (54), (67), (70)-(71) and (79). Agreement of a boundary
metric alone would not establish these conditions.

For a pure tension wall of rest density $\sigma$, let its velocity and normal
at the slice be $U=\gamma(N+v s)$ and $n=\gamma(vN+s)$. Restricting the
spacetime wall delta to the slice supplies a factor $1/\gamma$. Projecting
$S_{\mu\nu}=-\sigma(g_{\mu\nu}-n_\mu n_\nu)$ therefore gives
$E_S=\sigma\gamma$ and $J_{S s}=\sigma\gamma v$, with zero angular momentum
density. In spherical proper distance, $H_S=2R'/R$ and
$K^\theta{}_\theta=-V/R$, where $V=N(R)$. The jump tests reduce to

$$
[R']=-4\pi R\sigma\gamma,\qquad [V]=4\pi R\sigma\gamma v.
$$

These are exactly the independent distributional tests required of the
connection-slice adapter. A moving wall need not move relative to every
possible slice; conversely a nonzero $V$ does not require matter in a vacuum
region. The complete $(h,K)$ and the declared observer frame determine the
constraint test.

After those body tests pass, the assembly has three steps:

1. Select a smooth vacuum overlap inside each end, and use Theorem 1 to retain
   it as an exact copy of **both** $h$ and $K$ in the common exterior.
2. Attach the unchanged inner body through that overlap, using its isometry.
   Cut only in the overlap; the body shell and its neighbourhood are untouched.
3. Check the constraints locally. Near each shell the old distributions are
   unchanged. In the connector the constraints vanish by the imported theorem.
   Across the overlap the data agree on an open set, so there is no new jump
   and no artificial surface source.

A weak test function can be partitioned over this open cover: every local
constraint residual is zero, hence the global residual is zero. This partitions
test functions, **not metrics**. An arbitrary cutoff interpolation of metrics
would generally violate the constraints. No low-regularity extension of the
gluing PDE theorem has been assumed, and no arbitrary smoothing of the wall
has introduced a new source model.

The reattached remainder need not be compact. Locality also retains another
unmodified end of an original body; Theorem 1's surrounding discussion
explicitly allows multiple original ends. Only the selected ends are combined
into the designated common exterior.

## Unequal bodies in the same theory

For two shell-bearing bodies, use the same interior $H$, tension
$k=4\pi\sigma$, vacuum exterior and wall equation of state. Different
positive Schwarzschild masses $M_A$ and $M_B$ are different initial states
in that theory. They need separate slice and shell checks. Changing $M$ is
not a change in a coupling; changing $H$ or $k$ between a failed and successful
test would be a different test.

The first explicit body uses $H\ell=1$, $k\ell=2$, $M_A/\ell=0.4$ and the
wall event $R_w/\ell=1.5$. A second mass may be admitted only with its own
recorded spacelike slice, regular cap and bridge, exact AE end and jump
residuals. The existence theorem imposes no equality of the masses. It does
not let the atlas prescribe the final exterior ADM mass exactly as their
sum, select a particular numerical separation, or infer sibling signalling.

The smaller concrete assembly uses only **one** shell-bearing body. Choose
the checked body A with $M_A/\ell=0.4$, and a standard complete two-ended,
time-symmetric Schwarzschild bridge B with $M_B/\ell=0.6$. Use the right AE
end of each as a theorem input. Retain A through $R_{\rm cut}/\ell=4$ and an
open overlap beyond it; retain B's chosen compact collar and attach its whole
unchanged remainder, including its other AE end. B is a vacuum body of the
same gravitational theory and does not introduce a second value of $H$ or
$k$. It supplies no second wall, de Sitter cap or completed child attachment.

This manifold has **one designated common parent end and B's other AE end**,
not one total end. It has no unexplained inner boundary at B. Unequal masses
here compare two body states, not two already constructed cosmological
children. The relative placements remain an existence family.

| Map item | Permitted status | Evidence required |
| --- | --- | --- |
| One cap-wall-bridge-end body | Explicit initial-data connection when its adapter checks pass | Complete slice, bulk and shell constraints, orientation and centre/throat checks |
| That body plus a complete unequal Schwarzschild bridge B | Conditional shared initial-data existence family; B's other end remains | A's passing body record, two exact vacuum end overlaps; no second child claimed |
| Two such bodies in a common parent | Conditional existence family at sufficiently large separation | Two passing body records; smooth vacuum end and exact overlap certificates |
| Particular common-parent transition metric | Unresolved | An imported construction evaluated only if a later observable needs it |
| Future wall motion after exterior replacement | Unresolved beyond an established source evolution domain | Appropriate Einstein-matter-shell evolution and a domain containing that worldtube |
| Right exterior observer for a finite interval | Vacuum domain-of-dependence certificate below | Exact retained right exterior slice, fixed observer and outer boundary |

## A vacuum observer window that does not assume shell evolution

The source evolution question can be avoided for a deliberately narrower
observer. Suppose the retained slice contains the whole right exterior
$2M<R<R_{\rm cut}$ at future-oriented Schwarzschild time $T=0$, with $K=0$.
Take a static observer at $2M<R_o<R_{\rm cut}$. The right chart's future $T$
must not be confused with the oppositely oriented Killing coordinate used on
the wall's left exterior patch.

Import the tortoise/null coordinates from
[Kuchar, gr-qc/9403003v1](https://arxiv.org/html/gr-qc/9403003v1), section 2.1,
equations (11)-(12), in unscaled units:

$$
r_*(R)=R+2M\log\left(\frac R{2M}-1\right),\qquad
u=T-r_*,\quad v=T+r_*.
$$

The radial metric is $f(-dT^2+dr_*^2)$, $f=1-2M/R$. At $(T,R_o)$ the two
radial past null rays meet the initial slice at optical radii
$r_*(R_o)-T$ and $r_*(R_o)+T$. Since $r_*(2M)=-\infty$, the inner endpoint
lies strictly in the retained right exterior for every finite $T$. The
outer endpoint stays there precisely while

$$
0\leq T<T_{\rm safe}:=r_*(R_{\rm cut})-r_*(R_o),\qquad
0\leq\tau_o<\sqrt{f(R_o)}\,T_{\rm safe}.
$$

Nonradial causal curves cannot reach an initial radial edge faster. The
strict bound identifies an open vacuum domain whose smooth initial data are
unchanged. Smooth vacuum local existence and uniqueness are the relevant
input; see [Sbierski, 1309.7591v2](https://arxiv.org/html/1309.7591v2), section 2,
Theorem 2.6 and the stated smooth setting. No assertion about a global
Einstein-shell development is necessary to specify this local vacuum region.
Any compatible development containing this vacuum domain agrees there up
to isometry.

If only a finite annulus $R_{\rm in}<R<R_{\rm cut}$ is retained, with
$R_{\rm in}>2M$, the bound becomes the **minimum of both optical distances**
from the observer to its two edges. If the slice is not $T=0$, its actual
null-coordinate endpoints must be used. The formula above cannot be reused
on an arbitrary slice or inside a horizon.

Losing the certified domain at its outer edge is not a prediction that a
disturbance arrives then. This is not a parent-child travel time, a sibling
message, a cross-universe age, or a certificate for the future wall itself.
The wall would require the source's evolution and its own causal domain.

## Formation is a separate missing connection

The signed thin-wall branches and regular-looking slices have substantial
prior art. [Aguirre and Johnson, gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2),
section II, equations (12)-(20) and Figures 9-10, distinguish retained
regions by the signed wall factors. Their discussion of the Farhi-Guth
obstruction relates anti-trapped surfaces, the null energy condition and a
noncompact Cauchy surface to past null incompleteness. Footnote 6 explicitly
distinguishes incompleteness from a curvature singularity. Their exceptional
classically buildable inflating examples use nearly equal interior/exterior
vacuum energies; this is not a certificate for the present Schwarzschild
exterior. No formation result or perturbation stability result is imported.

[Ansoldi, Guendelman and Shilon, 0711.2198v1](https://arxiv.org/html/0711.2198v1),
sections 2-3.3, likewise separates the familiar spherical matching equations
from singularity, tunnelling, stability and mass-threshold questions. Its
survey points back to Blau-Guendelman-Guth and Farhi-Guth-Guven. Those original
papers' full proofs have not been independently audited here. A regular cap
on one initial slice does not establish a classically realizable history,
past completeness, generic stability, or nonsingular cosmology.

## Finite admission and remaining reading boundary

This note admits no gluing PDE solve or parameter sweep. It borrows the end
theorem, distributional constraint identities and Schwarzschild null
propagation. The finite decision is whether the fixed body data satisfy the
end and overlap conditions and whether a chosen observer stays in the
unchanged vacuum domain. That unlocks an initial-data attachment and a local
observer clock without pretending to solve the parent globally.

The connection-slice record owns its bounded body calculation. Any additional
observer evaluation uses fixed radii and at most 100 samples and 60 CPU
seconds. Stop if the data fail a shell constraint, lose the AE vacuum
overlap, require changing the source, or put the observer outside the stated
domain. Compare the null bound with regular Kruskal coordinates. A failure
leaves the already valid body or alternative parent branches visible.

The focused 10 October 2026 search covered N-body gluing and its references,
distributional hypersurface curvature, spherical shell initial-data
projections, smooth vacuum uniqueness, and false-vacuum formation challenges.
Reading was equation/theorem level for the sections cited above. The complete
weighted-space gluing proof, general low-regularity Einstein-shell
well-posedness, quantitative separation thresholds, global horizons of the
assembled parent, full bubble perturbation analysis and original historical
formation proofs remain outside this audit. This bounded comparison makes
no claim of priority for the construction.

The existing [reuse and admission rule](prior-art-and-open-claim.md)
remains the research owner. Imported formulas, their implementation and this
scoped assembly argument are not new-to-world gluing or thin-wall dynamics.
