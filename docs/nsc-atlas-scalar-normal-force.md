# Normal-force compatibility of the conformal probe

The same conformal scalar that crosses the existing wall adds no quadratic
regular normal-stress jump. This completes a local compatibility check left
open by the [surface and tangential-flux calculation](nsc-atlas-scalar-wall-consistency.md).
It gives a necessary condition on a later sourced metric and wall: the
quadratic correction to the mean extrinsic-curvature trace must vanish.
It does not require the wall displacement to vanish.

The [admission](prior-art-and-open-claim.md#scalar-normal-force-admission-13-october-2026)
was recorded before the reduction. The mathematics is borrowed from the
nonminimal Hilbert stress, scalar junction law and Gauss/Israel identities.
This finite application uses no new waveform or pressure curve and makes
no claim to a new general junction theorem. The
[record](../results/development/nsc-atlas-scalar-normal-force-v1/normal-force.json)
and [manifest](../results/development/nsc-atlas-scalar-normal-force-v1/manifest.json)
bind the derivation and exact checks to the original scalar-wall evidence.

## Fixed model and traces

Keep the original canonical real massless scalar, $\xi=1/6$, zero scalar
background and potential, scalar-independent vacuum terms, constant material
tension and no additional surface interaction. The source is the original
$C^3$ compact preparation. The retained wall has $R/\ell\in[0.45,0.55]$,
$H\ell=1$, $k\ell=2$, $M/\ell=0.4$ and the same no-return condition.
The causally excluded right-parent observer supplies no new data.

Use signature $(-+++)$, a common spacelike normal $n$ from the de Sitter
cap $(-)$ to the Schwarzschild side $(+)$, and
$K_{ab}=\tfrac12\mathcal L_n\gamma_{ab}$. Brackets mean plus trace minus
minus trace, and bars are arithmetic means. With $\kappa_4=8\pi G$,

$$
k=\frac{\kappa_4\sigma}{2},\qquad
[K_{ab}]=-k\gamma_{ab},\qquad [K]=-3k.
$$

The common field trace and its intrinsic derivatives are continuous. Put
$p_\pm=(n\phi)_\pm$, $v_a=D_a\phi$ and $d=[p]$. The borrowed scalar law
is $d=-2\xi\phi[K]=6\xi k\phi$, hence $d=k\phi$ at the selected coupling.
Sufficient one-sided regularity is required for the following smooth traces.
No product of a delta distribution with a discontinuous normal derivative
is used. All scalar stresses below are the formal quadratic terms on the
reference geometry, with $\phi=\varepsilon\Phi$.

## Project the complete stress before taking its jump

The original action gives, with its Einstein-tensor term included once,

$$
T^\phi_{\mu\nu}
=\nabla_\mu\phi\nabla_\nu\phi
-\frac12 g_{\mu\nu}(\nabla\phi)^2
+\xi\left[G_{\mu\nu}\phi^2
+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2\right].
$$

In Gaussian normal distance $\eta$, the metric is
$d\eta^2+\gamma_{ab}dy^ady^b$. For a scalar $F$,
$\Box F=\partial_\eta^2F+K\partial_\eta F+D^2F$.
The second normal derivative therefore cancels in the normal projection:

$$
T^\phi_{nn}=\frac12p^2-\frac12v_av^a
+\xi\left(\phi^2G_{nn}+2K\phi p+D^2\phi^2\right).
$$

The two intrinsic terms have common traces. Using
$[p^2]=2\bar p\,d$ and $[Kp]=\bar Kd+[K]\bar p$ gives

$$
[T^\phi_{nn}]=\bar p\,d+\xi\phi^2[G_{nn}]
+2\xi\phi\left(\bar Kd+[K]\bar p\right)
=\xi\phi^2\left([G_{nn}]+12\xi k\bar K\right).
$$

Thus the jump needs no separately evaluated tangential second derivative,
even though either individual normal stress generally does. The stored
waveform samples and geometric interpolants are not differentiated.

## Geometric cancellation and an essential vacuum contribution

The normal Gauss identity, in these conventions, is

$$
G_{nn}=\frac12\left(K^2-K_{ab}K^{ab}-{}^{(3)}\mathcal R\right).
$$

The intrinsic scalar curvature is common. The existing curvature jumps
then give the kinematic relation

$$
[G_{nn}]=\bar K[K]-\bar K^{ab}[K_{ab}]=-2k\bar K.
$$

It is an identity of the reference geometry, not an assumption that the
unknown sourced Einstein solution has already been constructed. Consequently

$$
[T^\phi_{nn}]=2\xi(6\xi-1)k\bar K\phi^2=0,
\qquad \xi=\frac16.
$$

This cancellation holds for every sufficiently regular field trace obeying
the stated matching conditions, not just for a specially fitted waveform.
For symbolic $\xi$, the expression equals the already established scalar
surface coefficient $2\xi(6\xi-1)k\phi^2$ contracted with $\bar K$.
That comparison is a conservation control; it does not create an admitted
nonconformal solution with the same material wall and geometry.

In the actual reference regions, $G^+_{nn}=0$ and $G^-_{nn}=-3H^2$.
Hence $[G_{nn}]=3H^2$ and $\bar K=-3H^2/(2k)$. One can also check this
using the signed radial normals of the existing moving-wall solution:

$$
\beta_- = \frac{M}{kR^2}+\frac{k^2-H^2}{2k}R,\qquad
\beta_+ = \frac{M}{kR^2}-\frac{k^2+H^2}{2k}R.
$$

On the moving branch, $K_\pm=d\beta_\pm/dR+2\beta_\pm/R$, with its
regular continuation through a zero of $\beta_+$. This reproduces both
$[K]=-3k$ and the mean above without a new trajectory integration.
Omitting the cap's Einstein tensor would falsely leave
$[T^\phi_{nn}]=-H^2\phi^2/2$ at conformal coupling. The vacuum term
therefore cannot be dropped by treating both sides as Ricci-flat.

## Force sign and the remaining geometric condition

With our common normal, the material-shell identity has a plus sign:

$$
S^{ab}\bar K_{ab}=[T_{nn}].
$$

There are two independent sign checks. First, contract the Israel law
$[K_{ab}]-\gamma_{ab}[K]=-\kappa_4S_{ab}$ with $\bar K^{ab}$ and use
the Gauss jump above; this gives $[G_{nn}]=\kappa_4S^{ab}\bar K_{ab}$.
Second, normal distributional conservation has a bulk-step coefficient
$+[T_{nn}]$ and a shell coefficient $-\bar K_{ab}S^{ab}$, which must sum
to zero. Neither check requires an evaluated waveform.

[Lobo--Crawford, gr-qc/0507063v2](https://arxiv.org/html/gr-qc/0507063v2)
provides the Israel and Gauss equations (4)--(5). Its printed equation (6)
has the opposite force sign; copying that sign alongside (4)--(5) fails
the checks just given under the stated conventions. We use the identity
fixed by those checks. The library records this reading limitation rather
than treating a printed sign as authority over the geometric equations.

For the constant-tension material wall the complete normal condition is

$$
-\sigma\bar K=[T^{\rm vac}_{nn}]+[T^\phi_{nn}].
$$

The fixed vacuum contribution is $3H^2/\kappa_4$. Since the additional
quadratic scalar jump vanishes, a compatible perturbation must obey

$$
\delta\bar K=0\qquad\text{at order }\varepsilon^2.
$$

Here the variation is the mean curvature of the actual wall in the actual
metric, pulled back to the reference wall for comparison. It includes
metric, induced-metric and wall-displacement corrections. It does not
imply $\delta R=0$, zero one-sided pressure, no energy transfer, or a
vanishing bulk gravitational response. The
[initialized source geometry](nsc-atlas-characteristic-initialization.md)
and [initial mass response](nsc-atlas-source-mass.md) already show why the
unbackreacted geometry cannot be asserted to solve the full problem.

## Evidence, controls and stopping point

The [recorder](../scripts/record_nsc_atlas_scalar_normal_force.py) and
[tests](../tests/test_nsc_atlas_scalar_normal_force.py) compare the direct
one-sided projection, product/jump reduction, Gauss/Israel sign, explicit
vacuum reference and generic-coupling conservation check using exact
symbolic arithmetic. Zero field, quadratic amplitude dependence and
length homogeneity are controls, not numerical waveform experiments.
Under simultaneous normal reversal and side exchange, $\bar p$, $\bar K$
and $[G_{nn}]$ negate, so the normal-stress jump negates too; its vanishing
and the force identity remain invariant. Altered jumps, omitted improvement
terms and the wrong force sign are required to fail their corresponding checks.

Creation authenticates the original scalar-wall evidence and its current
bound sources, then runs the new exact checks within 60 process-CPU seconds.
It refuses an existing result directory. The ordinary check authenticates
source, prerequisite and record bytes without loading symbolic algebra;
the explicit algebra check authenticates first and then repeats the
identities. The earlier record's normal-force-open flag remains unchanged
as a statement about that dated calculation.

This is an exact local compatibility result within its stated assumptions,
not a perturbative remainder bound, pressure history, stability statement
or a solved sourced spacetime. The branch stops here as admitted: the
normal jump is automatic, so pressure curves would not close the next gap.
The next physical question is a closed finite initial/boundary problem for
transporting the initialized geometry toward the wall while respecting
$\delta\bar K=0$ and the remaining junction equations. The cap vacuum
term and coupled material action must survive that construction. Other
atlas routes remain available if this one proves incompatible.
