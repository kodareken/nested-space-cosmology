# One spatial seed and a repeated offset: a constructive geometric map

Douglas's starting proposal is one spatial seed plus an offset operation,
applied repeatedly to generate a circle, a sphere and continued structure.
Here the offset is the primitive update; no independent proper-time
coordinate is used to generate the spatial map. The construction below
realizes part of that proposal with an explicit metric and a proof.

The [calculation](../results/summaries/nsc-atlas-offset-generator-v1.json)
uses one continuous seed coordinate and one irrational offset. A history
rule supplies the two phase values rather than choosing two independent
initial phases. Retaining the earlier draws, choosing a cross-pass distance
law and taking its metric completion are explicit assumptions.

## Preparation and choice

We compared supplying three independent spatial coordinates, claiming a
single smooth instantaneous trajectory fills a volume, and retaining an
offset history whose distance completion has higher dimension. The third
formulates the proposed drawing without making the false instantaneous-rank
claim. The earlier [suspension and golden-map controls](nsc-atlas-golden-vortex.md)
remain valid in their scope; this entry replaces independently supplied
phase coordinates with a single deterministic history.

The most likely mistakes are hiding a second freely supplied phase,
confusing a dense history with its completed space, and using an embedding's
component count as the intrinsic spatial dimension. The formula, completion
proof and induced-metric calculation check those separately. Finite coverage
uses fixed counts 100,200,1000,4000 and at most 120 CPU seconds; no physical
target selects a phase, count or distance rule.

## 1. The same accumulation operation is applied twice

Let alpha be irrational, theta_0=phi_0=0, and set, modulo one,

$$
\theta_{n+1}=\theta_n+\alpha,\qquad
\phi_{n+1}=\phi_n+\theta_n.
$$

Equivalently, Delta theta=alpha and Delta phi=theta. One fixed offset feeds
the first accumulator, whose history feeds the second. Thus

$$
\theta_n=n\alpha\quad\bmod1,\qquad
\phi_n=\frac{n(n-1)}{2}\alpha\quad\bmod1.
$$

The index n counts applications. It is not a separately supplied clock in
seconds. Both phases are determined by that count and alpha; there is no
second independent phase input. The inverse update exists, so the map
alone does not select a thermodynamic arrow or explain why updates occur.

**Density proof.** For any nonzero integer pair (k,l), the corresponding
Fourier character is the polynomial phase

$$
\alpha\left[kn+\frac{l n(n-1)}{2}\right].
$$

If l=0, its average tends to zero by the geometric-sum formula. If l is
nonzero, subtracting the value at n from the value at n+h gives a linear
phase with irrational slope l h alpha for every positive integer h.
Its average again tends to zero. The van der Corput bound then gives an
upper limit of order 1/H for the squared original average for every H;
letting H grow makes it zero. Weyl's criterion proves joint equidistribution
on the two-phase torus. The general polynomial theorem is attributed to
Weyl; [Tao's notes](https://terrytao.wordpress.com/2010/03/28/254b-notes-1-equidistribution-of-polynomial-sequences-in-torii/)
give the character/differencing framework, and
[Yifrach](https://arxiv.org/abs/2201.07138) records the classical statement
before developing a further generalization. No generalization beyond this
quadratic special case is needed here.

A rational offset need not work: alpha=1/4 returns to the initial state
after eight updates and supplies only finitely many directions. A zero
offset leaves the original line. Nonresonant accumulation is an essential
condition, not a word added after the calculation.

## 2. Line, circular boundary, spherical boundary

Put Theta_n=2pi theta_n and Phi_n=2pi phi_n. A first rotation of the retained
seed interval r in [-R,R] gives lines (r cos Theta_n,r sin Theta_n). Their
completion is a disk with a circular boundary.

For the second stage define the direction correlation

$$
C_{nm}=\cos\Theta_n\cos\Theta_m
+\sin\Theta_n\sin\Theta_m\cos(\Phi_n-\Phi_m),
$$

and the distance rule

$$
d^2[(r,n),(s,m)]=r^2+s^2-2rsC_{nm}.
$$

This is a declared metric kernel. Its Gram factorization is

$$
v_n=(\sin\Theta_n\cos\Phi_n,\sin\Theta_n\sin\Phi_n,
\cos\Theta_n),\qquad C_{nm}=v_n\cdot v_m.
$$

The distance is therefore Euclidean distance between r v_n and s v_m;
coincident labels, including every representation of the origin, are
identified. Density of the phase pairs implies density of directions on
the sphere. Every point of the ball is a limit of retained seed points.
Conversely every seed point lies in the ball, which is complete and closed.
**The metric completion is exactly the three-dimensional closed ball.**

This is a conditional constructive theorem: one spatial interval and one
offset history suffice for the stated completion. Finite histories and
the countable union of lines have zero three-dimensional volume before
completion. If motion were restricted to walking along the spokes through
their shared origin, a different distance would result. The distance law,
memory and completion must not be hidden inside the phrase "one plus time."

## 3. The same history supplies the closed parent geometry

For a radial seed r in [0,pi R], choose the positive-curvature distance law

$$
d=R\arccos\left[
\cos(r/R)\cos(s/R)+\sin(r/R)\sin(s/R)C_{nm}\right].
$$

Its completion is the round closed three-dimensional space used by the
parent model. A four-component representation is
(R cos(r/R), R sin(r/R) v_n), with one fixed-norm constraint. The intrinsic
metric, calculated by differentiation, is

$$
d\ell^2=dr^2+R^2\sin^2(r/R)
\left(d\Theta^2+\sin^2\Theta\,d\Phi^2\right).
$$

There are three intrinsic spatial directions. Four embedding components
do not make four independent spatial dimensions. Likewise a separate
five-coordinate auxiliary fluid/gravity model is not the count of primitive
variables in this generator.

The radius is still a parameter of this geometric construction. The
existing constant-field parent at its retained minimum supplies a compatible
choice through its Hamiltonian constraint. For zero extrinsic curvature,

$$
{}^{(3)}R=\frac{6}{R^2}=2\Lambda_*,\qquad
R=\sqrt{\frac{3}{\Lambda_*}}=3.684522144189
$$

in the current model units, with Lambda_star=.220983024213. Momentum
constraints vanish for those constant fields. This joins the generated
spatial geometry to existing admissible parent initial data. It is not a
primordial formation calculation or an SI scale selection.
The earlier action's phase ratio and radius are compatibility inputs;
this construction has not derived that action or its auxiliary
internal-dimensional interpretation from the offset alone.

The scalar Laplacian then has eigenvalues l(l+2)/R^2. The first two positive
values are .220983024213 and .589288064567. Under radius rescaling,
R^2(-Delta) preserves its dimensionless spectrum. This is an explicit
normalized operator inheritance law across separately specified regional
metrics; it does not yet attach those regions physically.

## 4. Finite resolution and measure

For N complete radial lines, let h be the spherical covering radius of
their two endpoint directions. The largest distance from the unit ball
to that line family is exactly sin(h). For N oriented north-to-south
meridians in the closed parent, it is min(h_oriented,pi/2) times R.
Spherical Voronoi vertices and a nearest-neighbor check evaluate these
finite errors; all phases are calculated at 80-digit precision before
converting the direction coordinates.

| Offset | Passes | Ball distance / R | Closed-parent distance / R |
| --- | --- | --- | --- |
| Retained action ratio | 100 | .35147195 | .50145491 |
| Retained action ratio | 1000 | .16691160 | .19978947 |
| Retained action ratio | 4000 | .08015006 | .12788791 |
| Golden-angle fraction | 100 | .39247396 | .57017792 |
| Golden-angle fraction | 1000 | .15218010 | .22441514 |
| Golden-angle fraction | 4000 | .08859494 | .13519849 |

The retained action offset is alpha=.0573586470504441, from the previously
derived ratio of linear mode frequencies. The golden comparison uses
alpha=(3-sqrt(5))/2=.381966011250105. Both are irrational; neither was fitted
to these errors. Their physical identification as the fundamental offset
is still a hypothesis. This construction does not select a cosmic pass
count, the 10^60 hierarchy, or a golden branching average.

Uniform phase passes are also not uniform physical volume. Their direction
second moments tend to (1/4,1/4,1/2). The metric area weight |sin Theta|
instead gives (1/3,1/3,1/3). The closed-parent volume element is
R^2 sin^2(r/R)|sin Theta| dr dTheta dPhi; the ball has r^2 in place of
R^2 sin^2(r/R). This is the metric measure used by the covariant action,
not an identification of update frequency with energy or particle density.

## What this adds to the atlas

The generator is now explicit and does not require an independent second
starting phase. Its completed spatial geometry and a compatible parent
slice are proved under the stated construction rules. Existing spacetime
calculations can be read as tests on a generated geometry, rather than
being asserted as the primitive origin of it.

What remains is to derive the physical selection of this distance kernel,
the offset law, two accumulation levels followed by nesting, and a duration
or irreversible arrow. Further accumulation levels can be constructed, so
exactly three is not selected by irrationality alone. A coherent physical
theory must also supply regional stress, formation and observable predictions.
