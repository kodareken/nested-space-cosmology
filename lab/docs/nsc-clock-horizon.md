# Clock inheritance and the exterior–interior observer map

The existing smooth black-universe geometry realizes a finite exterior horizon,
a trapped region, a minimum sphere and an expanding child region in one metric.
It also gives a precise example of why unbounded coordinate time does not prove
infinite affine runtime. This is a fixed-background audit, extending the
[Dirac tetrad calculation](nsc-dirac-tetrad.md); it does not derive the geometry's
quantum source or an ancestral clock map.

## An invariant observer map

Use signature `+---`, units `c=L_throat=1`, and the existing background

\[
ds^2=d\tau^2-(d\rho+\beta d\tau)^2-r^2d\Omega_2^2,\qquad
r=\sqrt{1+\rho^2},\quad \beta=\sqrt{1-A},
\]
\[
A=1+3\rho+3(1+\rho^2)(\arctan\rho-\pi/2).
\]

The sphere area is \(4\pi r^2\). The two future null directions normal to
these spheres and their expansions are

\[
\ell_\pm=\partial_\tau+(-\beta\pm1)\partial_\rho,\qquad
\boxed{\theta_\pm=\frac{\ell_\pm(4\pi r^2)}{4\pi r^2}
=\frac{2\rho}{1+\rho^2}(-\beta\pm1).}
\]

Multiplication of a future null normal by a positive function preserves each
expansion's sign. Thus the sign map is independent of that normalization:

| Domain | θ+ | θ− | Spherical geometry |
|---|---:|---:|---|
| ρ > ρh | + | − | Exterior untrapped spheres |
| ρ = ρh | 0 | − | Outer horizon |
| 0 < ρ < ρh | − | − | Future-trapped spheres |
| ρ = 0 | 0 | 0 | Minimum areal radius |
| ρ < 0 | + | + | Expanding, anti-trapped spheres |

This sign classification has an analytic basis. Set
\(x=\pi/2-\arctan\rho\). Then \(A'=6(1-\rho x)>0\): it is immediate
for \(\rho\le0\), and for \(\rho>0\) follows from
\(x=\arctan(1/\rho)<1/\rho\). Together with \(A(0)=1-3\pi/2<0\)
and \(A(+\infty)=1\), this gives one positive horizon. Positivity of
\(\beta^2\) follows from \(\arctan y>y/(1+y^2)\) for \(y>0\);
both arctangent bounds follow by integrating \(1/(1+t^2)\) on \([0,y]\).
The runner checks the exact expansion identities and independent numerical
probes in every sign region.

The horizon is at \(\rho_h\simeq1.90069160547\), with invariant areal
radius \(r_h\simeq2.14770309380\). A distant instrument can fail to
resolve it. A coordinate change cannot turn its nonzero area into zero. A
particle interpretation additionally requires mass, charges, spin/statistics
and a measured scattering response; unresolved appearance alone establishes
none of those identities.

A horizon's spatial cross-section is a **two-dimensional sphere**. The local
spatial domain has **three dimensions**. A particular observer's surrounding
cosmological horizon is a causal construction using that observer's worldline;
it need not be the parent event horizon. Spherical symmetry about the parent's
center is not a proof of isotropy about every child observer.

## Interior expansion and its limit

For \(A<0\), use the original static coordinate \(t\) as a spatial coordinate
and let \(dT=-d\rho/\sqrt{-A}\), so decreasing \(\rho\) is future-directed.
The interior metric is

\[
ds^2=dT^2-a_\parallel(T)^2dt^2-r(T)^2d\Omega_2^2,\qquad
a_\parallel=\sqrt{-A}.
\]

These homogeneous spatial sections have topology \(\mathbb R\times S^2\).
The two independent expansion rates and shear are

\[
H_\parallel=\frac{A'}{2\sqrt{-A}},\qquad
H_\perp=-\sqrt{-A}\frac{\rho}{1+\rho^2},\qquad
\sigma^2=\frac{(H_\parallel-H_\perp)^2}{3}.
\]

At the neck \(H_\perp=0\) and the shear is nonzero. As \(\rho\to-\infty\),
\(A/\rho^2\to-3\pi\), \(A'/\rho\to-6\pi\), and both expansion rates
approach \(\sqrt{3\pi}\), with vanishing shear. This is a Kantowski–Sachs
interior approaching local de Sitter behavior, with \(\Lambda=9\pi\) in the
declared units. It is not an exactly isotropic FLRW geometry at finite radius,
nor proof that its three-dimensional spatial volume is compact.

This exterior-to-expanding-interior geometry is established prior art:
[Bronnikov and Fabris](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.251101)
and [Bronnikov, Dehnen and Melnikov](https://arxiv.org/abs/gr-qc/0611022)
construct regular black-universe families. Their phantom source is not adopted
here. Deriving the source from the common Dirac operator remains the project's
physical task.

## A finite affine past despite infinite coordinate time

The two-dimensional radial metric has determinant \(-1\), including at the
horizon. Its coefficients and their derivatives are smooth there; the horizon
is not a curvature singularity. Direct Christoffel calculation gives

\[
\Gamma^\tau_{\tau\tau}=\frac{\beta A'}2,\qquad
\Gamma^\rho_{\tau\tau}=\frac{A A'}2.
\]

On \(\rho=\rho_h\), \(k=\partial_\tau\) is a future null horizon tangent,
with

\[
\nabla_k k=\kappa_h k,\qquad
\kappa_h=\frac{A'(\rho_h)}2\simeq0.238325799634.
\]

It is a geodesic tangent with a non-affine parameter. An affine parameter is

\[
\boxed{\lambda=\lambda_0+C e^{\kappa_h\tau},\qquad C>0.}
\]

The runner checks both \(\lambda''-\kappa_h\lambda'=0\) and the inverse
affine geodesic equation. For \(C=1\), the entire interval
\(-\infty<\tau\le0\) has affine length one. Thus the **declared PG patch**
has a past-incomplete null generator, even though this generator exists at
every finite coordinate time.

This is compatible with the earlier complete characteristic flows and Cauchy
Dirac evolution: coordinate-time completeness is a different claim from
affine geodesic completeness. An extension beyond this PG patch remains to be
constructed and tested. The finite affine endpoint does not establish a
curvature singularity or a beginning of the complete spacetime. It supplies a
concrete boundary where the proposed chart-extension mechanism can be tested.

## An operational ancestry clock

Let room \(m=0\) be ours, with successive ancestors \(m=1,2,\ldots\).
Choose a specified spectral clock in each room and derive its boundary
transport. For positive factors constant over the intervals being compared,

\[
q_m=\frac{d\widehat\tau_{m-1}}{d\widehat\tau_m},\qquad
T_0=\sum_{m=0}^{\infty}\left(\prod_{j=1}^{m}q_j\right)
\Delta\widehat\tau_m.
\]

For varying factors, integrate the pointwise composite clock map. Infinite
ancestral duration requires this positive sum to diverge. Infinitely many
rooms alone is insufficient: unit intervals with every \(q_m=1/2\) total
two room-0 ticks; every \(q_m=1\) or every \(q_m=2\) gives a divergent
sum. These are exact counterexamples to inference from the number of rooms,
not proposed inheritance factors. The current operator has not selected them.
Even a divergent assigned-clock sum must be distinguished from actual proper
or affine completeness of trajectories in a consistently glued spacetime.

For ruler and clock maps \(d\ell_c=\alpha d\ell_p\) and
\(d\tau_c=\beta_{\rm clock}d\tau_p\), the dimensionless parent speed in
child units is

\[
v=\frac{\alpha}{\beta_{\rm clock}}\frac{c_p}{c_c}.
\]

Pure coordinate rescaling of the same null cone gives \(v=1\). A physical
speed difference needs a dimensionless cone or transmission measurement from
the operator. In matched units the characteristic crossing-duration ratio is
\((L_p/L_c)/(c_p/c_c)\); both ratios can exceed one while their quotient
is smaller than, equal to, or larger than one. Larger and faster does not
determine the clock ordering.

Elapsed proper time from a specified formation event along a specified local
observer remains a meaningful local age. A common global age requires a
defined causal gluing and clock convention. Scale-dependent clocks do not by
themselves rule it out.

The original weak-gradient parent thesis remains the constructive working
direction. An entropy ordering is unresolved until matter entropy,
gravitational/boundary entropy and coarse-graining are defined together. This
audit does not edit the protected creative authority or interchange its
minimum-entropy wording with a maximum-entropy assertion.

## Reproduction and claim boundary

Run `python3 scripts/check_nsc_clock_horizon.py --check`. Exact expressions,
sign labels, numerical probes, and script/document/upstream-source hashes in
[the record](../results/nsc-5-clock-horizon.json) are all compared. Floating
fields use the existing comparator's `rtol=1e-8, atol=1e-8`; expressions and
metadata match exactly. `--output PATH` creates a new file exclusively.

The result verifies the fixed-background observer map and exhibits a specific
affine boundary. It derives neither the quantum source, a physical ancestral
clock map, the global extension, nor a particle identity. The general rule is
to test curvature and tidal invariants, proper/affine completeness and
extendibility; a chart change cannot remove an invariant curvature divergence.
