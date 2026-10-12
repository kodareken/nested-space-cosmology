# Gravitational initial data for the prepared null pulse

The [prepared field](nsc-atlas-prepared-response.md) has an independently
specified source, and its [conformal surface stress](nsc-atlas-scalar-wall-consistency.md)
passes the singular-surface and tangential-flux checks. Neither result supplies
gravitational initial data. This successor applies the borrowed null Einstein
constraint to the same source. A nonzero pulse requires a radius correction
when the null normalization and corner geometry are held fixed.

The calculation is a single characteristic constraint at quadratic order in
the field amplitude. It does not evolve the coupled spacetime or the wall.
The [admission](prior-art-and-open-claim.md#characteristic-gravitational-initialization-admission-13-october-2026)
and [primary-source cards](arxiv-studies.md#characteristic-gravity-library)
record its domain, budget, imports and reading boundaries. A generic bulk
boundedness argument was deferred because it would not specify these missing
initial data.

## Borrowed equations and specified corner

In the smooth retained exterior use the same Jordan-frame action, signature,
canonical scalar normalization and coupling $\xi=1/6$. Write

$$
ds^2=-A(U,V)\,dU\,dV+R(U,V)^2d\Omega^2,\qquad \kappa_4=8\pi G.
$$

The spherical null Raychaudhuri equation gives the geometric identity

$$
G_{VV}=-\frac{2}{R}\left(R_{VV}-\frac{A_V}{A}R_V\right).
$$

An--Lim, [2005.04090v1](https://arxiv.org/html/2005.04090v1),
Eqs. (2.1), (2.6)--(2.7), supply this constraint with minimally coupled
matter. The geometric identity is reused; its matter side is replaced by
the already borrowed full conformal stress. The affine-null constraints and
corner data are also described by Madler--Gannouji--Gallo,
[2503.24162v3](https://arxiv.org/html/2503.24162v3),
Eqs. (16), (29a), (86)--(92). Neither paper supplies a theorem for our
moving material wall.

The two initial legs remain $U=U_0$, $V_a\le V\le V_b$ and
$V=V_a$, $U_0\le U\le U_a$. Their corner has $R_0=0.4\ell$;
the corner values of $A,R_U,R_V$ are their original Schwarzschild values.
There is no added corner mass or expansion perturbation. Keep the bottom
geometry at its reference values; its identically zero scalar data satisfy
its intrinsic null constraint. The original source has an open initial
zero interval, so the initial wall corner and cap data stay compatible with
that reference. This supplies no later wall compatibility condition.

Along the source choose $A=A_0$, where the subscript denotes the original
vacuum geometry. Its background constraint implies $R_{0,V}/A_0$ is a
positive constant. Hence

$$
x=R_0(U_0,V)/\ell
$$

is a dimensionless affine parameter on this leg. The gauge choice keeps it
affine for the perturbed source, anchored by $x_a=0.4$,
$R(x_a)=0.4\ell$ and $dR/dx(x_a)=\ell$. The other corner expansion and
the corner lapse are also fixed as above. One may instead hold a radius
coordinate fixed and solve for the lapse; holding both reference functions
fixed would overdetermine the sourced constraint.

## The unchanged preparation in affine coordinates

Use the original function $F(x)=(1-x/0.8)e^{x/0.8}$, so that
$V(x)=F(x)/U_0$. With the same $V_a,V_b$ and their midpoint, put

$$
b(x)=B\left(\frac{4[V(x)-V_{\rm mid}]}{V_b-V_a}\right),\qquad
p(x)=\frac{b(x)}{x},\qquad
B(z)=(1-z^2)^4\quad (|z|<1),
$$

and $B=0$ elsewhere. The compact function is $C^3$. No source is chosen
from the desired gravitational correction. If $\chi=\varepsilon b$ is
prescribed using the perturbed areal radius as well, its leading scalar
trace is unchanged:

$$
\phi=\frac{\varepsilon}{\ell}p+O(\varepsilon^3),\qquad
\alpha=\frac{\kappa_4\varepsilon^2}{\ell^2},\qquad
\frac{R}{\ell}=x+\alpha h(x)+O(\alpha^2).
$$

These are formal perturbative coefficients with the other dimensional
parameters fixed; no bound on the omitted terms is established.
$\alpha$ is dimensionless. The coefficient $h$ is independent of a chosen
finite amplitude. Changing $\ell$ at fixed $\varepsilon$ and $G$ changes
$\alpha$; a comparison at fixed $\alpha$ must scale the field amplitude
accordingly. This is distinct from a change of units.

## Forced radial correction

Primes below mean $d/dx$, not a wall proper-time derivative. The null
component of the full stress at order $\varepsilon^2$ is

$$
T_{xx}=\frac{\varepsilon^2}{\ell^2}
 \left(\frac23p'^2-\frac13pp''\right)+O(\varepsilon^4).
$$

The reference exterior has $G_{xx}=0$; the nonminimal $\xi G_{xx}\phi^2$
term enters only at the next perturbative order. It is not added twice or
discarded from the underlying theory. The null constraint therefore fixes

$$
h''=-\frac{x}{3}p'^2+\frac{x}{6}pp'',\qquad
h(x_a)=h'(x_a)=0.
$$

The initial scalar amplitude and slope vanish. Integration by parts gives
two useful representations:

$$
h'(x)=-\frac12\int_{x_a}^{x}s p'(s)^2\,ds
       +\frac{x p(x)p'(x)}6-\frac{p(x)^2}{12},
$$

$$
h(x)=-\frac12\int_{x_a}^{x}(x-s)s p'(s)^2\,ds
      +\frac{x p(x)^2}{12}-\frac16\int_{x_a}^{x}p(s)^2\,ds.
$$

They agree with direct integration of $h''$, with the same two zero corner
data. After the compact trace ends, the local terms vanish. For this
nonzero source, both $h$ and $h'$ are strictly negative. The radial-slope
correction becomes constant on the rest of this finite source leg. Thus
keeping both $A=A_0$ and $R=R_0$ would violate the sourced constraint.

This sign does not assume pointwise positive conformal stress. The improved
stress may have either sign during the pulse; its boundary terms cancel in
the post-input slope. Nor is $h'$ itself the null expansion: the latter is
$2R_x/R$ and has a separate radius-denominator correction. The reported
quantity is the radial derivative in the declared affine normalization.

## Independent smooth-exterior frame check

Aviles--Maeda--Martinez, [1910.07534v2](https://arxiv.org/html/1910.07534v2),
Eqs. (2.21)--(2.24), supply the field and metric change. In this smooth
exterior with zero other matter, set

$$
\mathcal F=1-\kappa_4\phi^2/6>0,\qquad
\bar g=\mathcal F g,\qquad \bar R=\sqrt{\mathcal F}R,
\qquad \frac{d\psi}{d\phi}=\mathcal F^{-1}.
$$

The equivalent Einstein-frame scalar is minimally coupled. The affine
parameter must change too: $d\bar x=\mathcal F\,dx$, with the same corner
origin and scale. To the retained order,

$$
\bar x=x-\frac{\alpha}{6}\int_{x_a}^{x}p^2\,ds,\qquad
\frac{\bar R}{\ell}=x+\alpha\left(h-\frac{x p^2}{12}\right).
$$

Reexpressing the radius in $\bar x$ gives

$$
h_E=h-\frac{x p^2}{12}+\frac16\int_{x_a}^{x}p^2\,ds
    =-\frac12\int_{x_a}^{x}(x-s)s p'(s)^2\,ds,
\qquad h_E''=-\frac{x}{2}p'^2.
$$

This is precisely the borrowed minimal-scalar affine focusing equation.
Comparing the two frames at the same untransformed parameter would fail
this check. They are descriptions of the same exterior theory, not two
candidate geometries. Extending the comparison through the material wall
or into the de Sitter cap requires transforming those source actions too;
that extension is not used here.

## Finite evaluation and evidence

The creation-only [initialization record](../results/development/nsc-atlas-characteristic-initialization-v1/initialization.json)
and its [manifest](../results/development/nsc-atlas-characteristic-initialization-v1/manifest.json)
bind the unchanged preparation, response and surface-stress evidence to
this owner, the [adapter](../src/recursive_horizons/nsc_atlas_characteristic_initialization.py),
[producer](../scripts/record_nsc_atlas_characteristic_initialization.py) and
[tests](../tests/test_nsc_atlas_characteristic_initialization.py).
The calculation uses analytic source derivatives, not differences of the
saved waveform. It evaluates the direct stress integral and the integrated
expressions above independently, using Gauss quadrature only within the
compact support.

| Quantity on the source leg | Recorded value |
| --- | ---: |
| Last reference radius $x_b$ | 0.4815042082576364 |
| Compact support in $x$ | [0.4226917905945489, 0.4631925883576189] |
| $h(x_b)$ | -6.12283354819419 |
| $h'(x_b)$ | -162.18055362491776 |
| $\int p^2\,dx$ | 0.06157963228952755 |
| $\int x p'^2\,dx$ | 324.3611072498355 |

These are dimensionless coefficients, not finite radius displacements or
fractional mass changes. For example the formal radius shift is
$\alpha\ell h$, and its radial-slope shift is $\alpha\ell h'$.
No nonzero amplitude or finite-amplitude validity range has been selected.
The sign after the input follows from the positive-integral formulas, not
from numerical roundoff or the sampled extrema.

The constraint quadrature orders 32 and 64 differ by at most $2.54\times10^{-14}$
in $h$ and $1.37\times10^{-12}$ in $h'$. Independently changing the geometric
adapter from order 32 to 64 at fixed constraint order 64 changes them by
at most $4.62\times10^{-14}$ and $1.67\times10^{-12}$.
Direct integration and integration by parts differ by at most
$1.16\times10^{-14}$ and $9.38\times10^{-13}$. The affine frame identity
has a maximum numerical difference below $4.45\times10^{-16}$ and an
independent exact symbolic check. These are refinement diagnostics, not
rigorous decimal enclosures or a remainder bound for gravity.

The complete calculation has 129 displayed source samples, 15360 scalar
constraint integrand evaluations and 1440 nodes in the reused geometric
quadratures. The total 16800 also fits the 20000-evaluation cap; the creation
is guarded by the admitted 60-CPU-second limit. It performs no new wave
solve or bulk Einstein evolution. Controls
check the analytic derivative chain, exact constraint and frame identities,
zero input, corner compatibility, compact-support edges, scaling, domain
and theory changes, budget limits, source drift and evidence tampering.

## Knowledge boundary and next connection

The result supplies the intrinsic source-leg constraint correction and
compatible zero-bottom/corner choices. Mixed Einstein equations still have
to determine transverse derivatives and evolution. The moving wall and
interior must eventually obey the coupled constraints and matching laws.
Neither the full coupled characteristic problem nor its solution is
established by this one constraint.

After $b$ ends, $\phi=\phi_V=0$ on the source leg and $T_{VV}=0$ there.
This is not a vacuum-spacetime statement. For example, the original field
equation gives a generally nonzero transverse derivative

$$
\chi_U(U_0,V)=-\varepsilon\int_{V_a}^{V}K(U_0,q)b(q)\,dq.
$$

Other stress components may persist. No Schwarzschild mass change, total
pulse energy, wall displacement, corrected central waveform, formation
history or stability follows from the radial-slope result. The amplitude
has not been fitted or selected to make the geometry correction small.

The next useful question is which transverse constraint or propagation
equation connects these initialized data to an invariant exterior response
or to the wall. Its admission must name the remaining gravitational datum,
finite domain and required field derivatives before computing them. The
original linear waveform remains evidence in its own approximation, and
other atlas attachments remain open.
