# An invariant mass response on the prepared null surface

The [initialized radius](nsc-atlas-characteristic-initialization.md) satisfies
one intrinsic Einstein constraint. This successor closes the transverse
radius equation on that same initial null leg and computes its geometric
Misner--Sharp mass at quadratic order. It retains the fixed preparation,
canonical conformal action, corner data and source-leg lapse. It does not
evolve the wall or the interior. The [admission](prior-art-and-open-claim.md#invariant-source-leg-mass-admission-13-october-2026)
and [literature library](arxiv-studies.md#source-leg-mass-library) specify
the borrowed laws, finite domain and reading boundaries.

## Coordinates, field and invariant being compared

Keep the earlier notation $m=0.4$, $a=2m=0.8$ and $x=R_0/\ell$ on
$U=U_0$. Change only the second null coordinate from $V$ to the monotone
function $x(V)$. On this leg the dimensionless metric coefficient is

$$
ds^2/\ell^2=-C\,dU\,dx+r^2d\Omega^2,\qquad
C=-4a/U_0>0,\qquad r=R/\ell=x+\alpha h+O(\alpha^2).
$$

$C$ stays at its constant reference value on the source leg, as already
prescribed. Off the leg the coefficient is not assumed constant.
The reference transverse derivative obeys

$$
d(x)=\frac{r_{0,U}}{C}=\frac{a/x-1}{4}>0,
\qquad (xd)'=-\frac14.
$$

Both original future radial derivatives are positive in this sector-IV
interval. Exterior-region monotonicity statements that assume a negative
transverse derivative do not apply here.

Use the geometrically defined mass length

$$
\mathfrak{m}=\frac{R}{2}\left(1-g^{ab}\partial_aR\partial_bR\right),
\qquad \frac{\mathfrak{m}_J}{\ell}=m+\alpha\mu_J+O(\alpha^2),
\qquad \alpha=\frac{\kappa_4\varepsilon^2}{\ell^2}.
$$

With $c=1$, dividing this length by $G$ gives a mass unit. The defined
geometric quantity is a scalar under coordinate changes. It is not in
general invariant under a conformal change of metric, nor automatically
an asymptotic energy or the mass of the whole exterior. The reference
mass is constant, so changing the identification of cuts by an
order-$\alpha$ coordinate displacement does not change its first-order
coefficient. The physical preparation and corner are still held fixed.

## The transverse datum is already determined

Write $\phi=(\varepsilon/\ell)p+O(\varepsilon^3)$ and define the normalized
leading transverse scalar derivative

$$
s(x)=\frac{\ell\phi_U}{C\varepsilon}\bigg|_{\text{leading order}}.
$$

The earlier reduced-field equation and zero bottom data give
$\chi_U/\varepsilon=-\int_{V_a}^{V}K(U_0,q)b(q)\,dq$.
Equivalently, the scalar wave equation along the source says

$$
s'=-\frac{s+d p'}{x},\qquad s(0.4)=0,
\qquad
s=-\frac{d p}{x}-\frac{a}{4x}\int_{0.4}^{x}\frac{p(t)}{t^2}\,dt.
$$

This is determined data, not an independently chosen outgoing profile.
No additional gravitational datum is inserted: the corner transverse
radius derivative remains its original value. A nonzero $s$ can persist
after the compact source trace ends.

## Borrowed mass flux and complete frame conversion

In the smooth exterior let $\mathcal F=1-\kappa_4\phi^2/6>0$ and
$\bar g=\mathcal Fg$. The already imported scalar redefinition gives a
minimal Einstein-frame scalar with the same leading derivatives.
An--Lim's mass transport (2.15), restored to canonical real-field
normalization, supplies

$$
\mu_E'=-x^2d\,p'^2,
\qquad
\mu_E(x)=-\frac14\int_{0.4}^{x}t(a-t)p'(t)^2\,dt.
$$

The corner mass coefficient is zero. Leading background coefficients
suffice in this quadratic flux. This sign follows from $0<t<a$ in the
admitted strip. It is compatible with a positive minimal null stress;
the orientation is that of a white-hole sector, not the usual untrapped
exterior. It does not establish radiation escaping to infinity.

For $\bar g=\omega^2g$, Faraoni--Vitagliano's geometric transformation
gives

$$
\bar{\mathfrak{m}}=\omega\mathfrak{m}
-\frac{R^3}{2\omega}(\nabla\omega)^2
-R^2\nabla\omega\cdot\nabla R.
$$

Take $\omega=\sqrt{\mathcal F}$ and include both null derivatives. At the
retained order this fixes the Jordan-frame answer:

$$
\mu_J=\mu_E+\frac{m p^2}{12}
+\frac{x^2p}{3}\left(s+d p'\right).
$$

Dropping the transverse derivative or identifying the two frame masses
during the input would change this result. This frame comparison is
confined to the smooth exterior; crossing the material wall requires
transforming its source action too. The two frames describe the same
smooth-exterior theory and supply an independent consistency check.

## Independent mixed-geometry and improved-stress checks

Write $r_U/C=d+\alpha k+O(\alpha^2)$ and
$I_0(x)=\int_{0.4}^{x}p^2dt$. In the original $(U,x)$ coordinates put
$H=h-xp^2/12$ and $k_E=k-dp^2/12-xps/6$.
The borrowed minimal mixed equation becomes

$$
(xk_E+dH)'=\frac{p^2}{24},\qquad
k=-\frac{dh}{x}+\frac{I_0}{24x}+\frac{dp^2}{6}+\frac{xps}{6}.
$$

Its integration constant is fixed at the original corner. Substitution
into the Jordan geometric mass definition gives independently

$$
\mu_J=\frac{mh}{x}+2x\left(k+dh'\right).
$$

The full conformal stress on the vacuum reference yields

$$
\frac{\ell^2 T_{xx}}{\varepsilon^2}
=\frac{2p'^2-pp''}{3},\qquad
\frac{\ell^2 T_{Ux}}{C\varepsilon^2}
=\frac{p's-ps'}{3}.
$$

The general-matter mixed equation also checks the transverse construction
directly in the Jordan frame:

$$
(xk+dh)'=\frac{x^2}{6}(p's-ps').
$$

Hayward's general-matter spherical mass law, with this normalization,
then gives

$$
\mu_J'=x^2\left[\frac{p's-ps'}{3}
-d\frac{2p'^2-pp''}{3}\right].
$$

Differentiating the frame expression using the scalar transport equation
agrees exactly with this law. The improved mixed stress is retained;
setting it to zero as for a minimally coupled real scalar would fail.

## Post-input consequence and knowledge boundary

After the compact trace, $p=p'=0$ and both derivatives of $\omega$ vanish,
even though $\phi_U$ need not vanish. Therefore

$$
\mu_J=\mu_E=-\frac14\int_{\text{support}}t(a-t)p'(t)^2\,dt<0.
$$

The coefficient is constant along the remaining finite source leg.
Its negative sign means a correction relative to the reference, not
a negative total mass or a negative-energy field.
The equivalent expression from the earlier radius data is
$\mu_J=[h+(a-x)h']/2+I_0/12$. This is a conditional local mass
response to the specified data. It supplies neither a vacuum future nor
a global mass-loss history. No nonzero amplitude or bound on the
perturbative remainder has been selected.

The next physical connection is to transport compatible sourced geometry
toward the moving wall and test its dynamical junction conditions, or
identify which further field and gravitational data that requires. The
existing linear observer response retains its original approximation;
it has not yet been corrected for gravitational backreaction. Other
atlas geometries remain separate candidates.

## Finite evaluation

The [mass record](../results/development/nsc-atlas-source-mass-v1/mass.json)
and its [manifest](../results/development/nsc-atlas-source-mass-v1/manifest.json)
bind the [adapter](../src/recursive_horizons/nsc_atlas_source_mass.py),
[producer](../scripts/record_nsc_atlas_source_mass.py),
[tests](../tests/test_nsc_atlas_source_mass.py), this derivation and the
unaltered initialization evidence. The producer is creation-only. Its
ordinary check authenticates stored bytes without numerical integration;
explicit numerical replay repeats the admitted finite calculation.

| Recorded quantity | Value |
| --- | ---: |
| Last source radius $x_b$ | 0.4815042082576364 |
| Post-input $\mu_J=\mu_E$ | -28.88319705406449 |
| Remaining transverse scalar $s(x_b)$ | -0.07835520871821607 |
| Largest sampled $|\mu_J-\mu_E|$ during input | 4.265459834842238 |

The last difference is a sampled frame comparison, not an extremum bound.
The dimensionful formal correction is
$\delta\mathfrak{m}_J=\alpha\ell\mu_J$; a dimensionless coefficient alone
does not select an amplitude or give a finite physical mass change.

The independent original-$V$ integral and affine-$x$ primitive for $s$
agree to $8.9\times10^{-16}$ at the displayed points. The mixed metric and
frame mass expressions differ by at most $1.1\times10^{-14}$; the analytic
frame derivative and improved-stress mass law differ by at most
$6.9\times10^{-13}$. Ten exact symbolic identities check the geometry,
frame and stress equations, with omitted-term negative controls.

The 32/64 mass-quadrature comparison changes the mass coefficients by
less than $2.5\times10^{-13}$; changing geometric order independently at
fixed mass quadrature changes them by less than $4.6\times10^{-13}$.
These are floating-point diagnostics, not decimal enclosures or bounds on
the gravitational perturbation remainder.

The 129 displayed samples require 47520 counted quadrature nodes: 1440
in the reused geometry, 15360 in the radius constraint, 15360 in shared
affine mass/transverse moments, and 15360 in the independent original-$V$
kernel integral. Multiple integrands evaluated at one shared node count
as one node. The complete creation has a 100000-node and 60-process-CPU-second
guard. It performs no field re-solve or bulk Einstein evolution.
