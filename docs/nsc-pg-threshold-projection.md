# Resolving the massive threshold in the PG covariance

The new PG mode fields determine the source map. The remaining numerical
problem in the first compact packet projection is concentrated just below
the parent mass threshold. For group 13, comparing 16 and 32 Gauss nodes per
existing energy interval gives

| Energy interval | Gram difference, operator norm | Centered covariance difference |
|---|---:|---:|
| $[0,0.25]$ | $2.6\times10^{-15}$ | $9.8\times10^{-15}$ |
| $[0.25,1]$ | $4.8\times10^{-12}$ | $1.9\times10^{-12}$ |
| $[1,\pi/2]$ | $6.28\times10^{-2}$ | $1.85\times10^{-2}$ |
| $[\pi/2,4]$ | $9.0\times10^{-12}$ | $5.2\times10^{-12}$ |
| $[16,40]$ | $1.33\times10^{-4}$ | $1.31\times10^{-4}$ |

The other finite intervals are retained as resolved controls at their stated
accuracy. This comparison identifies the numerical owner of the large error;
it is not a new state or a refit of a physical parameter.

## Exact subgap split

Let the finite-horizon basis be projected onto the four existing real radial
packet profiles, with both spin components retained. Arrange its columns as
eight-component vectors $u,v,w$: $u$ is the exterior outgoing-horizon part,
$v$ is the common ingoing part, and $w$ is the interior horizon partner.
Below $E=m$, the physical source map has the form

$$
F_J(E)=\begin{pmatrix}u+R(E)v&w&0\end{pmatrix},
\qquad |R(E)|=1.
$$

The third column is absent because parent infinity has no propagating channel
there. The horizon covariance remains
$C_H=\left(\begin{smallmatrix}f&-is\\is&1-f\end{smallmatrix}\right)$.
With $d=f-1/2$, the centered covariance integrand is

$$
K(E)=K_0(E)+R(E)K_1(E)+R(E)^*K_1(E)^\dagger,
$$

$$
K_0=d(uu^\dagger+vv^\dagger-ww^\dagger)
+is(wu^\dagger-uw^\dagger),
\qquad K_1=dvu^\dagger-isvw^\dagger.
$$

The Gram integrand has the same split with
$G_0=uu^\dagger+vv^\dagger+ww^\dagger$ and $G_1=vu^\dagger$.
These identities reconstruct the physical packet projections from the saved
mode fields; no algebraic completion of an observed covariance block is used.

## Analytic frequency integration

The finite-horizon basis is analytic in frequency inside its Frobenius strip.
For complex $z$, a product is continued as

$$
u(E)v(E)^\dagger\longmapsto
u(z)v(\bar z)^\dagger.
$$

Using $u(z)v(z)^\dagger$ would destroy analyticity. The source functions
$f(z)$ and $s(z)$ are analytic for $|\Im z|<\kappa_h/2$.
Only the term $R(z)K_1(z)$ is moved off the real energy axis. For a subgap
interval $[a,m]$, Cauchy's theorem gives

$$
\int_a^m R(E)K_1(E)\,dE
=\int_{a\to a+ih\to m+ih\to m}R(z)K_1(z)\,dz,
\qquad 0<h<\kappa_h/2.
$$

Adding the adjoint recovers the conjugate term in the real physical integral.
The contour height is a numerical integration parameter; it specifies no
transport time, state, metric history or regulator. The source remains the
same Lorentzian horizon/infinity preparation.

## Why the threshold corner does not add a source

The required endpoint bound follows from the owned exterior Dirac current.
Choose a fixed regular exterior cut only for this argument. For the retarded
Jost solution and $\Im z>0$,

$$
\frac{d}{dx}(\psi^\dagger\sigma_3\psi)
=-2\Im z\,\psi^\dagger\psi,
\qquad
j(x_b)=2\Im z\int_{x_b}^{\infty}|\psi|^2\,dx>0.
$$

Thus the ratio $r_b=\psi_2/\psi_1$ lies in the unit disk. Let the analytic
horizon basis at that cut be
$F_h=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)$.
Its horizon reflection coefficient is

$$
R(z)=\frac{-c(z)+a(z)r_b(z)}{d(z)-b(z)r_b(z)}.
$$

At real $z=m$, current normalization gives $|d|^2-|b|^2=1$. Hence the
denominator is bounded away from zero uniformly for $|r_b|\le1$ in a
sufficiently small upper-half-plane neighborhood of $m$. The finite-horizon
basis and the packet coefficients are bounded there, so $R(z)K_1(z)$ is
bounded and the omitted small corner arc vanishes. This bound does not
require a pointwise real-axis limit of the rapidly winding phase.

Inside the contour, a pole corresponding to a Jost solution outgoing at
infinity and ingoing at the horizon would give a non-real eigenvalue after
continuation into the child. The already owned self-adjoint whole-line
operator excludes it. The Frobenius/source poles are avoided by the chosen
strip. A cut in this proof is an evaluation surface, not a new reflecting
boundary or an added term in the action.

## Numerical gate

The real split must reconstruct the saved physical field projection. The
complex extension must reproduce the existing real-frequency owner and use
the conjugate-energy dual. The physical subgap integral is then checked under
both quadrature refinement and a contour-height change within the same
analytic strip. Neither a uniform phase average nor a chosen covariance tail
is an acceptable substitute.

This calculation addresses the finite subgap panel. The high-energy packet
tail remains a separate part of the complete covariance integral. A complete
spatial $C_{\mathrm{PG}}$, transmitting endpoint jets and metric evolution are
not assigned by the contour identity alone.
