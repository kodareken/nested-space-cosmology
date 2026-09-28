# Continuous KS residual to local matter error

The numerical state comparison is replaced by a bound only when its
continuous residual norms are enclosed. The implemented arithmetic is exact
rational arithmetic with rational upper square-root enclosures. Missing
inputs are rejected. Solver tolerances and sampled residuals are not accepted
as established continuous bounds.

The general reconstruction/residual energy method is established numerical
analysis; see [Irene Kyza, ESAIM M2AN 45 (2011), 761–778](https://doi.org/10.1051/m2an/2010101),
especially the reconstructed error equation and energy argument in Section 3.
That paper treats a Schrödinger operator and Crank–Nicolson reconstruction.
The application below uses the existing first-order KS Dirac operator and
checks its own derivative commutators. No claim of a new general error method
or of transferring that paper's discretization theorem is made.

## Norm propagation for this operator

The source envelope obeys `X_rho=A X_z+B X`, where
`A=sigma3/a^2` is Hermitian and independent of z, and B is anti-Hermitian.
On a numerical period of length L, integration by parts gives a skew
generator. For a continuous reconstructed field Xhat and its actual PDE
residual `R=Xhat_rho-A Xhat_z-B Xhat`, set `e=X-Xhat`. With absolute d-rho
along the decreasing-rho evolution,

\[
d\|e\|_2\le\|R\|_2|d\rho|,\quad
d\|e_z\|_2\le(\|R_z\|_2+b_1\|e\|_2)|d\rho|,
\]
\[
d\|e_{zz}\|_2\le(\|R_{zz}\|_2+2b_1\|e_z\|_2+b_2\|e\|_2)|d\rho|,
\quad b_j=\|\partial_z^jB\|_\infty.
\]

Norms include all quadrature-weighted spinor columns. If `e_j` bound the
initial norms, `R_j` bound their residual integrals, and `K_j` bound the
integrals of b_j on a slab, valid endpoint and slab-sup bounds are

\[
E_0=e_0+R_0,\quad E_1=e_1+R_1+K_1E_0,
\quad E_2=e_2+R_2+2K_1(e_1+R_1)+(K_1^2+K_2)E_0.
\]

The square term uses the ordered double integral of b1, equal to half its
full squared integral. No exponential in the large source energy enters
this estimate. Source/preparation errors remain initial errors; they are not
reset to zero at each step.

For a periodic Hilbert-valued function, its mean has norm at most
`L^(-1/2)||e||_2`. The zero-mean primitive kernel has L2 norm `sqrt(L/12)`.
Since e_z has zero mean,

\[
\epsilon_F\le L^{-1/2}E_0+\sqrt{L/12}E_1,\qquad
\epsilon_{F_z}\le\sqrt{L/12}E_2+E_{max}\epsilon_F.
\]

The last term comes from the actual carrier derivative
`F_z=exp(-iEz)(X_z-iEX)`. It is not an outgoing-momentum substitution.

## Radius-dependent commutator constants

For the common normal window with support sigma, write Wj,Uj for global
sup bounds of the two axial functions' derivatives and rmin for a positive
radius bound. Since `|d rho|=a |dT|` and `|chi|<=1`,

\[
K_1\le\frac{|\ell|}{r_{min}^2}
 (\sigma^2W_1/2+\sigma^4U_1/24),
\]
\[
K_2\le|\ell|\left[
 \frac{\sigma^2W_2/2+\sigma^4U_2/24}{r_{min}^2}
 +\frac{2}{r_{min}^3}
 (\sigma^3W_1^2/3+\sigma^5W_1U_1/15+\sigma^7U_1^2/252)\right].
\]

These are consequences of the same reciprocal-radius potential. They are
zero for ell=0. The axial derivative sup bounds must be supplied with their
own proof; sampled maxima are insufficient.

## Coherent contraction and local scope

Let f,z bound the reconstructed weighted Frobenius norms of F,Fz, and let
e,ez bound their errors. The exact source has `0<=C_src<=I`; an optional
operator error ec for its numerical representation gives

\[
\epsilon_D=2fe+e^2+f^2e_c,\qquad
\epsilon_J=ze+fe_z+ee_z+fze_c.
\]

For the original signed-family multiplicity mu,

\[
\epsilon_N\le\mu[(m+|\ell|/r_{min})\epsilon_D+
\epsilon_J/a_{min}],\qquad \epsilon_\beta\le\mu\epsilon_J.
\]

This retains the full coherent source matrix. Numerical contraction
roundoff, geometry/coefficient error and spectral quadrature remain separate.
Both current and reference errors must enter a difference unless a proved
correlated estimate is supplied.

The period is a comparison domain. Its continuum solution equals the
physical solution on I only when the histories and upstream fields agree
throughout I's backward characteristic cone and that cone meets no periodic
copy. The owned speed integral and causal padding check establish the
geometric separation for the declared compact profiles; arbitrary callback
agreement at sample nodes does not establish the needed history equality.
Finite Fourier grids themselves have no exact causal-support theorem.

The missing implementation input is a directed enclosure of the reconstructed
PDE residual over each continuous rho/z cell, including spatial truncation,
coefficient evaluation and arithmetic. The error propagation above alone
does not close the physical local gate.
