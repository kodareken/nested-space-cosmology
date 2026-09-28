# Exact incoming covariance matching at a spectral pair

For one unchanged angular/compact block, the linearized fixed-C0 matching
kernel is exactly

\[
\boxed{\delta\widehat C(E_o,E_i)
=\delta A(E_o,E_i)C_i f_i^\dagger
+f_o C_o\,\delta A(E_i,E_o)^\dagger.}
\]

The corresponding canonical CAR identity is

\[
\boxed{\delta A(E_o,E_i)P_i f_i^\dagger
+f_o P_o\,\delta A(E_i,E_o)^\dagger=0.}
\]

These close a missing NSC connection: a prepared incoming response can be
tested at an individual spectral pair without replacing the full covariance
by a finite weighted energy sum. The first-order background delta functions
perform that energy integration exactly. This is a linear matching test,
not a full action variation, nonlinear matching proof or stationary history.

## Reuse and exact normalization

Reuse the [common KS restriction](nsc-common-ks-trace.md), the
[retarded preparation family](nsc-retarded-compatible-preparation.md), and
the existing full source fibers. No state law, field equation or geometry
parameter changes. At the translation-invariant reference surface, write

\[
F_0(z,E)=f(E)e^{-iEz},\qquad
C_0(z,z')=\int_{\mathbb R}\frac{dE}{2\pi}
F_0(z,E)C_{\rm src}(E)F_0(z',E)^\dagger.
\]

Here `Ci=C_src(Ei)`, `Pi=P_src(Ei)`, and similarly for the output label.
The integral is over signed real source energy, with the existing `k=-E`
convention. Define the transforms without prefactors:

\[
\delta A(E_o,E_i)=\int_{\mathbb R}dz\,e^{iE_o z}\delta F(z,E_i),
\]
\[
\delta\widehat C(E_o,E_i)=\int_{\mathbb R^2}dz\,dz'\,
e^{iE_o z}\delta C(z,z')e^{-iE_i z'}.
\]

For the fixed retarded source law, `delta C_src=0`. In each differentiated
term, the unvaried background plane wave gives `2*pi*delta(E-Ei)` or
`2*pi*delta(E-Eo)`, cancelling the source measure. No factor of `2*pi`,
inherited energy weight, positive-energy folding or group degeneracy is
added to the displayed pair formula. As a normalization check, the
unvaried double transform is
`2*pi*delta(Eo-Ei)*f(Ei)*C_src(Ei)*f(Ei)^dagger`.

Only the source covariance's energy diagonality is used; its matrix within
each source fiber need not be diagonal. If the source covariance or source
projector itself varies, its differentiated contribution must be added.
It is absent here because the entire earlier preparation is unchanged.

## Actual trace, analytic background and source fibers

Use the actual, unweighted rho1 tangent trace. The
[restriction owner](../src/recursive_horizons/nsc_common_ks_trace.py) and
[prepared-history adapter](nsc-compatible-prepared-history.md) fix its
normal/half-density map. With `z=tau+S1`, let
`g(E)=R_Sigma*Phi_reference(E,rho1)` and define the Fourier transform of the
time-dependent tangent in tau. Then

\[
f(E)=e^{iES_1}g(E),\qquad
\delta A_z(E_o,E_i)=e^{iE_oS_1}\delta A_\tau(E_o,E_i).
\]

The entire pair kernel therefore differs from its tau representation by
the common phase `exp(i*(Eo-Ei)*S1)`. Norms are unchanged. Do not reapply a
stationary phase to an already time-dependent trace, and do not use the
drifting unperturbed SBP trace as the translation-invariant background.
Use the authenticated stationary mode amplitude and its analytic phase.

The background transform is over the whole line, analytically. Only the
compact tangent trace is integrated numerically over the finite recorded
window. A double transform of the old finite two-position matrix would
introduce window sinc factors instead of delta functions. Although
`delta F` has compact incoming support, `delta C(z,z')` need not have compact
support in both coordinates. The old four-energy weighted position kernel
and its projector diagnostic do not implement the identities above.

Keep all three source columns and the full coherent `C_src` and `P_src`.
For the inherited group14 below-threshold controls, the infinity column is
closed and the coherent horizon block remains active. The
[mode owner](../src/recursive_horizons/nsc_pg_massive_modes.py) fixes
`P(E)=diag(1,1,0)` below `abs(E)=m`; the incoming label of each leg determines
its projector. Thus `delta A(Eo,Ei)=delta A(Eo,Ei)*P(Ei)`, even when the
output Fourier label is different. Work away from zero and threshold
energies unless their separate limiting/distributional treatment is supplied.

Do not reduce the source covariance to two local moments or drop
source-complement correlations. In an open fiber,
`Q=P-f^dagger*f` can be nonzero: replacing `delta A*C*f^dagger` by a
closed two-by-two contraction can lose `delta A*Q*C*f^dagger`.
The full-column formula retains this term automatically. The fixed radial
geometry and pure-radius variation do not mix angular blocks, so an
individual block is a legitimate operator test; other orthogonal blocks
cannot cancel its nonzero kernel.

## Exact causal coverage of the recorded trace

For the chosen compact variation, the axial support is
`z in [S1+0.09,S1+0.21]`, and the normal window has outer width0.03.
The pure-radius variation leaves the principal PG speeds
`d rho/d tau=+/-1-beta0` unchanged. Both point toward decreasing rho.
Using `S'=beta0/a0^2` gives `d z/d rho=+/-1/a0^2` along the two
characteristics. Only the upstream portion `1<=rho<=rho_plus` can affect
the incoming surface. Therefore

\[
\operatorname{supp}_{\tau}\delta F_\Sigma
\subset[0.09-X,0.21+X],\qquad
X=\int_1^{\rho_+}\frac{d\rho}{a_0(\rho)^2}.
\]

An exact coarse bound suffices; no numerical root or optical quadrature is
needed. The [fixed profile](../src/recursive_horizons/nsc_lorentzian.py) gives

\[
a_0^2=3[(1+\rho^2)\arctan(1/\rho)-\rho]-1,\qquad
(a_0^2)'=6[\rho\arctan(1/\rho)-1]<0\quad(\rho>0).
\]

Since `a0(1)<1`, the normal width
`integral_1^rho_plus d rho/a0=3/100` implies `rho_plus<103/100`.
At `rho=103/100`, use
`atan(100/103)=pi/4-atan(3/203)>pi/4-3/203` and `pi>314159/100000`.
Exact rational arithmetic gives

\[
a_0(103/100)^2>
\frac{547699824079}{812000000000}
=\frac{16}{25}+\frac{28019824079}{812000000000}
>\frac{16}{25}.
\]

Thus `a0>4/5` throughout that upstream interval. With normal coordinate
`u=T(1)-T(rho)`,

\[
X=\int_0^{3/100}\frac{du}{a_0}<\frac{3}{80},\qquad
\operatorname{supp}_{\tau}\delta F_\Sigma
\subset(21/400,99/400)=(0.0525,0.2475).
\]

The recorded `[0,0.3]` trace therefore covers the exact retarded response,
independently of the downstream half of the radius variation. This is a
continuum support statement, not a certificate for numerical SBP tails,
time quadrature or the accuracy of the computed tangent.

## What the pair test decides

The CAR equation follows from complete canonical resolution on the fixed
intrinsic/normal surface. It is an independent normalization and response
consistency check. Never subtract its numerical residual from the covariance
response or enforce it by projecting/normalizing the sampled fields.

Both terms of the covariance pair are essential: a nonzero field tangent
can cancel from the state tangent. Covariance tangents need not be positive;
the pair satisfies `delta Chat(Eo,Ei)=delta Chat(Ei,Eo)^dagger`. Positivity
belongs to the underlying covariance family, not to its derivative.

For a regular continuous pair kernel away from thresholds, a rigorously
resolved nonzero covariance pair disproves fixed-C0 linear matching for
this particular variation and block. A zero pair does not prove matching
at all pairs. Numerical resolution differences and a small CAR residual
are indicators, not error bounds. The
[response-resolution owner](nsc-retarded-response-resolution.md) retains
that numerical scope; this analytic reduction does not upgrade its evidence.

Reuse decision: apply these identities to saved tangent traces and exact
stationary source fibers; no old field, scattering, source or coefficient
producer is required for the reduction. Stop at a correctly normalized
pair test and its honest error account. No new state law, full4D Hadamard
gate, action term, nonlinear matching result, physical initial tuple or
metric timestep follows from it.
