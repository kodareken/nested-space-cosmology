# Conditional sign of one exact diagonal preparation response

For the selected positive energy `E=0.5562120090641313` in group14_1,
the following implication concerns the **full coherent source law**:

\[
\boxed{\operatorname{Re}(P_{v,\Sigma})_{01}>\frac14
\quad\Longrightarrow\quad
\delta\widehat C(E,E)_{00}
<-\frac{69}{875000000}\approx-7.885714\times10^{-8}.}
\]

This is a conditional analytic certificate. The missing input is a directed
bound on the exact selected-energy affine-vacuum projector in the stated
incoming canonical basis. A stored finite-frame mode or a tolerance
comparison does not provide that input. No ODE or source producer is run
for the argument below.

## Reuse and the source uncertainty being bounded

Reuse the [exact Fourier pair equations](nsc-incoming-fourier-matching.md),
the [short-shell radial response](nsc-retarded-radial-response.md), and the
unchanged compact positive radius pulse. The source state, occupations,
normalization and action remain unchanged. The construction isolates a
scattering-phase-independent dominant projector while retaining the entire
coherent remainder as a bound. It does not set that remainder to zero.

For positive E below the fixed mass threshold, the
[physical sewing owner](../src/recursive_horizons/nsc_pg_massive_modes.py)
gives

\[
S=\begin{pmatrix}0&1&0\\ R&0&0\end{pmatrix},\quad |R|=1,
\qquad P_v=\operatorname{diag}(0,1,0).
\]

Thus `S Pv S†=diag(1,0)`: the dominant column is the first horizon
fundamental column, independent of the Jost reflection phase. The exact
[source covariance](../src/recursive_horizons/nsc_paired_horizon_preparation.py)
has horizon block

\[
C_{\rm src}=\begin{pmatrix}f_H&-is\\ is&1-f_H\end{pmatrix},
\qquad s^2=f_H(1-f_H),
\qquad f_H=(1+e^{2\pi E/\kappa_{\rm src}})^{-1}.
\]

The third column remains closed. On the active block,
`(C_src-Pv)^2=f_H I`, so its operator norm is exactly `sqrt(f_H)`.
Canonical coisometry transports this to
`norm(C_Sigma-Pv_Sigma)=sqrt(f_H)`. All horizon coherence and relative
phases remain in this exact difference. No numerical column is normalized.

At the diagonal Fourier pair, `delta A=K F` and the exact CAR identity
gives `K†=-K`. Consequently

\[
\delta\widehat C(E,E)=[K,C_\Sigma],\qquad
\|[K,C_\Sigma]-[K,P_{v,\Sigma}]\|_2
\le2\|K\|_2\sqrt{f_H}.
\]

The background symbol `C_Sigma` is distinct from the distributional
unvaried double-Fourier kernel; no `delta(0)` is evaluated here. No energy
quadrature weight, extra2pi or angular multiplicity enters this pair.

## Orientation and short-shell generator bound

In the owned unitary rho frame,

\[
G_E=\frac{i}{a}H_E,\qquad
H_E=-m\sigma_1+\frac\ell r\sigma_2-\frac E a\sigma_3,
\qquad r^2=1+\rho^2.
\]

Let `s=T(rho)-T(1)<0` upstream, and let `U(1,rho)` be the exact unitary
reference transport to the incoming slice. The radius direction has U=0
in the compatible ansatz and `delta r=s chi(s) w(z)`. The
[radial response equation](../src/recursive_horizons/nsc_retarded_radial_response.py)
has forcing `-i ell s chi(s) w_hat(0) sigma2/(a r²)` and is integrated
from upstream rho to1. Reversing those integration limits therefore gives

\[
K=-i\int_1^{\rho_+}b(\rho)U(1,\rho)\sigma_2U(1,\rho)^\dagger\,d\rho,
\qquad b=\frac{\ell(-s)\chi(s)\widehat w(0)}{a r^2}\ge0,
\quad I=\int_1^{\rho_+}b\,d\rho>0.
\]

The weight vanishes at endpoints and is strictly positive on an interior
subinterval. K is independent of preparation and reflection because this
is the full local Duhamel insertion, not merely an inference from CAR.
The unchanged past/inflow and fixed incoming restriction supply no extra
preparation or frame derivative.

The [causal-cover proof](nsc-incoming-fourier-matching.md) establishes
`rho_plus<103/100` and `a>4/5`. These remain conservative for the actual
binary outer radius0.03, which is smaller than `3/100`. The selected binary
inputs satisfy `m<8/5`, `0<E<3/5`, and `ell²<128/25`; since `r²>=2`,
`ell/r<8/5`. Hence

\[
\|G_E\|_2^2
<\frac{2(8/5)^2+(3/4)^2}{(4/5)^2}
=\frac{2273}{256}<9.
\]

Unitary conjugation gives
`norm(U sigma2 U†-sigma2)<=2 integral norm(G_E) d rho<.18`.
With `K0=-i I sigma2`, it follows that

\[
\|K-K_0\|_2<\frac9{50}I,\qquad \|K\|_2\le I,
\qquad \|[K-K_0,P_{v,\Sigma}]\|_2<\frac9{25}I.
\]

In the fixed Pauli convention `sigma2=[[0,-i],[i,0]]`,

\[
[K_0,P_{v,\Sigma}]_{00}=-2I\operatorname{Re}(P_{v,\Sigma})_{01}.
\]

This fixes the response sign; an oppositely oriented radial integral would
give the wrong sign.

## Coherent-source margin and positive insertion weight

The actual frozen source inputs obey `E>55/100` and
`0<kappa_src<6/25`. With `pi>314/100`,

\[
2\pi E/\kappa_{\rm src}>\frac{1727}{120}>14.
\]

Since `e>27/10` and `(27/10)^14>10^6`, the exact occupation satisfies
`f_H<10^-6`. The full coherent correction is therefore less than `.002 I`.
Under the stated projector condition,

\[
\delta\widehat C(E,E)_{00}
<-\tfrac12I+\tfrac9{25}I+\tfrac1{500}I
=-\tfrac{69}{500}I.
\]

For an explicit lower bound on I, use normal coordinate `u=-s`. Its
orientation is `d rho=a du` (not a derivative of the areal radius).
On `0<=u<=3/500`, the inherited binary inner radius0.007 makes `chi=1`.
Also `r²<(103/100)²+1<21/10` and the positive angular label obeys `ell>2`.
The normalized [axial bump](../src/recursive_horizons/nsc_transmitting_resolvent.py)
is `exp(1-1/(1-x²))`; on `abs(x)<=1/2` it is at least `exp(-1/3)>2/3`.
Its binary halfwidth0.06 is greater than `1/20`, so
`w_hat(0)>1/30`. No Fourier quadrature is required for that lower bound.
Therefore

\[
I>\frac{20}{21}\frac1{30}\int_0^{3/500}u\,du
=\frac1{1750000},
\]

which yields the displayed response bound. All constants are conservative
for the stored binary controls; no parameter was rounded into a more
favorable physical value.

## Remaining certificate and interpretation

The required projector certificate must enclose the exact affine-vacuum
`Re(Pv_Sigma[0,1])` at this selected energy, in this canonical Pauli basis,
including its initialization, transport and arithmetic uncertainties. The
group32 high-energy certificate does not supply it. A precise Jost phase is
unnecessary for this sufficient condition because every coherent phase is
already covered by the remainder above.

Stop this analytic reduction at the conditional inequality and its single
missing projector bound. Numerical short-shell agreement is useful evidence
but is not substituted for that bound. Once supplied, the result would
exclude fixed-C0 linear matching for this specified compact variation and
angular block. It would not exclude other histories, preparations, surface
extensions or NSC families, nor establish full action stationarity or a
metric solution. The full research goal remains OPEN. No new source state,
coupling, Hadamard gate, action term or metric timestep is introduced.
