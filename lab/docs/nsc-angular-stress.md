# Full angular reference stress at finite child radius

The finite-radius canonical massless reference tensor is now evaluated on
the actual child geometry. Its covariant normalization, unitary evolution,
local conformal map and energy/work identity have independent controls.
Angular, momentum and time resolution have been varied. The reported
numerical sensitivities are measured convergence evidence, not rigorous
infinite-tail bounds.

The target is the canonical massless reference tensor specified in
[the reference construction](nsc-massless-reference.md), on the actual
child geometry. The physical parent-state correction, massive compact
sector and finite-Lambda causal conversion remain separate required
contributions. No reference width is changed to obtain a stress sign.

## Reuse the direct-product geometry

Set q=arctan(u), where u=-1/rho_coordinate. Then

\[
\bar g=g/r^2=\frac{dq^2}{W(q)}-W(q)dz^2-d\Omega_2^2,
\quad r=\csc q,
\quad W=3(\pi-q)+3\sin q\cos q-\sin^2q.
\]

Conformal proper time obeys d eta/dq=-1/sqrt(W). The metric is a direct
product M2 x S2, even when M2 depends on time. The published S2 Dirac
spectrum therefore reduces the massless four-dimensional field to
two-dimensional Dirac channels with m=kappa/R, kappa=1,2,... and
both-sign degeneracy 4 kappa. This channel mass is angular momentum;
it is not a four-dimensional rest mass.

The 2D canonical Hamiltonian is h=m sigma1+p sigma3, p=k/b,
b=sqrt(W). The specified auxiliary ultrastatic ground state is evolved
back into the actual child. We use the instantaneous eigenbasis to
avoid subtracting large sea energies numerically. If its spinor is
(B,A), then its excess energy is 2 omega |B|^2, omega=sqrt(m^2+p^2).
The pressure also retains the off-diagonal coherence; replacing it by
particle occupations alone would lose part of the source.

The native C++ transport and Python transport implement the same
fourth-order two-exponential unitary method. A small direct comparison
agrees to rounding accuracy. The native path accelerates unresolved modes;
it changes neither the reference nor the evolution equation.

## Covariant subtraction and the angular sum

The Dirac adiabatic expansion is a unitary spinor expansion, not the
scalar WKB ansatz. See
[Barbero–Ferreiro–Navarro-Salas–Villaseñor, section IV](https://arxiv.org/html/1805.05107).
For H=dot(b)/b and constant channel mass, the second-order subtractions are

\[
E_2=\frac{m^2p^2H^2}{8\omega^5},\qquad
P_2=\frac{m^2p^2(H^2-\dot H)}{4\omega^5}
       -\frac{5m^2p^4H^2}{8\omega^7}.
\]

The implementation obtains the fourth-order terms by the same normalized
spinor recursion. Their continuous momentum integrals agree with the
two-dimensional Dirac heat coefficient

\[
a_4^{(2)}=-R_2^2/120-\Box R_2/60,
\qquad
\Gamma_{2,m}^{(4)}=-\frac{1}{960\pi m^2}\int\sqrt{|g_2|}R_2^2
\]

for compactly supported bulk metric variations. In q coordinates
R2=-W'', and the R2-squared variational tensor has components

\[
H^{(1)q}{}_q=\tfrac12(W'')^2-W'W''',\qquad
H^{(1)z}{}_z=H^{(1)q}{}_q-2W W''''.
\]

Consequently int dp E4/(2 pi)=H1_q/(480 pi m^2), while the pressure
integral is -H1_z/(480 pi m^2). These were checked by independent
momentum quadrature before the large angular calculation.

Subtract orders zero, two and four before the angular sum. The fourth
order is a convergence subtraction for that sum, and must be restored
with the four-dimensional normalization. For the joint zeta function,
the relevant factor is z R^(2z+2) zeta(1+2z). Taking its derivative
before setting z=0, and converting zeta normalization to the existing
proper-time subtraction, gives the finite harmonic coefficient

\[
L_R=\ln(\mu R)+\gamma_E/2.
\]

It cannot be replaced by an arbitrary value for a divergent sum. At
R=1, the restored density and longitudinal pressure are respectively
L_R H1_q/(480 pi^2) and -L_R H1_z/(480 pi^2).

The static part of that same four-dimensional proper-time convention is

\[
\rho_{\rm cyl}(R;\mu)=
-\frac{2\zeta'(-3)+[\ln(\mu^2R^2)+1-\gamma_E]/120}
       {4\pi^2R^4}.
\]

This follows from the full R2 x S2 spectral zeta function, with the
4-component spin trace counted once. The existing zero-winding proper-time
owner approaches it on cutoff removal. Its static pressure and anomaly
also agree. The integrated R2 term is topological for the base-metric
variations; its role in the sphere response is retained by the exact
four-dimensional trace relation.

After restoring these pieces, the two independent base-metric responses
give bar rho and bar p_z. Spherical symmetry and the massless trace give
bar p_sphere=(bar rho-bar p_z-anomaly)/2, where

\[
\mathcal A_4=\frac{1}{16\pi^2}
\left[-\frac{(R_2-2)^2}{60}-\frac{11R_2}{90}
                         -\frac{\Box R_2}{30}\right],
\quad \Box R_2=-W W''''-W'W'''.
\]

Longitudinal momentum cutoffs used after subtraction are numerical
integration boundaries. They are not replacements for the physical
spacetime proper-time regulator.

## Transport the computed tensor to the physical child

The conformal factor is sigma=-ln(sin q). The implementation integrates
the known anomaly on g_t=exp(2t sigma) bar g, keeping the two base metric
functions N(q),a(q) independent during variation. The local functional is

\[
F_{\rm WZ}=\frac{1}{16\pi^2}\int_0^1dt\int\sqrt{|g_t|}
\,\sigma\,[\alpha C_t^2+\beta E_t+\gamma\Box_tR_t],
\]

with alpha=-1/20, beta=11/360, gamma=-1/30. Compact support of the
source variation allows integration by parts of the box-R term. No
physical boundary stress is inferred by dropping an endpoint variation.
The construction retains the nonzero Weyl contribution; sigma is the
specified geometry, not a new matter field.

The two bulk variations of F_WZ are added to the computed barred rho,p_z
and multiplied by exp(-4 sigma). The physical sphere pressure is fixed
by the physical four-dimensional anomaly. This local map approaches the
previously recorded massless-reference asymptotic tensor. Finite-difference
checks of the density derivative satisfy both the barred and physical
energy/work equations, including at the neck.

## Finite-radius result and numerical sensitivity

The specified smooth auxiliary continuation creates high-frequency
coherences. Their oscillatory pressure makes hard numerical boundaries
inefficient. The computation retains the original reference and uses
smooth windows on the already-subtracted, convergent remainder. The
counterterms remain unchanged. Windows are fixed in comoving labels and
are moved outward to assess the limit; they are not a physical regulator
or a change of state.

The selected calculation retains 2048 angular levels and longitudinal
momentum through 8192, with windows equal to one through angular level
1536 and momentum 6144 before tapering. Lower angular, momentum and time
resolutions are recorded. Already computed mode bands are reused when
extending the sum.

| Areal radius r | rho | p_parallel | p_sphere | rho+p_parallel |
|---|---:|---:|---:|---:|
| 10.0167 | 0.1022194 | -0.1027418 | -0.1028342 | -0.0005224 |
| 5.0335 | 0.0998904 | -0.1014998 | -0.1018347 | -0.0016094 |
| 2.0858 | 0.0988723 | -0.0904788 | -0.0888546 | 0.0083935 |
| 1.1884 | 0.2382098 | 0.0554015 | -0.0380241 | 0.2936113 |
| 1, the neck | 0.7101684 | 0.5026775 | 0.0376558 | 1.2128459 |

Units are hbar=c=L_throat=1 and mu=1. These are the specified auxiliary
reference values, not measurements or an identification of the physical
parent state. At the neck the resolved-window spreads are approximately
0.000755 in rho, 0.000482 in p_parallel, 0.000506 in p_sphere and
0.000990 in the radial null combination. They are sensitivity ranges,
not statistical error bars. The momentum-quadrature probe changes the
physical components by at most about 4.3e-8; the time-step probe changes
them by at most about 1.1e-4.

The independent density-derivative test gives barred and physical Ward
residuals below 1e-9 on its stated probe. Mode unitarity is preserved to
about 2e-13 in the large runs. No previously recorded sphere, compact or
geometry generator is rerun. The fast full replay requires a C++17 compiler;
the Python transport is retained as an independent implementation control.

## Bind the reference to the common source budget

After extracting one healthy Einstein term from the common functional,
write a_EH=A_EH L_throat^2>0. The existing neck geometry requires

\[
L_{\rm throat}^4 T_{\rm total,kk}
=4a_{\rm EH}(1-3\pi/2)=-14.8495559\,a_{\rm EH}<0.
\]

The calculated reference gives +1.2128459. It therefore cannot supply
the neck's null requirement by itself. No Einstein coefficient was fitted.
The required remaining null contribution is
-14.8495559 a_EH-1.2128459 in these units. This is a residual requirement,
not a value inserted into an unevaluated source. A pure vacuum term
proportional to g_mu_nu cannot change this null component.

The binding target remains one functional:

\[
\Gamma_{\rm one}^{\rm CTP}
=\Gamma_{\rm ref,ren}^{\rm CTP}
+\big(\Gamma_{H,\nu}^{\rm CTP}+\mathcal C_{\nu,\mu}^{\rm CTP}\big)
+\Delta\Gamma_{\rm state}^{\rm CTP}
+\Gamma_{\rm remaining}^{\rm CTP},\qquad
\frac{\delta\Gamma_{\rm one}^{\rm CTP}}{\delta g^{\mu\nu}}=0.
\]

This calculation supplies the first reference variation. The state term
must match the actual parent covariance and flux; it is not a fitted
fluid. The compact sector already contained in Gamma_H is counted there
once, and its causal evaluation is still needed. Remaining interactions,
boundary and compensator terms must come from the same action. The
Einstein term must not be counted again on the source side after extraction.

This is a completed reference-source calculation and a conditional
source-budget obstruction. It is not a full finite-Lambda tensor,
parent-selected state, self-sourced solution or cosmological prediction.

Evidence: [immutable result](../results/development/angular-stress.json),
reproduced by scripts/check_nsc_angular_stress.py --check. All source and
control fields use the declared comparison policy. Raw development outputs
and integrands are preserved in the private runs store. Two bounded Grok
investigations returned no usable reports and are not counted as proof;
the normalization and subtraction were checked directly against the
product heat/zeta identities, independent quadrature and the work equation.
