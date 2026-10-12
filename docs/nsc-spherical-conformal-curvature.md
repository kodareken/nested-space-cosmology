# Actual invariant curvature of the conformal trajectory

The saved conformal $T=0.05$ to $0.3$ trajectory now has a metric-curvature
measurement in
[the curvature record](../results/development/nsc-spherical-conformal-curvature-v1.json).
This is postprocessing of the sealed Cauchy frames. No geometric state is
re-solved or evolved. The calculation uses the actual projected ODE rate;
$\chi+2$ is compared afterward and is not substituted for metric curvature.

## The realized time jet

The conformal gauge has $L=Q$, $\beta=0$, so
$h=Q^2(dt^2-dx^2)$. Let $A_g$ denote prolongation and $P_g$ the existing
geometry pullback. The actual finite rate is

$$
\dot Q_g=P_g\left[\frac{Q_f p_{\chi,f}}{2F_\chi}\right],
\qquad Q_f=A_gQ_g.
$$

Because $F_\chi$ is constant and nonzero, its exact directional Jacobian
along the actual ODE rate gives

$$
\ddot Q_f=A_gP_g\left[
\frac{\dot Q_f p_{\chi,f}+Q_f\dot p_{\chi,f}}{2F_\chi}
\right],
\qquad
\dot Q_f=A_g\dot Q_g,\quad
\dot p_{\chi,f}=A_g\dot p_{\chi,g}.
$$

Both projections are retained. $\dot p_{\chi,g}$ is read from the actual
projected rate, rather than replacing the metric time jet by an unprojected
Euler equation. Independent directional differences of `rates(...).Q`
validate the Jacobian. A manufactured top-band example also detects the
wrong result obtained by omitting the final projection.

The lapse derivative is $\dot L=\dot Q_f$, and $\dot\beta=0$.
The existing direct metric helper constructs $R_h$ from these jets, with the
owned periodic spatial derivative. In this gauge the corresponding formula
is

$$
R_h=\frac{2}{Q^2}\left[
\frac{Q_{xx}}Q-\frac{Q_x^2}{Q^2}
-\frac{\ddot Q}Q+\frac{\dot Q^2}{Q^2}
\right].
$$

An independent Christoffel/Ricci contraction of the two-metric jets checks
the sign and factors. The spatial FFT helper and the dense differentiation
matrix represent the same trigonometric derivative in exact arithmetic;
their finite precision difference is recorded separately as a conditioning
indicator. It is not silently called a physical curvature difference.

## Invariant and measured domain

The four-dimensional spherical Weyl invariant follows the established
conformal product identity from
[the owned action chart](nsc-spherical-feedback-action.md):

$$
C^2=\frac{(R_h-2)^2}{3r^4}.
$$

Here $R_h$ comes from the metric jet calculation. The auxiliary proxy
$\chi^2/(3r^4)$ is retained as a comparison. Agreement of the two is an
observed property of these frames and rates; it is not an assumption in the
calculation.

At $T=0.3$ on the fine $dt=0.0005$ trajectory,

$$
R_h\in[11.45342664,15.71659763],\qquad
\max C^2\simeq0.11331608321.
$$

The maximum difference $|R_h-(\chi+2)|$ at that endpoint is about
$7.88\times10^{-8}$. The maximum pointwise difference between the actual
Weyl invariant and the auxiliary proxy is about $1.29\times10^{-9}$.
The record reports these gaps over all 51 frames and on both resolutions
and timesteps.

The endpoint invariant maximum differs by about $5.13\times10^{-11}$ under
timestep refinement and $1.25\times10^{-6}$ between the coarse and fine
spatial resolutions. The latter is about $0.00110\%$ of the measured maximum.
The endpoint dense-versus-FFT metric formula gap is about
$2.64\times10^{-6}$, approximately $2\times10^{-7}$ of the curvature change
from the initial $R_h=2$ value. The five sampled Christoffel
contractions agree with their direct metric-jet formula at about
$2\times10^{-15}$ on the fine endpoint. These remain finite numerical
indicators; the record does not claim a continuum curvature enclosure.

This measures invariant curvature on the realized metric, alongside the
existing regional normal-energy and local-observer measurements. It does
not establish a global horizon, a child region, renewal, or an infinite
extension.

## Reproduction

The producer is
[`derive_nsc_spherical_conformal_curvature.py`](../scripts/derive_nsc_spherical_conformal_curvature.py).
It reads only the sealed v2 trajectory, records every frame's scalar
summaries and stores the initial and final invariant/time-jet profiles for
each case. The postprocessing used about $12.94$ CPU seconds and the payload
is about $669$ kB. Its input and producer hashes are bound.

[`test_nsc_spherical_conformal_curvature.py`](../tests/test_nsc_spherical_conformal_curvature.py)
checks the actual projected directional Jacobian, detects an omitted
projection, verifies an independent Christoffel contraction and binds the
saved invariant profiles to their source trajectory. The calculation does
not change the sealed producers or Cauchy data.
