# Five-dimensional ultraviolet matching for the compact interaction

The leading proper-time bulk functional fixes the same factor-nine torsion
stiffness in five dimensions. Its quadratic derivative correction has now
been matched in that dimension. These results supply the cutoff-dependent
part of the [compact interaction](nsc-compact-interaction.md); a finite
relative coefficient, the boundary completion and the quantum state remain
to be specified by the full common functional.

## Reuse the dimension-independent result

[Pfäffle–Stephan, Proposition5.4 and equations41–43](https://arxiv.org/html/1101.1424v3)
give the Bochner formula and a0,a2 in arbitrary dimension. Their simplified
curved a4 in Proposition5.5 is explicitly four-dimensional and is not reused
as a five-dimensional formula. For the needed flat quadratic derivative
piece, use [Vassilevich, equation4.28](https://arxiv.org/html/hep-th/0306138v3).

Let K be the skew connection difference, so physical torsion is 2K and
the component norm is ||K||²_comp=6|K|²_form. For one complex
four-component spinor in five Euclidean dimensions, the bulk heat trace is

\[
\operatorname{Tr}e^{-tD_K^2}\sim(4\pi t)^{-5/2}\int\sqrt g\,
\left[4+t\left(-\frac R3+3\|K\|_{\rm comp}^2\right)
+t^2 a_4(x)+\cdots\right].
\]

The quoted a2 formula already establishes the ratio 9 at equal Einstein
coefficient. This is a published result applicable to our dimension, not a
new torsion theorem. These are local interior coefficients. The actual
compact endpoint conditions can add independent boundary contributions.

## Separate the matched determinant part from the remaining fermions

Introduce a proper-time matching scale 0<nu_match<Lambda, distinct from
both the physical cutoff Lambda and determinant normalization M:

\[
\Gamma_{\Lambda,\nu}=
\frac12\int_{\Lambda^{-2}}^{\nu^{-2}}\frac{dt}{t}
\operatorname{Tr}e^{-tD_K^2}.
\]

The complementary one-loop modulus uses proper times above nu^-2.
This split prevents adding the full determinant again after taking its
ultraviolet coefficients as an induced action. It is a one-loop bookkeeping
identity; a complete causal coarse-graining and measure are still required
when retaining interacting low-energy fermions. Proper time is not a sharp
partition of KK modes.

The corresponding weights of a0,a2,a4 are
(Lambda^5-nu^5)/5, (Lambda^3-nu^3)/3 and Lambda-nu. In particular,

\[
\Gamma_{\Lambda,\nu}\supset\int\sqrt g\left[
\frac{4(\Lambda^5-\nu^5)}{5(4\pi)^{5/2}}
-A_5 R+9A_5\|K\|_{\rm comp}^2\right],\qquad
A_5=\frac{\Lambda^3-\nu^3}{9(4\pi)^{5/2}}.
\]

The volume term is retained. Omitting it in an unconstrained metric variation
would change the source problem. This calculation adds no independently
weighted Einstein action.

If only this leading matched part is kept, kappa5²=1/(2A5) and s_T=9.
The common prefactor of the previously evaluated compact vertex is then

\[
\left.\frac{\kappa_{4,\rm bulk}^2}{32s_T}\right|_{\rm leading}
=\frac{(4\pi)^{5/2}}{64(\Lambda^3-\nu^3)I_3},
\qquad I_3=\int e^{3\sigma}dY.
\]

This identifies the leading contribution. It does not assign that value to
the full renormalized interaction in the previous record, where s_T remains
unresolved.

## The derivative correction that controls the local contact approximation

Around flat space and K=0, the integrated quadratic derivative part of a4
comes from tr(E1²)/2+tr(Omega1,ij Omega1,ij)/12. E1 and Omega1 follow
directly from the published Bochner connection. Two representative
polarizations with derivative along x0 fix the two possible parity-even
quadratic derivative invariants:

| 3-form profile | E1 contribution | Connection-curvature contribution | Sum |
|---|---:|---:|---:|
| tau(x0) dx0 wedge dx1 wedge dx2 | 0 | -3(tau')² | -3(tau')² |
| tau(x0) dx1 wedge dx2 wedge dx3 | 9(tau')²/2 | -9(tau')²/2 | 0 |

The Euclidean spinor calculation therefore gives the integrated bulk term

\[
a_4^{(K^2,\,\partial^2)}=-3|\delta K|_{\rm form}^2,
\qquad
\Gamma_{\Lambda,\nu}\supset
-\frac{3(\Lambda-\nu)}{(4\pi)^{5/2}}\int|\delta K|_{\rm form}^2.
\]

Total divergences are excluded only for this integrated interior test.
Curvature–torsion couplings, higher powers of K, higher heat terms and
boundary terms are not evaluated by this two-polarization calculation.

For the first polarization, the displayed quadratic coefficients have ratio

\[
\mathcal M_{\rm coeff}^2
=2(\Lambda^2+\Lambda\nu+\nu^2),\qquad
\frac{\text{displayed derivative contribution}}{\text{algebraic contribution}}
=-\frac{q_E^2}{\mathcal M_{\rm coeff}^2}.
\]

This is a coefficient scale, **not a physical torsion mass** or a pole of
the full kernel. The formal ultraviolet limit nu/Lambda→0 gives
2 Lambda². A zero extrapolated from this truncated Euclidean polynomial
does not establish Lorentzian propagation or an instability.

For a controlled local expansion of the matched window, the energy scales
q_E, |K|, sqrt(||nabla K||) and sqrt(||Riemann||) must be small compared
with nu_match; the small displayed ratio alone is insufficient. In the unwarped compact
reference, the first massive mode's current contains a harmonic q=pi/L_star.
Its displayed relative correction is

\[
\frac{\pi^2}{2\zeta(1+u+u^2)},\qquad
\zeta=(\Lambda L_\star)^2,\quad u=\nu/\Lambda.
\]

This is a scale-separation check, not an equation selecting zeta. Proper
momenta and derivative terms on the warped carrier require the actual
metric and boundary conditions. No old stationary-scale candidate is used.

## The finite datum that the leading calculation cannot fix

Write any remaining permitted local finite contributions to these channels
as -c_R R+c_T||K||²_comp. The complete functional must determine or exclude
these contributions; they are not set to zero here. The resulting local
ratio would be

\[
s_T^{\rm full}=\frac{9A_5+c_T}{A_5+c_R},\qquad
s_T^{\rm full}-9=\frac{c_T-9c_R}{A_5+c_R}.
\]

Thus matching the Einstein coefficient alone does not fix the torsion
coupling. Leading cutoff compensation changes c_T and c_R in the ratio9,
leaving c_T-9c_R unchanged. This is not a solved renormalization flow or a
recursive fixed point. Repeating the leading heat calculation cannot fix
this finite relative datum.

The next source calculation must therefore supply the finite matching and
the compatible low-energy state, using the existing coupled-mode vertex.
It must retain the complementary determinant contribution, compact boundary
action and regulator/normalization matching without double-counting fermions.
The current bulk calculation does not determine an absolute vacuum stress,
the full quantum measure, or a self-sourced transition.

The [exact development record](../results/development/torsion-uv-map.json)
is generated and authenticated by
[the focused runner](../scripts/check_nsc_torsion_uv_map.py). It uses the
existing Euclidean spinor matrices and checks both polarizations, proper-time
coefficient ratios and the finite-matching identity. All recorded scientific
values are exact strings or integers. Run
`python3 scripts/check_nsc_torsion_uv_map.py --check`.
This subsequent development evidence does not modify the v0.3.0 release.
