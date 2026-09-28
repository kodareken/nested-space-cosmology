# Carry the curvature correction into the same source action

The known local curvature-squared contribution can be retained in the
spherical action as a perturbative source interaction. A checked metric
field redefinition gives its coefficient and keeps the original physical
metric explicit. The resulting canonical constraints close through the
same first order. No independently weighted interaction or new fundamental
scalar is introduced.

This is a four-dimensional effective-field-theory (EFT) construction. It is
not a replacement of the full nonlocal spectral action, a nonperturbative
fourth-order gravity theory, or a proof that the expansion is valid at the
imposed throat.

## Reused methods and the new interface

The [compact/light map](nsc-compact-matching.md) already gives the Dirac
contributions to the local a4 coefficients. The [canonical spherical
action](nsc-spherical-action.md) supplies the leading metric/source system.
The missing local step is to retain the curvature correction in that same
system without assigning additional initial data to an EFT truncation.

Reuse [Parker–Simon](https://arxiv.org/html/gr-qc/9211002v1), especially
section5, for perturbative order reduction. Reuse
[Ruhdorfer–Serra–Weiler, sectionIII.1](https://arxiv.org/html/1908.08050v2)
for removing redundant curvature operators by metric redefinitions.
The calculation here checks the coefficients, matter terms, metric
interpretation and constraint map for this source. It does not claim a
new general order-reduction theorem.

## An off-shell bulk identity fixes the interaction

Use the laboratory conventions

\[
S_0=\int\sqrt{|g|}[-A R_L+\mathcal L_m],\qquad
T_{\mu\nu}=\frac{2}{\sqrt{|g|}}\frac{\delta S_m}{\delta g^{\mu\nu}},
\qquad E_{\mu\nu}=G_{\mu\nu}-\frac{T_{\mu\nu}}{2A}.
\]

T includes the vacuum term as well as the electromagnetic and collective
source. The full coefficients remain symbolic. In particular the zero
Dirac contribution to R-squared is not an assignment of the complete R-squared
coefficient to zero.

Write the original local correction as -cW C^2-cR R^2-cEuler E4-cBox boxR.
In four dimensions,

\[
C^2=E_4+2(R_{\mu\nu}R^{\mu\nu}-R^2/3).
\]

To first order in cW,cR, the inverse-metric map is

\[
g_{\rm old}^{\mu\nu}=g_{\rm new}^{\mu\nu}+\Delta g^{\mu\nu},
\]
\[
\Delta g^{\mu\nu}=-\frac{2c_W}{A}
 \left[R^{\mu\nu}-\frac R6g^{\mu\nu}
       +\frac{T^{\mu\nu}-Tg^{\mu\nu}/3}{2A}\right]
 +\frac{c_R}{A}\left(R-\frac{T}{2A}\right)g^{\mu\nu}.
\]

The variation of S0 is -A integral sqrt|g| E_mu_nu Delta g^mu_nu, up to its
boundary variation. Direct contraction with arbitrary symmetric E and T
tensors verifies that the transformed bulk correction is

\[
\boxed{\Delta\mathcal L_m=
 -\frac{c_W}{2A^2}\left(T_{\mu\nu}T^{\mu\nu}-\frac{T^2}{3}\right)
 -\frac{c_R}{4A^2}T^2.}
\]

The Euler coefficient becomes cEuler+cW in the inherited static-energy
convention. Euler, box-R and the induced boundary variations are retained
separately. The Euler identity is not topological in five dimensions; this
step applies after the four-dimensional reduction.

This is an off-shell action identity through the declared order, modulo
boundary terms and the stated change of variables. Simply substituting the
Einstein equation into an action without carrying the metric map would not
establish the same result. Keeping both the original curvature term and its
replacement contact would double-count it.

## The charge source supplies a specific spherical contact

Let H=4 pi r^2 be the sphere area, K=(grad chi)^2 and W=mu^2 chi^2/2.
The [previous charge-source calculation](nsc-gauge-source.md) identifies
the radial electromagnetic energy

\[
\rho_{\rm EM}=\rho_B+W/H,\qquad
\rho_B=\frac{C_Fq^2}{2r^4},\qquad p_{\rm sphere}=\rho_{\rm EM}.
\]

The classical charge/Maxwell source is trace-free in four dimensions when
that angular pressure is included. With vacuum energy V the total trace
is 4V. Direct tensor contraction gives

\[
T_{\mu\nu}T^{\mu\nu}-T^2/3
 =\frac{K^2}{2H^2}+4\rho_{\rm EM}^2-\frac43V^2.
\]

Consequently the induced two-dimensional contact density is

\[
\Delta\mathcal L_2=
 b\left[\frac{K^2}{2H}+4H\rho_{\rm EM}^2-\frac43HV^2\right]
 +16b_RHV^2,
\qquad b=-\frac{c_W}{2A^2},\quad b_R=-\frac{c_R}{4A^2}.
\]

Its derivative interaction, charge self-interaction, magnetic cross term
and vacuum term are fixed together. Eliminating the electric field at its
leading stationary value is valid for this first-order insertion; its
first-order shift multiplies the already vanishing leading electric equation.

The displayed specialization contains one charged collective mode. If
additional matter or neutral CFT sectors are retained, they belong in the
total T before it is squared. Summing separate stress squares would omit
the cross interactions. Quantum use also needs the renormalized composite
products and connected correlations, not a substitution of products of
mean stresses or an unmodified classical trace-free formula.

## The correction fits the canonical constraints at EFT order

Use the previous first-order variables and a bookkeeping parameter eps
which multiplies the curvature correction. In the leading canonical
variables,

\[
K_0=\frac{(\chi')^2-\pi_\chi^2}{2e^-_1e^+_1}.
\]

The perturbative Legendre transform gives H_m,1=-L_m,1 evaluated at the
leading velocity. Calling the contact F(X,chi,K0), the constraint correction
is

\[
\Delta G_1=0,\qquad
\Delta G_2=e^+_1F,\qquad
\Delta G_3=-e^-_1F.
\]

The implementation checks the smeared brackets for the class
F=d(X)K0^2+F0(X,chi), then binds the actual contact to it. The actual
d(X)=-cW/(2 A X). All derivative-delta residuals and closure residuals
vanish at orders eps^0 and eps^1. The new contribution to the 23 structure
coefficient multiplying G1 is partial_X F; the other first-order structure
corrections vanish in this convention.

This does not assert closure of a resummed G0+eps G1 at arbitrary eps.
Terms of second and higher order, and the corresponding Legendre transform,
are outside the record. The perturbative branch uses the leading theory's
initial data rather than adding the extra modes of a fourth-order truncation.

The new derivative contact also changes the matter dependence on the
coframe. The old exact linear geometry integration cannot simply be reused
unchanged. It must be handled as an EFT insertion with its corresponding
source, state and measure transformations.

## A change of variables does not remove physical vacuum curvature

A useful control has T_mu_nu=V g_mu_nu and lambda=V/(2 A). The contact
appears to shift the vacuum term in the new metric variables. The metric
map simultaneously rescales those variables. Pulling the curvature back
to the original physical metric cancels that apparent shift exactly through
first order. Thus no cosmological constant cancellation or new vacuum
prediction follows from the field redefinition.

The same requirement applies to the areal radius, lapse, observer clock,
boundary data and matter embedding. X in the corrected canonical action is
the area of the redefined metric. It is not automatically the original
operator's measured area.

## Quantum and validity requirements remain explicit

At quantum level the variable change carries

\[
\mathcal J_F=\det\!\left[\frac{\delta g_{\rm old}}
                                  {\delta g_{\rm new}}\right],
\]

together with any accompanying field, ghost, cutoff and state changes.
This finite-regulator Jacobian is not set to one. The induced-gravity
power counting must determine its order; an argument using dimensional
regularization alone cannot identify this project's finite prescription.
Boundary and zero-mode conditions must be transformed with the same map.

The construction requires |c_i R|/A and |c_i T|/A^2 to be small, and external
scales below the matched local expansion scales. That regime is not yet
established at the imposed throat. The full nonlocal h(D4^2) remainder still
belongs to the theory. Order reduction does not prove its physical-mode
health or justify hiding an instability within the regime being claimed.

No physical coupling, quantum state or self-sourced geometry is selected
here. The result is the explicit local source correction and its consistent
canonical/observable map at first EFT order.

The record `results/development/curvature-eft.json` is reproduced by
`scripts/check_nsc_curvature_eft.py --check`, using every formula, residual,
scope field and authenticated dependency. Existing generators are not
rerun. A bounded Grok review attempt returned no usable report and is not
counted as evidence; the checks use the primary equations and the direct
symbolic action, stress, Legendre and constraint calculations.
