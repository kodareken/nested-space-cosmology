# Matching static geometry to causal Dirac response

The same smooth Dirac geometry now has a direct Euclidean-frequency versus
Hamiltonian-response comparison. On a finite spatial grid, independent matrix
frequency integration agrees with the vacuum Kubo response. The
periodic-minus-antiperiodic response then supplies a finite relative control
of spatial, angular and regulator limits. No finite gravitational coefficient
is chosen in this calculation.

The [record](../results/nsc-11-response-matching.json) follows
[the covariant source](nsc-covariant-source.md) and
[normalized histories](nsc-influence.md). Its owners are
`src/recursive_horizons/nsc_response_matching.py` and
`scripts/check_nsc_response_matching.py --check`.

## One radius perturbation in the same operator

Use the existing smooth axial circle, R=2, a=1, with N=q=1 and the same
localized periodic profile s(x) used in the vacuum-work experiment. Set

\[
r(t,x)=r_0(x)/(1+J(t)s(x)),\qquad H[J]=H_0+J V,
\quad D_E[J]=D_E+J V_D,\quad V_D=i\beta V.
\]

The warped two-sphere contributes c(d log r) to the base Dirac operator.
Conjugating with r cancels this term, since
r(partial_a+partial_a log r)r^(-1)=partial_a for both time and radial
directions. Thus the time-dependent spinor rescaling cancels the radius
derivative in this fixed-lapse/radial-metric channel. The full
independent angular multiplicity of a representative radial block is 4kappa.
The calculation uses the negative-energy vacuum in each global spin
structure. It does not select one of them dynamically.

For occupied i and empty j, let Delta=E_j-E_i>0 and w=|V_ij|². With
H_int=+J V and Fourier convention exp(+izt),

\[
\chi_R(z)=\sum_{i-,j+}w_{ij}
\left[(z-\Delta_{ij})^{-1}-(z+\Delta_{ij})^{-1}\right],
\qquad \operatorname{Im}z>0.
\]

In particular chi_R(i nu)=-2 sum Delta w/(Delta²+nu²). The independent
Euclidean calculation uses matrix solves for

\[
\int_{-\infty}^{\infty}\frac{d\omega}{2\pi}
\operatorname{tr}\left[D_{\omega-\nu/2}^{-1}V_D
D_{\omega+\nu/2}^{-1}V_D\right].
\]

Its internal-frequency integral runs to infinity with an independent
quadrature transformation. It agrees with the spectral formula for both
P and AP finite grids. This is the ordinary finite Dirac determinant
identity; it does not equate an unregulated loop to a finite heat cutoff.

## The proper-time Hessian and its cutoff

The prescribed action profile is g(d)=E1(d²/Lambda²)/2. Its quadratic
Euclidean radius kernel is evaluated with divided differences

\[
K_\Lambda(\nu)=\int\frac{d\omega}{2\pi}
\sum_{ij}\frac{g'(d_i^-)-g'(d_j^+)}{d_i^--d_j^+}
|\langle i,-|V_D|j,+\rangle|^2,
\qquad g'(d)=-e^{-d^2/\Lambda^2}/d.
\]

The coincident-eigenvalue limit is g''; d^- and d^+ belong to frequencies
omega-nu/2 and omega+nu/2. At zero external frequency, independent analytic
integration yields the finite spatial functional

\[
E_\Lambda=\operatorname{Tr}F_\Lambda(H),\qquad
F_\Lambda(E)=\tfrac12\left[\frac{\Lambda}{\sqrt\pi}e^{-E^2/\Lambda^2}
-|E|\operatorname{erfc}(|E|/\Lambda)\right].
\]

Its spectral Hessian agrees with the frequency calculation. The log-radius
version also agrees with the previous covariant metric-variation owner,
including its full 4kappa weight. Independent second differences of the
energy check both coordinate choices.

## Which local terms cancel

The P-minus-AP difference uses identical local geometry and radius probes.
State-independent local metric counterterms therefore cancel in the
continuum relative response. This supplies a matching test before their
absolute coefficients are fixed. The spin structures and spectra remain
different global data.

An explicit operator contact term is different. For r=r0 exp(Js),
H'=-V and H''=rho1 kappa s²/r0. Then

\[
\frac{d^2 E}{dJ^2}=\chi_R(0)+\langle H''\rangle.
\]

That contact expectation is state-dependent and its P-minus-AP difference
need not vanish. It is retained in the record. Omitting it would produce a
false mismatch between radius coordinates. At finite proper-time cutoff,
the corresponding contact is Tr[F'_Lambda(H)H''].

## Numerical continuation control

The retarded spectral sum is independently Abel-regularized by the positive
weight exp[-(E_i²+E_j²)/K²]. This is a summation device, distinct from the
finite proper-time action. Quadratic Richardson extrapolation in 1/K²
removes the first two analytic cutoff corrections. Spatial resolution is
also varied at fixed K; angular sectors and frequency quadrature are varied
separately. The recorded remainder estimates are observed convergence,
not a rigorous infinite-dimensional bound.

The full angular sum through kappa=8 gives the following relative responses
at proper-time cutoff4, with96 internal-frequency nodes and33 P/32 AP
spatial points. The retarded comparison uses257 P/256 AP points and Abel
cutoff48 with its two smaller extrapolation cutoffs.

| External imaginary frequency | Euclidean relative kernel | Retarded extrapolation minus Euclidean |
|---:|---:|---:|
| 0 | 0.18577154999 | -1.49e-9 |
| 0.6 | 0.17252544942 | -1.41e-9 |
| 1.4 | 0.13062541226 | -1.14e-9 |

The last included angular increment is below7.4e-8; this is an observed
increment, not a bound on the omitted tower. In the dominant kappa=1
channel, increasing33/32 to49/48 spatial points changes the nu=1.4 result
by less than2e-15. The64-to96 frequency-node change is about5e-9; the coarse
32-node value differs by about6.3e-5 and is retained as a resolution control.
The relative log-radius contact for kappa=1 is -0.02201411785 before its
4kappa multiplicity, explicitly demonstrating that it cannot be dropped.

The relative proper-time kernels approach the corresponding extrapolated
vacuum response. Each finite
proper-time value, each Abel value and their differences remain visible.
Equality is a cutoff-removal control; no identification of the two finite
regulators is imposed. In particular, this does not establish a causal
Lorentzian continuation of every finite-cutoff heat form factor.

These are Dirac matter response kernels for a specified radius channel.
They are not a constrained graviton Hessian or a backreacted solution.
The full functional still needs its finite completion, compensator,
physical state/domain selection and joint field/metric equations. The
calculation strengthens that next assembly by making the static-to-causal
matching and its required contact term explicit.

The influence-functional use of quantum stress response has established
prior work, including [Martín and Verdaguer](https://arxiv.org/abs/gr-qc/9904021).
The present project contribution is this executable match on its smooth
Dirac geometry; the Kubo and analytic-continuation identities are prior art.
