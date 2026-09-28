# Which finite terms affect the absolute metric equation?

The [calculation](../scripts/check_nsc_finite_terms.py) and
[record](../results/nsc-8-finite-terms.json) distinguish **four independent local
bulk coefficients** on each actual smooth R=2 and R=4 profile. The Euler term
and the Laplacian of curvature have zero bulk variation on the declared closed
cell. Radial coordinate freedom produces a relation between metric equations;
it does not make one of the four coefficients disappear.

This moves beyond the earlier constant-cylinder Weyl-squared example. The
constant cylinder has only rank two in the probes used here. On the varying
neck, the Weyl-squared term contributes to the radial-null stress and remains
independent of the scalar-curvature-squared term. Their coefficients have not
been determined by the currently specified measure, so the absolute-stress
stationarity problem remains underdetermined.

These are sensitivities of the same action's unresolved finite local part.
No independent gravitational action, scalar field, fitted counterterm, or
stationary root is introduced. The positive result is an explicit map from
missing measure data to the remaining physical equations.

## Functional and curvature convention

Retain every metric function before variation:

\[
ds^2=N(x)^2dt^2-q(x)^2dx^2-r(x)^2d\Omega_2^2,\qquad x\sim x+L.
\]

Here x and t have length units, N and q are dimensionless, r has length,
and M has inverse length. M is the normalization mass, not a fermion mass.
Use the same Riemann convention as the [smooth-geometry
calculation](nsc-smooth-geometry.md). In particular, the Lorentzian scalar
curvature of the constant cylinder is -2/a². Its positive Wick-rotated metric
has the opposite scalar curvature and identical squared curvature invariants.

Define the coefficient convention through the static energy functional

\[
F_{\rm fin}=\sum_i c_i F_i,\quad
F_i=4\pi\int_0^L Nqr^2 I_i\,dx,\quad
I_i=(M^4,M^2\mathcal R,C^2,\mathcal R^2,E_4,\Box\mathcal R),
\qquad S_{\rm fin}=-\int F_{\rm fin}\,dt.
\]

The M²R coefficient therefore includes this static-action and curvature sign
convention. If an Einstein term is isolated from the completed common action,
its coefficient is counted once. The six symbols describe unresolved finite
data; the calculation never assigns them values. Fixed-M variations differ
from changes of units, which also scale M inversely.

Set D=q⁻¹∂x and introduce four curvature quantities

\[
A=\frac{D^2N}{N},\qquad B=\frac{DN\,Dr}{Nr},\qquad
C=\frac{D^2r}{r},\qquad D_s=\frac{1-(Dr)^2}{r^2}.
\]

Direct Christoffel and Riemann contraction of the general four-dimensional
metric, independently of these reduced expressions, verifies

\[
\begin{aligned}
\mathcal R&=2(A+2B+2C-D_s),\\
R_{\mu\nu}^2&=(A+2B)^2+(A+2C)^2+2(B+C-D_s)^2,\\
R_{\mu\nu\rho\sigma}^2&=4(A^2+2B^2+2C^2+D_s^2),\\
C_{\mu\nu\rho\sigma}^2&=\frac43(A-B-C-D_s)^2,\\
E_4&=8(-AD_s+2BC).
\end{aligned}
\]

All three direct-contraction residuals and both contracted curvature identities
reduce symbolically to zero for arbitrary N, q and r. The symbol C without
indices in this display is the radial curvature quantity; C² with tensor
indices is the Weyl invariant.

## Independent variations and exact null directions

For the reduced density L_i=Nqr²I_i, the implementation differentiates its
metric jets and forms the Euler–Lagrange densities

\[
e_{if}=4\pi\left[
\frac{\partial L_i}{\partial f}
-\partial_x\frac{\partial L_i}{\partial f'}
+\partial_x^2\frac{\partial L_i}{\partial f''}\right],
\qquad f=N,q,r.
\]

The local stress sensitivities, per unit c_i, are

\[
\rho_i=\frac{e_{iN}}{4\pi qr^2},\qquad
p_{xi}=-\frac{e_{iq}}{4\pi Nr^2},\qquad
p_{\perp i}=-\frac{e_{ir}}{8\pi Nqr}.
\]

The exact reduced boundary identities are

\[
Nqr^2E_4=\partial_x\left[-8\frac{N'}q
\left(1-\frac{r'^2}{q^2}\right)\right],\qquad
Nqr^2\Box\mathcal R=-\partial_x\left[\frac{Nr^2}{q}\mathcal R'\right].
\]

Periodic variations of this fixed-topology cell therefore annihilate their
bulk integrals. The Euler density is nevertheless differentiated before
numerical cancellation as a control. The box term's zero gradient follows from
its exact divergence. These conclusions need not survive physical boundaries
or topology changes without the corresponding boundary terms.

For every coefficient, radial coordinate transformations satisfy

\[
\delta N=\xi N',\quad \delta q=\xi q'+q\xi',\quad \delta r=\xi r',
\qquad N'e_N+r'e_r-q\partial_xe_q=0.
\]

The implementation checks the local identity and integrated gauge probe on
nonconstant N and q. An additional finite coordinate pullback uses the analytic
profile evaluated at x+ε sin(2πx/L), with q multiplied by its Jacobian; all six
integrals remain equal. No division by r' is used at the neck.

The Weyl identities remain different from this coordinate identity.
For δ(N,q,r)=σ(N,q,r), C² has exactly zero pointwise density variation.
The local integrated-variation densities are respectively 4L_M4, 2L_M2R,
zero, -12L_boxR, zero and zero, including the common angular factor. This
does not declare common Weyl rescaling a gauge symmetry of the incomplete
full measure, or remove the physical relative-scale/link equation.

## Response matrix on the actual smooth neck

Use the already declared smooth family

\[
r(x)=\sqrt{a^2+[\sin(kx)/k]^2},\quad
k=\pi/(2R),\quad x\in[-R,R),\quad a=1.
\]

After taking variations, evaluate N=q=1. Nine probes comprise three lapse
perturbations δN=Nf, three radius perturbations δr=rf, and three spatial
volume-preserving shape perturbations δq=qf, δr=-rf/2. In each group,
f=1, cos(2πx/L), cos(4πx/L). The last group preserves qr² pointwise to first
order. The matrix is

\[
J_{Ai}=\int dx\,\sum_{f=N,q,r}\delta_A f\,e_{if}.
\]

Every entry is stored. For rank assessment only, the first four columns are
divided by their Euclidean norms; the boundary columns have unit scale.
The absolute singular-value threshold is 10⁻⁸.

| Background | Nonzero normalized singular values | Rank |
|---|---|---:|
| Constant cylinder a=M=1, L=4 | Two nonzero values; full matrix in record | 2 |
| Smooth R=2 | 1.39721723, 1.01589563, 0.77305076, 0.64663174 | 4 |
| Smooth R=4 | 1.38077045, 0.99766263, 0.89610691, 0.54326295 | 4 |

The smooth matrices have only the two exact right null vectors
(0,0,0,0,1,0) and (0,0,0,0,0,1). These are the Euler and total-derivative
coefficients. Gauge directions act on metric probes on the other side of the
matrix. Generic linear dependencies among nine probes in a four-dimensional
response space are not all claimed to be gauge transformations.

At a=M=1 the cylinder has two additional right null vectors:
(1,1,0,1/4,0,0) and (0,0,1,-1/3,0,0). They arise from the restricted constant
geometry. Varying the radius along the cell resolves both degeneracies.
The cylinder control recovers 16παL/(3a²) for the C² energy and zero radial-null
response for every local term on that product.

## Exact throat equations that were previously missing

At the smooth neck, r'=r'''=0, r''=1/a, and
r''''=-3/a³-4k²/a. Substituting these jets **after** general variation gives

\[
\begin{aligned}
(\rho+p_x)_{M^4}&=0,\\
(\rho+p_x)_{M^2\mathcal R}&=4M^2/a^2,\\
(\rho+p_x)_{C^2}&=\frac{32(1+a^2k^2)}{3a^4},\\
(\rho+p_x)_{\mathcal R^2}&=-\frac{16(1+4a^2k^2)}{a^4}.
\end{aligned}
\]

The record includes exact expressions for density and both pressures, together
with all 128 spatial samples and all six coefficient sensitivities. In the
declared units, the neck null sensitivities are

| R | M²R | C² | R² |
|---:|---:|---:|---:|
| 2 | 4 | 17.2464029 | -55.4784176 |
| 4 | 4 | 12.3116007 | -25.8696044 |

Thus an unspecified finite C² part changes the very null residual needed for
the smooth self-sourcing equation, while preserving its Weyl identity. The
opposite signs of C² and R² sensitivities are not a prescription to tune their
coefficients against each other. The four coefficients must come from the
specified common measure or independent matching data.

## Reproduction and matching conditions

The grid calculation uses 32, 64 and 128 points. Real Fourier first derivatives
are antisymmetric, and the second derivative is their square; both Nyquist
symbols vanish. The symbolic-jet variation uses these exact discrete adjoints.
Pointwise continuum product rules are checked under refinement, not silently
assumed at finite resolution. At R=4 on 32 points, the raw Euler column has a
spurious singular value of 3.41×10⁻⁸, producing numerical rank five at the
declared threshold. It falls below 10⁻¹³ on 64 and 128 points. The exact
continuum divergence establishes why this coarse-grid value is an error,
rather than another physical coefficient. The integrated matrices stabilize much earlier
than fourth-derivative local residuals, whose roundoff is explicitly retained.

Controls include direct full-metric curvature contraction, exact Weyl and
Euler identities, fixed-M spatial scaling and joint unit scaling, a finite
local conformal transformation, finite coordinate pullback, and centered
energy changes of every N, q and r component on a nonconstant N,q profile.
The independent scalar-energy route contracts curvature invariants and never
uses the analytic response gradient. The perturbation steps are halved
separately from the spatial grid. No spectral cutoff of a spatial Hamiltonian
is called a physical covariant regulator in this calculation.

A separate 64-point check varies the individual throat node of N, q and r.
Its finite-energy derivatives agree with the local analytic lattice gradient
to 5×10⁻⁶ absolute error at step 10⁻⁶, with larger steps retained to show
convergence. This checks local derivatives of the discrete functional; the
single-node perturbation is not asserted to be a smooth continuum metric.

Run:

```sh
python3 scripts/check_nsc_finite_terms.py --check
PYTHONPATH=src:. python3 -m unittest discover -s tests -p test_nsc_finite_terms.py -v
```

The all-field comparison includes every sample, matrix, exact expression,
scope flag and source hash. Float tolerances are atol=2×10⁻⁶ and rtol=2×10⁻⁸
to accommodate local fourth-derivative roundoff. The calculation separately
checks integrated matrix refinement at 2×10⁻⁷ and exact throat matching at
10⁻⁷ absolute tolerance. Existing outputs cannot be overwritten by `--output`,
and `--output` and `--check` are mutually exclusive.

The next equations now have explicit coefficient dependence:

\[
\mathcal R_f(x)=\mathcal R_f^{\rm nonlocal,state,link}(x)
+\sum_{i=0,R,C,R^2}c_i e_{if}(x)=0,\qquad f=N,q,r,
\]

or b_A+J_Ai c_i=0 for the integrated probes. The missing b_A values must be
calculated independently from the same covariant functional, state and link;
setting b_A to obtain the imposed neck would merely fit a target. Four
independent finite matching conditions are needed to fix the four bulk
coefficients unless the completed measure supplies them directly. Euler and
box coefficients require additional boundary/topology data if those sectors
are used. Link, compensator and state equations remain necessary.

If a completed action has an invertible Hessian H in physical metric
coordinates at a stationary solution, its local sensitivity is
du_star^A/dc_i=-(H⁻¹)^AB J_Bi. This calculation supplies J, not that Hessian
or a stationary solution. It does not construct a complete ultraviolet
functional or an in-in measure, determine absolute AP stress, or close the
recursive tower.

The finite/anomalous factor distinction comes from
[Andrianov, Kurkov and Lizzi](https://arxiv.org/html/1106.3263v1), while
[Dappiaggi, Hack and Pinamonti](https://arxiv.org/abs/0904.0612) develop
covariantly conserved Dirac stress for backreaction. Those established results
motivate the missing data; the contribution here is the reproducible response
rank, nullspace and throat sensitivities for the project's smooth geometry.
