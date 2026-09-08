# A covariant Dirac source on the smooth geometry

The static four-dimensional Dirac determinant and its metric variations are
now evaluated on the smooth varying throat profile. Time frequency, both
angular signs and the correct spacetime measure are retained together. The
calculation supplies an explicit determinant contribution to the source
equation that previously contained an unevaluated term.

Three independent connections close in this realization: the frequency
integral reproduces the known cylinder vacuum response; the varying-metric
stress obeys the conservation and scaling identities; and its short-distance
trace reproduces the covariant Dirac heat coefficients. A prescribed
ultraviolet subtraction yields a finite determinant remainder with negative
neck null stress. Its finite completion, common compensator and physical
scale still have to be fixed by the same action. The source is not fitted to
hold the neck open.

The [record](../results/nsc-9-covariant-source.json) is generated and checked by
`scripts/check_nsc_covariant_source.py`. The numerical owner is
`src/recursive_horizons/nsc_covariant_operator.py`; an independent symbolic
owner, `src/recursive_horizons/nsc_covariant_identities.py`, checks the full
spin connection, Lichnerowicz identity, metric variation and heat coefficients.

## One spacetime operator before variation

Use the static positive Euclidean metric

\[
ds_E^2=N(x)^2d\tau^2+q(x)^2dx^2+r(x)^2d\Omega_2^2,
\qquad x\sim x+L.
\]

N, q and r are smooth and positive. The spatial circle has periodic (P) or
antiperiodic (AP) spin structure. This is the compact static development
geometry, not a continuation through the shifted trapped PG region.

The torsionless Dirac connection contains the radial term
q^-1(partial_x+N'/(2N)+r'/r). With the four-dimensional unitary rescaling

\[
\chi=r\sqrt{Nq}\,\psi,\qquad
p_q=-i\left(q^{-1}\partial_x+\tfrac12(q^{-1})'\right),
\]

the measure becomes flat in x and tau after angular separation. The lapse
connection must be included before rescaling; omitting it produces a
nonzero defect in the independent symbolic check.

In the paired angular representation used by the preceding boundary result,

\[
\alpha=I\otimes\rho_2,\quad \beta=I\otimes\rho_3,\quad
H_0=\alpha p_q+(\eta_3\otimes\rho_1)\kappa/r,
\]

and the self-adjoint Fourier-frequency operator is

\[
\boxed{D_\omega=\beta\omega/N+i\beta H_0.}
\]

The physical four-dimensional chirality anticommutes with D_omega. Its square
retains an explicit lapse-gradient contribution:

\[
\boxed{D_\omega^2=H_0^2+\omega^2/N^2+
                     \alpha\,\omega N'/(qN^2).}
\]

Only at constant unit lapse does this reduce to omega²+H0². The independent
unseparated calculation verifies every matrix entry of

\[
D_E^2=\nabla^*\nabla+R_E/4.
\]

The curvature term is already present after squaring; it is not added again.

## Cutoff functional, normalization and physical counting

For one complex four-component massless Dirac field, define the modulus per
coordinate Euclidean time by

\[
\boxed{E_\Lambda=
\frac12\sum_{\kappa\ge1}2\kappa
\int_{\mathbb R}\frac{d\omega}{2\pi}
\operatorname{Tr}_{x,4}E_1(D_\omega^2/\Lambda^2).}
\]

The factor 2 kappa multiplies the full paired four-component block. Static
real metric functions and symmetric signed Fourier momenta permit a reduction
to one angular-sign two-component block and positive frequency with prefactor
4 kappa/pi multiplying the half trace. Both symmetries are tested directly.

Periodic grids have odd size and include symmetric integer momenta; AP grids
have even size and symmetric half-integer momenta. The derivative matrix is
Hermitian. This new Fourier representation is independently compared with
the preceding staggered regulator. Finite matrices do not acquire a continuum
Leibniz rule by declaration: geometric Weyl variation and a direct matrix
sandwich can differ as operators, while their cutoff-weighted agreement is
checked separately.

The integral uses a physical proper-time cutoff Lambda and independent
spatial, angular and frequency refinements. A finite-operator lower estimate

\[
D_\omega^2\ge\omega^2/N_{\max}^2
 -|\omega|\,\|[p_q,N^{-1}]\|
\]

sets the frequency endpoint. E1(u)<=exp(-u)/u then gives the reported analytic
upper bound on omitted energy for the retained finite angular/spatial
operator. Derivative tails and other truncations require their own checks;
binary64 evaluation is not an interval certificate.

The old finite-rank normalization proportional to
rank(D)[log(M/Lambda)+gamma_E/2] is not integrated over unrestricted frequency
and angular labels: that would introduce an unregulated trace of the identity.
M enters the explicit subtraction below, not the Dirac spectrum as a fermion
mass. A global determinant phase, zero-mode law and complete invariant
functional remain separate parts of the common-action definition.

## Metric source and the clock test

Each frequency fiber is varied before integration. Spectral functional
calculus gives

\[
\delta E_\Lambda=-\sum_\kappa2\kappa\int\frac{d\omega}{2\pi}
\operatorname{Tr}\left[e^{-D_\omega^2/\Lambda^2}
D_\omega^{-1}\delta D_\omega\right].
\]

The implementation evaluates the inverse through its eigenvalues and retains
the noncommuting second variation. Local energy-density and pressure
conventions follow the earlier static source calculation:

\[
\rho=\frac{E_N}{4\pi qr^2},\qquad
p_x=-\frac{E_q}{4\pi Nr^2},\qquad
p_\perp=-\frac{E_r}{8\pi Nqr}.
\]

Here E_f denotes the functional derivative density per dx; the finite nodal
gradient is divided by the grid spacing. The conservation identity is

\[
p_x'+\frac{N'}N(\rho+p_x)+2\frac{r'}r(p_x-p_\perp)=0.
\]

For the nonconstant N,q control on the smooth R=2 profile, the maximum
conservation residual decreases from approximately 6.31e-6 at 16 points to
4.58e-8 at 24, 2.37e-10 at 32, and below 5e-14 at 48. Independent finite
variations of all three metric fields agree with the spectral gradients;
the noncommuting Hessian agrees with separate finite energy differences.

Uniform lapse scaling is an especially direct physical check:

\[
D_\omega[cN]=D_{\omega/c}[N],\qquad E_\Lambda[cN]=cE_\Lambda[N].
\]

At Lambda=2, the unit-lapse smooth profile has energy about 8.6421970.
Multiplying the lapse by 1.2 gives about 10.3706364. Applying the same fixed
cutoff to the spatial canonical Hamiltonian instead gives about 4.8197489.
The latter is a different regulator. This explicitly identifies why a spatial
cutoff cannot supply the missing spacetime lapse variation.

## Independent vacuum and ultraviolet connections

At unit lapse, integration over time frequency produces, for each spatial
eigenvalue h,

\[
\frac{\Lambda}{2\sqrt\pi}e^{-h^2/\Lambda^2}
-\frac{|h|}{2}\operatorname{erfc}(|h|/\Lambda).
\]

This agrees with the new frequency quadrature. On the constant cylinder,
subtracting the zero-winding contribution reproduces the independent
proper-time winding sum for both P and AP. At L=4,a=1,Lambda=2, the AP null
stress is about -0.00380414334 in both calculations.

On the same varying R=2 geometry, the covariant calculation at Lambda=3 gives
P-minus-AP energy approximately 0.158657420675 and neck null-stress difference
0.0186087268245. The prior staggered 256/512-point values extrapolate to the
same values within 2e-8. The extrapolation is a numerical control, not a
rigorous continuum bound. It connects the new covariant source to the existing
finite stress difference without modifying either operator to force agreement.

For D²=-(nabla²+E), E=-R_E/4, the computed spin-curvature trace gives

\[
a_0=4,\qquad a_2=-R_E/3,\qquad
a_4=(-18C^2+11E_4-12\Box_E R_E)/360.
\]

These use the universal Laplace heat formula; see
[Vassilevich, equations 4.26–4.28](https://arxiv.org/pdf/hep-th/0306138).
They are established coefficients, independently reconstructed here to bind
the actual operator and metric-response conventions. Numerical frequency and
angular integration on nonconstant N,q,r reproduces the expansion through a4
with a remainder decreasing as the next proper-time order. No extra R/4 or
independently weighted gravitational action is inserted.

## Computed subtraction remainder and the next matching equation

Write A_i for the integrated heat coefficients per coordinate time. The
explicitly defined determinant remainder is

\[
E_{\rm sub}(M)=E_\Lambda-
\frac{A_0\Lambda^4/2+A_2\Lambda^2+
A_4\log(\Lambda^2/M^2)}{32\pi^2}.
\]

The subtraction and its metric gradients use the independently calculated
finite-term owner. R_E=-R_L, while the squared invariants and static
Box_E R_E=Box_L R_L agree. M is held fixed during physical metric variation.

For a=1,R=2,M=1, the computed remainder is:

| Lambda | Remainder energy | Neck rho+p_x |
|---|---:|---:|
| 1 | -0.08044646 | -0.00879397 |
| 2 | -0.09076188 | -0.00970810 |
| 3 | -0.09057287 | -0.00928097 |
| 4 | -0.09049882 | -0.00910693 |
| 6 | -0.09044317 | -0.00897192 |

This is a calculated negative null contribution on the varying neck under
the stated subtraction. It is not a scheme-independent complete source or a
choice of finite coefficients that solves the metric equation. Normalization
changes satisfy dE_sub/dlogM=A4/(16pi²); compensating changes of the complete
action's finite coefficients must preserve physical predictions.

In the declared +--- stress convention the limiting trace is
rho-p_x-2p_perp=a4/(16pi²). This has the opposite sign to a trace convention
defined by deltaW=-int sigma T. The numerical trace approaches this local
coefficient as Lambda increases. Finite R² terms can alter the removable
Laplacian contribution; the full invariant functional is still to be matched.

The new record now supplies actual determinant source projections b_A on the
same lapse/radius/shape probes used in the finite-term calculation:

\[
b_A^{\rm determinant}+J_{Ai}c_i+
b_A^{\rm compensator,link,state}=0.
\]

All b_A and J_Ai entries are stored. No coefficient is fitted to the imposed
neck. This advances the missing source calculation while preserving the next
task: obtain the remaining finite completion and coupled state/geometry from
the same action. Anomaly matching alone does not select that invariant part;
see [Andrianov, Kurkov and Lizzi](https://arxiv.org/html/1106.3263v1).

## What the Schur signs do and do not imply

The child inverse is a reciprocal operator, not automatically a second minus
sign. At real energy above or below a child's spectrum, its inverse can have
different signs; at complex energy it is generally complex. A consistent
change of child spin frame transports both K_c and B and leaves
B K_c^-1 B† unchanged.

The record checks these statements with a noncommuting complex example. The
same coupled Hamiltonian has a nonzero Schur self-energy and a nonzero current
operator, yet its stationary vacuum has zero net regional power. Thus the
Schur subtraction is not an energy drain by itself. The physical sign of an
expansion response must follow from the derived stress and metric dynamics.
For a flat isotropic Einstein background, acceleration involves rho+3p;
the present static anisotropic neck is not that cosmological solution.

This preserves the constructive boundary hypothesis while replacing a sign
analogy with an explicit route: operator, state, metric variation, conserved
source, and coupled evolution. The present result computes the static
determinant part of that route. It does not claim the full compensator,
retarded metric response, physical mass/scale selection or final cosmology.
