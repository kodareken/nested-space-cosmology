# Continuous source-error bounds from the exact local propagator

This owner carries a bounded preparation covariance error into the actual
history-minus-reference N/beta action. It uses the existing pure-radius
evolution equation and analytic radius bounds. It does not require a numerical
field trajectory or assume that the inhomogeneous propagator is pointwise
unitary. It covers only the source rows whose covariance errors are supplied.

## Characteristic row bounds

Factor the incoming carrier as $F=e^{-iEz}G$. The implemented envelope equation
has principal part $\sigma_3\partial_z/a^2$, a diagonal imaginary energy term,
and off-diagonal potential

$$\frac{i}{a}V,\qquad V=-m\sigma_1+\frac{\ell}{r}\sigma_2.$$

For the canonical two-column propagator the initial matrix is I. Along each
characteristic the energy term is a phase. The two row norms therefore obey a
coupled integral inequality with off-diagonal coefficient at most
$b=(m+|\ell|/r_{\min})/a_{\min}$. If $B$ bounds its time integral, each row is
bounded by $e^B$. This is a supremum over space and over the preceding time
interval. The homogeneous reference propagator is unitary; its row norms are 1.

The existing [radius coupling bounds](nsc-ks-radius-coupling-bounds.md) enclose

$$M\ge\int\|\delta V/a\|\,|d\rho|,
\qquad M_z\ge\int\|\partial_z V/a\|\,|d\rho|.$$

Writing $D=G-G_0$, its forcing is $(i\delta V/a)G_0$. Gronwall gives a row
bound $e^B M$. Differentiating G in z gives forcing $(iV_z/a)G$ and zero
initial derivative. The two growth factors combine to $e^B$, giving a row
bound $e^B M_z$. Thus, in Frobenius norm,

$$f=\sqrt2,\qquad d=\sqrt2 e^B M,\qquad d_{z,0}=\sqrt2 e^B M_z.$$

The physical derivative includes the carrier:

$$\|F_0{}_z\|_F=|E|f,\qquad
\|(F-F_0)_z\|_F\le |E|d+d_{z,0}.$$

The second term is essential; it is not a replacement of evolved momentum
by the source energy. The source lies outside the changed history: on the
owned slab a<1, and rho_up-1>=normal_support guarantees the preparation is
before the support. Positive radius is checked by the existing exact-rational
history bounds. These inequalities apply over the whole spatial line, hence
continuously on the incoming interval I.

## From covariance errors to action errors

Let $\epsilon_E$ bound the operator error of the unweighted two-spinor
preparation covariance at one energy. Field linearity gives the error
$\operatorname{Tr}(\delta Q\,\Delta H_B)$, with the original local vertices.
The trace-norm inequality $\|X^\dagger Y\|_1\le\|X\|_F\|Y\|_F$ applies
after forming the history-minus-reference difference. Define

$$A=2fd+d^2,\qquad J_0=(f+d)d_{z,0},$$

$$S_0=\sum_E\frac{w_E}{2\pi}\epsilon_E,
\qquad S_1=\sum_E\frac{w_E}{2\pi}|E|\epsilon_E.$$

For the declared multiplicity mu, the continuous source-error contributions
are bounded by

$$\epsilon_\beta\le\mu(AS_1+J_0S_0),$$

$$\epsilon_N\le\mu\left[\sqrt{m^2+\ell^2/r_\Sigma^2}\,AS_0
+\frac{AS_1+J_0S_0}{a_\Sigma}\right].$$

Both energy signs are summed explicitly with their own covariance bound.
The original frequency weight is applied once; no extra three-column factor
is present. Intrinsic radius and axial scale use certified lower bounds.

These are bounds for the exact response to the preparation discrepancy.
Numerical field error at the archived preparation belongs to its separate
component. In contrast, an insertion based on numerical endpoint norms must
enlarge those norms by their field errors, as the
[endpoint insertion owner](nsc-ks-covariance-insertion.md) does.

## Verification and remaining scope

An independent matrix exponential of a spatially varying Dirac envelope checks
the row/difference bounds and the N/beta contraction. Dropping the envelope
derivative fails that control. A zero radius history gives zero source error
in the difference; it does not erase baseline source uncertainty.

The record aggregates only authenticated preparation bounds for selected
group-14 rows. Missing rows, the other positive angular family, other groups,
energy quadrature and the remaining physical error components stay open.
The aggregate is not a complete upstream-budget certificate and does not
change the candidate residual or the physical source.

Owners: `src/recursive_horizons/nsc_ks_source_operator_majorant.py`, its tests,
and `scripts/derive_nsc_ks_source_operator_majorant.py`.
