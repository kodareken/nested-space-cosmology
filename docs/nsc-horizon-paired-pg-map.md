# Horizon preparation in the magnetic paired basis

This calculation supplies the magnetic angular normalization, determines the
angular cross block of the declared boundary state, and evaluates signed
massive horizon-to-seed mode maps without using seed covariance as input.
It also identifies a finite-collar normalization error in the compact seed
generator. The global C1b-H gate remains **OPEN**; the old records and PDF
v0.25.0 are preserved.

## 1. Magnetic eigenspaces, rather than multiplicity as a map

Import the standard spin-c monopole harmonic construction and apply it to
the locked flux $q=4$. Its component weights and spin are

$$
s_a=\frac{q-1}{2},\qquad s_b=\frac{q+1}{2},\qquad
j=\frac{|q|-1}{2}+n,\qquad \lambda_n^2=n(n+|q|).
$$

The implemented north/south sections use the same normalized Wigner matrix
elements. Their two-spinor transition is

$$
\chi_N=e^{iq\phi}
\mathrm{diag}(e^{-i\phi},e^{i\phi})\chi_S.
$$

The first factor is the magnetic gauge transition; the second is the spin
frame transition. Direct sphere quadrature checks the paired spinors
$(Y_{s_a},\eta Y_{s_b})/\sqrt2$, their eigenvalue $\eta\lambda_n$ and their
orthonormality. The zero angular level keeps only its surviving component.
No angular pair is invented for that level.

The conjugate charge bundle has Chern number $-q$ in the same fixed magnetic
background. With the recorded phase convention for its angular eigenspinors,
angular $\sigma_1K$ and radial $\sigma_3$ reproduce the locked
$\eta_1\otimes\sigma_3K$ map up to a known magnetic-index phase. This is
an explicit current/gauge identification, not a change of flux sector.

## 2. What fixes Z

The already adopted affine-horizon vacuum is positive affine frequency along
the null generators with the identity kernel on their angular spinor bundle.
The inherited incoming occupation is also angular diagonal. Projection into
the orthonormal magnetic eigenspaces therefore gives

$$
C_{\mathrm{boundary}}(E)
=I_{\mathrm{angular}}\otimes\bigl(C_H(E)\oplus n_{\mathrm{in}}(E)\bigr).
$$

This is a boundary preparation statement, not a product state on a PG slice.
Its nonzero horizon-partner correlations remain inside each $C_H$ block.
The spherical magnetic operator and transparent radial domain conserve
$\eta$. Hence distinct angular eigenmodes have no preparation cross term,
and propagation gives $Z(E)=0$ in this declared realization. Direct angular
overlaps verify the zero. It is not selected by CAR or adjusted to a stress
target. A different angular boundary state would be a different preparation.

## 3. Y comes from its own massive mode map

For each angular sign, the original scattering owner supplies complex
$R_\eta(E)$ and transmission probability $T_\eta(E)$. The input map is

$$
M_\eta=
\begin{pmatrix}0&1&0\\R_\eta&0&\sqrt{T_\eta}\end{pmatrix},
\qquad M_\eta M_\eta^\dagger=I_2.
$$

The positive square root fixes an incoming-basis phase. The incoming state
is diagonal and uncorrelated with the horizon pair, so that phase does not
select a new state. The complex reflection phase is retained.

The new owner propagates both signs with the existing massive Dirac equation
on the fixed geometry. In the stored logarithmic collar coordinate
$y=\log\delta_q$,

$$
H_y=\frac{\delta_q}{\sqrt W}
\bigl(-m r\sigma_1+\eta\lambda\sigma_2\bigr)
-\frac{E\delta_q}{W}\sigma_3.
$$

No position-dependent mass-axis rotation is used. Both noncommuting mass
terms remain in the evolution. The initial horizon frame, its constant
working-frame rotation and $M_\eta$ give a mode map $F_\eta$. Then

$$
X=F_+\bigl(C_H\oplus n_{\mathrm{in}}\bigr)F_+^\dagger,
\qquad
Y=F_-\bigl(C_H\oplus n_{\mathrm{in}}\bigr)F_-^\dagger.
$$

The source matrix is formed from the original occupation laws, including the
locked compact incoming threshold. The recorded seed matrices are read only
after this calculation, for comparison. The signed-frequency source relation
is checked directly on those occupation laws and combined with the magnetic
charge map; it is not inferred from an arbitrary completion of X.

The same-surface T map and the static-coordinate relation
$t=\tau-F(\rho)$, $F'=\beta/A$, convert the mode columns into the PG current
frame. Its differentiated expression agrees with the locked PG Dirac
generator, including the coframe/half-density derivative. The numerical
mode chart evaluated here is the trapped probe interval $[-1,1]$.
It is not a completed global spectral resolution or a spatial Cauchy covariance.

## 4. The finite compact collar does not use the same distance

The inherited coordinate conversion is

$$
\epsilon_\rho=|\rho-\rho_h|=r_h^2\delta_q+O(\delta_q^2),
\qquad W=2\kappa_h\delta_q+O(\delta_q^2).
$$

Consequently the dimensionless massive Frobenius argument is

$$
z_H=\sqrt{\lambda^2+(m r_h)^2}\sqrt{\frac{2\delta_q}{\kappa_h}}
=\frac{\sqrt{\lambda^2+(m r_h)^2}}{r_h}
\sqrt{\frac{2\epsilon_\rho}{\kappa_h}}+O(\delta_q^{3/2}).
$$

The old compact **interior** initializer uses the right-hand $1/r_h$ factor
but supplies $\delta_q$ rather than $\epsilon_\rho$. Its initial argument is
therefore too small by $1/r_h$. The angular initializer and the earlier
`ParentDirac.initial_covariance` use the matched interior convention.
The exterior reflection initializer correctly uses a radial offset and its
$1/r_h$ factor is retained.

The new record tests $\partial_yF+iH_yF$ and separately propagates matched
and legacy compact frames with the same scattering and source data. This
isolates the finite-X discrepancy from integrator or occupation changes.
Both collar conventions approach the same affine-horizon limit as
$\delta_q\to0$; this is a finite-collar error, not a nonexistence result.
No old generator or seed record is overwritten.

## 5. Gate and scope

The magnetic basis and boundary angular pairing have evaluated residuals.
The record computes matched-collar X and Y, signed partners and a quadrature
check of Z at one existing frequency per massive group. This is a scoped
modal evaluation; it is not an all-frequency spatial covariance.

For the compact groups, the corrected normalization does not reproduce the
current finite seed X. The legacy frame does reproduce it, locating the
responsible owner precisely. Those corrected modal results are retained as
development counterparts; they are not silently substituted into the locked
retained-state checkpoint.

**C1b-H remains OPEN** on consistent adoption/reproduction of the compact
finite-collar preparation and on the globally normalized PG spectral mode
resolution. C1b-M2, full retained C1b and C2--C4 remain closed to execution.
The physical Weyl mismatch, stress, nulls and metric trajectory are not
assigned. The LLL and signed-preparation certificates remain intact.

See the [record](../results/development/nsc-horizon-paired-pg-map.json) for
per-group results, the exact normalization distinction and source hashes.

```sh
python3 scripts/derive_nsc_horizon_paired_pg_map.py --check
python3 -m pytest -q tests/test_nsc_monopole_paired_basis.py tests/test_nsc_paired_horizon_preparation.py
```
