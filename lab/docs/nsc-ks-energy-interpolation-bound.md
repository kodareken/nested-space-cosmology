# Conditional energy interpolation bound for the KS column propagator

This is a remainder bound for interpolating the existing KS envelope in the
source-energy label. It does not replace the 49,372 retained positive-energy
rows, change `A_up`, `C_src`, weights, the action, or `Gamma_rest`, or close
the local gate on `I=S(1)+[.12,.18]`. It does not apply one `W` across mixed
energy labels; per-label `X=W(E)A_up(E)` is a separate owner.

## Reuse

| Input | Owner |
|---|---|
| Envelope law `F=e^{-iEz}X`, `X_rho=L_E X` | [KS source envelope](nsc-ks-source-envelope.md) |
| Joint reference plus difference | [difference envelope](nsc-ks-difference-envelope.md) |
| Off-diagonal growth `K0`, `B_z` integral `K1`, and `sqrt(2)` row-to-Frobenius | [characteristic error](nsc-ks-characteristic-error.md) |
| Absolute `d rho` orientation | [residual error](nsc-ks-residual-error.md) |
| Original `A_up`, `C_src`, weights | [source inventory](nsc-ks-source-inventory.md) |
| Directed dyadic uppers | `nsc_ks_ball_trajectory.exact_upper` |

The missing connection is a certified interpolation remainder in the energy
label of the **operator basis** `W`, in the same row-l2/max comparison as
the characteristic majorant. Interpolate `W` only.

## Operator basis, not a new state

`W` is the fundamental 2x2 solution of the same real-energy generator

\[
L_E=\frac{S_3}{a^2}\partial_z+\frac{i}{a}\Bigl(-m S_1+\frac{\ell}{r}S_2-\frac{E}{a}S_3\Bigr),
\]

with `W(\rho_{\mathrm{up}})=I`. The identity is an operator basis. It is not
a covariance and not a replacement of any source column. After a
single-energy evaluation,

\[
X=WA_{\mathrm{up}},\qquad
F=e^{-iEz}WA_{\mathrm{up}},\qquad
F_z=e^{-iEz}(W_z-iEW)A_{\mathrm{up}},
\]

with the original `A_up` at that same energy. History tangents
`Y_\alpha=\partial_\alpha W` convert by the same columns. This module does
not reconstruct mixed-label tables; a single `W` applied to every energy,
or an extra energy axis on `W`, is not a general source map.

The generator is affine in real `E`. Second energy derivatives of `L_E`
vanish. The `E` term is anti-Hermitian for real `E` and does not enter the
row-norm growth; complex energy is outside the scope.

## Characteristic row-l2/max norm

The comparison is the same as the local characteristic helper:

\[
\|W\|_{\mathrm{row}}
=\max_{z,\,\mathrm{spin}\,s}
\bigl\|W_{s,:}(z)\bigr\|_2,
\]

the Euclidean l2 of each two-entry spin row, then the maximum over the two
spins and over `z` on the backward cone of `I`. The two rows propagate on
opposite `S_3` characteristics (speeds `\pm 1/a^2` in `rho`). They are not
one `C^2` vector at a common `z`, so this is **not** the induced 2x2
spectral norm `\|W\|_{2\to 2}`.

The matrix with both rows `(1,0)` has row-norm `1` and spectral norm
`sqrt(2)`. Stacking opposite-characteristic snapshots from two different
`z` into one 2x2 likewise inflates the spectral norm while the row-l2/max_z
comparison stays `1`. Homogeneous 2x2 ODE toys cannot identify the two
norms; this bound does not claim a tighter induced-2 estimate for
two-speed transport.

Row-norm still satisfies `\|S_3\|_{\mathrm{row}}=1` and
`\|-m S_1+(\ell/r)S_2\|_{\mathrm{row}}=\sqrt{m^2+\ell^2/r^2}`. All integrals
use absolute `|d\rho|` from `\rho_{\mathrm{up}}` down to `\rho_\sigma=1`:

\[
D=\int\frac{|d\rho|}{a^2},\qquad
K_0=\int\frac{\sqrt{m^2+\ell^2/r^2}}{a}\,|d\rho|,
\]
\[
K_1=\int\|B_z\|_{\mathrm{row}}\,|d\rho|,\qquad
J_{0,\alpha}=\int\|\delta B_\alpha\|_{\mathrm{row}}\,|d\rho|,\qquad
J_{1,\alpha}=\int\|\delta B_{\alpha,z}\|_{\mathrm{row}}\,|d\rho|.
\]

These must be supplied as already certified nonnegative uppers. Passing them
into this helper does not certify them. `D` is the envelope continuum speed
distance. `K_0` and `K_1` are the characteristic growth integrals.

## Derivative bounds: simplex cancellation and the rectangle

Let `U(\rho,\sigma)` be the propagator of `L_E` acting on spin rows. The
characteristic comparison gives `\|U(\rho,\sigma)\|_{\mathrm{row}}\le\exp(\int_\rho^\sigma c_0)`
along the cone, independently of real `E`. Concatenating
`U(\rho,\sigma)W(\sigma)` therefore produces one factor `e^{K_0}`, not
`e^{2K_0}`.

**Energy derivatives of `W`.** Variation of constants / Duhamel with the
affine generator yields

\[
\partial_E^n W(\rho)
=n\int_{\rho_{\mathrm{up}}}^\rho U(\rho,\sigma)\,(\partial_E L)(\sigma)\,
\partial_E^{n-1}W(\sigma)\,d\sigma,
\]

and `\|\partial_E L\|_{\mathrm{row}}=1/a^2`. Absolute integrals and the
comparison `u_n(t)\le n\int_0^t e^{K(t)-K(s)}\mu(s)u_{n-1}(s)\,ds` with
`u_0\le e^{K}` and `\mu=1/a^2` give, by induction,

\[
u_n(t)\le e^{K(t)}\,n\int_0^{D(t)} u^{n-1}\,du
=e^{K(t)}D(t)^n.
\]

The Leibniz factor `n` cancels the ordered-simplex volume `D^n/n`. The bound
is `D^n`, not `n!\,D^n` and not `D^n/n!`.

**Axial derivatives.** `W(\rho_{\mathrm{up}})` is independent of `z`, so
`W_z` starts at 0. Differentiating in `z` inserts one `B_z` vertex (`a` is
independent of `z`). Running integrals satisfy

\[
\frac{d}{dt}\bigl(K_1 D^n\bigr)
=\|B_z\|_{\mathrm{row}}D^n+n K_1 D^{n-1}\mu,
\]

which is the Duhamel integrand after the growth factor is removed. Bounding
each source term by the global `K_1` and `D` *before* adding would
manufacture a spurious factor 2. The `n` energy insertions and the one
`B_z` insertion fill a single rectangle of volume `K_1 D^n`.

**History tangents.** `Y_\alpha` is the same Duhamel source with
`\delta B_\alpha` in place of `B_z`, and `Y(\rho_{\mathrm{up}})=0`. Mixed
`z`/history insertions are one `\delta B` and one `B_z` together with `n`
energy vertices. Their unordered volume is `J_0 K_1 D^n`. The direct vertex
`\delta B_z` supplies `J_1 D^n`.

Thus, in the row-l2/max_z norm, with a common `e^{K_0}`,

\[
\|\partial_E^n W\|_{\mathrm{row}}\le e^{K_0}D^n,
\qquad
\|\partial_z\partial_E^n W\|_{\mathrm{row}}\le e^{K_0}K_1 D^n,
\]
\[
\|\partial_E^n Y_\alpha\|_{\mathrm{row}}\le e^{K_0}J_{0,\alpha}D^n,
\qquad
\|\partial_z\partial_E^n Y_\alpha\|_{\mathrm{row}}\le e^{K_0}(J_{1,\alpha}+J_{0,\alpha}K_1)D^n.
\]

These are majorants for the characteristic comparison system, not identities
for a unitary `L^2` or spectral norm. The commuting exponential with
`K_0=0` saturates the `D^n` row-norm bound.

## Chebyshev-root interpolation remainder

Take `N+1` Chebyshev **roots** of `T_{N+1}` on `[E_c-h,E_c+h]`,

\[
E_j=E_c+h\cos\frac{(2j+1)\pi}{2(N+1)},\qquad j=0,\ldots,N.
\]

These are not Chebyshev--Lobatto extrema. On `[-1,1]` the nodal polynomial
is the monic Chebyshev polynomial `T_{N+1}(x)/2^{N}`
([Trefethen, ATAP, Chapter 3](https://doi.org/10.1137/1.9781611975840),
after (3.10)): leading coefficient `2^{N}` and
`|T_{N+1}(\cos\theta)|=|\cos((N+1)\theta)|\le 1`. Hence on `[E_c-h,E_c+h]`

\[
\bigl|\omega(E)\bigr|
=\Bigl|\prod_{j=0}^{N}(E-E_j)\Bigr|
\le\frac{h^{N+1}}{2^{N}}.
\]

The interpolation remainder in a Banach space `B` is algebraic,

\[
f(E)-p_N(E)=f[E_0,\ldots,E_N,E]\,\omega(E).
\]

Here `B` is the space of 2x2 blocks with the row-l2/max_z norm. A scalar
mean-value representation `f^{(N+1)}(\xi)/(N+1)!` is not used. The
Hermite--Genocchi formula
([Hermite, J. Reine Angew. Math. 84 (1878), 70--79](https://doi.org/10.1515/crll.1878.84.70))
writes the divided difference as a Bochner integral of `f^{(N+1)}` over the
simplex of volume `1/(N+1)!`. Therefore

\[
\bigl\|f[E_0,\ldots,E_N,E]\bigr\|_{\mathrm{row}}
\le\frac{1}{(N+1)!}\sup\|f^{(N+1)}\|_{\mathrm{row}},
\]

provided the convex hull of the nodes and `E` lies in the real interval,
which it does for every `E\in[E_c-h,E_c+h]`. Combining with the derivative
majorants gives the ideal interpolation remainders

\[
\|W-p_N\|_{\mathrm{row}}\le e^{K_0}\frac{(Dh)^{N+1}}{2^{N}(N+1)!},
\]

and the same factor times `K_1`, `J_{0,\alpha}`, and
`J_{1,\alpha}+J_{0,\alpha}K_1` for `W_z`, `Y_\alpha`, and `Y_{\alpha,z}`.

This is the interpolation remainder of exact nodal values of `W` in the
row-l2/max_z norm. It is not a field-solve certificate. Rounded-node
polynomial factors are a separate owner.

## Field conversion: sqrt(2) times row times weighted A_up

Let `a` bound the original weighted Frobenius norm of `A_up` at one energy,
and let `E_\star` bound `|E|` on the interpolation interval. For a spin row
`r` of `Delta W` and a weighted column `w_c A_{:,c}`,

\[
\bigl|(Delta W\,A)_{s,c}w_c\bigr|
\le\|r_s\|_2\,w_c\|A_{:,c}\|_2.
\]

Summing columns and the two spins,

\[
\|Delta F\|_{F,w}
\le\sqrt{\|r_0\|_2^2+\|r_1\|_2^2}\,a
\le\sqrt2\,\|Delta W\|_{\mathrm{row}}\,a.
\]

The unitary phase `e^{-iEz}` does not change the norm. The envelope
derivative `W_z-iEW` is a triangle inequality on the same rows, so

\[
\|Delta F\|_{F,w}
\le\sqrt2\,a\,\|Delta W\|_{\mathrm{row}},
\]
\[
\|Delta F_z\|_{F,w}
\le\sqrt2\,a\bigl(\|Delta W_z\|_{\mathrm{row}}+E_\star\|Delta W\|_{\mathrm{row}}\bigr),
\]

and the identical formulae with `Y` in place of `W`. Omitting `sqrt(2)`
fails on the matrix whose rows are both `(1,0)`: the row-norm is `1` while
the weighted Frobenius of `Delta W A_up` saturates `sqrt(2)\,a`. That is
the same `sqrt(2)` the characteristic helper uses to pass from two row
bounds to a weighted Frobenius field.

The identity matrix must not be substituted for `A_up` or `C_src`. Node-solve
errors in `W` at the Chebyshev roots, and rounding of those nodes and
barycentric weights, are a separate missing budget and are not set to zero.

## Scope that remains OPEN

- The local incoming gate on `I=S(1)+[.12,.18]` stays OPEN.
- Finite real energy only. No UV tail, no high-energy contraction constant,
  and no physical PASS.
- Source parameters and the action are fixed. Conditional integral inputs
  do not become certified by being used here.
- Node field-solve errors and rounded node/weight errors are missing.
- Mixed-energy reconstruction is not performed here.
- A completed cone residual, source-preparation error, and continuum
  quadrature error remain separate, as in the characteristic helper.

The implementation uses python-flint directed balls and exports the same
dyadic `{mantissa, exponent}` records as the other KS validation owners.
