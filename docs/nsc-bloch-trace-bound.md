# Conditional Bloch/Sobolev covariance-to-local-stress error propagator

This helper converts certified Sobolev/nuclear seminorms of the source-fixed
defect `D=C-P` into pointwise density and current traces, then into `N,beta`
insertion error. It is a missing finite-history UV norm link for the actual
source-sum/field certificate. It is not a physical gate, an existence claim,
or a new source. The root certificate remains the actual-source sum and field
bound; this owner only propagates already certified uppers.

Reuse: owned Weyl Bloch normalization, KS characteristic commutators, the
source-fixed law `C_g=U_g C_{\rm up} U_g^\dagger`, and standard Schatten
Hölder/Parseval. An earlier symbol-`L^1` nuclear bound is not used. The
spinor factor `\sqrt{\mathrm{spin}}` in the Hilbert--Schmidt Hölder is kept.

`P` is only an auxiliary add/subtract for the **same** `C_g`. Physical
covariance contributions are reconstructed as `P_g-P_{\rm ref}` plus
`D_g-D_{\rm ref}`. This does not change the action and does not assert that
`P` generates `\Gamma_{\rm sub}`. Band-action closure belongs to a generating
subtraction functional, if that functional is used at all; it is not a
covariance nuclear seminorm and is not an input of this propagator.

## Bloch fiber and Weyl matrix

Period `L`, Bloch momentum `q\in[-\pi/L,\pi/L]`, `k_n=q+2\pi n/L`. The fiber
inner product is `L^{-1}\int_{\rm cell}` with orthonormal basis
`e^{ik_n z}`. Cell Fourier coefficients of a Weyl symbol `a(z,k)` are

\[
\hat a_l(k)=L^{-1}\int_{\rm cell}e^{-i\omega_l z}a(z,k)\,dz,
\qquad\omega_l=2\pi l/L.
\]

The owned Weyl matrix is `A_{mn}(q)=\hat a_{m-n}((k_m+k_n)/2)`. The kernel is
reconstructed by `\int_{\rm BZ} dq/(2\pi)\sum_n` with **no extra** `1/L`.
Diagnostic mode lists must be distinct integers; non-integer or repeated
indices are rejected so the orthonormal-fiber identities are not applied to
an invalid basis. There is no physical periodic boundary condition: the
periodically padded coefficient problem agrees with the physical problem on
`I` only when a causal embedding is separately established.

## Period `L` and `1/L` as one enclosing interval

Supplied `L` is kept as a strictly positive enclosing interval in every
formula that contains `L` or `1/L`. Norms collapse upward. Denominators are
not replaced by their upper endpoints: `1/L_{\mathrm{upper}}` can underbound
the true reciprocal. The same rule applies to `a` and `r` in
`j/a` and `\ell/r`. Collapsed dyadic-upper records are accepted for seminorm
inputs via `restored_upper`, and rejected for `L`, `a`, and `r`.

On each fiber

\[
S(q)^2=C_0(q)^2=\sum_n(1+k_n^2)^{-1},\qquad
C_1(q)^2=\sum_n\frac{k_n^2}{(1+k_n^2)^2}.
\]

The `n=0` terms are at most `1` and `1/4`. For `n\neq0`,
`|k_n|\ge 2\pi(|n|-1/2)/L`. The odd-Basel sum
`\sum_{m\ge1}(m-1/2)^{-2}=\pi^2/2` then yields the uniform computational-period
bounds

\[
C_0^2=S^2\le 1+\frac{L^2}{4},\qquad
C_1^2\le\frac14+\frac{L^2}{4}.
\]

These are **not** the real-line values `1/2` and `1/4`, and not the large-`L`
sums `L/2` and `L/4`. At `L=2/5`, `L/2=1/5` is smaller than the `q=0` mode
contribution `1`, so `\sqrt{L/2}` is false as an upper bound. The Brillouin-zone
measure under `dq/(2\pi)` is the reciprocal of the same enclosing `L`.

## Hilbert--Schmidt Hölder, with the accepted spin correction

Let `H(q)=\|\Lambda^2 A\Lambda^2\|_{\mathrm{HS}}` on the spinor fiber. Schatten
Hölder and `\|\Lambda^{-1}\|_{\mathrm{op}}\le1` give

\[
d_{11}(q)=\|\Lambda A\Lambda\|_1\le\sqrt{s}\,H(q)S(q),\qquad
d_{21}(q)=\|\Lambda^2 A\Lambda\|_1\le\sqrt{s}\,H(q)S(q),
\]

with `s=2` for the owned KS spinors. The scalar-mode formula `H S` fails for
the `2\times2` identity. Point evaluation does **not** get another spin factor:

\[
\|\rho(z)\|_1\le d_{11}C_0^2,\qquad
\|j(z)\|_1\le d_{21}C_0 C_1.
\]

For Hermitian `D` the symmetric derivative is `(X+X^\dagger)/2` with
`\|(X+X^\dagger)/2\|_1\le\|X\|_1`.

## Parseval bound on `H`

\[
\int_{\rm BZ}\frac{dq}{2\pi}H(q)^2
=\sum_l\int_{\mathbb R}\frac{dk}{2\pi}
\langle k+\omega_l/2\rangle^4\langle k-\omega_l/2\rangle^4\|\hat a_l(k)\|_F^2
\le\frac8L I[a],
\]

\[
I[a]=\int_{\rm cell}dz\int_{\mathbb R}\frac{dk}{2\pi}
\Bigl(\langle k\rangle^8\|a\|_F^2+\frac1{256}\|\partial_z^4 a\|_F^2\Bigr).
\]

Cauchy--Schwarz against the zone measure `1/L` yields
`\int_{\rm BZ}H\le\sqrt{8I}/L`, using the same enclosing `L`. Instantaneous
`I[a]` is optional symbol-route input. Time-integrated residual seminorms
`\int d_{11}(R)` and `\int d_{21}(R)` are different objects and are not filled
from `I[a]`.

The required mixed geometry order for `\partial_z^4 R` is **9**: `P_4` has
four physical derivatives, its time derivative adds one, then `z^4` adds four,
and the kinetic first-star term adds one spatial derivative. Spatial-only
bounds are insufficient for time-containing jets. Formal `P_4` does not
supply `I[a]`.

## Duhamel under the same KS `U`

`D=C-P` obeys Duhamel. `R=\dot P-[G,P]` is the residual of the **actual**
subtracted Weyl operator. Reuse the owned commutator integrals
`K_1=\int\|B_z\|`, `K_2=\int\|B_{zz}\|`. Then `M_1\le 1+K_1` and
`M_2\le 1+\sqrt{4K_1^2+(K_2+K_1^2)^2}`, and

\[
d_{11}^{\rm ev}\le M_1^2\bigl(d_{11}(D_{\rm up})+\textstyle\int d_{11}(R)\bigr),
\qquad
d_{21}^{\rm ev}\le M_2 M_1\bigl(d_{21}(D_{\rm up})+\textstyle\int d_{21}(R)\bigr).
\]

Those evolved covariance seminorms are the local-trace inputs. No band-action
remainder is added.

## Local traces

\[
\|\rho\|_1\le C_0^2 d_{11}^{\rm ev},
\qquad
\|j\|_1\le C_0 C_1 d_{21}^{\rm ev},
\]
\[
N\le\mu\Bigl(\sqrt{m^2+\ell^2/r^2}\,\|\rho\|_1+\|j\|_1/a\Bigr),
\qquad
\beta\le\mu\|j\|_1.
\]

Missing genuine norm inputs remain `None`. `physical_local_gate` is always
`OPEN`.

```sh
PYTHONPATH=src /tmp/nsc-validated-runtime.8pg_7t1z/bin/python -m pytest -q tests/test_nsc_bloch_trace_bound.py
```
