# The actual Dirac source columns acquire a density transport coefficient

For the fixed-source radius history, the occupied high-frequency source
column has relative total density coefficient

\[
\boxed{\delta n_{2,s}(z)=s\,\partial_z f_s(z)},\qquad
f_s(z)=\frac{\ell^2}{2}\int_1^{\rho_u}
\left(r_g^{-2}-r_{\rm ref}^{-2}\right)
\bigl(\rho,z-sD(\rho)\bigr)\,d\rho,
\quad D(\rho)=\int_1^\rho\frac{d\rho'}{a(\rho')^2}.
\]

Here `E=s*e`, `e>0`, `s=+/-1`; this is the original source label. It is
not identified with outgoing Weyl momentum. The result is a **formal
asymptotic coefficient** derived from the actual envelope PDE, not an
error estimate or a new stress term. A uniform remainder is still required
to use it in the ordered source/coincidence limit.

## Continuity fixes the amplitude that a leading phase alone misses

The [owned Dirac envelope](../src/recursive_horizons/nsc_ks_signed_state.py)
obeys

\[
X_\rho=\frac{S_3}{a^2}X_z+\frac{iV}{a}X
-\frac{iE S_3}{a^2}X,\qquad V=-mS_1+\frac\ell r S_2.
\]

Because V is Hermitian and a depends only on rho, the exact total density
and spin flux obey
`partial_rho(X†X)=a^-2 partial_z(X†S3 X)`. The background high-energy
source approaches the occupied affine-vacuum projector, of trace one,
with the already owned exponential coherent source remainder. This reuses
the source assumptions of the [fixed-transfer owner](nsc-incoming-fixed-transfer.md);
it does not normalize the evolved modes or replace the fixed source.

Remove the common phase whose rho derivative is `-i*e/a²`. The principal
occupied component has `S3=s`. Algebraic slaving in the other component
gives

\[
X_{-s}=-\frac{a}{2e}V_{-s,s}X_s+O(e^{-2}),\qquad
|X_{-s}|^2=\frac{a^2(m^2+\ell^2/r^2)}{4e^2}+O(e^{-3}).
\]

Writing `n=1+n2/e²+...` therefore gives

\[
X^\dagger S_3X=s\left[1+\frac{n_2-a^2V^2/2}{e^2}\right]+O(e^{-3}),
\qquad
\left(\partial_\rho-\frac{s}{a^2}\partial_z\right)n_2
=-\frac{s}{2}\partial_z V^2.
\]

The characteristic through the incoming point is `z-s*D(rho)`.
The reference potential and the unchanged upstream preparation are
spatially homogeneous. Integrating the difference of the last equation
from rho_u down to1 proves the boxed coefficient, with the displayed sign.
The same major-component recurrence gives relative phase `exp(i*f_s/e)`.
Leading transport has unit Jacobian, but this does not remove the
next-order total-density response.

The [executable check](../scripts/check_nsc_dirac_source_phase_transport.py)
also obtains the density equation from the next minor component and the
major-component recurrence, for both source signs. It verifies the
anti-Hermitian generator, minor coefficient and density transport exactly.

## Consequence for the unresolved source limit

At coincidence the inverse-energy terms in the momentum current cancel:

\[
\operatorname{Im}(X^\dagger X_z)-E n:
\qquad \frac{\partial_z f_s}{e}
-\frac{s\,\delta n_{2,s}}e=0.
\]

This cancellation does not justify interchanging the source cutoff and
the differentiated coincidence limit. At nonzero separation the phase is
`f_s(z+)-f_s(z-)`, not exactly `eta*f_s'(z)`. The endpoint amplitudes,
spinor terms and a uniform remainder must be carried through the common
kernel before that limit is taken. The
[globally unitary control](nsc-source-cutoff-unitary-control.md) explains
why unitarity alone cannot supply the missing exchange argument.

No correction is inserted into either constraint by this result. The
finite-source sum, all certified local/reference coefficients and the
locked frozen-C0 regressions retain their recorded values. A numeric tail
bound and an actual local gate solution remain separate requirements.

```sh
python scripts/check_nsc_dirac_source_phase_transport.py --check
```
