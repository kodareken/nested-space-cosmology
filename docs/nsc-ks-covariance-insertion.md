# Local source-covariance error insertion

This owner converts a supplied upstream covariance error into a bound on the
history-minus-reference $N,\beta$ difference. It does not change the action,
the archived source, or the actual axial derivative. It is not an inheritance
identity, an infinite-nest limit, or a closed physical source.

The committed records `nsc-subgap-upstream-covariance-v1` and
`nsc-vacuum-source-remainder-v1` already store covariance-operator bounds.
The first is the archived $\rho_{\rm up}$ comparison for one subgap row,
including the signed Gram and complement discrepancy. The second covers 768
high-energy signed rows as a preparation bound. This module does not open
those files and does not rerun their drivers. A caller may pass those
epsilons. Passing them does not extend the proof.

## One fiber

At fixed energy the homogeneous upstream matrix $A$ is $2\times 3$ and

$$
Q=A C_{\rm src} A^\dagger.
$$

The same retarded solution acts on every column of that fiber, so the weighted
fields are $F=KA$ and $F_z=LA$. The right inverse

$$
R=A^\dagger(AA^\dagger)^{-1},\qquad \delta C=R\,\delta Q\,R^\dagger
$$

reproduces $\delta Q=A\,\delta C\,A^\dagger$ and leaves the original $C_{\rm src}$
and $A$ untouched. No column phase is aligned. The Gram is invariant under
diagonal column phases. Q is invariant only when the phase also commutes
with the source covariance or is compensated there. A common phase cancels;
relative phases that change a coherent Q are retained.

If $\lambda_{\min}(AA^\dagger)\ge\gamma>0$ and $\|\delta Q\|_{\rm op}\le\epsilon$,
then $\|\delta C\|_{\rm op}\le\epsilon/\gamma$. The factor is the reciprocal of
a directed lower endpoint of the smaller eigenvalue of the original binary
$2\times 2$ Gram. Replacing $\gamma$ by 1 underestimates $\delta C$ when the
fiber is near rank failure. An enclosure that is not strictly positive, an
interval that straddles an invalid sign, a missing $\epsilon$, and an
incomplete fiber are rejected. Absent bounds are rejected before Arb conversion.

Fibers are complete three-column blocks of one nonzero energy. Covariance
between different energies is not an insertion fiber and energies are not
folded. For energy-diagonal blocks,

$$
\eta=\max_e\epsilon_e/\gamma_e.
$$

A negative energy uses its own $(\epsilon,\gamma)$. Copying the positive scalar
is rejected. That own bound has to include the numerical Gram and complement
discrepancy; a zero Bloch residual does not cancel a non-unitary numerical Gram.
The Gram is not assumed to be the identity. The later inhomogeneous $K$ is not
treated as pointwise unitary, so a negative field norm is not copied from the
positive field either.

## History minus reference

The vertices are those of `source_column_matter`:

$$
V=-m\sigma_1+\frac{\ell}{r}\sigma_2,
$$

$$
H_N=-\mu\left[K^\dagger VK+\frac{K^\dagger\sigma_3 L-L^\dagger\sigma_3 K}{2ia}\right],
\qquad
H_\beta=\mu\frac{K^\dagger L-L^\dagger K}{2i}.
$$

The scalar contraction is $\operatorname{Tr}(QH)$. The source error in the
difference of two histories is

$$
\operatorname{Tr}\bigl(\delta Q(H[g]-H[{\rm ref}])\bigr).
$$

The intrinsic $a,r$ in $H$ are the declared slice. The norm bound uses lower
endpoints of those same quantities and does not replace an evolved $F_z$ by
$-EF$. Column weights already contain $\sqrt{dE/(2\pi)}$ and are not applied
again.

Let $f,f_z,d,d_z$ be certified weighted Frobenius uppers for the reference
field, its actual axial derivative, the difference $D=F-F_{\rm ref}$, and
$D_z$. Certified field errors are added to those uppers before the inequality,
which keeps the source-field cross term. With the enlarged norms,

$$
\text{density}\le\eta(2fd+d^2),\qquad
\text{momentum}\le\eta(fd_z+f_zd+dd_z),
$$

$$
N\le\mu\left[\|V\|\,\text{density}+\text{momentum}/a_{\min}\right],
\qquad
\beta\le\mu\,\text{momentum}.
$$

$\|V\|=\sqrt{m^2+\ell^2/r_{\min}^2}$. The same polynomial is evaluated by
`finite_matter_error` with source norm $\eta$. A zero enlarged difference
cancels in the history-minus-reference source error even when $\eta$ and the
reference field are positive. A numerical reconstruction of $K,L$ and of
$\delta C$ is only a validation helper. It is not a certificate.

The endpoint adapter reads the named norm and error slots and an explicit
caller scope. It does not treat a source-norm slot as $\epsilon$, does not
build a campaign binding, and does not populate a physical budget.
`KSUpstreamBatch` supplies the fiber layout and the unchanged $A_{\rm up}$ only.
Its preparation is not rerun.

## Scope

The caller must authenticate each supplied covariance bound against the same
source row, energy, channel, preparation and slice. A matching panel name alone
is not evidence. Both the subgap and high-energy preparation records can supply
such bounds within their stated coverage. This insertion adds no rows,
families, quadrature or changed-history field remainder. The full upstream gate
remains OPEN, and no physical upstream budget component is filled.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_covariance_insertion.py
```
