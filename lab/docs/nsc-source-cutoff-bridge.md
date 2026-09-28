# Fixed-source cutoff and the owned differentiated kernel limit

For the declared smooth compact pure-radius histories, with unchanged
intrinsic incoming metric and upstream preparation, the complete source
change in the common-kernel prescription is

\[
\boxed{\Delta G_B^{\rm owned}[g]
=\Delta G_B^{\rm raw,\Lambda}[g]+\Delta G_B^{\rm edge}[g]
+\mathcal T_{B,\Lambda}[g].}
\]

The first term is the actual finite-source, coincidence-first column
contraction. For each original angular family of multiplicity mu, the
source-cutoff edge term is

\[
\boxed{\Delta G_N^{\rm edge}
=\frac{\mu}{2\pi a_\Sigma}(f'_+-f'_-),\qquad
\Delta G_\beta^{\rm edge}
=-\frac{\mu}{2\pi}(f'_++f'_-).}
\]

The functions f_s are the
[actual Dirac phase coefficients](nsc-dirac-source-phase-transport.md),
with E=s*e and e>0. Each angular sign has its original ledger multiplicity.
The source measure is de/(2*pi). The omitted coincidence-first tail
`T_{B,Lambda}` is O(1/Lambda); its coefficient at the actual panel ends
160/320 has **not** been numerically enclosed. Thus the identity supplies
no finite-cutoff accuracy certificate and no local gate result.

This term follows from evaluating the same source-fixed state in the owned
order of limits. It is not a new action, source, coupling or `Gamma_rest`.
The already certified local/reference coefficients and paired unweighted
band-bulk cancellation are unchanged. Historical finite-source records
retain their original raw values.

## Scope and reused inputs

The history has `a=a_ref(rho)>0`, real `r_g>0`, and a smooth compact radius
perturbation, flat near the unchanged upstream slice. On Sigma,
`r_g=r_ref` as a function of z. The fixed incoming identification, N and
beta are unchanged. The local backward cone fits in the already verified
[periodic embedding](nsc-ks-local-embedding.md). No global parent/child
matching is imposed.

The [fixed-transfer source owner](nsc-incoming-fixed-transfer.md) supplies
the affine-vacuum projector expansion through more orders than needed here
and the exponential full-coherent-source difference. Only its unchanged
upstream expansion is reused. Its affine error constants are not asserted
to be the error constants of a changed history.

## A finite expansion gives a uniform source-frequency remainder

Set `V=-m*S1+ell*S2/r_g`, `Pi_s=(1+s*S3)/2`, and

\[
L_0=\partial_\rho-\frac{S_3}{a^2}\partial_z-\frac{iV}{a},
\quad \Theta'=a^{-2},\quad
X^{(M)}=e^{-ie\Theta}\sum_{j=0}^M e^{-j}A_j.
\]

After removing the common phase, the exact envelope operator is
`L0-2*i*e*Pi_{-s}/a²`. Construct the coefficients using

\[
\Pi_{-s}A_0=0,\qquad \Pi_sL_0A_j=0,\qquad
\Pi_{-s}A_{j+1}=\frac{a^2}{2i}\Pi_{-s}L_0A_j.
\]

The minor coefficient is algebraic in the previous coefficient and its
derivatives. The major coefficient solves transport on the fixed
characteristics. Its upstream value is taken from the unchanged normalized
affine-vacuum expansion. Since g equals the reference on an upstream
neighborhood, the minor compatibility conditions there are the reference
conditions. Normalizing this mathematical vacuum column fixes its phase
convention; it does not normalize or replace the physical source columns.

For every fixed finite M, the smooth coefficients and their needed
derivatives have finite norms on the compact preparation slab. Direct
substitution telescopes exactly:

\[
\left(L_0-\frac{2ie}{a^2}\Pi_{-s}\right)
\sum_{j=0}^M e^{-j}A_j=e^{-M}L_0A_M.
\]

The original envelope generator is skew-adjoint in the periodic L2 norm.
Its energy-dependent term is spatially constant and skew-Hermitian, so it
contributes nothing to the H2 energy inequalities either. Commutators
involve only V_z,V_zz and the fixed a. The existing
[H2 error owner](../src/recursive_horizons/nsc_ks_residual_error.py)
therefore propagates the initial truncation and the displayed defect with
constants independent of e. There is a finite C_M such that

\[
\|X-X^{(M)}\|_{H^2}\le C_M e^{-M},\qquad e\ge e_*.
\]

Sobolev evaluation gives the same power for the envelope and its first
axial derivative. Restoring the source carrier `exp(-i*s*e*z)` costs at
most one power of e. The differentiated two-point kernel error is
O(e^{1-M}), uniformly on the compact pair domain. M>=3 is sufficient for
an integrable error; M=2 alone is not. The exponentially small coherent
source remainder remains integrable after the same propagation.

This proves finiteness of the uniform constant needed for the limit
identity. It does not evaluate C_M numerically or make the omitted tail
smaller than the local gate tolerance.

## Retain the actual two-point phase before taking coincidence

Write `z_+=Z+eta/2`, `z_-=Z-eta/2`, `Delta f=f_s(z_+)-f_s(z_-)`.
At Sigma the first minor spinor coefficient is the same constant as in
the reference. The trace of the changed occupied-column kernel, through
inverse-energy order two, is consequently

\[
e^{-ise\eta}\left[
\frac{i\Delta f}{e}
+\frac{\frac{s}{2}(f'_s(z_+)+f'_s(z_-))
-\frac12(\Delta f)^2+i\Delta h_s}{e^2}
\right].
\]

Here h_s is the imaginary second major coefficient; its precise value is
not needed for the exchange. The same trace with S3 inserted is s times
this expression at these orders. The remaining differentiated kernel is
integrable by the finite expansion above.

For positive e_*, define the oscillatory integrals
`I_j(eta)=int_{e_*}^infinity exp(-i*s*e*eta)/e^j de`. At nonzero eta,
I1 has its usual convergent oscillatory meaning. Differentiation is in
the off-diagonal kernel/distribution sense; Abel damping gives the same
derivative. No physical spectral cutoff is added. Integration by parts gives

\[
I_2=\frac{e^{-ise_*\eta}}{e_*}-is\eta I_1.
\]

The potentially noncommuting leading terms therefore combine as

\[
i\eta f'_s(Z)I_1+s f'_s(Z)I_2
=\frac{s f'_s(Z)}{e_*}e^{-ise_*\eta}.
\]

Their differentiated coincidence current is `-f'_s(Z)/(2*pi)`, while
their current at coincidence before integration is zero pointwise in e.
This is the finite edge term. It is independent of the arbitrary split e_*.

For smooth f, `Delta f-eta*f'(Z)=O(eta^3)` and the averaged derivative
minus f'(Z) is O(eta^2). Since `I1=O(1+|log|eta||)` and
`I1'=O(1/|eta|)`, these differences contribute zero to the limiting
current. The `(Delta f)^2 I2` term also contributes zero. The
`i*Delta h*I2` term commutes with the limit: its surviving derivative is
the ordinary integrable `h'(Z)/e²` current. Thus linearizing the phase
only after retaining these terms loses no finite exchange.

The bounded off-diagonal N vertex has no changed order-e^-1 contribution
at Sigma, because the intrinsic radius is fixed there. Its remaining
kernel is integrable without a differentiated coincidence exchange.
The owned momentum vertices are `S3/a` for N and `-I` for beta; action
gradients carry the common minus sign. Contracting the current edge term
therefore gives the boxed N,beta signs.

Finally, after the order-e^-1 coincidence cancellation, the raw changed
current is O(e^-2). Its tail after a finite source cutoff is O(1/Lambda).
Computing a rigorous coefficient for that tail, and its history derivative,
is the remaining numerical task before this bridge can certify a gate.

The [executable algebra check](../scripts/check_nsc_source_cutoff_bridge.py)
verifies the finite recurrence at orders3 and4, the kernel primitive, the
two limit orders and both action-vertex signs. Its record binds the reused
source and embedding owners. The analytic uniformity argument above is
part of the proof; an exact symbolic identity alone would not supply it.

```sh
python scripts/check_nsc_source_cutoff_bridge.py --check
```
