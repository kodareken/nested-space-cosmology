# Fourth-order reference symbol on the transmitting history

The new owner retains spatial operator composition in the reference
subtraction. It extends the existing homogeneous Bloch representation to a
matrix-valued phase-space symbol. The initial physical covariance $C_0$ is
unchanged; this object is a formal subtraction symbol, not a new vacuum state.

The projector construction is imported from
[Panati, Spohn and Teufel, sections 2 and 4.4](https://arxiv.org/abs/math-ph/0201055).
The new calculation applies it to the project's canonical half-density
Hamiltonian, its paired angular/compact channels, and the already specified
KS/PG history. No general adiabatic theorem is claimed as an NSC result.

## Fixed operator and clock

In the same Pauli frame as the existing fourth-order reference,

$$
H_T=-\beta_K k I+N_K\left[-m\sigma_1+
\frac{\lambda}{r_K}\sigma_2+\frac{k}{a_K}\sigma_3\right].
$$

Symmetric ordering of the spatial differential operator is already included
in its Weyl symbol. Adding a separate imaginary derivative term would count
that ordering twice. The negative-normal-energy band has gap
$2\epsilon=2N_K\sqrt{m^2+\lambda^2/r_K^2+k^2/a_K^2}$; it is not defined by
the sign of the shift-displaced coordinate energy. The gap-zero point is
outside this expansion and is rejected without a denominator floor.

The formal relations are

$$
P\star P=P,\qquad i\varepsilon\partial_TP=[H_T,P]_\star,
\qquad P=\sum_{j=0}^4\varepsilon^jP_j+O(\varepsilon^5).
$$

$\varepsilon$ counts derivatives; it is not a fitted scale or parameter.
The star product includes axial-coordinate and canonical-momentum
derivatives. Coordinate $\partial_T$ is used with the lapse-weighted
Hamiltonian. Its homogeneous limit yields the already owned
$N_K^{-1}\partial_T$ Bloch recursion. Spatial shift gradients also generate
the $k(\partial_z\beta_K)\partial_k$ term, which cannot be obtained by
substituting a scalar normal derivative alone.

The supplied spacetime metric uses the fixed chart $T=T(\rho)$,
$z=\tau+S(\rho)$. Its Taylor data follow $\partial_T\rho=-a_0$ and
$\partial_TS=-\beta_0/a_0$. Thus KS time is not replaced with PG time.

## What is computed

The executable recurrence solves the diagonal part from idempotence and
the off-diagonal part from transport. It retains the off-diagonal
idempotence defect and diagonal transport defect as compatibility residuals;
neither is discarded by projection or symmetrization. Formal Hermiticity,
idempotence and transport are checked order by order.

The spatial symbol need not have pointwise trace one. Already

$$
\mathrm{Tr}\,P_1=-\frac12
\widehat h\cdot(\partial_z\widehat h\times\partial_k\widehat h)
=-\frac{m\lambda\,\partial_zr_K}
{2a_Kr_K^2\left(m^2+\lambda^2/r_K^2+k^2/a_K^2\right)^{3/2}}.
$$

The old Bloch-only form misses this scalar component on the varying
geometry. A truncated pointwise symbol is not tested as a Gaussian
covariance with eigenvalues in $[0,1]$; that condition belongs to the physical
state, which this calculation does not alter.

The homogeneous restriction is compared with all 1,904 authenticated stored
reference blocks: full compact/LLL matrices and the angular energy/pressure
projections actually recorded. The old state generators are not rerun.
The spatial calculation covers all 63 signed families with both signs of
one inherited control momentum per family. This is a local symbol check,
not a spectral integral or a new global history.

## Remaining same-action connection

The [record](../results/development/nsc-spatial-reference-symbol.json) separates
the formal-symbol gate from the still-open generating subtraction functional
and its finite regulator/local-allocation identity. Those are necessary
before combining this reference with the physical Gaussian variation and
the [computed local action](nsc-spherical-local-history.md).
No finite matching coefficient is set to zero, no extra
$\Gamma_{\mathrm{rest}}$ is introduced, and no physical stress or stationary
solution is inferred from formal projector residuals.

```sh
python3 scripts/derive_nsc_spatial_reference_symbol.py --prepare
python3 scripts/derive_nsc_spatial_reference_symbol.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_spatial_reference_symbol.py
```
