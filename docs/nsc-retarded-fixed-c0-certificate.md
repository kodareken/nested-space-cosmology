# One retarded parent extension fails fixed-C0 linear matching

For the declared positive radius pulse and group14_1 at
`E=0.5562120090641313`, the combined certificate establishes

\[
\boxed{\delta\widehat C(E,E)_{00}
<-\frac{69}{875000000}<-7.8857\times10^{-8}.}
\]

This is a nonzero **full coherent** incoming covariance response. The
specified tangent therefore fails the fixed-C0 linear matching condition.
It does not exclude other parent extensions or establish non-existence for
the NSC architecture. The [record](../results/development/nsc-retarded-fixed-c0-certificate.json)
keeps that scope separate from global constraints and full stationarity.

## The previously missing bound

The [directed low-energy projector certificate](nsc-incoming-low-energy-vacuum-enclosure.md)
encloses the exact affine-vacuum component at rho1, with initialization,
whole-step transport, arithmetic and endpoint uncertainties included.
Its lower bound on `Re Pv01` exceeds1/4. This is the selected low energy's
own certificate; no group32 high-energy error was transferred to it.

The [analytic sign bound](nsc-retarded-diagonal-sign-bound.md) then applies.
The actual source remains the original coherent horizon pair with a closed
infinity channel. Its distance from the dominant vacuum component is exactly
`sqrt(f_H)` in operator norm. That entire remainder, including any allowed
reflection/relative phase, is included in the response allowance. The
vacuum component is a mathematical decomposition, not a replacement state.

The argument uses the existing unitary Dirac generator and the positive
short-shell insertion. Frozen input inequalities imply `norm(G)<3`,
`norm(K-K0)<.18 I`, and a full coherent allowance below`.002 I`.
Together with `Re Pv01>1/4`, the diagonal response is below`-.138 I`.
The original pulse has the exact positive insertion lower bound
`I>1/1750000`, yielding the displayed rational result.

## Why this is a complete selected-pair test

The [Fourier matching identity](nsc-incoming-fourier-matching.md) integrates
the stationary source-energy delta functions analytically at first order.
This selected pair therefore has no omitted energy-quadrature contribution.
All three source columns and horizon coherence are retained. Different
orthogonal angular blocks cannot cancel this block's nonzero kernel.

The earlier spacetime and short-shell numerical calculations locate the
response near`-1.05548e-5`. Their tolerance/refinement comparisons are not
used as proof inputs for this sign. The exact inequality instead combines
the new directed projector enclosure with conservative analytic bounds.
The compilation replays authenticated accounting and rational inequalities;
it performs no additional field propagation or source integration.

## Consequence for the constructive model

Holding the intrinsic surface and normal fixed does not by itself hold the
prepared incoming covariance fixed. This particular retarded extension of
the compatible normal data changes that state already at first order.
A parent extension intended to realize the conditional fixed-C0 constraints
must satisfy the quantum matching equation as well as those local geometric
equations. A different extension changes the question and must be tested
with the same source law and action, not supplied with a cancelling force.

The result is local in variation space. It does not claim matching fails
for every direction, that no isolated finite-amplitude return is possible,
or that a global constraint/metric solution has been ruled out. Full
extended stationarity remains OPEN; no metric equation is stepped.

```sh
python3 scripts/derive_nsc_retarded_fixed_c0_certificate.py --check
python3 -m pytest -q tests/test_nsc_retarded_fixed_c0_certificate.py
```

## Retained response methods and their distinct evidence

The spacetime Fourier-pair pilot kept the complete coherent source in
$\delta\widehat C_{oi}=\delta A_{oi}C_i f_i^\dagger
+f_oC_o\delta A_{io}^\dagger$ and analytically integrated the stationary
source-energy delta functions. At the selected diagonal energy it measured
`-1.0554769134063148e-5`, with spatial/time differences
`6.91800e-9`, `2.46272e-11` and CAR-pair residual `3.58011e-11`.
These are numerical diagnostics; the directed projector bound above closes
the later sign certificate rather than promoting those differences to proof.

The independent short-shell check uses

$$
A_i'=G_iA_i,\qquad
B_{oi}'=G_oB_{oi}
-\frac{i\ell s\chi(s)\widehat w(E_o-E_i)}{ar^2}\sigma_2A_i,
\qquad
G_E=\frac{i}{a}(-m\sigma_1+\ell\sigma_2/r-E\sigma_3/a).
$$

Here $s=T(\rho)-T(1)$ and $\widehat w(\omega)=\int e^{i\omega z}w(z)dz$.
Its upstream tangent is zero; all three source columns and four signed
energies remain. No $2\pi$, energy weight or angular multiplicity enters
this pair ODE. The conditional sign proof uses
$\|C_\Sigma-P_{v,\Sigma}\|=\sqrt{f_H}$,
$\|K-K_0\|<9I/50$, $K_0=-iI\sigma_2$ and
$[K_0,P_v]_{00}=-2I\operatorname{Re}(P_v)_{01}$.
Thus the full coherent remainder, including its phase, is bounded without
replacing the physical state. The selected-pair exclusion remains local;
it cannot exclude all histories or prove full action stationarity.
