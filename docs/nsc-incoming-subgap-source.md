# Group13 local subgap source on the incoming KS slice

The refined group13 contribution on $\rho=1$, in order
$(\rho,p_\parallel,T_{01},p_\perp)$, is approximately

$$
(1.19132077\times10^{-5},\;-5.12197861\times10^{-5},\;0,\;0).
$$

Relative to the inherited eight-node source panel, the density and axial
pressure change by $+3.33505\times10^{-11}$ and $-7.60184\times10^{-10}$.
The largest nested-nine versus seventeen-node difference in this physical
contribution is $4.67\times10^{-14}$. This is a refinement of one finite
spectral panel, not a complete 33-group source or constraint solution.

The [record](../results/development/nsc-incoming-subgap-source.json) retains
the full numerical vectors, raw current residual, independent analytic
control, quadrature comparisons and authenticated payload.

## Existing physical input, new local observable

Reuse the [group13 history-response input](nsc-subgap-history-response.md):
seventeen real Chebyshev unsewn basis nodes on $1\le E\le\pi/2$, the
separately prepared complex-frequency control and its dual, the normalization
control, and all 360 exact contour reflection samples. The new work continues
these columns only from $\rho=0$ to $\rho=1$ using the existing stationary
radial Dirac equation. It performs nineteen short solves. It performs no
horizon, scattering or history solve.

The exact trace inverse and clock phase of the common KS owner give the
canonical basis $(v,w)$ at the actual incoming slice. On the real axis the
physical source columns are $(Rv,w,0)$, with the unchanged horizon covariance
and closed incoming infinity channel. No integrated packet covariance or
seed matrix is used as a field.

The four local vertices are

$$
V_\rho=-m\sigma_1-\frac E{a_1}\sigma_3,\quad
V_\parallel=-\frac E{a_1}\sigma_3,\quad
V_{01}=\frac E{a_1}I,\quad V_\perp=0,
\qquad m=\pi/2.
$$

At complex $z$, the bilinear is
$B_{ab}(z)=\Phi_a(\bar z)^\dagger V(z)\Phi_b(z)$: its bra is obtained
from the independently archived conjugate-frequency basis. These analytic
columns are never passed to a Gaussian state constructor.

## Integrate the local source with its own normalization

With $d=f_H-1/2$ and the inherited horizon coherence $s$, each positive-energy
moment is

$$
d(B_{vv}-B_{ww})+2\operatorname{Re}[-isR B_{wv}]
+\frac12\operatorname{Tr}V-\operatorname{Tr}(P_{\rm ad4}V).
$$

Only the reflected coherence is integrated along the existing pole-free
contour. The smooth bilinear is interpolated; the reflection coefficient is
always an exact archived-node lookup. The half-identity and existing
fourth-order subtraction remain on the real axis.

The result is multiplied once by `incoming_group_factor`: copies and group
degeneracy are included, together with the negative-frequency counterpart.
The earlier history-action factor $-1/\pi$ does not belong to this local
stress calculation.

The record compares seventeen and nested-nine interpolation nodes, contour
orders 24 and 48, and heights $\kappa/3$ and $\kappa/4$. It also compares the
local formula against the archived real eight-node physical columns and
their existing incoming kernels. Raw analytic-bilinear residuals have their
own tolerance, separate from physical integrated-moment comparisons.

## Current identity and scope

The already verified sewing identity gives
$\operatorname{Tr}(S C_{\rm src}S^\dagger)-1
=\mathcal T(n_{\rm in}-f_H)$. Here $\mathcal T=0$ throughout the open subgap
interval. Therefore this panel's physical $T_{01}$ is exactly zero. The raw
integral gives about $9.85\times10^{-16}$ after the physical group factor;
that value remains in the payload and record as a numerical zero-identity
residual. It is not relabeled as physical flux or edited to zero.

No light restoration or compact induced term is included in this source
panel. The other 37 compact signed-family subgap refinements and the complete
source convergence remain separate. The measured controls do not constitute
a rigorous uniform spectral remainder bound.

```sh
python3 scripts/derive_nsc_incoming_subgap_source.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_subgap_source.py
```

Both commands read the prepared content-addressed artifact. Neither repeats
the nineteen short continuations or invokes a historical generator. The
verifier compares every recorded field, hashes and scope declaration.
