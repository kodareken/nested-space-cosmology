# Both energy measures in the linear large-transfer bound

The owned Fourier response has source E_i, output E_o and omega=E_o-E_i.
Its inverse kernel uses dE_i dE_o/(2pi)^2, and the symmetric momentum
insertion is -(E_o+E_i)/2. Counting only the omega integral gives the wrong
UV power. This owner connects the existing short-shell response and raw
N,beta vertices with both measures retained; it is not a full-source bound.

For the reference radius family, unitarity of the owned G_E gives

$$
\|B_{oi}\|\le |\ell|\left(I_w|\widehat w(\omega)|+
I_U|\widehat U(\omega)|\right),\quad
I_w\le\frac{\sigma^2}{2r_{\min}^2},\quad
I_U\le\frac{\sigma^4}{24r_{\min}^2}.
$$

At sigma=.03 and reference r^2>=2, these bounds are exactly9/40000 and
27/1600000000. Canonical reference fibers have operator norm one and the
source covariance has norm at most one. Retaining both coherent response
terms gives ||delta Chat||<=2|ell|(I_w|w_hat|+I_U|U_hat|). The Pauli trace
bound adds the spin dimension2; it does not add an angular or energy factor.

For p>3 and compact profiles with certified pth-derivative L1 norms, let
D_p=I_w||w^(p)||_1+I_U||U^(p)||_1. Integration by parts gives the
omega^-p bound. On |E_i|>=E0, |omega|>=|E_i|/2, define

$$
A_p=\frac{2^{p+1}}{(p-1)(p-2)},\qquad
B_p=\frac{2^p}{p-3}\left(\frac4{p-1}+\frac1{p-2}\right).
$$

Both E_i signs and both omega tails are included. The density integral is
bounded by A_p E0^(2-p), and the |E_o+E_i| integral by B_p E0^(3-p).
In particular B4=88/3 and B6=112/5. W^(4,1) therefore gives E0^-1 for
the linear current tail; E0^-3 requires W^(6,1). An upper-bound power
count is not a lower-bound counterexample for any particular smooth profile.

The resulting raw linear matter bounds are

$$
|\delta G_N|_{\cal T}\le\frac{\mu|\ell|D_p}{(2\pi)^2}
\left[4(m+|\ell|/r)A_pE_0^{2-p}
+\frac2a B_pE_0^{3-p}\right],
$$
$$
|\delta G_\beta|_{\cal T}\le
\frac{2\mu|\ell|D_p}{(2\pi)^2} B_pE_0^{3-p}.
$$

The implementation uses exact rational arithmetic and the conservative
pi>3 to enclose the normalization. Its inputs must be actual upper/lower
bounds; it does not certify profile norms from sample values. The same
signed-family multiplicity mu is used once.

This concerns the linear response about the reference. It neither bounds
the complementary frequency region nor the finite-history covariance
remainder. The latter must preserve cancellations among all terms of
F_g C_src F_gdagger; a separate wave-amplitude norm can lose phase
cancellations. No old frozen-C0 result, full-source accuracy or physical
local gate is claimed. Reference/action subtraction remains unchanged.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_local_transfer_bound.py
```
