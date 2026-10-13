# Source-preserving finite-energy KS propagator

`nsc_ks_energy_propagator.py` reuses the existing joint reference/difference
Dirac evolution on two basis columns per energy. The resulting operator can
be applied to several disjoint source batches without repeating their field
evolution. The physical incoming law remains
`C_Sigma[g] = U_g C_up U_g†`.

For the canonical envelope, let `W_g(E,z)` evolve from the identity matrix
on the fixed upstream slice. Linearity gives

\[
 X_g(E,z)=W_g(E,z)A_{\rm up}(E),\qquad
 F_g(E,z)=e^{-iEz}W_g(E,z)A_{\rm up}(E).
\]

The identity is an operator basis. The original physical `A_up`, `C_src`,
source energies and quadrature weights are retained on application. Their
values are never interpolated. All three coherent columns at each original
energy are reconstructed, and their covariance is contracted by the existing
matter owner. The baseline and channel multiplicities still belong to the
constraint accumulator.

`evolve_energy_propagator(interval, degree, family, grid, target_z, mass,
angular, rho_up, **solver_options)` returns a `KSEnergyPropagator` containing
the reference, evolved difference, its axial derivative, and retarded history
tangents at the operator's energy nodes. The existing evolution is used without
changing its generator. Its fingerprint binds the sampled history and mesh;
the operator digest additionally binds all operator values, target positions,
energy nodes and channel parameters.

`operator.apply(source, A_up)` interpolates those operators and constructs a
`KSDifferenceIncoming` bound to that exact source and initial field. It forms

\[
 F_z=e^{-iEz}(W_z-iEW)A_{\rm up},\qquad
 \delta F_z=e^{-iEz}(\delta W_z-iE\delta W)A_{\rm up}.
\]

The input energy remains a label. The `W_z` terms are retained. Retarded
variations use the same evolution and fixed `A_up`; no source derivative is
introduced. Positive and negative energies can be handled separately by the
existing signed-state map.

The numerical interpolant passes through the actual rounded Chebyshev nodes.
Barycentric weights are calculated from those nodes. Extrapolation outside
the declared finite energy interval is rejected. Ideal Chebyshev remainder,
node evolution error, node displacement, interpolation arithmetic and source
preparation errors are separate requirements. Recovery of unphased tangents
from the existing field interface is also floating-point arithmetic and is
included in the unresolved reconstruction error. No formal polynomial order
or requested ODE tolerance supplies those bounds.

Focused synthetic tests compare all four fields and both matter insertions,
including the full retarded derivatives, against independent direct-column
evolution. They retain a non-diagonal source, exercise original interval
endpoints and detect dropped coherence. Different upstream columns reuse the
same operator while yielding different physical states. These are numerical
implementation controls. Full finite-source accuracy, the changed-history
UV remainder and the physical local incoming gate remain OPEN.


## Conditional energy interpolation and all-momentum auxiliary bounds

The continuum comparison norm is the maximum over space and spin of each
two-entry row's Euclidean norm, not the induced matrix spectral norm.
For already certified cone integrals

$$D=\int a^{-2}|d\rho|,\quad
K_0=\int\sqrt{m^2+\ell^2/r^2}\,a^{-1}|d\rho|,\quad
K_1=\int\|B_z\|_{\rm row}|d\rho|,$$

and $J_{0,\alpha}=\int\|\delta B_\alpha\|_{\rm row}|d\rho|$,
$J_{1,\alpha}=\int\|\delta B_{\alpha,z}\|_{\rm row}|d\rho|$,
Duhamel for the affine real-energy generator yields

$$\begin{aligned}
\|\partial_E^nW\|_{\rm row}&\le e^{K_0}D^n,&
\|\partial_z\partial_E^nW\|_{\rm row}&\le e^{K_0}K_1D^n,\\
\|\partial_E^nY_\alpha\|_{\rm row}&\le e^{K_0}J_{0,\alpha}D^n,&
\|\partial_z\partial_E^nY_\alpha\|_{\rm row}
&\le e^{K_0}(J_{1,\alpha}+J_{0,\alpha}K_1)D^n.
\end{aligned}$$

The ordered-simplex factor cancels the energy Leibniz factor; no $n!$ or
second growth exponential enters these derivative majorants. At $N+1$
Chebyshev roots on $[E_c-h,E_c+h]$, Hermite--Genocchi gives the ideal
operator interpolation remainder

$$\|W-p_N\|_{\rm row}\le
e^{K_0}\frac{(Dh)^{N+1}}{2^N(N+1)!}.$$

Multiply this factor by $K_1$, $J_{0,\alpha}$ or
$J_{1,\alpha}+J_{0,\alpha}K_1$ for the other fields. These are roots,
not Lobatto extrema, and this exact-node bound does not include node
evolution, node displacement or barycentric arithmetic. Applying the
operator to original weighted columns of Frobenius norm $a_{\rm up}$
requires the two-row factor:

$$\|\Delta F\|_{F,w}\le\sqrt2a_{\rm up}\|\Delta W\|_{\rm row},\qquad
\|\Delta F_z\|_{F,w}\le\sqrt2a_{\rm up}
\left(\|\Delta W_z\|_{\rm row}+E_*\|\Delta W\|_{\rm row}\right).$$

The same holds for history tangents. Interpolate the operator only; the
49,372 retained positive-energy rows keep their individual original
$A_{\rm up},C_{\rm src}$ and weights. Finite real-energy interpolation
does not supply the source UV tail or physical gate certificate.

Large axial transfers also need auxiliary-projector bounds valid for all
real canonical $k$, including shifted arguments crossing zero. On a
positive geometry box let
$c_{\min}=\inf1/a$ and
$g_{\min}=\sqrt{m_{\min}^2+(|\ell|_{\min}/r_{\max})^2}>0$.
For $H=ck\sigma_3+V$, $V=-m\sigma_1+(\ell/r)\sigma_2$ and
$g^2=c^2k^2+m^2+\ell^2/r^2$, a nonnegative Taylor jet $X$ majorizes
$g^2/g_0^2-1$, with zero constant. Its ordinary-$k$ coefficients are
bounded by $\|(c^2)_\alpha\|/c_{\min}^2+
\|(v^2)_\alpha\|/g_{\min}^2$ at $k$-order zero and $\alpha\ne0$,
$\|(c^2)_\alpha\|/(c_{\min}g_{\min})$ at order one,
$\|(c^2)_\alpha\|/g_{\min}^2$ at order two, and zero above two.
Formal compositions $(1-X)^{-1/2}$ and $(1-X)^{-1}$ bound finite
derivatives without a convergence-disk claim. The normalized Hamiltonian
majorant $\mathcal A$ has constant one, and

$$\|P_0\|\prec\tfrac12\mathcal A(1-X)^{-1/2}+\tfrac12,
\qquad
\|P_n^{\rm off}\|\prec\frac{\mathcal A F_n(1-X)^{-1}}{2g_{\min}},
\qquad
\|P_n^{\rm diag}\|\prec2P_0G_nP_0.$$

Here $G_n$ and $F_n$ are the absolute idempotence and transport recurrence
majorants, retaining mixed derivative factorials and star factors
$\binom lh/(2^ll!)$. Physical order9 and ordinary-$k$ order5, trimmed
to $(9-j,5-j)$ for $P_j$, cover $\partial_z^4P_4$. Multiply Taylor
coefficients by $t!z!k!$ for derivative bounds, and use explicit
$\sqrt2$ or $2$ only when converting to Frobenius or nuclear norms.
Zero rest gap is rejected. These all-$k$ bounds concern $P_0,\ldots,P_4$,
not physical $C_\Sigma$; initial covariance, true Weyl-defect and UV
remainder bounds remain separate.
