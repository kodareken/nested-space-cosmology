# The PG source variation in raw KS metric coordinates

The transmitting source-mode variation can now be expressed in the raw KS
metric coordinates used by the history action. The conversion is evaluated
throughout the same compact region, with the same physical mode columns.
Its spatial dependence matters: multiplying an integrated PG derivative by
the neck's single Jacobian gives a different result.

## Fixed reference chart

Use the already owned future-child coframe. Write
$a_0(\rho)=\sqrt{\beta_0(\rho)^2-1}$ in the strictly trapped region. The
coordinate differentials are

$$
dT=-\frac{d\rho}{a_0},\qquad
dz=d\tau+\frac{\beta_0}{a_0^2}d\rho.
$$

They are integrable functions of the reference radial coordinate. This chart
is held fixed while varying the metric. It does not set a physical duration
or move the neck. In these coordinates the general raw KS orbit metric is

$$
ds_2^2=N_K^2dT^2-a_K^2(dz+\beta_K dT)^2.
$$

Here $a_K=q_{\mathrm{ADM}}$ is distinct from the PG radial scale. The
background values $N_K=1,\beta_K=0$ do not remove their independent variations.
Define

$$
s_K=\frac{\beta_0}{a_0^2}-\frac{\beta_K}{a_0},\qquad
F=a_K^2s_K^2-\frac{N_K^2}{a_0^2}.
$$

On the connected chart with $F>0,s_K>0$, the same metric has PG components

$$
q_P=\sqrt F,\qquad
N_P=\frac{N_Ka_K}{a_0q_P},\qquad
\beta_P=\frac{a_K^2s_K}{F},\qquad r_P=r_K.
$$

This supplies the actual Jacobian
$J^A{}_B(\rho)=\partial(\log N_P,\beta_P,\log q_P,\log r_P)^A/
\partial(N_K,\beta_K,a_K,r_K)^B$.
At the neck it is nonsingular with determinant one. The implementation
verifies the transformed two-metric, inverse field map and independently
differenced Jacobian.

## Pullback of the transmitting field variation

The existing compact shape is applied to the raw KS coefficients. Its
converted PG direction contains the full $\rho$-dependent Jacobian:

$$
M^K_B(E_o,E_i)=\int d\rho\,s(\rho)
\sum_A J^A{}_B(\rho)
\langle\Phi_{E_o},\delta H_A\Phi_{E_i}\rangle_\rho.
$$

The weak evaluation includes this Jacobian inside the integral. The
independent direct calculation transforms the perturbed KS metric and its
radial derivatives into PG variables, then differentiates the existing
Hamiltonian. Thus derivatives of the coordinate conversion are accounted
for, rather than only multiplying constant endpoint coefficients.

The previously computed mode columns are reused for every signed family;
there is no new horizon, scattering, covariance or stress integration.
The raw KS mode kernels also contract through the same finite source-fiber
CTP differential used by the PG calculation.

| Check | Maximum | Tolerance |
|---|---:|---:|
| Raw KS weak form versus transformed Hamiltonian variation | $2.18\times10^{-9}$ | $3\times10^{-8}$ |
| Hermitian energy exchange | $1.84\times10^{-15}$ | $3\times10^{-11}$ |
| CTP tangent condition | $2.09\times10^{-17}$ | $3\times10^{-11}$ |
| CTP differential versus source trace | $3.32\times10^{-18}$ | $3\times10^{-11}$ |
| Coordinate metric reconstruction | $1.78\times10^{-15}$ | $3\times10^{-11}$ |

These are derivative checks on the inherited quadrature fields, rather than
a new bound on the complete spectral integral. The constant-neck-Jacobian
control differs from the spatially pulled-back vertex by Frobenius norm
$0.476$–$7.256$ across the retained families; that shortcut is not used.

## Fixed normal and remaining endpoint map

At the fixed surface, $n_K=(1/N_K,-\beta_K/N_K)$ transforms by the inverse
coordinate Jacobian. Its four raw-field derivatives are explicit and agree
with the existing transmitting-domain unit normal. This is a metric normal
calculation, not a pointwise boost promoted to a Cauchy-state isometry.

**Raw KS metric-profile pullback PASS; full B1 OPEN.** The [record](../results/development/nsc-pg-ks-metric-pullback.json)
contains the coordinate residuals, fixed normal and all 63 mode-family
comparisons. An endpoint in the history action additionally needs the
specified spacetime support of its variation and its Cauchy embedding:
$\delta g_P^A(\tau,\rho)/\delta g^B_{K,e}$. The current compact direction
does not set those two endpoint profiles or select their history. Consequently
no complete `EndpointBranchJets`, two-sided Weyl mismatch or stationary
solution is inferred here. Physical scales and previous certificates remain
unchanged, and metric stepping remains closed.

```sh
python3 scripts/derive_nsc_pg_ks_metric_pullback.py --check
python3 -m pytest -q tests/test_nsc_pg_ks_metric_pullback.py
```

## Spacetime, nonlinear source and endpoint chain rules

A raw KS spacetime direction $\xi^B(T,z)$ pulls back as
$\delta x_P^A=J^A{}_B\xi^B(T(\rho),\tau+S(\rho))$, with

$$
\partial_\rho\delta x_P=J'\xi
+J(-\partial_T\xi/a_0+\beta_0\partial_z\xi/a_0^2).
$$

For $\xi=s(\rho)e^{-i\omega z}$ the spatial factor
$e^{-i\omega S(\rho)}$ and derivative $s'-i\omega S's$ are required.
The canonical first variation is
$\delta H=-i(\delta v\partial_\rho+\frac12\partial_\rho\delta v)+\delta M$;
its half-density transformation already includes the measure contribution.
No separate measure force is added. The PG source-mode kernel
$M_A(E_o,E_i)=\int\Phi_{E_o}^\dagger\delta H_A\Phi_{E_i},d\rho$
gives insertion
$-i e^{-iE_o(t_f-\tau)}M_Ae^{-iE_i(\tau-t_i)}$ with spectral measure
$dE_i/(2\pi)$. Sampled CTP weights are $\sqrt{w_ow_i}/(2\pi)$; source
covariance retains the coherent open fibers. The 33-group/63-family controls
do not complete the unsampled infinite trace or select a duration.

For a supplied nonlinear history with both PG speeds negative, restriction
to actual KS events uses
$\Psi_{\rm KS}=M\chi_{\rm PG}$,
$M=\operatorname{diag}(\sqrt{-v_+},\sqrt{-v_-})$ and
$\partial_z\Psi=(\partial_\tau M)\chi+M\partial_\tau\chi$.
Already evolved columns carry their phases; neither another $e^{iES}$ nor
the stationary replacement $\partial_z\Psi=-iE\Psi$ is admissible.
The symmetric insertion is
$K_\sigma=-\frac i2[\Psi^\dagger\sigma\partial_z\Psi
-(\partial_z\Psi)^\dagger\sigma\Psi]$,
with fixed coordinate measure $d\rho,d\tau/a_0$. This restriction is a
normal/half-density frame map on resolved solutions, not a Cauchy isometry.
Finite source checks still need spectral completion and reference allocation
before an absolute stress or stationarity claim.

The endpoint interface requires
$J^\pm_{A,e}=\partial V_\pm/\partial g^A_{\Delta,e}$ of shape
`(4, 2, d_c, d_c)` in the original eight-coefficient basis. Its algebraic
chain-rule tests do not supply physical branch jets. Full matching is
$R_{\partial,A,e}=w_{A,e}+b_{A,e}$, where $b$ is the same-action boundary
remainder's geometric pullback. Renaming `q` to `q_ADM` supplies none of
those data. The physical jets, remainder and two-sided mismatch remain
`null`; the diagnostic `93.54264532195464` is retained without a fitted
counterforce. The Gaussian action is counted once, and metric stepping
remains outside this retained calculation.
