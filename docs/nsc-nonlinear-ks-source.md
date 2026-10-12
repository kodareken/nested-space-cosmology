# Nonlinear Gaussian source in the common KS representation

The new owner restricts the already evolved PG solutions to actual KS
events and integrates their weak metric source. The comparison is against
the existing nonlinear endpoint-jet Gaussian derivative. It uses the same
metric, initial horizon/infinity preparation, source frequencies and spatial
Dirac discretization.

## Actual metric, actual field derivative

In the PG characteristic frame the canonical KS restriction is particularly
simple. With $v_\pm=\pm N_{\mathrm{PG}}/q_{\mathrm{PG}}-\beta_{\mathrm{PG}}<0$,

$$
\Psi_{\mathrm{KS}}=M\chi_{\mathrm{PG}},\qquad
M=\mathrm{diag}(\sqrt{-v_+},\sqrt{-v_-}).
$$

This is the existing normal/half-density frame on the **actual** geometry.
It is applied to resolved solutions, not called a unitary Cauchy propagator.
The evolved columns already include their source phases; no extra
$e^{iES}$ is multiplied onto them.

Since equal KS time means fixed $\rho$,

$$
\partial_z\Psi_{\mathrm{KS}}
=(\partial_\tau M)\chi_{\mathrm{PG}}
+M\partial_\tau\chi_{\mathrm{PG}}.
$$

The second term uses the actual PG equation, including its fixed exterior
inflow. Incoming $E$ is a source label after evolution, not the local
momentum. The calculation does not replace this derivative by $-iE\Psi$.

The symmetric kinetic insertion for each matrix $\sigma$ is

$$
K_\sigma=-\frac{i}{2}
\left(\Psi^\dagger\sigma\partial_z\Psi
-(\partial_z\Psi)^\dagger\sigma\Psi\right).
$$

Together with the mass/angular bilinears it forms the four existing raw KS
metric vertices. Their integration retains the fixed-chart Jacobian
$dT\,dz=d\rho\,d\tau/a_0(\rho)$. The instantaneous perturbed axial scale
does not replace $a_0$ in that coordinate Jacobian.

## What is recorded and compared

The spatial bulk is propagated before any output is projected. Within-step
quadrature samples use the same midpoint-frozen exponential as the prior
propagator. Continuous source coefficients use the actual supplied metric
at each quadrature event. Thus its finite-step comparison with the old
endpoint derivative is a measured numerical check, not an assumed identity.

The source is evaluated independently in PG and KS forms and compared with
$i\mathcal A_B$, where $\mathcal A_B=U_g^\dagger\delta_B U_g$ is stored by
the previous nonlinear jet owner. Contracting $-i\mathcal A_B$ through the
unchanged source covariance retains the existing Gaussian sign convention.

The fine output also records total fields and their KS spatial derivatives
on $\rho=0$, at the actual temporal quadrature nodes. Those data allow
equal-KS-time, unequal-PG-time kernel queries. They are observations of the
evolved full spatial field, not a few-mode replacement for that field.
Opposite-angular families use the already verified signed source/operator
intertwiner and their own covariance, rather than a second propagation.

The [record](../results/development/nsc-nonlinear-ks-source.json) distinguishes
coordinate-trace, endpoint-source, spatial and temporal errors. The physical
source-frequency sample is the inherited derivative control. Full spectral
completion and the common regulated reference trace remain required before
an absolute stress or stationarity claim. The selected canonical prescription
from [the common-trace scope](nsc-common-ks-trace.md) remains in force; no
raw-heat equivalence gate or additional $\Gamma_{\mathrm{rest}}$ is introduced.

```sh
python3 scripts/derive_nsc_nonlinear_ks_source.py --prepare
python3 scripts/derive_nsc_nonlinear_ks_source.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_nonlinear_ks_source.py
```
