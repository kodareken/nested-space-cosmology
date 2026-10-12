# The local action on the same spherical spacetime history

The local contribution now uses both PG coordinates, matching the domain of
the transmitting field/history jets. The old homogeneous KS implementation
is retained. Its identically zero shift variation is not imposed on the
spatially varying history.

The [new owner](../src/recursive_horizons/nsc_spherical_local_history.py) applies
the existing [finite-term curvature convention](../src/recursive_horizons/nsc_finite_terms.py)
and the [locked local ledger](../src/recursive_horizons/nsc_general_ks_local_history.py):

$$
S_{\mathrm{local}}=-4\pi\int d\tau\,d\rho\,\sqrt{-\det g_2}\,r^2
\left[A R+C_F\frac{q_{\mathrm{mag}}^2}{2r^4}
+C_W C^2+C_E E_4+C_{\Box}\Box R\right]+S_{\mathrm{GHY}}.
$$

The magnetic contraction uses the existing $F_{\theta\phi}
=q_{\mathrm{mag}}\sin\theta/2$ convention. There is no new multiplicity,
independent $R^2$ coefficient, cosmological term or interface force. All five
coefficients retain their signs and published values.

The supplied metric is the same compact raw-KS amplitude family already used
for the nonlinear field jets. It is pulled into the fixed PG chart using
$g_P=J^{\mathsf T}g_KJ$, with $T=T(\rho)$ and $z=\tau+S(\rho)$.
Second-order Taylor products give its exact coordinate derivatives, including
mixed time/space derivatives. These are derivatives of a specified metric
family, not new Einstein equations. Complex amplitude increments differentiate
the local action; they are not physical complex geometries.

The standard spherical curvature contractions are applied before integration.
Their static and homogeneous limits agree with the existing owners. Each of
the four amplitude derivatives acts on the full spacetime metric. The record
stores the separate Einstein, magnetic and Weyl action derivatives, their
sum, and its negative in the existing action-force convention.

## Boundary accounting

The metric perturbation and all its coordinate derivatives vanish on the
outer boundary of the integration domain. Consequently the existing GHY,
Euler and $\Box R$ endpoint variations are zero. The Euler density is also
integrated directly as a numerical cancellation check; its discretization
residual is not inserted as a source. $\Box R$ uses its already-owned boundary
identity instead of numerically differentiating flat endpoint tails four
times. The two orientations of the internal smooth $\rho=0$ cut still cancel.

Reported action differences subtract the fixed undeformed geometry only for
numerical presentation. Their derivatives equal those of the original action.
This fixed baseline does not replace the metric-dependent quantum reference
subtraction and does not set any physical vacuum energy.

## Composition status

The [result](../results/development/nsc-spherical-local-history.json) compares
time and radial quadrature separately, independently perturbed action values,
coordinate jets against the prior metric-value owner, and the inherited
curvature limits. This completes the local directional variation on the
specified family; it supplies no stationary history or absolute stress.

The full quantum spectral contribution and the fourth-order reference on the
spatially varying geometry remain to be joined. A homogeneous Bloch projector
cannot be copied into the initial PG source basis. The author's
[declared-action scope](nsc-declared-action-scope.md) remains in force:
no additional $\Gamma_{\mathrm{rest}}$ term is requested or invented, and
$93.54264532195464$ remains a discrete time-node diagnostic.

```sh
python3 scripts/derive_nsc_spherical_local_history.py
python3 scripts/derive_nsc_spherical_local_history.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_spherical_local_history.py
```
