# Common Cauchy state and Gaussian trace in KS coordinates

## Execute the selected canonical prescription

The adopted action is the canonical Gaussian CTP functional, the declared
fourth-order reference subtraction, and the locked local contribution once:

$$
\Gamma_{\mathrm{selected}}^{\mathrm{CTP}}
=\Gamma_G[C_0]+\Gamma_{\mathrm{sub,ad4}}^{\mathrm{CTP}}
+S_{\mathrm{local}}[g_+]-S_{\mathrm{local}}[g_-].
$$

This follows from the [canonical owner](nsc-causal-common-functional.md),
the [compact subtraction prescription](nsc-compact-ctp-completion.md), and
the author's [action inventory](nsc-declared-action-scope.md). The canonical
and finite proper-time kernels are explicitly distinguished in that choice.
Exact equality with the complete raw heat-vacuum functional is a stronger
UV-equivalence claim, not a prerequisite for executing this selected action.

The preceding reference-symbol/band-action records remain unchanged. Their
formal computations still apply, but this note corrects the subsequent work
gate: do not invent a finite conversion merely to force raw-heat equality.
The actual execution requirement is a common trace, measure and boundary
accounting for the three already specified contributions. This does not set
any missing source, cutoff boundary term or finite coefficient to zero.

## Restrict resolved solutions, rather than copying a seed

The supplied perturbation vanishes for $\rho\ge1$. Both Dirac characteristics
in the relevant trapped collar point toward decreasing $\rho$. The global
incoming modes therefore remain the reference modes on $\rho=1$, for every
axial coordinate $z$. That spacelike slice is a Cauchy surface for the future
interior collar; no statement about a Cauchy surface of the complete exterior
spacetime is required.

With the already fixed chart $T=T(\rho)$, $z=\tau+S(\rho)$ and the owned
canonical spin frame $F_{\mathrm{frame}}$, restriction gives

$$
k=-E,\qquad
\Psi_{\mathrm{KS}}(E,\rho)
=e^{iES(\rho)}F_{\mathrm{frame}}(\rho)^{-1}\Phi_{\mathrm{PG}}(E,\rho),
$$

$$
C_{\mathrm{KS}}(k)
=\Psi_{\mathrm{KS}}\,C_{\mathrm{source}}(E)\,
\Psi_{\mathrm{KS}}^\dagger.
$$

The input is the archived globally resolved field, together with the same
horizon pair and inherited incoming covariance. The scalar phase follows from
$e^{-iE\tau}=e^{-iEz}e^{iES}$; it is not an assigned transport duration.
The spin transformation is used after restricting actual solutions, not
declared to be $U_{L0}$.

Current-normalized scattering gives the coisometry
$\Psi_{\mathrm{KS}}\Psi_{\mathrm{KS}}^\dagger=I_2$. Combined with
$dk/(2\pi)=dE/(2\pi)$, this supplies the continuum Fourier CAR normalization.
The source covariance retains its horizon coherence. Its restriction is a
positive Gaussian covariance; source-complement correlations are recorded
separately and are not claimed to reconstruct an exterior joint state.
The evaluated samples do not replace the continuum by a closed finite system.

## The same metric variation in both representations

On the common reference geometry, the KS weak vertices act at midpoint
momentum $\bar k=-(E_i+E_j)/2$. They vary all four raw metric fields. The
spacetime measure is

$$
dT\,dz=\frac{d\rho\,d\tau}{a_0(\rho)}.
$$

The [new owner](../src/recursive_horizons/nsc_common_ks_trace.py) reconstructs
the complete source-space matrices in these coordinates and compares them
with the authenticated PG vertices. It retains the relative spatial phase,
the normal/half-density factor and the Jacobian before the source trace.
This checks off-diagonal frequency/source entries as well as diagonal ones.

The [record](../results/development/nsc-common-ks-trace.json) evaluates all
63 signed families on their inherited derivative-control frequencies. It
also constructs the incoming KS covariance directly from the stored global
fields at $\rho=1$. No scattering, horizon or old stress generator is rerun.

The remaining work is to use this common Cauchy representation for the
nonlinear combined Gaussian/reference trace, with convergent spectral tails
and explicit numerical cutoff terms, then add the locked local action once.
The local symbol cutoff terms cannot be discarded merely because the
metric perturbation is compact in spacetime. Physical source/endpoint
equations and stationarity remain open; no metric evolution is started.

```sh
python3 scripts/derive_nsc_common_ks_trace.py --prepare
python3 scripts/derive_nsc_common_ks_trace.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_common_ks_trace.py
```
