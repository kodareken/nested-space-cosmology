# Recovering the retained spectral fields for the common source trace

At each real frequency the existing stationary interior Dirac equation has
two independent solutions. The six interior packet overlaps stored by the
PG preparation therefore contain enough information to recover the full
stationary source columns when the observation map is numerically resolved.
This avoids repeating horizon and scattering calculations for those data.
It does not turn the packet observables into a closed bulk system.

Let $F_E(\rho,1)$ solve the already owned radial equation with the fundamental
matrix convention $F_E(1,1)=I_2$. The actual three-source field is
$\Phi_E(\rho)=F_E(\rho,1)X_E$. With the unchanged parent, child and child-bulk
packet profiles, define the six-by-two observation matrix

$$
\mathcal A_E=\int d\rho\,J_{\rm interior}^\dagger(\rho)F_E(\rho,1).
$$

The stored **frequency-resolved** complex columns satisfy

$$
\mathcal A_E X_E=\Phi_{J,E}^{\rm stored}[0:6,:],\qquad X_E=\Phi_E(1).
$$

The new [owner](../src/recursive_horizons/nsc_pg_mode_recovery.py) integrates
only this short interior fundamental matrix and its packet overlaps. A
resolved left inverse determines $X_E$, including the phases in the existing
source basis. Neither a covariance nor a seed is inverted. The same source
covariance, open-fiber projector, magnetic/angular labels and gauge convention
remain attached to the recovered columns.

## Acceptance and error propagation

Both singular values are recorded. No singular-value floor or regularized
solution is admitted. The declared sensitivity criterion is

$$
\frac{\epsilon_{\rm input}}{\sigma_{\min}(\mathcal A_E)}
\le\epsilon_{\rm field},\qquad
\epsilon_{\rm input}=3\times10^{-11},\quad
\epsilon_{\rm field}=3\times10^{-8}.
$$

This is propagation of a declared input accuracy, not a theorem making
upstream numerical errors rigorous. The verifier also checks the actual
projection residual, source-projector consistency and current normalization:

$$
X_E P_E X_E^\dagger=[-v(1)]^{-1}.
$$

The fundamental matrix preserves the same current between $\rho=1,0,-1$.
Independent controls project the authenticated old fields from $\rho=0$ and
recover their separately archived values at $\rho=1$. The frozen transmitting
domain, spin frame and source normalization are retained.

Only accepted fields are serialized as usable rows. Rejected rows contain
NaNs plus a boolean mask and an explicit reason; consumers must call
`require_all()` or account for every omitted panel. The upstream horizon and
infinity preparation is unchanged.

## Spectral scope

The [record](../results/development/nsc-pg-spectral-mode-recovery.json)
applies this inverse to the authenticated real low/middle panels for every
massive group, including the retained refinements and the separate group13
and group14 preparations. Base and refined panels remain separate; this
record does not add overlapping quadrature panels or silently change their
weights. Negative frequencies retain the already verified opposite-angular
signed map. The LLL preparation is unchanged and not assigned to these modes.

Integrated subgap and infinite-tail eight-by-eight covariance matrices are
not frequency-resolved observations and cannot enter this inverse. Their
packet error estimates likewise do not bound a differentiated local source.
Any unresolved high-frequency rows require the existing massive asymptotic
mode construction applied to the source observable, not guessed occupation
or an inverse of the integrated covariance.

The recovered fields provide the previously omitted stationary bulk input
for conditional evolution and the common equal-KS-time state/reference
kernel. Full spectral subtraction, reference-band boundary terms and
convergence still belong to that calculation. This record supplies no full
stress, new state, stationary metric or metric timestep. All locked scales,
earlier certificates, the PDF and the public repository are unchanged.

```sh
python3 scripts/derive_nsc_pg_spectral_mode_recovery.py --prepare
python3 scripts/derive_nsc_pg_spectral_mode_recovery.py --check
```
