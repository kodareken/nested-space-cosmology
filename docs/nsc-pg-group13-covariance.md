# Full-energy PG covariance of the first massive compact group

The retained group with $m=\pi/2$ and zero angular eigenvalue now has a
computed common-PG covariance, including the complete energy integral.
Its source is the same affine horizon covariance and inherited incoming
occupation. Neither historical seed generation is used as a spatial state.

The calculation compresses the full field onto eight orthonormal observables:
the original parent and child packets with both spin components, an orthogonal
child-bulk packet, and an exterior-bulk packet with both spin components.
These probes do not replace the eliminated bulk by an eight-mode Hamiltonian.

## Source and covariance

The field state remains the global mode construction

$$
C_{\mathrm{PG}}=\mathcal F(C_H\oplus n_{\mathrm{in}})\mathcal F^\dagger.
$$

For the original probe embedding $J$ and $Q=I-JJ^\dagger$, its unresolved-bulk
mode field and cross-covariance density are

$$
\Phi_Q(E,\rho)=\Phi(E,\rho)-J(\rho)F_J(E),
\qquad
C_{JQ}(E;\rho)=F_J(E)C_{\mathrm{src}}(E)\Phi_Q(E,\rho)^\dagger.
$$

The [state-block owner](../src/recursive_horizons/nsc_pg_state_covariance.py)
retains this spatial complement and its pair correlations. The numerical
eight-probe matrix includes four explicit probes of that complement.

## Complete energy account

The [reproducer](../scripts/derive_nsc_pg_group13_covariance.py) uses the
following connected parts of the same state:

1. Physical massive mode integration on the resolved low-energy panels.
2. The [causal subgap contour](nsc-pg-threshold-projection.md) on $[1,\pi/2]$.
3. Derived high-energy boundary modes followed by integration of the actual
   Dirac equation through the probe supports on $[16,160]$.
4. Massive mode endpoint Fourier series integrated from $160$ to infinity
   using generalized exponential integrals.

Mass, angular terms, the $A'/2$ connection term and the slow massive phase
remain in the high-energy calculation. The geometric characteristic
coordinates are reused; the LLL state or its covariance tail is not assigned
to this massive field. The approximation orders concern physical mode
evaluation and do not change the existing fourth-order stress reference or
local induced allocation.

The current fixes each branch normalization. Source terms omitted from the
vacuum approximation above $E=16$ have an explicit exponentially small
operator-norm bound from the existing horizon/incoming occupations. The
Fourier/normalization representation has a separate computable remainder
bound. Mode-order and spatial-grid changes assess the mode approximation;
they are reported as numerical convergence, not as a theorem about every
nonperturbative remainder.

## CAR is calculated, not filled

Both the spectral Gram matrix and the centered covariance are integrated.
The primary numerical covariance is assembled as

$$
C_{\mathrm{raw}}=\frac12G_{\mathrm{raw}}+
\int \Phi\left(C_{\mathrm{src}}-\frac12P_{\mathrm{src}}\right)
\Phi^\dagger\,\frac{dE}{2\pi}.
$$

The calculation does not insert $I-G_{\mathrm{raw}}$ to normalize a finite
frequency table. The independently computed $G_{\mathrm{raw}}$ is compared
with the canonical identity. Lesser, greater and Keldysh blocks retain this
same numerical Gram:

$$
G^<=iC_{\mathrm{raw}},\qquad
G^>=-i(G_{\mathrm{raw}}-C_{\mathrm{raw}}),\qquad
G^K=-i(G_{\mathrm{raw}}-2C_{\mathrm{raw}}).
$$

The [record](../results/development/nsc-pg-group13-covariance.json) contains
every matrix, the eigenvalues, the probe–bulk correlation norm and the error
account. It is a numerical PASS for this retained massive group at the
declared $3\times10^{-9}$ accuracy. Full retained C1b remains OPEN until the
other 31 massive groups have their corresponding covariance results.

| Quantity | Computed value |
|---|---:|
| Independently integrated CAR residual | $3.27\times10^{-11}$ |
| Covariance eigenvalue range | $[0.03570181,\,0.96429819]$ |
| Original-probe to bulk-probe correlation norm | $0.31816523$ |
| Finite-energy quadrature change in covariance | $2.81\times10^{-15}$ |
| Subgap contour-height change | $3.61\times10^{-12}$ |
| Tail mode-order change | $5.33\times10^{-15}$ |
| Tail spatial-grid change | $1.26\times10^{-12}$ |
| Fourier/normalization tail error bound for covariance | $5.70\times10^{-10}$ |

The action/scale ledger, old seeds, earlier certificates and PDF v0.25.0 are
unchanged. No stress tensor, selected metric history, transmitting endpoint
jet or metric timestep is inferred from this covariance result.

```sh
# Cold preparation of the already defined subgap controls, when absent.
python3 scripts/derive_nsc_pg_threshold_projection.py --points 48 --height-fraction 3
python3 scripts/derive_nsc_pg_threshold_projection.py --points 48 --height-fraction 4
python3 scripts/derive_nsc_pg_threshold_projection.py --points 24 --height-fraction 3

# Prepare the group's canonical inputs with the physical mode owners.
python3 scripts/derive_nsc_pg_group13_covariance.py --prepare-inputs --workers 2

# Recombine and compare all fields of the immutable record.
python3 scripts/derive_nsc_pg_group13_covariance.py --check
```
