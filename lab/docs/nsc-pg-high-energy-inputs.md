# Completing the archived high-energy mode inputs

The [spectral inverse](nsc-pg-spectral-mode-recovery.md) explicitly rejected
3,188 rows whose packet observations did not resolve the declared field
accuracy. Their original middle-panel producer already defines the boundary
mode columns, so these rows can be recovered directly from that authenticated
recipe. No singular-value regularization or new horizon calculation is needed.

`FastVacuumPacketProjector.project_many` initializes the interior field at
$\rho=1$ as

$$
\Phi_E(1)=\left(0,\;f_{+,E}(1),\;f_{-,E}(1)\right),
$$

where the two current-normalized vectors use its order-16 massive Riccati
series. The [new owner](../src/recursive_horizons/nsc_pg_archived_high_energy_modes.py)
exposes exactly these columns, including the original source order and phase
convention. It authenticates the middle-panel producer's dependency hashes
and reads its pinned expansion order. No extra clock phase is inserted.

## Checks and inherited approximation

The saved six-by-two observation matrix is sufficient to check the new
columns against the actual archived packet data:

$$
\mathcal A_E\Phi_E(1)-\Phi_{J,E}^{\rm stored}[0:6,:].
$$

This is multiplication by previously computed overlaps, not another packet
integral. The result also compares against accepted inverse rows in the
overlapping panels, checks current normalization, and evaluates the existing
Riccati equation residual. The order-12/order-16 difference is a numerical
approximation diagnostic; it is not promoted to a rigorous error bound.

The acceptance gate applies the equation/current/order checks to the **newly
filled indices**, while projection agreement and overlap checks cover their
whole original panels. A first implementation accidentally applied the new
$3\times10^{-11}$ threshold also to unchanged lower-energy rows. That OPEN
audit is preserved: the whole-panel Riccati maximum was
$1.288\times10^{-10}$ and the order difference $3.590\times10^{-10}$, near the
old lower panel edges. These remain explicit inherited-context diagnostics,
not passing residuals. No tolerance is raised and no accepted inverse field
is replaced. Input coverage does not certify uniform stress accuracy over
the original approximation domains.

The construction retains the original middle-panel **vacuum-mode
approximation**. In particular, the zero first interior source column is the
archived producer's convention, not a newly imposed incoming state or
counterflow. The physical source covariance remains the same horizon pair
plus inherited incoming occupation. Thermal/reflection corrections and
bounds for the differentiated stress still require their own spectral
accounting. Packet agreement alone cannot supply those bounds.

## Composable completion

The [record](../results/development/nsc-pg-high-energy-inputs.json) and its
artifact retain only the newly available rows and their indices. The original
inverse artifact, acceptance mask and OPEN rows remain byte-identical. The
consumer combines accepted inverse rows with these explicitly tagged direct
producer rows; it must cover each requested panel index exactly once.

After this gate all 54,548 archived real-frequency rows across the 138
base/refinement panels have input fields in their original approximation
domains. These are not 54,548 new particle modes or a combined spectral
quadrature: base and refinement panels remain separate. The integrated
subgap and infinite-tail matrices still do not provide frequency-resolved
fields for a local source.

The next calculation uses these inputs in the common KS state/reference
kernel, with source-specific subgap/tail accounting and convergence. This
record supplies no complete renormalized stress, stationary geometry, new
null contraction or metric timestep. The locked scales, source covariance,
PDF and public repository remain unchanged.

```sh
python3 scripts/derive_nsc_pg_high_energy_inputs.py --prepare
python3 scripts/derive_nsc_pg_high_energy_inputs.py --check
```
