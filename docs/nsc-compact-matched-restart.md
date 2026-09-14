# Versioned compact restart with the matched interior distance

This checkpoint reproduces the complete retained compact frequency grids from
the original horizon/incoming state. It keeps two explicitly named generations
side by side:

| Start convention | Version | Role |
|---|---|---|
| `historical_radial_r_h` | `compact-seed-v1-historical` | Regression against the immutable original seed |
| `matched_delta_q` | `compact-seed-v2-matched-delta-q` | New matched seed control, selected explicitly |

The historical matrices are not rewritten or made to agree with a different
start convention. The matched generation is a new control surface for the
same canonical seed-frequency representation. It is not yet a spatial PG
Cauchy covariance.

## 1. The changed normalization

Reuse the [verified collar identification](nsc-horizon-paired-pg-map.md).
With $\delta_q=q_{h,\mathrm{geom}}-q_{\mathrm{geom}}$ and
$|\rho-\rho_h|=r_h^2\delta_q+O(\delta_q^2)$, the matched argument is

$$
z_H=\sqrt{\lambda^2+(m r_h)^2}
\sqrt{\frac{2\delta_q}{\kappa_h}}.
$$

The historical interior argument divides this by $r_h$. That factor is
retained only in the named historical regression. The unchanged exterior
scattering owner still uses a radial offset and its correct inverse-radius
factor. No mass, flux, temperature, scale or physical interaction is adjusted.

## 2. Full retained compact grids, from physical source data

The calculation uses the existing 20 compact groups and their complete
64-node grids: 1,280 positive-frequency X blocks. For nonzero angular
eigenvalues it also computes the 1,152 opposite-angular Y blocks. Zero angular
modes are not artificially doubled. The already owned signed-frequency map
and boundary occupation laws produce the corresponding negative-frequency
blocks from source mode columns, rather than from stored seed matrices.

The grids, weights, complex reflection coefficients, transmission and incoming
threshold are those of the declared compact realization. Reflection on each
needed angular sign is computed once with the unchanged exterior owner and
saved in an authenticated scattering artifact. Modal verification reads this
artifact; it does not repeat the scattering calculation.

For each label, a single propagator U is evaluated with the same inherited
6,000-step commutator-free SU(2) method. It is applied to both initial frames:

$$
F_\mathrm{hist}=U F_{H,\mathrm{hist}}M,
\qquad F_\mathrm{match}=U F_{H,\mathrm{match}}M,
$$

$$
C_\mathrm{hist}=F_\mathrm{hist}(C_H\oplus n_\mathrm{in})F_\mathrm{hist}^\dagger,
\qquad
C_\mathrm{match}=F_\mathrm{match}(C_H\oplus n_\mathrm{in})F_\mathrm{match}^\dagger.
$$

Each frequency/angular label has its own U. Sharing it between start
conventions isolates the initial normalization change; U is not shared
between different angular signs or masses. The integration interval follows
from the fixed geometry and declared numerical collar, not a chosen physical
transport duration. No stress or reference subtraction is recomputed.

## 3. Two different comparisons

The historical generation must reproduce the original 1,280 X blocks within
$3\times10^{-11}$. The matched generation is compared with the separately
versioned matched checkpoint and with the preceding independent modal
controls where they exist. Its difference from the historical X remains a
reported diagnostic. Agreement between two different start conventions is
not manufactured.

The artifact preserves both covariances, mode maps, signed partners, labels
and source metadata. `CompactSeedGeneration.load` requires an explicit
start convention and artifact digest. It refuses an implicit/default
convention and rejects promotion of a seed-frequency table into a spatial
PG covariance. This prevents a quiet change of reference in downstream work.

The focused checks cover historical regression, CAR, modal normalization,
T feedback, signed partners and cold reproduction of the matched generation.
The unchanged magnetic and local horizon-frame certificates are reused.

## 4. What H1 closes, and what it does not

H1 closes when the named historical path reproduces its original control and
the named matched path reproduces its new control with the declared residuals.
The original seed remains permanently historical for this normalization.
The matched generation is the explicit control to select for subsequent
matched-convention calculations.

H2 still requires the complete massive mode transform into the common PG
Hilbert space, including its global current/Plancherel normalization and
appropriate contact contribution. Finite seed-frequency quadrature and
fiber coisometry do not establish that identity. Until it exists, no massive
spatial $C_{\mathrm{PG}}$, full C1b, transmitting endpoint jets or physical
Weyl mismatch is assigned. The LLL result, old certificates and PDF v0.25.0
remain unchanged; no metric timestep begins.

See the [restart record](../results/development/nsc-compact-matched-restart.json)
for the two immutable payloads, version tags and per-group residuals.

```sh
# One-time extraction of physical scattering inputs; not part of routine checks.
python3 scripts/derive_nsc_compact_matched_restart.py --prepare-scattering --workers 4

# Routine cold modal reproduction reads the authenticated scattering payload.
python3 scripts/derive_nsc_compact_matched_restart.py --check
python3 -m pytest -q tests/test_nsc_compact_matched_restart.py
```
