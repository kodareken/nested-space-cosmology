# Common-PG covariance of the retained charged field

The massive channels now use their own globally normalized horizon/infinity
mode fields on the same PG slice as the transmitting response. The calculation
retains both angular signs, the horizon-partner correlations and the inherited
incoming occupation. The matched start convention is `matched_delta_q`;
historical seed matrices remain immutable comparison data.

**Retained-state C1b PASS:** all 32 massive groups and the separate LLL
preparation are composed, covering 63 evaluated signed families. This closes
the common-PG Cauchy-state input to C1 at the declared numerical accuracy.

For each retained angular/compact family, the state is

$$
C_{\mathrm{PG}}=\mathcal F(C_H\oplus n_{\mathrm{in}})\mathcal F^\dagger.
$$

The magnetic basis makes the covariance between different angular signs zero.
This does not discard the correlations inside the two-channel horizon block
$C_H$. Each signed angular family has its own massive radial mode map. Negative
frequency is supplied by the previously verified paired-spinor identity.

## One connected energy integral

The low-energy mode equations are evaluated in independent vectorized batches.
Vectorization introduces no interaction between frequencies. The original
Dirac connection, finite horizon matching, Jost phase and current normalization
are retained. The new quadrature resolves the angular-barrier transition; old
coarse integrations are controls, rather than sources of spatial covariance.

For compact mass $m$, the real panels omit $[1,m]$, which is supplied by the
already computed causal subgap contour. The middle window continues to the
selected numerical endpoint $L$, and the massive endpoint series integrates
the tail to infinity. $L=320$ is used for groups 10–12 and 31–32; otherwise
$L=160$. These are integration splits, with no change to physical scales.
Groups 13 and 14 retain their already evaluated lower transition at 16.

All positive-energy components use the same measure $dE/(2\pi)$. Write their
integrated Gram and centered covariance as $(G_s,K_s)$. With the locked signed
packet map $S$, both angular families enter the result:

$$
G_{\rm raw,s}=G_s+S G_{-s}^{*}S,
\qquad
C_{\rm raw,s}=\frac12G_{\rm raw,s}+K_s-S K_{-s}^{*}S.
$$

Angular-zero groups use their own opposite-frequency counterpart without
duplicating their physical angular multiplicity. The Gram is integrated
independently; no missing tail is replaced by the identity.

## CTP data and the bulk complement

Each massive group is evaluated on eight orthonormal test coordinates: the
four original parent/child coordinates, an orthogonal child packet with two
spin components, and an exterior packet with two spin components. The separate
LLL result retains its seven-coordinate test. The original four-coordinate
probe embedding agrees in every group.

The whole-field state remains a source-mode operator. For $Q=I-JJ^\dagger$,
the existing state-block owner evaluates

$$
\Phi_Q=\Phi-JF_J,
\qquad C_{JQ}(E;\rho)=F_J C_{\rm src}\Phi_Q^\dagger.
$$

The additional packets inspect this complement; they do not exhaust it or
turn the field into a closed eight-mode Hamiltonian. Equal-time contour data
retain the actually integrated Gram:

$$
G^<=iC_{\rm raw},\qquad
G^>=-i(G_{\rm raw}-C_{\rm raw}),\qquad
G^K=-i(G_{\rm raw}-2C_{\rm raw}).
$$

## Evidence and numerical scope

The [record](../results/development/nsc-pg-retained-covariance.json) contains
every matrix and a residual entry for every signed family. Its numerical
tolerance is $3\times10^{-9}$. Covariance-specific low-energy refinement is
recorded for groups 6, 22 and 32; middle-band refinement for 12, 22 and 32.
Group 13 imports its completed certificate; group 14 retains its independently
computed refinement in both angular signs. Subgap and tail inputs retain
their separate authenticated records.

| Check | Maximum or range | Tolerance |
|---|---:|---:|
| Independently integrated CAR residual | $1.40\times10^{-10}$ | $3\times10^{-9}$ |
| Covariance eigenvalues | $[0.0027381,\,0.9972619]$ | $[0,1]$ within numerical tolerance |
| Low-energy covariance refinement | $2.89\times10^{-11}$ | $3\times10^{-9}$ |
| Middle-band covariance refinement | $1.18\times10^{-14}$ | $3\times10^{-9}$ |
| Batched field against independent massive mode owner | $2.94\times10^{-10}$ | $3\times10^{-9}$ |
| Fourier/normalization covariance-tail bound | $9.73\times10^{-10}$ | $3\times10^{-9}$ |
| Physical/middle-band join, Gram density | $4.88\times10^{-12}$ | $3\times10^{-9}$ |

The original probes have nonzero correlations with the sampled bulk in every
family; the smallest recorded operator norm is $0.10258$.

The Fourier/normalization tail bound and numerical mode-approximation checks
are distinct. Refinement measurements do not constitute a uniform theorem
about all spectral remainders. Source occupation corrections omitted from
the middle/tail vacuum approximation decrease exponentially with energy
under the same declared horizon and incoming distributions.
The [join audit](nsc-pg-retained-join-audit.json) retains physical and
approximate packet columns at existing frequencies immediately below 40 for
the largest angular and largest mixed compact groups. The old angular-only
seed grid ended below 16, so its new physical low-energy input supplies the
appropriate join control. The omitted source correction above 40 has the
conservative log-norm bound
$\log 6-40\min(\pi/\kappa_h,2\pi/(\Omega\kappa_h))$.

The current gate status is stored explicitly in the record. A completed PG
preparation supplies Cauchy-state data to C1. Full transmitting endpoint jets
still require unitary history derivatives and the KS endpoint pullback.
The remaining boundary derivative, two-sided Weyl mismatch and extended
stationarity are subsequent calculations. This construction assigns no
finite stress, selects no history and starts no metric timestep.

```sh
# Recombine authenticated inputs; does not rerun old mode or stress generators.
python3 scripts/derive_nsc_pg_retained_covariance.py --check
python3 -m pytest -q tests/test_nsc_pg_retained_covariance.py

# New input preparation accepts explicit paths to computed components.
python3 scripts/derive_nsc_pg_retained_covariance.py --prepare COMPONENTS.json
```

The component manifest lists `name` and `path` for each `low/G_S`,
`mid/G_S`, `low_ref/G_S`, `mid_ref/G_S` and optional `low_coarse/G_S` artifact,
plus `group14_input` with its repository path and SHA-256. Low inputs can be
recomputed by `derive_nsc_pg_batched_low.py`; middle inputs by
`derive_nsc_pg_mid_inputs.py`. The assembly authenticates the existing subgap,
tail, LLL and group-13 records. Seed covariances are excluded from its input
contract.

The `mode_control` component compares newly vectorized physical columns with
the authenticated independent mode-resolution artifact at four existing
frequencies of the largest mixed group. Its field values and indices are
retained so the residual is recomputed during record verification.
