# The one-way enclosure on the RN radial parent

Owning note for [nsc-rn-oneway-v2.json](../results/development/nsc-rn-oneway-v2.json).
One finite question in one named domain: do there exist modes of the
recorded state's own Dirac operator that carry outward current from
inside the measured trap?

## Domain and method

Areal radius on the recorded radial grid of `nsc-rn-neutral-family-v1`.
The operator is the instantaneous rate matrix of `nsc_rn_pg.dirac_rates`
on the reconstructed coupled lapse and shift — the same route the march
uses, guarded by the mass-reconstruction check of the cone map. A mode is
interior when at least half its weight sits below the state's own measured
`trap_radius`. The sample is the final station of each case, where the
shell has sunk farthest in.

The direction statistic is the **signed** probability current of the
chart,

    S = sum( v_out |chi_+|^2 + v_in |chi_-|^2 )
      / sum( |v_out| |chi_+|^2 + |v_in| |chi_-|^2 ).

Inside the trap both characteristic speeds are negative, so channel
content there is inward weight; the signed form is the one that reads
"can the inside emit outward". v1 used speed magnitudes, which counts a
trapped mode's chi_plus content as outward weight and mislabels the
same physics; v1 is sealed and named in the v2 `supersedes` block.

Each dense spectrum runs in a worker under a 240-second wall guard; a
station that does not return records `eigen_timeout` with no partial
spectrum. Every case also carries an eig-free `state_witness`: the signed
current and interior weight of the recorded state itself.

## What the record says

| Case | interior max S (bands) | verdict | state witness S (W_in) |
|---|---|---|---|
| μ=1.002, ε=0.001 | — | eigen_timeout | −0.840 (0.048) |
| μ=1.002, ε=0.010 | — | eigen_timeout | −0.841 (0.084) |
| μ=1.010, ε=0.001 | +0.846 (band |Im λ|<0.5) | interior_cardinal_outward | +0.278 (0.029) |
| μ=1.010, ε=0.010 | +0.823 | interior_cardinal_outward | +0.301 (0.041) |
| μ=1.040, ε=0.001 | +0.111 (band 0.5–1.0) | interior_cardinal_outward | +0.423 (0.022) |
| μ=1.040, ε=0.010 | −0.257 in every band | strictly one-way | +0.433 (0.033) |

Only the deep, thick case (μ=1.04, ε=0.010) is strictly one-way in every
resolved frequency band. The interior-supported transients that do carry
outward current live in the low band: the slow modes that straddle the
measured trap and ride the exterior outgoing channel. The measured trap
sits outside the source-free horizon (see the cone map), so "interior"
here includes an annulus outside r_plus where outward travel is allowed,
and the escape capability measured is exactly the content of that
annulus plus the slow straddling modes.

The recorded state's own signed current is outward for μ ≥ 1.01 and
inward for μ = 1.002: the lightest-mass shell is already beyond the
point of return at its final station, while the heavier shells are still
on their way in.

## Horizon-reaching domain (v3)

The bounded reduction now runs over the inward-march campaign as well
(`--family`/`--record`, with the superseded record named on the command line).
Record: [nsc-rn-oneway-v3.json](../results/development/nsc-rn-oneway-v3.json),
sample = the t=20 final station of each ε=0.001 case.

| Case | interior max signed current | reading |
|---|---:|---|
| μ=1.002 | — | `eigen_timeout` (240 s worker guard) |
| μ=1.010 | −0.356 | **strictly one-way** |
| μ=1.040 | +0.283 | interior transient remains |

The intermediate case reverses relative to the family sample, where its
interior maximum was +0.823: on the domain that reaches the inner horizon,
μ=1.010 carries no interior-supported outward current in any resolved band. The
heaviest case keeps a positive transient (+0.283 against +0.111 on the family
sample), consistent with the trap-offset ordering — the case with the smallest
measured offset retains the most escape capability. The lightest case did not
close its dense spectrum and is recorded as such, beside its eig-free witness.

## What this does not say

Version two does **not** falsify the 2005 reading: it localizes it. The
Bronnikov–Fabris one-way statement is about the classical trapped region;
this record measures the instantaneous operator of the recorded sourced
chart, one trap offset outside r_plus, and finds the outward capability
confined to the low band and absent in the deep thick case. It is an
instantaneous spectral statement, not a retarded response across time,
and it does not measure the actual universe.

## Owners

Module [`nsc_rn_oneway.py`](../src/recursive_horizons/nsc_rn_oneway.py),
driver [`derive_nsc_rn_oneway.py`](../scripts/derive_nsc_rn_oneway.py),
tests [`test_nsc_rn_oneway.py`](../tests/test_nsc_rn_oneway.py). The trap
radius comes from [`nsc-rn-cone-map.md`](nsc-rn-cone-map.md). Sealed v1
remains preserved at
[`nsc-rn-oneway-v1.json`](../results/development/nsc-rn-oneway-v1.json).