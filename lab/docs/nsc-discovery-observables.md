# Stage-1 discovery observer

`observe(pair, state, time)` in
[nsc_discovery_observables.py](../src/recursive_horizons/nsc_discovery_observables.py)
is the compact-row API for a later campaign worker. This note freezes that
row. It does not evolve a state, and the production atlas file was not
written in this wave.

The call is

```python
from recursive_horizons.nsc_discovery_observables import observe

row = observe(pair, state, time)
```

`pair` is a conformal `NestedPair`. `state` is a `NestedState`. `time` is
the caller's coordinate time. The return is a JSON-ready scalar row of
version `nsc-discovery-observe-v1`. Optional keywords are `capture_profiles`,
`bundle`, and `nodal_rate`. A bundle is accepted only together with the
coarse Cauchy rate of that same frozen-`W` image.

## What one row contains

The source probability is the occupation-weighted nodal mass. Its
coordinate density is that mass divided by the node spacing, and the
density per proper radial length divides again by `r Q`. Each fixed window
reports the probability, that probability divided by the window's proper
length, and the nodal peak of the proper-length density.

The windows are the coordinate intervals child `(1, 3)` and parent
`(0, 4)`. Normal energy, pressure work, lapse-gradient work, coordinate
metric work `F_Q Q̇ + F_L L̇`, and the three flux projections are all
nodal summands. Each window integral uses the conformal frame weights:
the trigonometric interpolant integrated and then divided by `dx`. On a
full period that equals the raw nodal sum. Boundary flux keeps the
collocation formula `(sample[a] - sample[b]) / dx`. The derivative matrix
still has a zero Nyquist symbol; the gap between that matrix and the
collocation flux is reported and is not removed by changing the matrix.

Proper lengths, proper means of `r` and `Q`, and the instantaneous clock
rates `r Q` at the pair's clock locations are copied from the existing
nested metric. Those rates are not integrated in the row.

Localization uses the existing proper-arc packet geometry. The weight is
the probability mass; `q` is not multiplied into it. The proper width is
the circular standard deviation of that packet times the carrier's proper
length. It is not a time.

Proper gradients are `(1/q) ∂_x` of the proper-length probability density
and of the areal radius, using the owned periodic derivative. Their values
at the four fixed cuts are part of the row.

## Windows are not mode projectors

The child and parent windows are coordinate cuts. The mode block is the
fixed original observer: geometry-child and geometry-parent coefficient
norms, plus the overlap of the live columns with the original child
columns and with the complementary columns. Those numbers answer a
different question from the probability inside `(1, 3)` or `(0, 4)`.

## Characteristics, clocks, and packet width

On the conformal null rays of `Q²(dt² − dx²)`, the coordinate speed is
`dx/dt = ±1`. The coordinate times required to cross the period-8 carrier,
the parent width 4, and the child width 2 are therefore 8, 4, and 2. The
row also carries the eight candidate rays that leave the fixed cuts
`0, 1, 3, 4` at those speeds. A ray that lands on a different fixed cut is
recorded as an ambiguity. None of these rays is a measured boundary.

Proper-clock integration is a separate accumulator. This row does not form
it, and it does not convert the packet width into a dispersal time.

## Moving windows

For the normal shell, `E = ∫ η_N dx` with `η_N` the summand divided by
`dx`. On a moving artificial cut

```
d/dt ∫_a^b η_N dx = ∫_a^b S dx + F(a) − F(b) + η_N(b) ḃ − η_N(a) ȧ.
```

The last two terms are the Reynolds term. `S` already contains the
volumetric pressure and lapse work. An artificial gradient cut does not
add `−p V`. `surface_action_work` is the separate surface-traction
evaluation and is not folded into that balance.

`control_mode="coupled"` uses the supplied coarse rate and bundle.
`control_mode="frozen_geometry"` keeps the current fields and sets the
geometry jets, including `Q̇` and `Q̈`, to zero, so the metric work and
the expansion pressure work vanish while the spatial lapse work remains.
The default is `coupled`. Naming `preparation` selects the saved-control
label; a live call does not assume `new_pair`.

The finite-difference slope of the normal shell is not part of this row.
`time_derivative_measured` stays false.

## Live curvature

The nested state is converted by `reconstruct_state` through the frozen
read-only geometry map. In coupled mode, `actual_q_second_rate` builds the projected `Q` jet.
In frozen-geometry mode that jet is zero. `direct_rh_grid` builds `R_h`
with `L = Q` and `β = 0`.
`W = (R_h − 2)² / (3 r⁴)` is the metric Weyl scalar. `chi + 2` is compared
afterward as `chi_shell_gap_max = max |R_h − chi − 2|` and is not written
back into `R_h`. A live call does not apply the saved replay-basis control;
that comparison belongs to the sealed records below. A new pair is flagged
`missing`.

## When to call it

`observation_due(step, every=..., on_rk_stage=True)` is false. Campaign
code should call `observe` on accepted output steps, not on the four
Runge--Kutta stages. The coarse Cauchy rate is `bundle["nodal_rate"]`
from `rates(pair, state, return_bundle=True)`. The nested coefficient
rate is a different vector and is rejected. `capture_profiles=True`
attaches the nodal arrays listed in `PROFILE_KEYS`. The default row does
not.

`assert_compact_row` checks the frozen keys and the refusals above:
renewal is not asserted, the continuum is not certified, and `chi` is not
substituted.

## Atlas check

[derive_nsc_discovery_atlas.py](../scripts/derive_nsc_discovery_atlas.py)
reads four sealed preparations. `--check` also recomputes the nf256
`T=0.3` conformal frame against the raw sums. `--write` evaluates the
saved baseline frames — original conformal `nf256`/`nf512` at `dt=0.0005`,
separated pair v1 `nf256`, and separated pair v2 `nf512` — through the
frozen replay-basis `W`. Original Cauchy momenta use `encode_state`
(`π = dx Wᵀ p`). Nested frames are already canonical and are not encoded
again. No trajectory is created and no initial radius is solved. The
successor is one new JSON+NPZ and is refused if either file exists:

| Preparation | Record |
|---|---|
| original conformal v2 | `results/development/nsc-spherical-conformal-episode-v2.json` |
| separated pair v1 | `results/development/nsc-nested-parent-child-v1.json` |
| separated pair v2 | `results/development/nsc-nested-parent-child-confirmation-v2.json` |
| frozen basis | `results/development/nsc-nested-parent-child-replay-basis-v1.json` |

Their schemas and top-level fingerprints differ. The frozen geometry
control is the saved replay basis, and only for separated pair v1, separated
pair v2, and the basis record itself. The check compares published hashes.
It does not rebuild `W`. The original conformal episode and any new pair
are marked missing.

```sh
python scripts/lab.py scripts/derive_nsc_discovery_atlas.py --check
```

`--write` creates the successor JSON and NPZ of those saved-frame series
and final profiles, and refuses to replace an existing file or a sealed
input. The record binds the producing commit, any dirty bound files, the
producer/observer/backend/core/consumer hashes, and every read sealed
JSON and NPZ. `--check --successor <file>` authenticates that record and
its stored time, curvature, localisation, and energy rows without
rewriting. The production files were not written here.

## Not claimed

No renewal, no continuum certificate, and no dispersal time. An imposed
window trajectory is not evidence that a boundary was regenerated. The
root index and claim ledger were not edited here; they remain with the
integrating checkout.
