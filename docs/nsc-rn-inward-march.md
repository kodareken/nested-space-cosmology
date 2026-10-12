# The inward march to the inner horizon

Owning note for the campaign directory
[nsc-rn-inward-march-v1](../results/development/nsc-rn-inward-march-v1/) and
its assessment record. The question is the one the
[gradient-bubble note](nsc-rn-gradient-bubble.md) left open: the recorded
family stops at the horizon midpoint (the classical inner mass), outside
the RN inner horizon in every case, so the Hayward regularity criterion
(g ≥ −1 on the compressive force as r → 0) was undecidable. This campaign
evaluates the region the old grid never reached.

## What changed

The areal grid is extended inward from the midpoint to one cell above
r_minus: same spacing, same outer boundary at 32 r_m, same packet, same
coupled/frozen two-arm march, and the same `nsc_rn_pg.stage_rates` route.
Just above r_minus the outgoing speed is still negative, so the chart
remains pure outflow at the new inner edge and the extension is admitted
(for the μ=1.04 seed, `v_out[0] = -0.015`). Below r_minus the outgoing
speed is positive and the excision would emit; that region belongs to a
new construction and is not opened here.

## Runner discipline

Each case runs one process and writes creation-only stations in the
family payload format (state arrays, radius, coupled mass, and the three
source pressures), so the cone-map, one-way and gradient-bubble readers
consume the new stations unchanged. Any exception writes a structured
`station-failure.json` with the blocker, stage, and clock; no station is
committed from a partial step. The runner pins its own sha256 and the
payload hashes in every record.

## What the march is for

The assessment successor
([`assess_nsc_rn_inward_march.py`](../scripts/assess_nsc_rn_inward_march.py))
reads the new stations with the decision helper
([`nsc_rn_inward.py`](../src/recursive_horizons/nsc_rn_inward.py)): the
compressive-force log-slope is measured on the inner ten percent of
cells, which now sit one cell above r_minus, and the Hayward bound is
decided on evaluated cells with no further extrapolation. The verdict
`meets_hayward_on_evaluated_domain` is a measured property of the
extended domain, not a claim about the actual universe.

## Status

All three ε=0.001 cases reached all four stations (t = 5, 10, 15, 20)
with no failure file. The sealed assessment record is
[nsc-rn-inward-march-assessment-v1.json](../results/development/nsc-rn-inward-march-assessment-v1.json).

Sealed verdict, with the inner edge one cell above r_minus: the
inner-region log-slope minimum over the inner ten percent of cells at
the last station is −171.6 (μ=1.002), −80.2 (μ=1.010) and −62.3
(μ=1.040); `meets_hayward_on_evaluated_domain` is false in all three
cases. The slope at the second cell oscillates in sign at the horizon
corner (positive for μ=1.002 and μ=1.010), which is the boundary-stencil
corner, not the verdict; the verdict rests on the region minimum, which
uses the centered interior stencil.

The compressive force does not flatten at the inner horizon on this
branch. The Hayward criterion that the gradient-bubble note could not
decide is now decided on evaluated cells: the neutral-Dirac-sourced
interior is not Hayward-regular in this family, and the programme's next
named question is which source branch (charged lowest-Landau Gauss,
torsion, or a temperature/sea-corrected state) brings the slope above
−1 — that question is a branch choice, not another extrapolation.

## Owners

[Runner](../scripts/derive_nsc_rn_inward_march.py),
[decision helper](../src/recursive_horizons/nsc_rn_inward.py),
[assessor](../scripts/assess_nsc_rn_inward_march.py),
[tests](../tests/test_nsc_rn_inward.py).