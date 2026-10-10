> **Withdrawn from publication, 10 October 2026.** This historical draft is
> retained to document the interpretation error. Pressure equality does not
> establish mechanical balance, Hayward does not supply its claimed slope
> threshold, and reaching a positive inner horizon cannot determine centre
> regularity. Use the [corrected owner](../../lab/docs/nsc-rn-gradient-bubble.md).
> The text below is not an active publication proposal and must not be merged
> into a paper as a measured force-balance result.

# Pending section: two-force gradient balance on the RN radial parent

Status: draft for review. Not merged into `paper/nested-space-cosmology.md`,
because that source is sha-pinned in `paper/catalog.json` and the PDF builder
needs ReportLab and the fonts the pinned PDF was made with. Merging this text
requires a rebuild and a re-pin on the build machine.

Owner note: `lab/docs/nsc-rn-gradient-bubble.md`
Record: `lab/results/development/nsc-rn-gradient-bubble-v5.json`
Module: `lab/src/recursive_horizons/nsc_rn_gradient_bubble.py`
Figure: `lab/results/development/nsc-rn-gradient-bubble-v5.png`

---

## Purpose

The foundation paper carries the periodic spherical chart. The RN radial
parent is a second, nonperiodic chart, and its gradient balance is a measured
result that appears in no paper. This is the text that would carry it.

## What the section would say

The question is finite and stated in one domain: on the areal radius of the
nonperiodic horizon-regular parent, do the enclosing and the enclosed form a
balanced pair, or does one force run away and leave the other as decoration.

Two force densities are read from the same recorded station payload and
neither is fitted. The compressive side is the radial pressure with the
circulating term removed, `(force_N - angular_force) / (4 pi volume)`. The
resisting side is the transverse pressure, `-force_r / (2 lapse radial_density
r 4 pi)`. Both come from one angular bilinear in the same spinor with
`kappa = 1`. Balance is `B(r) = |F_in| / (|F_in| + |F_out|)`, scale free;
`B = 1/2` is equal forces. The gradient is `d ln|F| / d ln r`.

The outside acts in every recorded case. The integrated outer force is 15.6 to
82.0 percent of the inner force, and no case falls below a tenth. The balance
crosses one half at 49.7 percent of stations, so the two forces genuinely reach
equality.

The crossing does not hold a radius. Its persistence falls as the exterior
compression strengthens, from about 65 percent of stations at `mu = 1.002` to
about 31 percent at `mu = 1.040`, and the inner excursion is pushed from 1.04
to 1.29 in units of the mass scale. A stronger exterior dissolves the crossing.
The balance is therefore real and transient, and a sustained gradient bubble is
not established by this family.

The member that reaches the verdict `two_force_balance` is
`mu = 1.040, eps = 0.001`, at balance 0.915 and throat slope -7.37. It is
recorded as corroboration of a known result, not as a discovery.

## Two withdrawals this section must carry

An earlier version judged the outcome against a chosen throat-slope limit of
eight. The measured slopes straddled it at -7.37 and -8.44, so the verdict
pointed the wrong direction. That limit is withdrawn and the verdict is now
the measured balance structure, which needs no chosen number.

An earlier version also recorded a hard `interior_not_regular` dead end. That
over-claimed. The grid's inner edge is the current inner mass and lies outside
the RN inner horizon in every case: `r_minus` is 0.9387, 0.8682 and 0.7543 for
`mu = 1.002`, `1.010` and `1.040`, leaving an uncovered gap of 0.063 to
0.286 mass units, and the centre is never reached. What can be said is bounded:
the fitted throat power law extrapolated across that gap would grow the force
by about `10^0.9` to `10^4.3`. The flag is now
`throat_extrapolation_diverges` and it does not decide the verdict.

## The 2005 criterion

Douglas's framing pointed at Hayward, `gr-qc/0506126`, "Formation and
evaporation of non-singular black holes", whose criterion is `g >= -1` across
a trapping horizon. That criterion is sourced in this work. For the member
that reaches the verdict `two_force_balance`, the balance at the throat is
0.915 and the exterior carries 15.6 percent of the inner force.

That figure is single-station. Across the family the balance at the throat is
not one number: its per-case mean runs 0.339 to 0.549 and its per-station range
runs from 0.0002 to 0.988. Reporting a single value without the range would
overstate what the reading supports.

It is recorded as corroboration of a known result and not as a discovery,
because this chart's domain does not reach the horizon. The grid's inner edge
sits outside `r_minus` in every case, and the trapping character is read at
the boundary. That is a limit of the instrument, and it is why the balance at
the horizon cannot be promoted to a test the parent passes.

## Proposed claim marker

The text would carry a new inline marker beside the others in the paper, so
that `scripts/check_claim_bindings.py` covers it:

```
<!-- nsc-claim:two-force-gradient-balance -->
```

with the record `../lab/results/development/nsc-rn-gradient-bubble-v5.json`
and the note `../lab/docs/nsc-rn-gradient-bubble.md` linked from the
paragraph, matching how the existing markers carry their backing.

## A number that did not survive checking

An earlier draft of this section stated the balance at the trapping horizon as
0.5407. That value could not be sourced. It came from a working summary, and a
search of the repository finds it only inside an unrelated constraint-germ
artifact, where it is a different quantity. The recorded throat balances are
the per-case means and per-station ranges quoted above. The unsourced figure
is removed rather than adjusted.

## Scope

Draft text for review. It measures nothing new. Every number above is taken
from the recorded note and its record, and the two withdrawals are stated
because the record retains them.