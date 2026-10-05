# The parent/child cut carries nothing

The construction nests a parent and a child on one chart, with the cut at
`r_plus + gap/2`. The only physical interface the chart contains is at
`r_plus`, which lies **below** the child rather than at its boundary. So the
question is whether the boundary between the two spaces carries any
signature at all.

Every quantity the chart carries was tested for one: the lapse, the shift,
the redshift and both characteristic speeds, at three mass ratios, with the
cut located to one part in `40000`.

## The statistic, and one that was wrong first

The first attempt compared each field's value at the cut against the median
of a local window. It returned exactly zero for every field, which looked
like a clean result and was not.

Over a locally linear stretch, the median of an odd window of nodes **is**
the centre node. A smooth field scores zero — and so would a symmetric bump
sitting exactly on the cut, which is the thing most worth detecting. The
statistic was blind in precisely the case it was built for.

The second difference is the statistic that can see a kink, and it is what
the record uses.

## What it found

| mu | cut | fields with a feature |
|---|---|---|
| 1.002 | `1.1286` | **none** |
| 1.010 | `1.2935` | **none** |
| 1.040 | `1.6113` | **none** |

Curvature excess at the cut runs `0.29` to `0.83`. That is the local
curvature divided by the chart's own maximum, so every field is *less*
curved at the cut than somewhere else on the chart, which is what a smooth
field looks like and not what a boundary looks like.

**There is no signature at the parent/child cut.** No kink, no change of
character, no condition, no membrane.

## What that means

In this construction the parent/child nesting is **bookkeeping**. It names
two intervals on one chart and draws a line between them, and the geometry
does not know the line is there.

That is consistent with everything measured along the way. The child and the
parent sit on the same smooth metric, with the lapse pinned at one and the
two characteristic speeds summing to two at every radius. The only place the
geometry changes character is `r_plus`, and that is the membrane between the
child and the trapped region.

## What it does not mean

It does not show that a parent/child boundary is unreal. It shows that
**this chart carries nothing there**, and that if the nesting is meant to be
physical then the boundary condition it needs has not been supplied.

Which is the sharpest statement available about the open question: not "we
looked and found nothing" but "the construction asserts a nesting and
supplies no physics for it, and the one interface it does supply is in the
wrong place for that reading."

## Scope

A construction check on the reference chart. **It does not measure the
actual universe.** Finding no signature is a statement about this chart. A
different construction could put a real boundary at `r_plus + gap/2`; nothing
here shows that one cannot exist.

Three mass ratios, `40000` samples per entry, no station payloads, no
marching, no construction code changed.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_child_cut.py`
- records: `nsc-rn-child-cut-v2.json` is current; `v1` used the degenerate
  median statistic and is retained
