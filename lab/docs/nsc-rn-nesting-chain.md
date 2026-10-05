# Walking the nested-space question link by link

The question, in its own terms: there is a local nestation; the parent is
obtained by stepping out of it; the same math is then applied to the parent;
and the question is whether the parent produces a path into the child that is
a black hole.

This walks that chain and marks each link as measured, unanswerable here, or
refuted. It is assembled from measurements already made in this line. Nothing
new is computed for it.

## Link 1 — stepping out of the local nest: what is the parent?

**Measured.** The chart names three radial regions for `mu = 1.04`:

| region | interval |
|---|---|
| unnamed interior | `[1.0400, 1.3614)` |
| child | `[1.3614, 1.6113]` |
| parent | `(1.6113, 32]` |

Stepping out of the local nest is moving outward along one continuous chart
from `1.04` to `32`. There is no separate space to jump to.

## Link 2 — applying the same math at the parent

**Measured.** The characteristic speeds are smooth functions of radius across
the whole chart, with the lapse identically `1`:

| `r` | outgoing | ingoing |
|---|---|---|
| `0.80` | `−0.0185` | `−2.019` |
| `1.361` | `+0.0058` | `−1.994` |
| `1.611` | `+0.0482` | `−1.952` |
| `2.00` | `+0.1112` | `−1.889` |
| `7.90` | `+0.5028` | `−1.497` |

**This link is where the picture breaks.** Across the child/parent cut at
`1.6113` the step is `3.1e-04`. The child range `[0.006, 0.048]` runs
directly into the parent range `[0.049, 0.506]`. There is no discontinuity,
no kink, no change of character.

The child interval is **a label on a smooth metric**, not a region with its
own propagation physics. Applying the local math to the parent gives the same
smooth functions continued inward. The parent does not "produce" anything at
the cut, because there is nothing there for it to produce.

## Link 3 — does the parent have a path inward?

**Measured: yes.** The ingoing speed is negative at every radius sampled in
every run, from `0.76` to `32`. A path inward exists everywhere. This is not
the interesting part.

## Link 4 — is that path a black hole?

**Measured: yes, and this is the one real interface.** The outgoing speed
crosses zero at

> `r = 1.3257`, which is `r_plus`.

It is the only zero of the outgoing speed in every station of every run. The
region below it has **negative outgoing speed**: light cannot move outward
there at all. The causal ordering is one-way in all nine runs checked —
sealed inside, the zero count never changes, and the zero never reverses at
any resolved scale.

So the answer to the question as posed is **yes**: the parent has a path
inward, and that path crosses a genuine one-way membrane at `1.33` into a
sealed region.

## Link 5 — what is in there?

**Measured, and this is where the genuinely new content is.** The 2005 study
gives a criterion on the metric. It does not say what falls through. This line
says what falls through:

- a coherent shell falls at `−1.49` per unit `t/r_m`;
- it compresses **×1.58** against a deep boundary, saturating there;
- it then **relaxes to `4%`** of that peak, nearly completely;
- it **retains `94%`** of its accreted mass below the horizon;
- the layer **narrows by `2.27×`** as the boundary deepens, amplitude up and
  width down together;
- the horizon sits **above** every boundary this programme has placed, and
  every one of those boundaries is already in the trapped region.

The last point is the one that reframes the earlier work. The `×1.58` is not
matter piling against a wall that was chosen. It is matter compressing in a
region it cannot leave.

## Link 6 — is the sealed region "the child"?

**Unanswerable from the states, and it is a question about labels.**

The chart's labelled `child` is `[1.3614, 1.6113]`. The sealed region is
`r < 1.3257`. **The labelled child lies outside the sealed region.** The two
are disjoint, and the horizon sits in the unnamed interior below both.

Either the sealed region is the child and the chart's labels do not mean what
they appear to, or the child is not the sealed region and the black hole lives
somewhere else in it. The data cannot say which, and the answer decides what
the whole construction is about.

## Link 7 — nest to nest to nest?

**Refuted as posed for this construction.** It defines one child interval and
one parent. There is no second nest inside the first for time to flow
*between*. One-way flow is demonstrated for a single step, from parent into
the sealed interior, and for a single nest. Iteration is not represented at
all.

## What survives

The chain holds from link 3 onward and fails only at link 2: **the parent does
not do anything special at the cut where the child begins.** The one real
interface is the horizon at `1.3257`, and it is not a labelled boundary.

That is not the nested picture as described. It is a smooth charged geometry
with one causal membrane in it, and a region below that membrane whose
contents are now measured.

## Scope

This document is an assembly of measurements already in the records. It is
**not** a scientific result and it does not measure the actual universe. It
marks what is measured, what is not, and where the picture fails.

## Records

- trapping horizon and causal ordering: `nsc-rn-trapping-horizon-v1.json`,
  `nsc-rn-causal-order-v4.json`
- shell motion and compression: `nsc-rn-mass-shell-v6.json`,
  `nsc-rn-packet-sweep-settle-v3.json`
- layer shape: `nsc-rn-layer-shape-v4.json`
- balance against a known criterion: `nsc-rn-balance-at-horizon-v1.json`