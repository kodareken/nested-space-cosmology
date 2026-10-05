# The outward region is algebraic, not exponential

One claim about this picture had never been tested here: that the outward
direction is far more sensitive than the inward one, exponentially so.

The outward region is everything above the outer horizon. What says whether
it bends gently or violently is how the outgoing characteristic speed
approaches its asymptote. An exponential approach shows a local exponent
that grows without bound as radius grows; a power law shows one that
settles.

| r | local exponent |
|---|---|
| 1.5 | `-0.2639` |
| 10 | `-0.4747` |
| 100 | `-0.4976` |
| 1000 | `-0.4998` |
| 10000 | `-0.5000` |

It settles, on **`-0.5` to four figures**, at every mass ratio tried.

That is the closed form doing what a closed form does: the shift goes as
`sqrt(2M/r - Q^2/r^2)`, which at large radius is `sqrt(2M/r)`, so the
outgoing speed approaches one as `r^(-1/2)`. Nothing exponential appears in
the outward direction of this geometry.

The same measure taken inward over the trapped region does not settle on a
single exponent, which is what a region containing a turning point in the
shift should look like.

## What this is and is not

It is a negative result about **this chart**. The chart has no expansion, no
cosmological horizon and no dynamics at all beyond what was marched into
it, so it cannot speak to whether an expanding universe is exponentially
sensitive in the outward direction. It does not speak to dark energy, to
curvature at cosmological scale, or to anything in the parent beyond this
one radial coordinate.

It does say that the thing which has been called the outward direction all
along — the region above `r_plus` — is algebraic in this construction, and
that the outward/ inward contrast which motivated pushing the edge deeper
was not, here, an exponential against an algebraic.

## What that leaves open

The contrast between outward and inward sensitivity is still worth having,
just not with these two labels. In this geometry the inward region is where
the interesting structure is: the trapped band, the turn in the shift, the
factor-of-two compression, the three-sided shell, the floor. The outward
region is close to empty by comparison.

That is consistent with everything measured so far, and it is not the same
claim as the one this note set out to test.

## Scope

A construction check on the reference chart. **It does not measure the
actual universe** and it tests no cosmological claim.

Three mass ratios, `20000` samples per entry, no station payloads, no
marching.

## Provenance

- reduction: `lab/scripts/derive_nsc_rn_outward_scaling.py`
- record: `lab/results/development/nsc-rn-outward-scaling-v2.json`; v1 is
  retained and carries the same numbers under a key that did not match
  the prose
- tests: `lab/tests/test_nsc_rn_outward_scaling.py`
