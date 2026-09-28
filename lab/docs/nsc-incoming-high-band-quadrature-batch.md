# Direct quadrature of every retained high-energy source region

The numerical source used here is the complete **order24-minus-ad4 vacuum**
integral on each authenticated high-energy region. Its directed quadrature
replaces the original high-region source; it does not append an order24-minus16
correction to floating baseline quadrature.

## Decision and fixed execution budget

Reuse the [certified quadrature primitive](nsc-incoming-source-quadrature-bound.md),
its interval Riccati coefficients, exact Gauss-root brackets and analytic
disk bounds. The missing connection is the rigorous integral of the entire
high-region approximation. Group6 `[32,160]` is the first gate: if its
lapse-projected numerical enclosure meets `3e-11`, evaluate the remaining31
channels once at the same48-node,80-digit resolution. Use32 circle arcs,
radius-eight disks and eight-unit real cells, with at most three CPU workers.
Stop after this fixed batch and report PASS or OPEN against the unchanged
aggregate `3e-11` tolerance. Do not refine failed channels automatically.

The primary risks are carrying an unbounded old baseline through a certified
correction, overlapping the early-middle intervals, and accumulating parameter
enclosure widths across channels. The direct integral removes the first
issue, original metadata determine the exact disjoint unions, and the final
aggregate is summed with directed arithmetic before testing the tolerance.

## Authenticated region and source convention

Groups13 and14 own `[16,160]`, already a single middle interval. Every other
non-LLL group owns the union of the selected LOW `[32,40]` cell and its middle
interval `[40,Emax]`; `Emax=320` for groups10,11,12,31,32 and160 otherwise.
There are32 disjoint channel regions. Actual angular signs, their shared
multiplicity, archived source arrays and original producer metadata are
authenticated. No overlapping LOW interval is added to groups13 or14.

The frozen primitive's historical name `low_vacuum` denotes its **direct
order24-minus-ad4 source formula**, not a restriction to LOW panel labels.
This owner uses that same formula over each complete declared high region.
The previous pilot's middle correction and all earlier source receipts stay
unchanged. Point coefficients, fixed geometry and the exact interval Gauss
rule are cached within each worker; no radial or physical-mode solve occurs.

The original high source includes its archived thermal insertion. The new
direct value supplies the vacuum approximation; physical thermal remains a
nonzero separately bounded remainder. It is never declared zero. The record
retains original high-source values and the explicit replacement delta for
root composition. Those original values precede later source corrections:
they must be replaced once on the original finite source, with overlapping
later corrections accounted for, rather than added on top of them.

## Receipts and numerical proof

Content-addressed per-group artifacts store every interval kernel sample,
exact-root bracket, weight enclosure, analytic circle bound, provenance and
checked source integral. Resumable pointers authenticate each artifact before
reuse and never relaunch a completed channel. Replay contracts the saved
intervals and verifies the cell coverage without evaluating the source.
The original quadrature proof supplies the analytic and arithmetic guarantees;
this batch adds the full-region source composition and directed32-group sum.

The [batch record](../results/development/nsc-incoming-high-band-quadrature-batch.json)
keeps its direct source totals, aggregate quadrature/arithmetic enclosure and
explicit replacement delta separate from physical projector/thermal bounds.
Low-energy regions below these unions and the infinite tails retain their
own source and error owners. No physical IV, global surface solution,
constraint root or metric timestep is selected. No scale, seed or action
coefficient is altered.

```sh
python3 scripts/derive_nsc_incoming_high_band_quadrature_batch.py --pilot
python3 scripts/derive_nsc_incoming_high_band_quadrature_batch.py --prepare --workers 3
python3 scripts/derive_nsc_incoming_high_band_quadrature_batch.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_high_band_quadrature_batch.py
```
