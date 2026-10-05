# Paper coverage audit

What the two papers own, what the laboratory has measured since, and where
the light result sits. Written because the papers are dated 15 September and
1 October while the ledger runs to 4 October, so a coverage question is a
real question with a checkable answer.

## The two documents

| | foundation | companion |
|---|---|---|
| title | Nested-Space Cosmology | Regional Transfer, Geometric Feedback, and Local Memory |
| version | 0.26.0 | 0.27.0 |
| date | 15 September 2026 | 1 October 2026 |
| source | `paper/nested-space-cosmology.md` | `paper/finite-regeneration/main.tex` |
| pages | 1612 source lines | 14 |
| status | WORKING PREPRINT | MEASURED_FINITE_REALIZATION |
| claims | 27 inline `nsc-claim` tags | — |

The foundation paper states its own evidence base: 100 preprint records and
51 post-preprint development records. The laboratory now holds 159 result
JSONs, 553 lab notes, 343 source modules, 424 tests and 489 scripts across
234 commits.

## The light result: yes, it is in the paper

This was the question, so this is the answer with the location.

**The claim lives in §2.1, "Constants and physical resolution".** "The causal
cone follows from the physical principal symbol. A dimensional speed is
expressed using the room's clock and ruler conventions."

That is the statement: the speed is a number in the room's own units, and the
room's units are the room's. Read locally, `c = 1` in every room. Read from
one room to the next, the same `c = 1` is a different dimensional speed
because the clock and ruler conventions differ.

**The mechanism lives in §2.2 and two records.** The horizon-penetrating
Hamiltonian retains the shift and the spin connection,

```
H_k = -i(sigma_2 - beta I) d/drho + (i/2) beta' I + (kappa/r) sigma_1
```

and `nsc-4-lorentzian-transport.json` gives the characteristic coordinate
speeds as `1 - beta` and `-1 - beta`. At the throat those are about `-1.1708`
and `-3.1708`: both travel toward the child side even when the local
normal-frame light direction is outward. The local speed stays one. The
coordinate speed carries the shift. Light rides the space.

The room-to-room part is in the inheritance theorem, `nsc-nested-qualities.md`:
regional eigenvalues obey `omega_n,k = Lambda_n * hat_omega_k(Theta)`. The
dimensionless spectrum is inherited unchanged; the per-room scale is not.
Same law in every room, different register.

So all three pieces you remembered are real, are recorded, and are in the
paper: local constancy (§2.1), the shift carrying the light (§2.2, transport
record), and the per-room scale (inheritance theorem).

## Coverage: what the foundation paper owns

Word-boundary counts on `paper/nested-space-cosmology.md`, after correcting a
scan that had matched `traction` inside "subtraction" and `RN` inside "return":

| in the paper | count |
|---|---|
| gradient | 6 |
| balance | 4 |
| shift | 7 |
| memory | 6 |
| retarded | 11 |
| parent-child | 2 |

## Coverage: what it does not contain

| not in the paper | count | lives in the lab as |
|---|---|---|
| two-force | 0 | `nsc-rn-gradient-bubble.md`, 4 Oct |
| collar | 0 | `nsc-discovery-parent-traction.md` |
| nested pair | 0 | `nsc-nested-parent-child.md`, 1 Oct |
| null expansion | 0 | `nsc-spherical-null-expansion.md`, 1 Oct |
| turnover | 0 | ledger, 30 Sept |
| Reissner / Reissner-Nordstrom | 0 | the whole RN radial programme, 3-4 Oct |

## The gap, by category

**Post-1-October, in neither document.**

- the finite parent-child pair, `MEASURED_CONFIRMED_NESTED_PAIR`: four
  reciprocal physical effects and a same-trajectory local reduction meeting
  one-percent targets under one evolving metric;
- the two-force gradient balance and its verdict, including the withdrawal
  of the invented throat-slope limit and the balanced-pair replacement;
- the RN radial chart: mass shell, compression ceiling of `2.05` at `97.3`
  percent of the trapped width, the simple pole at `r_minus`, the geometry
  map and its three regions;
- the local-response successors beyond v1.

**The companion paper's own domain.** The v0.27.0 companion owns the periodic
spherical realization: six Gaussian columns, conformal gauge, four matched
runs to `T = 0.3`, regional transfer, memory and drive. That is the 30
September to 1 October chain and it is covered.

## What this audit changed

The first pass reported `traction: 16` and `RN: 31` and would have concluded
that the collar and the RN programme were already in the paper. Both were
substring artefacts: `contractions`, `subtraction`, `extraction`, and
`return`. With word boundaries the counts are zero. The corrected reading is
that neither is present.

This is the fifth time in this stretch that a statistic was read off the
wrong matching rule. It is recorded here rather than quietly fixed.

## To do, concrete

1. Decide the next foundation version and its date. Everything in the gap is
   post-preprint and needs a home.
2. The gradient-balance record is the smallest complete addition: it already
   has a note, a module, a figure and a verdict, and its verdict is negative,
   which is a result.
3. The nested pair is the largest claim not yet in a document and is the one
   the programme calls confirmed. It needs a section, not a footnote.
4. The RN radial programme is a different construction from the periodic
   chart. It should be added as a second benchmark, clearly marked as a
   separate chart, or left in the lab and named as a limit.

## Scope

A coverage audit of two documents against a laboratory. It measures nothing
about the universe and adds no claim to either paper.
