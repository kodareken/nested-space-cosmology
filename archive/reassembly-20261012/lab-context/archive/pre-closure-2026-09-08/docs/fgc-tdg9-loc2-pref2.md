# FGC-1-TDG9-LOC2-PREF2: independent LOC2 terminal binder

`FGC-1-TDG9-LOC2-PREF2` independently binds the completed, two-leaf LOC2
diagnostic produced under immutable authority commit
`1cd78ee20100ae4ae6495161146c6f2fcf7de1bf`. The one-time construction path
uses a no-following, race-checked reader for the canonical raw manifest and
terminal, verifies their exact SHA-256 identities, authenticates the authority
lineage and six committed implementation blobs, and checks the read-only
campaign store before and after against its sealed 115-leaf snapshot.

The binder imports no LOC2 runner decision code. It reconstructs the three
historical TDG6 proposal surfaces from the authenticated store, selects exactly
the ten PREF1 failure occurrences, rebuilds their `122,640` base and `367,920`
value/slope/complete cubics, and independently applies the primary and v2 exact
localizers. For every one of the 60 occurrence/level/component evaluations it
requires:

- the same polynomial and candidate counts;
- the same route classification;
- the same per-polynomial stationary-count digest;
- identical survivor keys;
- overlapping survivor and global exact intervals;
- the exact raw primary and v2 evidence;
- the exact maximizer-input decomposition;
- the exact contraction classification; and
- the independently recomputed component-ownership reduction.

All 60 route checks pass. The candidate totals are `245,280` for `value_V`,
`391,791` for `slope_S`, and `326,164` for `complete_C`, or `963,235` in all,
below the frozen per-family and aggregate ceilings. Twenty-nine localizations
have a unique maximum, while 31 conservative nonunique-or-interval-
inconclusive routes retain multiple degenerate exact-point candidates. Across
the 3,062 retained survivor candidates, v2 uses 3,013 exact endpoints, 28 exact
endpoint deflations, and 21 adaptive nonsquare isolations. No ULP tolerance or
boundary clipping enters the classification.

The bounded outcome is narrow but useful: all ten historical sufficient
contraction failures remain exact, and in every occurrence both the endpoint-
state component and the width-scaled endpoint-RHS component independently
fail. None is a mixed-only cancellation artifact on these frozen samples.
This authenticates the repaired localization instrument and the component
ownership of the sampled discrete obstruction. It does not select the next
temporal remedy, prove continuum-order loss, finish a common event, calibrate
GR-0, or say anything about FGC-QR dynamics or physics.

Ordinary verification reads only the tracked compact certificate. The raw
namespace and campaign store are opened only by the explicit one-time live
construction command.

```bash
make fgc-tdg9-loc2-pref2
make verify-fgc-tdg9-loc2-pref2

# Explicit one-time construction only; not an ordinary verifier:
python3 scripts/reproduce_fgc_tdg9_loc2_pref2.py --live
```
