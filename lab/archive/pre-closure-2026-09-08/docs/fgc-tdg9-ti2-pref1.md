# FGC-1-TDG9-TI2-PREF1: independent TI2 terminal binder

`FGC-1-TDG9-TI2-PREF1` independently binds the completed two-leaf TI2
diagnostic produced under immutable authority commit
`963f772c7fd5de676c37c3ee385a930701df6341`. The one-time construction path
uses no-following, race-checked reads for the canonical raw manifest and
terminal, verifies their exact SHA-256 identities, authenticates the authority
lineage and six committed authority blobs, and checks the scientific campaign
store before and after against its sealed 115-leaf identity.

The binder imports no TI2 runner decision code. It restores the exact three
RK4-2049 retry predecessors, verifies their state, time, grid, SBP4 operator,
projector, transaction, tracers, temporal ledger, journal, widths, and
transaction hashes, then calls the lower-level TDG6 compositor with only the
tableau selector changed to SSPRK3. It constructs 21 read-only shadow proposals
and 84 SSPRK3 stage/endpoint records. It never calls temporal admission,
commits a fine path, advances a state, writes the campaign store, or opens a
candidate branch.

For the ten sealed occurrence selections, the binder rebuilds all `122,640`
base and `367,920` value/slope/complete cubics and independently repeats both
exact localization routes. All 60 occurrence/component/level checks require:

- exact polynomial and candidate counts;
- identical route classifications and stationary-count digests;
- identical survivor keys with overlapping exact intervals;
- byte-identical reproduction of the corresponding raw primary and v2 exact
  evidence;
- independently recomputed contraction and component ownership; and
- exact agreement with every raw reduction and aggregate count.

All 60 route checks pass. The candidate totals are `245,280` for `value_V`,
`397,721` for `slope_S`, and `332,438` for `complete_C`, or `975,439` total,
inside the frozen ceilings. Fifty routes have a unique maximum; ten retain a
conservative nonunique-or-interval-inconclusive survivor set. The 93 retained
survivor candidates comprise 60 exact endpoints, 17 exact endpoint
deflations, and 16 adaptive nonsquare isolations.

The bounded outcome is mixed and exact. Four original complete-C failures
clear under the tableau-only shadow: retry-3 `u:alpha`, retry-3 `u:R`, retry-5
`u:v`, and retry-5 `q:R`. Six persist. Retry-4 `u:alpha` and retry-4
`u:lambda` become endpoint-state-owned; retry-4 `u:R`, retry-5 `u:alpha`,
retry-5 `u:lambda`, and retry-5 `u:R` retain independent value and slope
failure. There are no threshold-inconclusive occurrences.

This proves that the frozen discrete obstruction is partly tableau-sensitive
and partly persistent across the original RK4 and this SSPRK3-tableau shadow
on the exact same restored state/grid/SBP4 operator. It is not a production
SSPRK3 comparator, independent-method
agreement, a temporal remedy, state progression, common-event completion,
GR-0 calibration, candidate dynamics, mechanism evidence, or physics.
Changing ownership alone is not clearance, and this binder selects no
successor.

Ordinary verification reads only the tracked compact certificate. The raw
namespace, campaign store, restored predecessors, and shadow proposals are
opened only by the explicit one-time live construction command.

```bash
make fgc-tdg9-ti2-pref1
make verify-fgc-tdg9-ti2-pref1

# Explicit one-time construction only; not an ordinary verifier:
python3 scripts/reproduce_fgc_tdg9_ti2_pref1.py --live
```
