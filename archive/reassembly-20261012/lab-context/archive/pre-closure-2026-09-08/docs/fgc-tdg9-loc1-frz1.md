# FGC-1-TDG9-LOC1-FRZ1: exact failed-cubic localization freeze

`FGC-1-TDG9-LOC1-FRZ1` is the prospective, premise-only successor to sealed
[`FGC-1-TDG9-AR1-PREF1`](fgc-tdg9-ar1-pref1.md) commit
`e6052fe62c4e8eb54b64f92fe501a923f3d0e087`. PREF1 proved ten exact-discrete
contraction failures on RCV3 retry widths `3`, `4`, and `5`; it did not say
where those failures attain their extrema or which Hermite inputs own them.

LOC1 authorizes one read-only answer to that narrower question. It consumes
exactly the ten sealed failure occurrences for `RK4-2049`, reconstructs the
same binary64 proposal surfaces independently of the AR1 runner and TDG8
persisted-retry wrapper, and inspects `4,088` D01 plus `8,176` D12 base cubics
per occurrence: `122,640` base cubics total. Each base cubic is split exactly
into endpoint-state `V`, width-scaled-RHS `S`, and complete `C=V+S`, so the
executed explanatory workload is `122,640` cubics per family and `367,920`
component cubics total. Candidate ceilings are separately frozen at `490,560`
for each of `V`, `S`, and `C`, or `1,471,680` in aggregate.

For each D01 and D12 maximum, two standard-library exact routes must preserve
the complete co-maximizer set or return a nonunique/interval-inconclusive
classification. A result records retry, channel, subinterval, owned row,
global index, exact radius identity, row region, endpoint or stationary-root
interval, exact absolute-value enclosure, the separately localized `V`, `S`,
and `C` maxima and threshold classes, cancellation/reinforcement and
value-only/slope-only/mixed-only ownership classes, restricted-parent versus
child endpoint data, endpoint updates, and the upward/downward binary64 ULPs
of contributors. ULPs are explanatory evidence; they never enter a tolerance
or pass rule. The original TDG6 `3/2` admission rule and its exact equality
precedence remain unchanged.

| Owned rows | Global rows | Region |
|---:|---:|---|
| `0` | `1` | centre-parity derivative row |
| `1..2` | `2..3` | inner stencil transition |
| `3..2040` | `4..2041` | smooth interior |
| `2041` | `2042` | outer dissipation transition |
| `2042..2043` | `2043..2044` | outer-projector-influenced |

The direct centre row `0` and fixed outer rows `2045..2048` are excluded from
the owned surface. Stage-2 pre/post-projector source decomposition and an
SSPRK3 comparator are not authorized by this freeze.

The raw namespace
`runs/fgc-2-sf1/tdg9-loc1/rcv3-retries-3-5` must be absent before the authority
commit and may later contain only canonical `manifest.json` and
`terminal.json`. Ordinary verification is compact and store-blind. No fourth
width, resource escalation, threshold change, state commit, continuation,
calibration, candidate branch, mechanism result, or physical inference follows.

```bash
make fgc-tdg9-loc1-frz1
make verify-fgc-tdg9-loc1-prelaunch
python3 scripts/run_fgc_tdg9_loc1.py --authority-commit <commit>
```
