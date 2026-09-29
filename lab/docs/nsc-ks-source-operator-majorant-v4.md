# Source insertion of the replayed direct-vacuum window

The v4 record replays the v3 aggregate and the saved direct-window traces, then
inserts those 87 unweighted signed pairs once. It does not solve, apply a
quadrature weight twice, or close the local gate. The independent read-only
review `82ecd706` was still running and had returned no terminal finding when
this page was written, so the bound is not an accepted final result.

Family `(14, 1)` now has 1042 authenticated signed rows: the previous 868 plus
174. Subtracting those keys from
`RetainedUpstreamArchive.family_entries((14, 1))` leaves 126 signed rows.
The family census is 1168. The v3 keys and the 174 direct identities are
disjoint. The upstream-budget component stays null, quadrature and numerical
field error stay out, and the local gate stays OPEN.

## Partial N and beta

Each archived weight is applied once to the unweighted signed total. The
continuous insertion of the covered rows is the stored dyadic upper bound

$$
N\le 444149617721196854374302166047853765198706369052164520085\times 2^{-227},
$$

$$
\beta\le 5645659491320673260464406121917161460134265034465380648797\times 2^{-231}.
$$

Read as decimals, those same bounds are about \(2.0593031173461625\times 10^{-12}\)
and \(1.6360089773354265\times 10^{-12}\). They are larger than the v3 partial
bounds because 174 previously omitted rows are now inside the sum. They are not
a bound on the remaining 126 rows, on every family, or on the full upstream
budget.

The stored weighted moments are mantissa
`5179819082872626942325109018578143493572420048632658129773` at exponent `-226`
and mantissa `443149945646076338904207074131538877453057508320262558083` at
exponent `-219`.

## Direct-window input

`results/development/nsc-direct-source-window-v1` holds 87 row directories and
174 files. Replay `--check` authenticated all 87 positive rows and their
negative partners, captured no new row, and left the checkpoint bytes
unchanged. Its process time was 115.087875 seconds. The maximum completed
signed total upper is mantissa `570057742100263` at exponent `-85`. The
production resume that wrote the rows used row budget 87 and CPU budget 600;
its capture process time was 147.245454 seconds, and its row records match the
replay.

The added positive rows are `group14/low16_1` rows 16 through 47 and
`group14/low32_1` rows 41 through 95.

## Remaining 126 signed rows

Every range below is a contiguous archived row interval. Each positive row
contributes its archived negative partner at the exact negation of that
energy. No new solver or pilot was run.

The existing subgap phase-tube owner
`subgap_mixed_phase_transport_bound` requires the archived mass to exceed the
positive energy. The saved low16 row 15 witness is that transport plus the
Bloch witness in `nsc_subgap_row_witness`. These remaining rows reuse that
method only where the mass hypothesis holds. They do not yet have production
witnesses.

| Panel | Positive rows | Signed rows | Archived energy hex, first to last |
|---|---:|---:|---|
| `group14/low16_1` | 1–14 | 28 | `0x1.c60a99e906500p-8` to `0x1.f1cfab30b7cd8p-3` |
| `group14/low32_1` | 0–31 | 64 | `0x1.010cf92bcf8b8p-2` to `0x1.ff79836a183a4p-1` |
| `subgap8/14_1` | 0–7 | 16 | `0x1.02e6bb940d290p+0` to `0x1.8f38f9b035a88p+0` |

That is 54 positive rows and 108 signed rows, all strictly below
\(\pi/2=0x1.921fb54442d18p+0\). Low16 rows 0 and 15 stay inside the previous
868 and are not in this remainder. The `subgap8/14_1` panel is stored as eight
inherited real subgap rows whose accuracy is OPEN. Its row 7 has only the
existing wider-tube pilot, tube \(10^{-3}\) and 16 subdivisions, not a saved
production witness.

Nine further positive rows are outside both that phase-tube hypothesis and the
completed direct window. `group14/low32_1` rows 32–40, and their negative
partners, are 18 signed rows. Their energies run from `0x1.92f9814c33825p+0`
to `0x1.fae56f00f6a85p+0`: each is at least \(\pi/2\) and strictly below 2.
The phase-tube owner rejects this side of the mass. The recorded direct window
starts at energy 2. No saved witness is claimed for these 18 rows.

## Files and checks

The aggregate owner is `scripts/derive_nsc_ks_source_operator_majorant_v4.py`.
The result is `results/development/nsc-ks-source-operator-majorant-v4.json`,
sha256 `4b2462987f7784ee77b58687521226f274400ac6b0924bdf493428037aabbdf1`.
The direct-window content sha256, over each relative path and file bytes, is
`1a4f1b592931a7a2804a6a502e90483f306672af4be5384f9ebd78f19a7b27b9`.

From the repository root:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_direct_source_window.py --check --output results/development/nsc-direct-source-window-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_source_operator_majorant_v4.py --check
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_source_operator_majorant_v4.py
```

In this session the direct-window check, v4 `--record`, v4 `--check`, and the
two v4 tests all exited 0. The tests passed in 123.96 seconds. Checkpoint
bytes were not rewritten.
