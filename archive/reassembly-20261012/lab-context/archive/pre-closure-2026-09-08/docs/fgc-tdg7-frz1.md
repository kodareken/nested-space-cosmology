# FGC-1-TDG7-FRZ1: exact binary64 subdivision-lattice design freeze

`FGC-1-TDG7-FRZ1` preserves CAL11's terminal invalid-runtime result while
freezing a prospective arithmetic diagnosis of its pre-shadow failure.  It
does not resume the PROTO14 checkpoint, read a physical-state array, construct
a TDG6 shadow path, mutate a runtime, create an output namespace, or execute a
calibration.

## Historical boundary

CAL11 binds the `RK4-2049` abort at common event `24`, after completed event
`23` at coordinate time `23/16`.  The historical TDG6 code formed, separately
for counts one, two, and four,

```text
step = width / count
boundaries[i] = current + i * step
```

and required each rounded adjacent subtraction to be bitwise equal to its
independently rounded `step`.  With CAL11's frozen requested cap, all three
old guards fail before source, projector, tracer, monitor, or state work.  The
TDG7 reproducer reads that source only from CAL11's immutable authorization
lineage; it does not call the historical campaign.

## Selected prospective route

For an increasing PROTO14 coordinate interval, TDG7 takes the larger binary64
ULP of current time and event target as `Q`, then uses the **stage-safe** macro
lattice `G = 8Q`, not merely endpoint-safe `4Q`.  A `4Q` macro may have an
odd-`Q` fine width.  Its c=`1/2` fine-stage time is then a half-`Q` coordinate,
which binary64 rounds or aliases rather than representing exactly.  Both RK4
and SSPRK3 actually evaluate at c in `{0, 1/2, 1}`, so endpoint equality alone
is insufficient.

Exact dyadic arithmetic floors the requested cap to a positive `8Q` multiple
and constructs one shared five-boundary partition.  The one-step path selects
boundaries `(0,4)`, the two-step path selects `(0,2,4)`, and the four-step
path selects `(0,1,2,3,4)`.  All path endpoints and every distinct c=`1/2`
stage coordinate are therefore exact, shared binary64 lattice coordinates.

For the frozen CAL11 event-24 witness,

```text
Q = 2^-52
G = 8Q = 2^-49
8Q ticks = 3664984285035
W = 3664984285035 / 562949953421312 = 0x1.aaa90b0fb5800p-8
fine = 3664984285035 / 2251799813685248 = 0x1.aaa90b0fb5800p-10
fine / Q = 7329968570070
reduction = 531 / 576460752303423488 = 0x1.0980000000000p-50
```

The selected five boundaries are

```text
0x1.7000000000000p+0
0x1.706aaa42c3ed6p+0
0x1.70d5548587dacp+0
0x1.713ffec84bc82p+0
0x1.71aaa90b0fb58p+0
```

and the four exact fine c=`1/2` coordinates are

```text
0x1.7035552161f6bp+0
0x1.709fff6425e41p+0
0x1.710aa9a6e9d17p+0
0x1.717553e9adbedp+0
```

The old independent four-path construction is preserved as historical
evidence, not reused as the selected partition.  Its old expected first macro
endpoint was `0x1.71aaa90b0fb5cp+0`; the conservative stage-safe first macro
endpoint is intentionally `0x1.71aaa90b0fb58p+0`.  The invariant is the
unchanged event target `24/16`, not equality of those two first-step endpoints.

The historical cap is rounded down, never up.  An unaligned target, a cap too
small to contain one macro quantum, or a minimum-width violation is a typed
`coordinate_lattice_limit_reached` numerical stop before any shadow or real
state operation.  It is not a temporal-admission, constraint, source, health,
or physical classification.

## Immutable scope

The binder hashes CAL11’s config, compact result, reproducer, document, and
diagnosis source at checkpoint `b59d214…`, requiring the current worktree to
match.  It separately reads the historical TDG6 runtime, PROTO14 runner,
run-plan, and authorization result from authorization commit `7534a16…`.

When the raw checkpoint exists, TDG7 verifies its whole-file SHA-256 before
reading only the ZIP member `metadata_utf8.npy`; it does not open state,
tracer, or debit arrays.  Complete absence is portable in a clean clone and
does not change the canonical certificate.

## Nonclaims and sequencing

This is a design freeze, not a runtime repair or an authority to patch
PROTO14.  The terminal checkpoint may not resume.  Runtime implementation and
qualification require a separately committed TDG7 binder; only afterwards can
a new protocol, new namespace, and fresh runtime authorization be considered.
No GR-0 eligibility or calibration, SGB-L/FGC-QR/DEF1 execution, retained-EFT,
transition, mechanism, or physical claim follows.
