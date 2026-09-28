# FGC-1-TDG7-PREF23: independent stage-safe lattice binder

`FGC-1-TDG7-PREF23` independently rederives the stage-safe exact-binary64
shared-lattice theorem frozen by `FGC-1-TDG7-FRZ1`.  It binds the sealed TDG7
freeze and the historical TDG6/PROTO14 sources at their immutable commits. It
does not import the TDG7 design helper as its theorem.

## Independent result

The independently reconstructed/prospective **`4Q` endpoint-only
alternative** is insufficient as stated for actual RK4 and SSPRK3 stage
execution: an odd-`Q` fine width leaves a c=`1/2` time at a half-`Q` coordinate
that binary64 cannot represent exactly. The selected fixed `8Q` quantum is a
simple sufficient choice: it makes every fine width an even multiple of `Q`,
so every stage time at c in `{0, 1/2, 1}` is an exact shared coordinate and
each positive fine midpoint is strictly interior. This does not claim that no
other explicitly parity-constrained construction could be safe.

For CAL11 event 24, the independent witness reconstructs `Q=2^-52`,
`G=8Q=2^-49`, the floor-selected macro width
`3664984285035/562949953421312`, and its four exact fine intervals. It also
reconstructs the historical one/two/four old guard failures from the frozen
expression rather than treating the freeze helper as the proof. The immutable
historical TDG6 source is also parsed directly: the binder requires the exact
`step=width/count`, `start+index*step`, `expected_final=start+width`, endpoint
guard, adjacent-pair guard, `(left,right)` loop target, and both historical
`ValueError` actions/messages before that reconstruction is promoted as a
source-bound result.

## Exact theorem domain and proof

The theorem is deliberately limited to one forward PROTO14 coordinate step:

```text
23/16 <= x < t <= 32,
Q = max(ulp(x), ulp(t)),
x mod Q = 0,
t mod Q = 0,
(t-x) mod (8Q) = 0.
```

With requested cap `C`, define the non-upward selection budget and macro step

```text
B = min(C, t-x),
W = floor(B/(8Q)) * 8Q,
h = W/4.
```

For positive selected width, `W = 8nQ` for an integer `n >= 1`; hence

```text
h = 2nQ,
h/2 = nQ.
```

Every shared boundary and every c=`1/2` stage coordinate is therefore a
`Q`-lattice coordinate. On this positive interval, `Q` dominates the binary64
ULP at every intermediate coordinate. Each exact endpoint difference is also
a `Q`-multiple no larger than the target, so it is itself representable;
correctly rounded binary64 subtraction therefore returns that exact
difference. Division by two retains a `Q`-multiple, and the resulting lattice
sums are likewise representable. Thus the actual float expression used by the
engine,

```text
left + (right - left) / 2.0,
```

equals the exact interior midpoint for all one-, two-, and four-path
substeps. The result is a sufficient fixed-quantum construction, not a
uniqueness claim: a different construction with an explicit equivalent parity
condition could also be safe.

The theorem parses the frozen historical numerical-engine AST, rather than
assuming its abscissae. It requires the live method branches to retain RK4 at
`start`, `start+dt/2`, `start+dt` and SSPRK3 at `start`, `start+dt`,
`start+dt/2`, with a unique post-branch candidate-endpoint evaluation at
`final_time`. Duplicate/dead-code stage decoys, missing branches, moved
midpoints, or a changed endpoint call are rejected.

The deterministic sweep of 3,014 successful plans and 21,098 stage triplets
is regression evidence across the declared event lattice and binade controls.
It is not an exhaustive computational proof over every binary64 input; the
universal conclusion comes from the exact `Fraction` argument above together
with the enforced domain and alignment preconditions.

Typed controls reject nonfinite coordinates, overflowed integers, nonpositive
caps, targets that are not ahead or not stage-lattice aligned, off-`Q`
endpoints, too-small caps, and below-minimum aligned widths. Mutation controls
also reject a forged out-of-envelope plan, odd-`Q` fine widths, cap rounding
up, endpoint/midpoint changes, historical guard/source-stage changes, and
claim promotion. These are executed during canonical certificate construction:
the record does not serialize a mutation-control success flag without running
the corresponding lineage, `8Q`/`4Q`, witness, engine-stage, TDG6 source-body,
claim-boundary, or typed-cap attack.

The historical first macro endpoint and the selected conservative first macro
endpoint intentionally differ. The event target remains unchanged. CAL11
therefore remains an invalid numerical-runtime result, not a temporal
admission, constraint, source, health, GR-0, candidate, or physical result.

## Narrow authorization

This independent theorem/binder authorizes only implementation and synthetic
qualification of a future TDG7 stage-safe runtime repair. That implementation
must use one shared plan for one/two/four paths and actual RK4/SSPRK3 stage
times, preserve the physical equations, source, projector, spatial operator,
and TDG6 threshold, and issue any coordinate-lattice stop before shadow or
real-state work.

It does not implement a repair, mutate or resume historical PROTO14, create a
new protocol or output namespace, authorize a fresh GR-0 calibration, select
an eligible case, open SGB-L/FGC-QR/DEF1, or establish an EFT, transition, or
physical claim. No campaign checkpoint/state/tracer/debit arrays are loaded or
resumed.
