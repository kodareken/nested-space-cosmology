# FGC-1-SRC4-VEC1 — tensor-contracted SRC3 source evaluator

## Result

`FGC-1-SRC4-VEC1` derives and machine-checks a faster binary64 instrument for
the unchanged GR-0 specialization of the ACT1/VAR1/REF1 equations. It replaces
SRC3's explicit Python loops over spacetime indices with NumPy tensor
contractions. It does **not** replace SRC3 as historical evidence: SRC3 remains
the independent, immutable differential oracle, and RSP1 remains reproducible
from the exact evaluator that generated it.

The result is deliberately narrow:

> SRC4 evaluates the same reference-covariant connection-difference and Ricci
> identities, preserves the strict raw residual and kinetic gates, agrees with
> exact-dyadic and independent SRC3 controls within frozen binary64 bounds, and
> may be consumed only by a separately frozen successor protocol.

It freezes no PROTO12 contract, opens no trajectory, and says nothing about
collapse or the FGC-QR mechanism.

## Why this gate exists

RSP1 required about `7571.79` wall seconds to move six members only from
`t=0` to `t=1/16` on the `4097 -> 8193 -> 16385` ladder. A naïve extrapolation
of that high-grid diagnostic to the old `t=32` calibration endpoint would
measure the same local tensor identities for weeks. That is not stronger
science; it is avoidable interpreter overhead.

Profiling localized most of one source call to SRC3's auditable nested loops in
the construction of

\[
C^a{}_{bc}
=\frac12 g^{ad}
\left(\bar\nabla_b h_{dc}+\bar\nabla_c h_{db}-\bar\nabla_d h_{bc}\right)
\]

and

\[
R_{bd}
=\partial_a C^a{}_{db}-\partial_d C^a{}_{ab}
+\bar\Gamma C+C\bar\Gamma+CC.
\]

The point axis and the seven affine-seed states are independent. Contracting
only the four-dimensional tensor indices therefore changes the execution
graph, not the mathematical map.

## Frozen evaluator boundary

SRC4 preserves all of the following:

- the unredefined ACT1/VAR1 GR-0 REF1 equations;
- the flat spherical reference geometry;
- PROTO11's reference-balanced radial derivative map;
- the six-field ADM state and equation order;
- the auxiliary-cone factors `4` and `9`;
- the affine zero-plus-unit seed construction;
- the strict complete-residual gate `||R||_infinity < 10^-12`;
- the kinetic condition ceiling `10^10`;
- the maximum sixteen complete residual-refinement iterations;
- the accepted root branch and fail-closed error behavior.

It adds no epsilon floor, tolerance, fitted coefficient, higher-precision
fallback, damping term, constraint cleaning, or physical assumption.

## Independent controls

### Exact spherical reference

At every frozen reference radius, all complete metric/scalar residuals,
constraints, invariants, and accelerations remain bitwise zero.

### PREF11 exact-dyadic point

The compact point that diagnosed the original false rejection has SRC4
complete residual

\[
2.8272774843121063\times10^{-27},
\]

and remains within sixteen ULP of the rounded exact-dyadic affine root. Its
exact-dyadic residual remains below `10^-24`.

### Independent nontrivial controls

The two pre-existing `DYADIC-A` and `DYADIC-B` controls are evaluated three
ways:

1. by the exact rational oracle;
2. by SRC3's explicit-index binary64 implementation;
3. by SRC4's tensor-contracted binary64 implementation.

The normalized residual differences must remain below `2^-46`. The connection
difference, its derivative, and Ricci tensor also have separately frozen
differential bounds.

### CAP1 full-grid control

When the optional raw CAP1 bundle is present, both implementations evaluate
all 8,192 positive-radius points. SRC4 returns

- complete residual infinity `2.4197549346183234e-14`;
- kinetic condition maximum `2859227.196975944`;
- acceleration difference from SRC3 no larger than `2e-15`;
- the same strict source and condition-gate decisions.

The raw bundle remains optional for a clean clone. Its expected compact hashes,
strict-gate decisions, SRC3/SRC4 array hashes, and exact cross-backend deltas are
frozen in the tracked configuration. When the bundle is present the reproducer
recomputes and checks all of them; when it is absent the canonical certificate
remains byte-identical and binds those values to the immutable raw-fixture hash.
Partial or hash-mismatched evidence fails closed.

### Adversarial controls

- one binary64 bit changed in the compact lapse changes the accepted root;
- nonpositive radius fails;
- a deliberately impossible condition ceiling fails;
- invalid auxiliary-cone ordering fails.

## Performance evidence is not a physics result

Local measurements on this Mac found about `3.2x` to `3.4x` faster complete
evolved-state source calls and about `4.7x` on the 8,192-point CAP1 fixture.
Those timings are platform and load dependent. They are contextual engineering
evidence only and are intentionally absent from the scientific pass/fail gate.

The authorization rests on algebraic identity, exact/differential controls,
strict gate preservation, and adversarial behavior—not on a stopwatch.

## Machine contract

The canonical result is generated by:

```bash
python3 scripts/reproduce_fgc_src4_vec1.py
```

and checked by:

```bash
python3 scripts/reproduce_fgc_src4_vec1.py --check
```

The result binds the configuration, this derivation, the new implementation,
the immutable SRC3 and RSP1 results, the SRC3 source, and both exact-control
fixtures.

## Epistemic boundary

This certificate establishes an equivalent numerical instrument under its
declared controls. It does not establish that every possible binary64 state
will round identically to SRC3. Bitwise identity is neither expected nor used:
different valid contraction orders may differ at roundoff scale. What must be
preserved is the unchanged equation, branch, strict gate, exact-reference
behavior, exact-oracle agreement, and bounded differential behavior.

In particular, SRC4 does not imply:

- PROTO12 is not frozen;
- that the generic RSP1 all-field raw tail-ratio aggregate passes;
- that a GR-0 amplitude is eligible;
- that SGB-L or FGC-QR may run;
- that collapse defocuses;
- that retained-EFT validity, singularity resolution, a child domain, a dark
  sector, or a variable locally measured speed of light has been derived.

The next gate is `FGC-1-PRO12-FRZ1`: a prospective protocol decision that must
separately own the spectral-conditioning interpretation exposed by RSP1.
