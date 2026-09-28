# FGC-1-SRC2-PREF11: exact diagnosis of the captured affine-source wall

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** machine-reproduced, point-local arithmetic diagnosis; no repaired
runtime source, successor protocol, eligible GR-0 case, candidate trajectory,
or physical mechanism result

```text
PROTO11_source_wall_replayed = true
CAP1_frozen_arithmetic_routes_cleared_the_wall = false
captured_GR0_point0_exact_affine_system_nonsingular = true
captured_GR0_point0_binary64_residual_evaluator_incomplete = true
PROTO12_frozen = false
FGCQR_holdout_execution_authorized = false
```

## What the expensive replay established

[SRC2/CAP1](fgc-src2-cap1.md) was run from clean immutable commit
`241e9b99a44985202ad35bdfaf9e6e1cf2d2c3d6`. It rebuilt only the implicated
GR-0 amplitude-`5/2`, `RK4-8193` member. The replay reproduced the first
evolved common event, all 164 canonical source-only rejection records, the
terminal accepted state, and the typed terminal exception. It made 2,130
source-operator calls, of which 687 failed, and matched the target source call
twice. It did not read an SGB-L or FGC-QR outcome.

The raw bundle is immutable by hash but remains under the ignored run
namespace. PREF11 tracks a compact canonical point fixture containing the
literal binary64 bits of the first positive-radius lower jet, the PROTO11
acceleration and residual, and the captured forward affine coefficients. A
clean clone can therefore reproduce the point-level arithmetic conclusion;
when the full raw bundle is present, PREF11 additionally proves that every
compact value is bitwise the value extracted from point zero of that bundle.

## What CAP1's frozen attacks did not establish

At coordinate radius `1/64`, the unchanged complete direct REF1 evaluator
reported the baseline infinity residual

```text
1.0659145473163184e-12
```

against the unchanged strict `<1e-12` rule. The captured affine
reconstruction residual was about `3.55e-15`, and the worst grid-wide captured
infinity-condition estimate was about `2.86e6`, below the frozen `1e10` stop.
Nevertheless, none of CAP1's prospective routes cleared the complete
evaluator:

- power-of-two row and column equilibration;
- symmetric coefficient probes at scales `1`, `16`, and `256`; or
- selected exact-rational correction solves of the captured binary64
  coefficients.

The PROTO11 baseline remained the best frozen CAP1 candidate. This is an
important negative result about those declared repairs. PREF11 does not
retroactively add a probe, relabel one of them as passing, or weaken the raw
gate.

## Exact dyadic oracle

CAP1 captured a binary64 lower jet. Every finite binary64 number is an exact
dyadic rational. PREF11 therefore reconstructs those literal bits as
`Fraction` values and sends them through the already certified exact
two-jet/tensor REF1 implementation. It evaluates the same GR-0 action,
modified-harmonic reference connection, and six unredefined equations with
no discretization or floating arithmetic inside the point evaluator.

Because the six rows are affine in the six coordinate-time accelerations,
PREF11 evaluates exact zero and unit seeds, forms the exact six-by-six matrix
and constant term, and solves them by exact Gaussian elimination. The exact
matrix is nonsingular, and the resulting rational acceleration zeros all six
exact residual rows identically.

The exact root is then rounded once, component by component, to binary64. The
rounded acceleration is approximately

```text
(-2.0722242346857163e-14,
  2.8695350241251600e-13,
  2.0939050943367178e-13,
 -3.0918808477578940e-15,
  0,
  0)
```

When those rounded binary64 values are reinterpreted as exact dyadic
rationals, the complete exact residual infinity is approximately
`5.00e-28`. The same binary64 acceleration is falsely rejected by both
floating implementations: the fast direct evaluator reports approximately
`7.47e-12`, and the generic vectorized evaluator reports approximately
`6.69e-12`.

The comparison is deliberately asymmetric: the exact tensor evaluator is the
oracle for this captured dyadic input. A floating acceleration that merely
pleases one floating evaluator while worsening the exact residual is not an
admissible repair.

## The bounded conclusion

For this captured GR-0 point and these unchanged equations:

1. the affine source system is nonsingular;
2. an exact acceleration root exists;
3. its binary64 rounding satisfies the exact captured equations by more than
   fifteen orders of magnitude below the declared `1e-12` gate; and
4. both current binary64 residual paths mismeasure that root above the gate.

The specific point-zero wall is therefore a floating cancellation and
evaluator-instrumentation failure, not evidence of a physical or structural
breakdown of the GR-0 equations there. CAP1 initially distinguished the
linear-solve route from the complete evaluator; PREF11 now localizes the
failure to the latter. This is stronger than saying merely that the solver is
ill-conditioned, and narrower than saying that the entire calibration would
otherwise pass.

## Why this does not freeze PROTO12

PREF11 diagnoses one point on one rejected evolved source call. It does not
yet provide a cancellation-resistant runtime evaluator across the full grid,
all Runge--Kutta stages, both methods, or independent nontrivial controls. It
also does not clear the separate amplitude-`3` direct coarse-to-medium
`phi/Lambda` derivative-tail veto.

The next gate is therefore **FGC-1-SRC3**. It must prospectively derive a
stable evaluation order, compensated formulation, or bounded higher-precision
fallback that:

- converges to the exact oracle on this compact fixture and independent
  controls;
- preserves the unredefined equations, reference derivative map, branch,
  source threshold, grids, methods, CFL/retry rules, and physical inputs;
- retains typed failure outside its proven numerical envelope; and
- is frozen before any fresh GR-0 campaign is inspected.

Only after SRC3 passes may a separate protocol consider PROTO12. The
resolution-spectrum obstruction requires its own prospective study and may
not be tuned away in the source repair.

## Scientific boundary

PREF11 does not complete GR-0 calibration, establish a trapped interval,
open SGB-L or FGC-QR, test regulator activation, compute a complete affine-null
Raychaudhuri margin, validate the retained EFT, resolve a singularity, derive
a daughter domain or dark sector, vary a locally measured speed of light, or
confirm or reject the general gradient programme. It repairs the epistemic
status of one numerical wall: the wall belongs to the current measuring
instrument, not to the exact captured equations.

## Reproduction

```bash
python3 scripts/reproduce_fgc_src2_pref11.py --check
```

The command does not launch the forty-eight-minute CAP1 replay. If the raw
bundle is available, it verifies and cross-checks it; otherwise it reproduces
the exact point result from the tracked compact bit fixture and immutable hash
bindings.
