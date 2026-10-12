# FGC-1-PRO11-FRZ1: outcome-neutral well-balanced map freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** immutable premise-only protocol boundary; no runtime, trajectory,
calibration, candidate, or mechanism result

```text
PROTO11_frozen = true
PROTO11_successor_runtime_implemented = false
PROTO11_fresh_GR0_dynamic_calibration_authorized = false
FGCQR_holdout_execution_authorized = false
```

## Why a successor is justified

The completed PROTO10 campaign did not reach its first evolved common event.
Both frozen matter amplitudes stopped in `RK4-8193` on the same source-only
residual floor, with exact rollback and no non-source failure. CAL7/PREF9 then
localized that floor to the first positive `alpha` row in the exact centre
vacuum buffer. Its magnitude scales as binary64 epsilon divided by radius
squared and is independent of the matter amplitude.

CAL7 also evaluated an outcome-neutral control on all twelve initial inputs.
Differentiating departures from the exact Minkowski spherical reference
reduced the largest fine primary residual from about `7.86e-13` to about
`7.58e-14`, while retaining the unchanged raw `1e-12` source gate. That was a
`t=0` control, not a trajectory. PROTO11 freezes the smallest separate
numerical premise warranted by that evidence before any successor runtime or
outcome exists.

## The only permitted numerical-map change

In canonical field order

```text
(alpha, shift, lambda, R, phi, chi),
```

the fixed reference is

```text
u_ref = (1, 0, 1, r, 0, 0)
q_ref = (0, 0, 0, 1, 0, 0).
```

PROTO11 freezes exactly three algebraically equivalent semidiscrete rules:

```text
initial q             = D_h (u - u_ref) + q_ref
source q_r            = D_h (q - q_ref)
reduction constraint  = q - (D_h (u - u_ref) + q_ref).
```

The common-event constraint calculation must use the same representation.
The reference has no fitted, evolved, amplitude-dependent, or
outcome-dependent parameters.

This is a well-balanced discretization of the same derivatives. It is not a
new term in the continuum equations. It is not damping, constraint cleaning,
or subtraction from a physical observable. It does not redefine the state:
`q` remains the independently evolved approximation to the physical radial
derivative of `u`.

In particular, the evolution remains

```text
partial_t u = p
partial_t q = D_h p.
```

No interior `q` value may be reprojected after initialization or after a
Runge--Kutta stage. A reduction defect must remain visible rather than being
erased by the map.

## What remains byte-for-byte or semantically unchanged

The immutable PROTO10 hash continues to own every unlisted rule, including:

- the ACT1/REF1 equations and the complete unredefined source residual;
- the `chi` profiles, amplitude order, `phi` seed, couplings, and grids;
- the `[2049, 4097, 8193]` ladder and both numerical methods;
- the raw `1e-12` source ceiling, affine source solve, kinetic gate, and
  monotonic-residual requirement;
- CFL, dissipation, half-step retries, retry counts, and minimum step;
- constraint magnitudes and orders;
- the direct `1/4` spectral ceiling, direct coarse veto, absolute spectral
  budgets, and complete PROTO10 conditioning guards;
- boundary, health, scale, trapped-sign, temporal, affine-null, and
  robustness rules;
- all candidate, retained-EFT, and physical nonclaims.

No epsilon floor or new tolerance is introduced. The PROTO10 campaign remains
immutable diagnostic history and cannot be relabelled as a PROTO11 result.

## Runtime burden deliberately left open

HLT9 must independently implement and attack the frozen map before any fresh
run may be authorized. It must at least:

- rebuild and serialize all twelve projected input hashes;
- prove bitwise exact Minkowski preservation on every frozen grid and method;
- cross-check the map against the CAL7 control without importing an outcome;
- show that nontrivial generic REF1, source, and reduction defects remain
  visible;
- apply the same derivative representation in initialization, source
  evaluation, and common-event constraints;
- prove that the interior evolved `q` field is never silently reprojected;
- reconstruct every unchanged physical input and every PROTO10 spectral gate;
- use fresh absent `proto11` calibration and holdout namespaces.

Until every HLT9 predicate passes, PROTO11 cannot launch even GR-0.

## Scientific boundary

PROTO11 repairs a test representation, not the FGC-QR model. It does not show
that GR-0 will calibrate, that a trapped sphere will form, that FGC-QR is
healthy or activates, or that metric-null defocusing occurs. It does not
validate a retained EFT, resolve a singularity, create a child domain, derive
a dark sector, or vary a locally measured speed of light.

The freeze is valuable because it removes one identified numerical ambiguity
without moving any physical goalpost. A fresh failure under this map must
still be classified by its actual first failed premise; a fresh success must
still earn all later authorization gates.

## Reproduction

Before the fresh PROTO11 namespace is ever used:

```bash
python3 scripts/reproduce_fgc_pro11_frz1.py --check
```
