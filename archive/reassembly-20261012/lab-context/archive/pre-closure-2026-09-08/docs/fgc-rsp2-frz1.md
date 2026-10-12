# FGC-1-RSP2-FRZ1: late-event one-finer-grid constraint freeze

**Authors:** Douglas Ek & ChatGPT 5.6 Sol

**Status:** prospective GR-0 numerical authorization; no new trajectory or
mechanism result

```text
RSP2_runtime_implemented = true
RSP2_execution_authorized = true
RSP2_outcome_read = false
PROTO13_frozen = false
FGCQR_holdout_execution_authorized = false
```

## Frozen question

CAL9/PREF13 reached late evolved GR-0 common events with both declared
amplitudes.  Its sole terminal comparator veto was the SSPRK3
`radial_momentum` finest-pair order.  For amplitude `3`, the effective
`4097 -> 8193` value was `1.4990947384363291`, below the unchanged required
`3/2` by `0.0009052615636708783`.

RSP2 does not reinterpret that value.  It prospectively asks whether the same
component passes or remains below `3/2` on the genuinely finer
`8193 -> 16385` pair at the same exact endpoint `t=23/16`.

Only one new member is advanced:

```text
branch:                 GR-0
amplitude:              3
method:                 SSPRK3 / D2-1 SBP
read-only predecessors: 4097, 8193
new member:             16385
endpoint:               23/16
target:                 radial_momentum
order gate:             p >= 3/2
```

The physical input, GR-0 equations, PROTO11 reference-balanced map, SRC4
source backend, CFL/retry transaction, boundary budget, evolution-owned
constraint domain, raw guards, and roundoff enclosure are unchanged.

## Interpretation fixed before launch

- A target pass means only that CAL9's miss was pre-asymptotic on the tested
  ladder.
- A target failure means only that it persists through the tested
  `8193 -> 16385` pair.
- Complete constraint admission is reported separately from the target.
- A runtime stop, failed premise, or invalid/nonconverged run is not classified
  as persistence.
- No fitted asymptote, threshold reduction, saturation rule, or post-outcome
  parameter change is permitted.
- Every result requires a separate post-run PREF14 reduction before PROTO13
  can be designed.

## Temporal provenance

The authorization binds the immutable CAL9 checkpoint at commit
`da1750ffe4139e806e151846fc9945c07b30b449`, reconstructs both predecessor
states, constructs and source-checks the new `16385` input, and proves that
the RSP2 namespace is unused.  That namespace observation belongs to the
prelaunch commit.  After execution, verification consumes the immutable
authorization rather than rerunning the historical absence test.

The input state itself is bound by an exact array-content hash.  The source
precheck is bound as a strict pass against the frozen `1e-12` maximum, while
the raw binary64 infinity norm is left to the runtime diagnostics.  Its last
bits can vary with the platform linear-algebra reduction and therefore do not
mutate the study identity; a nonfinite, negative, or above-threshold value
still fails closed.

Raw checkpoints remain outside Git.  Their hashes and the reduced evidence
needed to audit the conclusion are tracked.

## Scientific boundary

RSP2 is neither a collapse calibration nor an FGC-QR trajectory.  It cannot
select a GR-0 case, test regulator activation, establish defocusing, reject
FGC-QR, validate retained EFT, resolve a singularity, derive a child domain,
identify a dark sector, or vary a locally measured speed of light.

## Reproduction

Before the raw namespace is created:

```bash
python3 scripts/reproduce_fgc_rsp2_frz1.py --check
python3 scripts/run_fgc_rsp2_constraint_study.py --check-authorization-only
```

The long run is launched separately after the authorization commit:

```bash
python3 scripts/run_fgc_rsp2_constraint_study.py
```
