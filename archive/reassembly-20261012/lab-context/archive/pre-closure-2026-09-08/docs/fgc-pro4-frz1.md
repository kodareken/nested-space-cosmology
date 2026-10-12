# FGC-1-PRO4-FRZ1: outcome-neutral PROTO4 freeze

## Decision

`FGC-1-PRO4-FRZ1` freezes **FGC-2-SF1-PROTO4** as the premise-only successor
to the numerically obstructed PROTO3 contract. It proves that the successor is
cryptographically tied to immutable checkpoint
`b5aa4934b58afd505b22e008202285c475942892`, and that the checkpoint contains
the exact PROTO3 protocol and CAL0 diagnosis named by the amendment.

The machine decision is:

```text
PROTO4_outcome_neutral_protocol_frozen = true
PROTO4_immutable_lineage_verified = true
PROTO4_branch_specific_stop_partition_verified = true
PROTO4_convergent_admission_contract_frozen = true
PROTO4_resolved_holdout_manifest_authorized = false
FGCQR_holdout_execution_authorized = false
```

This is a repaired test contract, not a mechanism result.

## What changed

PROTO4 replaces only the three numerical-contract defects established by
CAL0:

1. all 26 historical HLT1 stops are partitioned exactly once into nine
   universal runtime stops, eleven candidate-branch-only stops requiring their
   own SGB-L or FGC-QR definitions, and six stops replaced by convergence-aware
   admission rules;
2. spectral admission uses top-band field power, derivative-weighted power,
   nested tail convergence, and RMS proper scales, while the last occupied
   Fourier bin is retained only as a diagnostic; and
3. physical, gauge, and reduction constraints are judged at common events on
   three nested grids, with a declared normalization and Richardson error that
   enters the observable error budget conservatively.

The calibration and holdout output namespaces are versioned under `proto4`.
Their emptiness is a future manifest precondition, not a live filesystem fact
on which this permanent certificate depends.

## What did not change

Every unlisted field is inherited by the exact SHA-256 digest of PROTO3. The
action, couplings, scalar profiles, amplitude order, held-out cases, physical
equations, affine metric-null observable, outcome labels, robustness box, and
negative-result burden are unchanged. The reproducer reads the predecessor
from the historical Git object, validates it as PROTO3, and serializes semantic
hashes for each inherited physical section.

The earlier exploratory GR-0 outcomes are disclosed. They are development
observations and are forbidden from entering the calibration certificate.
Fresh primary and comparator runs are required. No FGC-QR or SGB-L evolution
outcome was opened to design PROTO4.

## Why this does not open the holdout

PROTO4 names **FGC-1-HLT2-MON2** as the owner of the executable spectral and
constraint formulas, but HLT2 does not yet exist. The all-case static ledger,
fresh nested GR-0 calibration, selected-case comparator result, and resolved
`FGC-1-PRO4-HLD1` manifest also do not yet exist. RUN1 therefore remains
historically fail-closed under PROTO3, and no successor authorization is
inferred from this freeze.

## Reproduction and nonclaims

Run:

```bash
python3 scripts/reproduce_fgc_pro4_frz1.py \
  --output results/fgc-1-pro4-frz1.json
python3 -m unittest tests/test_fgc_protocol_v4.py -v
```

The certificate does not read either new output namespace. It derives no
collapse, trapped interval, activation, defocusing, transition, singularity
resolution, child domain, dark sector, varying locally measured light speed,
or rejection of the FGC-QR action or the general gradient route.
