# FGC-1-SGB1-CTL1 — prospective product-box continuous-no-trap resource ladder

This owner is a separate prospective instrument. It does not change the
`(C,k)` certificate chart, the matched family, the compact bump, the
`lambda`/`k`/`C` physical domains, or the strict `C<1` gate. It is not
`SGBL_branch_owned_and_healthy`, not `FRZ1`/`PREF1`, and not holdout.

## Claim

On the existing product-box graph `(lambda,k,C)`, freeze a finite
interval-Picard *resource* ladder and then evaluate
`sgbl_validated_lambda_k_constraint_ode` at every declared level.

A local pass may say only that one prospectively declared level proves
continuous no-initial-trap on that product-box graph. Sampled RK4/SSPRK3
nodes are not the proof. The instrument does not promote the product-box
chart to the `(C,k)` certificate owner.

## Frozen ladder, declared before evaluation

The ladder is a module-level constant. Caps are the combinatorial worst
case of the declared 16-base-cell binary tree, not a depth chosen after a
pass/fail, and not a fitted amplitude or `C` threshold.

```text
base cells              = 16
Picard iterations       = 12
lambda domain           = [1/8, 16]
k domain                = [-4, 4]
C resource domain       = [-16, 16]   (contains a neighbourhood of 0; upper > 1)
scientific gate         = C < 1, unchanged
affine evals / RHS call = 4
RHS calls / partition   = 1 + 12
cell cap(D)             = 16 * (2^{D+1} - 1)
RHS cap(D)              = cell_cap(D) * 13 * 4
bit cap(D)              = 16384 * 2^{(D-8)/2}
```

| Level | Depth | Cell cap | RHS cap | Bit cap |
|---|---:|---:|---:|---:|
| `depth_8` | 8 | 8176 | 425152 | 16384 |
| `depth_10` | 10 | 32752 | 1703104 | 32768 |
| `depth_12` | 12 | 131056 | 6814912 | 65536 |

Depth 8 keeps the existing owner's bit cap. Later levels only enlarge
predeclared cell/RHS/bit budgets. The first level's cell/RHS caps are
larger than the certificate owner's `1024/8192` so a full depth-8 tree can
finish; that is a resource declaration, not a change of the existing chart
owner.

The evaluator runs every level even if one passes. If a level exceeds its
own cap, the result is typed `resource_limit` and later independently
declared levels still run. Caps are not raised after seeing an outcome.

## Nominal `A_chi=3` result

Recomputed on the corrected gap-free partition lineage. The three declared
levels all finish inside their frozen caps. Each inventory is a strictly
ordered prefix from `r=10` through the last returned cell; none tiles
`[10, 14]`, so none is a continuous-no-trap proof. The obstruction at every
level is `compactness_enclosure_not_below_one`. Shared overlapping cells
contract, and the reported `C` upper bounds are nonincreasing. Affine
evaluation counts are unchanged from the dropped-sibling bug; only the kept
prefix grew. Aggregate health remains false.

```text
depth 8:  C_upper = 18453842205463973737 / 2^64
          margin  = -7098131754422121 / 2^64
          cells=18, RHS=384, C_r=96, bits=76
          coverage=[10, 727/64] ≈ [10, 11.359375], tiling=false
depth 10: C_upper = 9223921512318141519 / 2^63
          margin  = -549475463365711 / 2^63
          cells=21, RHS=464, C_r=116, bits=76
          coverage=[10, 23265/2048] ≈ [10, 11.35986328125], tiling=false
depth 12: C_upper = 4611781872439090751 / 2^62
          margin  = -95854011702847 / 2^62
          cells=25, RHS=536, C_r=134, bits=76
          coverage=[10, 186123/16384] ≈ [10, 11.36004638671875], tiling=false
```

The depth-8 `C` upper bound remains the existing product-box regression
`18453842205463973737/2^64`. The scientific gate was not moved.

Compact payload SHA-256:

```text
nominal A_chi=3     a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8
named A_chi=1/8     c3cca0ef56f4f1e257bd0d2b3e0ec088505aaae2c4e782950d4c632102afd299
```

The named `A_chi=1/8` control closes `C<1` on all 16 base cells at every
declared level, tiling `[10, 14]` exactly
(`C_upper = 31242091780885377 / 2^64`). That local pass still leaves
`SGBL_branch_owned_and_healthy`, `FRZ1`, `PREF1`, and holdout false.

## Coverage, obstruction, and exact margin

Each level records:

- whether the returned Picard cells are a gap-free cover of the compact support;
- the ODE obstruction, or a typed resource/domain stop;
- exact rational `C` upper bound and `C<1` margin on produced cells;
- exterior `C=2M/r` only when the ODE finished a complete cover;
- cells, affine RHS evaluations, compactness RHS evaluations, and observed bits;
- whether overlapping enclosures on shared coverage contract monotonically.

Coverage is the validated inventory of cells returned by
`sgbl_validated_lambda_k_constraint_ode`, checked with
`sgbl_validate_ode_cell_inventory`. A dropped left sibling, overlap, or
radial gap is refused. A gap-free prefix is not a complete tiling. A `C`
upper bound on a prefix is not a no-trap certificate.

A level proves continuous no-initial-trap only if the inventory tiles the
whole compact support, every cell is a strict Picard self-map with residual
and invariant containment, support and exterior `C` are strictly below 1,
and shared coverage with the other declared levels is monotone. That local
bit still leaves aggregate health false.

## Named small-amplitude control

`A_chi=1/8` is the existing named positive control. It uses the same frozen
ladder, domains, bump, equations, and `C<1` gate. It is not a retry that
changes the nominal amplitude after seeing `A_chi=3`.

## Attacks

The contract refuses:

- insufficient resources presented as a compactness nonpass;
- a changed ladder (caps, domains, or undeclared depth);
- reordered or skipped declared levels;
- a nonmonotone overlapping enclosure treated as contraction;
- a dropped or gapped cell inventory, or a forged complete tiling on a prefix;
- a forged pass, promoting `SGBL_branch_owned_and_healthy`/`FRZ1`/`PREF1`/holdout,
  or a tampered payload/hash.

Undeclared amplitudes are not accepted as informal retries.

## Aggregate flags

The compact contract payload hashes the frozen ladder and the exact level
outcomes. Wall-clock times are diagnostic only and are not in the hash.

The following remain false:

- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- execution and holdout authorization

The product-box graph remains not the continuum certificate chart.
