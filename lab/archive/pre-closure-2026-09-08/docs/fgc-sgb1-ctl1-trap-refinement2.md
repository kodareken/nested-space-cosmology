# FGC-1-SGB1-CTL1-TRAP-REFINEMENT2 — successor product-box resource ladder

This owner is a separately named prospective successor. It does not edit
the depth-8/10/12 instrument. It is not `SGBL_branch_owned_and_healthy`,
not `FRZ1`/`PREF1`, and not holdout.

## Predecessor bind

The preserved first-ladder nominal contract hash is required before any
successor Picard run:

```text
predecessor = FGC-1-SGB1-CTL1-TRAP-REFINEMENT
depths      = 8, 10, 12
hash        = a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8
```

A different predecessor digest is `wrong_predecessor`. The first ladder's
caps, files, and nonpass are not rewritten here.

## Frozen successor ladder, declared before evaluation

Unchanged: 16 base cells, Picard iterations 12, `lambda∈[1/8,16]`,
`k∈[-4,4]`, `C∈[-16,16]`, `A_chi=3`, bump/equations, scientific gate
`C<1`. Resource caps are explicit practical budgets, not combinatorial
worst-case formulas and not fitted after a pass/fail.

| Level | Depth | Cell cap | RHS cap | Bit cap |
|---|---:|---:|---:|---:|
| `depth_14` | 14 | 4096 | 131072 | 65536 |
| `depth_16` | 16 | 8192 | 262144 | 65536 |
| `depth_18` | 18 | 16384 | 524288 | 65536 |

Every level runs even if one passes. A cap miss is typed `resource_limit`.
No informal depth beyond 18.

## Nominal `A_chi=3` result

All three levels finish inside the frozen caps. Inventories are gap-free
prefixes from `r=10`. None tiles `[10, 14]`. None has `C<1`. The
obstruction is `compactness_enclosure_not_below_one`. Shared overlapping
cells contract. The reported support `C` upper bounds are **not**
nonincreasing: later last cells wrap slightly above the earlier uppers.
That is recorded, not repaired by raising caps.

```text
depth 14: C_upper = 2305843524457714205 / 2^61
          margin  = -515244020253 / 2^61
          cells=28, RHS=584, C_r=146, bits=76
          coverage=[10, 372247/32768], tiling=false
depth 16: C_upper = 18446753041991488113 / 2^64
          margin  = -8968281936497 / 2^64
          cells=32, RHS=648, C_r=162, bits=76
          coverage=[10, 744495/65536], tiling=false
depth 18: C_upper = 18446754249224814419 / 2^64
          margin  = -10175515262803 / 2^64
          cells=37, RHS=736, C_r=184, bits=76
          coverage=[10, 2977981/262144], tiling=false
```

Compact payload SHA-256:

```text
nominal A_chi=3     2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d
named A_chi=1/8     9a1f2c91b2b3682b88cdad9aedf121ee5bbcd35a13d1d8bfdd0443bd8cb4df1a
```

The named `A_chi=1/8` control tiles `[10, 14]` with `C<1` at every
successor level (`C_upper = 31242091780885377 / 2^64`). Aggregate flags
remain false.

## Attacks

The contract refuses a wrong predecessor hash, changed resource caps,
reordered or skipped levels, insufficient resources presented as a
compactness nonpass, a forged pass or complete tiling, and promoting
`SGBL_branch_owned_and_healthy`/`FRZ1`/`PREF1`/holdout.

## Remaining gap

Depth 18 is still a prefix with `C>1`. Closing nominal `A_chi=3` would
need a later predeclared instrument, not a post-result cap raise and not
an informal depth past 18.
