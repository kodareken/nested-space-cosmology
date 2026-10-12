# FGC-1-SGB1-CTL1-TRAP-TAYLOR — prospective Taylor-model graph enclosure

This owner is a separate prospective instrument. It does not edit the
Picard `(lambda,k,C)` ladders, the `(C,k)` certificate chart, the matched
family, the compact bump, the physical domains, or the strict `C<1` gate.
It is not `SGBL_branch_owned_and_healthy`, not `FRZ1`/`PREF1`, and not
holdout.

## Claim

On the existing product-box graph `(lambda,k,C)`, freeze a finite Taylor
order and a subdivision/resource/bit policy, then evaluate the whole
policy. Every required `r`/state derivative, including compact-bump
derivatives, is produced by exact interval automatic differentiation.
The Lagrange remainder of Taylor's theorem is the high-order remainder;
sampled finite differences are not a proof.

A local pass may say only that one prospectively declared level proves
continuous no-initial-trap on that graph. The instrument does not promote
the product-box chart to the `(C,k)` certificate owner.

## Taylor theorem

If `f` is `C^{n+1}` on a convex cell containing the expansion point `c`
and the evaluation point `x`, then

```text
f(x) = sum_{k=0}^{n} f^{(k)}(c)/k! (x-c)^k
     + f^{(n+1)}(ξ)/(n+1)! (x-c)^{n+1}
```

for some `ξ` between `c` and `x`. The module encloses `f^{(k)}(c)` by a
univariate interval jet and encloses `f^{(n+1)}` on the whole cell. The
pair `(polynomial, remainder interval)` is a Taylor model of order `n`.
The ODE recurrence `y^{(m+1)} = f^{(m)}` is the same theorem applied
along the affine constraint vector field.

Compact-bump derivatives of `B(x)=exp(1-1/(1-x^2))` use the same jets on
cells where `1-x^2>0`. Cells that touch the support edge use proved
global bounds obtained from `B^{(n)}=Q_n B` with
`Q_n=Q_{n-1}'+Q_{n-1}Q_1` and `e^{-u} u^k <= k!`.

## Frozen policy, declared before evaluation

Unchanged: 16 base cells, `A_chi=3`, support `[10,14]`,
`lambda∈[1/8,16]`, `k∈[-4,4]`, `C∈[-16,16]`, bump/equations, scientific
gate `C<1`. Taylor order, depths, and caps are module-level constants.
They are not fitted after a pass/fail.

```text
Taylor order            = 4
base cells              = 16
Taylor-Picard iterates  = 8
lambda domain           = [1/8, 16]
k domain                = [-4, 4]
C resource domain       = [-16, 16]
scientific gate         = C < 1, unchanged
jet calls / partition   = 2 + 8
cell cap(D)             = 16 * (2^{D+1} - 1)
jet cap(D)              = cell_cap(D) * 10 * (4 + 2)
bit cap                 = 16384
```

| Level | Order | Depth | Cell cap | Jet cap | Bit cap |
|---|---:|---:|---:|---:|---:|
| `order_4_depth_2` | 4 | 2 | 112 | 6720 | 16384 |
| `order_4_depth_4` | 4 | 4 | 496 | 29760 | 16384 |

Every level runs even if one passes. A cap miss is typed `resource_limit`.
Caps are not raised after seeing an outcome.

## Predecessor hashes (regression only)

The depth-ladder Picard instruments are not re-run here and are not the
Taylor proof:

```text
TRAP-REFINEMENT   a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8
TRAP-REFINEMENT2  2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d
```

A different digest is `wrong_predecessor`.

## Nominal `A_chi=3` result

Both declared levels finish inside the frozen caps. Inventories are
gap-free prefixes from `r=10`. None tiles `[10, 14]`. None has `C<1`.
The obstruction is `compactness_enclosure_not_below_one`. Shared
overlapping cells contract. Reported support `C` upper bounds are
nonincreasing.

```text
order 4 depth 2: C_upper = 9979474004043034349 / 2^63
                 margin  = -756101967188258541 / 2^63
                 cells=11, jets=190, C_r=80, bits=105
                 coverage=[10, 183/16] = [10, 11.4375], tiling=false
order 4 depth 4: C_upper = 9307090083172832251 / 2^63
                 margin  = -83718046318056443 / 2^63
                 cells=13, jets=250, C_r=90, bits=105
                 coverage=[10, 731/64] = [10, 11.421875], tiling=false
```

Compact payload SHA-256:

```text
nominal A_chi=3     7de1421d200f9d276b95a7e541c0dd8b4b75eddddbe966db6611316b969f820b
named A_chi=1/8     d3cb9b3dd90c267e9851516025376eccce41453bdb6b853fdfcf49e630852d5a
```

The named `A_chi=1/8` control tiles `[10, 14]` with `C<1` at every
declared level (`C_upper = 28607003284473089 / 2^64`). That local pass
still leaves `SGBL_branch_owned_and_healthy`, `FRZ1`, `PREF1`, and
holdout false.

## Existence, residuals, tiling, compactness

A cell is accepted only when:

1. a first-order interval self-map of `(lambda,k,C)` is strict (existence);
2. the order-4 Taylor model remainder is the Lagrange enclosure, not an
   injected zero, and intersects that existence box;
3. affine `H` and `M` residuals contain the origin;
4. algebraic `C=1+r^2 k^2-lambda^{-2}` overlaps the propagated `C` with
   invariant residual containing zero;
5. the returned inventory is strictly ordered, nonoverlapping, and
   gap-free from `r=10` through the last returned right boundary.

A gap-free prefix is not a complete tiling. A complete tiling with every
cell and the analytic exterior `C=2M/r` strictly below 1 is the only
local no-trap bit. That bit still leaves aggregate health false.

## Named controls

Exact polynomial IVPs `y=r^d` for `d<=4` have vanishing remainder
derivative and reconstruct the endpoint exactly. The known ODE `y'=2r`
reconstructs `y=r^2`. Those controls certify the jet/remainder calculus
independently of SGB-L.

`A_chi=1/8` is the existing named positive family control. It uses the
same frozen Taylor policy, domains, bump, equations, and `C<1` gate.

## Attacks

The contract refuses:

- a remainder that does not contain the Lagrange enclosure (`wrong_remainder`);
- a jet that drops a required derivative (`omitted_derivative`);
- insufficient resources presented as a compactness nonpass (`resource_limit`);
- a changed order, cap, domain, or undeclared depth (`changed_policy`);
- reordered or skipped declared levels;
- a dropped or gapped cell inventory;
- a wrong predecessor hash;
- a forged pass, promoting `SGBL_branch_owned_and_healthy`/`FRZ1`/`PREF1`/holdout,
  or a tampered payload/hash.

Undeclared amplitudes are not accepted as informal retries.

## Aggregate flags

The compact contract payload hashes the frozen policy, predecessor
regression hashes, and the exact level outcomes. Wall-clock times are
diagnostic only and are not in the hash.

The following remain false:

- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- execution and holdout authorization

The product-box graph remains not the continuum certificate chart.
FGC-QR health evidence is not imported.
