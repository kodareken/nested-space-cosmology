# FGC-1-SGB1-CTL1 — continuous initial health, centre monitor, and controls (prospective)

This owner is an implementation map for a validated continuous
initial-slice compactness proof, a fail-closed evolving-centre monitor
contract, and the required branch-owned algebraic controls. It is not
`SGBL_branch_owned_and_healthy` and not a holdout pass.

## Continuous compactness

Compactness is the polar-areal identity

```text
C = 2m/R = 1 + R^2 k^2 - lambda^{-2}.
```

Sampled RK4/SSPRK3 constraint-integrator nodes are not a continuum
certificate. The family property `continuous_no_initial_trapped_sphere`
remains unqualified there.

### Validated constraint ODE

The certificate state is the Misner–Sharp chart `(C,k)`:

```text
C = 1 + r^2 k^2 - lambda^{-2}
D = 1 + r^2 k^2 - C = lambda^{-2}
lambda = D^{-1/2}   (positive square-root branch, D>0)
```

`k` is a chart coordinate, not reconstructed from `J=r^3 k`. The affine
pair is unchanged, `lambda_r=-H0/H_L`, `k_r=-M0/M_k`. The state ODEs
are `k_r` and

```text
C_r = 2 r k^2 + 2 r^2 k k_r + 2 lambda_r / lambda^3.
```

At singletons the `(C,k)` roundtrip `(lambda,k)->(C,k)->lambda` has
residual `{0}` when `D` is a perfect square, and otherwise the residual
enclosure contains zero. A lost or sign-changing `D`, or the negative
square-root branch, is fail-closed. Initial data at the Minkowski-buffer
boundary are exact: `(lambda,k,C)=(1,0,0)`, `D=1`.

The declared interval-Picard integrator owns the continuum graph:

1. Picard self-maps in `(C,k)`. Reconstructed `lambda` is used for the
   affine H/M residual and for `H_L`, `M_k` denominators. Physical
   domains remain `lambda in [1/8, 16]`, `k in [-4, 4]`. Compactness
   chart `C in [-16, 16]` is a resource domain containing a
   neighbourhood of `0`; it is not the scientific `C<1` gate. There is
   no `J` Picard state.
2. The same declared tree: 16 base cells, max bisection depth 8, max
   1024 cells, max 8192 affine RHS evaluations, max 12 Picard
   iterations, bit cap 16384, Euler inflation 2, outward `2^{-64}`
   dyadic rounding. The compact bump rule is unchanged. None of these
   is fitted after seeing pass/fail.
3. On each accepted cell the Picard image
   `(C,k)(r0)+[0,h] f(R,Y)` lies strictly inside the `(C,k)` graph box,
   `D` is strictly positive, both Jacobian diagonals exclude zero, and
   affine residuals contain the origin.
4. Endpoints propagate in `(C,k)`. `lambda` at an endpoint is
   reconstructed from `(C,k)` there.
5. The product-box `(lambda,k,C)` Picard and the `(C,J)` chart Picard
   are regression only. They are not the certificate.
6. Bisection concatenates left and right children. The returned
   inventory is strictly ordered, nonoverlapping, and gap-free from the
   requested left boundary through the last returned right boundary. A
   complete record tiles `[r_support min, r_support max]` exactly. A
   failed covering keeps every proved prefix cell and must not claim
   completeness.

There is no `GB_EFFECTIVE_DENSITY_BOUND` and no integrated energy/mass
proxy in the nominal proof.

### Exterior

At the enclosed support-boundary state, Misner–Sharp mass satisfies
`M = r C / 2`. The analytic vacuum continuation has `C=2M/r`, which is
monotone for `M` of one sign, so its maximum on `[r_support, r_outer]`
is the boundary compactness (or zero if `C_end<=0`). The proof requires
that maximum strictly below 1.

### Minkowski buffer

`r` from the regular centre through the compact-support minimum has
`lambda=1`, `k=0`, so `C=0` and `D=1`.

### Nominal `chi=3` result

Recomputed on the same declared tree after fixing partition lineage so
left siblings are retained. The local continuum fact is
`interval_inconclusive` / `picard_strict_self_map_failed`. The covering
is a gap-free prefix of the compact support, not a complete tiling, so
it is not a no-trap certificate.

```text
(C,k) certificate
  obstruction = picard_strict_self_map_failed
  C_upper = 17053332284082708107 / 2^64  ≈ 0.924463
  margin  = 1393411789626843509 / 2^64   ≈ 0.075537
  cells   = 9 (was 1 before lineage fix)
  affine RHS = 492 evaluations
  C_r evals  = 123
  coverage = [10, 10749/1024] ≈ [10, 10.4970703125]
  last cell = [2687/256, 10749/1024] at depth 8
  D_margin  = 326806051576253549 / 2^62  ≈ 0.070865
```

Comparison with prior charts after the same lineage fix:

```text
(lambda,k,C) cells = 18, compactness_enclosure_not_below_one
  C_upper = 18453842205463973737 / 2^64 ≈ 1.000385
  coverage = [10, 727/64] ≈ [10, 11.359375]
  last cell = [11631/1024, 727/64] at depth 8
  affine RHS = 384
(C,J) cells = 12, picard_strict_self_map_failed
  C_upper = 17069501262171800641 / 2^64 ≈ 0.925340
  coverage = [10, 10757/1024] ≈ [10, 10.5048828125]
  last cell = [2689/256, 10757/1024] at depth 8
  affine RHS = 640
```

All three inventories abut and start at `r=10`. None tiles through
`r=14`. Affine evaluation counts are unchanged from the dropped-sibling
bug; only the kept prefix grew. The tree, `A_chi=3`, domains, and `C<1`
threshold were not changed. A smaller amplitude `chi=1/8` still closes
`C<1` on all 16 base cells on every chart, tiling `[10, 14]` exactly
(`(C,k)` `C_upper = 31314122837433715 / 2^64 ≈ 1.698e-3`;
`(lambda,k,C)` `31242091780885377 / 2^64 ≈ 1.694e-3`;
`(C,J)` `31411545268290931 / 2^64 ≈ 1.703e-3`).

Lost `D`, a non-positive square-root argument, and the negative
square-root branch are typed failures. A flipped `lambda_r/lambda^3`
sign is not `C_r`. Forging a pass from the product-box or `(C,J)` graph
is refused.

### When it does not pass

If the declared tree cannot produce a strict Picard self-map, a
zero-free diagonal, residual containment, overlapping algebraic and
propagated `C` enclosures, or correlated `C<1` on the graph, the record
is typed `interval_inconclusive` with the precise obstruction
(`picard_strict_self_map_failed`, `jacobian_diagonal_contains_zero`,
`physical_domain_escape`, `compactness_enclosure_not_below_one`,
`residual_enclosure_misses_origin`, `exterior_compactness_not_below_one`,
`zero_in_interval_reciprocal`, `mean_value_inconsistent`,
`compactness_mean_value_inconsistent`,
`algebraic_compactness_mean_value_inconsistent`,
`compactness_invariant_misses_origin`,
`chart_denominator_not_strictly_positive`, `sqrt_branch_not_positive`, or
`chart_roundtrip_misses_origin`) and the matching missing theorem
name. Sampled nodes cannot fill the gap. Resource and domain failures
are typed stops, not compactness nonpasses. A compactness chart with
upper bound `<=1` or that does not contain a neighbourhood of `C=0` is
a domain stop: the chart must not assume the scientific gate.

RK4/SSPRK3 containment is a regression check only.

## Evolving-centre monitor

The monitor is a contract for later trajectories, grounded in:

- `sgb1_ctl1_center.sgbl_initial_center_series`
- `validate_sgbl_initial_center_profile`
- `regular_center.LaurentSeries` / `SeriesJet2`

Required later checks: even parity, elementary flatness `A=lambda` in
value/dt/dtt, certified negative Laurent powers absent, and emptiness
not assumed after matter arrives.

This slice has no trajectory owner. The contract is therefore
`monitor_not_executed`. An unexecuted monitor is not a pre-holdout pass.
`sgbl_execute_evolving_center_monitor` raises `monitor_not_executed` or
`missing_trajectory_owner`. Forging `executed=True` is refused.

## Branch-owned controls

All four controls remain local algebraic identities in
`sgb1_ctl1_controls.py`. Aggregate health stays false. None is a full
backreacted SGB-L solution. The established control is the regular
decoupling-limit linear-GB scalar on fixed Schwarzschild, not a
backreacted solution.

## Aggregate flags

`continuous_no_initial_trapped_sphere` may be true as a local continuum
fact on a slice that closes. The monitor never contributes a pass here.
`SGBL_branch_owned_and_healthy`, `FRZ1`, `PREF1`, execution, and holdout
remain false. No FGC-QR or GR-0 health certificate is imported.
