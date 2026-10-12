# FGC-1-SGB1-CTL1-SOL1 — orbit-local Picard-Lindelof continuation (prospective)

This owner is a separately named prospective instrument. It does not edit
the Picard or Taylor ladders, the `(C,k)` certificate chart, the matched
family, cell admission, principal feeder, continuity, or trap
refinement/barrier owners. It is not `SGBL_branch_owned_and_healthy`, not
`FRZ1`/`PREF1`, and not holdout.

## Claim

Start from the authenticated Misner-Sharp `(C,k)` prefix of the affine
SGB-L constraint IVP on the declared `A_chi=3` family, ending at

```text
r0 = 10749/1024.
```

On the remaining declared support `[r0, 14]`, freeze a tube policy and ask
whether interval Jacobian Lipschitz bounds prove a unique local affine
`(lambda,k)` orbit on each tube. Concatenate only exact neighboring
endpoints. Evaluate compactness on the certified orbit image, not the
whole product box `lambda∈[1/8,16] × k∈[-4,4]`.

The instrument proves or refuses. It does not retune `rho`, step count,
amplitude, thresholds, or charts after the outcome.

Classification of the frozen nominal policy:

```text
interval_inconclusive / prefix_endpoint_not_inside_declared_tube
```

That is wrapping of the authenticated prefix endpoint, not an on-orbit
compactness nonpass and not a Jacobian nonpass. It is evidence that a
later matched-data redesign is required if a unique remaining-support
orbit is to be certified. It is not a trapped-orbit claim and does not
reject SGB-L.

## Affine IVP

The polar-areal constraint pair remains the inherited affine diagonal
system

```text
lambda_r = -H0 / H_L
k_r     = -M0 / M_k
```

whenever `H_L` and `M_k` exclude zero. Compactness is the inherited
identity

```text
C = 1 + r^2 k^2 - lambda^{-2}
```

evaluated on the Picard image of each certified tube. Sampled RK4/SSPRK3
nodes are regression controls, never the certificate. This is not a
spacetime or quasilinear IBVP theorem.

## Frozen tube policy, declared before evaluation

Unchanged family and support: `A_chi=3`, compact support `[10,14]`,
`lambda∈[1/8,16]`, `k∈[-4,4]`, `C∈[-16,16]`, scientific gate `C<1`,
inherited bump and H/M equations. Tube constants are module-level. They
are not fitted after a pass/fail.

```text
continuation steps      = 16
remaining support       = [10749/1024, 14]
step width              = 3587/16384
tube radius rho         = 1/8
bisection depth         = 0   (no adaptive refinement)
rational-bit cap        = 16384
infinity-norm Lipschitz = ||Df/d(lambda,k)||_inf
Picard operator         = P(y)=y0+[0,h] f(T)
```

Each tube is the closed infinity-norm ball of radius `rho` about the
midpoint of the current endpoint enclosure. The current endpoint must lie
inside that tube. The Picard image `y0+[0,h]f(T)` must lie strictly
inside the same tube. Uniqueness is the bounded interval Jacobian on `T`
with `H_L` and `M_k` excluding zero. A failed tube is not bisected.

## Predecessor identities

Bound before evaluation, not re-run as this certificate:

```text
(C,k) prefix
  obstruction = picard_strict_self_map_failed
  cells       = 9
  coverage    = [10, 10749/1024]
  last cell   = [2687/256, 10749/1024]
  C_upper     = 17053332284082708107 / 2^64
  margin      = 1393411789626843509 / 2^64
  D_margin    = 326806051576253549 / 2^62
  affine RHS  = 492
  C_r evals   = 123

TRAP-REFINEMENT   a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8
TRAP-REFINEMENT2  2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d
TRAP-TAYLOR       7de1421d200f9d276b95a7e541c0dd8b4b75eddddbe966db6611316b969f820b
TRAP-BARRIER      a75d564c0ea0658802c138e0575edf386375ebba35728661bd8a11c579e7d689
```

A different digest or prefix identity is `wrong_predecessor`. The
whole-domain Nagumo barrier nonpass is historical evidence only. SOL1
tests unique on-orbit continuation.

## Picard-Lindelof predicates

On a tube `T` of radius `rho` about the current `(lambda,k)` seed:

1. The initial enclosure lies inside `T`.
2. Affine `H_L` and `M_k` exclude zero on `T`.
3. The interval Jacobian is bounded; `L=||Df||_inf` is the Lipschitz
   constant.
4. The Picard image lies strictly inside `T`.
5. Affine `H` and `M` residuals contain the origin.
6. Algebraic `C` on the image overlaps any propagated `C` with invariant
   residual containing zero.

A unique local affine orbit is (1)--(6). Compactness `C<1` is then
evaluated on that image. If (1)--(6) hold and `C>=1`, the outcome is
on-orbit compactness nonpass. If `H_L` or `M_k` contains zero, the
outcome is on-orbit Jacobian nonpass. Wrapping, failed self-map, or an
endpoint that does not fit in `T` is interval-inconclusive. Domain and
bit-cap failures are typed stops, not compactness nonpasses.

Constructor replacements of the frozen policy, predecessor hashes, or
aggregate claims fail closed.

## Nominal `A_chi=3` result

The authenticated prefix endpoint at `r0` has `lambda` width about `2.69`,
which exceeds `2 rho = 1/4`. It is not contained in the declared tube.
The frozen policy therefore refuses unique continuation on the first
step. Caps are not raised. `rho` is not reduced. Amplitude is not
changed.

```text
classification = interval_inconclusive
obstruction    = prefix_endpoint_not_inside_declared_tube
unique tubes   = 0
coverage       = [10749/1024, 175571/16384]  (refusing first step)
tiling         = false
```

Compact payload SHA-256:

```text
nominal A_chi=3     5fc44f7023eb82632309c6c15120f2e62e5e6483911e18589579576b36654f4b
named A_chi=1/8     661f06f9e6d6fa91d64f1767cf8b8b810c6ce17a4e9bd7ccf8a293b72ed946b9
```

This obstruction should redirect to a later matched-data redesign if a
tighter authenticated endpoint is required. SOL1 does not perform that
redesign.

## Named controls

Exact `y'=0` has a unique constant orbit in the frozen tube, independently
of SGB-L wrapping. The inherited Minkowski-buffer reduction
`(lambda_r,k_r,C)=(0,0,0)` is the same control on the exact annular kernel.

`A_chi=1/8` uses the same remaining-support grid, `rho`, bit cap, and
depth 0. Its state at `r0` fits in the declared tube. The first tube is a
unique local affine orbit with on-orbit `C<1`. The second tube fails a
strict Picard self-map. That local control never replaces the nominal
family and leaves every aggregate flag false.

## Typed outcomes

| Outcome | Meaning |
|---|---|
| `unique_local_affine_orbit` | Unique affine orbit on the remaining support, `C<1` along it. |
| `on_orbit_compactness_not_below_one` | Unique orbit; on-orbit `C` is not strictly below 1. |
| `jacobian_diagonal_contains_zero` | `H_L` or `M_k` contains zero on a tube. |
| `interval_inconclusive` | Wrapping, failed self-map, or prefix endpoint outside the tube. |
| `domain_error` | Tube left the declared physical chart. |
| `resource_limit` | Rational-bit cap. |

## Attacks

The contract refuses:

- a mutated `rho`, step count, bit cap, depth, or remaining-support grid
  (`changed_policy`);
- a wrong prefix identity or trap-instrument hash (`wrong_predecessor`);
- a dropped or gapped tube inventory, or concatenated endpoints that are
  not exact (`gapped_inventory`);
- injected Jacobian-zero presented as a unique orbit;
- compactness crossing `C>=1` presented as a unique-orbit pass;
- wrapping presented as a compactness or Jacobian nonpass;
- insufficient bits presented as a compactness nonpass (`resource_limit`);
- `A_chi=1/8` or any undeclared amplitude presented as the nominal family
  (`family_mutation`);
- promoting `SGBL_branch_owned_and_healthy`, `FRZ1`, `PREF1`, holdout,
  execution, trapping, model rejection, IBVP existence, or a retuned
  tube (`claim_mutation`);
- a tampered payload or hash (`tampered_contract`).

## Aggregate flags

The following remain false:

- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- execution and holdout authorization
- `whole_domain_nagumo_invariant`
- `trajectory_trapping_claimed`, `actual_orbit_traps`
- `sgbl_model_rejected`
- `replaced_nominal_family_by_small_amplitude`
- `spacetime_ibvp_theorem`, `quasilinear_existence`

FGC-QR health evidence, campaign runners, PRO20, and shared pass bits are
not imported. Final FRZ1/PREF1 is a later independent binder.
