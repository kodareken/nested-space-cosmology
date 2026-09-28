# FGC-1-SGB1-CTL1-TRAP-BARRIER — trajectory-free `D=1-C` Nagumo nonpass

This owner is a separately named prospective instrument. It does not edit
the Picard or Taylor ladders, the `(C,k)` certificate chart, the matched
family, the compact bump, or the initial-health ODE. It is not
`SGBL_branch_owned_and_healthy`, not `FRZ1`/`PREF1`, and not holdout.

## Claim

On the affine SGB-L constraint vector field, freeze the barrier

```text
D = 1 - C = lambda^{-2} - r^2 k^2
```

and ask whether `{D>0}` is a whole-domain Nagumo invariant of that field.
The instrument binds one exact outward witness on `D=0` for nominal
`A_chi=3`. That is enough to reject a trajectory-free invariant-barrier
proof. It does not claim that the actual constraint orbit traps, and it
does not reject SGB-L.

Classification: `barrier_nonpass_not_trajectory_trapping`.

## Affine substitution

The polar-areal compactness and its chain-rule derivative are the existing
identities

```text
C   = 1 + r^2 k^2 - lambda^{-2}
C_r = 2 r k^2 + 2 r^2 k k_r + 2 lambda_r / lambda^3.
```

The annular kernel remains affine and diagonal, `lambda_r=-H0/H_L` and
`k_r=-M0/M_k`, whenever `H_L` and `M_k` exclude zero. Substituting those
roots gives

```text
D_r = -C_r
    = -2 r k^2 + 2 r^2 k (M0/M_k) + 2 H0 / (H_L lambda^3).
```

A mutated formula, a dropped `H_L M_k` denominator, or a chain-rule
mismatch is refused.

## Barrier parameterization, both signs

`D=0` is the surface

```text
k = sigma / (r lambda),    sigma in {+1, -1}.
```

Both signs are exact rational identities at the frozen witness radius and
lambda. The outward witness uses `sigma=+1`. The opposite sign is recorded;
it is not an outward witness at this point.

## Null-expansion equivalence

On the unit-lapse, zero-shift, `R=r` slice

```text
theta_+ = 2 (-k + 1/(r lambda))
theta_- = 2 (-k - 1/(r lambda))
theta_+ theta_- = -(4/r^2) D.
```

Vanishing `D` is therefore equivalent to a vanishing future null expansion
of the 2-sphere. At the outward witness, `theta_+=0` and
`theta_-=-80/209`.

## Frozen outward witness

Declared before evaluation, not fitted after a pass/fail:

```text
r      = 209/20
lambda = 1
k      = 20/209
A_chi  = 3
```

Exact rational identities: `D=0`, `C=1`. Interval affine evaluation:
`H_L` and `M_k` exclude zero, and

```text
D_r < -3/8.
```

The vector field on the barrier points toward `D<0` (`C>1`). A whole-domain
Nagumo proof of `{C<1}` therefore fails. The actual radial orbit from the
Minkowski buffer `(lambda,k,C)=(1,0,0)` is not computed here and is not
claimed to hit this point.

## Open neighborhood

On the open radius interval `(52/5, 21/2)`, with the same `lambda=1` and
the `+` barrier parameterization, the closed interval enclosure satisfies

```text
D_r < -1/20.
```

The witness `r=209/20` is strictly inside. Both Jacobian diagonals still
exclude zero. The bound on the closed interval implies the open interval.

## Named controls

Exact vacuum, vanishing scalars, both `k` signs on `D=0`:

```text
D_r = 1/r = 20/209.
```

Vacuum therefore points inward on the barrier. The obstruction is the
nominal family, not the substitution algebra.

Named `A_chi=1/8` at the same geometric witness has `D_r>0`. That local
inward control uses the same formula, denominator, and `C<1` barrier. It
is not a retry that changes the nominal amplitude after seeing `A_chi=3`,
and it still leaves every aggregate flag false.

## Scalar-mass majorant nonclosure

The declared `phi` amplitude is `1/131072`. Dropping `phi` entirely, or
replacing the frozen scalar mass `mu=3` by the a priori majorant `mu=6`,
leaves `D_r < -3/8`. No scalar-mass majorant restores an inward barrier
field at the witness, because the obstruction is the `chi` kinetic stress
of nominal `A_chi=3`.

## Attacks

The contract refuses:

- a mutated `D_r` formula (`formula_mutation`);
- the opposite `k` sign, a boolean sign, or a `k` that is not the `D=0`
  parameterization (`sign_mutation`);
- a dropped or replaced `H_L M_k` denominator (`denominator_mutation`);
- a different radius, lambda, or witness point (`witness_mutation`);
- `A_chi=1/8` or any undeclared amplitude presented as the obstruction
  (`family_mutation`);
- promoting `SGBL_branch_owned_and_healthy`, `FRZ1`, `PREF1`, holdout,
  execution, trajectory trapping, or SGB-L rejection (`claim_mutation`);
- a tampered payload or hash (`tampered_contract`).

## Aggregate flags

Compact payload SHA-256:

```text
a75d564c0ea0658802c138e0575edf386375ebba35728661bd8a11c579e7d689
```

The payload hashes the frozen witness, both `k` signs, the neighborhood
enclosure, vacuum `1/r`, the `A_chi=1/8` control, and the scalar-mass
majorant nonclosure. The following remain false:

- `SGBL_branch_owned_and_healthy`
- `FRZ1`, `PREF1`
- execution and holdout authorization
- `whole_domain_nagumo_invariant`
- `trajectory_trapping_claimed`, `actual_orbit_traps`
- `sgbl_model_rejected`

Sampled constraint-integrator nodes are not the certificate. FGC-QR health
evidence is not imported.
