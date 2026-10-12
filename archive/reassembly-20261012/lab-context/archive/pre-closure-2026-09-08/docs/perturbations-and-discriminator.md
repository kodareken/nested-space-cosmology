# Closed de Sitter core: perturbations and discriminator

This note defines a deliberately narrow calculation for the selected **closed
de Sitter core**.  It is a regular, homogeneous benchmark for part of a
transition calculation.  It is not yet a covariant parent-black-hole-to-child
solution, and does not establish that a physical transition takes this form.

## Background and domain

In dimensionless conformal time

```text
-pi/2 < eta < pi/2,
```

use

```text
a(eta) = L sec(eta),
a''/a = 2 sec^2(eta) - 1.
```

`L` is the curvature radius.  The code evolves a finite symmetric interval
`[-eta_core, +eta_core]`, strictly inside this coordinate patch.  It does not
interpret either endpoint as an asymptotic in/out vacuum, a horizon, a junction
surface, or a reheating hypersurface.

## Gauge-invariant variables in this benchmark

For a spectator/order-parameter field with constant homogeneous background
`Φ' = 0`, the field fluctuation is gauge invariant at linear order because its
gauge transformation is proportional to the background derivative.  Expanding
in scalar harmonics on the unit three-sphere, with `n >= 0`, the canonical mode
obeys

```text
u_n'' + [(n + 1)^2 + ((m L)^2 - 2) sec^2(eta)] u_n = 0.
```

Here `m L` is dimensionless.  This is a spectator/order-parameter calculation.
It is **not** the Mukhanov--Sasaki curvature perturbation: when `Phi' = 0`, the
usual single-clock curvature variable is not supplied by this field and its
standard `z=a Phi'/H` normalization is singular.  A completed model must derive
its actual constrained scalar sector from its covariant action.

For tensor harmonics, using the project convention `n >= 3`, the canonical
mode is

```text
mu_n'' + [n^2 - 2 sec^2(eta)] mu_n = 0.
```

The convention is stated explicitly because closed-universe harmonic labels
are not universal across the literature.  Bonga, Gupt, and Yokomizo derive the
same `S^3` tensor eigenvalue and canonical tensor equation in this convention;
their scalar label is shifted by one relative to the `n >= 0` scalar label used
here [arXiv:1612.07281](https://arxiv.org/abs/1612.07281).

Positive canonical kinetic terms and unit principal propagation speeds exclude
ghost and gradient instabilities in this benchmark.  They do not imply that
every infrared solution is bounded.  In particular,

```text
omega_n^2 = (n + 1)^2 + ((m L)^2 - 2) sec^2(eta)
```

can be negative for sufficiently light, low-`n` scalar modes.  That is a
tachyonic/infrared growth question, not a ghost or gradient instability.  A
bounded-evolution claim would require a declared mass, harmonic, core interval,
and amplitude criterion; the present code reports transfer and convergence
only for its explicit benchmark inputs.

Thus this calculation is also the recovered common-cone baseline for
[DSF-1](domain-settings.md), not a variable-`c` result. A future domain-setting
extension must derive a modified principal symbol from its covariant action
and report an invariant relative characteristic speed such as
`c_mode/c_gamma-1`; changing coordinates or rescaling the metric conformally
would not supply such a signal.

## What the code calculates

[`perturbations.py`](../src/recursive_horizons/perturbations.py) integrates
each real second-order equation as two independent phase-space solutions.  The
result is a real fundamental transfer matrix

```text
(y, y')_out = M (y, y')_in,
M = [[M11, M12], [M21, M22]].
```

For a real linear oscillator without dissipation, exact evolution preserves the
Wronskian:

```text
det(M) = 1.
```

The implementation uses deterministic classical RK4 and reports both the
determinant error and changes under step doubling.  The unit test with `mL =
sqrt(2)` is a controlled exact case: the scalar equation becomes the constant
frequency oscillator `u''+(n+1)^2 u=0`, whose transfer matrix is known
analytically.

This transfer is **not** yet any of the following:

- a Mukhanov--Sasaki curvature spectrum;
- a specified adiabatic quantum state or Bogoliubov particle-production result;
- a prediction of `A_s`, `n_s`, running, tensors-to-scalars, non-Gaussianity,
  BBN, or CMB bandpowers;
- an information-preserving black-hole transition map;
- evidence that our universe emerged from a black-hole interior.

Those require a completed covariant action, constraint analysis, endpoint and
state prescription, coupling to the hot Big Bang, and an observational
likelihood.  A determinant close to one only establishes the expected
classical symplectic property of this specified linear ODE.

## Single independent discriminator: closed spatial curvature

The sign test requires an additional global postulate: the observable child
must remain on the strict `S^3`, `K=+1` branch after the controlled core.  CCT-1
alone is a local/global-de-Sitter core laboratory and does not derive that
post-core matching.  Conditional on the added branch postulate, the model
predicts the **sign**

```text
Omega_K < 0.
```

With the standard FLRW convention

```text
Omega_K = -K c^2 / (a_0^2 H_0^2),
```

positive closed curvature `K>0` has negative `Omega_K`.  This sign statement is
independent of H-SAT: it does not use the parent entropy relation, an inferred
parent mass, or the observed target value of `Lambda`.  It can be tested using
geometric cosmological curvature constraints rather than fitting an H-SAT
parameter back to the same target.

It is deliberately modest.  Negative `Omega_K` (positive spatial curvature) is not unique to a
black-hole-to-child model, so confirming it would not establish the proposed
genealogy.  Conversely, a robust positive-`Omega_K` inference that survives a
stated family of nonflat cosmological models and systematic checks would rule
out the strict closed branch.  Finally, an inflationary or otherwise prolonged
post-core expansion can make the magnitude `|Omega_K|` arbitrarily small.  An
observational result statistically consistent with zero therefore does not
confirm or falsify the branch unless the full transition fixes the amount of
dilution and predicts a nonzero lower bound.

## Reproducibility contract

Run the focused calculation with:

```bash
python -m unittest tests.test_perturbations
```

Before treating a numerical transfer as stable, require at least:

1. agreement under step refinement;
2. determinant/Wronskian error near numerical precision;
3. agreement with the constant-frequency `mL=sqrt(2)` regression case;
4. all endpoints remain strictly inside `(-pi/2, pi/2)`.

These are necessary numerical checks, not sufficient physics checks.
