# Static periodic leading-sector virial

The [read-only assessment](../scripts/assess_nsc_discovery_leading_virial.py)
authenticates the actual `e6e784b` leading evaluator and preparation inputs.
It derives the continuum boundary sign and measures the exact finite SBP
identity on an explicitly declared NF32 positive-spectrum source. It does
not solve initial constraints, evolve a trajectory, or construct a boundary
pressure model. Existing leading and exact auxiliary histories remain intact.

## Exact periodic identity

Write $g=8\pi A$, $F=-gr^2/2$, $Z=-3g$,
$V=gr^2-\mathfrak m$, with $\mathfrak m=2\pi C_F\mathrm{flux}^2$.
For static geometry in $L=Q$, beta zero, positive r and Q, the invertible
leading velocities force both canonical momenta to vanish.

The massless four-dimensional Dirac source is independent of r at fixed
canonical half-density and Q. Its reduced angular gap $\kappa/r$ is not a
fundamental four-dimensional mass. Maxwell is trace-free as well. The
Einstein bulk Hamiltonian is quadratic in r under these conditions.

The actual sampled derivative D and its adjoint therefore give

$$
\boxed{-\langle Q\rho\rangle-\mathfrak m\langle Q^2\rangle
=-\tfrac12\langle r\dot p_r\rangle-\langle QC\rangle.}
$$

The fine pairing uses $\Delta x_q$; the geometry pullback gives precisely
the same coarse pairing with $\Delta x_g$. This finite identity does not
replace $D(r^2)$ by $2rDr$ or assume a continuum spectral product rule.

For the actual retained Dirac operator, a positive-spectrum covariance has

$$
E_D=\langle Q\rho\rangle=M\operatorname{Tr}(C_{\rm source}H_Q)>0,
\qquad M=4\kappa\text{ once}.
$$

Thus static $\dot p_r=C=0$ on a periodic positive geometry would require
$E_D+\mathfrak m\langle Q^2\rangle=0$. Positive source energy, or nonzero
positive magnetic energy, excludes that **static periodic leading class**.
This is the trace/topology blocker: the closed carrier has no boundary
traction term to supply the static balance. It removes a static search in
that class. It does not exclude dynamic, oscillatory or other-source regimes,
and it is not a verdict on NSC.

## What changes on a finite interval

In the smooth continuum interval $[a,b]$, direct integration gives

$$
\boxed{
E_D+\mathfrak m\int_a^bQ^2dx
=\int_a^b QC\,dx+\tfrac12\int_a^b r\dot p_r\,dx
+g\,[rr_x+r^2(\log Q)_x]_a^b.}
$$

Static balance therefore requires the signed exterior/boundary charge
$g[rr_x+r^2(\log Q)_x]$. The assessment verifies its sign as a local total
derivative, separately from the periodic finite SBP identity.

The electrovac product $r=1$, $Q=L=\sec x$, beta zero, on
$|x|<\pi/2$, has $R_h=2$. The audited magnetic-radius ratio gives
$\mathfrak m=g$, and

$$
g[\tan x]_a^b=\mathfrak m\int_a^b\sec^2x\,dx.
$$

It satisfies the leading source-free constraints and radial force exactly.
Adding positive Dirac energy to that unchanged product gives $C=\rho$,
so the background alone is not a self-sourced static construction. A static
Dirichlet value for r is not evidence that a physical parent supplies the
needed traction. Its domain, geometric boundary variation, matching and work
must come from the same action and actual exterior. No wall or exterior
solver is introduced by this note.

## Prototype and interpretation

The finite control uses the existing NF32/NQ128 carrier at period 8,
$Q=.6+.03\cos(2\pi x/8)$ and $r=1.3+.08\sin(2\pi x/8)$, zero momenta,
and six real positive eigenmodes of its actual retained Hamiltonian. Fixed
occupations are $(.75,.75,.5,.5,.25,.25)$, with trace 3 and M=4.
This is a constructor-prepared eigenstate probe, **not source-solved Cauchy
data or a static candidate**. Its residuals deliberately remain nonzero.
It measures the SBP/adjoint identity, Gram/CAR and equality of spectral and
owned source energy. The finite arithmetic differences are indicators, not
continuum or rigorous roundoff-error certificates.

A genuine fundamental Dirac mass from a compact/chiral operator would need
its four-dimensional trace and corresponding radial force derived first.
The angular gap alone does not break this trace-free virial. The separate
[curvature EFT branch](nsc-curvature-eft.md) retains its stated perturbative
scope; the present diagnostic does not restore an independent Weyl mode or
change any coefficient.

## Reproduction

Default calculation is pure and read-only. The explicit writer requires a
frozen producing commit, authenticates source blobs and coefficient/data
inputs, and creates a new record exclusively. The checker repeats the
symbolic and finite control without a trajectory or evidence write.

```sh
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_leading_virial.py
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_leading_virial.py --producer-commit HEAD --write
.venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_leading_virial.py --check
```

The [tests](../tests/test_nsc_discovery_leading_virial.py) cover the boundary
sign/BR checks, actual source/SBP adjoints, narrow static scope and exclusive
record replay. Root owns record generation and any subsequent physical
construction.
