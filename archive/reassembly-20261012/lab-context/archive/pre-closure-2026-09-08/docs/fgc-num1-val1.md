# FGC-1-NUM1-VAL1: independent numerical and affine-measurement validation

**FGC-1-NUM1-VAL1** is the implementation-validation gate immediately before
the resolved holdout manifest. It tests the numerical machinery on problems
whose answers are known independently of an FGC-QR collapse outcome. It does
not evolve the FGC-QR holdout, observe regulator activation, form a trapped
surface, or claim defocusing.

In the canonical lower-case contract language, this gate **does not evolve the FGC-QR holdout**.
Each time integrator constructs its own candidate endpoint,
which must pass the composed health and causal checks before it can become the
next accepted state.

The primary method is the frozen diagonal-norm fourth-order SBP operator with
second-order boundary closure and classical RK4. The comparator is a separate
diagonal-norm second-order SBP operator with SSPRK3. Both use a regular-centre
parity extension and independently evaluate the final Runge--Kutta
combination; the last internal stage is never substituted for the candidate
endpoint. Checkpoints encode little-endian binary64 arrays with per-array
checksums and a canonical outer hash.

## Outcome-blind controls

The certificate requires all of the following:

1. the unprojected derivative matrices satisfy their exact diagonal-norm SBP
   identities to the declared binary64 tolerance, and the Kreiss--Oliger
   filters decrease the energy of an injected highest-frequency mode;
2. all six ADM fields preserve exact spherical Minkowski data bitwise under
   both time integrators;
3. an even, centre-regular manufactured wave converges under both methods;
4. the odd field `psi=r*chi` obeys the independent GR-0 spherical scalar-wave
   control with the declared solution and reduction-constraint orders;
5. a parity-regular annular acceleration limit converges to the analytic
   centre value without evaluating the singular spherical reference at
   `r=0`;
6. the new batched source sends the complete grid through the unchanged REF1
   tensor evaluator, recovers the flat zero root, and agrees with the existing
   scalar and exact-rational activated roots after all routes are expressed in
   the same ADM acceleration variables;
7. the two aligned outer domains agree identically inside the retained region
   before either numerical boundary can enter its domain of dependence;
8. every HLT1 typed stop remains present, while a composed HLT1/BND2 runtime
   transaction commits all stages atomically and rolls back both a late-stage
   health failure and a causal-buffer failure;
9. a split checkpoint/restart run is bitwise identical to its uninterrupted
   counterpart under both methods; and
10. an affinely normalized outgoing physical-metric null generator in the
    exact polynomial flat-FLRW control remains null and affine, while the
    direct derivative of its expansion agrees with the complete spherical
    Raychaudhuri assembly. Spherical radial shear and twist are explicitly
    zero by the symmetry and hypersurface-orthogonality specialization already
    proved by DEF0; they are not silently generalized beyond spherical
    symmetry.

The FLRW control is useful because it is time dependent and curved while
remaining analytically specified. It therefore exercises interpolation,
geodesic transport, affine normalization, expansion differentiation, the
warped-product Ricci contraction, and both integrators without importing any
FGC mechanism outcome.

## Scientific boundary

A positive NUM1 result means only that the two declared numerical methods,
the vectorized unchanged-equation source adapter, the fail-closed stage
transaction, restart format, centre limit, and affine observable pass their
independent controls. It is necessary evidence that a later signal is not an
obvious implementation artifact. It is not evidence that such a signal
exists.

The canonical machine decision is
`independent_solver_validation_and_measurement_contract_passed = true`.

The retained-EFT authorization remains false. NUM1 supplies no nonlinear FGC
trajectory, no collapse, no trapped interval, no positive Raychaudhuri margin,
no finite transition, no singularity resolution, no child domain, no dark
sector, and no varying locally measured speed of light. With the canonical
result now present, the preserved PROTO3 RUN1 structure advances from `6/8`
to `7/8`. NUM1 tests low- and near-Nyquist estimator controls; it does not test
whether PROTO3's declared compact pulse passes HLT1. The later
[`FGC-1-CAL0-PREF2`](fgc-cal0-pref2.md) composition test fills that gap and
closes PROTO3 as an executable contract before any holdout outcome.

Reproduce with:

```bash
python3 scripts/reproduce_fgc_num1_val1.py
```
