# Explicit order24 correction to the group14 middle vacuum approximation

## Decision and stopping condition

Reuse the authenticated group14 signed middle panels, their order16 vacuum
source and existing `[16,160]` interval, and the local recurrence
`nsc_incoming_source_tail._riccati_at_one(mass, angular, order=N)`.
The [middle-bound record](nsc-incoming-middle-bound.md) shows that the
ordinary order16 absolute-defect enclosure cannot certify the existing
`3e-11` lapse tolerance and group14 dominates that enclosure. This is the
specific reason to change the **numerical Riccati order** to24. The physical
operator, horizon source, state parameters, geometry and fourth-order
subtraction remain unchanged.

The missing quantity here is the explicit finite-band vacuum-source
correction. Evaluate `order24-order16` on the existing24-node panels and an
independent48-node positive quadrature over the same interval, then check
arithmetic precision. Stop after this correction record and focused replay
verification. No archived source is overwritten, no mode or scattering
equation is solved, and no physical source convergence follows from an
order difference. The independent physical order24 defect bound remains
OPEN in this record.

Three likely issues are cancellation in nearly equal projectors, missing
angular-sign or multiplicity factors, and confusing numerical quadrature
agreement with an exact-mode error certificate. Use stable projector
differences, verify the original order16 source convention, retain both
actual angular signs and the owned factor once, and record these error
categories separately.

## Same subtraction cancels exactly

For the same raw vertex `V_A` and the unchanged formal fourth-order projector,

$$
\Delta K_A=\operatorname{Tr}[(P_{24}-P_{16})V_A].
$$

The reference subtraction cancels algebraically. Writing
`S24=S16+dS`, `u=a^2|S16|^2`, `D=1+u`, and
`du=a^2(2 Re(conj(S16)*dS)+|dS|^2)`, the stable Bloch difference is

$$
\Delta b_x=\frac{2a(D\Im dS-\Im S_{16}\,du)}{D(D+du)},\quad
\Delta b_y=-\frac{2a(D\Re dS-\Re S_{16}\,du)}{D(D+du)},\quad
\Delta b_z=-\frac{2du}{D(D+du)}.
$$

The four correction kernels in the existing order are

$$
\left(-m\Delta b_x+\lambda\Delta b_y/r-E\Delta b_z/a,\;
-E\Delta b_z/a,\;0,\;\lambda\Delta b_y/(2r)\right).
$$

The exact zero current follows from equal unit projector traces; no numerical
small-value threshold is used. The two angular signs use the unchanged
`incoming_group_factor` once. The source16 baseline is verified against the
archived vacuum kernels with the original ad4 owner, while the correction
itself does not reevaluate that subtraction.

## Additive numerical correction and remaining error

This is a vacuum-only correction. The already evaluated thermal insertion
stays as archived, with its conservative exact-versus-approximate difference
bound retained from the middle-bound record. The record contains the old16
group14 middle source, the new additive vacuum correction and the corrected
approximation separately. It does not relabel the old artifact as order24.
The rest of the finite source and the already certified order16 tail above160
remain separate numerical contributions to the same physical source.

Stored-grid versus refined-grid differences and precision differences are
numerical controls. They are not rigorous quadrature or exact-mode error
bounds. A physical source certificate still requires the independent
order24 transport-defect bound, any necessary rigorous correction-integral
accuracy, and the other unresolved finite spectral contributions.

```sh
python3 scripts/derive_nsc_incoming_middle_order_correction.py --prepare
python3 scripts/derive_nsc_incoming_middle_order_correction.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_middle_order_correction.py
```

Preparation uses only local coefficient and positive quadrature algebra.
Replay authenticates the coefficient/quadrature artifact and checks its
contractions without regenerating coefficients or quadrature nodes.
