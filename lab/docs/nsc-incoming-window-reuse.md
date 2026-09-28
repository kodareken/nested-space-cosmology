# Two source windows using already-certified order24 radial coefficients

## Decision and stopping condition

Reuse the group12 order24 radial certificate from `94dba54` and the group22
certificate from `bed8b1d`. Their positive radial defect integrals have no
energy dependence. A separately recorded interval-binding view therefore
permits using the group12 certificate on LOW `[32,40]` and the group22
certificate on middle `[40,160]`, with matching order24 source approximations.
The original endpoint labels, coefficients and result bytes remain intact.

The missing connection is an explicit source update on these two intervals:
direct order24 vacuum minus the selected archived LOW contribution for
group12; order24-minus16 vacuum correction with the old thermal insertion
retained for group22. Reuse the actual selected panels and signed factors.
Use positive48/64 quadrature and50/70-digit precision controls. Compare each
physical lapse/shift component bound to the unchanged `3e-11` tolerance.
Stop after these TWO source updates and replay/tests/local commit. No radial
recurrence, mode/scattering solve, or full-group sweep belongs to this record.

Likely issues are mismatched endpoint metadata, accidentally changing the
physical source while improving its approximation, and treating numerical
quadrature differences as rigorous error bounds. Reject metadata mismatch,
retain original values and explicit thermal policies, and keep each error
category separate.

## Reusable window API

`nsc_incoming_window_refinement.py` provides authenticated panel selection,
the local order24 source/correction functions, positive source quadratures,
and explicit reuse of a saved radial certificate. It supports the common
retained LOW and middle families, including later group32 work; only groups12
and22 are calculated in this record.

The local `rho=1` coefficient algebra reuses `_riccati_at_one`. It does not
reconstruct radial defect integrals. For the direct LOW source, the original
`_mp_bloch` fourth-order subtraction remains unchanged. For middle corrections,
the same subtraction cancels in `Tr[(P24-P16)V_A]`; stable Bloch differences
avoid subtracting unit traces. The old16 source convention is independently
checked against the stored middle vacuum kernels. Both angular signs and
`incoming_group_factor` enter once.

## Physical source and error accounting

The declared affine-horizon source law, incoming occupation, metric, normal,
mass and angular labels are unchanged. Group12's new numerical insertion
contains its vacuum term; the thermal term remains a nonzero bounded remainder.
Group22 retains its archived thermal insertion. In both cases the existing
conservative `2*b(E)` trace-norm estimate covers the thermal/scattering
remainder, including possible horizon coherence.

For each window, a fresh bound-consumer view changes only the endpoint tag
required by the integration API and records the original tag and requested
interval. This is numerical interval bookkeeping, not a physical cutoff or
a mutation of the old certificate. The saved order24 radial integrals and
local cross-product coefficients are used without regeneration.

The rank-one energy bound certifies vacuum density/lapse and vacuum current;
thermal remainder bounds are included in the displayed component budget.
Pressure-mode error, rigorous numerical source quadrature error, other
windows and full source/constraint stationarity remain separate. A numerical
PASS does not close those gaps or authorize metric evolution.

```sh
python3 scripts/derive_nsc_incoming_window_reuse.py --prepare
python3 scripts/derive_nsc_incoming_window_reuse.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_window_refinement.py
```

The [result](../results/development/nsc-incoming-window-reuse.json) retains
before/after values, signed contributions, explicit additive corrections,
reused-certificate provenance, physical bounds and numerical indicators.
Replay contracts the authenticated artifact without preparing coefficients,
radial recurrences, modes or quadrature nodes.
