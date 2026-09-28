# Group22 source on the existing LOW interval [32,40]

## Decision and stopping condition

The authenticated group22 LOW24-to-LOW32 density change on `[32,40]` is
`+9.03092560547e-11`, larger than its net full-LOW change. Those values use
packet-inverted physical fields. Reuse both archived quadratures, the same
affine-horizon source law, the same ad4 subtraction, and the local Riccati
recurrence at numerical order24. The exact gap is an accurate source
approximation on this one existing interval, with a physical vacuum-energy
certificate in the existing lapse-action tolerance `3e-11`.

Use the frozen `trimmed_defect_jets` helper with128 radial cells, centered
Taylor depth4, and the existing rank-one energy bound. Evaluate positive
48- and64-node source quadratures with independent50/70-digit controls.
Stop after one group22 record and focused verification, whether the energy
enclosure fits the tolerance or remains OPEN. No old field/scattering
generator, full-group sweep or physical-state selection is authorized here.

Likely obstacles are excessive defect-enclosure width, conflating numerical
quadrature agreement with rigorous integration error, and silently replacing
the physical covariance. Keep all three explicit: retain an OPEN result if
the bound is broad; separate quadrature indicators from physical mode-error
bounds; preserve the original LOW sources and record the additive change.

## Same affine-horizon source, improved numerical representation

The [spatial mode owner](nsc-pg-massive-mode-resolution.md) defines the source
as the affine-horizon covariance and inherited incoming occupation. Its
interior vacuum is the normalized partner projector. The finite matched
collar and later packet inverse are numerical representations of that source,
not different free physical states.

At `rho=1`, the new local order24 polynomial gives the partner Bloch vector.
The original generated fourth-order adiabatic expressions are subtracted
at the same incoming geometry. Mass, signed angular labels, state scales and
the negative-frequency/group factor are unchanged. Both angular signs enter
exactly once. The existing LOW24/LOW32 kernels, weights and selected intervals
are authenticated and retained; the new source is an explicitly different
numerical approximation with its own error account.

The new numerical insertion contains the vacuum contribution only. The
thermal source is **not declared physically zero**: the existing horizon
coherence and incoming occupations give the nonzero trace-norm bound
`b(E)=2 exp(-pi E/kappa)+exp(-2 pi E/(Omega kappa))`. This record conservatively
reuses the already-owned `2*b(E)` finite-interval bound as an omitted thermal
remainder. It neither overwrites the source covariance nor discards its
coherence as a physical statement.

The finite-offset initialization uncertainty belongs to the old LOW method.
The new certificate instead starts from the same declared affine-horizon
projector limit and does not assign a zero error to the old finite-offset
archive. Differences from the old sources are not identified as recovered
"noise" or as bounds on their errors.

## Matching order24 source and order24 certificate

The frozen trimmed recurrence encloses `f24,...,f48` over each complete radial
cell. Centered Taylor remainders are intersected with the natural interval
rectangles, then integrated with the existing positive horizon weight.
The [rank-one energy owner](../src/recursive_horizons/nsc_incoming_projector_energy_bound.py) uses
the matching order24 local cross-product coefficients:

$$
|\Delta\rho|\le2|h\times n_{24}|\,d+2|h\cdot n_{24}|\,d^2,
\qquad d=\|P_{\rm affine}-P_{24}\|.
$$

Positive energy-power primitives integrate this bound on `[32,40]`.
The lapse conversion is `4*pi*a*r^2`, with `a^2=3*pi/2-4`, `r^2=2`.
The exact vacuum current difference is zero because both projectors have
unit trace; the physical current still has its bounded thermal remainder.
Pressure-source values are recorded, but the rank-one certificate here
specifically concerns density/lapse and vacuum current.

The source quadrature and precision differences remain numerical controls.
A passing vacuum-plus-thermal physical error budget is not, by itself,
a rigorous bound on numerical source integration or full-source convergence.
Remaining groups, other LOW/subgap intervals, altered normal jets and
extended stationarity are not certified by this record.

```sh
python3 scripts/derive_nsc_incoming_low_high_source.py --prepare
python3 scripts/derive_nsc_incoming_low_high_source.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_low_high_source.py
```

Replay authenticates and contracts the stored source/coefficient artifact;
it does not rebuild interval recurrence jets, modes or quadrature nodes.
