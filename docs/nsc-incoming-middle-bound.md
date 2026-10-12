# Incoming middle-band error from the stored Riccati defect enclosure

The retained middle-band source error has the directed upper enclosure

$$
(\epsilon_\rho,\epsilon_\parallel,\epsilon_{01},\epsilon_\perp)
\le(1.970589\times10^{-6},\;1.770801\times10^{-6},\;
2.154307\times10^{-46},\;7.290434\times10^{-8}).
$$

In lapse/shift action units its bounds are approximately
`(4.1801744e-5, 3.8571317e-45)`. The lapse enclosure exceeds the existing
`3e-11` tolerance, so the result is **OPEN: a valid but too-broad bound**.
It does not establish that the actual physical error is large. Group14's
`[16,160]` interval dominates the density enclosure at `1.8595814e-6`,
followed by groups12 and32. The thermal/scattering contribution is separately
bounded at density `2.386710e-46`; vacuum approximation controls this result.

This calculation reuses the directed radial coefficient bounds in
`9655f6d` and integrates them over the **existing finite middle bands**.
Groups 13 and 14 use `[16,160]`; groups 10–12 and 31–32 use `[40,320]`;
the remaining retained non-LLL groups use `[40,160]`. The metadata of the
actually selected archived panels owns these endpoints.

## Decision and stopping condition

The missing quantity is the integrated physical mode error of the archived
order16 middle recipe. Reuse the stored positive radial defect bounds,
the authenticated channel inventory and the already-owned source trace-norm
bound. Compute finite-energy primitives and project the aggregate density
and current errors into the raw lapse/shift action equations. Compare to the
existing stationarity tolerance `3e-11`; do not enlarge it. Stop after this
one record, its replay and focused tests, even if the bound is broad. No
interval coefficient preparation, mode/scattering solve, physical geometry
selection, source refit, or metric step belongs to this calculation.

Likely obstacles are broad interval enclosures at the lower energies, stale
or inconsistent panel metadata, and conflation of exact-vacuum error with
the physical thermal/reflection correction. These are addressed respectively
by reporting OPEN and group rankings, rejecting provenance mismatches, and
retaining a separate conservative source-state bound.

## Positive finite-energy primitives

For each saved radial coefficient upper bound `I_n`, the projector transport
estimate from the [vacuum-tail owner](nsc-incoming-vacuum-tail-bound.md) uses

$$
J_s(n;L,R)=\int_L^R E^{s-n}\,dE
=\frac{L^{s+1-n}-R^{s+1-n}}{n-s-1},\qquad s=0,1.
$$

The four stress bounds use the same vertices and group factor as that owner:
`2*(J1/a + M*J0)`, `2*J1/a`, zero vacuum current, and
`|lambda|*J0/r`, each multiplying `I_n/2^n`. Directed arithmetic preserves
the stored coefficient bounds and rounds the resulting positive sums
outward. No quadrature or order difference is promoted to an error bound.

The exact source and the archived vacuum-mode recipe each differ from their
own vacuum projector by a trace-norm bound inherited from the horizon and
incoming state law,

$$
b(E)=2e^{-\pi E/\kappa_h}+e^{-2\pi E/(\Omega\kappa_h)}.
$$

For current-normalized mode compression, the difference between these two
thermal corrections is therefore bounded conservatively by `2*b(E)`.
This retains possible reflected horizon coherence rather than treating the
recipe's zero exterior-horizon column as a physical zero. Elementary directed
exponential integrals give the finite-band stress bound. This also bounds
the otherwise omitted thermal/scattering contribution; it does not approximate
reflection coefficients or invent an independent state.

## Constraint units and interpretation

On the unchanged incoming surface `a^2=3*pi/2-4`, `r^2=2`, the existing
source-action map gives

$$
|\Delta\mathcal E_N|\le4\pi ar^2\,\epsilon_\rho,
\qquad
|\Delta\mathcal E_\beta|\le4\pi a^2r^2\,\epsilon_{01}.
$$

Each bound is compared to `3e-11`. A valid enclosure above that tolerance
means this error estimate is too broad to certify the approximation. It is
not a measured residual, a lower bound on physical error, or NON-EXISTENCE.

The local middle recipe is evaluated directly at `rho=1`; it does not use
the finite-offset `working_frame` initialization. Thus no infinity or middle
finite-start term is introduced here. Finite-offset modal accuracy belongs
to the separately resolved low/subgap panels. Low-panel quadrature, remaining
subgap refinements, the angular/compact complement, changed-normal-jet source
matching, and extended stationarity are outside this record.

The [result](../results/development/nsc-incoming-middle-bound.json) includes
the component bounds, action-unit comparison and density ranking. `--check`
authenticates endpoints and inputs and replays only the finite primitives.

```sh
python3 scripts/derive_nsc_incoming_middle_bound.py --write
python3 scripts/derive_nsc_incoming_middle_bound.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_middle_bound.py
```
