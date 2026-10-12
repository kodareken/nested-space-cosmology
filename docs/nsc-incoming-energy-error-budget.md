# Projector geometry in the incoming lapse-error budget

The same stored eight-cell coefficients now give a total middle-band
lapse-error enclosure of `1.15107e-7`, reduced from `4.18017e-5` by the
energy-specific projector geometry. Group14's128-cell bounds are
`3.24091e-10` at numerical order16 and `2.51409e-10` at order24. The latter
still exceeds `3e-11`; the rigorous component gate remains OPEN.

Combining that new order24 group14 bound conditionally with the other31
groups and the retained thermal bound gives approximately `4.998e-9`.
Using the order24 approximation additionally requires the explicit
[source correction](nsc-incoming-middle-order-correction.md), whose lapse
change is `-3.30349e-12`. The old source receipt is preserved. The remaining
uncertainty is now dominated by groups12 and32 in the unchanged coarse
enclosures, with group14 still requiring a tighter component certificate.

For the existing normalized vacuum projector `P=(I+n.sigma)/2`, exact vacuum
`P*`, and `d=||P*-P||op`, the same rank-one geometry gives

$$
|\operatorname{tr}[H(P_*-P)]|
\le 2|h\times n|d+2|h\cdot n|d^2,
\qquad H=h\cdot\sigma.
$$

The parallel change is quadratic: `n.(n*-n)=-2*d^2`. This imports the
rank-one projector identity into the **actual incoming energy vertex**,
instead of bounding every direction by `2*||H||*d`. Both vacuum traces are
one, so their current error is exactly zero. Thermal contributions retain
their independently bounded source-law correction.

At the incoming surface, `h=(-m,lambda/r,-E/a)` and the owned Riccati mode
has `n=(2*a*Im(S),-2*a*Re(S),1-u)/(1+u)`, with `u=a^2*|S|^2`.
Writing `x=1/(2E)`, the inherited `c1=lambda/r+i*m` cancels the constant
cross-product terms exactly. The remaining numerator coefficients begin at
positive powers of x, giving a directed polynomial upper bound `B(E)` for
`|h cross n|`. The existing radial defect coefficients bound `d` by `D(E)`.
Consequently the new energy bound is

$$
2B(E)D(E)+2(E/a+M)D(E)^2,
\qquad M=\sqrt{m^2+(\lambda/r)^2}.
$$

Positive coefficient products are integrated analytically on the same
authenticated finite middle intervals, with the same signed multiplicity.
No negative-frequency factor is added again. The source-action projection
remains `4*pi*a*r^2` for the lapse equation. This improves the energy bound;
it does not silently improve pressure bounds.

## Reuse, intervention and stopping condition

Reuse the committed radial defect certificate and its finite-band mapping,
the exact incoming Hamiltonian, and the existing Riccati recurrence. The
missing connection is a useful error bound on the lapse source at `3e-11`.
The original eight-cell middle enclosure is valid but too broad. A group14
sixteen-cell control and a centered Taylor enclosure test whether interval
dependency alone explains the gap. A 128-cell energy-specific control and
an explicitly separate numerical order24 representation then test the
remaining approximation error. They change neither the action nor the
physical state law, radius, flux or other scale.

The order16 source and all old artifacts are preserved. A numerical-order24
source is usable only with its own defect bound and the separately recorded
additive source correction; it is never pasted into an order16 receipt.
The fourth-order adiabatic reference is unchanged. Stop each experiment when
its bound resolves the stated budget decision or identifies the obstruction;
do not repeat an unchanged enclosure merely to continue.

Material risks are interval dependency through many derivative recurrences,
loss of oscillatory cancellation in the positive defect estimate, and
asymptotic-series growth at higher numerical order. Each remains visible in
the comparison rather than being absorbed into a new tolerance.

## Centered interval control

The [centered owner](../src/recursive_horizons/nsc_incoming_defect_taylor_bound.py)
retains exactly sixteen inverse-energy terms. Only derivative bookkeeping
depth increases. For a cell I and a center c it intersects the natural
rectangular complex enclosure with

$$
\sum_{j=0}^{p-1}\frac{f_n^{(j)}(c)}{j!}(I-c)^j
+\frac{f_n^{(p)}(I)}{p!}(I-c)^p.
$$

The stored Taylor jets already include `1/j!`. Tests retain that convention,
the existing coefficient recursion, and full directed collar coverage.
The depth4, sixteen-cell control barely improves the raw group14 lapse bound
(`4.37e-6` to `4.35e-6`), so it is retained as an unsuccessful enclosure
strategy at that resolution. This is information about the estimate, not
a physical exclusion. Point-sample diagnostics used to choose the next
numerical experiment are explicitly not certified bounds.

## Evidence and scope

The [record](../results/development/nsc-incoming-energy-error-budget.json)
keeps the old and new numerical representations separate, with directed
linear and quadratic error contributions, source hashes and raw radial
coefficient inputs. Its replay does not rerun radial coefficient generation,
mode equations or scattering. The original source is unchanged until its
explicit correction is composed by a later owner.

Low-panel quadrature/modal accuracy, finite-collar errors, remaining subgap
refinement, the retained inventory's scope, full incoming constraints and
extended stationarity remain separate. No physical initial geometry,
endpoint family, extra action term or metric timestep is selected here.

```sh
python3 scripts/derive_nsc_incoming_energy_error_budget.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_incoming_defect_taylor_bound.py tests/test_nsc_incoming_projector_energy_bound.py tests/test_nsc_incoming_higher_order_defect.py
```
