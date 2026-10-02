# Ambient-extent intervention

This owner tests whether candidate source return timing depends on the periodic
carrier. High tagged return alone is not holding: the period-eight candidate
can simultaneously have almost collapsed proper child length and large actual
tides. Probability, geometry, constraints and energy remain separate readouts.
The code is [nsc_discovery_extent.py](../src/recursive_horizons/nsc_discovery_extent.py),
with [driver](../scripts/derive_nsc_discovery_extent.py) and
[focused tests](../tests/test_nsc_discovery_extent.py). Existing producers and
evidence are preserved.

## Native physical domains

The actual Galerkin constructor already supports `length=L`, including periodic
geometry derivatives, half-integer AP momentum, interpolation and half-density
maps. The adapter calls `build_grid(nf, quadrature=4*nf, length=L, gauge="conformal")`.
It does not patch period-eight operators onto longer coordinate arrays.

| Period | Fermion count | Quadrature | dx_f | dx_q |
|---:|---:|---:|---:|---:|
| 8 | 256 | 1024 | 1/32 | 1/128 |
| 12 | 384 | 1536 | 1/32 | 1/128 |
| 16, conditional confirmation | 512 | 2048 | 1/32 | 1/128 |

Geometry has `ng=nf-1`; its spacing and maximum wavenumber change slightly.
These are nominally matched resolutions, not identical geometric bands. Those
spacings are serialized. The fermion spectrum remains
`k_m=2*pi*(m+1/2)/L`; a full-period spinor translation is minus the identity.

The separated source and fixed rank-six observer retain their physical widths,
occupations, actual owning carrier parameter and local phase origins. They are
not dilated with L. A full spectral/AP translation by **two units** moves the
source, observer and columns of the complete geometry basis W. The parent is
`(2,6)`, child `(3,5)`, source supports `[2,3],[3,5],[5,6]`, and clock worldlines
are `3,4,5`. They stay away from the seam. Translation is on the fermion and
quadrature lattices; the odd geometry grid is translated spectrally. New nodal
and fine detail tails are measured rather than retaining an old zero-tail label.

The period-eight adapter audits source covariance against an independent
integer AP node roll, and checks source forces against the existing nested
constructor after a periodic fine-grid roll. The actual carrier parameter is
read from `coupling.CARRIER_K`, not remembered from a review. This audit concerns
the source/operator construction; it does not substitute an old radius for the
new initial solve.

## Own preparation and protected mismatch

Each domain gets its own source-consistent dense initial radius solve, with
constant initial Q0, zero chi and zero geometric momenta. Full current and shift
residuals remain reported. No mean current is deleted and no old radius is
padded into a larger domain.

The periodic elliptic constraint sees a changed source mean as L changes.
Therefore exact local matching is not guaranteed. The saved protected neighborhood
is the parent plus a collar, **`[1,7]`**. Common physical quadrature points retain
r, Q, physical N and q (`r*Q` on this chart), first and second spatial jets,
source rho, constraints and initial actual tidal components. Both T0 and
handoff protected comparisons report absolute and scaled mismatches. A mismatch
beyond the declared numerical-equivalence tolerance is labelled ambiguous for
pure causal attribution. That is a preparation limitation, not a continuation
veto or an invented physical error bound.

The two control branches of each domain begin with that domain's identical
full **T=0.3** state, field, momenta and integrated normal clocks. There is no
source reset or second solve at handoff. Coupled and frozen geometry then reuse
the existing owned RK4 steppers. Normal clocks use the same per-step endpoint
trapezoid protocol as the episode controller.

Conformal frozen-field H has kinetic factor one and potential **kappa times
its own handoff Q(x)**. It is independent of r directly. It need not equal the
common constant initial Q0 after a coupled prefix. The protected handoff Q
comparison makes this qualification explicit. Frozen return timing is useful,
but these differently prepared potentials are not an exact matched-locality
proof.

## Observation, checkpoints and budget

The period-eight discovery atlas rejects other periods and fixed-window
translations. This adapter instead uses underlying period-aware metric,
interval, matter-energy and tidal helpers. It reports tagged columns 2 and 3
inside the common shifted child `(3,5)`, total probability/CAR/Gram, proper
lengths and r/Q/N, source current, Noether flux, constraints, full field and
geometry energy, and actual scalar/radial/angular tides. Chi is not substituted
for curvature. Tagged ancestry is a computational label, not an independent
energy parcel. No probability-return threshold certifies holding or renewal. The final
assessment retains sampled post-parent-crossing peaks (including their brackets),
flags maxima at the recorded boundary, and reports period12-minus-period8 values
at common stations plus coupled-minus-frozen differences within each extent.
It leaves `extent_discriminates` for root review rather than asserting an echo.

Diagnostics occur at cadence **at most 0.1** and at stations **8 and 12**.
Each station or stop commits an immutable full checkpoint and its observations.
Budget and chart exits preserve the last admissible arrays without clamping or
non-existence claims. Checkpoints serialize actual L as both metadata and a
hashed array; reconstruction uses L and nq explicitly, never `nf` alone or the
old period-eight loader. Source, observer, W, weights and extent pins are checked
before continuation commits. JSON/NPZ chunks are capped at 64 MiB.

Preparation and the short coupled prefix are serial. Long cases use exactly
one root-owned process pool, at most six workers, with one-thread math scopes.
The aggregate budget is six CPU-hours, including measured preparation, prefix,
child CPU and coordinator CPU. Immutable station records retain cumulative
child CPU for recovery; a small coordinator/commit reserve is kept before
allocating the remaining child allowances. Each case gets a disjoint remaining allowance;
the forecast multiplier is **1.5**. The existing field plus geometric principal
step restriction is reevaluated on the actual state and L. It is an admission
restriction, not a stability or continuum certificate.

Default execution compares period eight and twelve in both modes: **four cases**.
A read-only `audit_existing_period8` API can assess an explicitly identified old
checkpoint against the newly prepared comparator. Its default reads the latest
checkpoint; `ordinal=0` can select a handoff. It audits translated source covariance,
physical state, geometric means and physical rates. Admission also requires the
candidate's declared time to equal the legacy checkpoint time. It does not automatically import
or extend old evidence. Reusing old period-eight continuations requires matching
inputs, source and physical states, plus explicit root integration of the audited
comparator. Otherwise the new four-case batch remains the declared experiment. The known
late comparator is `nsc-discovery-crossing-v1`, cases
`nf256_coupled_dt0.0005` and `nf256_frozen_geometry_dt0.0005`; upstream handoffs are
in `nsc-discovery-episode-v1`. These references do not authorize automatic reuse.

Period sixteen is not part of the default batch. Its preparation requires a
root-reviewed admission JSON containing `extent_discriminates:true`; this decision
must follow the period-eight/twelve comparison. The adapter never automatically
opens it because a return fraction is high.

## Commands

No flags prints the interface preview and writes nothing. After root freezes
and commits the producer, run from the repository root. The scientific launcher
changes cwd to `lab`, so relative evidence paths start with `results/`:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --prepare --output results/development/nsc-discovery-extent-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --materialize --production --output results/development/nsc-discovery-extent-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --run --production --workers 6 --output results/development/nsc-discovery-extent-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --check --output results/development/nsc-discovery-extent-v1
```

The focused tests use smaller matched grids and coordinate times at most 0.01:

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_extent.py -q
```

## Real periodic observation successor

The frozen v1 batch was produced by commit
`697efe0155cdbf796c0174468c90a4c368bb7dbe`. L8 coupled and both frozen branches
reached T12. L12 coupled retained its exact last state at
`T=6.699999999999557` after the normal-energy observation raised
`interval integral did not preserve a real density`. This was an observation
failure, not a chart-exit classification. All original v1 manifests and chunks
remain immutable. `--check` authenticates their historical producer Git blobs;
it does not require today's observer to have the old producer's bytes.

The offending normal-energy density had maximum magnitude about `4.12e8`.
Contributions with absolute sum about `2.34e8` cancelled to a child integral
near `-0.221255`; a complex FFT sum left an imaginary rounding remainder about
`1.92e-8`. The cut endpoints are even grid nodes, so the lone Nyquist mode's
imaginary contribution was negligible in this specific failure. This does not
justify relaxing an imaginary-output threshold. It calls for the explicitly
real interpolant of real nodal data.

`real_interval_integral(grid, values, interval, return_indicator=False)` and
`real_periodic_values(grid, values, coordinates)` now form the real conjugate
pairs and the even-N Nyquist cosine. If `c=RFFT(values)/N`, `omega_k=2*pi*k/L`,
the integral is

\[
c_0(b-a)+2\operatorname{Re}\sum_{k=1}^{N/2-1}
c_k\frac{e^{i\omega_kb}-e^{i\omega_ka}}{i\omega_k}
+c_{N/2}\frac{\sin(\omega_{N/2}b)-\sin(\omega_{N/2}a)}{\omega_{N/2}}.
\]

Inputs must be finite real nq-vectors. A conjugacy check uses coefficient-scale
roundoff rather than a tiny cancelled integral as its scale. Integration uses
the sinc form of the primitive and compensated real summation. No mode is
discarded and the owning derivative's zero Nyquist symbol remains unchanged.
Point evaluation uses that same cosine convention. The integral indicator
retains the conjugacy gap/tolerance, absolute term sum, cancellation factor and
`eps*absolute_term_sum`. That last number is a sensitivity indicator, not a
Fourier interpolation, propagated-state, or total arithmetic error enclosure.
The extent's physical metric lengths and clock-rate observations use these
real primitives; Hamiltonian rates and RK4 are unchanged.

The observer also retains the constraint pair
`h_c=Q*(C_g+rho)` and `D=D_g+current`, alongside raw `C_g+rho`, projected and
held-out h_c, and the like-unit source scales `Q*rho` and `current`. A large raw
C when Q is tiny is not by itself a chart failure or proof of physical
instability. Ratios to the corresponding source density are descriptive; no
bound maps the sampled constraint pair or RK-stage residual into the claimed
tagged probability or proper-length effect. `observable_error_bound` and the
stage constraint bound remain null. The assessment explicitly leaves
constraint propagation and late curvature resolution unresolved.

## Exact L12-only recovery

`prepare_resume(output, reference)` is creation-only and refuses an existing
exact output directory or an output inside its reference. It accepts only the
completed frozen v1 batch and its latest named L12 coupled T6.7 observation
stop. It authenticates historical Git producer bytes, every reference JSON/NPZ,
the initial solver, prefix-state hash, source/W/observer/weights, physical L/nq,
windows, and saved continuous clocks. Current dynamical/source owners must
still match the frozen producer; only this observer adapter and its driver may
have changed. The new handoff contains every old array byte-for-byte, including
canonical momenta, spinors and clocks. It calls no initial solve and resets no
prefix, physical time, or step count.

The prepared successor has one runnable case, `L12_coupled`, with stations 8 and
12. Its run uses the existing RK4, controls, restrictions and aggregate budget,
with one executor worker. Historical aggregate CPU and the saved case's CPU
are counted once. Original L8 coupled and both frozen branches are authenticated
summary references and are never resubmitted. Assessment merges those old
observations with the new L12 continuation. Period16 is not admitted by this
recovery: the observed frozen return shifts from eight to twelve already, while
the roughly twenty-percent protected initial-radius mismatch remains an
attribution limitation.

After root freezes and commits the corrected producer, choose the exact new
directory and run:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --check --output results/development/nsc-discovery-extent-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --prepare-resume --reference results/development/nsc-discovery-extent-v1 \
  --output results/development/nsc-discovery-extent-v2
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_extent.py \
  --run --production --workers 1 --output results/development/nsc-discovery-extent-v2
```

Development tests authenticate and observe the frozen state, check manufactured
Nyquist/cancellation cases, and prove exact array handoff. The recovery test
uses `max_steps=0` and forbids the production stepper and initializer; it creates
no new production trajectory.
