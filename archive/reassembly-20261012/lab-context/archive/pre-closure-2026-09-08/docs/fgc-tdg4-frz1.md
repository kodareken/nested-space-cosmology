# FGC-1-TDG4-FRZ1: finite-sampling identifiability theorem freeze

`FGC-1-TDG4-FRZ1` freezes the next outcome-neutral mathematical question after
TDG3-PREF18. It consumes only the tracked compact PREF18 record at immutable
commit `63f5278…`; it does not load the CAL10 checkpoint or history arrays,
resume a trajectory, or read SGB-L or FGC-QR data.

## Question being frozen

Let

```text
S(f) = (f(x_0), ..., f(x_63))
```

be the finite sampling map and let `T(f)` be the five continuum Fourier
coefficients in bins `28..32`. TDG4 asks whether `S(f)` alone can place a
finite upper bound on either of the continuum powers used by the old temporal
gate:

```text
P_field      = sum_k 2 |T_k(f)|^2
P_derivative = sum_k 2 (2 pi k)^2 |T_k(f)|^2.
```

This is not a third numerical estimator. It asks whether the target is
identifiable from the information supplied to every such estimator.

## Frozen linear criterion

For a linear sampling operator `S` and target operator `T`, exact
identifiability from the samples is equivalent to

```text
T(h) = 0 for every h in ker(S).
```

Equivalently, `T` must factor through the sample vector: there must be a map
`A` on `range(S)` such that `T = A S`. If instead one direction `h` satisfies

```text
S(h) = 0,     T(h) != 0,
```

then every member of

```text
f_A = f_0 + A h
```

has exactly the same samples as `f_0`, while the target grows without bound as
the unrestricted amplitude `A` grows. The exact rational finite-dimensional
preflight checks this criterion through row-space ranks in both directions: a
factorable target and a target that sees a sampling-kernel direction.

## Frozen smooth witness

For any finite strictly ordered set of sample times, choose an open interval
strictly inside one adjacent-sample gap. Let `phi` be a nonnegative, nonzero
smooth bump compactly supported in that interval, and for any declared target
bin `k` define

```text
h_A,k(x) = A phi(x) cos(2 pi k x).
```

The witness is zero in a neighborhood of every sample, so its value and every
derivative vanish at every sample node. Its coefficient at `k` is

```text
hat(h_A,k) = A (I_0 + I_2k)/2,
```

where `I_0 = integral phi > 0`. The Fourier phase is nonconstant on a
positive-measure interval, so strict triangle inequality gives

```text
|I_2k| < I_0.
```

The coefficient is therefore nonzero for every nonzero `A`. Because the bump
is compactly supported away from the endpoints, integration by parts also
gives

```text
hat(h'_A,k) = 2 pi i k hat(h_A,k).
```

The field and derivative powers both scale as `A^2`. The prospective theorem
must validate this construction independently for all five frozen bins.

## Exact 64-node alias control

TDG4 also freezes a separate exact control for the nominal uniform nodes
`x_j=j/63`:

```text
h_k(x) = sin(2 pi 63 x) sin(2 pi (63-k) x).
```

The first factor is exactly zero at all 64 nodes. Product-to-sum gives complex
coefficients `1/4` at `+/-k` and `-1/4` at `+/-(126-k)`. The declared positive-
bin field power at unit amplitude is therefore exactly `1/8` for every
`k=28..32`. This control is not used as the universal proof; it attacks the
same conclusion through an independent exact construction tied to the frozen
sample count.

## What would be sufficient

The theorem freezes three legitimate routes by which a future admission could
recover a continuum bound:

1. A declared finite-dimensional function class on which sampling is
   injective, together with a quantitative stability bound for the inverse.
2. An independently proved continuum derivative or Sobolev norm. For example,
   if `|f'| <= L`, the piecewise-linear coefficient error obeys

   ```text
   |hat(f)_k - hat(f_PL)_k|
       <= (L/2) sum_j Delta_j^2
       <= L h_max/2.
   ```

3. A continuum evolution-residual theorem and stability estimate strong enough
   to produce a direct coefficient error or one of the required regularity
   bounds.

Every constant must be fixed independently of the observed top-band margin.
The sample values do not supply `L`, a bandlimit, or a continuum stability
constant by themselves.

## What is not sufficient

The frozen attacks explicitly reject the following promotions:

- mere smoothness with no quantitative norm bound;
- agreement of two or more estimators built from the same finite samples;
- adding any finite number of samples without regularity control;
- a small discrete residual without a continuum stability estimate;
- treating a piecewise-linear or quadrature surrogate as an enclosure of all
  smooth histories through the nodes.

Finite densification shrinks the largest gap, which improves a bound once a
derivative constant exists. Without such a constant, another smooth bump can
always be placed inside a remaining gap and scaled independently.

## Freeze boundary

At this checkpoint the theorem design and its exact no-history controls are
frozen, but the separately hash-bound theorem result has not yet executed.

```text
TDG4_sampling_identifiability_theorem_design_frozen = true
TDG4_exact_no_history_controls_pass = true
TDG4_sampling_identifiability_theorem_completed = false
finite_samples_alone_proved_insufficient_for_continuum_top_band_bound = false
current_temporal_spectral_rule_retired_as_continuum_admission = false
replacement_temporal_gate_design_authorized = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
GR0_case_eligible = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

After this freeze is committed, only the read-only `FGC-1-TDG4-PREF19`
theorem binder may execute. A passing theorem would retire the present sampled
spectral rule as a continuum admission and authorize design—not acceptance—of
a replacement premise. It would not make GR-0 eligible or answer any candidate
or physical question.
