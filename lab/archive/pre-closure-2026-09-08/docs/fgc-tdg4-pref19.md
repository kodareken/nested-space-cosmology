# FGC-1-TDG4-PREF19: finite-sampling identifiability theorem

`FGC-1-TDG4-PREF19` independently executes the theorem frozen by
`FGC-1-TDG4-FRZ1` at immutable commit `8d0df53…`. It reads only tracked
configs and compact canonical records. It does not load a campaign checkpoint,
inspect a terminal history array, resume or advance a state, or read an SGB-L
or FGC-QR trajectory.

## Result

Let `S` map a continuum history to its values at finitely many sample nodes,
and let `T` map that history to the continuum Fourier coefficients used by the
declared top-band field and phase-derivative powers. On
`C_c^infinity(0,1)`, finite point samples alone do not place a finite upper
bound on either power.

The exact linear criterion is

```text
T is determined by S  iff  T(h)=0 for every h in ker(S).
```

Necessity is immediate: if `T=A S`, every sampling-kernel direction is also
annihilated by `T`. For sufficiency, define `A(S(f))=T(f)`. This is
well-defined because two histories with the same samples differ by an element
of `ker(S)`. The certificate attacks both directions with exact rational
matrix controls: injective sampling, noninjective but factorable sampling, and
a target that detects a sampling-kernel direction.

For any finite ordered sample set there is an open gap between adjacent
samples. Choose a nonnegative, nonzero smooth bump `phi` supported strictly in
that gap and define, for each target bin `k`,

```text
h_A,k(x) = A phi(x) cos(2 pi k x).
```

The witness and all of its sample-node derivatives vanish, yet

```text
hat(h_A,k) = A (I_0 + I_2k)/2 != 0
```

for every nonzero `A`. The inequality is strict because the Fourier phase is
nonconstant on the positive-measure support of `phi`, so
`|I_2k| < I_0`. The derivative coefficient is
`2 pi i k hat(h_A,k)`. Both declared powers therefore grow as `A^2` while the
complete sample vector remains unchanged. The construction is checked for all
five bins `28..32`.

An independent uniform-grid alias supplies a second exact route. At the 64
nodes `x_j=j/63`,

```text
sin(2 pi 63 x) sin(2 pi (63-k) x)
```

vanishes at every node but has target coefficient `1/4` and unit-amplitude
positive-bin field power `1/8` for every declared bin. These results match the
prospectively frozen FRZ1 controls exactly.

## What has been retired

PROTO4 introduced a temporal rule that accepted or rejected a continuum
trajectory using 64 sampled history values, fixed top-band budgets, and a
nested sampled-tail ratio. PROTO13 inherited that rule unchanged, and CAL10
preserved its historical failure at the first 64-sample audit.

PREF19 now retires that rule **only as a continuum-admission test**. The reason
is not that CAL10 should have passed. The reason is that, on the function space
the rule implicitly left unrestricted, identical finite samples permit
arbitrarily different continuum top-band power. A sampled statistic therefore
cannot by itself certify the required continuum claim.

The rule remains available as a descriptive sampled diagnostic. The immutable
PROTO13/CAL10 record remains exactly what occurred under its then-frozen
contract:

```text
historical_PROTO13_stop_preserved = true
historical_CAL10_result_reclassified = false
historical_thresholds_changed = false
actual_history_classifications_changed = false
```

PREF19 does not infer that the actual histories contained any particular
between-sample oscillation. It proves an information boundary, not a hidden
trajectory.

## What the theorem does not cover

The theorem is universal on the declared unrestricted smooth function space.
It does not prove that the actual PDE solution manifold contains the scaled
bump or alias families. A future gate may succeed precisely by proving that
the equations, initial data, constraints, and health domain restrict histories
to a quantitatively smaller class.

That replacement must independently establish at least one sufficient route:

1. a stable injective finite-dimensional class;
2. a quantitative derivative or Sobolev bound; or
3. a continuum evolution-residual plus stability estimate.

Merely taking more finite samples, assuming smoothness, comparing sampled
surrogates, or observing a small discrete residual does not provide such a
bound. Every replacement constant must be fixed independently of the observed
target margin.

## Scientific boundary

The theorem closes a numerical-admission design question and authorizes the
next prospective design step. It does not accept a replacement, make GR-0
eligible, open a candidate branch, or answer the mechanism question.

```text
TDG4_sampling_identifiability_theorem_completed = true
finite_samples_alone_proved_insufficient_for_continuum_top_band_bound = true
current_temporal_spectral_rule_retired_as_continuum_admission = true
current_temporal_spectral_rule_retained_as_sampled_diagnostic = true
actual_unsampled_power_of_PROTO13_histories_inferred = false
PDE_solution_manifold_restriction_proved = false
replacement_temporal_gate_design_authorized = true
replacement_temporal_admission_defined = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
classical_spherical_diagnostic_authorized = false
SGBL_execution_authorized = false
FGCQR_holdout_execution_authorized = false
```

The next task is prospective replacement-gate design. It must choose and
quantify a sufficient continuum-control route before any new calibration or
candidate trajectory is authorized.
