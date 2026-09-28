# FGC-1-TDG1-FRZ1: prospective temporal-gate diagnosis

`FGC-1-TDG1-FRZ1` freezes the first diagnosis permitted by
FGC-1-CAL10-PREF15. It advances no state, loads no terminal history array,
changes no historical threshold, defines no replacement admission, and opens
neither SGB-L nor FGC-QR.

## Why this backtrack is necessary

PROTO13 stopped at its first 64-sample temporal-spectrum assessment. Its
spatial and constraint gates still passed, but both methods failed a temporal
rule that compares normalized top-band power fractions across *spatial*
resolutions.

The synthetic preflight exposes a structural problem with that quantity. The
historical rule is applied to six controls on three spatial resolutions:

```text
zero
constant
affine
one low-frequency sine
one smooth pulse
one low-frequency sine whose amplitude contracts by four per refinement
```

The zero history passes. Every nonzero control passes every individual
bandwidth budget but fails the nested-ratio layer. When a nonzero history is
identical at all spatial resolutions, every normalized field- and
derivative-tail ratio is exactly `1`. The same remains true when the entire
signal amplitude contracts by four at each refinement—and its power therefore
contracts by sixteen—because normalizing each tail by its own total power
cancels the amplitude convergence.

Therefore:

> Normalized raw temporal tail fractions are not a general spatial-convergence
> observable for a nonzero continuum temporal history.

This conclusion is independent of the PROTO13 history values. It does not say
that all temporal checks are unnecessary, that the individual bandwidth
budgets pass, or that the GR-0 trajectory is valid. It says only that this
cross-resolution ratio cannot carry the meaning assigned to it.

## Frozen raw-history diagnosis

The separately authorized TDG1 execution will consume the immutable CAL10
checkpoint once and record, without admission promotion:

- the inherited compact-window budget for every tracer and field;
- rectangular, mean-centered, and affine-detrended diagnostic views;
- the absolute FFT amplitude and the already declared `1e-30` amplitude floor;
- exact-zero, subnormal, below-`epsilon^2`, below-`epsilon`, and resolved-power
  bins;
- 16-, 32-, and 64-sample resampling-sensitivity views, explicitly not called
  independent temporal-convergence evidence;
- per-signal proper-time interpolation round-trip debit;
- three-resolution history differences on the common proper-time
  intersection; and
- their raw observed order against the already declared `3/2` diagnostic
  threshold.

These views separate questions rather than choosing the most favorable
answer. In particular, a detrended or downsampled pass cannot retroactively
pass PROTO13. A raw history-difference order cannot become a replacement gate
without a later prospective protocol.

## Immutable boundary

The freeze binds:

```text
CAL10/PREF15 commit:  fa56acf4d8c6297eb6f9b5dc8afdb5e2b5ce9020
terminal event:       63
terminal time:        63/16
history samples:      64
tracers:              48
fields:               6
historical ratio:     1/4, unchanged
amplitude floor:      1e-30, inherited
raw-order diagnostic: 3/2
```

The checkpoint is read-only and the terminal PROTO13 campaign may not resume.
TDG1 may inspect the frozen arrays, but it may not evolve them.
The ignored raw checkpoint may be absent in a clean clone. Its hash is then
verified through the immutable CAL10 certificate and commit lineage; whenever
the local raw file is present, TDG1 verifies its bytes directly and fails
closed on any mismatch.

## Claim boundary

At this freeze:

```text
TDG1_diagnostic_contract_frozen = true
TDG1_synthetic_counterexample_confirmed = true
TDG1_terminal_histories_diagnosed = false
replacement_temporal_admission_defined = false
PROTO14_frozen = false
fresh_GR0_dynamic_calibration_completed = false
GR0_case_eligible = false
```

The structural counterexample rejects one numerical observable, not GR,
FGC-QR, Finite Gradient Closure, thermodynamics, collapse, or any physical
mechanism. The actual-history result must remain equally narrow.

Reproduce the prospective freeze with:

```bash
make fgc-tdg1-frz1
```
