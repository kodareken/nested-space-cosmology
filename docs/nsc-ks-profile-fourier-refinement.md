# Sharper bounds for the same axial profiles

The initial p=8 profile table left a first-derivative tail of about
`2.76e-8` for w. That conservative geometry error is too large to allocate
across the complete source inventory. This refinement retains the same
history, period and profile law, and tightens the existing rigorous
integration-by-parts/alias bounds.

Use derivative order 16, 1024 cells per transition, 65536 ball-arithmetic
samples and retained Fourier index 8192. The larger retained band is a
numerical representation with a bounded remainder, not a physical cutoff.
The prior coefficient intervals remain a compatibility control; the new
enclosures must overlap them on the shared indices. No source or Dirac
evolution is repeated. The preparation is capped at 120 CPU seconds.

The exact dyadic balls and remainder bounds are authoritative; binary64
numbers in the compact report are display approximations. This is a
geometric axial-profile certificate, not the source-energy UV certificate
or a completed field-error budget. The physical local incoming gate is OPEN.
