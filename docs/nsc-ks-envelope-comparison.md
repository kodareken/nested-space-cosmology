# Independent KS-envelope comparison with the stored PG state

The PG CF4 comparison stopped at its declared three cases with a local
120/240-step matter difference about2.3e-9. This experiment uses the same
continuum Dirac law in canonical KS envelopes. It reuses the stored PG
fields and never reruns a PG case or a horizon/scattering/source generator.

Initial canonical amplitudes are the fixed7601 reference columns at the
first saved node at or above rho1.03, transformed by the owned
`archived_amplitudes`. Actual rho is1.0300000000000002. The Gram defect is
about4.92e-13 and is retained, not normalized away. The source covariance,
four real signed energies, weights and multiplicity12 are unchanged.
The same alpha=.001,T1*plateau,U=0 history is integrated to rho1.

The computational envelope period is .4 with support padding .14 on each
side. The continuum characteristic length is about.0432803, leaving a
conservative double-travel margin about.0534394. This is numerical domain
padding, not a physical periodic condition on the incoming constraints.

The bounded cases are reference64; nonzero64/128; nonzero256 only if the
64/128 local matter change exceeds1e-11; then a tighter DOP853 comparison
at the finest completed grid. Maximum5 cases and120 CPU seconds. Base
rtol/atol/max_step are2e-12/2e-14/.001 and tight options2e-13/2e-15/.0005.
Targets are the21 nodes of the approved I=S(1)+[.12,.18].

Reference drift is tested against3e-11. Spatial and temporal differences
are indicators with target1e-11. Agreement with the stored PG240 matter
uses5e-9, reflecting that path's still-open2.3e-9 time indicator; it is not
a rigorous combined error bound. Field, axial and retarded-tangent
differences are reported separately. No comparison changes the state law,
and no measured reference drift is subtracted.

Full source, continuum, spectral truncation and between-node errors remain
OPEN. This is not a root or a physical local-gate certificate.

```sh
python3 scripts/derive_nsc_ks_envelope_comparison.py --run
python3 scripts/derive_nsc_ks_envelope_comparison.py --check
```
