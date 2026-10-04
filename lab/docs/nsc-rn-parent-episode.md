# RN radial parent episode

Observers, the boundary ledger, and the campaign CLI live here. Vacuum
geometry, horizons, SBP weights, the weighted Gram and curvature invariants
belong to
[`nsc_rn_reference.py`](../src/recursive_horizons/nsc_rn_reference.py).
Source-free PG reconstruction, rates, metric jets and RK4 belong to
`nsc_rn_pg`. Neutral preparation belongs to `nsc_rn_source`. This owner
does not rebuild those formulas and does not turn a first-derivative
product \(D@D\) into a metric jet.

The driver is
[`derive_nsc_rn_parent.py`](../scripts/derive_nsc_rn_parent.py).
Diagnostics are
[`nsc_rn_observables.py`](../src/recursive_horizons/nsc_rn_observables.py).
Checks are
[`test_nsc_rn_observables.py`](../tests/test_nsc_rn_observables.py) and
[`test_nsc_rn_parent_driver.py`](../tests/test_nsc_rn_parent_driver.py).

## Magnetic family

\(r_m=\sqrt{P^2}\) is the magnetic radius. The initial reference sets
\(P^2=r_m^2\) with fixed magnetic \(q=1\) and \(V_4=0\). It is not
\(r_-\). The mass ratios are

\[
M/r_m\in\{1.002,\,1.01,\,1.04\}.
\]

Horizons are the reference roots of \(r^2-2Mr+P^2=0\). Excision is the
midpoint \((r_-+r_+)/2\), where both radial null speeds leave the domain.
The child collar is \((r_++\Delta/16,\,r_++\Delta/2)\) with
\(\Delta=r_+-r_-\). The parent is the disjoint exterior from that collar
to the member outer boundary. The default outer boundary is \(32 r_m\) and
the control outer boundary is \(64 r_m\). Point floors are 32 and 64. The
spacing is the physical throat rule

\[
\min(\lambda_{\min}/24,\,\Delta/16,\,\sigma/24),
\]

using only the scales the source owner has supplied. A \(D(2,1)\) product
\(D@D\) is not that grid: its endpoint second derivative is half the
consistent jet, and a sample with that bias is rejected rather than doubled.

## Matter and the ledger

The shell is closed and neutral: \(\kappa=1\), multiplicity \(4\) once, no
filled sea, no electric current. Canonical \(\phi\) has shape `(2, N, rank)`.
The Gram and the Hamiltonian sample use the reference SBP weights. The
angular CAR inserts \(\kappa/r\) once and does not multiply by the
multiplicity again. Ledger accounts stay separate: probability remaining,
energy remaining, outer outflow, excision outflow, SAT debit, and the
moving-horizon Reynolds term \(\rho(r_h)\dot r_h\).

A station is compared with the initial reference RN and with one
enclosed-mass RN built by `RNReference.from_action`. The sample lapse is
not replaced by 1. There is no per-radius refit. On a stationary source-free
slice the event horizon is the trapping horizon. On a dynamical slice it
stays unresolved.

## Calibration and admission

The default CLI is a read-only preview. `--calibrate` calls, when present,

- `nsc_rn_pg.sourcefree_rates`
- `nsc_rn_pg.metric_jets` on that reconstruction
- `nsc_rn_pg.rk4_step` for one step
- `nsc_rn_pg.build_grid` with the magnetic radius and the route outer boundary

and compares them with the reference jets, charged mass, the trapping
horizon and the curvature scalars. \(R_4\) and Ricci² are computed from
the numerical lapse, shift and their radial jets, then compared with the
independent reference values. A non-negligible \(\partial_t N\) or
\(\partial_t\beta\) leaves the curvature unresolved rather than inserting
the static formula. The tolerance is \(10^{-4}\). Injected analytic jets
are not a calibration. A pass boolean, including `root_scientific_calibration_passed`,
does not grant admission. Admission reads the measured residuals and the
producer-hash binding. Root review is a separate decision and does not
override a failed measurement.

The route spacing is \(\min(\lambda_{\min}/24,\Delta/16,\sigma/24)\).
\(\sigma=2r_m\) and \(\lambda_{\min}=2\pi\min(|a_+|,|a_-|)/E_{\max}\) come
from the positive-frequency window at \(12r_m\). The march step is
\(\min(\Delta t_{\mathrm{CFL}}/2,\,0.001\,r_m)\). On the first pair the cap
is the smaller one: half the full CFL is about \(0.0087\,r_m\) and is not
the step that the forecast counts.

The first pilot uses nine positive-frequency samples. Rank 5 on the
869-point route splits into lobes near \(12r_m\) and \(30r_m\). Rank 9 has
one dominant lobe. The measured RMS width is kept as measured, about
\(2.94r_m\) on that grid, and is not refitted to the request \(2r_m\).

The matched pair shares one prepared column and the metric reconstructed
from that column at \(t=0\). The coupled arm calls `rk4_step`, so the source
is rebuilt at every stage. The fixed arm is that same initial lapse and
shift held frozen, advanced only by `dirac_rates`. It is not a second
coupled evolution and it is not the vacuum RN. Probability stocks, mass-flux
stocks and the SAT debit stay separate. A sourced station whose lapse or
shift rate exceeds \(10^{-5}\) stores null \(R_4\) and Ricci². The missing
primitive is the radial derivative of those time rates, \(N_{tr}\) and
\(\beta_{tr}\).

`forecast_first_pair` times one coupled step and one frozen Dirac step and
projects each station with the factor 1.5. It does not march and does not
write. `--run` performs the march only after the measurement binding and a
frozen commit:

```
python scripts/lab.py scripts/derive_nsc_rn_parent.py --run \
  --output <unsealed-directory> \
  --producer-commit <40-hex> \
  --calibration-record <calibration.json> \
  --mass-over-rm 1.04 --energy-fraction 0.001 \
  --station-limit 5 --cpu-seconds 21600
```

`--station-limit` accepts 5, then 10, 20 or 40. A CPU stop leaves the last
completed step in memory; a reached station is a creation-only checkpoint
under 64 MiB. `resume` reloads that checkpoint and does not evolve.
This follow-up does not write a scientific record.

A positive frequency passed into the PG state is numeric metadata only.
The reference does not provide `require_exterior_killing_frequency`, and
the calibrator does not call it or treat the scalar as a Killing-mode
proof. `--calibrate` stores the measured source-free residuals. It does
not substitute analytic jets.
