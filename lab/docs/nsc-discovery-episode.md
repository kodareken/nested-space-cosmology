# Discovery episode: stage 0 handoff and stage 2 continuation

This note owns the first continuation batch after the saved separated
parent–child pair. It does not replace that pair, the frozen replay
basis, or any sealed trajectory. No batch was launched while this
interface was frozen.

The saved inputs are the matched \(dt=0.0005\) baseline states at
\(T=0.3\) for \(n_f=128\) and \(n_f=256\), in
`results/development/nsc-nested-parent-child-v1.json` and its payload,
together with the frozen replay basis. The continuation reads those
final arrays and the stored \(W\). It does not call `initial_state`,
does not repeat the radius solve, and does not rebuild \(W\). The dense
grid is only the initializer. Evolution requests carrier backend `auto`.
In this environment that recorded 9.55× FFT comparison selects FFT.

## What is continued

The first batch has six cases, all taken from the same \(dt=0.0005\)
state at \(T=0.3\):

| Case | Step cap | Controller |
|---|---:|---|
| `nf128_coupled_dt0.001` | 0.001 | coupled |
| `nf128_coupled_dt0.0005` | 0.0005 | coupled |
| `nf256_coupled_dt0.001` | 0.001 | coupled |
| `nf256_coupled_dt0.0005` | 0.0005 | coupled |
| `nf128_frozen_geometry_dt0.0005` | 0.0005 | frozen geometry |
| `nf256_frozen_geometry_dt0.0005` | 0.0005 | frozen geometry |

The default stations are \(T=1\) and \(T=3\). Stations \(8\), \(16\), and
\(24\) remain accepted when requested. There is no extra maximum-time
stop and no zero-residual gate.

The frozen controller holds \(Q,r,\chi\) and their canonical momenta at
the handoff values and steps only the carrier columns. Its observation
and stability row receive that actual rate with
`control_mode="frozen_geometry"`. \(\dot Q\), \(\ddot Q\), \(\dot L\),
and the metric coordinate work are zero. Field rates, source forces, and
spatial lapse work stay. That row is not a coupled stability classification.

`nsc_discovery_observables.observe(pair, state, time, *, nodal_rate, control_mode, ...)`
receives those names when its current signature has them. Older observers
are still called positionally. A resolver that changes \(W\) or the saved
state bytes is not accepted.

## Output

New immutable chunks may be created under
`lab/results/development/nsc-discovery-...`, and tests may use a temporary
directory. An existing sealed file is never overwritten. A later commit
uses a new chunk name. `--check` reloads hashes and source pins and does
not evolve. Source-column, observer, weight, and \(W\) pins are the same
before and after a step. The state \(\Phi\) is not stored as those source
columns, and clock rates are not stored as the source weights.

## Step restriction

The field frequency remains `subspace_frequency`. It is not used alone.
The geometric restriction comes from the owned conformal principal part.
With the coefficient identity \(Z=3a\), that second-order symbol is
\(\partial_{tt}=\partial_{xx}\) on \((Q,r,\chi)\), so the identity speed is
the coordinate speed \(1\), widened by the light-cone speed \(|L|/Q+|\beta|\).
The wavenumber is the largest nonzero symbol of the owned periodic
derivative, whose Nyquist entry is zero. The step is

\[
\Delta t=\min\left(\Delta t_{\mathrm{cap}},\frac{1.4}{\max(\omega_{\mathrm{field}},\omega_{\mathrm{geometric}},1)}\right).
\]

If \(Z=3a\) is not present to rounding tolerance, the step uses the
absolute majorant of the uncancelled symbol instead. This is a sufficient
frequency restriction. It is not a stability certificate:
`stability_certificate` stays false. Lower-order geometric terms, time
discretization, and projection defects are outside this bound.

## Diagnostics

These quantities are stored separately. They are not collapsed into one
score and none of them vetoes the continuation:

| Record | Content |
|---|---|
| Actual metric curvature | \(R_h\) and Weyl \(C^2=(R_h-2)^2/(3r^4)\) from the metric jets |
| Chi-shell | \(\chi^2/(3r^4)\) and \(R_h-\chi-2\), not substituted for \(R_h\) |
| Projected constraints | Projected Hamilton and momentum residuals, with held-out residuals beside them |
| CAR | Preserved carrier frame: column Gram and occupation-weighted eigenvalues |
| Constraint-algebra residual | Sampled defect of \(\partial_t h_c=\partial_x D\) and \(\partial_t D=\partial_x h_c\) |
| Energy, work, chart | Total, field, and gravity energy; fieldwork power; chart admissibility |

## Chart exit

`PositiveChartExit` ends the case. The committed state is the last state
that was still inside the chart. The event bracket is the step that left,
\([t,t+\Delta t]\). The failed values are not stored, not clamped back into
the chart, and not reported as a physical instability.

## Checkpoints and budget

A commit is an immutable `.npz` plus its `.json`, each written once and
then marked read-only. A chunk larger than 64 MiB is rejected. The arrays
are the native fields \(Q,r,\chi\) and canonical \(\pi\), or an explicit
nodal-momentum representation, plus \(\Phi\), source weights, \(W\),
clocks, the case, the step count, and source pins. Reloading those bytes
is the continuation. `--check` verifies the hashes, the read-only mode,
and the pins. It does not evolve and it does not rewrite a chunk.

Root launches the batch. These commands do not run it from this note:

```sh
python scripts/lab.py scripts/derive_nsc_discovery_episode.py --prepare --output results/development/nsc-discovery-episode-v1
python scripts/lab.py scripts/derive_nsc_discovery_episode.py --check --output results/development/nsc-discovery-episode-v1
python scripts/lab.py scripts/derive_nsc_discovery_episode.py --run --output results/development/nsc-discovery-episode-v1
```

`lab.py` uses `lab/` as the working directory, so that output path is the
new successor `lab/results/development/nsc-discovery-episode-v1`. `--run`
uses one process pool. The defaults are backend
`auto`, 6 workers, 21600 seconds of summed child-process CPU, forecast
factor 1.5, and an 8 GiB resident ceiling. A case whose forecast does not
fit, or whose resident peak passes that ceiling, stops on the last
admissible commit. `--dense` or `--fft` forces the carrier. The saved
state is still the dense initializer, not a new solve.

## Samples, snapshots, and quadrature

Scalar rows are stored in `<case>-observations.jsonl`, separate from the
state chunks. The default cadence is \(0.05\) in coordinate time. A row is
also written at each station and at the last admissible event. Full states
are committed at the handoff, at each station, and at that event, so a run
that reaches \(T=3\) keeps the \(T=1\) state. Chunks stay at most 64 MiB.
No dense propagator history is stored.

Clock integrals use the trapezoid of the normal rates \(rQ\) at \(x=1,2,3\)
on every accepted step:

\[
\Delta \tau = \tfrac12 \Delta t\,(\dot\tau_- + \dot\tau_+).
\]

Fieldwork, pressure, lapse, and the child and parent boundary fluxes are
spatial integrals of the actual control rate. Frozen control samples have
zero coordinate fieldwork. Their time integrals use the same trapezoid
between successive samples, not a second coupled reconstruction. Restart
continues the stored sample and the stored integral.

Read these records with `read_observations`, `read_station_states`,
`read_case_ledger`, and `read_physical_audit`. The audit reports the
observation and step process times.

## Not claimed

The v1 record's own status is unchanged. Reaching a station is not a
continuum statement, not a certificate that the step restriction is sharp,
and not evidence that a chart exit is a physical instability. This note
does not authorize root's later launch by itself.
