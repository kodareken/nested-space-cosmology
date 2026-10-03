# Parent episode adapter

This note owns the checkpointed parent episode around the existing leading
carrier. It does not replace the leading rates, the analytic Jacobian action,
or RK4, and it does not reopen a cosmological fit, an infinite-time claim, or
a bounce search.

The public reconstruction entry point is
`reconstruct_parent_pair(arrays, record)` in
[the adapter](../src/recursive_horizons/nsc_discovery_parent_episode.py).
The [driver](../scripts/derive_nsc_discovery_parent_episode.py) is a pure plan
unless prepare, run, or check is requested. Focused tests live in
[the test module](../tests/test_nsc_discovery_parent_episode.py).

## State and regions

The carried state is \(Q,r,p_Q,p_r,\phi_0,\phi_1\). Canonical momenta are
\(\pi=\Delta x_g W^{\mathsf T}p\). \(W\) may be the identity. Either way it is
a full \(n_g\times n_g\) frame, so the geometry band stays active. The fermion
columns stay on the full antiperiodic band. Their rank is the length of the
stored frame weights, not a fixed column count. The three populations are
weight contents of one source geometry: `child_heavy` (population 0,
\(+\alpha\)), `balanced` (population 1, no extra imbalance), and `parent_heavy`
(population 2, \(-\alpha\)). They are not three different source geometries.

The centre is \(L/2\) (so \(4\) when \(L=8\)). \(s\) is the signed distance from
that centre. The child is \(|s|\le 0.5\), the protected collar is \(|s|\le 1\),
the parent is \(|s|\le 3\), and the parent annulus is \(1.2\le|s|\le 3\).
Quadrature occupations are the stored weights, installed by replacing the fine
system dataclass. The historical six-column pair constructor is not used.

## Frozen parent record

Production input is one `NSC-DISCOVERY-PARENT-v1` directory per
\((n_f,\mathrm{population},\mathrm{sign})\), written by the parent producer as
`parent.json` and `parent.npz`. The episode reads that record; it does not
invent a second preparation schema.

`prepare_parent(nf, population, sign, k_override, profile, cpu_limit)` returns
`(pair, leading.State, report)`. `population` is `0`, `1`, or `2`. `sign` is
the geometric momentum sign `+1` or `-1`. The accepted NPZ arrays are `Q`,
`r`, `pi_Q`, `pi_r`, `phi0`, `phi1`, `W`, `weights`, `source_columns`, and
`reference_columns`. JSON carries `schema`, `nf`, `population`, `sign`, `k`,
`k_common`, `intervals`, `clock_locations`, `payload_sha256`, and
`array_sha256`, as well as `input_hashes`, `producers`, and the immutable
`producing_commit`. The native loader authenticates those pins before the
episode adopts any arrays. Manufactured loader-boundary tests are not frozen
scientific records.

Phi is copied byte for byte. A record whose momenta were already solved at
`sign` is not signed again. If a test field is only the positive solution,
the minus case negates both `pi_Q` and `pi_r`. Conjugation is not a sign and
is not an exchange reversal.

The declared offset comes from the authenticated parent producer's `k`, carried
as `parent_k` on episode checkpoints. These authoritative values precede legacy
`common_k`/`k_common` labels, and a null label never masks an available offset.
The six magnetic NF128 originals have chosen `k=0.002621213673327129` while the
old suggested `k_common` was about1e-8. Old episode labels remain untouched;
read-only checking reports the mismatch and binds the canonical initial state
to the original authenticated parent record. New headers retain the true offset.

The separately measured kinetic mean is `mean(p_Q^2/r^3)` on the actual lifted
state. It agrees with the chosen offset in the uniform magnetic branch; a
general collar contains `J` and can have a different mean. No equality is
imposed outside the uniform branch. Source-free controls use the original
chosen offset, not the obsolete suggested minimum.

The primary episode cases are \(n_f=128\) with step cap \(0.001\). Confirmation
is \(n_f=256\) with step cap \(0.0005\). The active first batch records stations
\(1,3\) and runs each case uninterrupted through both before physical review.
The historical default remains \(1,3,8,16,24\). A stop
before the next station is kept as an event checkpoint.

## Numerics

The timestep delegates to `nsc_discovery_parent_step_control.step_restriction`
and passes `control_mode`. The returned record is pinned to the live source
weights and source columns. The episode does not keep a second `4+4R` probe.
The older raw coordinate norm is not installed. Native FFT carriers are built
before the actual occupation width replaces the factory weights, then limited
to one thread. Hooks live only inside the scoped adapter and the spawn
initializer, and both paths restore the previous callables.

Every RK4 stage evaluates the source from the stage state. Frozen geometry
keeps the same columns and sets the geometry jets identically to zero.
Source-free control empties the source on a copy, keeps the weights and the
source metadata, and solves both momenta against the new constraints. Zero
weights are not treated as a solved constraint. The caller's arrays are not
written. The source-free copy retains its column slots and weights but its
actual covariance rank is zero; its column Gram distance from the identity
is one. It retains the original `pair.common_k` when that offset is admissible
for the empty-source square root. An inadmissible offset is explicitly OPEN
and refused; the control does not select a different motion. The record stores
the actual control offset, kinetic anchor, momentum correction and own C/D.
Frozen geometry is a prescribed control and does not demonstrate autonomy.
Resume reloads the stored state, normal clocks, and work-ledger
stocks. It does not project or replace the source.

Scalar rows use the leading rates, constraints, curvature jets, clocks, and
ledger at the stored rank. They do not call a fixed six-column Gram identity.
They add the source density, source current, and clock samples. When
`nsc_discovery_parent_observer.observe(pair, state, time, control_mode="coupled")`
is importable, its result is attached and a hook refusal is recorded on the
row. No held-out population or gradient measurement is evaluated before a
later prediction.

Every scalar row also stores weighted probability in child, parent, protected
collar, parent annulus and the full carrier for each original column index.
The tags use evolved columns with their original weights; no observer or
source frame is rebased. They describe source-column ancestry, not independent
particle energies. Explicit stations `1,1.25,2.25,3` retain full-state chunks at
the two additional ancestry checkpoints as well as1 and3.

## Budget and commands

One process pool runs at most six one-thread cases. The aggregate CPU budget,
including child time, is \(21600\) seconds. The forecast factor is \(1.5\), the
resident ceiling is \(8\) GiB, and each immutable chunk stays within \(64\) MiB.
Hitting the budget stops the case and keeps the last checkpoint. Only `run`
launches that pool. The root executor freezes the producer commit before
production; preparation rejects a source closure that is not that commit.
Analytic fixtures are for tests and cannot be prepared as a campaign.

~~~sh
python scripts/lab.py scripts/derive_nsc_discovery_parent_episode.py
python scripts/lab.py scripts/derive_nsc_discovery_parent_episode.py --prepare --stations 1,3 --source <owner1-output> --producer-commit <frozen> --cpu-budget <remaining-global-CPU>
python scripts/lab.py scripts/derive_nsc_discovery_parent_episode.py --prepare --confirm --source <owner1-output> --producer-commit <frozen>
python scripts/lab.py scripts/derive_nsc_discovery_parent_episode.py --run --workers 6 --cpu-budget <remaining-global-CPU>
python scripts/lab.py scripts/derive_nsc_discovery_parent_episode.py --check
python scripts/lab.py -m pytest tests/test_nsc_discovery_parent_episode.py -q
~~~

The default command prints the plan and writes nothing. No regeneration proof,
stability certificate, or held-out prediction is claimed.

`prepare(..., stations=(1, 3))` and CLI `--stations 1,3` record the same positive,
finite, strictly increasing targets in the manifest, every case and its
numerical binding. The last stored target controls the admission cost forecast;
the active first batch forecasts through3, not24. Run and check read the frozen
targets and reject a CLI station override. Source preparation and authentication
are unchanged; there is no global station mutation.

Preparation/preview also accepts `control_mode="coupled"|"frozen_geometry"|
"source_free"` and an optional positive finite `step_cap`. CLI forms are
`--control-mode frozen_geometry --step-cap 0.0005` or `--control-mode source_free`.
Each output remains creation-only; root supplies the selected one- or two-record
source directory. Run/check consume the recorded control and cap. Future run
metadata reports whether a pool actually launched and counts the final manifest
rewrite bytes explicitly; worker chunks and observation streams are excluded
from that byte counter. Older records are not rewritten to add these fields.

When the live source closure differs, read-only check authenticates the original
producer at its frozen commit and qualifies the result as historical source and
stored-array structural authentication. It does not replay current dynamics or
claim current producer compatibility.

The generic station-resume helper remains separate. It copies state, clocks and
work stocks but creates a fresh ledger and omits the parent manifest's top-level
source closure. A later root-owned resume must preserve that binding and deduct
measured prior preparation/episode CPU from the global allowance; a new output
directory does not renew the global budget.
