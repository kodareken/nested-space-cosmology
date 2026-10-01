# Variational Fourier–Galerkin coupling

## Current execution

`nsc_spherical_galerkin_coupling.evolve` is the state-and-geometry step.
[derive_nsc_spherical_galerkin_refinement_v5.py](../scripts/derive_nsc_spherical_galerkin_refinement_v5.py)
owns the self-contained saved diagnostic and preparation. It is not a
permanent production-driver requirement. The episode driver
`scripts/derive_nsc_spherical_feedback_episode.py` reuses that preparation.
The saved \(T=0.05\) record is
`results/development/nsc-spherical-feedback-episode-v1.json`, with payload
`results/development/nsc-spherical-feedback-episode-v1.npz`. Verdict
`MEASURED_FEEDBACK_UNRESOLVED_CONSTRAINT_CONTROL`. Four runs reach
\(T=0.05\). All 62 physical refinement rows meet one percent. Renewal is
false, and a constraint-consistent solution is not validated. This note's
v1, v2, and v3 records remain the earlier failed controls. They did not
produce that episode.

Each rate evaluation prolongs the state and recomputes the column source
with `source_from_columns`. Geometric rates come from the summation-by-parts
chart. With the matter force on,
\(\dot p_Q\leftarrow\dot p_Q-F_Q/\Delta x_Q\) before the adjoint pullback.
The initial radius uses \(\rho=F_L/\Delta x_Q\) from that source.
`solve_initial_radius` does this when \(\rho\) is omitted. The v5 replay's
`full_source_rho` does it and holds \(\rho\) fixed while Newton updates
\(r\). \(L\) and \(\beta\) stay the sampled calibration gauge. \(Q\), \(r\),
\(\chi\), the momenta, and the columns evolve. That gauge is not frozen
physical geometry: the shift still drags \(Q\), and \(F_Q\dot Q\) is the
coordinate work of that drag.

v5 already measured \(T=0.005\) at \(n_f=512\), \(n_q=2048\), ten steps of
\(dt=5\times 10^{-4}\). The record verdict is `DIAGNOSTIC_MEASURED`.
Renewal is false. `time_extension_T_0_05` is false in this record.
Window maxima of absolute values are Hamilton
`9.981468739539423e-06`, momentum `1.3753934746951746e-06`, chart proper
speed `0.007950515372629119`, and lifted proper speed
`0.007950515381253666`. Those speed fields are max-abs. They are not
signed radial velocities. The initial Hamilton
residual has that same maximum and is held-out
(`9.978646438460075e-06`). The initial momentum residual is
`9.184208948907546e-12`. `measured_within_window_1e-3` is true.
`labeled_as_window_pass` is false.

Full and projected strong residuals both remain. The initial `1e-8` flag
is false and did not stop the run. It is not an automatic veto. A
weak-form error statement still needs the unresolved complement, a
quadrature comparison, roundoff, and observable convergence. The projected
residual is not that statement. At \(n_f=256\), doubling quadrature moves
full Hamilton by `4.177484925094177e-06`. No doubled-quadrature pair is
stored at \(n_f=512\).

The stored radius is the retained best projected iterate. The last
accepted Newton trial is smaller by `1.610240673480347e-09` and is not
the evolved state. Both values sit above the Newton `1e-10` line. The
independent file now has ten tests, sha256
`dde0ca7dcd8d870f772aff3f6069286cddb8493c291bddb715232b493b44f55f`.

The chart sample in `normal_velocities` is
\((\dot r_{\mathrm{unprojected}}-\beta r_x)/(r L)\) from `geometric_rates`,
which includes \(\Pi\). The returned `chart_proper_max` and
`lifted_proper_max` are maxima of absolute values. The chart sample is zero
on the initial slice because \(\Pi=0\). The function comment that the chart
velocity uses \(p=0\) matches that slice and not a later state whose
\(\Pi\) is nonzero. At the end of the \(n_f=512\) window the two max-abs
speeds agree near `0.00795051538`.

The sections below are the v2 control and the v3 initial solve. Their
verdicts and JSON bytes stay historical.

## v2 control

This is a fixed subspace of the existing spherical Hamiltonian. It is not a
filter of the v1 radius, not a constraint projection, and not a renewal.
The checkpointed v1 record
`lab/results/development/nsc-spherical-coupling-control-v1.json` is not
edited. The coupling source and its tests now carry the minus-column phase
correction. This record reads that corrected preparation and does not
rewrite the v1 JSON. No force, state reset, prescribed radius, or physical
constant was added. Signs, occupations, multiplicity `4κ`, and `D(F)` are the existing
action and source owners, evaluated on the quadrature grid.

The packet is not rebuilt here. Physical columns come from
`nsc_spherical_coupling.prepare_rank6`: real envelopes, local carriers
`exp(i k (x-n*ell))`, and the recorded phase on the minus column.
`phase_applied_to_odd_lobe` is false and
`phase_on_minus_spinor_column` is true. The old link `B` is compared and
not copied. `B_phys` is not the old `B`
(maximum link difference `8.950690833082303`). Old circulation is not
transferred.

| Path | sha256 |
|---|---|
| `lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py` | `85d1f3dbbd82a38fe14d2405ad9782af6de68095052ff522ac6602e5d50a14ab` |
| `lab/tests/test_nsc_spherical_galerkin_coupling.py` | `e928b0d9b9d1162e3f665681adf1d4c9a3ff00c2941189d20bdda937c4b51a19` |
| `lab/src/recursive_horizons/nsc_spherical_coupling.py` | `64e066b9b3210121461f902e748d4fc8a9cfc5978af2af47ab4b12fd8797e4f5` |
| `lab/tests/test_nsc_spherical_galerkin_independent.py` at the v2 note (6 tests) | `01afcc16c4dda9447ab7072a87dc50c3c5a74eb5470518793c5b00054a82abbb` |
| `lab/results/development/nsc-spherical-coupling-control-v2.json` | `26030900805ebe27553e1b8d1cd490c6ca1a79c4631b65dad5483db3c27b370c` |

Ten tests in `test_nsc_spherical_galerkin_coupling.py` passed in 0.35 s. The control took 0.820 s wall and 0.843 s CPU.
Verdict `FAIL_HELD_OUT_CONSTRAINT`. This v2 control did not run `T=0.05`. Renewal is false.
The saved record still has `independent_review` false. This note does not
rewrite that JSON.

Independent review run `9be42d26` passed 15 tests. The source symplectic
pullbacks and the lifted work identity are valid. The verdict
`FAIL_HELD_OUT_CONSTRAINT` stands. A projected residual near `1e-10` is
not a physical pass.

## Representation

Geometry uses an odd count `ng = nf - 1`, so the Nyquist degree is absent.
Fermions stay on the even antiperiodic count `nf`. The Hamiltonian is the
existing fine-grid operator on a quadrature grid `nq`. The default is
`nq = 4 nf`; the record also evaluates `nq = 8 nf`.

`A_g` interpolates geometry nodes `ng → nq` and
`(ng/nq) A_g.T A_g = I` (defect `8e-15` at `ng=63`). `A_f` interpolates the
fixed half-integer band. Canonical columns use `U_f = sqrt(nf/nq) A_f`, and
`U_f† U_f = I`. Occupation eigenvalues of the six-mode Gaussian are
unchanged. Column rates are `U_f† (-i H_fine U_f Φ)`, not `-M i H Φ`.
Geometry rates, for both coordinates and momentum densities, are
`(ng/nq) A_g.T` times the fine rates. `L` and `β` are the same static gauge
functions sampled on the quadrature nodes. Source weights include `4κ` once.
Only densities are divided by the quadrature spacing. Work uses the lifted
coarse `Q̇ = A_g Q̇_g`. On the physical slice that identity holds to
`5.7e-10`; the unprojected fine rate misses it by `1.7e-5`.

The initial radius is a fresh Newton solve inside the geometry subspace,
with the actual Gaussian, all momenta zero, and `χ = 0`. It is not a
filtered v1 profile. `r` and `Q` stay positive on the quadrature grid
(`r` about `4.210` to `4.877` at `nf=128`). The chart proper normal velocity
is `0` because the momenta vanish. The lifted proper velocity, which uses
the actual coarse `ṙ`, is `3.52e-4` at `nf=64` and `1.22e-6` at `nf=128`.

Directional derivatives of the fine energy against the pulled-back
right-hand side have relative error `1.65e-10`.

## Alias test

Low-band geometry, chart `F`, and matter-symbol gaps agree at
`2.06e-13`. A pure mode-31 square is a different statement. On 64
collocation nodes its product-rule gap has maximum `8π = 25.132741228718498`
and Euclidean norm over `2π` equal to `22.62741699796956`. Those nodes are
identical to the folded partner mode `-33`. On 128 nodes the same true mode
has gap `9.6e-13`. The Galerkin quadrature (`ng=63`, `nq=256`) evaluates
that true mode, not the folded samples, and the gap is `2.6e-12`.

## Constraints, by axis

Projected and full residuals are both reported. A projected residual near
`1e-10` is not a pass. Tolerances are unchanged: initial `1e-8`, window
`1e-3`.

| Slice | Projected Hamilton | Full Hamilton | Held-out Hamilton | Full momentum | Doubled-quadrature full Hamilton |
|---|---:|---:|---:|---:|---:|
| `nf=64`, `ng=63`, `nq=256` | `4.77e-11` | `2.375108194941006` | `2.375108194930269` | `5.4e-14` | `2.446258099704416` at `nq=512` |
| `nf=128`, `ng=127`, `nq=512` | `8.47e-11` | `0.034135561038210` | `0.034135561030981` | `1.5e-13` | `0.034135686719168` at `nq=1024` |

The held-out fraction of the Hamilton residual is the whole residual: the
subspace piece was solved, and the pointwise failure sits outside
`|mode| ≤ (ng-1)/2`. Doubling `nq` at `nf=128` changes the full Hamilton
maximum by about `1.3e-7`. Quadrature is not the blocker. Space refinement
cuts it from `2.375` to `0.0341`, which is still above `1e-3`.

Subspace frequencies are `211.583` and `424.915`. Quadrature Nyquist
frequencies are `100.531` and `201.062` and are not used as the timestep.
`dt = 0.0005`, and the halved step is a separate axis.

| Run | Full Hamilton | Projected Hamilton | Full momentum | Energy drift | Unitarity |
|---|---:|---:|---:|---:|---:|
| `nf=64`, `T=0.005` | `2.375108194940754` | `0.074127075589877` | `0.193282408341238` | `2.56e-9` | `5.25e-11` |
| `nf=64`, `dt/2` | `2.375275703620485` | `0.074127086411347` | `0.193282413851183` | `9.3e-11` | `1.6e-12` |
| `nf=128`, `T=0.005` | `0.034267189761761` | `0.010409885819528` | `0.006228653251853` | `3.55e-9` | `7.16e-11` |
| `nf=128`, `dt/2` | `0.034271614544118` | `0.010409906456804` | `0.006228665766706` | `1.26e-10` | `2.2e-12` |

No positive-chart stop. `r` stays near `4.21`. Halving `dt` does not move
the constraint maxima. The held-out Hamilton tail is already present at
`t=0` and does not grow. The projected Hamilton piece and the momentum
residual do grow, and both exceed `1e-3` inside the short window.

Recorded v1, which used the old odd-lobe packet and collocation, reached
Hamilton `36.98886722265645` at `N=64` and `2.330762882869135` at `N=128`
by `T=0.005`. Those numbers are the saved v1 diagnosis, not a rerun. This
window is smaller (`2.375` and `0.0341`) and does not show the old
shift-driven jump from a solved collocation residual. It still fails
`1e-3`. The packet correction and the subspace are both different from v1,
so the drop is not a single-cause claim.

## Blocker

The physical constraint is the full quadrature residual, including modes
outside the odd subspace. That residual is `2.375108194941006` at
`ng=63` and `0.034135561038210` at `ng=127`. Doubled quadrature does not
remove it. This v2 control did not start `T=0.05`. A short step with a tiny projected
residual is not a renewal.

## Matched refinement v3

This successor asks whether matched refinement `ng = nf - 1`, `nq = 4 nf`
brings the full quadrature constraints under the declared initial and
window tolerances. It does not. The driver and record below are immutable
failed evidence. They do not rewrite the module or the v1 and v2 JSON
files. `module_unchanged`, `v1_bytes_preserved`, and `v2_bytes_preserved`
are true. Renewal is false. This v3 record did not run `T=0.05`.

| Path | sha256 |
|---|---|
| `lab/scripts/derive_nsc_spherical_galerkin_refinement.py` | `2033c6ec99143ccc04666717533c8f67748570f403e7d8228fda5e7e80936b3c` |
| `lab/results/development/nsc-spherical-coupling-refinement-v3.json` | `5f324e51f87ccc54781e6d4142cb394fb4f32a5bc4505b8791900dcc5de27a97` |

At `nf=256`, `ng=255`, `nq=1024`, internal Newton converged at scales
`0.125`, `0.25`, and `0.375`. It stopped at scale `0.5` with residual
`1.2876060991167562e-10`, above the Newton tolerance `1e-10`. The recorded
blocker is that Newton left the positive quadrature chart or stalled
before full rho. The verdict is `FAIL_INITIAL_SOLVE`. The record's
top-level `status` remains `PROVISIONAL`. Full Hamilton, momentum,
held-out, and projected diagnostics are null. No full-source result was
produced. The saved v1 sha256 remains
`88bea1c962478a2ec57e0c0f3844ab6ee9237d148dc1863f10122ee8fa5a9b60` and the
saved v2 sha256 remains
`26030900805ebe27553e1b8d1cd490c6ca1a79c4631b65dad5483db3c27b370c`.
