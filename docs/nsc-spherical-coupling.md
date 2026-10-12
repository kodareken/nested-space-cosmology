# Provisional spherical coupling, v1 constraint diagnosis

This note records the failed v1 control and one diagnostic pass. It is not a
renewal, not a constraint pass, and not a change of the equations. The v1
source and the v1 JSON were not edited. No new result record was written.

| Path | sha256 at this diagnosis |
|---|---|
| `lab/src/recursive_horizons/nsc_spherical_coupling.py` | `8727a46aa7785a326e301eeeee370a77aeaa066ddb1dbaee022a255616887dbf` |
| `lab/tests/test_nsc_spherical_coupling.py` | `379fc239a1a3279a842db48e92d41c5be326ac85f616c773f420cd94ffa28083` |
| `lab/results/development/nsc-spherical-coupling-control-v1.json` | `88bea1c962478a2ec57e0c0f3844ab6ee9237d148dc1863f10122ee8fa5a9b60` |

The saved verdict is `PROVISIONAL_CONSTRAINT_DRIFT`. Seven tests of that
source passed in 0.45 s. The control write took 0.044 s wall and 0.043 s CPU.

## What v1 actually integrated

The state is six columns with occupations `(0.75, 0.75, 0.5, 0.5, 0.25, 0.25)`.
Empty complement occupation is 0 and is not a vacuum. The column equation is
`Phi_dot = -i H Phi`, not `-M i H Phi`, with `M = 4` once. Geometry `Q`, `r`,
`chi` and their momenta evolve from the discrete summation-by-parts
Hamiltonian. `L` and `beta` stay the static gauge.

The initial radius is the Newton solution of the discrete Hamilton constraint
with `rho = F_L/dx` from that state. It is not a prescribed profile. At N=64,
`r` runs from `3.816474920558259` to `4.7962017918006286`. The magnetic
bracket minimum is `0.0631642220827372`. `r_mag^2 = 0.9999999999999984`.

On that slice the field energy and the gravitational energy cancel:

| N | Field energy | Gravitational energy | Total |
|---:|---:|---:|---:|
| 64 | `24.7195086` | `-24.7195086` | `-1.4e-14` |
| 128 | `24.7498820` | `-24.7498820` | `8.2e-14` |
| 256 | `24.7499999` | (cancels the field energy) | `-2.1e-14` |

Initial N=64 residuals, zero momenta and `chi = 0`:

| Quantity | Value |
|---|---:|
| Hamilton residual max | `9.556799795973347e-12` |
| Momentum residual max | `6.661338147750939e-15` |
| Shift current max | `6.661338147750939e-15` |
| Proper normal radial velocity max | `0` |
| Coordinate radial velocity max | `0.1864899169029376` |

The coordinate velocity is the shift piece. Directional Hamiltonian
finite-difference relative error on a nonzero-momentum slice is
`3.7537523222797665e-10`. Dropping the matter force changes the power from
`8.88e-9` to `-7.971510461501907`, which matches the fieldwork
`-7.9715104692307985`.

## Short validation, not a renewal

`T = 0.05` was not run. The v1 driver also did not run N=128, because N=64
had already failed its tolerance. This diagnosis ran the same preparation
farther. `dt` stayed at the v1 cap `0.0005`. Estimated frequencies were
`211.58`, `424.90`, and `851.53`, so `omega dt` stayed below the RK4
stability limit. No positive-chart stop occurred.

| Grid | Hamilton max | Momentum max | Energy drift | Unitarity | Number drift |
|---:|---:|---:|---:|---:|---:|
| 64 | `36.98886722265645` | `0.293922324095206` | `1.053e-8` | `5.25e-11` | `1.61e-11` |
| 128 | `2.330762882869135` | `0.013190147175404` | `3.28e-9` | `7.16e-11` | `2.18e-11` |
| 256 | `0.575679312602035` | `0.001291419430788` | `3.61e-9` | `7.39e-11` | `2.25e-11` |

Halving `dt` at N=64 does not change the constraint maxima
(`36.988867326189634` and `0.2939223454447635`). The growth is the
semi-discrete trajectory. Fieldwork on N=64 is `-0.040021316927000936` and
the geometric counterwork is the opposite. Final N=64 `chi` max is
`1.0701968477796657`. Its neighbor product is `-0.1658934550814095`, so
adjacent samples often have opposite sign. The first saved samples are
`0.0087, 0.0054, -0.0121, 0.0330, -0.0398, 0.0669, ...`.

N=256 needs one qualification. Unmodified `solve_initial_radius` stalled at
source scale `0.75` with residual `1.353689282533704e-10`, just above its
`1e-10` line-search test, and returned `converged: false`. The N=256 row
above uses the same residual operator, started from the nested N=128 radius,
and reached `2.6204904912674465e-10` on the full source. It is not a second
equation and it is not the unmodified driver return.

## Operators

Periodic `D` and antiperiodic `P` have the expected symbols.

| Check | N=64 | N=128 |
|---|---:|---:|
| `D sin(2πx/L)` max error | `3.8e-14` | `1.9e-13` |
| `D sin(4πx/L)` max error | `6.6e-14` | `2.0e-13` |
| `D` on the checkerboard | `2.2e-13` | `5.9e-13` |
| `D + Dᵀ` | `0` | `0` |
| `P` eigenvalue at `k = π/L` | error `6e-16` | error `1e-15` |
| `P` eigenvalue at `k = 5π/L` | error `2e-16` | error `5e-16` |

An integer wave is not an antiperiodic eigenvector. That is the boundary
condition, not a wrong sign. Nyquist differentiation is zero because the
real antisymmetric symbol sets that mode to zero. The saved `chi` is not
that single mode: at N=64 only `0.59%` of its power is the Nyquist bin and
`0.42%` has `|mode index| ≤ 4`. The rest is broadband high mode content.
At N=128, `89%` of `chi` power is in `|mode index| ≤ 4` and the neighbor
product is no longer negative (`0.00107` against mean absolute value
`0.024`).

## Initial constraint derivatives

Centered differences along the actual right-hand side, at three steps
`1e-4`, `1e-5`, and `1e-6`, agree. At zero time `chi_dot = 0`. The Hamilton
residual still has a large derivative because `p_Q` and `p_r` are already
moving under the shift.

| Grid | `max \|dH/dt\|` | its `\|k\|≤4` part | `max \|dM/dt\|` | its `\|k\|≤4` part |
|---:|---:|---:|---:|---:|
| 64 | `7345.604` | `36.700` | `46.067` | `27.626` |
| 128 | `459.754` | `0.01491` | `2.603` | `0.01124` |

`7345.604 × 0.005 = 36.728`, matching the N=64 trajectory maximum. The same
identification holds at N=128 (`459.754 × 0.005 = 2.299`). Setting `beta = 0`
and keeping the same `r`, `Q`, `chi`, momenta, and columns drops
`max |dH/dt|` to `1.8e-10` at N=64 and `3.6e-10` at N=128. The Hamilton
propagation failure is the shift coupling. It is not a sign error in `D` or
`P`, and it is not present in the shift-free chart.

The momentum residual does not vanish at `beta = 0`. Its derivative is the
incomplete cancellation of two large pieces: the geometric momentum density
and the shift-current density. At N=64 those pieces have maxima `187.258`
and `172.780`, and their sum has maximum `46.067`. At N=128 the pieces are
`174.867` and `176.793`, and the sum is `2.603`. A wrong relative sign would
add them and would not shrink under refinement. Nyquist is a small fraction
of either derivative (`7.95` against `7346` for `dH/dt` at N=64).

## Phase mistake, separate from the propagation failure

The calibration number `phase = -0.842669616901` is the constant relative
phase of the minus spinor column. The checkpointed v1 source multiplied the
odd spatial lobe by that phase and then took the minus column to be the
conjugate. Its preparation flag
`phase_applied_to_odd_lobe_after_real_normalization` is true in the v1 JSON.
That changes the packet. Its overlap with the real-envelope packet is
`0.8327 + 0.3732 i`, not 1. The corrected source is described below and does
not rewrite this JSON.

Continuum packet on `[0, 2]`, same geometry coefficients:

| Construction | Off-diagonal | Absolute error against `i/5` |
|---|---:|---:|
| v1 odd-lobe phase | `-0.3318151250545991 - 0.12703307964704 i` | `0.4658883046379534` |
| Real envelope, no extra phase | `-0.14928446688289454 + 0.13309450757810007 i` | `0.16359155530994546` |
| Real envelope, phase on the minus column | `7.3e-14 + 0.20000000000000 i` | `1.1e-13` |

The recorded phase equals the argument that rotates the real-envelope
matrix element onto `i/5`, to `3.6e-13`. Diagonals are already `1` and `2`
to `1e-12` in all three constructions. A column phase does not change
`C = Φ diag(c) Φ†`. The lobe phase does. This mistake does not select the
constraint failure: with the v1 state held fixed, `dH/dt` still disappears
at `beta = 0` and still falls under refinement.

## Smallest justified correction

Do not change the Hamiltonian, project the constraints, reset the fields, or
insert a force. `D` and `P` already have the right signs and units. The N=64
Hamilton residual is unresolved high-mode aliasing of the static shift, and
the momentum residual is the truncation of a cancellation that is already
the right way up. Both fall from N=64 to N=256 (`37.0` to `0.576`, and
`0.294` to `0.00129`). N=256 at `T = 0.005` is still above the v1 validation
tolerance `1e-3` on the Hamilton residual, so the interval remains a failed
short validation.

## Preparation correction after the v1 checkpoint

The v1 JSON is unchanged and still records the odd-lobe phase. After that
file was checkpointed, the source was corrected without touching the
Hamiltonian or the derivative. Real lobes stay real. The plus carrier on
region `n` is `exp(+i k (x - n ell))`. The minus column is
`exp(i phase)` times the conjugate, which is the carrier
`exp(-i k (x - n ell))` with the recorded constant phase. That column phase
does not change `C`. The old hand-set link is still not copied.

Continuum onsite error against `Ω^n H`, three regions, maximum absolute
entry error `3.2e-12`. Region-0 off-diagonal error is
`1.099741404857301e-13`. Executed weight residual is `1.05e-13`.

Discrete onsite off-diagonal error against `Ω^n i/5`:

| N | Region 0 | Region 1 | Region 2 |
|---:|---:|---:|---:|
| 64 | `1.90e-4` | `2.86e-4` | `4.28e-4` |
| 128 | `2.91e-7` | `4.37e-7` | `6.56e-7` |
| 256 | `9.75e-11` | `1.46e-10` | `2.19e-10` |

Link difference from the old `B`, N=128: `5.9845119231619615` and
`8.97676788474295`. No successor result record is written. The v1 entry
point refuses to overwrite `nsc-spherical-coupling-control-v1.json`.
