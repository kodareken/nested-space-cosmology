# Reproducing the root project

The [remaining Lambda mechanisms](schur-memory.md#remaining-lambda-mechanisms-and-the-matched-nesting-test)
have a small paired one-level/two-level calculation:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_lambda_mechanisms.py
python scripts/run_python.py -m pytest tests/test_nsc_atlas_gravitational_reduction.py -q
```

It reuses the existing fixed quadratures, keeps scalar density and pressure
alongside their trace, and directly compares the common quiet exterior
Einstein tensors. Ten focused controls pass. No evolution, physical grid
scan, fitted source, result file or retained arrays are needed. The
nonzero scalar increment fails vacuum form; the exterior Lambda difference
is exactly zero. The actual collective parent remains unspecified.

The [inherited local Dirac limit](nsc-atlas-dirac-probe.md#inheritance-and-the-late-local-cap-state-19-october-2026)
has an analytic calculator, without a new parent benchmark or field grid:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_inheritance.py -q
python scripts/run_python.py scripts/check_nsc_atlas_dirac_inheritance.py
```

Seven controls check the actual quiet wall's regular horizon clock, full
four-component conformal field/energy, inherited remainder bound, central
KMS kernel, detector spectrum by an independent contour quadrature, units
and the conserved late vacuum source. Together with the existing clock
and bulk-trace controls,20 tests pass in2.625 process-CPU seconds. The
shared new-control allowance is30 seconds; all development runs stay below
it. The calculator prints conditional scales for a supplied10km ruler.
It does not save arrays, assign a finite-time ancestor covariance, derive
the cap's classical H or compute compact-source Psi/Phi.

The [infinite inherited-law proof](nsc-infinite-inheritance.md) includes
its domain, convergence bound and spectral dimension. Its small calculator
can be reproduced without an exponentially large tree or any P001 grid:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_infinite_inheritance.py -q
python scripts/run_python.py scripts/check_nsc_infinite_inheritance.py
```

The nine controls compare full finite trees, a closed-form infinite
response and separately assembled length/volume-weighted edge actions.
They check physical dilation, lapse variation and root work as well as the
root compliance, depth bound and calibrated dynamic response. The latter
compares the complete static-profile inertia, differentiated Schur recursion
and disjoint finite-band intervals for two nesting ratios. The combined controls and illustration
use .9375 process-CPU seconds, without physical grids or an observational fit. The analytic depth
bound excludes arithmetic and physical-model uncertainty. No full arrays
or result JSON file is created.

The [common-cap causal support check](nsc-atlas-dirac-probe.md#causal-support-before-another-finite-part-inference)
compares the existing plus-source rows before and after projection:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py::test_common_cap_plus_columns_have_exact_causal_support_before_projection -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --diagnose-causal-support
```

Freeze the method first. The C1 control uses1.296875 CPU seconds of a
30-second shared allowance. The physical comparison uses the central
N128 nested/cap maps, with field, PDE-spatial and finite-parameter source
rows. It retains no arrays on disk and adds no vacuum contraction or grid
resolution. A forbidden-support identity is not a continuum error bound.

The [finite algorithm derivative](nsc-atlas-dirac-probe.md#differentiating-the-existing-finite-cell-map)
uses the same Cayley field kernel, normalized source knots and constraints:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_mean.py::test_parameter_tangent_differentiates_the_same_rectangular_gate_map_and_flux_weights tests/test_nsc_atlas_dirac_mean.py::test_same_CAR_transfer_flux_weights_and_bounded_evaluation tests/test_nsc_atlas_dirac_retarded.py::test_finite_cell_projection_tangent_includes_density_scale_and_moving_constraints tests/test_nsc_atlas_dirac_retarded.py::test_discrete_retarded_adapter_matches_parameter_difference_at_fixed_source -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --diagnose-discrete-derivative
```

Freeze the method first. These controls use6.203125 CPU seconds; the
shared allowance is90 seconds. The same N128/two-split diagnostic tests
the finite reconstruction derivative. It does not establish continuum
accuracy, evaluate a vacuum contraction or add a P001 grid.

The [differentiated source projection](nsc-atlas-dirac-probe.md#differentiating-the-constrained-source-projection)
uses analytic second potential jets and the existing two boundary rays:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py::test_projection_derivative_recovers_cubic_finite_hilbert_moment_and_mean tests/test_nsc_atlas_dirac_retarded.py::test_projection_derivative_matches_fixed_source_differentiation_of_smooth_kernel tests/test_nsc_atlas_dirac_retarded.py::test_volterra_source_observation_targets_match_separable_potential tests/test_nsc_atlas_dirac_preparation.py::test_auxiliary_second_null_jets_match_independent_optical_taylor_algebra -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --diagnose-projection-derivative
```

Freeze the method first. The diagnostic compares the independently
propagated spatial row with the continuum projection derivative evaluated
on cell rows. Its nonconstant residuals and finite endpoint Hilbert moment
are compatibility checks, not a vacuum stress. Controls use8.3125 CPU
seconds of the90-second shared allowance. No new field resolution or
smaller split is admitted.

The [nonconstant observer-moment check](nsc-atlas-dirac-probe.md#nonconstant-source-moment-compatibility-check)
uses the same N128 maps and two existing separations:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py::test_observation_moments_keep_moving_boundary_signs_and_exact_quadrature tests/test_nsc_atlas_dirac_retarded.py::test_nonconstant_probe_detects_projection_commutator_hidden_from_mean -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --diagnose-observer-moments
```

Freeze the method before the diagnosis. It prints only residual summaries
for fixed source probes, including nested-minus-cap and the correction's
contribution. Two controls use 1.09375 CPU seconds; the shared allowance
is 60 CPU seconds. No finer grid, smaller split or vacuum contraction runs.

The [source-compatible reconstruction](nsc-atlas-dirac-probe.md#enforcing-the-volterra-source-derivative-before-recomputation)
preserves the kernel mean and exact endpoint value/source derivative:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py::test_source_constraints_preserve_mean_and_are_the_unique_minimum_L2_correction tests/test_nsc_atlas_dirac_retarded.py::test_source_compatible_reconstruction_retains_the_exact_cylinder_kernel_and_log -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --stage source-jets --output .build/nsc-atlas-dirac-root-response-source-jets.json
```

Freeze the method before the second command. The 11.75-CPU-second controls
check weighted optimality and the exact cylinder kernel/derivative. The
physical comparison keeps N128 and the same two splits, uses eight
covariance processes, and reproduces the four v2 legacy cases within 1e-9.
No finite coincident stress or full Psi/Phi follows from log cancellation alone.

The [frozen v3 result](../results/summaries/nsc-atlas-dirac-root-response-v3.json)
uses method `a097d42`. All four old cases reproduce exactly, the imposed
logarithmic identity passes, and the corrected order48 split variation
falls to .385%. The 26.8% smaller-split polynomial variation remains
unresolved. The eight-process run uses 92.296875 CPU/17.42970 wall seconds
and retains 20761 bytes; no finer field grid or smaller split is used.

The [source-endpoint logarithmic check](nsc-atlas-dirac-probe.md#source-endpoint-jets-and-a-relative-logarithmic-diagnostic)
uses analytic Volterra source jets and the same central detector:

```powershell
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --diagnose-source-logs
```

Freeze the method first. This prints a small coefficient table within the
30-CPU-second shared diagnostic allowance; it creates no data file and
does not evaluate a coincident finite stress. Three independent algebraic
controls use 1.109375 process-CPU seconds before the physical diagnosis.

The corrected diagnosis is pinned to `8bad871`. It prints nonzero legacy
relative pressure-log residues at both orders; substituting the exact
source jets cancels those coefficients. Calculation uses 2.75 CPU/1.57345
wall seconds, or 3.59375 CPU seconds including imports. The owner retains
the small coefficient table and its limits; no extra result JSON is saved.

The [fixed-grid polynomial comparison](nsc-atlas-dirac-probe.md#fixed-grid-polynomial-comparison-and-useful-process-parallelism)
keeps N=128 and compares orders 24/48. Its order-48 cylinder check uses the
existing N=64 grid and 8.125 process-CPU seconds:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py::test_order_48_uses_the_same_existing_grid_and_exact_cylinder_kernel -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --stage polynomial --output .build/nsc-atlas-dirac-root-response-polynomial.json
```

Freeze the method before the second command. Field maps are prepared once,
then independent covariance contractions run in eight processes. The old
order-24 values must agree with the frozen v1 baseline within 1e-9 before
the new order is interpreted. No additional field grid is admitted here.

The [v2 fixed-grid result](../results/summaries/nsc-atlas-dirac-root-response-v2.json)
uses method `e6bd1bb`, 54.625 total CPU seconds and 12.77376 wall seconds.
The old cases reproduce, but the 24/48 spread grows rather than contracts.
Keep both pilot amplitudes provisional; the owner records the next
endpoint-identity diagnosis without admitting a finer field grid.

The [retarded-map extraction](nsc-atlas-dirac-probe.md#extracting-retarded-maps-from-the-existing-field-kernel)
uses the existing characteristic solver and analytic derivatives. Its
32/64 cylinder control checks both polynomial orders and contracts the
complete radial kernel and derivative errors by about four:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_retarded.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_root_response.py --output .build/nsc-atlas-dirac-root-response.json
```

The second command requires a frozen method and an unused output path.
It compares the actual nested/cap state on two nonzero splits, retaining
the first signed angular sector only. Neither the coincidence limit nor
the angular sum or full Psi/Phi is supplied. Control development used
46.84375 process-CPU seconds; the owner specifies the physical run's
separate budget and stopping conditions.

The [frozen v1 full-root result](../results/summaries/nsc-atlas-dirac-root-response-v1.json)
uses method `dbbfba7`: 47 CPU seconds, 42.47898 wall seconds and 22899
retained bytes. It supplies a finite-split first-sector vacuum difference,
with polynomial and separation effects larger than the field-grid spread.
No extra quadrature or finer P001 grid is needed for that diagnosis.

The [complete root-covariance pairing](nsc-atlas-dirac-probe.md#pairing-the-complete-root-covariance-in-null-coordinates)
has four controls against the exact cylinder vacuum kernel and its spatial
derivative, including its one-sided endpoint and Hilbert terms:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_root_pairing.py -q
```

Development and derivative regression used 13.09375 process-CPU seconds
with no physical field/P001 grid or retained arrays. This integrates the
full radial root covariance on supplied polynomial maps. Their numerical
representation, the coincidence limit and the angular sum remain external.

The [global Dirac bounds and real boundary symbol](nsc-atlas-dirac-probe.md#real-frequency-boundary-ellipticity-for-the-next-regularity-estimate)
have six controls, including sharp balanced endpoints, the explicit
decaying Poisson lift and an exact
constant-mass line resolvent:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_boundary_ellipticity.py tests/test_nsc_atlas_dirac_bounds.py -q
```

They use 1.953125 process-CPU seconds across development and create no
reference data. Calculator inputs must be genuine uniform geometry bounds;
the outputs are not a bound on the full evolved stress.

The [correlated parent and auxiliary preparation](nsc-atlas-dirac-probe.md#exact-correlated-root-data-and-its-finite-compressions)
have seven controls, using exact root kernels, manufactured plane modes
and the common equal-time kernel:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_cylinder.py -q
```

The original six used 3.546875 process-CPU seconds across development;
the recovered common-kernel control used another 1.28125. They retain no
new reference JSON or arrays. The full continuum state is prescribed; its
propagated local stress remainder is not inferred from these controls.

The [radar subtraction-symbol adapter](nsc-atlas-dirac-probe.md#reusing-the-fourth-order-symbol-on-the-actual-radar-geometry)
uses analytic one-sided metric jets through order four and the existing
Weyl projector/transport core. It runs no physical field reference:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_symbol.py tests/test_nsc_spatial_reference_symbol.py -q
```

Seven independent controls pass (7.875 process-CPU seconds across development).
The local symbol is not a state or a numerical UV-remainder estimate.

The [finite clock-energy preparation](nsc-atlas-dirac-probe.md#actual-geometry-clock-energy-candidate-admission)
has thirteen preparation and existing-field controls, including the exact
continuum-versus-finite Hardy projector comparison. The actual-geometry
reference uses sixteen independent cases on eight process workers, with
one BLAS thread each. It stores spectra, summaries and matrix hashes only:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_preparation.py tests/test_nsc_atlas_dirac_mean.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_preparation.py --output .build/nsc-atlas-dirac-preparation.json
```

Freeze the source before the second command and use an unused output path.
This calculates a finite state candidate, not a continuum vacuum or full
Psi/Phi, and runs no P001 grid. Physical preparation duration is held fixed
between field resolutions; changing mode count is a separate cutoff check.
The [frozen candidate result](../results/summaries/nsc-atlas-dirac-preparation-v1.json)
uses method `0892dbc`. It resolves the finite field eigenspaces but fails
to give a mode-count-stable tidal value; the owner records that conclusion
without interpreting the two/four-mode difference as an observational error.

The conserved compact-source tangent has five independent nonlinear
wall/ray and metric controls. It supplies a mass-dependent input to a
later Dirac calculation, with no new physical field or P001 evolution:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_compact_source.py -q
```

The [source owner](nsc-atlas-dirac-probe.md#conserved-compact-source-mass-tangent)
states the fixed initial speed, required parent-mass variation, annular
domain, analytic optical derivatives and parametric rather than retarded
interpretation. Development controls use 3.84375 process-CPU seconds.
The extension now has eight source-family controls, plus the five existing
GR controls on the changed reconstruction (19.15625 CPU seconds shared
across development). Three independent local curvature controls use
another 1.8125 CPU seconds. The bounded physical source-response command is

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_tidal.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_mass_response.py --output .build/nsc-atlas-dirac-mass-response.json
```

This command requires frozen source and an unused output path. Its
[admission](nsc-atlas-dirac-probe.md#prepared-compact-source-response-admission)
fixes the three masses, two field resolutions, two quadrature orders,
eight-process pool and 480-CPU-second limit. It returns only summaries;
the local tidal ratio is distinct from the full potential ratio Psi/Phi.
The [frozen v3 result](../results/summaries/nsc-atlas-dirac-mass-response-v3.json)
uses method `358bf25` and the exact triangular wall primitive, jointly
checking wall/bulk Gauss orders 8/12. It passes in 6.97495 calculation
wall seconds and 50.71875 charged CPU seconds. The failed v1 derivative,
its same-field precision diagnosis and v2 setup exit are retained in
the owner; reproducing v3 does not run them or create another field grid.

The prepared Dirac mean stress and spherical metric/clock reconstruction
use the same quiet geometry and existing wall-mass coefficients:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_mean.py tests/test_nsc_atlas_dirac_backreaction.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_mean_clock.py --output .build/nsc-atlas-dirac-mean-clock.json
```

The recorder requires frozen source and an unused output path. It runs
only 64/128 field cells with fixed quadrature orders 8/12, using up to
eight independent process workers and one BLAS/Numba thread per worker.
The [owner](nsc-atlas-dirac-probe.md#mean-stress-and-a-direct-spherical-metric-reconstruction)
specifies the state difference, shared observer motion, corner data,
proper clock, uncertainty and scope. No P001 evolution is invoked.
The [frozen mean-clock reference](../results/summaries/nsc-atlas-dirac-mean-clock-v1.json)
uses method `be1a11e` and occupies 13053 bytes. It completes in 11.3722
wall seconds and 40.265625 observed process-CPU seconds, with a conservative
48.265625 total charge. Mass-path gaps contract and every admitted
validation condition passes. The scale and numerical spreads are defined
in the owner; this is a conditional finite mean response, not a baseline
vacuum or full Psi/Phi.

The actual Dirac gravitational/observer bridge is an equation interface,
with independent four-dimensional curvature, constraint propagation,
moving-wall and proper-clock controls:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_gravity.py -q
```

These six controls use 4.390625 process-CPU seconds across development
executions, including one scalar-array fixture correction. They generate
no physical evolution or result JSON. The
[owner](nsc-atlas-dirac-probe.md#einstein-moving-wall-and-proper-clock-response)
retains the complete conserved source and constrained initial data as
requirements for the next gravitational solve.

Reproduce the actual curved Dirac transfer and ordered metric response:

```powershell
$env:NUMBA_NUM_THREADS="8"
$env:OPENBLAS_NUM_THREADS="1"
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_scattering.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_scattering.py --output .build/nsc-atlas-dirac-scattering.json
```

The recorder requires frozen source and an unused output path. It evaluates
only the admitted 64/128 characteristic Dirac grids, not P001 grids.
Positive-k multiplets are separate prepared occupation differences, not
a vacuum sum. The
[owner](nsc-atlas-dirac-probe.md#causal-transfer-inherited-correlations-and-nonlocal-metric-response)
defines the supports, conditional inheritance and remaining slip inputs.
The committed [comparison](../results/summaries/nsc-atlas-dirac-scattering-v1.json)
uses method `7567c2a` and occupies 22783 bytes. Its first save exceeded
the summary size gate after completing both resolutions. The successful
completion charges 30 prior CPU seconds with `--prior-cpu-charge 30`,
then uses 10.234375 CPU seconds. Both stay within the original allowance;
the current recorder's default zero prior charge is suitable for a fresh
independent reproduction.

The actual wall-clock radar construction supplies a finite causal diamond
for the quiet Dirac geometry, using its proper-wall null labels:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_radar.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_radar.py --output .build/nsc-atlas-dirac-radar.json
```

Commit the method first and use an unused output path. Optional
`--figure .build/nsc-atlas-dirac-radar.png` draws the same geometry.
The [owner](nsc-atlas-dirac-probe.md#proper-wall-radar-coordinates-give-a-finite-causal-propagation-chart)
keeps analytic domain bounds, sampled checks, the physical clock protocol
and state limitations separate. It does not run a field or a P001 grid.
The [saved reference](../results/summaries/nsc-atlas-dirac-radar-v1.json)
and [figure](../paper/figures/dirac-wall-radar.png) retain the declared
background protocol and producer binding; full metric arrays regenerate.


The fixed quiet time-chart calculation preserves a coordinate-domain exit
as an unsuccessful full-domain result:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_time_chart.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_time_chart.py --output .build/nsc-atlas-dirac-time-chart.json
```

The method must be committed and the output path unused. It uses the
[predeclared gauge and cuts](nsc-atlas-dirac-probe.md#finite-joined-time-chart),
stops on the first failed bulk metric cut and keeps its location. Controls
do not call the nine-cut report. Earlier unfrozen implementation previews
are disclosed in the admission; no profile is tuned to hide their fold.
The [saved reference](../results/summaries/nsc-atlas-dirac-time-chart-v1.json)
is an unsuccessful full-domain chart. A successful recorder exit means the
outcome was saved; `status=coordinate_domain_exit_recorded` does not mean
that the full time chart passed.

The shared original action, material junction and local Legendre-map
identities reproduce without a trajectory or generated result file:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_conformal_action.py tests/test_nsc_atlas_collar_field.py -q
```

The [owner](nsc-atlas-conformal-action.md) retains the boundary divergence,
physical radius kink and canonical/velocity distinction. This is not a
new P001 evolution or a correction to its archived mass values.

After committing the method, one full-strip manufactured scalar control
uses the exact quiet normal collar, a shared six-coordinate field, both
outer GHY terms and the restricted ADM momentum:

```sh
python scripts/run_python.py scripts/measure_nsc_atlas_collar_field.py --output .build/nsc-atlas-collar-field.json
```

The output path must be unused. The
[finite-action owner](nsc-atlas-conformal-action.md#a-shared-finite-field-action-on-a-quiet-wall-normal-collar)
states the continuum test source and energy/work identity. This is a
method control on prescribed geometry, not a P001 grid or physical signal.
The [saved application](../results/summaries/nsc-atlas-collar-field-v1.json)
contains only the summary, limits and frozen source hashes; full trajectories
are regenerated and are not retained.

The local conformal Dirac interface and bulk curvature calculation reproduce
with the committed method and an unused output path:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_conformal.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_conformal.py --output .build/nsc-atlas-dirac-conformal-interface.json
```

The [owner](nsc-atlas-dirac-probe.md#conformal-calculation-variables-remove-the-singular-wall-jet)
specifies the conformal weights and their limits. This calculation takes no
evolution steps and does not select a quantum state or compute a loop.

The same committed producer can separately evaluate the local one-loop
Dirac bulk trace in the original quiet physical metric:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_trace.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_conformal.py --bulk-trace --output .build/nsc-atlas-dirac-bulk-trace.json
```

This imports free Dirac quantization, calculates its spin/statistics heat
coefficient and uses the existing one-sided curvature invariants. The
`10 km` length is supplied independently. It computes no nonlocal loop
polarization, singular-wall term or slip.
The two compact reference records are
[conformal matching](../results/summaries/nsc-atlas-dirac-conformal-interface-v1.json)
and [original-metric bulk trace](../results/summaries/nsc-atlas-dirac-bulk-trace-v1.json).

The same producer has a separate Euler interface calculation in its
[declared conformal smoothing class](nsc-atlas-dirac-probe.md#euler-curvature-contribution-at-the-material-interface):

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_euler_wall.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_conformal.py --euler-wall --output .build/nsc-atlas-dirac-euler-wall.json
```

Commit the method first and use an unused output path. The exact normal
curvature and Chern checks require no smoothing scan. The result is the
Euler part of the wall trace, excluding `Box R`, finite interface terms
and the state-dependent nonlocal response. It is not a slip calculation.
The [frozen surface result](../results/summaries/nsc-atlas-dirac-euler-wall-v1.json)
retains the arithmetic, stated smoothing class, excluded terms and source
hashes only.

The same-action Dirac stress and bare metric contacts reproduce without a
reference trajectory or new state:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_stress.py -q
```

The [owner](nsc-atlas-dirac-probe.md#dirac-stress-and-the-actual-gravitational-metric-vertices)
retains the angular normalization, normal-frame energy/flux, full covariant
conservation and open-port term. This does not evaluate a loop or a slip.

The physical initial Dirac time generator and its moving-wall compatibility
reproduce with:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_evolution.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_wall_time.py --output .build/nsc-atlas-dirac-wall-time.json
```

The method must be committed and the output path unused. It evaluates the
same quiet source anchor algebraically and takes no field evolution steps.
The [owner](nsc-atlas-dirac-probe.md#physical-time-generator-at-the-moving-interface)
states the common initial lapse, its necessary jump, and the frame change;
a full foliation and quantum state remain unselected.
The [saved reference](../results/summaries/nsc-atlas-dirac-wall-time-v1.json)
contains the initial lapse and moving-domain coefficients only.

The source-anchor joined initial geometry and its intrinsic spatial Dirac
operator reproduce with:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_slice.py tests/test_nsc_atlas_dirac_probe.py -q
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_slice.py --output .build/nsc-atlas-dirac-source-slice.json
```

The producer requires its method to be committed and its output path unused.
It uses the same quiet sources at `Rwall/ell=.45`, a regular outgoing
Eddington--Finkelstein frame, and fixed 32/64 quadrature on the slicing taper.
No P001 time step, new evolution grid, quantum state or loop polarization is
introduced. The [owner](nsc-atlas-dirac-probe.md#joined-initial-surface-at-the-source-anchor)
distinguishes the intrinsic spatial operator from four-dimensional dynamics.
The [saved scalar record](../results/summaries/nsc-atlas-dirac-source-slice-v1.json)
binds the frozen method and contains no spinor or metric arrays.

The actual archived cap RHS check reproduces with:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_wall_rhs.py --output .build/nsc-wall-pointwise-rhs.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_wall_rhs.py --analyze-existing results/summaries/nsc-atlas-wall-pointwise-rhs-v1.json --output .build/nsc-wall-boundary-lift.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_pointwise_rhs.py -q
```

It requires the original ignored 3333-byte level-0 v3 seed; the producer
checks its binding and never starts an evolution to recover a missing seed.
It rebuilds the original history once, evaluates one unchanged coupled
stage and projects its pointwise field equations on the existing panels.
The [saved result](../results/summaries/nsc-atlas-wall-pointwise-rhs-v1.json)
contains no raw arrays beyond seven-component load summaries. Hashes and
dimensions are exact; stage recovery uses `1e-10*max(1,value)` and the
projected load identity uses the recorded `1e-8` scale. Compare local rates
within `1e-8` absolute/relative tolerance, not byte equality of runtime or
commit metadata. The diagnostic does not advance or correct P001.
The analysis mode reads the committed saved point and computes only the
boundary Riesz gains of the existing three spaces. The
[lift record](../results/summaries/nsc-atlas-wall-boundary-lift-v1.json)
separates the actual level-0 density force from its boundary projection;
it supplies no evolved force at levels 1 or 2. Compare gains within `1e-9`
relative tolerance and the lifted consistency residual within `1e-6`
absolute tolerance; neither is a physical uncertainty bound.

The actual cap--wall Dirac inputs and P001 active-region diagnosis reproduce
without an evolution or new grid:

```sh
python scripts/run_python.py scripts/measure_nsc_atlas_dirac_probe.py --output .build/nsc-atlas-dirac-probe-reference.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_resolution.py --output .build/nsc-atlas-endpoint-active-resolution.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_dirac_probe.py tests/test_diagnose_nsc_atlas_endpoint_resolution.py -q
```

Use new output paths. The first command requires its producer/component to
match the current committed source. Compare the algebraic geometry values
within `1e-12` relative/absolute tolerance and exact domain/status labels;
runtime and source-commit metadata may differ. It computes no nested slip.
For the second, receipt hashes and support-cell indices must match exactly;
fronts and widths can be compared within `1e-12` in the declared units.
The input comparison authenticates the original v4 endpoint receipts. Its
incoming support concerns the retained background, not unarchived departures.
The [operator record](../results/summaries/nsc-atlas-dirac-probe-reference-v1.json)
and [resolution record](../results/summaries/nsc-atlas-endpoint-active-resolution-v1.json)
remain compact; no raw array or previous evolution is regenerated.

For the later same-background incoming diagnostic:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_projection.py --output .build/nsc-endpoint-frozen-projection.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_projection.py --analyze-existing results/summaries/nsc-atlas-endpoint-frozen-projection-v1.json --output .build/nsc-endpoint-incoming-control.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_endpoint_incoming_control.py -q
```

The first command regenerates the original prepared history once and takes
three instantaneous cap solves on the existing spaces, with zero departures
and diagnostic rates. It is not a P001 trajectory replay. Compare mass,
lapse and gap values within `1e-9` absolute/relative tolerance; the original
prepared-input hash must match exactly. The second command uses the already
committed control and performs only algebra; compare its derived gaps within
`1e-12`. It copies measured incoming mismatches and is not a prediction of
an independent source. The first attempt's NumPy/JSON format failure is
retained in the frozen control; runtime metadata are not equality targets.

The current Step-7 exact calculations and causal exterior operator reproduce with:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_effective_vacuum.py --output results/summaries/nsc-atlas-effective-vacuum-v1.json
python scripts/run_python.py scripts/record_nsc_atlas_geometric_average.py --output results/summaries/nsc-atlas-geometric-average-v1.json
python scripts/run_python.py scripts/record_nsc_atlas_exterior_schur.py --output results/summaries/nsc-atlas-exterior-schur-metric-v1.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_effective_vacuum.py tests/test_nsc_atlas_exterior_dtn.py tests/test_nsc_atlas_geometric_average.py -q
```

The compact record contains the classical and conditional quantum trace
coefficients, fiducial SI inputs, interface polynomial, zero residuals and
source hashes. The response tests include the existing prepared pulse and
an exact constant physical field. They use the quiet reference geometry;
they do not rerun the finite-amplitude P001 endpoint trajectories.
The exterior Schur record evaluates that actual prepared pulse's conformal
Robin response and Weyl metric check on one 129-node grid. Its additional
local classical bulk Lambda is zero for variations supported strictly inside
the retained cap at fixed interface/exterior data. The nonzero boundary
pairing is recorded separately. The combined 120 CPU-second / 30-million-value
budget covers both prepared and prescribed-wall solves and their queries.
No quantum determinant or covariance is inferred from a retarded pairing.

The quiet-cap distance comparison reproduces with:

```sh
python scripts/run_python.py scripts/measure_nsc_atlas_distance_contrast.py --fetch --output .build/nsc-distance-contrast-reproduced.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cosmological_distance.py -q
```

Use a new output filename. The script retrieves two pinned publisher files
into `.build/data/pantheonplus-c447f0fea703` (about 34 MB), checks their Git
object hashes and dimensions, and uses the full joint covariance. Raw input
files are ignored; the [saved summary](../results/summaries/nsc-atlas-distance-contrast-v1.json)
contains the result, uncertainty, selection, hashes and exact zero quiet-cap
prediction. Compare contrast and standard uncertainty within `1e-10 mag`,
with input hashes and row counts exact; runtime/commit metadata and tiny solve
residuals can differ. Five synthetic controls need no downloads. The observed
contrast is a quiet-limit consistency check, not an induced-Lambda fit or a
withheld observation for a later model calibrated from these data.

The geometric record uses one existing initial slice through the right-end
cutoff `R/ell=4`, retaining the wall curvature distribution and boosted stress.
It records three quadrature methods at that same domain, their differences,
and the conditional parent-volume test. Its geometric density and pressure
fail the vacuum relation; the trace coefficient is not an observed Lambda.

The [parent-side metric calculation](schur-memory.md#whole-body-gravitational-reduction-at-the-parent-end)
reuses the quiet body's mass and the same outer cut:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_gravitational_reduction.py --output .build/nsc-atlas-gravitational-reduction-reproduced.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_gravitational_reduction.py -q
```

Its [summary](../results/summaries/nsc-atlas-gravitational-reduction-v1.json)
contains the boundary energy, surface pressure, exact action identities and
source hashes. Numerical tube values can be compared within `1e-12` relative
tolerance; the exact symbolic checks are zero. This is a 30-CPU-second
algebraic/reference calculation without field evolution. The coefficients
multiplying number density are recorded without inventing that density.
There is no complete parent solution or physical error estimate for the
omitted interactions.

The one-stage compiled-coefficient diagnostic also reports reconstructed
cap-rate Galerkin norms in the model's `dz` measure and physical wall traces.
Those norms are invariant under a consistent basis rescaling; they are not
proper-volume energy norms or a trajectory error bound. The original
coefficient-rate comparison remains visible with its original tolerance.
The diagnostic does not change the backend used by the running endpoints.

There is one package, `nested-space-cosmology`, importing `recursive_horizons`
from root `src`. Use the pinned versions in [pyproject.toml](../pyproject.toml).
Scientific conclusions retain the domains in their owning notes.
The current paper and test dependencies use Matplotlib 3.11.2; the scientific
stack includes python-flint 0.9.0 and Numba 0.68.0.

The [same-family formation diagnostic](nsc-atlas-finite-amplitude.md#formation-diagnostic-from-the-quiet-sphere-and-material-wall)
has four exact controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_formation_diagnostic.py -q -s
```

They check area/boost/mass consistency, affine focusing with a regular
caustic control, agreement of surface-stress and Israel impulses, and the
conformal affine relation directly from the metric connection. The shared
allowance is twenty active CPU seconds and zero physical samples; the
recorded run uses 0.296875 seconds. No global formation simulation is run.

The [earlier spacelike realization](nsc-atlas-finite-amplitude.md#earlier-spacelike-realization-of-the-prescribed-source)
has four exact controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_past_preparation.py -q -s
```

They check the past-directed primitives, source transverse ODE comparison,
Jordan Cauchy frame and scalar derivatives, and extrinsic curvature from
the metric connection. The shared allowance is twenty active CPU seconds
and zero physical samples; the recorded run uses 0.765625 seconds. The
analytic collar width follows from the owner's range and contraction
inequalities, not a numerical scan or a physical field integration.

The [exact null-preparation normalization](nsc-atlas-finite-amplitude.md#exact-null-stress-normalization-of-the-prescribed-preparation)
has four symbolic controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_null_normalization.py -q -s
```

They check the full-stress boundary identity, compact-support weak-field
limit, affine/dimensional scaling and the limit of the weight comparison.
The admitted shared budget is ten active CPU seconds and zero field samples;
the recorded run uses 0.515625 seconds. No original pulse integral, physical
source solve or numerical quadrature is part of this command.

The cap response solution also has `wall_sample(t)`, which returns its
wall-radius and polar-slice clock coefficients at fixed exterior label
without evaluating an arrival ray. The [coordinate/source bridge](nsc-atlas-finite-amplitude.md#connecting-first-coefficients-at-a-fixed-source-label)
converts original affine source label xi to the existing normalized query
t. Numerical source, geometry, field and clock-quadrature errors still
require their own bounds. Its manufactured controls run with:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_wall_coefficient_pullback.py -q -s
```

This command performs no physical source solve, evolution or grid scan.
The cap module's registered implementation successor authenticates its
original result-producing bytes; historical source-context commands below
retain their declared snapshots.

The direct onset comparison has independent positive-kernel controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_onset_endpoint_bound.py -q -s
```

After freezing its producer sources, the one-event application is:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_onset_endpoint_bound.py --output results/summaries/nsc-atlas-onset-endpoint-bound-v1.json
```

Use a new output path when reproducing an existing receipt; it will not be
overwritten. The [first receipt](../results/summaries/nsc-atlas-onset-endpoint-bound-v1.json)
binds frozen producer `2052a1d44da698f20b1c4180317d005e2596d4d1` and its exact
source bytes. The command authenticates the original source-cell and reference
image, then encloses the leading wall field and its first affine-source
derivative at the original right endpoint. It evaluates no trajectory or
sampled field and performs no grid scan. Radius and clock coefficient
integration is handled by the following separate component; their
finite-amplitude remainder constant remains unevaluated.

The first geometric coefficient enclosure has independent exact controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_onset_geometry_bound.py -q -s
```

After its sources are committed, one fixed 128-panel interval pass runs with:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_onset_geometry_bound.py --output <new-output-path>
```

The output parent must exist. This command authenticates the original source,
reference and scalar endpoint receipts. It bounds the continuum mass,
Einstein radius and anchored polar-clock first coefficients at the original
right source label. The panel ranges are enclosures of the analytic source;
there is no evolved field grid, geometric ODE, adaptive refinement or retry.
Its independent controls use 19,048 charged values / 0.65625 active CPU
seconds across the initial execution and one corrected reference-identity
assertion. The one physical application has a separate 200,000-value /
30-CPU-second limit and retains only endpoint summaries and hashes.

The [first geometric receipt](../results/summaries/nsc-atlas-onset-geometry-bound-v1.json)
uses producer `7a38bc965c822bf63fd66913a2ed38fb7356dc79`. Its original
128-panel application uses 140,288 charged values / 0.03125 CPU seconds;
the separate authentication/display audit uses 1,024 / 0.046875.
Both fixed-source radius and clock first coefficients are strictly
negative. Reproduction requires a new output path and the frozen source
version; later finite-amplitude errors are outside this certificate.

The amplitude-proportional geometric comparison has independent controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_geometric_departure_bound.py -q -s
```

After freezing the producer, its one algebraic application is:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_geometric_departure_bound.py --output <new-output-path>
```

The parent directory must exist. It authenticates the original source,
reference and strict coupled-tube receipts, extracts four affine directions
with those coefficients held fixed, and solves one rational three-by-three
comparison. It performs no source query, field evolution, physical time
step, parameter scan or retry. Its 32,768-value / 10-CPU-second allowance
is separate from the seven manufactured controls, which use 11,264 values
and active CPU time below timer resolution. This supplies a geometric
first difference, not a scalar difference or quadratic remainder constant.

The [first departure receipt](../results/summaries/nsc-atlas-geometric-departure-bound-v1.json)
uses frozen producer `459b9fdf6ad8fde233ce5f7137179eeda5a3c8bd` and
28,672 charged values / 0.03125 CPU seconds. The separate rational audit
uses 1,024 values / 0.09375 CPU seconds. The matrix norm is below 0.146;
the exact solution and dependent rows have zero rational residual.
Use a fresh output path when reproducing this retained certificate.

The canonical scalar comparison uses the same original onset and geometric
certificate. Its ten manufactured controls use 6,656 charged values and
0.421875 active CPU seconds in the first execution:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_departure_bound.py -q -s
```

The bounded original application requires the producer committed first and
a new output path in an existing directory:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_scalar_departure_bound.py --output <new-output-path>
```

It authenticates the source, reference, strict tube and geometric departure
inputs, then performs algebra only under a 24,576-value / 10-CPU-second
allowance. It saves normalized scalar difference bounds, comparison
constants, hashes and the command, with no source or field history. The
mathematical scope is in the [scalar comparison](nsc-atlas-finite-amplitude.md#amplitude-proportional-canonical-scalar-departure).
The [first scalar receipt](../results/summaries/nsc-atlas-scalar-departure-bound-v1.json)
uses producer `a0fce1a9a77b0a5c70463dd656ba2288a29161d2` and completes
18,944 charged values / 0.03125 CPU seconds. Its separate exact audit uses
2,048 values / 0.109375 CPU seconds. Both loops and all provenance, rounding
and alpha-scaling checks pass. It is a normalized scalar difference bound;
no numerical quadratic geometric remainder is supplied.

The radial cap consumer takes the same tube and completed scalar certificate.
Its ten independent controls initially use 6,144 values / 0.625 active CPU
seconds. Run them with:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_departure_bound.py -q -s
```

The original-data application requires a committed producer and new output
path, and uses at most 12,288 charged values / 10 CPU seconds:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_cap_departure_bound.py --output <new-output-path>
```

It evaluates positive algebraic cap bounds without radial integration,
source queries or field evolution. Its [mathematical scope](nsc-atlas-finite-amplitude.md#normalized-radial-cap-departure-and-vacuum-subtracted-residuals)
distinguishes normalized N/D differences and actual-radius vacuum-subtracted
residuals from the complete geometric Taylor remainder.
The [first radial receipt](../results/summaries/nsc-atlas-cap-departure-bound-v1.json)
uses producer `d5b7794f40884d5506c639622a3e71a5a4c13660`, completing
7,680 values / 0.015625 CPU seconds. The independent positive-product
audit uses 2,048 values / 0.09375 CPU seconds. Exact inequalities,
rounding, alpha powers and provenance all pass.

The quadratic source consumer reuses the original source-forcing receipt.
Its ten controls initially use 6,656 values / 0.40625 active CPU seconds:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_source_remainder_bound.py -q -s
```

After committing its producer, run one exact-rational application to a new
output path under the 8,192-value / 5-CPU-second allowance:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_source_remainder_bound.py --output <new-output-path>
```

It bounds four [source-geometry remainders](nsc-atlas-finite-amplitude.md#quadratic-remainders-of-the-original-source-geometry)
without source integration, field evaluation or pointwise h0 reconstruction.
The coupled radius and clock remainder remains separate.
The [first source remainder receipt](../results/summaries/nsc-atlas-source-remainder-bound-v1.json)
uses producer `7aa098a4dec0fbdc0ac08fb43f3cc86d462ac57b`, completing
3,584 values with CPU below timer resolution. The independent rational
product audit uses 1,024 values / 0.109375 CPU seconds; constants,
rounding, alpha powers, displays and provenance all pass.

The full onset geometric remainder combines the source/scalar/cap receipts
with interval second jets and the existing geometric matrix. Its ten
manufactured controls initially use 13,312 values / 5.609375 active CPU
seconds:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_geometric_remainder_bound.py -q -s
```

After committing the producer, the single original application requires a
new output path and uses at most 49,152 values / 15 CPU seconds:

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_geometric_remainder_bound.py --output <new-output-path>
```

It authenticates compact certificates and historical input hashes without
sampling historical arrays. The [comparison owner](nsc-atlas-finite-amplitude.md#uniform-quadratic-geometric-comparison-on-the-original-onset)
specifies the original finite domain, half-Hessian factors, reused matrix
and strict endpoint threshold. No source/field evolution or grid scan is
part of this calculation. The first application stopped on an indeterminate
Arb power of a zero-crossing derivative ball. The corrected implementation
uses interval multiplication for those squares; eleven controls pass in
14,336 values / 5.609375 active CPU seconds. The stopped v1 receipt is
preserved. The [corrected v2 receipt](../results/summaries/nsc-atlas-geometric-remainder-bound-v2.json)
uses producer `7e4734bc75ce04f9ba3afef3addcdfe700f5e8d5`, 37,376 values /
0.03125 CPU seconds. The independent 4,096-value / 0.09375-CPU-second audit
passes. Both fixed-label endpoint signs hold throughout the admitted
positive-amplitude range; all nine uniform remainder coefficients are
retained with the exact inputs and source bindings.

The material-radius conversion keeps the same original inputs and evaluates
one positive half-Hessian bound for the exact canonical-to-Jordan frame map:

```powershell
python scripts/run_python.py -m pytest tests/test_nsc_atlas_jordan_radius_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_jordan_radius_bound.py --output <new-output-path>
```

The six manufactured controls use 15,872 values / 0.3125 active CPU seconds
across the initial five-test execution and the added radius-floor guard.
The admitted original application has a 12,288-value / 10-CPU-second ceiling;
its output path must be new and its producer committed. The
[frame owner](nsc-atlas-finite-amplitude.md#jordan-material-radius-and-the-retained-clock-convention)
keeps the canonical scalar, invariant Jordan trace, source-label event shift
and clock interpretation explicit. No source or field integration is involved.

The [first material-radius receipt](../results/summaries/nsc-atlas-jordan-radius-bound-v1.json)
uses producer `d897d7969b7238d71c03948bdc83da5fa3f24779`: 7,680 values /
0.03125 CPU seconds in one application. The independent rational audit uses
2,048 values / 0.109375 CPU seconds and passes, including a positive-series
cosh enclosure independent of Arb. The material-radius sign is strict on
the entire admitted positive-amplitude interval; the clock receipt is
retained unchanged with its original slicing and quiet-centre anchor.

Saved evidence bytes, input identities, array shapes, discrete counts and
symbolic identities retain exact checks. Recomputed floating-point arithmetic
uses the numerical allowance of its validator. The clock boundary checker
uses relative `2e-11` / absolute `8e-12`; the rectangle diagnostic comparison
uses relative `2e-12` / absolute `2e-13`. Compact-source checks propagate
binary64 coordinate rounding through the source polynomial and its derivatives.
Source-frame checks also propagate Hermite reconstruction roundoff through
derived constraint residuals; cancellation is assessed against its operands.
Reference-geometry derivative allowances account for the actual spline knot
spacing, while retained coordinates and copied inputs remain exact where declared.
Restored forcing grids allow `2e-13` absolute drift in the independently rebuilt
wall coordinate. Their saved affine node spacing and fixed source endpoint
still use the original eight-epsilon check; saved arrays are never moved.
Restored reference splines evaluate a saved endpoint just outside their live
domain at the live endpoint, within the already admitted eight-epsilon bound.
Interior queries are unchanged; larger excursions and nonfinite values fail.
Departure-history constraint replay propagates rounding through the analytic
vacuum, terminal subtraction, centred stencils and moving-grid coefficients.
Only `E_Ubar`, `E_x` and their numeric summary leaves receive those additional
operand-dependent allowances. Mass, source constraints, counts, masks and
other diagnostics keep their existing checks. The eight-unit elementary-function
assumption is a replay policy, not a certified SciPy accuracy bound.
The compatible, background, SSPRK3 and midpoint-complete validators all
forward these transient bounds to their departure-history comparisons.
Midpoint profile brackets are checked against their own floating query and
knot coordinates. A changed integer bracket must be adjacent and correspond
to one unique common knot within sixteen epsilon units of both coordinates;
crossing counts must exactly follow those brackets. Saved-cut coordinates
and metadata remain exact. This handles a one-ULP knot-side difference seen
in GitHub's PDF authentication, without a general integer tolerance or any
change to saved data. The standalone midpoint check passes with 1,418,330
charged points, two final radial stages and no source solve or evolution.
Profile-fit replay similarly propagates rounding through its polynomial basis
conversion and differentiated residuals. Adaptive inverse iterations and their
work totals are checked for internal consistency and algorithmic bounds in each
run; floating-point stopping decisions need not match across platforms. Historical
budget receipts use the saved costs, while replay charges its actual work.
These are replay allowances, not physical uncertainty or convergence bounds.
One-ULP conditioning controls and deliberate larger corruptions test the
acceptance and rejection sides.

The [validation successor record](storage/validation-source-successors.json)
identifies the original Git bytes and the reviewed current checker bytes.
Old result manifests remain intact, current validation code runs, and newly
created records bind current code. Unknown source changes and all data-hash
changes still fail. Verified original Git bytes are cached within a process;
current files, registry entries and checker hashes are checked on every call.
Named `implementation-successor` entries also preserve the original producer
bytes when a calculation evolves. This classification does not assert
equivalent numerical results: the cap-field quiet branch and new history
adapter remain separate from the older evidence. The wall, overlap, cap and
scalar-boundary checkers authenticate saved inventories without rerunning
their producers; creation still binds current bytes. The October 17 update
passes 62 source-authentication and recorder tests. Check results report
`reviewed-successors` when stored producer bytes differ from current files;
`current` is reserved for exact current bindings. The degree-seven
characteristic profile also has an explicit implementation-successor entry.

Compact summaries use the same authentication for producer bindings while
keeping input files byte-exact. A successor may name a full `original_commit`
when its original source was introduced after the shared registry origin.
The cap third-moment helper uses its `bc70162` origin: its later extension
from quintic to degree-seven histories does not rewrite the two saved
mass-ledger summaries. The summary checker itself has a validation-only
successor entry. Check receipts list `original_producer_paths` and the
current/original context; creation still binds current files. Forty storage
and source controls pass, including the three real tracked summaries with
their ignored raw payloads treated as absent. All three standalone summary
checks also pass with local raw payloads present. Neither check regenerates
scientific data. This fixes the producer-binding failure seen on both
platforms in GitHub run `38977670691`; it was not a raw-file absence error.

An active GitHub PDF check is allowed to finish when later work is pushed;
the concurrency group keeps the newest pending revision for the next run.
The shared synthetic-test budgets count only CPU time
inside their tests, retaining cumulative point and active-CPU limits.

When a publication replay fails on another platform, the checker reports
the rejected numeric leaves or separates wall-coordinate drift from internal
forcing-grid drift. These diagnostics read the failed comparison's operands;
they do not rerun a producer or change acceptance. Investigate those values
before modifying a replay allowance.
Conditioned comparisons retain the actual rejected subtree for that report,
so accepted conditioned differences cannot obscure a later failing value.

The Verify workflow lets an active full-suite run finish while subsequent
commits are saved. It keeps the latest pending run in the same concurrency
group; intermediate pending runs can be replaced, as documented by
[GitHub's concurrency rules](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/control-workflow-concurrency).
A successful older run is evidence for that commit, not for a newer queued one.

## Current work

The next two-level composition uses the same quiet membrane law at both
boundaries. Its tests do not run the full physical record:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nested_slice.py tests/test_nsc_atlas_nested_spacetime.py -q
python scripts/run_python.py scripts/record_nsc_atlas_nested_slice.py --output .build/nsc-atlas-nested-slice-reproduction.json
python scripts/run_python.py scripts/plot_nsc_atlas_nested_slice.py --input .build/nsc-atlas-nested-slice-reproduction.json --output .build/nsc-atlas-nested-slice-reproduction.png
```

Freeze the producer and its listed sources before the application; both
output paths must be new. The twenty controls cover the junctions, bulk
curvature, horizons, units and cap-side detector clock. The complete scalar
record is limited to 30 KB and retains its one 32/64 comparison, source
hashes, inputs and finite domain. The figure regenerates the radial curve
with fixed quadrature and checks its endpoint against the compact record;
it is a radius profile rather than a Euclidean embedding or a causal diagram.
The [owner](nsc-atlas-connection-slice.md#two-successive-sourced-boundaries)
defines the model and central flash. Neither a new cosmological vacuum
value nor a complete two-wall clock atlas is supplied by these commands.

The [saved application](../results/summaries/nsc-atlas-nested-slice-v1.json)
binds method `acc5a176ac074999ed23af770045003c4123cd57`. It stores all
assumptions, milestones, moments and quadrature differences in 11,474 bytes.
Compare its numerical values using the declared tolerances, retaining
runtime and current-commit metadata as reproduction metadata. The flash's
80-digit comparison tests roundoff; it supplies no input/model uncertainty.

The two unequal sourced parent inputs have a compact scalar record:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_connection_slice.py tests/test_nsc_atlas_geometric_average.py -q
python scripts/run_python.py scripts/record_nsc_atlas_unequal_children.py --output .build/nsc-atlas-unequal-children-reproduction.json
```

Use a fresh output path. The [saved result](../results/summaries/nsc-atlas-unequal-children-v1.json)
binds frozen method `f904e0c5184a4a15a05d9d9e1184e8b4dfa46b02` and exact
source hashes. It reuses the first body's saved adaptive moments and
evaluates only the new second body's moments, comparing Gauss 64 with
Gauss 32. There are 42 focused checks. The application takes 0.078125 CPU
seconds / 13,932 pointwise geometric samples under its 30-second / 100,000
sample limit; runtime can differ on another computer. No field solve,
evolution, parameter scan or connector solve is part of this command.
Compare construction identities in ell units with absolute `1e-8`, and the
quadrature moments with `1e-7*max(1,abs(moment_64))`; preserve the reported
32/64 differences as indicators. Runtime and current-commit metadata are
not numerical predictions. The exact source bytes, assumptions and finite
domain must agree. A future changed method requires its own source version.

The [shadow-radial comparison](nsc-atlas-finite-amplitude.md#separating-radial-quadrature-from-reconstructed-field-transport)
has exact controls for its local identity, wall conversion, differentiated
panel-error sum and moving-interface signs:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_mass_transport_identity.py -q -s
```

These controls share 20 active process-CPU seconds and take no numerical
field samples. They do not evaluate a retained cap or authorize a new
trajectory. The numerical shadow-flow/error evaluator is still to be built.
Its selected reference uses the smooth analytic potential; the affine panel
component and cap composition are now implemented for manufactured controls.

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_radial_reference.py -q -s
```

These ten manufactured controls compare the potential and its small implementation
gap with 120-digit arithmetic, check exact moving affine maps and an independent
fixed-radius variation, and exercise the born-panel and budget failure paths.
They share 200,000 charged values / 20 active CPU seconds and use no retained
cap. The returned ODE tolerances and Taylor remainder estimates are separate
quantities; neither asserts a physical or uniform integration-error bound.
Development uses 9,356 charged values / 0.03125 active CPU seconds.

The read-only field adapter supplies fixed-radius and radial derivatives
of the same reconstructed fields, with explicit owner-cell traces:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_fields.py -q -s
```

Its 27 manufactured controls pass. Cumulative development uses
53,753 charged values / 0.734375 active CPU seconds of the
200,000-point / 20-second allowance. No retained cap or radial solve is used.
Native samples exclude axis jets and exact formula-switch queries. Reference
samples provide explicitly named small-radius limits, and a separate method
returns regular axis coefficients without claiming axis field jets.

The cap composer retains the producer's native panels and splits only the
reference at the fixed formula switch. Local map errors and their derivatives
are summed independently of the final mass/rate subtraction:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_radial_comparison.py -q -s
```

Seven controls check an exact constant-field cap, a nonconstant reference
against independent fixed-radius variation, moving switch coincidences,
frozen partitions and invalid native evidence. Cumulative development uses
691,326 charged values / 9.40625 active CPU seconds of the 1,000,000-point /
40-second allowance. A native born panel exactly on a formula switch remains
outside the existing tangent's domain. The terminal comparison below combines
these maps with the independent lapse and residual/jump integrals.

The next component integrates mass, lapse and the three field-residual channels
on one declared smooth piece, with no input from the native mass tangent:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_reference.py -q -s
```

Five manufactured checks pass using 9,828 charged values / 0.046875 active CPU
seconds of a 200,000-point / 20-second allowance. Exact polynomial integrals
and a 60-digit independent mass integral check the states and attenuation.
The full assembler below supplies the interface assembly. The expanded seven
transport-algebra controls also check wall/centre coincidences and separate
interface velocities; that run uses zero field samples / 0.8125 CPU seconds.
Neither component evaluates the retained cap or provides an error enclosure.

The full cap assembler now connects these components and every ordered
interface event, while keeping the affine and residual-ODE reference solutions
independent:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_transport_balance.py -q -s
```

Ten manufactured controls pass. Cumulative development, including failed
and corrected attempts, uses 1,327,411 charged values / 18.25 active CPU
seconds of a 2,000,000-point / 60-second allowance. One coupled-rate test uses
the original RHS and discrete ledger; the new channel sum differs by about
`8.1e-16`. The full-balance wall remains off the native formula switches, and
the centre lapse is independently checked against the frozen field. Reported
reference mass/rate/jump discrepancies still require a declared accuracy
budget before interpreting a retained-state result. That next fixed terminal
comparison is specified in the calculation owner. Its first attempt stopped
before integration at a limiter-branch mismatch. A bounded diagnosis identified
component-order cancellation; the corrected adapter now selects the native
branch before contracting reconstruction derivatives. The exact guard remains.
The [completed fixed comparison](current-result.md#terminal-cap-transport-attribution-16-october-2026)
passes the declared accuracy gates in both cases.
Creation plus reconstruction uses 23,516,078 points / 185.40625 CPU seconds;
including the initial failure and branch diagnosis gives 25,382,986 /
190.84375. The standalone reconstruction passes with 11,758,039 / 93.265625.

The retained-state recorder and its synthetic routing controls are:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_cap_transport_balance.py -q -s
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_transport_balance.py --output .build/results/nsc-atlas-cap-transport-balance-v1.json --summary results/summaries/nsc-atlas-cap-transport-balance-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_transport_balance.py --output .build/results/nsc-atlas-cap-transport-balance-v1.json --check
python scripts/run_python.py scripts/summarize_result.py --output results/summaries/nsc-atlas-cap-transport-balance-v1.json --check
```

Creation reconstructs the analysis once for independent validation. The
[fixed comparison](prior-art-and-open-claim.md#next-fixed-terminal-comparison)
sets its shared 48-million-point / 300-second allowance and the separate
24-million-point / 150-second check. Only terminal cut 128, 64 cells and native
order eight are used. Scalar summaries select `/analysis/summary`; full
piece arrays stay local. Twenty-one manufactured recorder controls pass.
Actual input and producer hashes remain exact. Numerical replay uses absolute
`1e-8` / relative `2e-12`, with all original gates also checked independently
on the saved detailed data. Adaptive ODE counts are bounded work receipts;
the two-metric/two-RHS/two-tangent/four-balance inventory remains exact.

The fixed physical-coordinate mass and initial-representation identities:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_fixed_space_identity.py -q -s
```

Four controls pass once with 20,000 charged algebra values / 0.265625 CPU
seconds of 60,000 / 20. They integrate the quadratic cap mass independently,
check a basis not closed under dilation, distinguish coefficient and moment
rates, and reconstruct an arbitrary local quintic exactly. The [fixed-space
owner](nsc-atlas-finite-amplitude.md#fixed-physical-coordinate-bulk-space)
records why this is the next bulk route and which numerical/physical checks
remain. No field samples, time stepping or raw files are produced.

The fixed local basis, exact initial transfer and independent comparison:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_fixed_basis.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_cap_fixed_transfer.py
```

The thirteen basis/transfer controls use 197,250 charged values / 0.28125
active CPU seconds across fifteen development executions, within 2,000,000 / 30.
The second command is a separately admitted single retained-input check
with 3,000,000 values / 30 CPU seconds. It authenticates the existing payload,
transfers the original exact-binary Hermite data, checks the polynomial
identities through independent Bernstein/de Casteljau evaluation, and
factors one clipped Gram matrix. It prints compact hashes and scalar
evidence without generating raw arrays. It performs no radial solve or
time evolution. The [calculation record](prior-art-and-open-claim.md#one-retained-fixed-basis-transfer-check)
separates coefficient rounding, sampled numerical consistency and the
remaining error analysis.
The retained check passes once with 936,018 charged values / 2.484375 CPU
seconds. Its [compact receipt](../results/summaries/nsc-atlas-cap-fixed-transfer-v1.json)
records exact identity counts, sampled errors and the scaled Gram condition.

The fixed-field metric-domain controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_fixed_bound.py -q -s
```

Five controls pass with 17,600 charged values of 500,000 / 20 active CPU
seconds; measured active time is below the reported timer granularity.
They check the exact Bernstein expansion and source envelope on existing
cells. This command does not run a retained-input certificate or radial ODE.

The separately admitted retained application:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_cap_fixed_domain.py
```

It passes once with 278,272 charged values / 2.296875 CPU seconds of
1,000,000 / 15. Its compact [receipt](../results/summaries/nsc-atlas-cap-fixed-domain-v1.json)
records the positive metric floor for the rounded coefficients. This
regenerates the initial transfer and performs one certificate, without
rerunning the polynomial comparison, Gram factor or a radial integration.

The fixed-space mass/gradient controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_fixed_reference.py -q -s
```

Six manufactured controls check independently integrated alpha=0 energy
and rates, finite coefficient/shape derivatives, stable centre limits and
the strong/weak wall distinction. The [admission record](prior-art-and-open-claim.md#fixed-local-field-mass-gradient-adapter)
preserves the initial conditioned-rate failure, one tighter-accuracy repeat,
and the subsequent physical-norm diagnosis. Aggregate development work is
5,226,742 values / 4.640625 active CPU seconds, including one standalone
diagnostic; a fresh suite has a 5,000,000 / 60 allowance. A separate
1,000,000 / 15 follow-up checked the stable Gram norm. Numerical tolerances
and same-metric forward/adjoint agreement are not uniform error bounds.

The separately admitted single retained reference command is
`python scripts/run_python.py scripts/check_nsc_atlas_cap_fixed_reference.py`.
Its first 25,000,000-value / 120-CPU-second attempt reaches adjoint
quadrature, then exhausts at 25,000,008 values / 7.515625 CPU seconds.
The forward integrations are complete, but there is no checked retained
mass/RHS. The owner bounds the adjoint/operator tail explicitly and admits
one completion replay with 30,000,000 / 120, at identical numerical settings.
No scan or time trajectory is involved.
The one completion replay passes with 25,463,040 values / 7.734375 CPU
seconds. Both retained attempts together use 50,463,048 values / 15.25 CPU
seconds. The [compact result](../results/summaries/nsc-atlas-cap-fixed-reference-v1.json)
records mass, lapse, weak balance, actual wall differences and propagated
gradient comparison. The script now uses the admitted 30,000,000 / 120
allowance. No further replay is required for this recorded result.

The moving-material-wall connection has exact algebra controls and a
separate initial-event audit:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_hybrid_compatibility.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_cap_initial_junction.py
```

The six controls use 40,000 algebra values / 0.484375 active CPU seconds
across their original executions. They distinguish weak, actual and matched
wall traces and include the physical mass dictionary. The audit uses 2,180
values / 2.1875 CPU seconds of 10,000 / 15. It authenticates the retained
input and cap receipt, records the terminal common-data differences, and
checks the existing quiet junction algebra. It runs no radial/prefix solve
or time step. The [compact audit](../results/summaries/nsc-atlas-cap-initial-junction-v1.json)
is saved; these commands are reproduction entry points, not an instruction
to repeat completed calculations before the next connection.

The shared-value instantaneous hybrid adapter has focused controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_fixed_hybrid.py -q -s
```

All four pass once with 375,431 charged values / 0.375 active CPU seconds
of the shared 1,500,000 / 45 allowance. The independent zero-coupling energy
oracle is included on the small test space. These controls call no stored
physical trajectory and take no time step. Historical cap-reference
producer hashes are authenticated at commit `8f6c73c`; current code changes
do not update the original numerical receipt.

The common fixed-cap/exterior stage and its zero-coupling time estimates:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_fixed_cap_stage.py -q -s
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_time_operator.py -q -s
```

Four stage controls pass, plus one focused repeat for the output names,
using 250,248 values / 0.265625 active CPU seconds of 2,000,000 / 60.
The five exact time controls use 76,000 algebra values / 7.984375 active
CPU seconds of 80,000 / 20 across their original executions. They verify
the frozen generator, an energy-norm time bound and the moving-cut energy
identity. These commands do not run a physical trajectory or a grid scan.
No unchanged suite repeat is required for the next connection.

The fixed-coefficient time adapter and original-source quiet-front audit:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_fixed_cap_time.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_quiet_source_front.py
```

The four time controls and one focused extension use 531,005 values /
0.53125 active CPU seconds in development. They include one manufactured
nonzero coupled step and an independently evaluated dense-path residual.
The source command prints a compact audit, authenticating the original
129-knot Hermite input without another source integration. Its saved
[receipt](../results/summaries/nsc-atlas-quiet-source-front-v1.json) records
the completed run and the first failed comparison's upper-bound cost.
The [owner](nsc-atlas-finite-amplitude.md#retained-source-at-the-first-scalar-front)
explains why the earlier quiet knot, rather than the scalar front, supplies
the candidate numerical starting values. No full retained trajectory or
uniform time-error bound is produced by either command.

The retained-source onset's exact time scale and boundary-order controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_linear_bound.py -q -s
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_onset_lag.py -q -s
python scripts/run_python.py -m pytest tests/test_nsc_atlas_hybrid_onset.py -q -s
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_characteristics.py -q -s
```

The two matrix controls use 38,144 values / 0.03125 active CPU seconds;
the two onset-order controls use 24,000 / 0.40625. The three source/hybrid
and incoming-cap controls use 38,000 / 1.28125. The six finite-coupling
characteristic/port controls pass in seven executions using 27,948 /
0.234375, including the affected port repeat after stable frame evaluation.
None replays fields.

The incoming-history extension of the same fixed radial mass adapter:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_fixed_background_reference.py -q -s
```

The four controls check total-field equivalence, independent explicit-time
mass sensitivity, the incoming quartic mass and common-wall flux, the
finite-coupling remainder, quiet/cancelling fields and domain rejection.
Including one execution of the six existing fixed-reference controls and
targeted repairs, development used 4,569,376 values / 5.28125 active CPU
seconds. The [owner](nsc-atlas-finite-amplitude.md#incoming-history-with-fixed-departures-and-an-actual-mass-derivative)
keeps numerical resolution distinct from the exact-input metric bound.
This command uses manufactured histories, without a physical source or
time-trajectory replay.

The source-history onset obstruction and its exact first-piece solution:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_history_onset.py -q -s
```

All three symbolic controls pass with zero field samples / 0.21875 active
CPU seconds. They test onset order and coefficient, not a bound on the
prepared-history interpolation error or a physical evolution.

The implemented first active piece, centre moments and mass/metric checks:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_retained_onset.py -q -s
```

Eight new and 65 affected existing controls pass across 75 development
executions with 787,030 charged values / 2.0625 active CPU seconds. The
existing suites cover characteristic background, retained background jets,
background bounds and cap mass derivatives. No physical source is replayed.

The prepared-source adapter's manufactured interface checks:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_prepared_background.py tests/test_nsc_atlas_background_profile.py -q -s
```

The 33 controls pass with 8,596 charged values / 0.09375 CPU seconds.
The separately admitted physical check is
`python scripts/run_python.py scripts/check_nsc_atlas_prepared_history.py --output results/summaries/nsc-atlas-prepared-history-v1.json`.
It refuses to overwrite a result, shares one 10-million-value / 60-second
budget, and saves only summaries and hashes, including a stopped outcome.
The retained v1 source completes but fails the exact quiet-jet check; it
must not be overwritten. After the causal quiet-branch correction, the
same command with output
`results/summaries/nsc-atlas-prepared-history-v2.json` completed 392
comparisons on 393 retained knots, with 1,362,327 values / 2.3125 CPU
seconds. The [receipt](../results/summaries/nsc-atlas-prepared-history-v2.json)
binds the source commit, inputs, comparisons and finite interpretation.
Both outputs are immutable; use a new output name for a justified replay.
This comparison is not a continuum error enclosure or a cap trajectory.

The retained background in the common stage has separate finite controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_fixed_background_stage.py tests/test_nsc_atlas_fixed_cap_stage.py tests/test_nsc_atlas_fixed_cap_time.py -q -s
```

Eleven controls pass across thirteen development executions with 2,375,777
values / 3.6875 active CPU seconds. They compare the same total polynomial,
its constrained mass direction and one common manufactured step at the
evolving cap clock, and verify invariance under a constant history-clock
translation. They do not integrate the prepared physical source. For the
prepared history, pass its inherited `physical_input_domain.Ta` as
`background_clock_offset`; the reference wall gauge uses T=0 at that anchor.

The exact quiet-prefix integration simplification is checked separately:
`python scripts/run_python.py -m pytest tests/test_nsc_atlas_quiet_history_partition.py -q -s`.
One control passes with 386,961 values / 0.609375 CPU seconds, using exact
rational projection integrals and an unchanged original-cut certificate.

The numerical attempt was produced with
`python scripts/run_python.py scripts/run_nsc_atlas_fixed_onset.py --output results/summaries/nsc-atlas-fixed-onset-v1.json`.
It stopped at the 180-second CPU check after 341/371 steps. Preserve that
receipt; any future rerun needs its own output and an admitted reason.
The [owner](nsc-atlas-finite-amplitude.md#first-retained-source-fixed-cap-interval)
records the residuals and the analytical next action. No replay is needed
to verify the saved bytes and producer hashes.

The prepared-history successor uses the same full interval and time
schedule, with a separately admitted budget:

```sh
python scripts/run_python.py scripts/run_nsc_atlas_fixed_onset.py --retained-history --output results/summaries/nsc-atlas-fixed-history-onset-v1.json
```

This is a creation command, not part of routine verification. The runner
requires committed producer sources and refuses to overwrite any existing
output. It regenerates the leading field once, retains the inherited clock
origin, and charges preparation and common evolution to the same overall
budget. The original left endpoint is retained because the source metric
interpolant has a pre-front constraint defect. Independent derived-residual
controls compare clock/radius/field perturbations with recomputed masses
and wall scalars; source-failure accounting and receipt preservation are
checked by `tests/test_run_nsc_atlas_fixed_onset.py` without a physical solve.

The [v1 history receipt](../results/summaries/nsc-atlas-fixed-history-onset-v1.json)
completed 371 steps with 1,176,614,645 values / 619.4375 CPU seconds. Its
producer is frozen at `2551c9f52703fb418bec5a5ebf3d0bde136cfcc5`; 70 source
hashes, the original input/geometry/schedule and the v2 prepared-history
identity have been authenticated without a replay. The committed receipt
is 18,513 bytes and contains no full history. Preserve it when regenerating
into a separately named local output. The one-sided source endpoint check
uses the saved input's `exact_rows[::2]`, knot 36, and the actual left/right
Hermite curvatures; its equations and values are in the
[owner](nsc-atlas-finite-amplitude.md#derived-time-residuals-and-the-full-retained-history-interval).

The subsequent diagnostic extension separates the wall residual into
weak-velocity and incoming-characteristic terms. Two exact controls use
14,000 algebra values / 0.390625 CPU seconds; run them with
`python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_hybrid_compatibility.py -k 'velocity_and_incoming or mass_sensitivity' -q -s`.
Future onset outputs include one small final state and its wall/gradient
diagnostics so a local attribution need not reconstruct the whole history.
That extension does not backfill or change the frozen v1 result. These algebraic controls
perform no further physical run.

The separately admitted single-event alpha-zero velocity control uses the
completed endpoint geometry and the same prepared history with zero cap
departures:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_fixed_history_endpoint.py --output results/summaries/nsc-atlas-fixed-history-endpoint-velocity-v1.json
```

This creation-only command shares its 10-million-value / 60-second budget
between preparation and one radial reference. It takes no time step and
does not reconstruct the completed finite-alpha endpoint state.

The separate leading-reference endpoint comparison checks the prepared
wall jets, inherited arrival clock, retained profile and one existing
common stage:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_prepared_endpoint_connection.py --output results/summaries/nsc-atlas-prepared-endpoint-connection-v1.json
```

This creation-only comparison has its own admitted 10-million-value /
60-second budget. It samples 63 original exterior interior positions from
one prepared field, takes no time step and reports an alpha-zero reference
event with zero cap departures. It does not recreate the evolved endpoint.

The separately admitted local dense-polynomial split uses the same event:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_prepared_endpoint_connection.py --audit-dense-piece --output results/summaries/nsc-atlas-prepared-endpoint-connection-v2.json
```

It analytically differentiates one existing dense polynomial and records
its difference from the prepared RHS jet and retained history derivative.
It retains those two small polynomial pieces, not a trajectory. The
physical derivative interface and evolution operator remain unchanged.

The first enriched reference comparison uses the unchanged source and knot
grid, with additional ordinary-piece value/first jets:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_prepared_endpoint_connection.py --enrich-history --audit-dense-piece --output results/summaries/nsc-atlas-prepared-endpoint-enriched-v1.json
```

It separately authenticates the original nodal history and records the new
input identity. Declared local decision limits are enforced in the receipt.
This command also has a 10-million-value / 60-second budget and takes no
time step. Existing creation-only output files are preserved.

The separately admitted **P001 endpoint comparison** uses three fixed
nonuniform cap meshes and three exterior resolutions. It starts at the quiet
`x=0.4` anchor and targets `x=0.4815042082576364` (reference radius `0.55`),
with the same finite amplitude and source on every mesh. The opening-phase
robustness draft is paused. The resource limits and interpretation are in the
[calculation admission](prior-art-and-open-claim.md#p001-full-endpoint-fixed-grid-comparison-17-october-2026).

From a committed source tree, create new compact receipts with:

```sh
python scripts/run_python.py scripts/run_nsc_atlas_fixed_endpoint.py --level 0 --output results/summaries/nsc-atlas-fixed-endpoint-level0-v4.json
python scripts/run_python.py scripts/run_nsc_atlas_fixed_endpoint.py --level 1 --output results/summaries/nsc-atlas-fixed-endpoint-level1-v4.json
python scripts/run_python.py scripts/run_nsc_atlas_fixed_endpoint.py --level 2 --output results/summaries/nsc-atlas-fixed-endpoint-level2-v4.json
```

The commands are creation-only and never overwrite a saved result. Each
receipt records completion or the reason for stopping, common-cut scalars,
residual summaries and hashes; it does not save full histories. The physical
observable conversion is in `nsc_atlas_endpoint_comparison`, with a fiducial
`ell=10000 metres` fixed before the runs. A three-grid resolution range is
reported separately from common-input and physical-model uncertainty.

The first middle/fine receipts stopped before accepting a step because of a
redundant zero-history panel boundary. Later coarse and middle runs stopped
on the radial gradient check; the old fine run was deliberately interrupted.
The v3 successors use the separately confirmed forward accuracy, preserving
every earlier receipt. The coarse v3 run nevertheless stops later, after
464 accepted steps. The middle v3 also stops. Both saved failure states pass the unchanged
gradient guard with quarter-panel forward radial steps. The v4 commands above
now select the panel's own history polynomial at its endpoints, at unchanged
physical coordinates. The one-sided diagnostics identify this as the source
of the forward/adjoint mismatch. The provisional radial step cap is not used.
Certified exactly quiet data retain their previous zero-field startup. The other solver settings,
physical grids and time schedules stay fixed. Once the three v4 endpoints
complete, compare their saved scalars without another physical evolution:

```sh
python scripts/run_python.py scripts/analyze_nsc_atlas_fixed_endpoint.py --inputs results/summaries/nsc-atlas-fixed-endpoint-level0-v4.json results/summaries/nsc-atlas-fixed-endpoint-level1-v4.json results/summaries/nsc-atlas-fixed-endpoint-level2-v4.json --output results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json
```

The analyzer rejects incomplete endpoints, changed sources and inconsistent
frame conversions. It records both pointwise residual maxima and accumulated
transport gaps. A nonshrinking result remains a nonshrinking result; it does
not become a converged prediction through the physical-unit conversion.

The later CPU-only resource amendment admits ceilings of 20,000 / 40,000 /
80,000 process-CPU seconds for these same three levels. Add
`--extended-cpu-allowance` to a fresh runner command to select those limits.
This changes no physical input, grid, numerical tolerance, charged-value
limit or accepted-step cap. The running v4 jobs retain their original
scientific source and, when amended, record the separate controller source,
old/new CPU limits and cumulative accounting in compact amendment receipts.
Do not replay the live-process amendment against a different process: a clean
reproduction selects its finite CPU allowance at startup.

Render the authenticated saved comparison with:

```sh
python scripts/run_python.py scripts/plot_nsc_atlas_fixed_endpoint.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json --output docs/figures/nsc-atlas-fixed-endpoint-comparison.png
```

The plotter checks receipt hashes and recomputes the saved comparison before
showing the three residual maxima and the wall-radius resolution envelope
in metres. Increasing and unchanged gaps remain visible. It performs no
physical evolution and refuses to overwrite an existing figure. Fixture
renders used to check its layout stay under ignored `.build`.

Inspect the mass-gap decomposition at the same saved cuts with:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_gaps.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json --output results/summaries/nsc-atlas-fixed-endpoint-gap-diagnostic-v1.json
```

This authenticates the completed comparison and its three inputs, then reports
the signed trace and time-projection contributions, arithmetic identity
residuals, and the coverage of the all-stage maxima by the nine saved cuts.
It performs no physical solve or replay. A maximum between those cuts has no
saved state, so this diagnostic cannot establish its cause. Small algebraic
identity residuals alone do not establish physical mass or wall closure.

The first coarse receipt stopped at the radial forward/adjoint check after
419 accepted steps. The bounded replay and subsequent one-stage comparison
have separate, predeclared budgets in the same admission:

To regenerate the original failed-stage seed, use the recorded replay
source commit `eb38441` in a clean checkout. The current runner has a
tighter radial tolerance and is not the original failing implementation.
Keep the generated receipt with its own seed hash; cross-platform
regeneration need not have identical state bytes.

```sh
python scripts/run_python.py scripts/run_nsc_atlas_fixed_endpoint.py --level 0 --diagnostic-step-limit 430 --output results/summaries/nsc-atlas-fixed-endpoint-level0-diagnostic-v2.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-diagnostic-v2.json --output results/summaries/nsc-atlas-fixed-endpoint-radial-diagnostic-v1.json
```

The second command requires the authenticated single-stage NPZ captured by
the first command on a stage failure. This small local seed stays in ignored
`.build/results`; its hash and byte count are recorded in the committed
summary. A clean checkout regenerates it with the bounded replay. The
radial comparison performs no time step and preserves the original
forward/adjoint discrepancy guard.

After the first diagnostic attributed the discrepancy to the forward solve,
the separately admitted two-case confirmation uses:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-diagnostic-v2.json --confirm-forward-accuracy --output results/summaries/nsc-atlas-fixed-endpoint-radial-confirmation-v1.json
```

For the later coarse v3 failure, its seed comes from the full level-0 command
at source commit `6b6e69c`. The later diagnostic commands are:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --confirm-forward-accuracy --output results/summaries/nsc-atlas-fixed-endpoint-radial-later-diagnostic-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --panel-step-check --output results/summaries/nsc-atlas-fixed-endpoint-radial-panel-diagnostic-v1.json
```

The panel comparison holds the floating-point tolerances fixed and restricts
the forward radial step to half/quarter of each existing smooth panel. It
changes neither the physical cap mesh nor the forward/adjoint acceptance
guard, and it performs no time evolution.

The subsequent endpoint-side diagnosis and public-interface confirmations use
the same authenticated local failure seeds, with no trajectory replay:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --one-sided-panel-check --output results/summaries/nsc-atlas-fixed-endpoint-radial-sided-diagnostic-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level1-v3.json --one-sided-panel-check --output results/summaries/nsc-atlas-fixed-endpoint-radial-middle-sided-diagnostic-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --confirm-panel-polynomials --output results/summaries/nsc-atlas-fixed-endpoint-radial-panel-polynomials-coarse-v2.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level1-v3.json --confirm-panel-polynomials --output results/summaries/nsc-atlas-fixed-endpoint-radial-panel-polynomials-middle-v2.json
```

The diagnosis moves only numerical endpoint samples inward by a bounded
roundoff distance. The final method instead uses the panel's own polynomial
at the unchanged query coordinate. Its 32/64-point confirmations enter via
the actual coupled-stage interface, without injecting a solver wrapper.

After that comparison passes, the separately admitted full enriched
finite-coupling onset interval is:

```sh
python scripts/run_python.py scripts/run_nsc_atlas_fixed_onset.py --retained-history --enrich-history --output results/summaries/nsc-atlas-fixed-enriched-history-onset-v1.json
```

The source, geometry and 371-step schedule stay those of the original
retained-history attempt. The 900-second / two-billion-value budget counts
one preparation and all stages. The command checks the committed enriched
reference and its input hashes, and saves one small final state plus
summaries. Its sampled scalar decision is distinct from interval completion
and from a continuum validity claim. Eight orchestration controls pass.
The frozen `42c7741` execution completes in 660.703125 CPU seconds with
1,244,922,621 charged values and passes the sampled scalar decision at
`2.75430e-13`. The 28,312-byte receipt contains the final state, its hash
and component maxima; the 371-step history is represented only by its hash.

GitHub Verify run `38979718935` confirms the storage fix on both Linux and
Windows, including all retained tests and the current-paper checks. Later
research commits have their own queued checks; this completed run applies
to commit `6239ee1`.
The later Verify run `38983430995` at `b422d9b` also passes both complete
Linux/Windows suites and current-paper checks. The PDF workflow at `0a2799c`
passes separately; newer commits retain their own queued verification.

The uniform reference gauge has three independent controls and one bounded
two-integral application on the already certified original source cell:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_reference_gauge_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_reference_gauge_bound.py --output results/summaries/nsc-atlas-reference-gauge-bound-v1.json
```

The application requires committed producers and a new output path. Its
100,000-value / 20-second budget covers two predeclared radius endpoints and
interval coefficient evaluation, with no root search or evolution. The
summary retains signed dyadic enclosures. Exact corner normalization is
explicit; numerical gauge and sourced-tube errors are separate. The controls
pass with 16,000 charged values / 0.09375 active CPU seconds.

Frozen producer `81b0170` completes the application with 5,248 charged
values and CPU below timer resolution. Its 13,822-byte summary proves the
strict reference-radius bracket `(0.474,0.477)` on the original cell.
All 17 producer bindings, signed displays and inverse inequalities pass a
separate 256-value rational audit. Grok independently confirms the formulas
and the distinction between this reference bound and a sourced tube.

The exact initial anchor uses the saved onset event and the preceding gauge
certificate, with no new source solve or trajectory:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_reference_initial_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_reference_initial_bound.py --output results/summaries/nsc-atlas-reference-initial-bound-v1.json
```

The runner authenticates the original source-geometry file as well as the
compact receipts. Three integrals at one saved seed have a shared
100,000-value / 20-second budget. It checks the exact inverse correction,
both clock corrections, quiet source/history margins and original cap cut.
Its four controls pass with 15,004 values / 0.40625 active CPU seconds.
The numerical source curvature defect remains recorded. Future residual
certification must evaluate the specified exact-reference reconstruction.

The frozen `33445f7` application passes with 7,872 charged values and CPU
below timer resolution. It saves a 12,790-byte summary with 43 producer
bindings. A separate 512-value rational audit checks all displays, input
hashes, the strict quiet margins and the original fixed-basis cut. The saved
physical initial radius/clock differences are below `1.125e-16`/`1.244e-16`;
the whole cap remains quiet with margin above `0.0014530219858557`.

The actual source forcing relative to the vacuum reference uses the same
original cell and pulse certificate:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_source_forcing_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_source_forcing_bound.py --output results/summaries/nsc-atlas-source-forcing-bound-v1.json
```

Its four controls pass with 21,192 values / 0.28125 active CPU seconds.
The physical application requires committed producers and a new output path;
its allowance is 4,096 values / 5 seconds. Exact rational calculation with
128-bit outward dyadic rounding replaces no source data and runs no source
integration or evolution. It bounds actual forcing, not interpolation error.

Frozen producer `7025f41` passes with 2,048 charged values and CPU below timer
resolution. The 6,365-byte summary retains exact rational bounds for the
whole admitted coupling interval. Its 19 producer hashes and input identity
authenticate; a 512-value rational audit verifies the displays, reported
ceilings and source/frame continuation margins. Grok independently confirms
the equations and outward rounding.

The fixed conditional coupled-tube check uses the authenticated source,
gauge and initial certificates above:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_coupled_tube_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_coupled_tube_bound.py --output results/summaries/nsc-atlas-coupled-tube-bound-v1.json
```

The four controls pass with 26,192 charged values and active CPU below timer
resolution. A prior control invocation's output was unavailable after
context transfer; its cache retained no failed nodes. Conservatively charge
that invocation's full 40,000-value / 20-second allowance in the development
ledger, for at most 66,192 values / 20 seconds across both invocations.
The physical application is separate: one fixed box evaluation, 8,192
charged values within 20,000 values / 10 CPU seconds, with committed producer
sources and a new output path. There is no source solve, trajectory or box
search. The receipt explicitly distinguishes a conditional a priori bound
from construction of a smooth solution.

Frozen producer `706fe56` returns `strict_conditional_bootstrap_image`:
all ten ratios are below `0.516` and all seven margins are positive.
The application uses 8,192 values and CPU below timer resolution. Its
29,878-byte receipt binds 24 producer files and retains 64 interval records.
The separately admitted 512-value / 0.09375-second rational audit passes
all hashes, exact ratios, outward displays and scope checks. Grok's completed
implementation review confirms the bounds; no physical rerun is needed.

The conditional first-derivative closure has separate symbolic controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_derivative_closure.py -q -s
```

All three pass with 24,000 values / 0.328125 active CPU seconds. They verify
the quotient/kernel differentiation, cap commutators and radial time
variation; they do not evaluate physical coefficients or construct a solution.

For the primitive-preserving null construction map:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_null_construction.py -q -s
```

Four exact controls pass with 34,000 cumulative values / 0.25 active CPU
seconds (the first three together, then only the added fourth control).
They include a fixed source, a spacelike initial curve, old and wall-born
points, the six incoming wall jets' common tangent identities, and a
counterexample to inferring derivative convergence from weak error bounds.
They verify the map's primitives and wall-value transport, not a contraction
theorem; no physical trajectory is calculated.

For the next null-jet and chart identities:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_null_continuation.py -q -s
```

Three exact controls pass with 24,000 charged values / 1.90625 active CPU
seconds. They use one nonlinear manufactured solution, a curved initial
line and moving wall, and independent symbolic chart variables. The checks
cover the second triangular derivative hierarchy, wall entry derivatives,
metric/scalar conversion and proper wall normalization. They do not evolve
the physical preparation or certify full-domain continuation.

The coupled comparison controls are:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_integrated_feedback.py -q -s
```

They pass with 29,000 charged values / 0.0625 active CPU seconds, including
the moving interval's complete boundary history and the production lapse
source in canonical variables. They do not evaluate a physical inverse
constant or run a trajectory.

For the connected scalar inverse:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_transmission_bound.py -q -s
```

Its seven controls pass within 99,024 cumulative development values / 0.796875
active CPU seconds, including one corrected clock-conversion control. It
performs no physical source solve or evolution. The comparison helper requires
proved coefficient/error bounds and does not certify them from samples.
The source-conversion helper maps invariant value/first-jet bounds into the
canonical scalar bounds using an explicitly supplied common source tube.
The admitted application authenticates the existing v2 certificate and its
proved tube, then uses exact rational arithmetic only:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_canonical_source_bound.py --output results/summaries/nsc-atlas-canonical-source-bound-v1.json
```

It requires committed sources and a new output path, with a 4,096-value /
5-second budget. Display floats are summaries; the rational strings carry
the error bounds. This command performs no source readout or evolution.
The frozen `d365a13` application completes with 1,024 charged values and CPU
below timer resolution, after a pre-calculation source-graph correction. Its
6,437-byte summary retains exact rational bounds; a 128-value rational audit
checks the outward displays `4.658e-38` and `2.576e-33`. All 13 producer
bindings authenticate.

The source representation certificate uses the same source cell and
its actual binary cubic coefficients:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_source_cell_bound.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_source_cell_bound.py --output results/summaries/nsc-atlas-source-cell-bound-v2.json
```

The seven controls pass. The separately admitted physical certificate has a
100,000-value / 30-second budget and requires committed calculation sources.
It encloses two original panels without sampling, source integration or
time evolution, and refuses to overwrite an existing receipt. Its frame
bounds concern the mathematical representation, with runtime floating
evaluation and coupled propagation errors still separate.
The frozen `ac98a98` v1 receipt bounds the closed left cubic. The v2 command
also includes the actual next-piece value and first jet selected at its
right endpoint; their coefficient-rounding jumps are added explicitly.
The `191a05b` v2 execution completes in 10,372 values / 2.5625 CPU seconds;
its 14,751-byte receipt retains all six frame bounds as outward dyadics.
The reported zero source readout calls exclude the constructor's 132
existing-knot/geometry validation values, which are included in the budget.

The internal wavefront identities and bounded position/amplitude stage:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_wavefront_identity.py tests/test_nsc_atlas_cap_wavefront.py -q -s
```

Five exact identities pass across six development executions, using 36,000
charged algebra values / 0.390625 CPU seconds of 50,000 / 20. The first
factored-versus-expanded assertion was corrected and rerun in isolation.
Five stage controls pass across six executions with 404 charged values /
0.015625 CPU seconds of 50,000 / 15; the focused repeat checks scalar-real
orientation rejection. No radial solve or numerical
time trajectory is run. The [owner](nsc-atlas-finite-amplitude.md#internal-wavefront-compatibility-and-evolution)
states the exact two-front interval, source-free interface hypotheses and
the distinction between this stage and a coupled cap evolution. A full
algebra-suite repeat exceeds its remaining original development allowance
and needs its own justified reproduction account.

The moving tangent-enrichment regularity control:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_background_enrichment.py -q -s
```

Three exact algebra controls pass across four development executions,
using 24,000 charged values / 0.765625 CPU seconds of a 50,000 / 15
allowance. They check the moving-interface mass term, loss of H1 time
differentiability and the fixed-space tangent error. They run no physical
trajectory and create no data file. The [owner](nsc-atlas-finite-amplitude.md#tangent-enrichment-and-its-moving-knot-obstruction)
records the first symbolic-form assertion correction and the scope.

The retained-background field jets, metric certificate and explicit-time
reference controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_background_jets.py tests/test_nsc_atlas_cap_background_bound.py tests/test_nsc_atlas_cap_background_reference.py -q -s
```

Fourteen jet controls pass in fifteen development executions, using 9,430
charged values / 0.34375 CPU seconds of 300,000 / 20. Seven certificate
controls use 7,213 / 0.203125 of 500,000 / 30. Seven reference controls use
378,959 / 1.0 of 4,000,000 / 60. The [admission record](prior-art-and-open-claim.md#retained-background-jets-and-local-metric-domain)
keeps these separate from historical controls. No scientific raw files are
created by these tests. After independent source review, the fixed-profile
certificate and two-reference commands have each run once:

```sh
python scripts/run_python.py scripts/check_nsc_atlas_cap_background_domain.py
python scripts/run_python.py scripts/check_nsc_atlas_cap_background_reference.py
```

The [fixed-input admission](prior-art-and-open-claim.md#fixed-retained-pulse-time-remainder-comparison)
allows one authenticated certificate and then two instantaneous references
under separate allowances. Both commands print compact values and hashes;
they create no files, re-fit no source and run no physical time evolution.
The certificate uses 78,388 charged values / 2.15625 CPU seconds of
100,000 / 10; the two references use 15,264,426 / 7.984375 of
20,000,000 / 120. The [recorded outcomes](nsc-atlas-finite-amplitude.md#fixed-retained-pulse-mass-remainder-diagnosis)
include the exact-input metric floor and the finite coupling remainder.
A full repeat would exceed the remaining original allowances; a justified
reproduction receives its own bounded work account.

The fixed background-profile projection diagnostic:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_profile_projection.py -q -s
python scripts/run_python.py scripts/check_nsc_atlas_cap_profile_projection.py
```

Eight controls pass with 234,742 charged values / 0.90625 CPU seconds of
their 2-million / 30-second allowance. The one fixed input projection uses
2,328,407 / 2.296875 of a separate 4-million / 30-second allowance. It
prints compact results and hashes; it creates no files and does not replay
an evolution. The existing zero-case terminal profile is not adequately
represented by sixteen global modes. The [result](nsc-atlas-finite-amplitude.md#fixed-profile-representation-diagnosis)
records value, derivative, wall and centre discrepancies. Keys ending in
`_xi` label derivative integrands, but returned norms include the coordinate
factor and are physical z-derivative norms.

The explicit-time mass connection for a retained background:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_affine_background_identity.py -q -s
```

Three exact controls pass with 6,000 algebra values / 0.21875 CPU seconds
of a 10,000 / 10-second allowance. No numerical fields are evaluated.
The [remainder and counterexample](nsc-atlas-finite-amplitude.md#mass-balance-with-a-retained-characteristic-background)
show the compatibility term that naive background subtraction would miss.

The orthogonal basis and its shared cap adapter:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_polynomial_basis.py tests/test_nsc_atlas_cap_legendre_reference.py -q -s
```

Twelve basis controls pass in thirteen development executions, using 34,140
charged values / 1.390625 CPU seconds. Five adapter controls use 1,835,038 /
0.6875 and compare the same field across bases, the exact quadratic limit
and one fixed sixteen-mode nonlinear reference. The [separate allowances](prior-art-and-open-claim.md#orthogonal-basis-cap-integration)
do not admit a grid scan or physical time evolution.

The numerical odd-polynomial reference and affected common-engine checks:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_variational_reference.py tests/test_nsc_atlas_cap_projection_bound.py tests/test_nsc_atlas_cap_polynomial_reference.py -q -s
```

All 58 controls pass. Monomial-reference development uses 1,252,570 charged values /
2.84375 active CPU seconds across 54 executions; the separately admitted extraction regressions
use 490,356 / 1.15625 for P1 and 4,620 / 0.5 for its projection bound.
The [record](prior-art-and-open-claim.md#numerical-regular-polynomial-cap-reference)
keeps those allowances separate from earlier development. One moving-domain
RHS and its radius derivative are evaluated; no physical evolution occurs.

The regular odd-polynomial space has an exact moving-domain mass identity:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_dilation_identity.py -q -s
```

Six controls pass in seven development executions, using zero field samples
and 1.296875 active CPU seconds, within the 20-second algebra allowance.
The independent quadratic check and general gradient pairing supply no
numerical radial or physical evolution result.

The moving-cap reference identity and its timelike counterexample are exact:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_moving_projection_identity.py -q -s
```

Seven controls pass using no field samples, radial solve or physical
evolution. The [record](prior-art-and-open-claim.md#moving-cap-geometric-transport-identity)
documents the measured algebra cost. These controls expose the additional
mass term in a proposed coordinate-advection rule; no solver is replaced.

The exact wall-momentum projection bound uses no radial ODE or retained data:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_projection_bound.py -q -s
```

Thirteen controls check the interpolation constant, mass-matrix margins,
constraint derivatives, affine-element normal-field bounds and independent
high-precision enclosure values. The [development receipt](prior-art-and-open-claim.md#exact-galerkin-wall-momentum-projection-bound)
records 8,872 charged values and the incomplete CPU receipt. Numerical radial
error and normal-trace compatibility remain separate from this bound.

The finite-backreaction fixed-cap reference evaluates radial mass gradients
and one weak RHS on manufactured primitive-conforming fields:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_variational_reference.py -q -s
```

Nineteen controls pass with 1,704,332 charged coefficient/derivative values /
3.859375 active CPU seconds across 71 executions, within 2 million / 60.
Forward sensitivities and attenuation-weighted quadrature must agree before
the weak RHS is constructed. The exact-model all-radius Arb bound is distinct
from numerical ODE error; weak and reconstructed boundary traces are reported
separately. No physical history or retained-field comparison is run.

The regular-centre fixed-radius reference implements one zero-backreaction
weak RHS and its mass/physical-boundary balance:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_geometric_reference.py -q -s
```

Twenty-six controls pass, including independent element matrices, direct
radial density integration, exact incoming data, symbolic regularity and
small-radius arithmetic. Cumulative development uses 10,600 charged values /
0.53125 active CPU seconds in 97 executions, within 200,000 / 20. The
component exposes the weak-versus-reconstructed wall discrepancy. It runs no
retained data, radial ODE or physical time evolution.

The mass-gradient and geometric-energy identities use no numerical fields,
reference solves or retained states:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_flux_projection_identity.py -q -s
```

Fifteen exact controls check the adjoint tail, reaction cancellation, mixed
reconstruction derivatives, upwind projection, geometric-coordinate energy,
moving cells and both boundaries. Development used zero field samples /
2.984375 active CPU seconds across 39 executions, within 20 seconds. These
are mathematical controls for the next discretization, not a changed native
evolution or another observational prediction.

The interface algebra uses no numerical fields or root search:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_interface_identity.py -q -s
```

All three exact controls pass in 1.71875 active CPU seconds within their
20-second allowance. They establish primitive continuity, formula-switch
orders and the necessity of the moving potential-jump term.

The discrete cap mass derivative and departure-flux ledger are implemented
in a diagnostic companion to the unchanged producer. Run their controls with:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_mass_derivatives.py tests/test_nsc_atlas_cap_mass_tangent.py tests/test_nsc_atlas_cap_mass_balance.py -q -s
```

The [mathematical owner](nsc-atlas-finite-amplitude.md#discrete-mass-derivative-and-departure-flux-balance)
states the branch-specific derivative and telescoping identity. The controls
have passed. The [completed four-case result](current-result.md#discrete-cap-mass-balance-16-october-2026)
accounts for the upwind contribution and retains the remaining discrepancy.
The [calculation record](prior-art-and-open-claim.md#discrete-cap-mass-balance-16-october-2026)
preserves the failed rounding check and its bounded diagnosis.

Full payloads belong in ignored `.build/results`. Creation refuses to
overwrite an existing output. Reconstruct a missing payload, or validate
an existing one, using the frozen producer and input bindings:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_balance.py --output .build/results/nsc-atlas-cap-mass-balance-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_balance.py --output .build/results/nsc-atlas-cap-mass-balance-v1.json --check
python scripts/run_python.py scripts/summarize_result.py --output results/summaries/nsc-atlas-cap-mass-balance-v1.json --check
```

The last command checks hashes and structure without a scientific replay.
The compact summary contains the exact regeneration argument vector.

The [first prediction-register entry](prediction-register.md#nsc-p001--perturbative-cap-and-moving-wall)
binds the closed perturbative cap–wall relation and four computed
coefficients to their original evidence. Authenticate its projections without
repeating a calculation:

```sh
python scripts/run_python.py scripts/check_prediction_register.py
```

The [terminal cap mass-functional derivative](current-result.md#terminal-cap-mass-functional-derivative-16-october-2026)
uses the [analyzer](../scripts/analyze_nsc_atlas_cap_mass_direction.py) and
[recorder](../scripts/diagnose_nsc_atlas_cap_mass_direction.py) on the authenticated
terminal cap states. Its 26 mathematical, orchestration and provenance controls
pass with cumulative usage of 171,957 points and 1.34375 process CPU seconds,
including the initial epsilon-tolerance failure and its correction:

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_cap_mass_direction.py tests/test_nsc_atlas_cap_mass_direction_orchestration.py tests/test_diagnose_nsc_atlas_cap_mass_direction.py -q
```

The [compact summary](../results/summaries/nsc-atlas-cap-mass-direction-v1.json)
uses 1,711,610 points and 6.6875 process CPU seconds including creation
validation. One standalone check passes with 855,805 points and 4.078125
CPU seconds. The 520,080-byte record has SHA256
`d0fec7048ab1cb476e15672ccf334c2c8e628c059c1fdec709e702c9dafd6f74`.
The original producer was frozen at `fc71709`. Full data remains local in
ignored `.build/results/`. A clean checkout can regenerate it when needed:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_mass_direction.py --output .build/results/nsc-atlas-cap-mass-direction-v1.json
```

Creation refuses to overwrite an existing record. Append `--check` to check
an existing raw record. For an inexpensive summary/input authentication that
does not repeat the numerical diagnostic, use:

```sh
python scripts/run_python.py scripts/summarize_result.py --output results/summaries/nsc-atlas-cap-mass-direction-v1.json --check
python scripts/run_python.py scripts/check_result_storage.py
```

The first command reports whether the raw file is present; its absence on
a clean checkout is expected. [Storage policy](nsc-storage.md#current-calculation-outputs-16-october-2026)
keeps future arrays local and their summaries, hashes and commands in Git.

Each invocation pins the prior cap audit, verifies its seven producer bindings
within the twelve audit bindings, and authenticates the 104 inherited bindings
once. Each analysis uses 20 cap metrics and four cap RHS evaluations at the
terminal cut, holding the directions fixed for backward probes. Creation
performs one analysis and one validation under a shared 60-CPU-second,
four-million-point allowance; a check performs one analysis under 30 CPU
seconds/two million points. No source, full stage or trajectory is replayed.
The three control suites share 30 CPU seconds/300,000 points. Fresh creation
requires a new output file without `--check`; first failures leave an ignored
`.build` receipt without retry.

The finest order-eight finite directional residual is `-0.6612504470`, compared
with `-0.6407403055` from saved history expressed with the same rates; its
order-eight minus order-four change is `1.09497e-5`. These are qualified
backward functional diagnostics, without a convergence-order claim. At zero
amplitude the normalized residual is not a physical gravitational defect.

The [selected cap mass and wall-curve audit](current-result.md#selected-cap-mass-and-wall-curve-audit-16-october-2026)
uses the [recorder](../scripts/diagnose_nsc_atlas_cap_consistency.py) and
[analyzer](../scripts/analyze_nsc_atlas_cap_consistency.py) on the fixed
operator-comparison histories. Its 34 manufactured and provenance controls
pass on their first attempt with 62,119 shared points and 0.34375 process
CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_cap_consistency.py tests/test_nsc_atlas_cap_consistency_orchestration.py tests/test_diagnose_nsc_atlas_cap_consistency.py -q
```

The [completed record](../results/development/nsc-atlas-cap-consistency-v1.json)
uses 4,393,156 points and 15.359375 process CPU seconds including creation
validation. One standalone check passes with 2,196,578 points and 8.453125
CPU seconds. The 820,093-byte record has SHA256
`79d1e02917ba7c1a7d3326cd80612b0467f19b0b5b7b1e69b6e1e2f1f3ae9728`.
To authenticate it:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_cap_consistency.py --output results/development/nsc-atlas-cap-consistency-v1.json --check
```

Each invocation authenticates the 104 inherited bindings once and checks
seven audit source/control bindings. Each analysis reconstructs 64 cap
metrics and 12 selected exterior prefixes and hybrid wall maps, comparing
quadrature orders 4/8 and radial/time widths one/two at cuts 63, 68 and 128.
Creation performs one analysis and one validation under a shared 120-CPU-second,
eight-million-point allowance; a standalone check performs one analysis under
60 CPU seconds/four million points. Neither invokes a source solve, full
coupled stage, cap RHS or trajectory replay. The three control suites share
30 CPU seconds/one million points. A fresh creation requires a new output
file without `--check`; the first failure leaves an ignored `.build` receipt
with its cut, order and phase when available, without retry or completed output.

Pulse residuals remain strongly sensitive to radial derivative width despite
negligible order-4/8 changes. At the final finite cut the temporal mass
residual is about `-8.32e-5`, while radial mass and lapse residuals are below
`1.7e-9` and `9.3e-9`. The actual wall-curve acceleration also differs from
the hybrid prediction. The terminal mass-functional calculation above now
tests the cap transport connection; neither audit supplies a new wall law or
a new evolution.

The [near-wall saved-row audit](current-result.md#saved-near-wall-row-audit-16-october-2026)
uses the [bounded diagnostic](../scripts/diagnose_nsc_atlas_boundary_rows.py)
and [eight synthetic controls](../tests/test_diagnose_nsc_atlas_boundary_rows.py).
Controls pass with 60,460 aggregate points/0.15625 CPU seconds. Creation
including validation uses 97,350 points/1.59375 CPU seconds; its one
standalone check uses 48,675 points/1.5625 CPU seconds. Each authenticates
the 104 input bindings once and evaluates three fixed cuts, without a
stage or trajectory. Both production allowances are 30 CPU seconds and
150,000 points. To authenticate the retained result:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_boundary_rows.py --output results/development/nsc-atlas-boundary-rows-audit-v1.json --check
```

Fresh reproduction selects a new output file without `--check`; existing
files are protected. The result reports stencil sensitivity and leaves
the origin of the boundary shape unresolved. It does not impose a wall
derivative; the selected cap audit above now evaluates the independent
mass and wall-curve consistency questions.

The [same-grid exterior-operator comparison](prior-art-and-open-claim.md#fixed-grid-exterior-operator-comparison-16-october-2026)
uses the [quadratic derivative](../src/recursive_horizons/nsc_atlas_departure_gradient.py),
[explicit coupled stage](../src/recursive_horizons/nsc_atlas_gradient_evolution.py)
and [callback SSPRK3 driver](../src/recursive_horizons/nsc_atlas_stage_ssprk3.py).
Its [recorder](../scripts/record_nsc_atlas_operator_comparison.py) and
[analyzer](../scripts/analyze_nsc_atlas_operator_comparison.py) reuse the
authenticated 64-cell baseline and unchanged source/gauge cache. Focused
controls are:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_departure_gradient.py tests/test_nsc_atlas_gradient_evolution.py tests/test_nsc_atlas_stage_ssprk3.py tests/test_record_nsc_atlas_operator_comparison.py tests/test_analyze_nsc_atlas_operator_comparison.py -q
```

All 89 controls pass before physical creation, reporting 2,047,112 points
and 3.84375 process CPU seconds across their admitted shared allowances.
The [completed comparison](current-result.md#fixed-grid-exterior-operator-comparison-16-october-2026)
uses 55,398,966 total points/256.25 process CPU seconds. Its one standalone
check passes with 4,040,624 points; the [manifest](../results/development/nsc-atlas-operator-comparison-v1/manifest.json)
binds 104 paths. To authenticate the retained result:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_operator_comparison.py --output results/development/nsc-atlas-operator-comparison-v1 --check
```

A fresh reproduction uses a new output directory without `--check`.
The numerical comparison uses 64 exterior/cap
cells, 512 time steps and both original amplitudes. It shares 450 CPU
seconds/80 million points including creation validation; its standalone
check has 60 CPU seconds/eight million points and no trajectory replay.
Old baseline arrays are referenced rather than duplicated. Same-grid cap
comparisons retain total physical cell integrals without adjacent summation.

The [completed fixed-input spatial comparison](current-result.md#fixed-input-spatial-refinement-16-october-2026)
uses the [recorder](../scripts/record_nsc_atlas_spatial_refinement.py) and
[saved-state analyzer](../scripts/analyze_nsc_atlas_spatial_refinement.py).
It adds zero and finite-amplitude 64-cell windows to the retained 32-cell
comparison, holding the 512 SSPRK3 steps, 1,025 exact cache rows, source
snapshot and computational profile fixed. Run its focused controls or
authenticate the [saved record](../results/development/nsc-atlas-spatial-refinement-v1/spatial.json):

```sh
python scripts/run_python.py -m pytest tests/test_analyze_nsc_atlas_spatial_refinement.py tests/test_record_nsc_atlas_spatial_refinement.py -q
python scripts/run_python.py scripts/record_nsc_atlas_spatial_refinement.py --output results/development/nsc-atlas-spatial-refinement-v1 --check
```

The original creation and validation use 60,087,926 points and 256.8125
process CPU seconds. One independent standalone check passes with
6,684,624 authentication points. The [manifest](../results/development/nsc-atlas-spatial-refinement-v1/manifest.json)
binds 92 files. Creation requires a fresh directory and shares 450 process
CPU seconds and 80,000,000 points across authentication, both trajectories, saved diagnostics
and creation validation. A standalone check has 60 CPU seconds/8,000,000
points; it reconstructs saved physical fields, final stages and diagnostics,
without source integration, cache reconstruction or trajectory replay. The
two new control suites together have 60 CPU seconds/2,000,000 points.
Failures leave an ignored `.build` receipt and no completed canonical
record. Comparisons sum adjacent total physical cap integrals and distinguish
native from restricted derivative stencils, retaining the nonzero source
constraint. Two spatial grids do not establish a convergence order or error
enclosure. The paired wall-scalar spatial difference is 0.774698% of the
fine response peak, while the paired cap `Iq` difference is about 7% of
the old 32-cell paired peak at common normalized cells. Wall curvature and
quiet-profile support remain unresolved. A fresh numerical repeat uses the same producer with
`--output` naming a fresh directory and without `--check`; the completed
record is never overwritten.

The [incoming wall-constraint audit](current-result.md#incoming-wall-curvature-and-transported-constraint-16-october-2026)
uses pinned histories, four one-sided spatial stencils and causal birth
interpolation. Authenticate the completed [record](../results/development/nsc-atlas-wall-constraint-audit-v1.json)
without replaying its source or trajectory:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_wall_constraint.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_wall_constraint.py --output results/development/nsc-atlas-wall-constraint-audit-v1.json --check
```

Creation requires a fresh output path; it never replaces the completed
record. Creation and standalone checking each have a shared allowance of
30 CPU seconds and 2,000,000 points. Creation uses 1,391,719 points and
1.046875 CPU seconds; one independent check passes with the same point
count. Both use 129 analytic vacuum evaluations and no source, profile,
cache, coupled-stage or time-evolution replay. The 19 controls share
30 CPU seconds/200,000 points and pass with 71,773 points and 0.0625 CPU
seconds. All four wall subtractions worsen the disjoint-subset RMS;
spatial/transport accuracy remains unresolved, with no junction-failure
or error-enclosure claim.

The [constraint transport audit](current-result.md#saved-constraint-transport-and-sampling-16-october-2026)
uses only pinned saved arrays, separating the inherited source residual,
diagnostic time-stencil sensitivity and incoming-characteristic origin:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_constraint_transport.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_constraint_transport.py --output results/development/nsc-atlas-constraint-transport-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_constraint_transport.py --output results/development/nsc-atlas-constraint-transport-audit-v1.json --check
```

Its 21 controls pass with aggregate 72,316 points and 0.359375 CPU seconds.
The completed audit uses 689,270 points and 0.921875 CPU seconds. Creation
and a standalone check each have 30 CPU seconds/1,000,000 points; controls
share 30 CPU seconds/300,000 points. Creation requires a fresh path and
checking repeats only saved arithmetic and hash authentication, with no
source/cache/profile/stage or time-evolution replay. The producer, controls
and input bindings remain fixed with the result.
One standalone check passes with 689,270 points.

The current cap representation carries an exact characteristic background and
evolves departures. Run its 72 manufactured, endpoint and independent radial controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_background.py tests/test_nsc_atlas_background_cap.py tests/test_nsc_atlas_background_profile.py tests/test_nsc_atlas_background_evolution.py tests/test_nsc_atlas_background_endpoint.py -q
```

The profile suite has one 30-CPU-second/20,000-point allowance; the cap
suite has 30 CPU seconds and 100,000 points for positive partial quadrature
and an independent narrow-pulse order-4/8/16 comparison. The input adapter
has 30 CPU seconds/20,000 points; the coupled stage has 30 CPU
seconds/100,000 points. The endpoint suite has 30 CPU seconds/20,000 points;
it accepts the exact terminal profile label and rejects an excursion beyond it.
The
[owner](nsc-atlas-finite-amplitude.md#carrying-a-characteristic-background-and-evolving-its-correction)
specifies the background's finite domain, centre limits and free incoming
correction. These tests do not replay the prepared physical source.

The [six-case recorder](../scripts/record_nsc_atlas_background_evolution.py)
keeps the raw cap departures and reconstructed physical integrals separately.
Its controls check profile fitting versus withheld rows, physical reconstruction,
paired comparisons, shared accounting, immutable output and failed-attempt receipts:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_background_evolution.py -q
```

The completed [first physical comparison](current-result.md#characteristic-background-coupled-evolution-16-october-2026) uses 228.9375 CPU seconds and 43,551,872 points. Its shared allowance is 600 CPU seconds and
100 million evaluated points, including authentication. It reuses the bound
source cache and extracts the computational profile once. Its [calculation
record](prior-art-and-open-claim.md#coupled-characteristic-background-comparison-16-october-2026)
fixes six cases and stopping conditions. Creation uses a new output directory;
`--check` authenticates saved data without replaying time evolution:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_background_evolution.py --output results/development/nsc-atlas-background-evolution-v1
python scripts/run_python.py scripts/record_nsc_atlas_background_evolution.py --output results/development/nsc-atlas-background-evolution-v1 --check
```

The [four-panel comparison](figures/nsc-atlas-background-evolution.png) is
drawn only from saved arrays and authenticated input hashes. Its original
render uses 1.15625 CPU seconds and 36,051 arithmetic points. To reproduce
it into a fresh local output (the renderer never overwrites an existing file):

```sh
python scripts/run_python.py scripts/plot_nsc_atlas_background_evolution.py --output .build/nsc-atlas-background-evolution.png
```

The [scalar-time diagnostic](../scripts/diagnose_nsc_atlas_background_time.py)
uses the retained reference profile and four frozen stages. Its 35 controls
exercise independent quadrature and chain identities, hash-only input
authentication and one shared allowance. It performs no time evolution:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_background_time.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_time.py --output results/development/nsc-atlas-background-time-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_time.py --output results/development/nsc-atlas-background-time-audit-v1.json --check
```

Each audit has 30 CPU seconds/150,000 points; controls have 30 CPU
seconds/100,000 points. Creation requires a new output path. Checking
repeats the same bounded diagnostic and excludes elapsed CPU time from
the numerical equality comparison.

The [saved audit](../results/development/nsc-atlas-background-time-audit-v1.json)
completes with 57,152 points and 1.0625 CPU seconds; its standalone check
passes. It exposes substantial cancellation and does not isolate the
accumulated scalar time error. The subsequent [local accepted-step comparison](prior-art-and-open-claim.md#actual-scalar-step-balance-16-october-2026)
uses identical saved starting states and the existing source/gauge cache:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_background_step.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_step.py --output results/development/nsc-atlas-background-step-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_background_step.py --output results/development/nsc-atlas-background-step-audit-v1.json --check
```

Its 22 controls use 5,312 points and 0.125 CPU seconds. The physical
diagnostic uses 138,891 points and 1.546875 CPU seconds, within its
30-second/250,000-point allowance. It takes six local SSPRK2 steps with
twelve coupled stages and saves an exact five-part scalar balance; it
replays no complete trajectory. Creation is exclusive and `--check`
repeats only those same local steps. The [outcome](current-result.md#actual-scalar-step-balance-16-october-2026)
retains the distinction between local increments and full-history error.

The [SSPRK3 operator](../src/recursive_horizons/nsc_atlas_background_ssprk3.py)
uses the existing coupled stage and exact supplied midpoint cache keys.
Its 20 controls take 27,607 points and 0.5625 CPU seconds, plus two
intentional budget-exhaustion charges. The 17 local-recorder controls use
3,808 points and 0.140625 CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_background_ssprk3.py tests/test_diagnose_nsc_atlas_ssprk3_step.py -q
python scripts/run_python.py scripts/diagnose_nsc_atlas_ssprk3_step.py --output results/development/nsc-atlas-ssprk3-step-audit-v1.json
python scripts/run_python.py scripts/diagnose_nsc_atlas_ssprk3_step.py --output results/development/nsc-atlas-ssprk3-step-audit-v1.json --check
```

The [saved local result](../results/development/nsc-atlas-ssprk3-step-audit-v1.json)
uses six stages, 70,118 points and 1.1875 CPU seconds, within its
30-second/150,000-point allowance. It reads the pinned earlier local
SSPRK2 record without replaying those steps. The weighted five-part
ledger, common starting state and stage keys are verified independently.
Creation requires a new path; `--check` repeats only the six new stages.
The [full original-window comparison](prior-art-and-open-claim.md#fixed-input-ssprk3-coupled-window-16-october-2026)
uses four fixed 128/256-step cases and no new source. Its recorder's 27
controls pass with 37,821 points and 1.796875 CPU seconds:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_ssprk3_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_evolution.py --output results/development/nsc-atlas-ssprk3-evolution-v1
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_evolution.py --output results/development/nsc-atlas-ssprk3-evolution-v1 --check
```

Creation has one 300-CPU-second/50-million-point allowance for four
evolutions, diagnostics and authentication together. The standalone
60-second/5-million-point check reconstructs saved fields, four final
stages and diagnostics without a trajectory or source replay. Creation
requires a fresh directory; failures retain their first context under
ignored `.build`. The fixed model and 80 bound paths remain unchanged.

The [completed four-case comparison](current-result.md#fixed-input-ssprk3-coupled-window-16-october-2026)
uses 130.734375 CPU seconds and 24,939,592 points, including validation.
Standalone checking passes with 956,575 points. Its scalar time difference
remains 21.54% of the fine response peak; closer cross-method agreement
does not replace a temporal error estimate. The saved `diagnostics`
contain `finite_minus_zero`, `paired_time`, `paired_versus_ssprk2`,
constraint histories and pre-support controls. The old SSPRK2 paired
histories supply the pointwise comparison `(4*paired_512-paired_256)/3`;
compare this array with the new paired 256-step history without a solve.

The [midpoint-complete comparison](prior-art-and-open-claim.md#midpoint-complete-ssprk3-comparison-16-october-2026)
adds only the missing cache labels and two 512-step cases. The original
source snapshot, profile and spatial operators stay fixed. Its 17 cache
and 19 recorder controls pass with 540,731 points and 2.890625 CPU seconds,
plus an intentional 513-point exhaustion attempt:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_refined_cache.py tests/test_record_nsc_atlas_ssprk3_refinement.py -q
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_refinement.py --output results/development/nsc-atlas-ssprk3-refinement-v1
python scripts/run_python.py scripts/record_nsc_atlas_ssprk3_refinement.py --output results/development/nsc-atlas-ssprk3-refinement-v1 --check
```

Creation shares 300 CPU seconds/50 million points across cache evaluation,
gauge-order comparison, evolution and validation. A standalone check has
60 CPU seconds/5 million points; it reproduces the midpoint cache and two
final radial stages, with no trajectory or source-ODE replay. The 86
input/source bindings remain frozen. Creation requires a new directory.

The [completed 512-step pair](current-result.md#midpoint-complete-ssprk3-comparison-16-october-2026)
uses 174.90625 CPU seconds and 33,604,170 points, including validation.
The standalone check uses 1,418,330 points. Its cache construction uses
575,112 points, including both gauge orders; all original 513 rows are
preserved. The paired scalar temporal difference is `5.84100e-7`, or
1.67616% of the fine peak. `diagnostics` retains full signed comparison
arrays and maximum locations, quiet channels and profile-knot brackets.

The [three-resolution figure](figures/nsc-atlas-ssprk3-refinement.png)
uses only authenticated saved arrays. Its original render uses 21,801
arithmetic points and 1.09375 CPU seconds within a 30-second/200,000-point
allowance. Reproduce it at a fresh path:

```sh
python scripts/run_python.py scripts/plot_nsc_atlas_ssprk3_refinement.py --output .build/nsc-atlas-ssprk3-refinement.png
```

The [compatible transport audit](../scripts/diagnose_nsc_atlas_compatible_transport.py)
uses the fixed [compatible result](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
without replaying its source or time evolution. Its 25 controls exercise
the exact incoming cap pullback, a closed polynomial prefix oracle,
derivative/integration separation, changed-input rejection and the shared
allowance:

```sh
python scripts/run_python.py -m pytest tests/test_diagnose_nsc_atlas_compatible_transport.py -q
```

Authenticate only the pinned input payload, manifest and its 51 source
bindings, without reconstructing a scientific state:

```sh
python scripts/run_python.py -c "import sys; sys.path.insert(0, 'scripts'); import diagnose_nsc_atlas_compatible_transport as audit; print(audit.authenticate()[2])"
```

Reproduce the saved-state audit into a file that does not exist:

```sh
python scripts/run_python.py scripts/diagnose_nsc_atlas_compatible_transport.py --output .build/compatible-transport-audit-review.json
```

Omitting `--output` prints the JSON without saving it. Each invocation has
one cumulative 30-process-CPU-second/100,000-point allowance, including
the before/after input hash checks, two zero-amplitude cap pullbacks and
four frozen finite-amplitude stages. Each stage reconstructs its cap once,
then crosses two departure derivative rules with two prefix quadratures;
it does not integrate a new trajectory. The
[completed audit](../results/development/nsc-atlas-compatible-transport-audit-v1.json)
uses 48,910 points and 0.859375 CPU seconds. Its input identity, audit
source hashes, runtime and counts are saved in the JSON. The
[outcome](current-result.md#compatible-transport-audit-16-october-2026)
distinguishes cap-profile smearing and derivative sensitivity from
quadrature sensitivity, without treating the postprocessed alternatives
as evolved solutions or error bounds.

The compatible successor evolves exact-vacuum geometric departures and
conserved cap integrals. Its single common scalar wall value comes from
the cap gradient primitive. Run its manufactured and analytic controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_vacuum_departure.py tests/test_nsc_atlas_reduced_cap.py tests/test_nsc_atlas_compatible_evolution.py tests/test_nsc_atlas_departure_diagnostics.py -q
```

The 94 controls include actual incoming-wave propagation, moving vacuum
preservation, finite-amplitude wall conversion and independent pure-null
constraint violations. Each regional/control suite shares its declared
30-CPU-second/20,000-point budget. The
[owner](nsc-atlas-finite-amplitude.md#conservative-reduced-variables-and-exact-vacuum-departures)
specifies the integral weights, regular-axis parity and background
subtraction. These controls do not replay the prepared physical source.

The [successor recorder](../scripts/record_nsc_atlas_compatible_evolution.py)
authenticates and reuses the first coupled run's source and reference caches.
Its [completed five-case comparison](../results/development/nsc-atlas-compatible-evolution-v1/compatible.json)
has the same grids, time steps and amplitude, with one shared
180-CPU-second/8,000,000-point allowance including inherited authentication
and final checks. Run the recorder controls and authenticate the canonical
result without a physical replay:

```sh
python scripts/run_python.py -m pytest tests/test_record_nsc_atlas_compatible_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_compatible_evolution.py --check --output results/development/nsc-atlas-compatible-evolution-v1
```

A fresh reproduction requires a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_compatible_evolution.py --output .build/compatible-evolution-review
```

Checking reconstructs raw-to-physical saved fields, final cap constraints,
wall/prefix stages and independent diagnostics. It performs no source,
reference or time evolution. The raw finite-difference and analytic-vacuum
departure diagnostics remain distinct. Fine cap integrals are restricted
by summing adjacent cells. The first failure is saved under the ignored
`.build/compatible-evolution-failures` area with no completed result
directory, automatic retry or change of numerical settings.

The saved creation and validation total 64.40625 process CPU seconds and
5,371,737 charged points. Standalone authentication uses 1,202,991 points;
it reuses the inherited source/reference caches and evolves no trajectory.
All five cases retain 129 slices. Exact vacuum and zero-outgoing recovery,
paired spatial responses, raw finite-case time changes and remaining cap
transport/constraint limits are reported in the
[current result](current-result.md#compatible-coupled-evolution-16-october-2026).

The simultaneous exterior, wall and polar-cap stage uses a source snapshot,
ordered exterior integrals, the constrained regular cap and the existing
hybrid wall law. Its fixed reference wall label evolves the actual wall
lapse through the acceleration equation. Run the analytic and manufactured
controls without a new physical source integration:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_trace.py tests/test_nsc_atlas_exterior_characteristic_slice.py tests/test_nsc_atlas_polar_cap.py tests/test_nsc_atlas_reference_wall_gauge.py tests/test_nsc_atlas_coupled_evolution.py -q
```

The [scientific owner](nsc-atlas-finite-amplitude.md#a-simultaneous-exterior-wall-and-cap-stage)
specifies the coordinate, axis, characteristic and averaging conventions.
The first-order spatial scheme and fixed-step SSPRK2 retain measured
residuals; they do not supply a continuum error certificate.

The [coupled recorder](../scripts/record_nsc_atlas_coupled_evolution.py)
completes five declared cases, with two spatial resolutions and one
separate time refinement. The independent saved-history diagnostics and
original-rectangle overlap have their own manufactured controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_coupled_diagnostics.py tests/test_record_nsc_atlas_coupled_evolution.py -q
python scripts/run_python.py scripts/record_nsc_atlas_coupled_evolution.py --check --output results/development/nsc-atlas-coupled-evolution-v1
```

For a fresh reproduction, use a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_coupled_evolution.py --output .build/coupled-evolution-review
```

Creation uses one cumulative 180-process-CPU-second/8,000,000-evaluation
budget including authentication, with the original bounded source children.
The record separates calculation and verification charges. It retains
129 accepted slices per case, exact-index cached source and reference
inputs, final radial/wall data and independent comparisons. Source hashes
are checked before and after calculation. A first failure is saved under
the ignored `.build/coupled-evolution-failures` area, without creating a
completed result directory or changing the declared grid.

Standalone checking performs no source integration, reference quadrature
or time evolution. It authenticates the original rectangle, reconstructs
the source interpolant and five final radial-cap/prefix states, and
recomputes the saved-array diagnostics. The creation-only
[result](../results/development/nsc-atlas-coupled-evolution-v1/coupled.json)
retains unresolved spatial/interface and cap field-compatibility errors;
see the [current outcome](current-result.md#first-coupled-exterior-wall-and-cap-window-16-october-2026).

The [nonlinear rectangle marcher](../src/recursive_horizons/nsc_atlas_nonlinear_rectangle.py)
assembles the exact source and quiet bottom in the original affine null
rectangle. Its independent nodal diagnostics compute pure null constraints
and geometric mass without using the solved mixed equations as derivatives.
The [four-grid recorder](../scripts/record_nsc_atlas_nonlinear_rectangle.py)
compares zero amplitude and half the source-only bound on 16/32 cells per
direction. Run its implementation controls and authenticate the saved result:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_rectangle.py tests/test_record_nsc_atlas_nonlinear_rectangle.py -q
python scripts/run_python.py scripts/record_nsc_atlas_nonlinear_rectangle.py --check --output results/development/nsc-atlas-nonlinear-rectangle-v1
```

For a fresh reproduction, select a directory that does not exist:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_nonlinear_rectangle.py --output .build/nonlinear-rectangle-review
```

Creation has one 120-process-CPU-second, 60,000-evaluation budget, including
initialization and all four marches. Source children retain their own
30-second/20,000-evaluation limits and debit the parent. Boundary queries
finish before marching. Failure records its phase and any failed cell under
`.build/nonlinear-rectangle-failures`, without writing a completed result
directory. Authentication checks source/output hashes, settings, boundaries,
counts and recomputed diagnostics with zero source integrations or cell
solves. The [owner](nsc-atlas-finite-amplitude.md#original-exterior-rectangle-and-independent-diagnostics)
and [admission](prior-art-and-open-claim.md#nonlinear-exterior-rectangle-16-october-2026)
define the finite scope; sampled differences do not certify a continuum
error bound or the later moving-wall and cap evolution.

The [exact source integrator](../src/recursive_horizons/nsc_atlas_nonlinear_source.py)
retains the original compact preparation and quiet corner. Its amplitude
control comes from a continuous source-only bound. The
[characteristic-cell solver](../src/recursive_horizons/nsc_atlas_characteristic_cell.py)
advances one positive-radius mixed cell with a bounded Newton iteration.
Their independent controls run without a prepared-source regional march:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_source.py tests/test_nsc_atlas_characteristic_cell.py -q -s
```

The source tolerance test prints its bounded comparison receipt. Source
queries retain the shared CPU/evaluation budget and reject extrapolation.
Cell diagnostics report the nonlinear residual, correction and iteration
counts; they do not certify the regional Einstein constraints or continuum
error. Both implementations and their physical domains are specified in
the [same nonlinear owner](nsc-atlas-finite-amplitude.md).

The [normalized bulk operator](../src/recursive_horizons/nsc_atlas_nonlinear_bulk.py)
and [hybrid wall map](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
connect the exact regional equations to a cap whose lapse is determined by
its central clock. Their [owner](nsc-atlas-finite-amplitude.md) derives the
regular axis system and identifies the source, moving-grid and constraint
requirements for a coupled evolution. Run the independent canonical-action,
vacuum, coordinate and mixed-chart controls together:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_bulk.py tests/test_nsc_atlas_nonlinear_boundary.py tests/test_nsc_atlas_hybrid_boundary.py -q
```

These pointwise helpers do not evolve a spacetime. In particular, bulk mass
rates labelled `on_constraints` require the null constraints; first jets
alone cannot measure those residuals. The hybrid map requires a constrained
untrapped cap and a positive incoming radial sum `H`.

The [exact boundary and central-clock helpers](../src/recursive_horizons/nsc_atlas_nonlinear_boundary.py)
evaluate the same-action nonlinear matching algebra without evolving a
spacetime. Their [owning derivation](nsc-atlas-finite-amplitude.md)
specifies compatible gauge inputs and the remaining coupled problem.
Run their independent vacuum, canonical-action and coordinate controls:

```sh
python scripts/run_python.py -m pytest tests/test_nsc_atlas_nonlinear_boundary.py -q
```

These tests create no physical evolution record and select no amplitude
interval for the full model.

The [profile-moment successor](nsc-atlas-profile-moment.md) compares the
pre-arrival mass identity against the original radial evaluator, then
compares local lifts on 16/32 active cells with 40 withheld queries:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_profile_moment.py --output .build/profile-moment-review
python scripts/run_python.py scripts/record_nsc_atlas_profile_moment.py --check --output results/development/nsc-atlas-profile-moment-v1
python scripts/run_python.py scripts/record_nsc_atlas_local_lift.py --output .build/local-lift-review
python scripts/run_python.py scripts/record_nsc_atlas_local_lift.py --check --output results/development/nsc-atlas-local-lift-v1
```

Each complete physical comparison has one shared 120-process-CPU-second,
30-million-kernel budget. The mass pair uses moment orders 4/8 and the
original radial/ray order 16 reference; the local lifts share 81 direct
moment-order-8 boundary jets. Both retain the same 129-node field and
257-knot geometry. Creation refuses an existing output directory. The
local-lift `--check` authenticates its predecessor and reconstructs only
the saved polynomials; it performs no physical solve. Numerical outcomes
are in the [current result](current-result.md#profile-moment-and-local-lift-16-october-2026).

The [conservative evolution implementation](nsc-atlas-scalar-evolution.md)
adds coefficient caching and the forced fixed-coordinate scalar equation.
Its manufactured tests do not themselves supply a physical waveform.

The first forcing-cache comparison shares 65 direct slices between 17/33
training times and 32 withheld times. The latter include radial midpoints
to distinguish time-only from combined time/radial discrepancies:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_forcing_cache.py --output .build/forcing-cache-review
python scripts/run_python.py scripts/record_nsc_atlas_forcing_cache.py --check --output results/development/nsc-atlas-forcing-cache-v1
```

Its 120-process-CPU-second, 30-million-kernel budget includes preparation.
Authentication reconstructs the interpolants from retained slices with no
physical fallback. A cached evaluation inside the scalar solver performs
no new field or metric integration.

The boundary-component diagnosis uses only the authenticated 81 saved
jets and reconstructs eight polynomial fits:

```sh
python scripts/run_python.py scripts/assess_nsc_atlas_lift_components.py --output .build/lift-components-review.json
python scripts/run_python.py scripts/assess_nsc_atlas_lift_components.py --check --output results/development/nsc-atlas-lift-components-v1.json
```

Its 30-process-CPU-second budget excludes physical solves entirely. The
record binds its producer and input, and `--check` repeats the inexpensive
component arithmetic and withheld comparisons.

The clock refinement reuses those outgoing cap jets, recomputes wall clocks
with 16/32 panels, and keeps the same boundary knots and holdouts:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_lift_clock.py --output .build/lift-clock-review
python scripts/run_python.py scripts/record_nsc_atlas_lift_clock.py --check --output results/development/nsc-atlas-lift-clock-v1
```

One field preparation and two wall solves share 120 process-CPU seconds
and 30 million kernels. The check restores saved first jets and clock
polynomials without new physical solves. The 32-panel result meets the
original sampled boundary target; the [current result](current-result.md#clock-resolved-scalar-boundary-16-october-2026)
records the numerical comparison and subsequent propagation status.

The first exterior scalar propagation restores these authenticated inputs
and runs four cases separating cache time, correction grid and ODE tolerance:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_evolution.py --output .build/scalar-evolution-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_evolution.py --check --output results/development/nsc-atlas-scalar-evolution-v1
```

The complete comparison shares 120 process-CPU seconds and 30 million
kernels, with at most 6,000 RHS calls per solve. It performs no new
leading-field, wall, forcing-slice or boundary solve. Authentication restores
the coefficient polynomials and reconstructs 68 discrete wall rates from
retained interior arrays, without replaying propagation. The
[result](current-result.md#exterior-scalar-correction-16-october-2026) records
separate sensitivities and the limit of the formal coefficient.

The central reconstruction reuses the earlier cap clock and current
exterior/wall records, and compares new reflected incoming-ray integrals
at radial orders 16/32 with ray order 16:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_central_scalar.py --output .build/central-scalar-review
python scripts/run_python.py scripts/record_nsc_atlas_central_scalar.py --check --output results/development/nsc-atlas-central-scalar-v1
```

The five common interior arrival cuts and one leading-field preparation
share 120 process-CPU seconds and 30 million kernels. It performs no new
wall, cap-clock or exterior evolution. The check authenticates saved inputs,
quadrature inventory and reconstruction algebra without repeating incoming
integration. The [central result](current-result.md#central-scalar-coefficient-16-october-2026)
retains all terms, the quiet control and separate numerical sensitivities.
Its independent incoming-ray successor retains that saved radial-32/ray-16
case and computes only radial-32/ray-32 at the same five cuts:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_central_ray.py --output .build/central-ray-review
python scripts/run_python.py scripts/record_nsc_atlas_central_ray.py --check --output results/development/nsc-atlas-central-ray-v1
```

One fresh leading preparation and five incoming integrals share a new
120-process-CPU-second, 30-million-kernel budget. Prior costs are recorded
separately. Its check authenticates the reused coarse side and new
reconstruction algebra without physical replay.

The [scalar-input successor](nsc-atlas-scalar-forcing.md) compares a 9/17-knot
Hermite boundary lift against 16 withheld queries, or cumulative exterior
forcing with radial orders 2/4 at nine times:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_inputs.py --comparison lift --output .build/scalar-input-lift-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_inputs.py --comparison forcing --output .build/scalar-input-forcing-review
```

Each complete pair shares a 120-process-CPU-second, 30-million-kernel
budget, including the one field and one wall preparation. Both hold 129
field nodes and 257 geometry knots fixed. Existing wall settings and cap
radial/ray orders 16 remain unchanged. The boundary fits share their 33
unique direct queries; withheld points never become fitted knots in that
comparison. The interpolation target comes from the saved field sensitivity
and is a numerical diagnostic. Exceeding it is recorded as an unresolved
discretization, rather than silently accepting the fit.

Creation refuses an existing output directory. Use `--check --output PATH`
to authenticate a saved `inputs.json` and its manifest without a physical
solve. These inputs do not yet propagate the scalar correction.

The saved pairs are under `results/development/nsc-atlas-scalar-inputs-v1`,
in `lift` and `forcing`. Both use current source bindings. The
[result entry](current-result.md#scalar-lift-and-exterior-forcing-16-october-2026)
reports the unresolved coarse interpolation and the separate exterior
quadrature/transport discrepancies.

The subsequent saved-data assessment changes local Hermite stencil size
while keeping the original training and withheld samples fixed:

```sh
python scripts/run_python.py scripts/assess_nsc_atlas_scalar_lift.py --output .build/scalar-lift-assessment-review
python scripts/run_python.py scripts/assess_nsc_atlas_scalar_lift.py --check --output results/development/nsc-atlas-scalar-lift-assessment-v1
```

This calculation uses no field or wall solve. Its separate creation-only
record binds the predecessor evidence and the assessment producer; `--check`
repeats only the inexpensive interpolation arithmetic. The results motivate
additional samples near the active pulse edges, rather than higher degree
on the unchanged knots.

The [scalar backscatter calculation](nsc-atlas-scalar-metric-response.md)
evaluates the outgoing correction and displaced exterior boundary derivative
from the same prepared field, cap metric and wall solution:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_scalar_boundary.py --comparison field --output .build/scalar-boundary-review
python scripts/run_python.py scripts/record_nsc_atlas_scalar_boundary.py --check --output results/development/nsc-atlas-scalar-boundary-v1/field
```

Its five comparisons are `field` (129/257), `radial` (8/16), `ray` (8/16),
`radial_fine` (16/32) and `ray_fine` (16/32). The finer radial/ray controls
hold 129 field nodes fixed; the other quadrature controls use 257. Geometry
knots remain 257, and cap radial/ray baseline orders are 16/16. Upstream
wall settings remain radial order 2, field-derivative order 32, 8 clock
panels, source order 64 and ODE tolerance $2\times10^{-9}$. The cap field
keeps derivative order 2, centre quadrature 8 and threshold $10^{-4}$.

Each pair shares one 120-process-CPU-second, 30-million-kernel budget,
including field/wall preparation and both independent integral coordinates.
The producer preserves existing output directories and binds source bytes
before and after the calculation. All five saved records use current
sources; `--check` performs zero physical solves. Nine samples, signed
contributions and domain margins document a boundary input, not a boundary
interpolant or propagated waveform. The [current result](current-result.md#scalar-backscatter-boundary-16-october-2026)
states all measured numerical sensitivities.

The [smooth cap calculation](nsc-atlas-cap-response.md) reconstructs the
interior metric and central null-arrival correction from the same prepared
field and wall solution. Its two clock coefficients keep proper wall time
and exterior null label fixed, respectively. Create a comparison, then
authenticate its saved settings, source bindings and bytes:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --comparison field --output .build/cap-review
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --check --output .build/cap-review
```

Each complete comparison shares one 120 process-CPU-second and
30-million-scalar-kernel budget across field preparation, wall evolution,
cap reconstruction, rays and diagnostics. Authentication performs no
physical solve. The baseline uses field/geometry 257/257, cap radial/ray
orders 16/16 and 16 cap time panels; upstream wall radial order 2 and 8
time panels retain the existing field-derivative order 32, source order 64
and ODE tolerance $2\times10^{-9}$.

| Comparison | Changed setting | Fixed field nodes |
| --- | --- | --- |
| `field` | 129 to 257 field nodes | Varied |
| `radial` | Cap radial order 8 to 16 | 257 |
| `ray` | Ray order 8 to 16 | 257 |
| `time` | Cap time panels 8 to 16 | 257 |
| `geometry` | Geometry knots 129 to 257 | 129 |
| `derivative` | Continuum second-derivative order 2 to 4 | 129 |
| `centre_order` | Small-radius quadrature order 8 to 16 | 129 |
| `centre_radius` | Small-radius threshold $10^{-4}$ to $10^{-3}$ | 129 |
| `ray_fine` | Ray order 16 to 32 | 129 |
| `radial_fine` | Cap radial order 16 to 32 | 129 |

The field-129 overrides are explicit in each record's `comparison_base`;
all other settings stay fixed within each pair. The default derivative
order is 2, small-radius order 8 and threshold $10^{-4}$. These comparisons
measure numerical sensitivity of formal metric and null-clock coefficients.
They do not compute the higher-order scalar amplitude or a finite-amplitude
error enclosure. Signed ray integrals and extrema refer to the quadrature
samples, with no claim of continuum extrema.

In the `time` comparison, the prepared field, wall and radial/ray rules are
identical. Only the wall-time coefficient's integration panels change. The
second solution reuses the first solution's nine reference-ray integrals;
`ray_nodes` counts fresh work and separate reused-ray counts preserve that
distinction. The shared budget still includes the first ray calculation
and both wall-time integrations.

The retained `field`, `radial`, `ray`, `centre_radius` and `ray_fine` records
keep the original source context `19983a7b4730d1fdc8e59ae4ba2f3efd560913a6`.
Their bytes are unchanged after the declaration checker gained a
$10^{-13}$ absolute tolerance for geometry-endpoint roundoff. Authenticate
those five records with an explicit context, for example:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_response.py --check --source-context 19983a7b4730d1fdc8e59ae4ba2f3efd560913a6 --output results/development/nsc-atlas-cap-response-v1/field
```

This reads the pinned source bytes from a full Git checkout, then checks
the saved record with the current structural validator. It does not run
the historical producer. Missing Git objects fail explicitly; there is no
substitution of current sources. The `geometry`, `derivative` and
`centre_order` records similarly retain context
`5e58c9fd2866f8adcc107ee22f13f389e7f15c77`; pass that full hash to
`--source-context` for their authentication. The `time` and `radial_fine`
records use current sources. The publication checker selects all three
contexts by explicit comparison name, and authentication reports the
selected context and zero cap solves.

The [moving-wall calculation](nsc-atlas-wall-evolution.md) extends the same
prepared field to its retained wedge and evaluates two mass paths, wall
displacement and the exterior null clock. Create or authenticate a comparison:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_wall_evolution.py --comparison field --output .build/wall-review
python scripts/run_python.py scripts/record_nsc_atlas_wall_evolution.py --check --output .build/wall-review
```

Its independent comparison names are `field`, `radial`, `derivative`, `time`,
`geometry`, `source` and `ode`. Each two-setting comparison includes its field
solves and queries in one 120 process-CPU-second, 30-million-kernel budget.
The saved summaries are formal quadratic coefficients and numerical
refinement diagnostics. The smooth-cap calculation above supplies the next
metric and central-clock connection.

The independent radial/mixed overlap uses the same field on the old rectangle
and a different gravitational integration path:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_wall_overlap.py --output .build/wall-overlap-review
python scripts/run_python.py scripts/record_nsc_atlas_wall_overlap.py --check --output results/development/nsc-atlas-wall-overlap-v1
```

It records nine common cuts, mass and lapse differences, the mixed solver's
own constraint diagnostics, and the cumulative compute budget. Authentication
performs no field, metric or wall evolution.

The [smooth exterior metric](nsc-atlas-bulk-response.md) has separate
finite comparisons for metric grid, leading field, geometry interpolation,
field derivative quadrature, source quadrature and time tolerance. For
example, create a fresh calculation or authenticate an existing one with:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --comparison metric --output .build/bulk-review
python scripts/run_python.py scripts/record_nsc_atlas_bulk_response.py --check --output .build/bulk-review
```

`--check` checks the recorded sources, runtime declaration, configuration and
output integrity; it does not repeat the metric solve. The six comparison
names are `metric`, `field`, `geometry`, `derivative`, `source` and `time`.
Each complete two-solve comparison has its own cumulative 120 process-CPU
second and 30-million-scalar-kernel ceiling. The extended prepared-field
interface leaves the original v8 evidence untouched; replay that earlier
record explicitly with `scripts/replay_history.py atlas-v8`.

The six retained bulk-v1 records were produced at `3470152e1766`, before
the wedge API was added. Authenticate them with their original recorder,
or explicitly repeat all six original numerical studies:

```sh
python scripts/run_python.py scripts/replay_history.py check-bulk-v1
python scripts/run_python.py scripts/replay_history.py bulk-v1
```

The first command checks retained bytes and declarations without a metric
solve. The second requires the recorded runtime and uses each original
comparison's compute budget. Both select the named source context explicitly;
they never replace a failing current-source check with historical code.

```sh
python -m pip install -e '.[paper,dev]'
python scripts/check_lab_snapshot.py
python scripts/check_publication.py
python scripts/check_claim_bindings.py
python -m pytest tests/
python scripts/build_paper.py --check
```

`make verify` runs current checks, every retained test and the deterministic
paper check. CI uses the same pinned stack on Ubuntu and Windows. Historical
replay is an explicit separate operation, selected through `make history-check`
or the commands below.

The history-evaluator test records CPU time as a diagnostic. Set
`NSC_ENFORCE_TEST_CPU_BUDGETS=1` to enforce its original 30-second reference
on a comparable, idle host, running that module serially. Selecting only the
timing test does not run the work it measures. Operational calculation budgets
remain enforced; parallel pytest workers do not provide comparable timings.
Both this control and the geometry replay carry the `cpu_budget` marker.
The geometry replay retains its producer's 20-second cap. For a parallel run,
`-m 'not cpu_budget'` optionally excludes those two whole tests, including the
geometry replay's numerical comparisons; the default serial suite includes them.
The [14 October Win11 follow-up](win11-verification-20261014.md) lists the
11 requested reruns. The [subsequent receipt update](receipt-rebinding-20261014.md)
records the metadata regeneration that resolved its four binding failures.

`check_lab_snapshot.py` authenticates the current paths, hashes and sizes in the
v2 [retained-file receipts](lab-snapshot.json), plus the declared object archive
and unchanged original pins in [source history](../archive/source-history/manifest.json).
The [history compaction record](repository-consolidation.md#explicit-historical-replay)
describes the sparse authentication transport retained after pruning. It
preserves the bytes required by current checks, without keeping full parent
ancestry. Complete old checkouts and campaigns require the Zenodo archive.
`check_publication.py` checks retained result identities, publication routing
and explicit DOI/hash receipts for archived evidence. It does not download or
verify a remote archive. Current checks do not require the retired reassembly
manifest or migration scripts; their carriers are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

Immutable records can bind earlier source bytes. Authentication then uses the
record's explicitly declared Git context, preserving its original numerical
claims and hashes. A pass authenticates that recorded calculation; it does not
identify its original source with the current implementation. Tests of active
APIs use this root package. Original-context authentication and replay are
never a silent fallback from a failed current-source check.

The threshold18 import has two normalized creator path labels. Its metadata
checker authenticates the exact old and current bytes and the single string
replacement, while numerical replay uses current local files. It records the
original and current tree hashes separately; any other change is rejected.

The bounded prepared-history receipt update has a separate exact
[ledger](../results/metadata/prepared-history/receipt-rebind.json). Verify its
current carriers with
`python scripts/run_python.py scripts/rebind_nsc_prepared_history_receipts.py --verify`.
The original execution signatures and numerical members remain fixed.

When changing the manuscript's figures, the retained figure producers are:

```sh
python scripts/build_geometry_overview.py
python scripts/build_characteristic_initialization_figure.py
python scripts/build_source_mass_figure.py
```

Run the paper check after figure changes and inspect the result visually.

The compact analytic successor authenticates the original v7/v8 evidence and
the proof of the fixed-model preparation inverse. Its checker records no new
waveform or numerical inverse:

```sh
python scripts/record_nsc_atlas_observer_inference.py --check
```

The [inference owner](nsc-atlas-observer-inference.md) supplies the mathematical
argument. The checker authenticates that argument's bytes and declaration;
passing the checker is not a machine proof. Creation requires a fresh output
directory, for example `--output .build/observer-inference-review`, and refuses
an existing directory. Its source bindings cover the proof, producer and tests.
Its pinned predecessor manifests and all stored output hashes are checked
without importing historical producers or substituting current source bytes.

The separate [scalar wall check](nsc-atlas-scalar-wall-consistency.md) records
the canonical probe's singular surface stress and regular tangential flux.
It authenticates the fixed preparation and includes exact symbolic controls;
it does not rerun the waveform or evolve the geometry:

```sh
python scripts/record_nsc_atlas_scalar_wall_consistency.py --check
python scripts/run_python.py scripts/record_nsc_atlas_scalar_wall_consistency.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_wall_consistency.py
```

Its owner, producer and tests are source-bound. As with the inverse, creating
another record requires a fresh output directory. Independent analytic review
and exact-algebra checks are distinct from provenance authentication. Bulk
backreaction, normal-force balance and amplitude bounds remain open.

The [characteristic initialization](nsc-atlas-characteristic-initialization.md)
evaluates one borrowed null constraint for the original compact source.
Its corner geometry and affine gauge are fixed; it runs no field or coupled
Einstein evolution:

```sh
python scripts/record_nsc_atlas_characteristic_initialization.py --check
python scripts/run_python.py scripts/record_nsc_atlas_characteristic_initialization.py --check-numerics
python scripts/run_python.py -m pytest tests/test_nsc_atlas_characteristic_initialization.py
```

Authentication checks the source-bound owner, implementation, tests and
prerequisites. Numerical replay separately checks exact identities and
the finite one-dimensional quadratures in the recorded budget. Creation
requires a fresh output directory. Refinement agreement is a numerical
diagnostic, not a rigorous enclosure or a finite-amplitude metric error bound.
That intrinsic constraint alone leaves the transverse geometry open.

The [source-leg mass successor](nsc-atlas-source-mass.md) closes the
transverse radius equation on the same initial null surface. It compares
the geometric mass using the full frame transformation, mixed geometry
and improved-stress transport:

```sh
python scripts/record_nsc_atlas_source_mass.py --check
python scripts/run_python.py scripts/record_nsc_atlas_source_mass.py --check-numerics
python scripts/run_python.py -m pytest tests/test_nsc_atlas_source_mass.py
```

Its numerical replay checks exact identities, independent transverse-field
integrals and separate mass/geometry refinements within the recorded finite
budget. Ordinary authentication performs no numerical integration. The
geometric mass is invariant under coordinate changes, but generally differs
between conformal frames during the input. The negative post-input
coefficient is relative to the fixed reference; it does not select a finite
amplitude or establish a vacuum future. Later wall matching, bulk evolution
and the gravitational correction to the observer signal remain open.

The [normal-force successor](nsc-atlas-scalar-normal-force.md) checks the
remaining regular scalar stress jump on the fixed wall. Its ordinary
check authenticates the new declaration and original scalar-wall evidence;
explicit algebra replay independently tests the projection and force sign:

```sh
python scripts/record_nsc_atlas_scalar_normal_force.py --check
python scripts/run_python.py scripts/record_nsc_atlas_scalar_normal_force.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_scalar_normal_force.py
```

Creation requires a fresh output directory and stays inside the declared
60-process-CPU-second budget. It runs no waveform, pressure history or
coupled metric evolution. Its necessary mean-curvature condition does not
set the wall displacement to zero; individual one-sided pressures and
finite-amplitude validity remain uncomputed.

The [wall-mass closure](nsc-atlas-wall-mass-closure.md) expresses the
exterior mass coefficient using the already fixed leading wall trace.
It checks both null fluxes, full frame conversion and original-frame stress
transport analytically, without a numerical mass curve:

```sh
python scripts/record_nsc_atlas_wall_mass_closure.py --check
python scripts/run_python.py scripts/record_nsc_atlas_wall_mass_closure.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_wall_mass_closure.py
```

Ordinary checking authenticates stored evidence; explicit algebra replay
checks the identities and omission controls. Creation requires a new output
directory within the 60-process-CPU-second budget. A numerical evaluation
needs its own quadrature, field/geometry error control and admission. The
existing waveform samples do not certify this new integral. Coupled cap/wall
geometry and the corrected observer signal were left to the successor below.

The [cap/wall response](nsc-atlas-cap-wall-response.md) derives the
vacuum-subtracted cap mass, wall-displacement functional and compatible
temporal/normal conditions. It records the common model specification with
the two-sided result:

```sh
python scripts/run_python.py scripts/record_nsc_atlas_cap_wall_response.py --check
python scripts/run_python.py scripts/record_nsc_atlas_cap_wall_response.py --check-algebra
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cap_wall_response.py tests/test_record_nsc_atlas_cap_wall_response.py tests/test_nsc_cosmological_atlas.py
```

Ordinary checking authenticates the stored identities, current sources and
unchanged exterior evidence. Explicit algebra replay also repeats the exact
stress, frame and junction calculations. Creating a record uses a fresh
directory and a 60-process-CPU-second limit. These commands integrate no
mass, trajectory or bulk metric. `field_mass_ingredients()` exposes the
first-derivative integrands and endpoints; `wall_response_coefficients()`
exposes the regular displacement equation. The corrected observer signal
requires the next bulk radius/lapse and field-response calculation.

Atlas v8 remains the finite exterior-to-interior waveform result. Its original
replay is explicitly selected below. The fixed input, numerical budget,
refinement diagnostics and independent analytic tail bound remain in the
[response owner](nsc-atlas-prepared-response.md). Refinement diagnostics are not
a rigorous numerical enclosure. The [v8 guide](nsc-cosmological-atlas-v8.md)
describes the outputs and causal limits.

The paper builder renders local figures, authenticates their bytes and compares
two builds in check mode. Visual inspection is still required for a manuscript
change. The software release version and the manuscript revision are independent.

Revision `working-2026-10-16.5` has passed visual inspection of the cover,
contents, abstract, new Section 13.18, quantitative table, conclusion and
evidence appendix. Its 37-page PDF is 2,810,596 bytes with SHA256
`5d28ba600c527daa00e8497203faf9aa001b17d0fc283c90191bed093fb95aaf`.
Two clean builds match the tracked PDF byte for byte. The current paper
and catalog bindings pass, and the source has 32 claims with 112 resolved
evidence links. The saved near-wall audit is documented in its research
owner; it does not change the manuscript's fixed-grid numerical result.

The canonical PDF uses the pinned Zopfli compressor for ReportLab streams, so
standard zlib and zlib-ng installations produce the same encoded document.
This controls compression only; fonts, images, equations and layout retain their
bound inputs. The builder checks its rendering dependencies and keeps strict
byte equality. A failed comparison saves diagnostic PDFs and runtime details
under ignored `.build/paper-check`; the paper workflow preserves those files.

## Preserved editions and source contexts

Routine authentication uses the retained receipts and explicitly pinned Git
objects. Historical publication tests run the original result contracts against
baseline `04a0cb9378977772243ebf7c9be3c74182b030f2`, without restoring deleted
producers into the active package. The original result manifest and numerical
claims remain unchanged. [Preservation and source history](repository-consolidation.md)
describes these receipts; retired migration and result carriers are cited at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

Scientific recomputation is separate. `replay_history.py` selects an original
commit and materializes its source layout under ignored `.build/replay`.
Its unchanged checkers enforce the original evidence hashes and comparison
policies. A replay requiring an archived carrier needs its exact hash-bound
bytes; missing evidence or source objects fail explicitly. Edited current
owners are never substituted for historical source bytes.

### Atlas v1–v6 and supporting editions

The v1–v6 manifests record Python 3.14.8 and Matplotlib 3.11.2. Their original
image checks require byte equality in the recorded Windows environment;
matching package versions alone does not guarantee identical images on another
machine. [requirements-history.txt](../requirements-history.txt) preserves
these scientific pins and supplies the replay tools. Pillow 12.3.0 is a replay
dependency; Pillow and pytest versions were not recorded in the atlas manifests.
For a deliberate original-context replay, use an isolated environment without
installing the current project package there. After history compaction, restore
the exact archived context from Zenodo before these full-context commands;
the sparse authentication store alone does not supply the complete old layout:

```powershell
py -3.14 -m venv .venv/history
.\.venv\history\Scripts\python.exe -m pip install -r requirements-history.txt
.\.venv\history\Scripts\python.exe scripts/replay_history.py prepare
.\.venv\history\Scripts\python.exe scripts/replay_history.py all
```

`prepare` selects the baseline context without scientific recomputation. `all`
recomputes the six atlas bundles within their recorded finite domains and
budgets, authenticates the finite and OPEN companions, and runs the original
publication tests. It does not restart their full research campaigns. The
individual selectors are `atlas`, `finite`, `companion` and `tests`.

With Make installed, the optional historical command is:

```sh
make history-check HISTORY_PYTHON=.venv/history/Scripts/python.exe
```

`HISTORY_PYTHON` defaults to `PYTHON`; it never selects another interpreter after
a failure. `make verify` checks the current project and does not run historical
replay. Current CI runs the retained suite on Linux and Windows with Matplotlib
3.11.2; the former separate historical rendering jobs are retired.

### Atlas v7: original root context

Atlas v7 uses commit `7fee894ad5a4da2ba378f724e45e926401f82670` under
`.build/replay/7fee894ad5a4`:

```sh
python scripts/replay_history.py prepare-v7
python scripts/replay_history.py atlas-v7
```

`prepare-v7` authenticates the original evidence and sources without a field
solve or renderer. `atlas-v7` runs the unchanged original checker, including
exact figure comparisons. That manifest requires Windows, Python 3.14.8,
NumPy 2.5.1, SciPy 1.17.1, SymPy 1.14.0 and Matplotlib 3.11.1. This immutable
historical requirement differs from the current 3.11.2 pin; the current
validation environment is not an exact-image replay environment for v7.

### Atlas v8: original prepared-response context

Atlas v8 uses commit `1154baf290d901e73058392e143d0acb37d22951` under
`.build/replay/1154baf290d9`:

```sh
python scripts/replay_history.py prepare-v8
python scripts/replay_history.py atlas-v8
```

`prepare-v8` authenticates the original source bindings and retained evidence
without running a field solve. `atlas-v8` runs the unchanged original producer
within its finite numerical budget, including exact image comparisons. It
requires the same recorded Windows/Python 3.14.8/Matplotlib 3.11.1 profile as v7.
Use these full replay selectors only with that original environment; the
current renderer does not promise historical PNG identity.

Neither v7 nor v8 is added to legacy `all`. Changed preserved outputs, missing
source objects and incompatible runtimes fail explicitly. Original image
comparisons retain byte equality; numerical comparisons retain their recorded
tolerances. On compatible x86-64 Windows hosts, the established replay setting
is `$env:OPENBLAS_CORETYPE='Haswell'`. CPU dispatch and BLAS configuration matter
in addition to version pins, and a later passing run does not excuse a mismatch.

The retired `make reproduce`, `make reproduce-development` and
`make reproduce-exact` routes stop with guidance. A new campaign needs a named
question, fresh admission and a finite budget. Its producer must run in an
explicitly selected source context; pinned hashes are never changed to fit the
current package.

## Windows runtime

The active Windows environment is the project-owned Conda prefix
`.venv/science-runtime`, using Python 3.14.8 and the pinned scientific stack.
This matches the Python version used by CI. Its interpreter and packages live
inside this checkout; Douglas's Python 3.12 and global Python configuration are
unchanged. Create it from the repository root with:

```powershell
conda create --prefix .venv/science-runtime --override-channels -c conda-forge --no-default-packages --yes python=3.14.8 pip
& ./.venv/science-runtime/python.exe -m pip install -e '.[paper,dev]'
python scripts/run_python.py --runtime-info
```

`python scripts/run_python.py <script-or-python-arguments>` selects this prefix,
clears inherited `PYTHONHOME`, disables the user site and supplies only this
checkout's sources and packages. For example,
`python scripts/run_python.py -m pytest tests/`. Its PATH changes apply only to
the child process. An incomplete new prefix is reported explicitly.

The older `.venv/validation` packages and declared base-interpreter fallback
remain available on checkouts without the new prefix. Both the moved and a
fresh standard venv hit `failed to locate pyvenv.cfg: Access is denied` on this
host; the complete Conda prefix runs successfully. The three already-running
P001 evolutions retain their original interpreter and packages. No live process
is migrated by this launcher change. Normal installations use their activated
interpreter. Root `scripts/lab.py` remains a launcher for the same root package.

The Ryzen 9 6900HX has 8 physical cores and 16 logical processors. Measure
useful throughput before increasing concurrency. The v4 endpoint solver
started three independent processes, each measured near one busy core; its adaptive
time/radial chain is serial and does not currently call Numba. Installed Numba
alone does not parallelize that chain. Compile independent kernels or batches
only after profiling and same-state numerical comparison, retaining ordered
reductions, tolerances and guards. Keep the total worker allocation appropriate
to the host rather than assigning all eight cores to every simultaneous run.

The optional compiled coefficient callback is now available through
`coupled_fixed_stage(..., compiled_coefficients=True)`. The default and live
v4 jobs retain the reference implementation. On the authenticated coarse v3
saved stage after aligning the reference operation grouping, the reference
took 1.9375 CPU seconds and the warm compiled call 0.515625 seconds (about
3.8 times faster). Both used 65 panels, 2647 function
evaluations, 151 accepted radial steps and the unchanged gradient guard.
Mass differed by `2.66e-15`, lapse by `8.67e-19`, and flux by `3.75e-12`.
The strict coefficient-rate comparison (`atol=1e-10, rtol=1e-11`) failed:
maximum differences were `8.04e-9` for cap chi and `8.15e-8` for canonical
density. The reconstructed rate differences in the model's `dz` Galerkin norm
are `6.33e-13` for chi and `2.47e-11` for canonical density, against reference
norms `0.02071195` and `1.47646770`. The scaled Gram condition is `4.70e6`.
The added norm report uses 8.40625 CPU seconds and 9,709,459 charged values
including preparation and all three stage evaluations. This diagnoses the
coordinate dependence of the coefficient comparison; it does not certify a
complete trajectory. The accelerator therefore remains optional; no trajectory equivalence
or live-run speedup is claimed. The batch API uses Numba `prange` for independent
nodes, but the physical stage comparison uses the scalar compiled callback.

Reproduce that performance check with
`python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --compare-compiled-coefficients --output .build/p001-compiled-stage-comparison.json`.

The independent full-mass/lapse diagnostic uses the same authenticated seed,
one reference stage and Gauss orders 32/64 on its existing smooth panels:
`python scripts/run_python.py scripts/diagnose_nsc_atlas_endpoint_radial.py --receipt results/summaries/nsc-atlas-fixed-endpoint-level0-v3.json --independent-mass-lapse --output .build/p001-independent-mass-lapse.json`.
Choose a new output path. It keeps failed comparisons, including a nonzero
temporal mass residual, rather than evolving or refining automatically.
The shared allowance is 60 CPU seconds / 40 million charged values.
Twelve focused diagnostic controls include manufactured mass/rate solutions
with nonzero entering data and a nonlinear radial-coordinate Jacobian.

The observer-clock scalar-average identities are checked with
`python scripts/run_python.py -m pytest tests/test_nsc_atlas_geometric_average.py -q`.
The nine controls include a static periodic lapse and a homogeneous changing
clock. These exact controls do not supply a new parent metric or raw geometry.

For the separate existing conformal common-action realization, reuse its
archived metric snapshots and actual projected rates:

```text
python scripts/run_python.py scripts/measure_nsc_conformal_geometric_average.py --restore-geometry --geometry .build/data/nsc-spherical-conformal-episode-v2.npz --output .build/conformal-geometric-average.json
```

The optional restoration reads only the required member of Zenodo record
23090015's ZIP archive, with a 60 MB transfer cap, and verifies the original
NPZ SHA256 before writing. An existing different cache is rejected. The
53,100,722-byte NPZ stays ignored. The historic `payload_bytes` value counts
JSON plus NPZ. Choose a new output path; no geometry, field or time step is
reintegrated. The 30-CPU-second calculation uses the four existing cases,
three declared times, original observer clock and one initial areal ruler.
Refinement differences remain indicators, and the trace is reported alongside
the vacuum-form defect. This realization is distinct from P001's parent.
The metric and second-rate controls run with
`python scripts/run_python.py -m pytest tests/test_nsc_atlas_geometric_average.py tests/test_nsc_conformal_geometric_average.py -q`.
The [recorded scalar result](../results/summaries/nsc-conformal-geometric-average-v1.json)
binds producing commit `8ecbeff`, the original geometry hash, metric identities
and four-case indicators. Its correction fails vacuum form on the declared
interval. This completes that one application; it does not authorize a new
trajectory or refinement campaign.

The physical-resolution comparison uses that same immutable NPZ. After the
method's tests pass and the source is committed, run once to a fresh path:

```text
python scripts/run_python.py -m pytest tests/test_nsc_conformal_resolution.py tests/test_nsc_conformal_vacuum_generator.py tests/test_nsc_conformal_resolution_correlation.py -q
python scripts/run_python.py scripts/measure_nsc_conformal_resolution.py --geometry .build/data/nsc-spherical-conformal-episode-v2.npz --output .build/conformal-resolution.json
```

If the cache is absent, the geometric-average command above supplies its
bounded recovery. The resolution producer verifies its original hash and
refuses uncommitted method changes. One reference `-Delta_h` eigenspace is
fixed at `T=.05`, then the original and filtered `T=.3` Cauchy states use
the same Hamiltonian and full fermion columns. The independent physical
ruler and cutoff are fixed in the admission. The corrected producer's 200-process-CPU-second
limit includes one eigensolve and two Hamiltonian/vector evaluations; there
is no evolution, grid scan, source fit or constraint reset. Eight BLAS
threads perform useful linear algebra without changing existing processes.
The compact output records the energy-projected vacuum-volume coefficient,
each retained-force residual and both constraint defects. It is not a
cosmological Lambda prediction or a retarded effective action. Full arrays
remain in the ignored cache and can be regenerated. Eigen/replay residuals
do not supply a propagated physical error bound.
The [recorded fixed-state result](../results/summaries/nsc-conformal-resolution-v1.json)
binds method `cd00b35`. It retains one geometric mode and records an energy
coefficient `Lambda_E*L_ref^2=-.100890404750047`, with substantial force
residuals. Its fifteen controls include the archived lifted-rate interface;
the initial shape failure and corrected one-shot application are described
in the admission. The finite probe is closed; reproducing it is not a new
evolution or an instruction to scan cutoffs.

The paired preparation uses the original `nf512_original_T0_*` arrays inside
the same authenticated cache. After freezing its producer, use a fresh output:

```text
python scripts/run_python.py -m pytest tests/test_nsc_conformal_preparation_tangent.py tests/test_nsc_conformal_preparation_fourier_control.py -q
python scripts/run_python.py scripts/measure_nsc_conformal_preparation_tangent.py --geometry .build/data/nsc-spherical-conformal-episode-v2.npz --output .build/conformal-preparation-tangent.json
```

This application holds the action and equal pair occupations fixed. A
specified unitary phase changes the source, and the existing radius Jacobian
determines the tangent to the saved, possibly nonzero constraint level.
The corrected application's budget is 140 process-CPU seconds, one
dense factorization and eight BLAS threads; checkpoints stop further matrix
work when the budget is exhausted. The summary includes the initial radius
and clock-rate derivatives, original constraint defect and unresolved
quadrature remainder. The original defect is preserved; initial constraints
are not solved. It supplies no evolution, light-cone signal or Lambda.
No arrays are written. Ten manufactured controls precede the application.
The [recorded response](../results/summaries/nsc-conformal-preparation-tangent-v2.json)
binds method `76cf805`, costs 9.984375 CPU seconds and gives the initial
fractional observer clock/radius derivative `-.272037713746`. It preserves
the nonzero original defect. The first method's incorrect zero-root premise
and stopped application are retained in the admission and v1 diagnosis.

The coupled derivative has a separate small-grid control command:

```text
python scripts/run_python.py -m pytest tests/test_nsc_conformal_tangent_evolution.py -q
```

It checks the analytic Jacobian-vector action against independent full-rate
differences and Hamiltonian identities. This command reads no archive and
takes no physical time step. The propagated photon response has its own
bounded producer and controls:

```text
python scripts/run_python.py -m pytest tests/test_measure_nsc_conformal_causal_response.py -q
python scripts/run_python.py scripts/measure_nsc_conformal_causal_response.py --geometry .build/data/nsc-spherical-conformal-episode-v2.npz --output .build/conformal-causal-response.json --raw .build/conformal-causal-response.npz
```

Use fresh output paths. The method is frozen before application. One initial
factorization and a directional check precede the original 300/600-step
cases, run in two processes with four BLAS threads each. The setup and case
limits are 140, 1000 and 2000 process-CPU seconds within a 3600-second total
envelope. Both saved baseline checkpoints and both photon events are
required. The JSON limit is 24,000 bytes; the ignored final-state cache limit
is 2,000,000 bytes. A stopped case preserves its last accepted state and
reports missing output. Time-step differences are empirical temporal
indicators, not continuum or physical error bounds. The first application stopped at its radial-force difference check before
evolution. The independent initial Fourier action now supplies that block's
check; its failed central difference remains recorded. Ten controls pass.
The corrected [application](../results/summaries/nsc-conformal-causal-response-v2.json)
uses frozen method `efc4af7` and retains the original stopped result. Both
time steps reach the endpoint and pass their two baseline comparisons.
The redshift derivatives are `-.0009145171828668826` and
`.0027489767089183288`, with empirical time-step differences below `1.8e-13`.
These are not continuum or physical error estimates. Reproduction may
change timing and tiny rounding residues; compare each coefficient to
`1e-9` absolute, with the original input, source hashes, domain, baseline
tolerance and endpoint conditions fixed. No value selected a parameter.
The canonical summary is 20665 bytes. The final-state cache is 477146 bytes,
ignored and hash-bound in that summary. Creation used 1736.1875 CPU seconds
including the stopped first attempt. This response is closed; no unchanged
trajectory replay is needed for the next calculation.

The constrained flat-vacuum exterior has exact controls and one bounded
scalar application:

```text
python scripts/run_python.py -m pytest tests/test_nsc_spherical_vacuum_response.py -q
python scripts/run_python.py scripts/measure_nsc_vacuum_gravitational_response.py --output .build/vacuum-gravity-reference.json
```

Use a fresh output path. No archived state, trajectory, root or grid is
needed. The application reads the original couplings and independently
fixed ruler from the two compact causal/resolution summaries, authenticates
its committed sources, and evaluates one declared Laplace point within
5 CPU seconds and 12000 output bytes. The symbolic controls share 60 CPU
seconds. Zero Maxwell and Dirac fields define this flat reference; it is not
an attachment to the charged periodic episode. The mass sign, kinetic
coefficient and valid retarded half-plane remain part of the result.

The [completed scalar record](../results/summaries/nsc-vacuum-gravitational-response-v1.json)
uses method `6c68d00`; compare its dimensionless coefficients within `1e-10`
absolute, with source/input hashes and the reference definition exact.
Timing and new commit metadata may differ. It records the failed stationary
stability of the flat reference; a negative result is not a producer error.
There is no raw archive to regenerate. The physical validity prescription
and actual nested exterior/state remain required for Lambda.

The smooth timelike boundary response of the owned auxiliary action has
independent exact controls:

```text
python scripts/run_python.py -m pytest tests/test_nsc_spherical_timelike_boundary.py -q
```

Six controls derive the Ricci density from Christoffel symbols and compare
the complete covariant boundary variation with the first-order action.
They also check the physical lapse, both orientations, the auxiliary
surface identity and the Einstein-sector Schwarzschild reference limit.
They read no archive and evaluate no physical boundary energy. The formula
is for a fixed-coordinate, orthogonal timelike face with variations away
from the time caps; a matter state, moving-face corner law and parent
normalization remain additional inputs. The shared symbolic allowance is
60 process-CPU seconds, with no numerical evolution.

The shared-action wall check is a direct symbolic calculation:
`python scripts/run_python.py -m pytest tests/test_nsc_atlas_gravitational_reduction.py -q`.
The seven controls include coordinate Riemann/Ricci, normal Weyl under a
radial boost, pure-trace singular curvature, units and the zero-mass limit.
`schwarzschild_normal_weyl(.4, 1.5)` returns the regular exterior tensor in
inverse-length-squared units. The conditional leading double-layer argument
is in the Schur owner; these tests do not implement a full quadratic junction.

The finite causal light-cone and radar comparison reuses the same cached
metric with its saved rates:

```text
python scripts/run_python.py scripts/measure_nsc_conformal_lightcone.py --restore-geometry --geometry .build/data/nsc-spherical-conformal-episode-v2.npz --output .build/conformal-lightcone.json
python scripts/run_python.py -m pytest tests/test_nsc_atlas_cosmological_distance.py tests/test_nsc_conformal_lightcone.py -q
```

The method uses reception `T=0.3,x=1`, the three predeclared ray lengths and
both radial directions. Radar return trips remain within `T>=0.05`.
Fixed-coordinate normal probes can be accelerated; their redshift is not
a cosmological Hubble inference. The eleven controls include a separate
Jacobi ODE and exact flat/de-Sitter optical and radar limits. New output
paths are required, and full data remain in the ignored input cache.
The [recorded result](../results/summaries/nsc-conformal-lightcone-v1.json)
and [figure](figures/nsc-conformal-lightcone.png) retain the actual predictions:

```text
python scripts/run_python.py scripts/plot_nsc_conformal_lightcone.py --input results/summaries/nsc-conformal-lightcone-v1.json --output docs/figures/nsc-conformal-lightcone.png
```

The plot is a saved-data rendering of the primary case. Its lines connect
three predeclared samples; no field, quadrature or cosmological fit is run.
It evaluates one frozen state, takes no time steps, and retains the existing
60 CPU-second / 40-million-value limit. Do not infer a new grid calculation
or a change in the physical model from this benchmark.

## Flat induced gravitational spectral response

The [profile connection](schur-memory.md#the-induced-spectral-profile-beyond-the-local-gravitational-polynomial)
has independent manufactured controls:

```text
python scripts/run_python.py -m pytest -q tests/test_nsc_spectral_gravitational_response.py
```

After freezing the method, its single admitted application uses:

```text
python scripts/run_python.py scripts/measure_nsc_spectral_gravitational_response.py --output .build/spectral-gravity-profile.json
```

The output must not exist. The producer reuses the saved coefficient moment,
cutoff and observer ruler, compares fixed 16/32 Jacobi rules, and runs no scale
root, full moment integration or trajectory. It bounds its process CPU at
120 seconds and its summary at 12000 bytes. Its Euclidean flat profile
response and empirical quadrature change do not establish a curved retarded
response or a physical Lambda error range.

The [completed record](../results/summaries/nsc-spectral-gravitational-response-v1.json)
uses frozen method `9b68122`. Compare its normalized response and scale
ratios within `1e-10` absolute, keeping the declared profile, input hashes
and reference exact. Runtime, destination and newer commit metadata may
differ. The rounded zero quadrature gap does not certify total accuracy.

## Covariant curvature channels and the static probe

The [covariant response](schur-memory.md#covariant-curvature-channels-and-a-static-gravitational-observable)
has independent curvature and source controls:

```text
python scripts/run_python.py -m pytest -q tests/test_nsc_covariant_spectral_response.py
```

After freezing the method, its one admitted application uses:

```text
python scripts/run_python.py scripts/measure_nsc_covariant_spectral_response.py --output .build/covariant-spectral-response.json
```

The destination must not exist. The producer checks the saved TT profile's
inputs and coefficients, evaluates the two curvature channels with fixed
16/32 Jacobi rules at the old reference scale, and compares its spin-two
factor with that saved result. It reruns no old producer, scale root or
trajectory. Limits are 120 process-CPU seconds and 12000 retained bytes.
The static potential and lensing ratios do not supply a physical length
calibration, a real-time exterior state or a cosmological Lambda.

The [completed summary](../results/summaries/nsc-covariant-spectral-response-v1.json)
uses method `9796add`. Compare dimensionless transfer ratios within `1e-10`
absolute, with the source/input hashes and declared fixed reference exact.
Runtime, destination and newer commit metadata may differ. Its rounded
zero rule-to-rule gaps are not certified zero errors.

## Scope

The main [paper](../paper/nested-space-cosmology.md) is a living atlas manuscript.
The [paper guide](papers.md) routes archived and supporting editions. Spherical
symmetry, specified states and fixed parameters are assumptions of individual
tests. No new constants, propagation or recurrence campaign was admitted during
repository reassembly. New calculations use the [existing admission rule](prior-art-and-open-claim.md).
