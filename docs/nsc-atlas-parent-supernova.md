# Inherited parent against supernova distance shapes

## Decision and preparation, 28 October 2026

The quiet cap's 7.278-sigma contrast and the controlled inherited matter
background answer different questions. The next calculation compares their
actual distance shapes, together with a LambdaCDM control, using the same
pinned Pantheon+ release and covariance. The latest request explicitly
authorizes that control; the earlier instruction excluding such comparisons
does not describe this task.

Three routes were considered: another local wall repair, adding an adjustable
dark component to the cap, and testing the existing inherited parent without
altering its source. The third connects the new source theorem to real data
and reveals whether the missing component is evolution, preparation, or
observer coverage. It introduces one fitted observer epoch, not a new dust
normalization. A magnitude intercept belongs to the data likelihood and is
profiled analytically. Controls introduce their own explicitly labelled
matter fraction. All state ownership stays local to the calculation.

Use `nsc_atlas_inherited_matter.exact_initial_geometry` and its uniform bound;
reuse `nsc_atlas_cosmological_distance.distance_contrast` and the original
publisher cache/verification. Extend the package with a small FLRW forward
map and common covariance likelihood. No new cosmology dependency or field
simulation is needed. Source parameters remain rho_v=9/32, E0=3/8 and the
exact junction curvature; the high-frequency member is omega=10^9.

1. Compare quiet de Sitter, fixed flat LambdaCDM fractions .3, .315 and .334,
   the inherited parent at A_o=2, and separately calibrated one-parameter
   LambdaCDM/parent histories on the existing 596-row selection.
2. Repeat the likelihood on non-calibrator rows with zHD>=.025 through the
   release maximum. Every inherited-parent emitter must have A_e>=1, so
   the fitted observer must satisfy A_o>=1+max(zHD). Do not extrapolate the
   protected future through its initial Cauchy slice.
3. Report both the unchanged one-dimensional contrast and the full selected
   covariance likelihood with a common magnitude intercept. Count fitted
   shape parameters and distinguish calibration from prediction.
4. Propagate the proved high-frequency density envelope to the light cone;
   check quadrature independently. A domain failure or a poor likelihood
   redirects the source/history construction instead of adding another
   component or silently modifying the data.

The calculation has a 120 process-CPU-second allowance, four BLAS threads,
at most 200 objective evaluations per fit, and zero field/grid evolutions.
Only compact summaries, hashes, producer and regeneration commands are saved.

Completion requires an identical data/zero-point/covariance protocol for both
models, a protected-domain check and a stated numerical/model error. A later
independent distance/growth comparison with the calibrated state held fixed
can falsify its physical interpretation; these already-seen data are a
retrospective comparison. Main risks are an inconsistent heliocentric
prefactor, extrapolation before A=1, and treating a calibrated epoch as a
parameter-free prediction. Each is an explicit code/test boundary.

## Remaining connected closures

The absolute SI ruler, inherited material-wall state, unconverged P001
endpoint, baryon sector and baryon asymmetry retain their separate owners.
A successful supernova shape cannot close them. Grok audits the retained
wall and endpoint equations without running new grids; their useful
consequences will be integrated with this observational result.
