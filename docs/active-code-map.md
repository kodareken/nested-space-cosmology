# Active code map

The active source relation is parent deformation -> inherited compression
offset -> child drawing. `inherited_compression_offset` reads corresponding
local deformation histories and computes additive log-length contrast.
Use `--inheritance-only` with the existing CLI for the three-level controls.
The constant-alpha connection and winding transport remain comparison
constructions; they do not select the enclosing region's history.
`intrinsic_phase_cost` computes the volume-weighted phase cost and directional
response; `generator_least_stretch_cost` and `least_resistance_compression_minimum`
derive and minimize the A0 1D generator least-stretch cost under enclosing
compression; `conservative_neck_shape` and `continue_across_neck_break`
evaluate the constrained neck motion and its throat-inversion continuation into
the child frame; `spiral_decay_to_turn_comparison` checks turning and opposite-lobe
asymmetry without imposing the golden spiral. Use `--intrinsic-only` for all
scale, shape, least-resistance and break-continuation comparisons.
`unified_least_resistance_cycle`, `nested_zero_crossing_map`,
`spectral_band_and_baryon_harmonic`, `branching_cosmos_chain_evolution`,
`cyclic_vortex_diffusion_ring`, `signed_circle_pendulum_duality`,
`vacuum_bubble_crystallization_and_mass_spectrum`,
`nonlocal_to_local_atlas_readout`, `two_force_zero_crossing`,
`spaghettification_to_spherical_winner`, `helical_buckling_under_compression`,
`damped_horizon_oscillation` and `buckling_to_goldilocks_calculation` connect
the single-cost 5-variable cycle, nested `+ -> 0 -> -` local-space readout,
moving Planck-to-Absolute-Zero spectral band, branching `Oooo -> oOOO...` chain
and cyclic vortex-diffusion ring, circle-pendulum `+ -> 0 -> -` duality and
non-local to local Atlas fixed-point and K=7 harmonic readout (`--atlas-only`),
two-force zero-crossing, axial compression buckling (omega^* > 0), conservative
`+`/`0`/`-` damped wave transfer (eta = tanh(pi*gamma/(2*omega))), and
Goldilocks child equilibration (`--buckling-only`).

| Current work | Implementation |
| --- | --- |
| Offset history, distance kernel and winding geometry | `src/recursive_horizons/nsc_atlas_offset_generator.py` |
| Independent geometric and symbolic controls | `tests/test_nsc_atlas_offset_generator.py` |
| Existing finite-coverage calculation | `scripts/check_nsc_atlas_offset_generator.py` |
| Isolated scientific runtime | `scripts/run_python.py` |
| Derivation and explicit assumptions | [Offset construction](nsc-atlas-offset-generator.md) |

Run `python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_offset_generator.py`.
The geometry-only reference calculation is
`python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --winding-only`;
it prints a small result and writes no data files.
Use `--transport-only` instead for the exact conservative evolution and its
open/closed comparison. `winding_transport` evaluates the phase, its first
derivatives, density and flux; `winding_transport_ledger` integrates them
at their exact moving breakpoints. No spatial grid scan is needed.
The old full-coverage calculation imports the conditional parent only when
requested; that parent is not the foundation's input.

The [previous code map](active-code-map-history-20261102.md) preserves
the other implementations. They remain executable at their original paths.
Their old next-step instructions do not restart a campaign.

Default ripgrep discovery follows `.rgignore` and shows the active desk.
For an explicit historical investigation, use `rg --no-ignore 'pattern'`
with the particular source or owner path. Test discovery and imports are
unaffected by these search exclusions.
