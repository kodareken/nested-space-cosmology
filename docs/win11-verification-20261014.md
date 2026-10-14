# Win11 follow-up verification, 14 October 2026

Baseline: `f2671a8c4c72add006ba746a555a601a311b9629`, local `main`.
Runtime: Python 3.14.8 with the root validation environment, Numba 0.68.0,
NumPy 2.5.1 and one BLAS/OpenMP thread. The requested 11 nodes ran serially,
without xdist: **7 passed, 4 failed in 184.01 seconds**, with no skips.

| Requested node | Win11 result |
| --- | --- |
| `tests/test_nsc_geometry_graded_window.py::test_record_replay_is_deterministic` | PASS |
| `tests/test_nsc_compatible_prepared_history.py::test_integrated_receipt_and_provenance` | PASS |
| `tests/test_nsc_middle_source_coverage.py::test_two_original_rows_then_replay_resume_and_refusals` | PASS |
| `tests/test_nsc_ks_history_evaluator.py::test_zzz_cpu_budget_is_control_only` | PASS |
| `tests/test_nsc_middle_source_coverage_v2.py::test_two_original_rows_then_replay_resume_and_refusals` | PASS |
| `tests/test_nsc_retarded_radial_response.py::test_radial_replay_rejects_errors_beyond_operand_envelopes` | FAIL |
| `tests/test_nsc_retarded_radial_response.py::test_saved_replay_has_no_propagation_and_preserves_full_response` | FAIL |
| `tests/test_nsc_retarded_spatial_diagnosis.py::test_saved_and_local_diagnosis_replays_without_any_evolution` | FAIL |
| `tests/test_nsc_rn_stencil_audit.py::RecordedVersusExactTests::test_extended_precision_is_unavailable_here` | PASS |
| `tests/test_nsc_ks_difference_receipt.py::test_exact_local_arrays_and_separate_phase_indicators` | PASS |
| `tests/test_nsc_upstream_compensator.py::test_saved_replay_retains_all_source_fibers_and_fixed_coefficients` | FAIL |

## Shared failure boundary

All four failures raise `response verification producer/input changed` in
`scripts/derive_nsc_retarded_compatible_response.py`. Its saved signature has
17 bindings. Exactly one differs:

- Path: `results/development/nsc-compatible-prepared-history.json`.
- Saved digest: `66cb7f086bb60bc646533ed09ca565cc02fedcafbc609d68b5dcf662a2e61318`.
- Current digest: `6eae36dbb91885caee51ef6d7e1ec0fa06983a70765f814761128b64c2d6514b`.

The only changed JSON leaf is the requested owner-document SHA correction,
to `870bac904d641a4ad1cfb96643faf1640f08bc4a82922ddb3561d8635b31db69`.
The document's earlier edit replaced an obsolete citation with the Zenodo DOI.
The downstream saved receipts still declare the original input identity.
Their numerical comparisons are not reached. This identifies the current
Windows blocker; it does not establish the causes of the earlier Linux run.
The remaining source bindings, saved arrays and numerical tolerances are unchanged.

## Timing controls

The history evaluator's 30-second reference is opt-in through
`NSC_ENFORCE_TEST_CPU_BUDGETS=1`; its default diagnostic remains finite and
nonnegative. Selecting just its timing node measures no preceding module work.
The additional full-module run exercises that work and records 3.734375 CPU
seconds at the diagnostic, below its 30-second reference.

The geometry deterministic replay and history timing control carry the
`cpu_budget` marker. Geometry retains its producer's 20-second cap and all
numerical comparisons; its marked rerun passed. The two middle-source tests
contain no 20-second test assertion: both passed serially, taking approximately
79 seconds each for their calls. Their operational budgets remain unchanged.
Optional `-m 'not cpu_budget'` excludes those two whole marked tests; default
serial verification includes them.

The Win11 POSIX adapter reports `resource`/`fcntl` stand-ins: OS CPU timers
and advisory run locks are unavailable. Explicit Python timing checks still run.

## Other changes and checks

- Numba 0.68.0 is pinned in both dependency files; the history requirements
  comment now names Matplotlib 3.11.2.
- Private Mac path labels are normalized in 24 results JSONs and the storage
  necessity audit. Numeric spellings and digest values in those 25 files are
  unchanged. Historical scratch/Trash labels are identified as historical.
- The two threshold18 creator-label changes have exact byte-transformation
  authentication and separately declared original/current tree hashes. The
  current local numerical inputs remain required. Original v6/v7 records and
  source-history pins remain unchanged.
- The obsolete restored-record private-path exemption is removed. Current
  normalized records receive the regular publication checks.
- 18 threshold18 metadata/majorant checks passed. A further 60 focused tests
  passed, including the full history evaluator, marked geometry replay,
  publication boundaries, catalog and historical replay controls.
- The five RN modules previously blocked by missing Numba collected all
  67 tests successfully. Collection alone does not execute those tests.
- Claim routing passed: 22 claims and 41 resolved record links. Retained-file
  verification passed, including the 25 unchanged pinned source objects.

An additional `build_finite_regeneration.py --check` could not complete because
the companion source archive `paper/finite-regeneration-arxiv.tar.gz` is absent
locally following the earlier cleanup. Its original snapshot authentication
passed in the catalog tests; full companion verification remains incomplete.
The existing `make finite-check` target also still names a retired verifier;
neither issue was changed in this dependency and test follow-up.

No full retained-suite run is claimed for this follow-up. Raw logs and JUnit
reports are local under `.build/ci-fix-20261014/`; the requested node list above
can be passed directly to `python scripts/run_python.py -m pytest -q`.
