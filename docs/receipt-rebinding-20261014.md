# Prepared-history receipt metadata, 14 October 2026

The four response replays requested after commit `88b2e7b1` now pass serially
on Win11: **4 passed in 20.48 seconds**. Their saved numerical values and
NPZ numerical members are unchanged.

The prepared-history control's corrected owner-document hash gives its current
JSON digest `6eae36dbb91885caee51ef6d7e1ec0fa06983a70765f814761128b64c2d6514b`.
Live downstream verification bindings formerly named
`66cb7f086bb60bc646533ed09ca565cc02fedcafbc609d68b5dcf662a2e61318`.
The update propagates the current digest through the live response chain.
It performs no field evolution, ODE solve or new scientific calculation.

## Scope and preservation

The update covers 13 JSON receipts and five NPZ metadata carriers:

| Receipt | Carrier updated |
| --- | --- |
| Retarded compatible response | JSON and NPZ metadata |
| Retarded response resolution | JSON and NPZ metadata |
| Retarded Fourier response | JSON |
| Retarded radial response | JSON and NPZ metadata |
| Upstream compensation rank | JSON |
| Upstream compensator | JSON and NPZ metadata |
| Retarded spatial diagnosis | JSON |
| Evolved incoming state | JSON and NPZ metadata |
| Evolved incoming constraints | JSON |
| Incoming endpoint obstruction | JSON |
| Incoming fixed transfer | JSON |
| Retarded fixed-C0 certificate | JSON |
| Incoming jet rigidity | JSON |

All 62 replacements are keyed digest strings. Every file keeps its byte count.
All numerical JSON token spellings, controls, domains, tolerances, statuses,
recorded runtimes and scalar measurements remain unchanged. All 193
nonmetadata NPY members remain byte-identical, including signed-zero bits and
three captured source members. Original execution signatures retain their
capture identities; current verification signatures name the rebound inputs.

The [rebinding ledger](../results/metadata/prepared-history/receipt-rebind.json)
enumerates exact old/current identities and permitted pointers. The maintenance
tool authenticates the original carriers at
`88b2e7b1bf265772eac7108827ccb5539889e5bd`, rebuilds the bounded transformation
and compares it with the current local files. Editing a number and refreshing
the ledger's current hash cannot approve that edit.

Historical campaigns, KS capture files, majorants, source-history objects and
frozen paper evidence retain their original identities. A strict KS consumer
may authenticate a listed current metadata successor against its original
capture through `nsc_receipt_metadata`; this API reports the two identities
separately. Numerical computations continue to use current local inputs.
Original Git bytes serve only the transformation authentication.

The prepared control itself declares current source and document bindings.
Its source authentication therefore uses current bytes even inside the
explicit original response-source context `329ef809e466a86ed4a290ee32c19485f85add9d`.
Other original producer-source authentication retains that explicit context.

## Verification commands

```sh
python scripts/run_python.py scripts/rebind_nsc_prepared_history_receipts.py --verify
python scripts/run_python.py -m pytest -q \
  tests/test_nsc_retarded_radial_response.py::test_radial_replay_rejects_errors_beyond_operand_envelopes \
  tests/test_nsc_retarded_radial_response.py::test_saved_replay_has_no_propagation_and_preserves_full_response \
  tests/test_nsc_retarded_spatial_diagnosis.py::test_saved_and_local_diagnosis_replays_without_any_evolution \
  tests/test_nsc_upstream_compensator.py::test_saved_replay_retains_all_source_fibers_and_fixed_coefficients
```

`--apply` accepts only the original, SHA-bound eighteen-file input state and
plans all replacements before writing. It is a bounded maintenance operation;
`--verify` authenticates the resulting current state without rewriting it.

The Windows PDF writer now stages beside its output so atomic replacement
works when the checkout and temporary output are on different drives. The
portability tests measure actual long-double capability, authenticate original
and current fingerprints separately, and validate a fresh proof's exported
cell census against its directed proof returns. The original witness hashes,
original 820/23 census and all numerical upper bounds remain asserted.

Local logs and JUnit reports are retained under `.build/receipt-rebind-20261014/`.
GitHub CI status is reported for the final pushed commit separately.

## Finite replay arithmetic on Linux, 15 October

The first full CI run passed on Windows (2,834 passed, 12 skipped) and exposed
two remaining Linux replay comparisons. Saved receipts and arrays did not
change while correcting these comparisons.

The `401_64` stationary-field drift maximum differed by one binary64 ulp:
`0.000431005732696872` versus the stored `0.00043100573269687207`.
That leaf had no allowance. The existing stationary-trace maximum rule now
enumerates field and trace drift maxima in all three saved cases. Its allowance
charges both evaluations of the uncancelled subtraction operands and complex
modulus, with the same operation counts. Types, peak identities, source hashes,
statuses and every unlisted field remain exact.

The compensator's six raw fields remain byte-exact in all five families.
Stationary maps differ across hosts by at most `1.3877787807814457e-16`, inside
the existing archived-map diagnostic allowance. Both saved and current maps
are checked independently against the same authenticated history arrays.
Each derived covariance is checked independently against its own exact
rational `F C F*` product, using the uncancelled operand masses and two
three-term complex contractions. Zero source mass has zero allowance;
gradual underflow has an explicit operation-count term. The complete
compensator JSON comparison remains exact.

The `3e-11` numerical control, calculation CPU caps, original execution
identities and OPEN physical qualifications remain unchanged. Negative tests
reject shared corruption, changes beyond the operand allowance, raw-input
byte or type changes, nonfinite values and altered metadata.

Linux diagnostics use this same checkout and pinned libraries in an ignored
project-local validation environment. Its Ubuntu WSL interpreter is Python
3.14.4; GitHub CI uses Python 3.14.8. Both Linux BLAS libraries use one Haswell
thread. The four requested tests passed again serially on Win11 in 26.98
seconds and on Linux in 75.25 seconds. The additional comparison and rejection
tests passed on both platforms (30 tests each, without repeating the two full
replay nodes). The independent eighteen-carrier preservation check and the
claim-binding check passed again after these code changes.
