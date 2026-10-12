# Local incoming-gate search v1

`nsc_ks_local_gate_search.py` owns the production search state machine for the
switched law

\[
C_\Sigma[g]=U_g C_{\rm up}U_g^\dagger.
\]

It is fixed to 32 Chebyshev coefficients for each of `w` and `U`, 129 nodes
on `I=S(1)+[0.12,0.18]`, 258 residual values and 64 retarded directions.  It
does not issue an existence or non-existence certificate.  Reaching the
`1e-11` nodal target only freezes a candidate for the independent 257-node,
continuous-enclosure and nine-component certification pass.

## Admission

The driver refuses to construct an evaluator unless
`NSC-KS-GATE-ERROR-BUDGET-v5` supplies all nine named componentwise bounds,
encloses their directed sum below `2e-11`, preserves at least `1e-11` for the
physical residual, and records all 60 positive and 120 signed source paths as
closed.  The current repository has only budget v4, with six missing
components, so `--audit-admission` presently reports `OPEN` and no scientific
evolution is launched.

The admitted manifest binds the budget, fixed upstream source, initial
profile, evaluator implementation, solver settings, physical metric and every
declared dependency by SHA-256.  Resume fails if any binding changes.

## Evaluations and acceptance

The concrete driver uses the existing value-only family cache for every trial
history.  Its fixed settings are DOP853, grid 1024, Fourier degree 48,
`rtol=5e-14`, `atol=5e-19`, `max_step=1/8192`, joint control and the exact
192-node source-phase contraction.  A full refresh separately evolves all 60
positive families, retains all 120 signed contributions, adds the three
analytic ell-zero responses as zero changes, and returns the complete
258-by-64 retarded Jacobian.

Each full refresh checks centered finite differences for `w` columns 0, 15
and 31, `U` columns 32, 47 and 63, plus the two remaining columns with largest
Jacobian norms.  Certified value uncertainty is included in each derivative
tolerance.  A mismatch aborts the campaign.

The Levenberg proposal uses the certified geometric principal matrix and the
`H^3(w)+H^1(U)` metric.  Radius positivity and preparation support are checked
before a trial is evolved.  Acceptance requires a fresh value record, positive
model reduction, measured merit reduction larger than the combined certified
evaluation uncertainty, and the fixed support checks.  Broyden updates are
used between full refreshes.  A refresh follows five accepted steps or two
poor model ratios.

## Journal and resume

`scripts/run_nsc_ks_local_gate_search_v1.py --run` writes an append-only
journal under
`results/development/artifacts/nsc-ks-local-gate-search-v1`.  Every trial event
contains the parent, candidate and accepted coefficient vectors; trust radii;
old, predicted and measured merit; reduction ratio; Jacobian and evaluator
identities; Broyden counters; support checks; and the fresh value identity.

Events are written through a flushed temporary file and an exclusive hard
link, so a complete destination appears atomically on Windows without
overwriting an earlier event.  Each event contains a body digest and the full
digest of the preceding event.  Resume authenticates the manifest and the
entire contiguous chain, restores the latest accepted coefficient vector, and
forces a fresh full-Jacobian evaluation before proposing another step.

Run the current admission audit with:

```powershell
$env:PYTHONPATH = (Join-Path (Get-Location) 'src')
C:\Users\swedi\AppData\Local\NSC\validation-venv\Scripts\python.exe `
  scripts\run_nsc_ks_local_gate_search_v1.py --audit-admission
```

An OPEN admission is the expected result until budget v5 exists.  Optimizer
failure remains OPEN and is never converted into NON-EXISTENCE.
