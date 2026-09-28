# Contributing

Recursive Horizons welcomes criticism, calculations, counterexamples, and clearer formulations. Contributions are evaluated by whether they make the project easier to falsify or reproduce—not by whether they protect its current narrative.

## Before opening a change

1. Read [docs/epistemic-status.md](docs/epistemic-status.md) and [docs/claim-ledger.md](docs/claim-ledger.md).
2. State whether the change concerns established physics, a repository derivation, a model postulate, an analogy, an observation, or an open problem.
3. Cite primary literature for technical factual claims.
4. For a numerical result, include inputs, units, uncertainty treatment, code path, package versions, and an explicit failure condition.
5. Do not multiply anomaly p-values or call correlated/post-selected evidence independent without a justified joint likelihood.

## Verification

Run:

```bash
python3 -m unittest discover -s tests -v
python3 scripts/reproduce_core.py --output results/core-identities.json
python3 scripts/reproduce_controlled_model.py --output results/controlled-model.json
python3 scripts/reproduce_fgc_action.py --output results/fgc-1-action-gate.json
python3 scripts/reproduce_fgc_metric_variation.py --output results/fgc-1-metric-variation.json
python3 scripts/reproduce_fgc_hyp1_reduction.py --output results/fgc-1-hyp1-reduction.json
python3 scripts/reproduce_fgc_hyp1_symbol.py --output results/fgc-1-hyp1-symbol.json
python3 scripts/reproduce_fgc_hyp1_modified_harmonic.py --output results/fgc-1-hyp1-modified-harmonic.json
python3 scripts/verify_data.py
python3 scripts/build_paper.py
python3 scripts/check_repo.py
```

For a fast check of the current extension surface, run:

```bash
make verify-current-development
```

Before closing a scientific change, run the complete sealed pre-result gate:

```bash
make verify-fgc-sf1-foundation
```

`make verify-fgc` is retained as an alias for that foundation gate. The old
result-rewriting recipe is retired, and completed one-shot campaign targets
now fail before Python. Historical source files remain available for compact
reproduction and provenance; they are not active execution invitations.

`FGC-1-ACT1` is an action and linear-background certificate;
`FGC-1-VAR1` is only its covariant bulk metric-variation certificate; and
`FGC-1-HYP1-RED1` is only an exact local spherical-reduction and uneliminated
second-jet certificate; `FGC-1-HYP1-SYM1` is only an exact pointwise covariant
scalar-quotient, ADM principal preflight, and rejected-formulation certificate.
`FGC-1-HYP1-MHG1` adds an exact pointwise modified-harmonic gauge/metric/scalar
principal system, gauge-constraint propagation cone, and complete frozen
radial eigenbases, but not complete lower-order sources or a uniform open
domain. A contribution must not turn the still-open lower-order formulation,
reduction-constraint propagation, open-domain full-system hyperbolicity,
collapse, or affine-Raychaudhuri flags true without supplying the corresponding
derivation and reproducible evidence.

A calculation may be valuable even when it falsifies a proposal. Negative results belong in the claim ledger and changelog, not in a private discard pile.

## Scientific writing rules

- Say “conditional on” when a result depends on a postulate.
- Say “inferred from” when the target observable is used to determine an upstream parameter.
- Reserve “prediction” for a quantity fixed without using the data later compared with it.
- Reserve “measurement” for a validated estimator with uncertainty and provenance.
- Separate a coordinate statement from an invariant or locally observed one.
- State the domain of a theorem or model result; symmetry-reduced bounces are not automatically generic collapse results.

## Archive policy

Historical files are immutable evidence of the project’s development. Correct canonical documents and add an archive note; do not silently rewrite an old snapshot to make it look prescient.

Default `rg`, pytest, and Ruff discovery exclude `archive/`. Use an explicit
path (and `rg --no-ignore` where needed) only when historical evidence is the
actual subject of the task.
