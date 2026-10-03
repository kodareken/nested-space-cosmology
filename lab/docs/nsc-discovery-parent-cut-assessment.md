# Completed exterior-population cut assessment

The [consumer](../scripts/assess_nsc_discovery_parent_cut_response.py) reads
`prepare`, the locked `prediction`, and both completed
`measurement-a±0.05-T2.25` records in
[the saved batch](../results/development/nsc-discovery-parent-cut-response-v1/).
Its domain is the NF128 strong balanced minus parent at fixed
$k=0.45388971484879426$, protected cuts $x=3,5$, child window $[3.5,4.5]$
and matched centre clock $\tau=2.8568786296681035$.
The [producer method](nsc-discovery-parent-cut-response.md) owns preparation
and evolution. This consumer performs neither.

The stages bind producer commit
`7004a7915a9375c2d5a420853158c4cf19bcdca1` and 39 source hashes through Git
provenance. The consumer authenticates payloads, the baseline's array hashes,
initial state/direction digests and forecast/measurement stage links. It
also records dtype/shape array hashes, actual current evaluator hashes and
unchanged input hashes. Historical source authentication does not require
a later note to match historical note bytes. Existing stages keep their
permissions and bytes; the assessment writer creates one new readonly file
outside the stage directory, refuses overwrite, and enforces 64 MiB.

```sh
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_parent_cut_response.py
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_parent_cut_response.py --output results/development/nsc-discovery-parent-cut-assessment-v1.json
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_parent_cut_assessment.py -q
```

At this clock, the $+5\%$ exterior occupation changes total child content by
$0.00284200211428$, versus locked $0.00285160805944$ (residual $0.3369\%$
of the predicted effect), and proper length by $0.00476723947267$, versus
$0.00477253592531$ ($0.1110\%$). The negative arm gives
$-0.00286118972482$ and $-0.00477780635435$. Both clocks match within
$1.78\times10^{-14}$; total energy is approximately $1.38\times10^{-11}$.
The full record recomputes all nine headline readouts, CAR eigenvalues and
both definitions of effect-relative residual.

The original child-column weight $c_0=0.3062058613$ stays fixed; only $c_1$
changes by $\pm5\%$. Initial $Q,r,\Phi$ are fixed by the authenticated
preparation owner, with momenta re-solved under its finite constraints.
Frame, reference and source columns remain unchanged. The assessment
recomputes each original column's matched-clock derivative as
$\int_C[\delta c_j|\Phi_j|^2+2c_j\Re(\Phi_j^\dagger\delta\Phi_j)]/\Delta x_q
-\dot P_j\delta\tau/N$, using the actual projected leading velocity.

| Original column | Baseline child content | Locked derivative | Held +5% change | Held −5% change |
|---|---:|---:|---:|---:|
| child, 0 | 0.0385173686504 | 0.00861542060057 | 0.000429207849261 | −0.000432368908844 |
| exterior, 1 | 0.0509827134758 | 0.0484167405883 | 0.00241279426502 | −0.00242882081598 |

The derivatives sum to $0.0570321611888$. Centered held slopes differ from
them by $0.0040274\%$ and $0.0012181\%$; even-amplitude remainders are
$-1.58052979116\times10^{-6}$ and $-8.01327548076\times10^{-6}$.
Column norms remain approximately one. The nonzero column-0 response is a
change of its regional localization, despite its unchanged occupation;
it is not simply the addition of exterior occupation inside the cut.
For an unchanged frozen generator and fixed $c_0$, that column's derivative
would be zero. This is an algebraic control statement; no frozen rerun is
performed. Column tags identify computational ancestry, not particles or
energy parcels.

Recorded stage CPU times sum to $126.950233$ seconds, below the 300-second
allowance. The final saved aggregate is $126.946791$ seconds; the 3.442 ms
difference comes from the separately sampled bookkeeping fields. Both are
reported. The baseline child content is still depleted from its initial
value. The source-to-incoming-to-child relation is a conditional finite-band
numerical result, with no continuous error certificate, strong-curvature
EFT guarantee or autonomous-renewal conclusion. Spatial confirmation has
no completed result in this assessment.

Does the original-child-column feedback survive spatial confirmation at the
same proper clock?
