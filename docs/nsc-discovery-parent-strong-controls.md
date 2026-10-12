# Matched strong parent controls

The [consumer](../scripts/assess_nsc_discovery_parent_strong_controls.py)
reads the completed strong balanced $T=1$ episode, its exact $T=1\to3$
successor, the new two-sign empty-prepared-source control at the same
$k=0.45388971484879426$, and the saved frozen parent-heavy trajectory.
It performs no preparation, evolution or source solve. The numerical
producer frozen for the new empty control is
`5ce6bd02e84c37764bb3df57d405b7c3839cbfd3`; historical strong and frozen
manifests authenticate their own producing commits through Git/source
hashes. Every snapshot NPZ and dtype/shape array hash is checked. Snapshot
JSON/NPZ files must be immutable; manifests, scalar streams, input files
and actual current evaluator source hashes must remain unchanged during
assessment. The strong successor's stale observation headers are never
used as measurements.

```sh
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_parent_strong_controls.py
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_parent_strong_controls.py --output results/development/nsc-discovery-parent-strong-controls-v1.json
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_parent_strong_controls.py -q
```

Stdout is the default. `--output` creates one new readonly record, bounded
by 64 MiB, outside the sealed input directories. It refuses an existing
path and does not mutate sealed evidence. The integration owner creates
the final record after freezing evaluator bytes.

## Attribution and frozen reuse

The empty preparation retains the strong initial $Q,r$, frame, observer
and occupation weights, selects the same authoritative `parent_k`, zeros
the spinor/source columns, and re-solves C/D. Both canonical momenta change.
It is a source-free preparation comparison, not force deletion on the
occupied Cauchy state. Empty matter probability is zero; its fractional
localization is undefined and is returned as null. This does not imply a
filled negative vacuum.

The frozen trajectory has exactly the strong initial $Q,r,\Phi,W$ and
original source/observer columns. The consumer rejects any mismatch. It
checks equal Dirac coefficients and Fourier momentum symbols, equal
initial generator actions on both columns, unchanged frozen canonical
geometry arrays, and zero saved geometric jets. In the conformal chart
the Dirac kinetic coefficient is $L/Q=1$ and the mass coefficient is
$\kappa Q$, so occupations and canonical momenta do not change this frozen
generator. Reusing its saved columns therefore requires no new trajectory.
The original columns are reweighted by
$(303.99138553987126,94.76215444295961)$ to the strong balanced occupations.
Direct NPZ child integrals agree with the scalar reweight at all five
stored states (recorded arithmetic gap zero).

| Frozen coordinate time | Strong-weight child fraction |
|---|---:|
| 1 | 0.3729217437 |
| 1.25 | 0.3762049855 |
| 1.5 (scalar only) | 0.3903296378 |
| 2.25 | 0.3942388995 |
| 3 | 0.3634687913 |

Only weighted regional probabilities are reused. Old aggregate work,
field/gravity energies and constraint readouts do not become strong
readouts after reweighting. Occupied and empty snapshot energies,
constraints, child length, centre radius and projected curvature are
recomputed from their own arrays.

## Centre-clock comparison

Saved NPZ normal clocks are authoritative at stations. Between stations,
the consumer integrates positive scalar-stream centre clock rates by
trapezoids and scales each interval to its exact saved clock increment.
The maximum unscaled interval discrepancy is $0.001066$ for occupied
minus, $0.000200$ for occupied plus, $0.000723$ for empty minus and
$0.000311$ for empty plus. These are sampling indicators, not bounds.
Interpolation is linear in proper clock and strictly refuses extrapolation.
Empty length/radius comparisons use the explicitly labeled station
interpolation; empty coordinate time and energy use the scalar stream.

| Occupied sign / coordinate time | Occupied centre $\tau$ | Occupied child fraction | Frozen fraction at same $\tau$ | Difference |
|---|---:|---:|---:|---:|
| minus / 1 | 1.6689173396 | 0.3403922161 | 0.3887899351 | -0.0483977190 |
| plus / 1 | 0.5587679859 | 0.3261680959 | 0.3729154737 | -0.0467473777 |
| plus / 3 | 0.7064165019 | 0.0596698351 | 0.3729153649 | -0.3132455298 |
| minus / 3 | 3.0507894717 | 0.0904964242 | outside saved range | undefined |

At minus $T=1$, the matched frozen coordinate time is $1.6691805542$.
Coarsening its scalar fraction sampling to coordinate spacing $0.1$
changes the matched fraction by $-0.0003794544$; the maximum such difference
among these occupied stations is $0.0006956021$. Neither number supplies
an interpolation error bound.

At coordinate $T=3$, occupied/empty child lengths are
$0.1493327575/0.1268874675$ for minus and
$0.0005310529/0.0004728301$ for plus. Their centre clocks are different:
$3.0507894717/1.9788561960$ and $0.7064165019/1.1882379468$ respectively.
At the plus occupied $T=3$ clock, the matched empty scalar time is only
$0.6905540716$; station interpolation gives empty length $0.8703609817$
and centre radius $0.8098710354$, compared with occupied radius
$0.5385445337$. These sparse geometry interpolations are finite indicators.
The minus occupied $T=3$ clock exceeds both the frozen endpoint
$2.9995269273$ and the empty endpoint $1.9788561960$, so all matched
control entries are null there.

Total energy stays near $10^{-11}$ in these recorded endpoint readouts,
while constraints and projected curvature grow: the plus occupied
$T=3$ raw C maximum is approximately $46.9561$ and projected $|R_4|$
maximum approximately $5.30275\times10^7$. The empty plus geometry also
reaches a small lapse and large curvature diagnostic. There is no
propagated observable error bound, continuum instability verdict or
autonomous-renewal proof. Proper-clock control differences establish finite
readout differences within the saved domain.

The next physical question is whether this source-dependent transfer
persists at matched proper clocks under resolved spatial/time refinement.
