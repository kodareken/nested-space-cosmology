# The isotropy trajectory and the registered relaxation test

Owning note for
[nsc-rn-isotropy-v1.json](../results/development/nsc-rn-isotropy-v1.json) (the
nine-case driven window) and its successor
[nsc-rn-isotropy-v2.json](../results/development/nsc-rn-isotropy-v2.json) (the
relaxation campaign, read against the registration v1 wrote first).
One finite question in one named domain: how far from its isotropic point does a
recorded state sit, and how does that departure move with time?

## What is measured

The two forces that hold a region are read from the recorded pressures. Their
balance is isotropic at `B(r) = 1/2`, and the reduced scalar is the mean
departure over the recorded cells,

    I = <|B(r) - 1/2|> .

Both recorded layouts are reduced: the 177 family stations and the 12
inward-march stations (grid one cell above the inner horizon) — 9 cases, 189
stations, from the same owning module that produced the two-force balance. The
successor record adds the 6 stations of the relaxation campaign
(`nsc-rn-relaxation-v1`, μ=1.040, ε=0.001, t = 5, 10, 15, 20, 40, 80), which is
the run the registration below named.

## What the recorded states say

| Layout | Case | I first | I last | registered verdict |
|---|---|---:|---:|---|
| family | μ=1.002, ε=0.001 | 0.4675 | 0.4290 | no relaxation |
| family | μ=1.010, ε=0.001 | 0.4675 | 0.4091 | no relaxation |
| family | μ=1.040, ε=0.001 | 0.4686 | 0.3939 | no relaxation |
| family | μ=1.040, ε=0.010 | 0.4686 | 0.3822 | no relaxation |
| march | μ=1.002, ε=0.001 | 0.4467 | 0.4060 | no relaxation |
| march | μ=1.040, ε=0.001 | 0.4447 | 0.3839 | no relaxation |
| relaxation | μ=1.040, ε=0.001 | 0.4447 | 0.3454 | no relaxation |

The full record carries all nine cases. Every one of them **falls
monotonically** from first to last station — the state does move toward the
isotropic point as the shell falls — but the fall is 10 to 18 percent, far short
of the registered half, so the registered verdict is `no_relaxation_observed`.
That is the honest reading of the driven window: a trend in the predicted
direction, not a verdict.

The relaxation campaign extends the same case to $t=80$ and **does not fall
monotonically**: the index runs 0.4447, 0.3966, 0.3934, 0.3839, then **rises to
0.4191 at $t=40$** before falling to 0.3454 at $t=80$. The first-to-last fall is
22.3 percent, larger than the driven-window trend, but the tail still sits above
the registered half of the driven maximum (0.3454 against a threshold of
0.2223), so the registered verdict is again `no_relaxation_observed` — with the
mid-run rise recorded as measured rather than smoothed away. This remains a
driven window; nothing here is an undriven late-time sector.

## The registered test

The relaxation test was written into the record before the run existed: a tail
index at or below **half** the driven-phase maximum, and a fall larger than
**three times** the largest recorded refinement indicator, is
`relaxation_observed`; a fall inside the refinement multiple is
`inconclusive_on_this_domain`; anything else is `no_relaxation_observed`. The
campaign `nsc-rn-relaxation-v1` (μ=1.04, ε=0.001, stations to t=80) has now been
run and read against that registration, and the verdict is
`no_relaxation_observed`. The registered refinement indicator is `0.0` on this
layout — no refinement channel was recorded — so the second clause is
quality-of-evidence only and the verdict turned on the first: the tail is not at
or below half the driven maximum.

## Cost, measured rather than assumed

One two-arm march step on the extended grid takes **0.12 s** single-threaded
(`driver.march_pair`, 876 points, dt = 0.001, measured on the recorded t=20
state), so stations to t = 80 cost about 2.7 h of wall clock. The first launch
of this campaign was stopped on a projected cost, not on a measured one: a
process-CPU snapshot was divided by an assumed elapsed wall time, giving about
0.6 s/step and roughly 14 CPU-hours to the endpoint, over the six-CPU-hour
batch rule. The direct benchmark then measured 0.12 s/step, inside the rule,
and the campaign was relaunched. The stopped launch had kept no station — only
its prepare marker, which was removed — so nothing measured was lost. The
relaunched campaign completed all six stations and is the run reduced into
[nsc-rn-isotropy-v2.json](../results/development/nsc-rn-isotropy-v2.json) and
read against the registration above.

## Arm control, measured

The frozen-metric arm is the same action without geometric feedback, and its
pressures can be recomputed from the stored `fixed_phi` on its own frozen lapse
and shift with the same statistic
([nsc-rn-isotropy-arm-control-v1.json](../results/development/nsc-rn-isotropy-arm-control-v1.json)):
the two arms agree to **0.0014** worst-case per station across all 189 stations,
and both read `no_relaxation_observed` on the driven window. The same control on
the relaxation campaign
([nsc-rn-isotropy-arm-control-v2.json](../results/development/nsc-rn-isotropy-arm-control-v2.json))
agrees to **0.00127** worst-case per station (largest gap at $t=80$), and it
reproduces the mid-run rise as well: the frozen arm runs 0.4447, 0.3966, 0.3936,
0.3843, 0.4189, 0.3467 against the coupled arm's 0.4447, 0.3966, 0.3934, 0.3839,
0.4191, 0.3454. The bump is therefore in the field, not in the feedback.

So the fall is a property of the evolving field, not of the geometry's response:
the feedback that distinguishes the arms moves the index by three orders of
magnitude less than the trend. That also fixes what this instrument can claim.
The RN chart is spherically symmetric by construction — its metric carries no
shape anisotropy at all — so the index measures the **matter** stress, which is
the only anisotropy these recorded states possess. A test of shape selection
itself needs a non-spherical sector these runs do not contain: $l\ge2$ metric
perturbations on the same background. That is the instrument the question now
names, and it is a new calculation rather than a reduction of what is on disk.

## Declared limits

The index measures the **matter** balance, not the metric's anisotropy; the
metric reading is the four-metric tidal comparison of the discovery family, and
this module does not replace it. Every recorded station is driven — the shell is
still falling — so the trend above is a statement about the driven phase only.
The record is a development record; it claims nothing about the actual universe.

## Owners

Module [`nsc_rn_isotropy.py`](../src/recursive_horizons/nsc_rn_isotropy.py),
driver [`derive_nsc_rn_isotropy.py`](../scripts/derive_nsc_rn_isotropy.py),
tests [`test_nsc_rn_isotropy.py`](../tests/test_nsc_rn_isotropy.py). The balance
profile comes from [`nsc-rn-gradient-bubble.md`](nsc-rn-gradient-bubble.md).