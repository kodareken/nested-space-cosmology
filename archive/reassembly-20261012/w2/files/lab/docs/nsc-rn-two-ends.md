# The two ends of a recorded room

Owning note for [nsc-rn-two-ends-v1.json](../results/development/nsc-rn-two-ends-v1.json).
One finite question in one named domain: are the bottom and top of a room
set globally, or does each recorded room carry its own two ends?

## Definitions

Read from the same recorded state's metric, so nothing is tuned:

| End | Definition | Meaning |
|---|---|---|
| Bottom | `omega_zero = N/r` | the angular ground of the Dirac shell, the cheapest note the room holds |
| Top | `omega_cutoff = (N + |beta|)/spacing` | the CFL resolution the march itself respects |
| Ratio | `Q = omega_zero / omega_cutoff` | dimensionless room quality |

`Q` is self-similar: under a coordinate doubling with the same
dimensionless lapse and shift profiles it is unchanged. Its closed form
is `N^2 spacing / (r (N + |beta|))`.

## What the six recorded cases say

| Case | Q inner (last) | Q inner / Q outer |
|---|---:|---:|
| μ=1.002, ε=0.001 | 3.9408e-03 | 19.909 |
| μ=1.002, ε=0.010 | 3.9299e-03 | 19.872 |
| μ=1.010, ε=0.001 | 8.7246e-03 | 19.687 |
| μ=1.010, ε=0.010 | 8.6871e-03 | 19.620 |
| μ=1.040, ε=0.001 | 1.6829e-02 | 18.920 |
| μ=1.040, ε=0.010 | 1.6758e-02 | 18.858 |

The room quality is not a single global number: the inner edge holds a
room about nineteen times sharper than the outer boundary in this
flat-slice gauge, and the inner value moves with the mass ratio from
`3.9e-3` to `1.7e-2` across μ. The energy fraction changes it only at the
fourth digit. The two ends are therefore set by the recorded room, not
read from a universal constant.

## Horizon-reaching domain (v2)

Reduced over the inward-march campaign (grid one cell above the inner horizon),
the room sharpens further:
[nsc-rn-two-ends-v2.json](../results/development/nsc-rn-two-ends-v2.json),
12 stations, ε=0.001.

| Case | Q inner (last) | Q inner / Q outer |
|---|---:|---:|
| μ=1.002 | 4.1734e-03 | 21.084 |
| μ=1.010 | 9.9791e-03 | 22.517 |
| μ=1.040 | 2.2380e-02 | 25.161 |

So the inner/outer ratio rises from the family's 18.9–19.9 to **21.1–25.2** once
the domain actually reaches the bottom, and the inner value itself grows (μ=1.04:
`2.24e-2` against `1.68e-2` on the old grid). The two ends are set by the local
space, and the sharper end sharpens as the measurement approaches the horizon.

## Declared limits

`Q` is a measured dimensionless invariant of the recorded states. It is
not Planck's constant and not absolute zero; the reading ("constants are
the other end of the same thing") stays a reading. The family never
reached a `summary.json`; coverage is reported per case.

## Owners

Module [`nsc_rn_two_ends.py`](../src/recursive_horizons/nsc_rn_two_ends.py),
driver [`derive_nsc_rn_two_ends.py`](../scripts/derive_nsc_rn_two_ends.py),
tests [`test_nsc_rn_two_ends.py`](../tests/test_nsc_rn_two_ends.py). The
metric reconstruction is shared with
[`nsc-rn-cone-map.md`](nsc-rn-cone-map.md).