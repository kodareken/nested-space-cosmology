# Working method

Repo note, not a manuscript. How this programme steers its own work. Papers
carry the science; this page carries the steering, so it stays visible instead
of leaking into the prose.

## Rules

- The question is owned; the conclusion is not. Choose what to compute, never
  what it came out to be.
- Frame changes are dated and attributed to a record.
- Every view carries its exit condition.

## View changelog

| Date | View | Moved by |
|---|---|---|
| 2026-09-15 | Spectrum-first (foundation v0.26.0) | The operator, source, and boundary calculations of the first three months |
| 2026-10-06 | Geometry-first: the enclosure is the object; the spectrum is what the shape makes inheritable | [Cone map](../lab/docs/nsc-rn-cone-map.md), [two-force balance](../lab/docs/nsc-rn-gradient-bubble.md), [inward march](../lab/docs/nsc-rn-inward-march.md) |

The second entry supersedes the emphasis of the first, not its mathematics.

## Exit conditions in force

| Reading | Ends if | State |
|---|---|---|
| Sphere as the least demanding equilibrium | A non-spherical perturbation does not decay under the same action | Open; the driven window falls monotonically 10–18 percent in the *matter* balance, and the relaxation campaign to $t=80$ falls 22.3 percent but **not monotonically** (a rise to 0.4191 at $t=40$); both read `no_relaxation_observed` under the thresholds registered in `nsc-rn-isotropy-v1.json`. The frozen-arm control reproduces both (0.0014 and 0.00127), so the trajectory is field-driven, and the spherical chart has no metric anisotropy to test. The shape test needs an $l\ge2$ perturbation sector |
| One-way enclosure | An interior mode carries net outward current in a band | Strict for μ=1.010 on the horizon-reaching domain; a +0.283 low-band transient remains for μ=1.040 |
| Hayward regularity, neutral branch | A measured inner-region slope reaches $g\ge-1$ | Failed on the evaluated domain in all three cases |
| Ratio-structure inheritance | A finite family whose normalized spectrum does not transport | Not triggered |

Physics for the first two readings is in
[the programme statement](../paper/nested-cosmology.md), §2 and §8.