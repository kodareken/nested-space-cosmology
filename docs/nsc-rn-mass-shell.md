# The mass shell of the RN radial parent

This is the first reduction in this line that measures the geometry rather
than a force balance.

Everything else in the RN family has been a force question: what the radial
and angular pressures do, where they are equal, how far each one reaches.
Those are real measurements, but they describe what the forces are doing,
not what the object looks like.

In the PG chart the geometry is fixed by the mass function `M(r)`. So the
shape of `dM/dr` **is** the shape of the object, and it can be read straight
off the recorded stations with no model, no fit, and no reference solution.

## What was found

A coherent shell of accreted mass, present in all six cases, which:

| | measured | spread across the family |
|---|---|---|
| starts at | `r = 12.010 .. 12.121` | `0.11` |
| falls at | `-1.4417 .. -1.4527` per unit `t/r_m` | `0.011`, or `0.8%` |
| loses prominence at | `t = 5.5` | identical |
| reaches the inner edge at | `t = 7.0` | identical |
| carries amplitude scaling with `eps` | factor `10.007` per decade | linear to `0.07%` |

The linear fits are good: `r^2` between `0.99865` and `0.99917`, with a
largest residual of `0.15` against a fall of `7.3`.

Within the falling phase the shell accelerates gently and monotonically,
from about `-1.344` to `-1.613`. It is not a straight constant-speed body;
the linear fit is a good summary, not an exact law.

## The two knobs do different things

The mass ratio `mu` and the source strength `eps` do not act on the same
quantity here.

`eps` scales the shell amplitude exactly, decade for decade. Doubling the
source buys twice the shell.

`mu` changes the mass geometry by almost nothing: the infall rate varies by
`0.8%` across a `1.002` to `1.040` range in `mu`. The shell arrives at the
inner edge at `t = 7.0` in every case.

**This does not mean `mu` is irrelevant.** `mu` moves the radial and angular
pressures, and therefore the balance wall and the force reach, both of which
were measured earlier and both of which do respond to it. What is insensitive
to it is the mass geometry alone. The family is separating a force question
from a mass question rather than sweeping one dial, and the two must not be
conflated when the wall and the shell are discussed together.

## What was corrected on the way

Three earlier statements from this same reduction were wrong and are recorded
here rather than quietly dropped.

**The first fit reported an infall rate that depended on `mu`.** It did not.
The fit ran through the point where the shell reached the inner edge and
reappeared later as a different, much weaker feature. Reading a mixture of
two things as one trajectory produced a `mu` dependence that is an artifact.
Restricting the fit to the contiguous falling phase removes it entirely.

**A bare sign change in `dM/dr` was reported as a shell.** Every recorded
reversal in the family sits at or below the median `|dM/dr|`. A shell is a
peak in the mass-per-radius with a width, not a sign flip in the noise. The
estimator now requires the excursion to clear the noise floor.

**`edge_arrival_t` was mislabelled.** It was recording when the shell lost
prominence, which is `t = 5.5`, not when it reached the inner edge, which is
`t = 7.0`. The field is now split so neither is called the other.

Records `v1`, `v2` and `v3` are retained and sealed. They contain those
errors. `v4` is the current record.

## Where this leaves the picture

The shell is the part of the object you can see. It has an envelope, a
width of roughly `4` to `4.6` at half maximum, and it falls coherently.

Where it goes is not measured, and that is the significant limit. The
computation's inner edge sits at `r = 1.04`, and at the end of every run
between `61%` and `99%` of the accreted mass lies **below** that edge. The
accumulation, which is the part of the picture that matters most, happens
outside the domain.

This also explains two earlier puzzles:

- the throat looks steep and never becomes regular, because the grid's inner
  boundary cuts across the middle of an inflow rather than sitting outside it;
- the oscillation audit found grid-scale ringing at large radius, because
  there is very little mass out there. The far field is nearly empty, so
  ringing there is numerical noise on a negligible amplitude.

Neither is a property of the object. Both are properties of where the domain
was cut.

## Scope

This is a construction check on recorded finite states in a named collar. It
does not measure the actual universe, and a coherent shell falling into a
boundary is not evidence that any real object has one.

The shell's motion was fitted over `t = 0` to `t = 5.5`, which is `66` of
`177` recorded stations. The remaining `111` stations are edge-pinned
boundary structure or the late weak feature, which is `1.7%` to `3.3%` of
the shell's amplitude and is a different thing.

## Provenance

- module: `lab/src/recursive_horizons/nsc_rn_mass_shell.py`
- driver: `lab/scripts/derive_nsc_rn_mass_shell.py`
- record: `lab/results/development/nsc-rn-mass-shell-v4.json`
- tests: `lab/tests/test_nsc_rn_mass_shell.py`, `24 passed`
- source family: `nsc-rn-neutral-family-v1`, `6` cases, `177` stations,
  reduced with no time step taken

## What would settle the open question

The shell arrives at the inner edge at `t = 7.0`. The domain ends at
`r = 1.04`. Nothing in this family says what the shell does next, and that
is exactly where the geometry of the picture lives.

To see it, the march has to start inboard of where the shell crosses,
with the inner boundary condition applied at a radius smaller than `1.04`,
and the shell released inside the domain rather than arriving from outside.
That is a different construction from the one recorded here, and it is the
next thing to build.