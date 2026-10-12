# Atlas v5: a causal map around an interior observer

## Route and admission

The requested outcome is one reproducible spacetime diagram and observer
access table, carried into the working paper while preserving the other atlas
routes. The v4 body, source and initial slice are fixed inputs. All v1–v4
source-bound files remain unchanged; this version supplies successor files.

Three routes were considered: import the isolated spherical wall's known
causal completion; solve a global many-body Einstein–shell evolution; or start
a new communicating-mouth model. The first closes the current observer
question with already available mathematics. The second needs unspecified
connector and evolution data. The third remains an independent atlas route.
No global multi-body evolution is admitted by this milestone.

The [classification owner](nsc-atlas-causal-classification.md) matches the
fixed signed branch to published thin-wall results. The
[observer owner](nsc-atlas-causal-observers.md) attaches central, wall and
right-parent worldlines, with event-specific access. A thin reference adapter
supplies regular Schwarzschild null coordinates and the v4 slice embedding.
The successor producer combines their records and draws matched coordinate
panels; the [v5 paper](../../paper/atlas-working-draft-v5.md) explains which
connections are established, blocked in this reference, or still unknown.

The root reference calculation borrows the signed junction equation from
[Aguirre–Johnson, gr-qc/0508093v2](https://arxiv.org/html/gr-qc/0508093v2),
Section II, and Schwarzschild null coordinates from
[Kuchar, gr-qc/9403003v1](https://arxiv.org/html/gr-qc/9403003v1), Section 2.1,
equations 9–13. The missing connection is the location of our existing slice
and observers within that known reference history. It decides allowed null
paths and unlocks a physically specified response question. It does not
rederive a new spacetime classification or formation mechanism.

Before compute, fix H ell=1, k ell=2, M/ell=0.4, R_wall/ell=1.5 and the expanding
branch. Use at most 300 plotted wall samples on 0.02<=R/ell<=6, fixed 32/64-node
quadrature, and a 60 CPU-second budget per adapter. Analytic signs and imported
end limits, rather than the finite plotting cutoff, own global statements.
Stop if the signed branch, regular null embedding or independent checks fail.
Do not change the source to rescue an access verdict. A failure leaves the
previous atlas and independent alternatives visible.

## Regular reference coordinates

Put x=R/ell. The imported signed trajectory gives

$$
\beta_-={1\over5x^2}+{3x\over4},\quad
\beta_+={1\over5x^2}-{5x\over4},\quad
s^2=\dot R^2={25x^2\over16}-1+{3\over10x}+{1\over25x^4}.
$$

The expanding branch has s>0. Its future-directed Schwarzschild advanced
coordinate v=T+r* has the horizon-regular derivative

$$
{dv\over dR}={\beta_++s\over f s}
=-{1\over s(s-\beta_+)},\qquad f=1-2M/R.
$$

The second expression extends through R=2M for the selected negative-beta
horizon crossing. Take V_K=-exp(v/(4M)) and
U_K=(1-R/(2M))exp(R/(2M))/V_K. The signs select the white-hole sector before
the crossing and the left exterior afterwards. Both Kruskal null coordinates
increase to the future. Their monotonicity, product identity and wall metric
normalization are independently checked. The integration origin is set by
the v4 slice's Schwarzschild T=0 bridge, not by pretending that all regional
clocks agree.

For a finite early-wall witness R/ell=0.5, U_K<0 and V_K<0. A constant-U_K
outgoing radial null path reaches the right static parent observer at
R/ell=1.6. A signal emitted by that parent has V_K>0 and cannot reach the wall
on this branch, whose V_K stays negative. These are reference-spacetime
statements. Post-cut and earlier-history rows are recorded separately.

## Completion criteria

The [recorded access table](../results/development/nsc-cosmological-atlas-v5/observer-access.md)
contains early centre/wall and centre/parent null witnesses, along with
exclusions for new emissions at or after the chosen slice. Earlier wall
signals can arrive at the centre arbitrarily late; emission and reception
times must not be conflated. The parent escape threshold is
t_C/ell=-1.2312944525 in the recorded convention. At equality, no finite
static-parent reception occurs. This escape horizon is distinct from the
central observer's incoming de Sitter horizon.

The [matched charts](../results/development/nsc-cosmological-atlas-v5/spacetime.svg)
clip horizon segments at their corresponding wall events. C's incoming
horizon continues into the Schwarzschild panel along constant V_K. The
uncut past R=0 chart edge is shown faint and dashed because part is excised;
finite wall endpoints are labelled rather than drawn as computed infinities.
The [map](../results/development/nsc-cosmological-atlas-v5/map.svg) keeps the
multiple-daughter, communicating-mouth and shared-parent routes visible.

Reproduce with `python lab/scripts/derive_nsc_cosmological_atlas_v5.py`
in the scientific environment described in the repository instructions.
The producer creates a new bundle and refuses to overwrite one. Add `--check`
to replay the existing record; use `--output` only for a fresh directory.
The manifest binds all imported v1-v4 sources plus this version's adapters,
tests, notes and paper. The focused test owners are
`test_nsc_atlas_causal_reference.py`, `test_nsc_atlas_causal_observers.py`
and `test_nsc_cosmological_atlas_v5.py`. Null-coordinate finite differences,
independent boost formulas, proper-clock agreement, quadrature order, scale
changes and domain exits provide integration checks. The prior atlas tests
and replay checks retain the earlier incompatible and unresolved branches.

- The known branch is selected by equations and signs, including the
  Schwarzschild exterior limit; no diagram is selected by visual resemblance.
- Named observer worldlines, emission events, clock origins, null paths and
  excluded or unknown domains accompany every access verdict.
- The diagram includes the initial slice and distinguishes a conformal chart
  from a Euclidean embedding or a computed many-body future.
- Independent null-coordinate and unit-rescaling checks pass; the observer
  table and paper agree with the recorded inputs and proof boundaries.
- Previous atlas replay, source hashes, publication links and affected checks
  pass; frozen sources and unrelated work are preserved.

## Failure forecast

1. Squared dynamics can select the wrong exterior. A mismatch in signed beta
   or Kruskal quadrant stops that attachment before it reaches the figure.
2. Coordinate singularities or shifted clock origins can create false access.
   Regular null coordinates, independent metric normalization and explicit
   event labels detect this; coordinate time is never a cross-universe age.
3. A reference diagram can be mistaken for the glued parent's global future.
   Records keep separate scope fields and unknown many-body rows. The v4
   protected smooth-vacuum domain remains its own finite certificate.
