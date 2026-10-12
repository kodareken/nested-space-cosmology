# Atlas v6: an early field imprint at the interior observer

## Decision and calculation admission

The requested outcome is a specified field, initial state and wall transmission
law, followed by a calculated interior measurement and an updated paper/map.
The common scientific question is how regional geometry, boundaries and state
shape local measurements. Spherical charts make this first connection
tractable; they do not select a globally spherical many-body cosmology.

Three routes were considered. An exact conformal test scalar closes a finite
wall-to-centre transfer. A minimally coupled scalar adds a bulk tail and a
different reduced wall jump, so it remains a separate candidate. A global
Schwarzschild scattering calculation needs additional exterior initial data
and cannot import the usual black-hole boundary conditions unexamined into
our white-hole branch. Geometric optics would add an unnecessary approximation
to the first route. No new PDE or quantum-vacuum calculation is admitted here.

The fixed four-dimensional classical field has action

$$
S_\phi=-\frac12\int\sqrt{-g}\,[(\nabla\phi)^2+\tfrac16\mathcal R\phi^2]d^4x.
$$

It is massless, spherical (angular mode zero), and treated at linear probe
order on the frozen v5 geometry. There is no extra scalar interaction on the
wall. The amplitude epsilon is a formal small parameter, not a fitted
cosmological energy density. Neither a quantum vacuum nor backreaction is
selected by this classical calculation.

- **Borrowed mathematics:** the distributional scalar jump in
  [1512.02120v2](https://arxiv.org/html/1512.02120v2), equations 31-34, and
  [1901.07524v1](https://arxiv.org/html/1901.07524v1), equations 1, 4 and 15;
  conformal field transformation in
  [gr-qc/0603031v1](https://arxiv.org/html/gr-qc/0603031v1), equations 21, 24,
  30 and 36; and the fixed v5 wall and observer charts. The imported papers'
  different static/cylindrical geometries and coupling conventions are not
  substituted for our moving spherical wall.
- **Missing connection:** a geometric null path did not specify which field,
  incident state or matching law supplies the measured signal.
- **Decision and unlock:** obtain the conditional waveform, clock distortion,
  finite support and any remaining central field. This separates a transmitted
  transient from a persistent local source and identifies the data still
  needed to infer an exterior from a measurement.
- **State and finite domain:** zero classical field and time derivative on the
  cap at the retained wall radius R/ell=0.05. Prescribe one incoming
  characteristic on the exterior side of the wall. Its compact profile is
  fixed in wall proper time, centred at the early R/ell=0.5 event, with
  half-width 0.03ell. Observe through central proper time t_C=0. The incoming
  history is a boundary-state assumption; its remote preparation and smooth
  exterior scattering are not derived.
- **Budget and stop:** at most 513 wall samples, paired 32/64-point quadrature,
  and 60 CPU seconds per adapter. Stop if the pulse leaves the retained
  history, reaches R>=ell, violates junction/centre regularity, or lacks the
  claimed finite no-return domain. Keep the model fixed; changing the pulse
  or coupling begins a separately labelled input or candidate.
- **Verification and attribution:** independent junction identities, radial
  wave and regular-centre checks, proper-clock change of variables, analytic
  compact-pulse integrals, unit rescaling, coverage errors and numerical
  refinement. This applies known field and matching mathematics to a project
  interface; it is not a new field equation or an established cosmic signal.

## Implementation and evidence owners

The [matching note](nsc-atlas-field-matching.md) owns source reading,
distributional signs and the limits of transmission. The
[field-transfer note](nsc-atlas-field-transfer.md) and
[adapter](../src/recursive_horizons/nsc_atlas_field_transfer.py) own the state,
arrival map and finite observer record. They reuse v5 wall clocks, with no
new solver dependency. The
[successor producer](../scripts/derive_nsc_cosmological_atlas_v6.py) preserves
all earlier candidate scenarios and adds this field-conditioned connection.
The [v6 paper](../../paper/atlas-working-draft-v6.md) connects it to the
foundation and the focused companion without identifying the companion's
unresolved mode-space complement with a geometric parent universe.

Write new versioned evidence under
`lab/results/development/nsc-cosmological-atlas-v6`. All v1-v5 bound sources,
outputs and papers remain unchanged. Current indexes and the literature
library point to the successor; their publication and import receipts are
refreshed after the scientific sources are frozen.

## Recorded result and reproduction

The fixed pulse occupies wall radii
`R/ell=[0.4749900536,0.5228504193]`. Its central arrival support is
`t_C/ell=[-0.5059339731,-0.3479958418]`: the duration is 2.6323021882
times the wall pulse duration. The normalized central extrema are
`ell*phi_C/epsilon=+23.8414479985,-24.4653347800`. Their unequal
magnitudes follow from the changing clock Jacobian, with no retuning of
the symmetric wall input. The field is exactly zero after the last arrival
through `t_C=0` for the stipulated incoming history. A separate no-return
bound leaves a comoving-radius margin of `0.1414661768ell` before any
reflected central pulse could revisit the wall during that window.

The adapter uses 440 distinct wall radii and 620,000 evaluations of the
borrowed clock integrands. Paired numerical rules and an independent radial
centre limit check agree; these diagnostics are not rigorous error bounds.
The 27 field tests and three successor-bundle tests pass. The first eight
atlas scenarios are preserved unchanged in the ninth-scenario successor.

With the scientific Python environment and `lab/src` and `lab/scripts` on
`PYTHONPATH`, generate a fresh directory with
`python lab/scripts/derive_nsc_cosmological_atlas_v6.py --output <fresh-directory>`.
Check the recorded bundle with
`python lab/scripts/derive_nsc_cosmological_atlas_v6.py --check`.
The [manifest](../results/development/nsc-cosmological-atlas-v6/manifest.json)
binds the sources, runtime versions, record, figures and narrative. Replaying
the known formulas is an integration check; the project's added result is
the explicitly conditional connection between the selected wall state and
this interior observer.

## Finish-line proof

- The action, orientation, surface law, classical initial data and the single
  incoming characteristic are explicit and mutually compatible.
- A reproducible waveform gives the named central observer's proper-time
  support and amplitude; finite late-time absence is distinguished from a
  statement about arbitrary fields or uncomputed exterior tails.
- Independent checks can falsify the transmission and clock normalization;
  numerical comparisons are not presented as rigorous error enclosures.
- The diagram keeps state preparation, wall matching and interior propagation
  separate. Other geometries and the minimally coupled field remain open.
- The paper's common thread is supported by the actual claims of the three
  manuscripts; no inference of universal topology, constants or dark-sector
  fractions is added.
- Source hashes, earlier atlas replay, local links and affected repository
  checks pass before a task-owned commit.

## Failure forecast

1. Omitting the wall's curvature distribution would give the wrong jump.
   Check the scalar and reduced-field conditions independently before drawing
   a transmission arrow.
2. A chosen incoming trace could be mistaken for a solved parent source.
   Keep the initial-boundary data in every record and leave exterior
   preparation/scattering visibly unresolved.
3. A coordinate amplitude or stored response could be called permanent
   cosmological memory. Record the observer's proper-time waveform and its
   post-pulse value separately from any future detector or state model.
