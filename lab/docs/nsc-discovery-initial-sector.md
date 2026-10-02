# Same-action initial expansion and auxiliary-sector controls

Owner: [module](../src/recursive_horizons/nsc_discovery_initial_sector.py),
[thin driver](../scripts/derive_nsc_discovery_initial_sector.py), and
[tests](../tests/test_nsc_discovery_initial_sector.py).
No new scientific trajectory is produced by the implementation or its tests.
The sealed v1 production has reached T=3 for all six cases; root owns those
observations. Its producer `cafb976cc2bbb3190ae8e258c481982a437d1d46`
and all v1 preparation/episode bytes remain unchanged by the successor below.

## Question and inherited action

All existing uniform-Q, chi-zero, zero-momentum preparations start with
Qdot=rdot=chidot=0 and Qddot=-Q^3. Their future contraction is evidence for
those initial sectors, not a broad negative result about other Cauchy data.
This bounded comparison changes existing initial geometric degrees of freedom
under the same Gamma_one action. It changes no coefficient, cosmological term,
mass, angular copy count, field covariance, observer, or physical force.

The input is the uniform case of the frozen PREPARED-v2 record,
`nsc-discovery-dynamic-preparation-v2.json/.npz`, produced by
`67ab5700229216612e7383dc23d0ae1ce6e479c0`. Its actual real standing
eigenvectors, order, signs, six weights, native W and observer are loaded;
no new degenerate eigenframe is selected. The completed uniform h=0/chi=0
cases in `nsc-discovery-dynamic-episode-v1` are authenticated read-only
comparisons. They are not copied into a new claim or evolved again.

For a=-8 pi A, f=-8 pi C_W/3 and Pi=p_r-ar p_chi/f, the family

\[
p_\chi=2fh,\qquad p_r=2arh,\qquad p_Q=0
\]

has Pi=0. Uniform geometry and the measured zero-current source make the
initial lapse and shift kinetic contributions zero. The rates are
Qdot=hQ, rdot=chidot=0. The normal radial expansion is
K_r=h/(rQ), with K_perp=0; changing h is not merely a lapse change with
the future normal and observers fixed. Opposite h may represent opposite
time directions of one unoriented solution for real time-reversal-invariant
matter; the future-directed experiments remain distinct.

The declared initial h values are -0.5 and +0.5. The h=0 comparison is the
saved case. h is not selected by a new law and is not held constant later.
Every evolution stage calls the existing full coupled rates. In particular,
homogeneous hdot=-Q^2(2+chi)/2 initially, so these are not ongoing pumps.
Nested momenta are pi=dx_g W^T p, not the nodal densities above.

## The chi=-2 control

Write alpha=-4 pi C_W/3 and mag=2 pi C_F flux^2. The initial constraint
with zero p_Q and Pi selects

\[
r_*^2=\frac{\rho/Q+\mathrm{mag}+\alpha(4\chi+\chi^2)}{8\pi A}.
\]

Thus chi=-2 gives r_*^2=(rho/Q+mag-4alpha)/(8 pi A). The locked coefficients
give mag-4alpha>0. This radius is computed from the actual frozen source;
the chi=0 radius is not copied. The existing general-Q sampled radial
Jacobian supplies a bounded positive finite correction. Full, projected,
and held-out constraints and the undeleted mean current are reported.

At chi=-2/h=0, p_chidot=Qddot=0 initially, but p_Qdot and p_rdot generally
do not vanish. In fact rddot=-rQ^2/3, so removing the forced Q turn does
not make proper radial length stationary. The actual metric has initial
R_h=0 and C^2=4/(3r^4), whereas the h controls retain R_h=2 and initial
C^2=0. These statements are checked with the actual projected time jets,
not by inserting chi into the curvature formula. h changes initial
extrinsic data; chi=-2 changes the invariant-curvature sector.

The source energy and geometric kinetic budget are recorded separately.
The angular multiplicity M=4 kappa is applied once inside the source;
it does not multiply the representative field generator. CAR, source
trace, current and source/observer/W hashes are retained. Algebraic
constraint compatibility and finite arithmetic residuals are not a continuum
initial-state certificate or evidence of holding.

## Bounded thin episode

There are three NEW preparations: h_minus, h_plus, chi_minus2. Each forks
its exact stored Cauchy state into caps 0.001 and 0.0005, yielding six
cases. Stations are common coordinate T=0.3,1,3. The source and observer
are never reset. Initial-clock differences and initial offsets remain
explicit; comparisons do not claim matched proper times.

The runner reuses the frozen period-aware dynamic-preparation adapter,
the existing field/geometric principal timestep restriction, coupled RK4,
actual curvature/tidal observer, disjoint instantaneous normal ledger,
proper clocks, checkpoint publication and locked aggregate CPU ledger.
One existing pool uses at most four workers, forecast factor 1.5, and a
fresh maximum 300 CPU seconds. Preparation and observations are charged.
Budget/chart/diagnostic stops retain the last admissible state. There is
no new framework, inhomogeneous h solver, restoring force, constraint reset,
or claim that a startup force supplies a maintained structure.

The primary source remains spatially uniform at preparation. This comparison
tests a restricted initial-rate/curvature assumption; it does not itself
prepare a localized parent-child structure. Later outcomes can support a
named finite question, not a broad conclusion about every allowed initial
state or an infinite future.

## Commands after the root freezes the producing source

Let `SCIENCE_SHA` denote the new commit containing all four owning files and
the observed dependency closure. The preserved producer 67ab570 remains the
input producer, not the new initial-sector producer.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_initial_sector.py -q
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --prepare --producer-commit SCIENCE_SHA --write results/development/nsc-discovery-initial-sector-preparation-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --check results/development/nsc-discovery-initial-sector-preparation-v1.json
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --episode results/development/nsc-discovery-initial-sector-preparation-v1.json --execute --producer-commit SCIENCE_SHA --write results/development/nsc-discovery-initial-sector-v1
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --run-episode results/development/nsc-discovery-initial-sector-v1 --workers 4 --episode-cpu 300
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --assess-episode results/development/nsc-discovery-initial-sector-v1
```

No production execution is authorized by importing the module or reading this
note. The root coordinates scientific execution after the source freeze.

## Conditional source-balanced curvature, balanced-v2

The authorized successor is one source-conditioned initial curvature choice,
not another arbitrary h/chi scan. It uses the SAME frozen uniform Q, actual
field columns, six occupations, observer, native geometry map and coefficients.
All initial momenta remain zero. The other two cases offset chi-star by exactly
plus/minus one percent of its magnitude; they are robustness comparisons and
are not called balanced preparations. Each source/chi choice has its own
constraint-selected positive radius and actual finite constraint trace.

For the uniform initial slice, with S the weighted canonical spin density
and rho=M(K/Q+kappa S), differentiating the owned chi rate gives

\[
\ddot\chi=\frac{2Q\rho-M\kappa QS+(8\pi A/3)r^2Q^2\chi}{4\alpha}.
\]

Substitution of the owned lapse constraint reduces chi-ddot=0 to

\[
P(\chi)=\alpha Q\chi^3+4\alpha Q\chi^2
 +(\rho+\mathrm{mag}Q)\chi+6\rho-3M\kappa S=0.
\]

The locked alpha is positive and mag>16alpha/3, so

\[
P'(\chi)=3\alpha Q(\chi+4/3)^2+
 \rho+Q(\mathrm{mag}-16\alpha/3)>0.
\]

This proves uniqueness under the measured positive-source hypotheses. The
vacuum limit rho=S=0 recovers chi=0. On the cached uniform source,
Q=0.463944714431183, rho=3.43433694084887 and S=0.203289532850862,
the unique chi-star is approximately -4.58327910764894, with radius
2.74786461178492. These are conditional preparation values, not a fitted
fundamental law. Q and zero initial momenta remain declared inputs.

The chi-zero matter preparation instead has initial chi-ddot approximately
198.6406658465. That was valid initial data with the vacuum auxiliary value,
not a demonstrated matter-selected curvature. The balanced-v2 rule cancels
ONE initial acceleration. It does not impose later curvature, match higher
adiabatic jets, select an entire dynamical solution, or establish holding.
No mass, source occupation, eta, action term, restoring force, or pump is added.

Every new case records the polynomial and coefficients, global positive
derivative lower bound, actual rho/S/M/kappa factors, its radius/constraint
correction, actual chi acceleration from the full analytic rate Jv, independent
metric tides/curvature, CAR, source trace, kinetic budget, canonical momenta and
normal clocks. The perturbation cases keep their measured nonzero initial
chi acceleration. Source selection and phase/order are not repeated after
handoff. Both caps fork exactly the same state in each case.

The v1 checker authenticates its recorded producer and executes the matching
Git module bytes when the live module has advanced. It checks that the
physical dependencies are still unchanged; it does not heal a v1 source
binding, rewrite a payload, or relabel v1 data with the successor producer.

The new immutable stems are
`nsc-discovery-initial-sector-preparation-balanced-v2` and
`nsc-discovery-initial-sector-balanced-v2`. One pool, at most four workers,
fresh 300 CPU seconds, forecast factor 1.5, caps 0.001/0.0005 and stations
T=0.3,1,3 are inherited without another runner. No production is performed
during development.

If all three cases show curvature departure growth or lose the positive chart,
report a limitation of this finite conditional initial-sector family and stop
this route. Do not pursue higher-jet matching, a ghost stable manifold, alpha
thresholds, or fine tuning. The consumer reports the actual retained chi
departure maxima and chart stops; a strictly increasing departure across
retained positive-time stations is a finite trend, not an infinite instability
proof. This route policy does not veto physical growth or make a broad theory
claim. Uniform-source comparisons do not themselves establish a localized
structure, a black-hole bounce, or an eternal engine.

After root pins the NEW source commit (not caf), execute explicitly:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --balanced
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --balanced --prepare --producer-commit SCIENCE_SHA --write results/development/nsc-discovery-initial-sector-preparation-balanced-v2
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --check results/development/nsc-discovery-initial-sector-preparation-balanced-v2.json
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --episode results/development/nsc-discovery-initial-sector-preparation-balanced-v2.json --execute --producer-commit SCIENCE_SHA --write results/development/nsc-discovery-initial-sector-balanced-v2
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --run-episode results/development/nsc-discovery-initial-sector-balanced-v2 --workers 4 --episode-cpu 300
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_initial_sector.py --assess-episode results/development/nsc-discovery-initial-sector-balanced-v2
```

The baseline criterion is initial chi-ddot=0; nearby chi values do not satisfy
that criterion. Subsequent growth remains an observed outcome of the same
action, not evidence that a reset or another source force should be inserted.

## Recorded outcome of the bounded comparison

The balanced successor is produced by `65ebad8634c2c486e67761c2161346e9030d5748`.
All six cases are retained in
`results/development/nsc-discovery-initial-sector-balanced-v2/manifest.json`.
Aggregate charged CPU is 137.43 seconds. The conditional candidate reaches
T=3 with child proper length 0.0023929, down from 2.5497, and mean areal
radius 9.8771. The upper-curvature control also contracts; the lower control
reaches the positive-radius chart edge near T=1.482. Its two step caps bracket
the same event within 0.001 and 0.0005 respectively and retain their last
admissible states. No variable is clamped and no source is reset.

These finite outcomes trigger the declared route stop. More curvature jets,
coefficient tuning or a finer root search will not be used to turn this
preparation into a holding claim. The action approximation and its additional
curvature sector require a separate assessment; the earlier transfer and
reduction results remain valid in their stated numerical domains.
