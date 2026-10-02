# Same-action initial expansion and auxiliary-sector controls

Owner: [module](../src/recursive_horizons/nsc_discovery_initial_sector.py),
[thin driver](../scripts/derive_nsc_discovery_initial_sector.py), and
[tests](../tests/test_nsc_discovery_initial_sector.py).
No new scientific trajectory is produced by the implementation or its tests.

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
