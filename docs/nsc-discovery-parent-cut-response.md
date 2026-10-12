# Exterior-population response at protected cuts and child proper clock

This thin adapter uses the authenticated NF128 strong balanced MINUS parent,
`k=.45388971484879426`, strength `173.16013550038755`, with constant magnetic
initial radius and no geometry fallback. It changes only the parent-annulus
column's occupation: `delta c=(0,c1)`. The child column weight, initial `Q,r`,
canonical fields, source identities and selected `k` remain fixed. This is a
new source state; weighted-Gram CAR is checked without clipping.

Initial momenta must change to satisfy the actual finite constraints. The
adapter uses the parent's analytic momentum Jacobian and hard scalar anchor.
Its source rows are the actual `Q delta rho`, `delta j`, and negative derivative
of the SAME projected fixed-`k` seed anchor. `delta J` is derived from the
actual finite periodic primitive; constant magnetic radius makes it tiny but
it is not set to zero. Fixed `k` does not mean fixed anchor. Two centered
momentum preparation widths check the implicit derivative. Initial source
fields and geometry are never rebuilt. A singular analytic rank or unresolved
perturbed momentum root is a preparation obstruction, not an exact tangent.
The old radius-only response preparation helper is not used.

Protected cuts are explicitly `x=3,5`, with outward normals `n=-1,+1`.
The distinct child window is `[3.5,4.5]`. Incoming geometric diagnostics are
`u_t+n u_x` for `u=Q,r`, using projected velocities; the radial derivative is
normalized by `N=rQ`, while `Q` is labeled a conformal-chart diagnostic.
The Dirac incoming projector is `(I-n sigma2)/2`; cut values use the native
half-integer AP Fourier carrier, in canonical density units, without a periodic
spinor wrap. Its weighted quadratic current derivative includes both `delta Phi`
and `delta c` and the proper-clock factor. Child spatial content, proper length,
center clock rate and cumulative proper clock variations are also reported.
The center-clock view changes only the clock list to `(4,)`.

Retarded evolution calls the existing analytic response owner and shared RK4
base/tangent stages repeatedly, with the rank-general admission cap `.001`.
The returned `tau` and `delta_tau` are INCREMENTS and are accumulated.
Stations are coordinate times `0,.25,.5,1,1.5,2.25`. Readouts match the center
proper clock with `delta O|tau=delta O|t-Odot delta_tau/tau_dot`.
Initial constraint-generated response is saved explicitly, and evolved change
is measured relative to it. Finite-band tails prevent an exact compact causal
delay claim; no expected-null result is presumed.

The forecast locks both `alpha=+.05` and `-.05` predictions at every attained
baseline proper clock BEFORE either held nonlinear arm. Root selects the
`T=2.25` baseline clock for the first comparison. Held arms independently
reprepare actual finite momenta at fixed `k` and then integrate to that proper
clock. Small amplitude controls may be run only after the forecast lock.
This is one source-to-incoming-to-child response comparison, not a imposed
traction, new dynamics, static wall, gate or full renewal claim.

After freezing the four adapter files and existing producer closure, root runs:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_response.py --prepare
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_response.py --predict
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_response.py --run --station 2.25 --amplitude .05
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_response.py --run --station 2.25 --amplitude -.05
```

Default operation authenticates and previews without preparation/evolution or
writing. Stages reuse the existing response owner's exclusive JSON/NPZ writer,
source closures and readers, with a fresh designated output prefix and 64 MiB
limit. The first batch has an aggregate 300 CPU second allowance, including
preparation, forecast and completed measurement arms. Budget stops retain finite
endpoint payloads and do not compare an unmatched clock as a measurement.
Numerical comparisons remain indicators: continuous error, physical vacuum
identification and strong-curvature EFT validity remain open.
