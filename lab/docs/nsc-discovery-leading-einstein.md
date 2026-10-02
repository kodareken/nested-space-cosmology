# Leading-Einstein diagnostic of the fixed curvature EFT

This is a named leading-order control, separate from the preserved exact
auxiliary-action histories. It measures generated state–geometry feedback after
removing the independent higher-derivative initial data responsible for the
extra pole of the resummed finite truncation. It does not implement the full
first-\(C_W\) correction, certify an EFT range, or declare a completed theory.
The frozen \(C_W\) value remains input history and is not retuned.

Owners are the [canonical evaluator](../src/recursive_horizons/nsc_discovery_leading_einstein.py),
[driver](../scripts/derive_nsc_discovery_leading_einstein.py), and
[focused tests](../tests/test_nsc_discovery_leading_einstein.py).
The only independent geometry fields are \(Q,r,p_Q,p_r\), together with the
unchanged six spinor columns. There is no independent \(\chi,p_\chi\).

## Finite canonical action

Let \(a=-8\pi A,\ b=ar,\ Z=3a,\ F=ar^2/2\), and
\(V=-ar^2-2\pi C_F f^2\). The pre-gauge constraints are varied before imposing
\(L=Q,\beta=0\). Their finite SBP realization is
\[
C_g=\frac{p_Qp_r}{2b}-\frac{ZQp_Q^2}{4b^2}
+\frac{Z(Dr)^2}{Q}-QV-2D[(DF)/Q],\qquad
D_g=p_rDr-QDp_Q .
\]
Matter adds the original \(F_L/\Delta x_q\) and \(F_\beta/\Delta x_q\).
The new gravitational Hamiltonian is
\[
H_g=\Delta x_q\sum\left[
\frac{Qp_Qp_r}{2b}-\frac{ZQ^2p_Q^2}{4b^2}
+Z(Dr)^2-Q^2V+2(DF)(DQ)/Q\right].
\]
It differentiates actual sampled \(D(F)\), rather than substituting a
continuum product rule. The geometry rates include
\[
\begin{aligned}
\dot Q&=Qp_r/(2b)-ZQ^2p_Q/(2b^2),&
\dot r&=Qp_Q/(2b),\\
\dot p_Q&=-p_Qp_r/(2b)+ZQp_Q^2/(2b^2)+2QV
+2D[(DF)/Q]+2(DF)(DQ)/Q^2-(F_Q+F_L)/\Delta x_q,\\
\dot p_r&=aQp_Qp_r/(2b^2)-aZQ^2p_Q^2/(2b^3)
+2ZD^2r+Q^2V_r+2bD[(DQ)/Q].
\end{aligned}
\]
The original representative field generator is unchanged:
\(\dot\Phi=-i(\sigma_2P+\kappa Q\sigma_1)\Phi\), with no multiplicity in
that generator. Source forces retain \(M=4\kappa\) once and the conformal
lapse chain \(F_Q+F_L\). Angular pressure and source cross terms are retained.

Odd real geometry/AP spinor carriers, adjoint pullbacks, and the frozen W
canonical frame are reused. Geometry momenta encode as
\(\pi=\Delta x_gW^\mathsf T p_{\rm nodal}\). The new analytic Jv differentiates
this new SBP rate, including field/source variations. Actual \(Q_{tt},r_{tt}\)
follow its directional action along the actual projected rate. Pure existing
metric contractions then produce normal tides and four-dimensional invariants.
No old auxiliary rate, old six-variable Jv, or old tidal acceleration helper
is called. Prescribed gauge is rejected.

## Principal restriction and finite controls

On a homogeneous frozen background, the two geometry blocks are
\[
A_{\rm kin}=\begin{pmatrix}-ZQ^2/(2b^2)&Q/(2b)\\Q/(2b)&0\end{pmatrix},
\quad
B_{\rm grad}=\begin{pmatrix}0&2b/Q\\2b/Q&2Z\end{pmatrix},
\quad A_{\rm kin}B_{\rm grad}=I .
\]
The leading principal coordinate speeds are one, with no independent Weyl
pole. The implementation combines the conservative quadrature principal
wavenumber, retained field wavenumber, local coupled reaction-Jacobian/source
majorants, and coefficient-gradient terms. This is a step restriction, not
a stability or continuum certificate.

Tests compare canonical action gradients, analytic Jv and projected jets,
dense/FFT carriers, the executed Fourier principal block, and a source-free
magnetic product. That product has \(r=1,\ Q=Q_0\operatorname{sech}(Q_0t)\),
\(p_Q=0,\ p_r=2ar\,[-Q_0\tanh(Q_0t)]\); it has no \(p_\chi\).
Only manufactured short steps below T0.01 are used in tests. Gauge and chart
failures are explicit, and historical auxiliary callbacks are blocked in tests.

## Initial data, runner, and scope

The frozen PREPARED-v2 uniform and coherent fields supply Q and six source
columns. Removed variables and all initial momenta must vanish. At nf256,
the same frozen AP columns and finite Q/r interpolants are extended without
reselecting the eigenframe or occupations. Each resolution prepares its own
source-consistent leading zero-momentum lapse constraint using the equivalent
general-Q C0 radius machinery; no old full-action preparation rate is used.

The six coupled cases are uniform and coherent at nf128, and coherent at
nf256, each with caps 0.001 and 0.0005 and stations T0.3/1/3. Two matched
nf128 coherent frozen-geometry arms use those same caps, for eight cases total.
They share the leading branch's own prepared Q, radius, source and constraint
solution with their coupled counterparts; no old exact-action T0.3 control is
substituted. A scoped adapter
reuses the existing episode runner, process pool, CPU admission ledger, cadence,
event retention, normal-clock trapezoids, immutable chunks and checksum reader.
There is one pool, at most four one-thread workers, a shared budget at most600
CPU seconds, and 64-MiB chunks. No second campaign engine is introduced.
Chart exits retain the last admissible state and never clamp the metric.

Observations report proper length, areal radius, child probability, contrast,
actual tides, normal clocks, new total/gravitational/field energy, instantaneous
matter work and boundary terms, Gram/CAR, and new lapse/shift residuals including
represented and unresolved scales. Space/time comparisons remain finite
indicators, not a holding claim or an error certificate. A small residual in
Hamiltonian units is not promoted to a radius error.

The frozen controller zeros all four geometry/momentum rates while evolving
the same field Hamiltonian on its initial metric. Its actual metric time jets
are zero, proper lengths and normal-clock rates stay constant, and coordinate
metric work vanishes. The new Jv differentiates this controller, rather than
assigning the coupled geometry acceleration afterward. Observations retain
the inactive leading geometry forcing and the evolving lapse/shift residuals.
The initial constraint solution is recorded; preservation of the coupled
geometry constraints under this external control is not claimed.

Coupled/frozen comparisons use matched coordinate times and retain separate
checkpoint-carried normal clocks. These are coordinate-controlled readouts,
not equal-proper-time measurements. A common-proper-time comparison would
require clock alignment and is not inferred from the coordinate tables.

Default invocation is a pure preflight. Preparation is explicit, requires a
frozen producer commit, and creates a new directory only. Run and read-only
check are separate. From the repository root, after source freezing:

~~~sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_einstein.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_einstein.py --prepare --producer-commit HEAD
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_einstein.py --run --workers 4 --cpu-budget 600
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_leading_einstein.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_leading_einstein.py -q
~~~

The new output owner is lab/results/development/nsc-discovery-leading-einstein-v1.
Historical exact-action files and evidence stay preserved. Generated leading
feedback is the scientific comparison; this document records no research run.
