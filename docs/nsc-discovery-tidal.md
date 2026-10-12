# Actual metric tidal consumer

The owner is
[`nsc_discovery_tidal.py`](../src/recursive_horizons/nsc_discovery_tidal.py),
with the creation/check driver
[`derive_nsc_discovery_tidal.py`](../scripts/derive_nsc_discovery_tidal.py)
and eight focused tests in
[`test_nsc_discovery_tidal.py`](../tests/test_nsc_discovery_tidal.py).
It reads authenticated saved discovery stations at coordinate times
0.3, 1 and 3. It takes no evolution step and writes no checkpoint.
The finite, positive chart at these stations is its measured domain;
neither a singularity nor an infinite-time conclusion follows.

## Metric and convention

The realized metric has signature \(+---\):

\[
g=N^2(dt^2-dx^2)-r^2d\Omega^2,\qquad N=rQ.
\]

The convention is
\(R^a{}_{bcd}=\partial_c\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{cb}
+\Gamma^a{}_{ce}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{cb}\),
with \(R_{bd}=R^a{}_{bad}\).
Put \(\sigma=\log N\), and define the warped base pieces

\[
\begin{aligned}
A&=-(\sigma_{tt}-\sigma_{xx})/N^2,\\
B_{00}&=(r_{tt}-\sigma_t r_t-\sigma_x r_x)/(rN^2),\\
B_{11}&=(r_{xx}-\sigma_t r_t-\sigma_x r_x)/(rN^2),\\
B_{01}&=(r_{tx}-\sigma_x r_t-\sigma_t r_x)/(rN^2),\\
S&=\{1+(r_t^2-r_x^2)/N^2\}/r^2.
\end{aligned}
\]

Contracting these pieces gives

\[
\begin{aligned}
R_4&=2A-4(B_{00}-B_{11})-2S,\\
R_{ab}R^{ab}&=(A-2B_{00})^2+(-A-2B_{11})^2-8B_{01}^2
 +2(S+B_{00}-B_{11})^2,\\
K&=4A^2+8(B_{00}^2+B_{11}^2-2B_{01}^2)+4S^2,\\
C^2&=\tfrac43(A+B_{00}-B_{11}-S)^2
 =\frac{(R_h-2)^2}{3r^4},\\
R_4&=\frac{R_h-2}{r^2}-\frac{6(r_{tt}-r_{xx})}{r^3Q^2},\\
R_h&=\frac2{Q^2}\left[\frac{Q_{xx}}Q-\frac{Q_x^2}{Q^2}
-\frac{Q_{tt}}Q+\frac{Q_t^2}{Q^2}\right].
\end{aligned}
\]

The module retains the two scalar routes, the factorized Weyl route and
the cancellation-sensitive identity \(C^2=K-2R_{ab}R^{ab}+R_4^2/3\)
as comparisons. Spatial derivatives use the owned Fourier derivative
with zero Nyquist symbol. The spatial derivative of the product \(N=rQ\)
is compared with its product-rule expansion, retaining their sampled gap.

## Observer and actual time jets

For the normal orthonormal frame
\(e_0=N^{-1}\partial_t\), \(e_1=N^{-1}\partial_x\), and the
usual radius-normalized angular vectors, the covariant tidal tensor is
\(R_{0i0j}=\operatorname{diag}(-A,B_{00},B_{00})\).
These are the instantaneous geodesic-deviation components for freefall
momentarily comoving with that normal frame. A fixed-\(x\) normal worldline
is generally accelerated: its signed orthonormal radial acceleration is
\(a^{\hat1}=N_x/N^2\). Relative acceleration of such supported worldlines
also contains the acceleration field; the normal tidal tensor alone does
not represent that supported motion.

The actual projected Hamilton rate supplies all first time derivatives.
Writing \(P\) for pull followed by prolongation to the fine grid,

\[
\begin{aligned}
Q_t&=P[Qp_\chi/(2F_\chi)],\\
Q_{tt}&=P[(Q_t p_\chi+Q(p_\chi)_t)/(2F_\chi)],\\
\pi&=p_r-F_r p_\chi/F_\chi,\qquad F_r=-8\pi_{\rm math}A_{\rm action}r,\\
r_t&=P[\pi/(2Z)],\\
r_{tt}&=P[((p_r)_t-\(F_r\)_t p_\chi/F_\chi
 -F_r(p_\chi)_t/F_\chi)/(2Z)].
\end{aligned}
\]

\(F_\chi\) and \(Z\) are constant on this chart. The action coefficient
\(A_{\rm action}\) in \(F_r\) is distinct from the curvature piece \(A\).
Both projections remain present. The full actual rate includes the
state-dependent forces; these formulas are the geometry components of
its analytic directional Jacobian along that rate. The tests compare
both accelerations with the independently implemented full-state analytic
Jacobian-vector product. The frozen controller has zero geometry time
jets while its field columns continue to flow. Auxiliary χ is compared
with the metric curvature afterward and is never substituted into it.

## Controls, refinement and provenance

Local analytic jets check spherical Minkowski \(r=x,Q=1/x\) and Milne
\(r=e^t\sinh x,Q=1/\sinh x\): all four-dimensional curvature and tidal
components vanish. The Bertotti–Robinson product with constant \(r\) and
\(Q=Q_0\operatorname{sech}(Q_0t)\) gives \(R_h=2,R_4=C^2=0\),
\(R_{ab}R^{ab}=4/r^4,K=8/r^4\), radial tide \(-1/r^2\), and zero angular
tide. Static sections and mixed Hessian contractions provide additional
sign checks.

The saved nf=256 half-step station at \(T=3,x=2\) gives radial tide about
\(-524203.28\), angular tide \(5632447.06\), Ricci square \(2.545895\times10^{14}\)
and Kretschmann \(5.091790\times10^{14}\). These values remain finite.
The report compares nf=128→256 and, when selected, nf=256→512,
the two saved timestep caps, fixed \(x=2\), each profile's angular-tide
peak, and stored normal proper clocks. Shape comparisons restrict the
fine profile to the shared periodic coarse nodes and report normalized
maximum and L2 differences. Refinement indicators are reported per
observable; ill-conditioned Weyl and scalar cancellations do not impose
a veto or a new acceptance tolerance on actual four-dimensional tides.
Conditioning factors retain the sum of the absolute contracted terms,
including a separate pointwise angular-tide and scalar condition at \(x=2\).

The baseline run envelope binds producing commit `1c9e770`; the nf=512
confirmation envelope binds `528efc7`. Their physics hashes are verified
against the live numerical owners before consumption. Current reader and
consumer SHA256 bindings remain separate from those saved producer
bindings. Checkpoint records and per-array hashes are authenticated by
the existing loader. Before/after directory hashes prove that this
consumer preserved its saved inputs. The nf=512 binding is required only
when selected. One numerical thread is used for the saved-data analysis.

From the repository root, after freezing the consumer source, create a new
record at the requested destination and check it without recomputation:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_tidal.py --write /tmp/nsc-tidal.json
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_tidal.py --check /tmp/nsc-tidal.json
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_tidal.py -q
```

nf=512 is selected by default; `--no-nf512` omits it. `--episode-dir` and
`--confirmation-dir` select other authenticated saved input directories.
Creation refuses an existing destination, including at file opening.
Checking only reads the record and its bound files. The library flag
`production_record_written=false` describes the measurement library;
an explicit CLI creation writes the requested consumer record. Focused
validation completed in about eight seconds, and the measured consumer
CPU time stayed below the assigned thirty-second budget.

## T8 crossing assessment

The successor consumer is
[assess_nsc_discovery_crossing.py](../scripts/assess_nsc_discovery_crossing.py),
with focused [authentication/replay tests](../tests/test_nsc_discovery_crossing.py).
It reuses this note's existing numerical owner without changing its v1 source
pins. Creation writes only the new crossing-assessment-v1 JSON owner. Checking
authenticates its repository-relative source/input hashes and recomputes actual
metric jets, with no evolution step, solve, or checkpoint reset. Both operations
retain a thirty CPU-second budget and one numerical thread.

The crossing campaign's observed-run-binding names producing commit
b7f0dab8de62d9a3472ac4e7b33351ab65d25847. It is a **post-run authenticated
coordinator observation**, not a pre-run signed attestation. The consumer checks
its artifact and producer hashes, separately verifies the original 1c9e770
physics envelope, and authenticates the exact T3 predecessor payloads and complete
per-array hash dictionaries. State, momenta, clocks, and basis are carried into
the crossing without a geometry or memory reset.

At T8, the available nf128→256 comparisons give maximum profile differences
of about 0.30% for the radial tide, 0.51% for the angular tide, and 0.99% for
Ricci square and Kretschmann. The two saved timestep caps change the tidal
profiles by less than \(4.4\times10^{-8}\) relative. At \(x=2\), nf256 half-step
gives radial tide \(-9.665013073\times10^{22}\) and angular tide
\(3.067404935\times10^{24}\). The angular component is well conditioned;
the radial condition is about 2111 and its tested refinement agreement survives
that moderate cancellation. Independent full analytic Jv accelerations agree
with the existing metric-jet consumer. These are measurements of the finite
discrete realization; they do not supply a continuum certificate or a
singularity claim.

The chart remains positive, with \(r\) about 61.9–62.3 and \(Q\) about
\(6.6\)–\(7.6\times10^{-15}\). Normal clock rates are about
\(4.1\)–\(4.7\times10^{-13}\). The occupied covariance spectrum remains within
its CAR range, with the nf256 half-step Gram gap about \(2.6\times10^{-10}\).
Returning child probability about 0.999687 coexists with child proper length
about \(9.14\times10^{-13}\); it does not demonstrate maintained physical size.
Frozen geometry also returns, to child probability about 0.994718 at proper
length 2.3752. The warranted label is **carrier-period return candidate**.
Ambient extent and an exterior echo are not established. Regional probability
and normal-energy stocks remain measurements, not a closed transport ledger.

\(R_4,R_h\), and Weyl remain unresolved cancellation diagnostics. They are kept
separate from the actual normal tidal agreement. The first retained
auxiliary-gap/\(R_h\)-maximum ratios above illustrative 0.1%, 1%, and 10% occur
at T2.40/2.70/2.95 for nf128 and T3.20/3.60/3.95 for nf256. Those markers describe
loss of shell agreement; they are not scientific acceptance thresholds or
physical curvature breakpoints.

Observation joins use timestamps, never row positions. The half-step and frozen
crossing ledgers retain cadence rows through T5.8 and then the T8 endpoint,
with no intervening cadence rows. The cap .001 coupled ledgers retain the full
cadence. The assessment reports these gaps and the unmatched timestamps
explicitly. A numerical replay tolerance applies to the stable endpoint
quantities; replay of cancellation-sensitive scalars is not promoted to a
physical certificate.

From the repository root, after freezing the new consumer:

~~~sh
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_crossing.py --write
.venv/validation/bin/python scripts/lab.py scripts/assess_nsc_discovery_crossing.py --check
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_crossing.py -q
~~~

The immutable output owner is
lab/results/development/nsc-discovery-crossing-assessment-v1.json.
Default invocation previews without writing. The existing tidal-v1 record and
its numerical source bindings remain preserved.
