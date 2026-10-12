# PG-normal observer clock

This note corrects a measurement label. It does not exclude a geometry, change the six-case family, or add a publication barrier. The field and metric owners stay as they are. The new adapter is `nsc_rn_observer_clock`. It reads a supplied metric and returns a clock. It does not advance φ.

The chart, shared with the RN source and the flat-slice producer, is

\[
ds^2=N^2dt^2-(dr+\beta\,dt)^2-r^2d\Omega^2.
\]

In the coordinate basis, \(g_{tt}=N^2-\beta^2\), \(g_{tr}=-\beta\) and \(g_{rr}=-1\). The unit normal to a \(t=\)constant slice is

\[
u=N^{-1}(\partial_t-\beta\partial_r),
\]

so its worldline and proper time are

\[
\frac{dr}{dt}=-\beta(t,r_{\mathrm{obs}}),\qquad \frac{d\tau}{dt}=N(t,r_{\mathrm{obs}}).
\]

A metric variation \((\delta N,\delta\beta)\) at fixed radius, together with a displacement \(\delta r\), has the tangent

\[
\frac{d(\delta r)}{dt}=-\delta\beta-\beta_r\,\delta r,\qquad
\frac{d(\delta\tau)}{dt}=\delta N+N_r\,\delta r.
\]

\(\beta_r\) and \(N_r\) are derivatives of the same interpolant that supplies \(N\) and \(\beta\).

## What the saved scalars are

`nsc_rn_pg.stage_rates` stores two passive numbers. The updates are the stage increment and the RK4 combination at the cited lines (1002–1005, 1052, 1072–1075, 1119–1128). Neither number is read by the Dirac rate or by the metric reconstruction.

The killing rate is \(\sqrt{N_{\mathrm{out}}^2-\beta_{\mathrm{out}}^2}\) at the fixed outer node when that \(g_{tt}>0\). That is the static proper-time rate of the outer node. When \(g_{tt}\le 0\) the producer stores 0. The zero is a clamp, not a record that the curve stopped being timelike.

The normal rate is \(N[0]\), the lapse at the fixed excision node. No observer radius enters. The saved normal clock is the integral of that lapse along the curve of constant areal radius at excision. On the vacuum RN chart this curve is spacelike: \(g_{tt}(r_{\mathrm{excision}})<0\), while \(N>0\). The integral is therefore not proper time along the PG normal, and it is not a child clock. The adapter leaves every stored value unchanged and only attaches this meaning.

## Transported clock

The primary normal observer starts at the parent packet centre, \(12 r_m\), at a declared initial coordinate time \(t_0\). Proper time is measured from a declared origin, not from the saved normal scalar. An optional second observer starts at the outer edge of the child interval. On this chart \(\beta>0\), so the normal moves toward smaller \(r\). The step whose stage would pass the excision node is refused, and the clock stops at the last state still inside the closed nodal interval. It is not moved onto the boundary and it is not continued by extrapolation.

A fixed-areal static clock has \(dr/dt=0\) and

\[
\frac{d\tau}{dt}=\sqrt{N^2-\beta^2}
\]

only while \(g_{tt}>0\). A null or spacelike sample, and any sample with \(N\le 0\), returns `static_not_timelike` or `lapse_not_positive`. The step is refused. The last finite \((t,r,\tau)\) is kept. The rate is not replaced by zero.

The radial domain is the closed interval of the supplied nodes. An RK4 stage whose radius lies outside that interval returns `domain_exit`, keeps the last valid state, and does not call for a metric value beyond the nodes. A later step on a stopped observer does not move it.

## Interpolation and the RK4 step

Radial samples use piecewise-linear interpolation on the supplied nodes. On an open piece the slope is that piece's secant. At an interior node the value is the nodal value and the slope used for the tangent is the right secant. The left secant is reported and is not averaged into the rate. A centered difference that straddles the node sees the average of the two secants, which is a different number when the nodal data kink. The left endpoint uses the first interval only. The right endpoint uses the last interval only.

The observer step is classical RK4. It consumes either the four coupled metric frames of the existing field step, in the same order as `rk4_step`, or one evaluator \((t,r)\mapsto(N,\beta,N_r,\beta_r)\). A `stage_rates` mapping supplies lapse and shift; the grid radius is passed once beside those four frames. The frame is sampled at the observer radius of that stage. No new field equation is solved. `initial_observers`, `rk4_normal_step`, `rk4_static_step` and `label_saved_clocks` are the calls a later response adapter can use.

## Saved stations

The six-case archive writes a station every \(0.5\,r_m\). Each station stores the historical normal and killing scalars, the radius, the coupled mass profile, and the frozen initial lapse and shift. It does not store the evolved lapse and shift, and it does not store the four RK4 metric stages inside a march step.

Those historical scalars are recoverable as stored, with the labels above. The proper time of the transported normal observer on the old coupled run is not recoverable: the stages that would determine \(N(t,r_{\mathrm{obs}})\) and \(\beta(t,r_{\mathrm{obs}})\) were not saved. Replaying the frozen initial metric is the clock of that frozen metric, not of the coupled run.

If lapse and shift samples are supplied on a common radial grid, the adapter can integrate an approximate worldline by piecewise-linear interpolation in time between those samples. Comparing stride \(0.5\) with stride \(1\) measures that sampling. It is not an error bar on the missing stages, and the result is not marked as the old run's proper clock.

## Checks

The focused tests cover a flat metric, the vacuum RN relation between the static rate \(\sqrt{f}\) and the normal rate \(N=1\) (Killing energy 1 along the normal), one manufactured polynomial trajectory, domain exit, static spacelike refusal, the tie and edge secants, and the observer tangent against a centered amplitude difference at two amplitudes. Producer files are hashed and left unchanged. Verification uses one numerical thread and stays inside a 60 s process-time budget. No family campaign and no result record is written.
