# Spherical episode assessment

Consumer for one stored episode from the driver job
`coupled-transfer-event-episode`. It does not evolve the state, change the
gauge, add a force, or reset a field. The action, geometry, observer, and
energy ledger stay with their owners:

| Role | Owner |
|---|---|
| Action | `nsc_spherical_feedback_action` |
| Geometry chart and periodic derivative | `nsc_spherical_coupling` |
| Observer clocks and null frame | `nsc_spherical_null_expansion` |
| Normal-energy ledger | `nsc_regional_energy_exchange` |

Positive chart means finite \(r>0\), \(Q>0\), and \(L>0\), the same exit as
`chart_failure`. Proper lapse and radial metric are \(N=rL\) and \(q=rQ\).

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_episode_assessment.py -q
```

The controls are analytic charts on at most 64 nodes, one Christoffel
contraction, and a header read of the saved episode. They do not integrate
a trajectory.

## Direct scalar

On \(h=L^2\mathrm dt^2-Q^2(\mathrm dx+\beta\mathrm dt)^2\),

\[
K=\frac{\dot Q-\partial_x(\beta Q)}{LQ},\qquad
R_h=\frac{2}{LQ}\partial_x\!\left(\frac{L_x}{Q}\right)
-2\left[\frac{\dot K-\beta K_x}{L}+K^2\right].
\]

\(\dot K\) is built from the stored realized \(\partial_t\dot Q\). A declared
`gauge_hold` sets \(\dot L=\dot\beta=0\) because those controls are not
evolved. A missing \(\dot Q\) is not replaced by zero.

The sign is the owned Ricci contraction

\[
\Gamma^a{}_{bc}=\frac12 g^{ad}(\partial_b g_{cd}+\partial_c g_{bd}-\partial_d g_{bc}),
\]
\[
R_{bd}=\partial_a\Gamma^a{}_{db}-\partial_d\Gamma^a{}_{ab}
+\Gamma^a{}_{ae}\Gamma^e{}_{db}-\Gamma^a{}_{de}\Gamma^e{}_{ab},
\qquad R_h=g^{bd}R_{bd}.
\]

It matches that contraction on a flat chart, on \(Q=\mathrm e^{\alpha t}\)
(\(R_h=-2\alpha^2\)), on a static \(L(x)\), and on a static shift. The same
values match `areal_curvature_scalar`. Substituting the Euler–Lagrange
\(\dot p_\chi\), or setting \(\chi=R_h-2\), is rejected before a curvature is
returned.

\[
C^2_{\mathrm metric}=\frac{(R_h-2)^2}{3r^4},\qquad
C^2_{\mathrm proxy}=\frac{\chi^2}{3r^4}.
\]

The proxy gap is reported and is not a metric check. The saved frame field
`weyl_C2` is the proxy.

## Goal

The goal is a renewed structure or a maintained structure with throughflow.
Both require a positive chart, admissible jets, a complete ledger, a
meaningful transfer episode, localization, resolved throughput, and resolved
pressure-plus-lapse work. A reversal of \(0.1\) and a leader share of \(0.5\)
are not gates.

The driver stores canonical Galerkin columns. The conformal ADM owner
writes the half-density \(u=r\sqrt{q}\,\psi\), and the coupling stores
\(\sqrt{\mathrm dx}\) samples whose flat sum \(\sum|\Phi|^2\) is the Hilbert
norm. That sum is the probability element. The proper radial element is
\(\mathrm ds=q\,\mathrm dx\), so the probability per proper length is
\(|\Phi|^2/q\). Weighting the same samples by another \(q\) would count the
radial metric twice. It is not particle or field content.

Shell energy \(F_L/r\) is the ledger's nodal normal energy. It is a
different measure. A caller may instead declare `positive_shell_energy`
nodal samples as a structure weight. Those samples are still not multiplied
by \(q\), and they are not labeled mode probability.

Localization uses the excess nodal concentration
\(n\sum p_i^2-1\) of the declared weight. The proper-arc phase
\(2\pi s/S\), \(s=\int q\,\mathrm dx\), carries the first circular moment.
When that moment cancels, the second-harmonic axis labels an opposing pair.
Equal clumps therefore stay localized. Basis occupations and global \(Q\) or
\(\chi\) maxima are not this measure.

Throughput uses the ledger windows that carry the leading initial nodal
measure, including every window tied with that leader:

\[
\frac{|\int F\,\mathrm dt|}{|\bar E|\,\tau},\qquad \tau=\int N\,\mathrm dt.
\]

The budget term stored by the driver is the signed flux contribution

\[
\frac{\mathrm dE}{\mathrm dt}=F_{\mathrm budget}+W_{\mathrm pressure}+W_{\mathrm lapse}.
\]

\(F_{\mathrm budget}\) is the window reduction of minus the owned flux
divergence from `proper_balance_terms`. It is not the content rate. The
residual of the identity widens the throughput uncertainty. Comparison
uncertainties are absolute differences against a paired timestep or
refinement, or a supplied nonnegative indicator. That indicator is not an
actual bound. `actual_bound` and `closed_regime` stay null. A location claim
that the packet stayed also requires the location indicator to lie inside the
packet's circular width.

An episode needs at least three times on a declared finite window. Two
samples are a cut of a series and do not meet the goal. A resolved arc or
bridge displacement is renewal. Absence of that displacement, with resolved
throughflow and work, is the maintained reading.

## What the driver stores

`assess` reads one dictionary. Arrays use the last axis as the periodic
coordinate. A length-\(n\) vector is spatial and is repeated in time. Time
dependence is shaped `(n_t, n)`.

- `binding`: the four owners above, `driver_job="coupled-transfer-event-episode"`, a nonempty `source_id`, and `phi_bound=true`.
- `chart`: `period`, `x`, `times`, `L`, `beta`, `Q`, `r`, optional auxiliary `chi`, and either `gauge_hold=true` or stored `L_dot` and `beta_dot`.
- `time_jet`: `origin` is `realized_increment` or `finite_ode_derivative`; `projection_kept=true` with the projection name; `Q_dot`; `Q_dot_rate`; `neighbor.dt` and `neighbor.Q_dot_next`. A realized increment also stores `neighbor.Q_next`. The neighbor quotient must match the stored rate. The helper does not difference coarse frames itself.
- `ledger`: `source_id`, `flux_is_budget_term=true`, `normal_energy`, `flux_budget_term`, `pressure_work`, `lapse_exchange`, `window_edges`, `packet_weight`, and `packet_weight_kind`. The kind is `canonical_half_density` for nodal \(|\Phi|^2\) of \(u=r\sqrt{q}\,\psi\), or `positive_shell_energy` for declared nonnegative \(F_L/r\)-class samples. Optional `phi_basis_occupations` must not have the spatial length.
- `comparison`: a paired episode or `refinement_indicator` with nonnegative `arc_location`, `bridge_location`, `packet_resultant`, `occupation_concentration`, `throughput_ratio`, `work_integral`, and `content`. A supplied `actual_bound` is ignored.

Missing `Q_dot`, `Q_dot_rate`, the neighbor increment, or the source binding
returns `completed=false` and does not emit \(R_h\).

## Saved four-case episode

`assess_saved_episode` reads
`results/development/nsc-regeneration-episode-v1.json` and the matching NPZ.
It does not rewrite them and it does not take another step.
The phase is \(T=0.05\) to \(T=0.085\) on
`nf256`/`nf512` at \(\Delta t\le 0.0005\) and \(0.00025\).

The canonical columns are \(\sqrt{\mathrm dx}\) samples. Summed with the six
weights they total occupation 3 on every frame and both fermion grids.
The quadrature image under the owned \(U_f\) has the same total. That nodal
element is the probability weight. It is not multiplied by \(q\).

`frame_quad_Q` matches \(A_g\) applied to `frame_coarse_Q`. The same map
lifts `frame_Q_dot` and `frame_indicator_Q_dot`. The indicator is the stored
neighbor increment divided by `frame_increment_dt`. `frame_p_chi_dot` is not
substituted into \(R_h\). \(L\) and \(\beta\) are the static calibration on
the quadrature nodes. The window budget is `observer_boundary` plus
`observer_pressure` plus `observer_lapse`, not the cancelling packet plus
reservoir flux.

On the finest run the direct \(R_h\) mean is about 8.146. The auxiliary
\(\chi^2/(3r^4)\) differs from \((R_h-2)^2/(3r^4)\) by as much as 0.028.
The proper-arc mean stays near \(x=1.58\) and drifts by 0.0234, inside a
proper width of 1.02. That drift is the same on all four runs. The
null-expansion event is separate: the bridge's positive arc enters the
packet edge at \(x=3.96875\).

The shell leader remains the window \([0,2]\). Its share moves from
0.60870 to 0.61653 and its content rises by 0.08117. The integrated
observer residual is \(3.7\times 10^{-8}\), under \(10^{-5}\). Proper time
is 0.10673 and the leader clock is 0.06438. The stored packet and reservoir
fluxes meet at magnitude 4.98756; that pair is balanced and is not the
content rate. The Gram gap is about \(7\times 10^{-10}\), inside the
\(10^{-8}\) admissibility margin, and \(r\), \(Q\), and \(L\) stay positive.

The reading is maintained structure with throughflow. It is not a renewal:
the probability did not leave its proper width, and no reversal, cyclic
return, or new matter is required. The largest recorded movement of a
compared change is 0.138 percent of that change, on \(K_{\perp,\min}\).
Those movements are refinement indicators. The propagated initial bound is
null, so `closed_regime` stays false. This consumer does not re-audit the
source constraints.

## Saved feedback frames

`inspect_saved_episode` reads headers of
`results/development/nsc-spherical-feedback-episode-v1.npz` only. On
`nf512_dt_0_0005` the frames are 11 by 2048 and contain `r`, `Q`, `chi`,
`proper`, `K_r`, `K_perp`, `rho`, `current`, and `weyl_C2`. The window ledger
is 101 by 4 and includes normal energy, proper flux, proper work, and lapse
work. `final_rate_Q_dot` and `final_p_chi` are one coarse slice of length
511, not the frame nodes, and there is no neighbor increment of that rate.

The metric assessment is not completed. The missing primitives are frame
\(\dot Q\), the realized rate of \(\dot Q\), the neighbor increment, frame
\(p_\chi\), frame \(L\) and \(\beta\), `projection_kept`, and a packet-arc
series. Frame-to-frame differences are not used as those rates.

## Not claimed

No closed regime, no constraint-consistent solution, and no continuum bound.
Renewal of the saved \(T=0.05\) episode is not decided here, because the
realized jet is absent. The driver that stores the schema above is the
campaign; this note is the consumption interface.
