# A fixed observer measures the forming concentration

## Decision and finite calculation

The [seed formation](nsc-atlas-seed-collapse.md) is complete in its finite
envelope domain. Compare a longer stability run, an instantaneous projected
lens, and a fixed-observer signal using the actual time-dependent source.
The third connects the result outward without installing another mass model.

Complete: the short source interval, weak metric, mass-current derivative
and fixed-observer null problem are calculated below. Image direction,
proper delay and frequency satisfy independently tested arrival-map
identities. The [NSC-P004 continuation](nsc-atlas-dynamics-lensing.md) now
uses freely falling motion to predict a withheld optical response,
keeping model and numerical errors separate.

Keep the formed 10% seed, inherited action, background and source unchanged.
Use a local averaged weak metric with equal leading scalar potentials and
the actual closed-parent spatial curvature. This is a declared metric
approximation, not a newly solved full Einstein-scalar evolution. The
emitter/receiver are held at comoving radius r=.25, at angles0 and pi-.08.
They are test clocks, not assumed to fall freely or move with the scalar
fluid. This sphere is not the previously tracked contracting mass shell
that started at r=.25 in the much earlier preparation.
Their clocks start together on the declared earlier slab slice. Match
emission proper age and receiver worldline in source and homogeneous control.

The saved three grids, time and outer controls remain source inputs.
Short backward/forward propagation uses the same field operators, and
time interpolation uses mass-current derivatives. Additional source-time
and ray-tolerance controls have a120 CPU-second postprocessing budget in
four workers. No full formation grid is rerun; missing raw data have an
explicit regeneration command. Full short-slab fields remain local.

The likely failures are an inaccurate weak-metric dictionary, numerical
cancellation in the small frequency/time differences, and mismatch of
observer endpoints. Keep approximation scope separate from numerical
error, integrate frequency contrast directly, and test the fixed-endpoint
arrival-time/redshift identity. A failed identity is diagnosed before
reporting an observer signal.

## The declared metric and its source dictionary

Retain the omega1e9 source, Lambda=9/32, D=3/8, C=.0227765023676067 and
L=.0004815623798469198. The new calculation does not rerun formation or
replace its field by a fitted halo profile. Every starting field is checked
against the preceding saved hash. Set N=ln A and B0=1-C L^2 r^2. The metric is

\[
ds^2=-(1+2\Phi)dt^2+A^2L^2(1-2\Phi)
\left(\frac{dr^2}{B_0}+r^2d\Omega^2\right),\qquad
\Phi=\frac{D L^2}{A}V.
\]

V is the existing boundary-clock potential, with Delta_r V=(|f|^2-1)/2
and V(Lbox)=0. Equal leading scalar potentials are part of the local
averaged weak-field dictionary. Null propagation is solved in this specified
metric; this does not supply the omitted higher-order Einstein-scalar
metric, anisotropic stress or rapid carrier corrections. The numerical
error estimates below do not enclose all those model errors.

The compact source interval is
ln(1.2)-.0002 <= N <= ln(1.2)+.0003. Symmetric field steps extend backward
and forward from the retained endpoint, using the same radial Poisson
and sine kinetic maps. Reversal returns the field at relative error around
1e-13 or smaller. No source parameter or phase is retuned.

## Mass current determines the potential's time derivative

Let delta=|f|^2-1, j=Im(f* f_r), m(r)=integral_0^r s^2 delta ds and
a=1/(omega L^2 A^2 H), as in the formation calculation. Continuity implies

\[
\delta_N=-\frac{a}{r^2}\partial_r(r^2j),\qquad
m_N=-a r^2j.
\]

The boundary-clock potential can be written as

\[
V(r)=-\frac{1}{2}\int_r^{L_{\rm box}}\frac{m(s)}{s^2}ds.
\]

Consequently its exact reduced-source derivatives are

\[
\Phi_r=\frac{D L^2}{2A}\frac{m(r)}{r^2},\qquad
\Phi_N=-\Phi+\frac{D L^2 a}{2A}\int_r^{L_{\rm box}}j(s)ds,
\]
\[
\Phi_{Nr}=-\Phi_r-\frac{D L^2 a}{2A}j(r).
\]

These identities connect the forming field's mass transport directly to
its photon frequency response. They include the outer contribution and
the boundary clock; no differencing of saved density plots is needed.
The field derivative is evaluated spectrally. Cubic radial interpolation
has exact polynomial mass/current integrals, and time-Hermite interpolation
uses the derived potential slopes. The metric and both derivatives then
belong to one differentiable interpolant. A pure-phase initial field gives
an independent analytic check: with zero density contrast and phase theta,
Phi_N=-(D L^2/A)(a/2)theta when theta vanishes at the boundary.

## Proper clocks and a receiver fixed before solving the ray

The emitter is at r_s=.25, phi=0; the receiver is held at the same comoving
radius and at Theta=pi-.08. Their common proper-clock rate is

\[
\frac{d\tau_s}{dN}=\frac{\sqrt{1+2\Phi(r_s,N)}}{H(N)}.
\]

Both clocks start at zero on the slab's initial slice. The specified
emission age equals the homogeneous clock's elapsed age from that slice
to A=1.2. Solve the emitter's actual N event from that same proper age.
This preserves the source clock and receiver worldline when comparing with
the homogeneous parent. Fixing the coordinate emission event instead would
be a different protocol.

The clocks require the acceleration needed to hold their comoving positions.
The scalar's local radial flow does not automatically give their velocity.
A freely falling clock or one tied to the scalar flow needs its own initial
velocity, trajectory, Doppler term and aberration. These can be larger than
the small gravitational frequency contrast here. No galaxy/star kinematics
or clock-detector microphysics is derived by this ideal test-clock model.

## Null evolution and cancellation-safe frequency

Use conformal path time dq=dt/(A L), with
F=sqrt[(1+2Phi)/(1-2Phi)]. After a constant momentum normalization, write

\[
E=\sqrt{B_0p^2+J^2/r^2},\qquad \mathcal H_q=F E.
\]

Hamilton's equations give r_q=F B0 p/E, phi_q=F J/(E r^2),
p_q=-F_r E-F[(B0)_r p^2-2J^2/r^3]/(2E), and N_q=A L H.
For inward source angle psi, take p_e=-cos(psi)/sqrt(B0(r_s)) and
J=r_s sin(psi). Shoot psi so the first outward recrossing of r_s reaches
the specified Theta. This solves one weak-image branch in the declared
angular bracket; it is not a global uniqueness or multiple-image proof.

At reception, the measured angle relative to the radial axis is
beta_o=atan2[J/r_s,sqrt(B0)p_o]. Normalize the emitted frequency to one.
The exact redshift formula for the specified metric is

\[
\log\frac{\nu_o}{\nu_e}=N_e-N_o+
\frac{1}{2}\log\frac{1+2\Phi_e}{1+2\Phi_o}
+\int\frac{2\Phi_N N_q}{1-4\Phi^2}dq.
\]

Integrating the small last term directly avoids subtracting two almost equal
photon energies. The endpoint Hamiltonian energy gives an independent
identity. The measured proper flight time integrates the receiver clock
rate along the arrival interval; a coordinate Shapiro delay is not that
same observable. The homogeneous curved-parent geodesic provides an exact
control: its conformal length is
2 asin[sqrt(C L^2) r_s sin(Theta/2)]/sqrt(C L^2).

## One arrival map generates the linked measurements

Let T(tau_e,Theta)=tau_o be the arrival map. Variation of a null geodesic
between the specified timelike worldlines gives

\[
\nu_o\,d\tau_o=\nu_e\,d\tau_e+J_{\rm physical}\,d\Theta.
\]

Therefore, with unit emitted frequency,

\[
\frac{\partial T}{\partial\tau_e}=\frac{1}{\nu_o},\qquad
\frac{\partial T}{\partial\Theta}=R_o\sin\beta_o,
\quad R_o=A_oLr_s\sqrt{1-2\Phi_o}.
\]

For the excess delay relative to the homogeneous reference, subtract the
corresponding derivatives. These are general endpoint identities, not a
new NSC-only law. The new result is their evaluated application to the
actual formed source. They show how a single arrival map constrains
direction, delay and frequency. Those readouts cannot be counted as three
independent physical discoveries when the full arrival map is already known.

The spacetime Fermat/endpoint framework is attributed to
[Perlick](https://arxiv.org/abs/1010.3416) and
[Fleury--Larena--Uzan,2011.04440v4](https://arxiv.org/html/2011.04440v4).
Perlick's abstract and the published review's Fermat discussion, and
Fleury et al.'s time-delay/Fermat discussion and Appendix A, were read.
Their specific lens approximations are not substituted for this ray
integration. The endpoint relation is derived directly above and checked
independently by changing emission age and receiver angle.

## Evaluated observer response

The [compact record](../results/summaries/nsc-atlas-formed-observer-v1.json)
contains the three retained source grids, source-time and outer controls,
a finer local source slab and tighter ray integration. Four workers use
about40 CPU seconds for the complete final producer. Including the first
bounded calculation, the turn remains within120 CPU seconds. All starting
field hashes match the preceding formation record.

For the fixed primary receiver:

| Quantity | Fine result and numerical comparison range |
| --- | ---: |
| Image-direction change | (1.027923 +/- .009616)e-7 rad |
| Proper flight-time excess | (1.347564 +/- .000582)e-11 atlas time units |
| Log received-frequency change from the matched homogeneous experiment | (-1.44088 +/- .01399)e-11 |
| Log frequency after removing its own A_e/A_o factor | (-7.61575 +/- .14212)e-12 |

The comparison ranges use twice the largest source-grid, source-time,
outer-domain, local-slab or ray-tolerance difference. They are not certified
full model-error intervals. The reader should use the exact saved values
when reproducing numbers rather than the rounded table.

The proper-age derivative of excess delay is1.44090e-11; its independently
predicted frequency value differs by1.36e-15. The angle derivative is
1.49205e-11 at the smaller angular step, differing from R_o sin(beta_o)
minus the reference by6.59e-16. These checks hold the same source age and
receiver definition; they do not compare photons arriving at different
unmatched observer positions.

Using the previously supplied Lambda=1e-52 m^-2 ruler gives an image shift
21.20 +/- .20 milliarcseconds and a proper delay27.5906 +/- .0120 days.
The angular unit conversion does not itself need an SI ruler; the delay
conversion does. An equivalent frequency-velocity notation is about
-4.32 mm/s, but this is a frequency contrast, not a predicted physical
observer velocity. The actual fractional frequency change is the table's
dimensionless quantity.

The same model also fixes the ruler-independent relation

\[
\Lambda_{\rm physical}(c\,\Delta\tau)^2
=(5.10730\pm.00441)\times10^{-23}.
\]

It connects the inherited vacuum to the local excess flight time under
the stated protocol. It does not select the observed absolute Lambda or
establish that the preparation describes a real lens.

## Why the evolving source matters

Freeze the comoving potential while preserving the same clock/receiver
protocol. The expansion-removed frequency term then vanishes exactly for
equal-radius endpoints. The total log-frequency contrast becomes
-6.79268e-12 instead of -1.44088e-11. This discrepancy is much larger than
the numerical comparison range. The direction and delay change much less.
A static snapshot can therefore conceal the formation dynamics that the
frequency channel measures. The frozen potential is an optical control,
not another self-consistent formation solution.

The figure's two additional receiver offsets illustrate the same finite
map. Only the central protocol has the full quoted source/control error
comparison; the extra points are not independent precision predictions.
An independent Plummer example checks the fixed-endpoint ray angle and
proper delay against their Born limits, and the no-potential curved-parent
case checks the exact background geometry. Pure-phase source, polynomial
mass integrals and spectral derivatives test the metric construction.

Grok independently recovered the weak-image scale and highlighted the
kinematic difference between held clocks and the flowing scalar. Its
frozen-coordinate-flight estimates use a different event/clock protocol;
the proper-clock and evolving-source values here are checked by the two
endpoint identities and the numerical controls rather than adopted from
that estimate.

This is a conditional forward prediction, not an empirical confirmation of
nesting or a unique inverse reconstruction of the exterior. The early seed,
boundary preparation, observer trajectories and dimensional ruler remain
specified. The next observational use must supply measured or physically
prepared emitter/receiver kinematics and assess the omitted metric terms;
it must not reuse these numbers for arbitrary moving astronomical objects.

```powershell
# Regenerate the preceding source first only if its ignored raw fields are absent.
python scripts/run_python.py scripts/check_nsc_atlas_formed_observer.py --output results/summaries/nsc-atlas-formed-observer-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_formed_observer.py
```
