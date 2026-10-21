# From inherited perturbations to a three-clock prediction

The [parent perturbation calculation](nsc-atlas-parent-perturbations.md)
left two independent amplitudes: initial density and scalar phase. The
present connection determines what freely falling clocks receive from
that same geometry and eliminates both amplitudes from a third response.
It is registered as NSC-P003 in the [prediction register](prediction-register.md).

For the fixed compact-profile family and operational protocol below,

\[
y_3=-2.88990465\,y_1+2.94232704\,y_2,
\qquad
y_i=\log\frac{(\nu_o/\nu_e)_i}{(a_e/a_o)_i}.
\]

The relation is at retained linear order. The amplitudes need not be known,
but the shape, background, initial source positions and motion rule must
be specified. No astronomical observations have been used.

## Why the next step is an observer relation

A new quantum-vacuum choice or a formation history would add a physical
preparation requirement that the present action does not yet supply.
Instead of inventing that state, use operational measurements to constrain
it and test another consequence. This reuses the two existing compact
perturbation columns and adds test emitters and receiver event definitions.
There is no new coupling, incoming waveform or cosmological fit.

Relativistic redshift calculations must retain source velocities, endpoint
potentials, propagation and the observer lapse. See
[Scaccabarozzi, Yoo and Biern, Section 2](https://arxiv.org/html/1807.09796)
for the general importance of observer terms. Their flat-background
survey model is not copied into this closed parent. Below, radial geodesics
derive the clock transfer directly in the actual curved background.

The calculation uses three fixed emitter separations, one receiver age,
two spectral/tolerance comparisons and independent finite metric
variations, within 120 CPU seconds. A rank defect, inconsistent geodesic
variation or source outside the retained light cone stops the inference.
Only compact responses and producer hashes are saved.

## The existing finite dynamics does not select a phase

Substituting the metric momentum equation into the existing linear
system gives a three-dimensional flow for z=(X,Y,Phi):

\[
\dot z=\begin{pmatrix}
0&1&0\\
-\omega^2-\lambda/a^2+2\kappa\Pi^2&-3H&-4H\Pi-2\omega^2\chi\\
\kappa\Pi/2&0&-H
\end{pmatrix}z.
\]

Its trace is -4H. Liouville's identity therefore gives the fundamental
matrix determinant

\[
\det U(t,t_0)=\exp[-4N(t)]>0
\]

at every finite expansion N. At one e-fold it is .018315638889.
Since C evolves as exp(-3N)C0, this invertible flow maps the initial
two-dimensional constraint space bijectively onto the later one.
The constraint covector never vanishes for n>=2: nonzero Pi gives a
nonzero Y coefficient, and at Pi=0 its Phi coefficient is (lambda-3)/a^2>0.
Distinct allowed states cannot become exactly the same full state at
finite time. Covariances likewise propagate by U Sigma0 U-transpose;
the dynamics does not select a unique Sigma0 from unrestricted inputs.

This is a finite deterministic statement about the retained protected
system. It does not exclude asymptotic attraction, tracing out an open
subsystem, a regularity boundary condition or a separately derived
formation mechanism. Those would be additional physical information.

## Fix physical events before calculating redshift

Keep the same action, initial parent and compact profile from the preceding
entry. Write the perturbation as alpha times its density column plus beta
times its phase column. The phase column has zero initial metric/density
and initial field velocity twice omega times the density column's initial
field displacement. Both are supported away from the material walls.

Set L_ref=a0=6.626077076100. Three ideal radial test clocks start at rest
relative to the initial slice normal, at proper separations
(.04,.08,.12)L_ref from the regular parent pole. Their own rest-frequency
emission is continuously available. They follow timelike geodesics and
do not backreact. The central receiver is freely falling by symmetry.
It samples all three signals at a fixed proper age
tau_o=2.251557477076, the background age at N=.8. This is a supplied
dimensionless clock protocol, not an SI ruler or a selected galaxy sample.

Use Newtonian-gauge metric

\[
ds^2=-(1+2\Phi)dt^2+a^2(1-2\Phi)
       (d\theta^2+\sin^2\theta\,d\Omega^2).
\]

For background source angle s, initial proper separation is fixed by

\[
\delta\theta_s(0)=\int_0^s\Phi(0,\theta)\,d\theta.
\]

Holding the source coordinate instead would change its initial physical
distance. Its first-order peculiar velocity and displacement obey

\[
\dot v_s+Hv_s=-\frac{\partial_\theta\Phi(t,s)}a,
\qquad \delta\dot\theta_s=v_s/a,\qquad v_s(0)=0.
\]

The receiver coordinate time shifts because its proper age is held fixed:

\[
\delta t_o=-\int_0^{t_o}\Phi(t,0)\,dt.
\]

The background emission time solves eta_o-eta_e=s. On its incoming
reference ray, theta(t)=eta_o-eta(t). The perturbed null travel equation
then fixes the emission-time shift:

\[
\delta t_e=a_e\left[\frac{\delta t_o}{a_o}
       +2\int_{t_e}^{t_o}\frac{\Phi(t,\theta(t))}{a(t)}dt
       -\delta\theta_s(t_e)\right].
\]

These event shifts are part of the measurement, not optional corrections.

## The derived received-frequency response

To first order, with the source's outward peculiar velocity taken positive,

\[
y=H_e\delta t_e-H_o\delta t_o
  +\Phi_e-\Phi_o-v_e
  +2\int_{t_e}^{t_o}\partial_t\Phi(t,\theta(t))dt.
\]

The integral contains the partial time derivative at fixed spatial
coordinate, not the total derivative along the photon. Field values and
first derivatives come from the evolved solution and its equations,
not differenced stored waveforms. The outward source Doppler sign is
negative because it recedes from the central receiver.

| Initial separation / L_ref | Background emission time | Background frequency ratio | Density response | Phase response |
| --- | ---: | ---: | ---: | ---: |
| .04 | 1.711177664367 | .840672581087 | .017861103048 | -.000731094396 |
| .08 | 1.251919826145 | .721722111018 | .030644152339 | -.001671366286 |
| .12 | .854677804399 | .628648035135 | .038548233209 | -.002804913120 |

Each response is y_i=R_i1 alpha+R_i2 beta. The rays and test bodies stay
well inside the protected parent. The curvature enters through the actual
background and spatial eigenfunctions; radial light propagation itself
uses the exact radial metric coefficient rather than a flat-parent model.

The receiver-age contribution to y is
(H_e a_e/a_o-H_o)delta t_o. Its density-column magnitude exceeds .01
for the furthest emitter. Dropping it would change the chosen physical
events and materially alter the result.

## Independent geodesic variation

Prescribe a finite metric deformation with
A=1+2Phi and B=a^2(1-2Phi). This finite metric is only a device to verify
the first variation; it is not claimed to solve the nonlinear field equations.
Its radial null Hamilton equations give

\[
\frac{d\theta}{dt}=-\sqrt{A/B},\qquad
\frac{d\log(-p_t)}{dt}=\frac{\partial_tA}{2A}
                        -\frac{\partial_tB}{2B}.
\]

For an exactly integrated source velocity v relative to the cosmological
orthonormal frame, the frequency ratio is

\[
\log\frac{\nu_o}{\nu_e}
=-\log(-p_{t,e})+\frac12\log\frac{A_e}{A_o}
 -\operatorname{atanh}v_e,
\]

with -p_t normalized to one at the observer. The code independently solves
the source timelike geodesic, initial proper-distance integral, receiver
proper-age integral and source/photon emission intersection.

The zero-deformation frequency agrees with a_e/a_o to 1.54e-13 in its
logarithm. Along amplitude direction (1,-.4), the analytic coefficient is
.0396701984571. Symmetric finite variations with steps .001 and .0005
differ from it by 2.10e-10 and 5.23e-11. Pure density and phase variations
are also checked, and the mixed-direction check covers all three emitter
positions. These are geodesic-variation checks, not nonlinear
Einstein-scalar validity bounds.

## Calibrate two clocks, predict the third

The first two response rows have determinant -7.44867743642e-6 and
condition number 169.341 in the declared density/phase normalization.
They determine the two amplitudes, after which the third row yields

\[
y_3+2.88990465\,y_1-2.94232704\,y_2=0.
\]

Both amplitudes cancel exactly in the symbolic matrix identity. The
numeric coefficients are conditional on the fixed background, shape and
clock protocol. Refining cutoff 128 to 256 and tightening time tolerance
changes response coefficients by at most 1.85e-12 and null weights by
1.35e-9. The retained null-weight comparison tolerance is 1e-7.

A synthetic check uses (alpha,beta)=(1e-5,-4e-6), giving frequency-log
departures (1.81535408,3.13126989,3.96701985)e-7. The first two recover
the supplied amplitudes and predict the withheld third value. These
generated values are not independent observations or fitted couplings.

For hypothetical independent equal measurement errors sigma_y, the null
residual uncertainty is 4.243682044 sigma_y. At sigma_y=1e-10 the inferred
amplitude standard deviations are 2.4491e-8 and 4.7618e-7. Changing only
the third synthetic value by 2.5e-9 makes a 5.89-sigma null residual under
that assumed noise model. No instrument accuracy or actual detection is
claimed. Correlated calibration errors require the corresponding full
covariance, and nonlinear corrections remain outside the linear claim.

![Three-clock geometry and withheld synthetic prediction](figures/nsc-atlas-parent-clocks.png)

## Reproduction and claim boundary

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_parent_clocks.py
python scripts/run_python.py scripts/plot_nsc_atlas_parent_clocks.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_parent_clocks.py
python scripts/run_python.py scripts/check_prediction_register.py
```

The [compact record](../results/summaries/nsc-atlas-parent-clocks-v1.json)
contains the two response columns, all event shifts, independent geodesic
comparisons, symbolic identities and synthetic noise calculation.
NSC-P003 is retrospective and conditional, with no observational validation.

Two independent Grok reviews confirmed finite state-flow invertibility and
the operational redshift derivation, including initial proper distances,
receiver age, source velocity and integrated-potential signs. The response
matrix's nonzero rank is calculated here rather than assumed from that proof.

The advance is an operational observable relation that survives the two
unselected amplitudes. An arbitrary spatial profile, unknown initial source
velocity or extra source component adds degrees of freedom and need not
obey this relation. Identical local initial data in an unnested model also
give the same protected-domain result. A physical formation/preparation
mechanism and observer tests of the assumed profile remain necessary for
a cosmological identification.
