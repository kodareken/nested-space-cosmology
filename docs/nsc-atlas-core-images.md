# One source, one receiver, two paths through the evolving parent

The [optical calculation](nsc-atlas-core-optics.md) maps an emitted ray to
its receiving location. This connection inverts that map for a fixed receiving
worldline. Both images belong to the same emitted pulse and reach the same
observer at different proper times. No new lens profile or scalar source
is fitted.

## Choice and finite calculation

The broader alternatives were a physical population law, nonspherical
persistence, and inversion of the existing optical map. The first still
needs a state-selection mechanism; the second requires a coupled perturbation
system on the time-dependent background. The inverse problem uses the
completed geometry to produce directly comparable image, clock and frequency
quantities. It is the next outward observational connection.

The only new physical datum is a receiver's angular location. Two numerical
angular patches and five nearby source ages resolve the inverse and its
derivative. They add no dynamical state: the unchanged field solver supplies
the geometry, and photons remain passive. Each of four parallel runs saves
its own completed result. An unresolved root, failed optical identity or
source-hash mismatch stops the interpretation locally without changing the
physical preparation. The initial budget is 4500 CPU seconds; fields are
regenerated rather than committed.

Implementation: `nsc_atlas_core_images.py` calls the existing `stream_rays`
and `evolve_arrays`. `check_nsc_atlas_core_images.py` authenticates the final
fields against the already retained evolution, inverts the angular patches,
and compares three resolutions and a separate optical-time control. No
retained evolution or optics source file is changed.

The general exact lens-map framework comes from
[Frittelli and Newman](https://arxiv.org/html/gr-qc/9810017); the
[spacetime lensing review](https://doi.org/10.12942/lrr-2004-9) supplies the
general Fermat, Jacobi and brightness context. Their generic identities are
reused here. The contribution is their explicit connection and calculation
on this inherited, forming scalar-parent solution, not invention of the
general lens equations. Existing real-scalar optical research, including
[Aimar et al.](https://arxiv.org/abs/2512.23800), remains relevant prior art.

## The fixed receiver

The source and receiver are held at areal radius Rs=.6. Their proper clocks
are zero on the initial slice. Put the source at azimuth zero and the receiver
at pi-b, with b=.0002 radians. Emit the central pulse at source proper age2.8.

Let Phi(psi,tau_e) be the unwrapped positive-orientation landing angle from
the old optical calculation. Two branches solve

\[
\Phi(\psi_+,\tau_e)=\pi-b,\qquad
\Phi(\psi_-,\tau_e)=\pi+b.
\]

The second photon actually travels in the negative azimuthal direction.
Its arrival coordinate is -(pi+b), the same receiver modulo 2pi. Spherical
reflection symmetry permits using the positive-orientation ray calculation
for its scalar observables. This identifies the receiver exactly; it does
not compare two nearby but different observers.

Search only the two declared patches near source angles .00652 and .00540.
Both are on decreasing parts of Phi. Existence in these patches follows
from bracketed roots of the calculated map; the finite search does not
exclude an additional central image or other branches.

The received sky angle has magnitude

\[
\sin|\xi_i|=\frac{\sin\psi_i}{\nu_i}.
\]

Choose the positive sky axis opposite increasing azimuth. Then xi_+ is
positive and xi_- negative. Here nu_i is received energy divided by the
common unit emitted energy. This relation follows from conserved angular
momentum, L=Rs sin(psi), in the local observer tetrad. The receiving clocks
and boosts are evaluated at each image's own arrival event.

## An exact connection between image geometry, delay and frequency

For a smooth family of null geodesics, variation of the null action gives
the endpoint identity

\[
p_a\delta x_o^a-p_a\delta x_e^a=0.
\]

At the held emitter, with fixed spatial position and emitted energy one,
the second term is -d tau_e. At the held receiving sphere, the first is
-nu_o d tau_o + L d phi_o. Thus

\[
d\tau_o=\frac{d\tau_e}{\nu_o}
+\frac{L}{\nu_o}\,d\varphi_o.
\]

This is a local branch identity; it requires neither a stationary metric
nor conserved photon energy. Angular momentum is conserved by spherical
symmetry. At fixed source event it predicts

\[
\frac{d\tau_o}{d\psi}
=\frac{R_s\sin\psi}{\nu_o}\frac{d\Phi}{d\psi}.
\]

At fixed receiving worldline each image instead obeys d tau_i/d tau_e=1/nu_i.
For the delay Delta=tau_--tau_+,

\[
\frac{d\Delta}{d\tau_e}=\frac{1}{\nu_-}-\frac{1}{\nu_+},\qquad
\frac{d\Delta}{d\tau_+}=\frac{\nu_+}{\nu_-}-1.
\]

The second relation uses only the receiver's clock and the images' frequency
ratio. Correlating the same source feature across two images removes the
need to know its absolute emitted frequency. This is a geometric consistency
relation; the particular values on its two sides test the selected source
and geometry, not nesting uniquely.

The calculation checks the first identity with independently interpolated
arrival times and landing angles. It checks the fixed-worldline derivative
by emitting at proper ages2.8 plus/minus .0002 and .0001, solving the inverse
again at each age, and comparing Richardson-combined delay differences with
the transported frequencies. Angles, delays and frequencies cannot be
adjusted separately without violating this relation.

## Pulse energy and brightness

For an infinitesimal isotropic emitted pulse energy E_e, the observed
bolometric fluence on a detector perpendicular to the ray is

\[
\mathcal F_i=\frac{E_e}{4\pi}
\frac{\nu_i}{|D_{\perp,i}D_{\parallel,i}|}.
\]

The Jacobi area is the received cross-section per emitted solid angle.
One redshift factor accounts for photon energy. Instantaneous bolometric
flux from a continuous source of luminosity L_e has an additional rate
factor, since d tau_o/d tau_e=1/nu_o:

\[
F_i=\frac{L_e}{4\pi}
\frac{\nu_i^2}{|D_{\perp,i}D_{\parallel,i}|}.
\]

Pulse fluence and continuous-source flux are distinct observables. Their
ratios cancel a common isotropic emission normalization. Anisotropic source
emission, absorption, finite duration, finite source size and wave interference
need their own model before these ratios can be compared with such data.
The signed determinant records parity, separately from positive brightness.

## What can be calibrated and what can be tested

This benchmark specifies a dimensionless geometry, emitting age and receiver
position. Angles, image frequency ratios and fluence ratios require no
choice of metres or seconds. If a dimensional embedding assigns length
L_unit to one atlas unit, a calculated delay delta corresponds to

\[
\Delta t=\frac{L_{\rm unit}}{c}\,\delta.
\]

A measured delay could calibrate L_unit=c Delta t/delta **conditional on this
geometry and preparation**. The angles and ratios would then be separate
tests, not additional quantities to fit by changing that same unit. Scaling
the entire solution also scales dimensional couplings; the calculation has
not selected a scalar mass or established that an observed lens is this
benchmark. A narrow-band common pulse and a sufficiently small, transparent
source are explicit measurement assumptions.

## Reproduction and current state

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_core_images.py --output results/summaries/nsc-atlas-core-images-v1.json
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_core_images.py
```

The inverse implementation and independent flat/closed-FLRW tests are in
place. The four-run calculation is active; no image numbers or successful
convergence are claimed before its checks finish. Numerical interpolation,
near-caustic conditioning and small delay differences are the principal
failure risks. Completed evidence determines whether the two-image result
is resolved or requires a specifically identified correction.
