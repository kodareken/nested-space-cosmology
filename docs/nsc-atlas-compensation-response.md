# From inherited supply to a spatial compensation signature

## Decision and calculation domain

The completed finite core perturbation test does not require another grid.
The external quadrupole amplitude remains unresolved in its own owner.
Three next routes were compared: a complete exterior memory boundary,
formation from an inherited reservoir, and a spatial signature of the
already constrained core plus reservoir. The third has the shortest
verified connection to a local-universe correlation: its source, support,
mass law and inherited vacuum already exist.

Active: derive the compensated spatial kernel and its initial envelope
error from `ScaleHierarchy`, calculate its turnover and vacuum/mass relation,
then bind the result to the atlas. Next: physical population/preparation and
time-dependent observables; an initial kernel alone does not supply them.

Reuse spherical mass constraints, Fourier multipole expansion and the
mass-conserving halo-model idea. Schmidt
[1511.02231v2](https://arxiv.org/html/1511.02231v2), Introduction and Section II,
distinguishes conserved stochastic power from the full matter spectrum.
Chen and Afshordi [1912.04872v5](https://arxiv.org/html/1912.04872v5),
Section III, equations (8)--(11), separate the large-scale field from
compensated profiles. Their fitted profile is not imported here.
The k-squared amplitude and k-fourth power are established mathematics.
Our contribution is their calculated kernel, scale and controlled source
dictionary for the existing inherited geometry, without fitting a new halo.

Keep the existing omega=1e9, zeta=.001, half-depletion initial preparation,
Pi=0 and unchanged inherited action. Compare quadrature orders and tighter
constraint tolerances; quarter/three-quarter depletion controls are distinct
previously declared preparations. At most 120 CPU seconds for the producer,
four useful independent workers, no new evolution or grid campaign.
Independent mass integration, integration by parts, an analytic compensated
shell and the positive gradient-energy remainder test the calculation.
Stop if the half-response root or error comparison is unresolved.

Likely failures are confusing geometric with proper mass, promoting a local
spatial transform to the closed parent's global spectrum, and erasing the
envelope's small uncompensated residual. Each receives an explicit equation
or retained diagnostic below. Population correlations, formation and
late-time lensing are not silently supplied by this initial calculation.

## The same geometry supplies the spatial kernel

Write M(R) for geometric excess mass above the inherited excited parent,
Mc=M(Rc), and Rd for the solved compensation radius. Pi=0 on this slice.
The already verified Einstein constraint is

\[
M_R=\frac{\kappa R^2}{2}\left[U-U_0+\frac{B F^2}{2}\right],
\qquad M(0)=M(R_d)=0.
\]

Here F is the radial scalar derivative and B is the positive radial chart
coefficient. This is an areal geometric mass measure. It is not the
proper-volume energy integral, which has a different volume factor.
Define s=R/Rd, K=kRd and the dimensionless spatial window

\[
W(K)=\frac{1}{M_c}\int_0^{R_d}j_0(kR)\,dM(R),\qquad
j_0(x)=\frac{\sin x}{x}.
\]

This transform is exact as a descriptor of the given areal mass profile.
Its local flat-space Fourier interpretation uses scales much smaller than
the parent's curvature radius. It is not a continuous k spectrum of the
whole compact parent. For that interpretation the overlap is
1/Lparent much less than k, and kRd much less than one; the mathematical
K-to-zero expansion does not erase the global topology.

An independent evaluation uses the solved mass primitive, without
differentiating stored data:

\[
W(K)=j_0(K)\frac{M(R_d)}{M_c}
+K\int_0^1\frac{M(sR_d)}{M_c}j_1(Ks)\,ds.
\]

Cancellation of the monopole and spherical cancellation of the dipole give

\[
W(K)=a_2K^2+\mathcal R_4(K),\qquad
a_2=-\frac{1}{6M_c}\int s^2\,dM,
\]
\[
|\mathcal R_4(K)|\leq
\frac{K^4}{120M_c}\int s^4|dM|.
\]

**Proof.** Use j0(x)=integral from zero to one of cos(xt) dt and Taylor's
remainder |cos y-1+y^2/2| <= y^4/24, then integrate over t and the signed
mass measure. Compact support makes its fourth absolute moment finite.
This is an analytic bound on the defined transform; numerical quadrature
error in its coefficients remains separately measured. In the retained
preparation a2=.0608058925116 and the fourth-order bound coefficient is
.00129586307998. The sign is set by the negative mass contrast farther
out, not by assigning negative total material energy to the reservoir.

## A uniform bridge to the initial envelope, without retuning its mass

The leading classical envelope uses the same initial field at its coherent
turning point. Its initial mass measure drops the gradient contribution:

\[
dM_{\rm env}=\frac{\kappa R^2}{2}(U-U_0)dR,
\qquad dM-dM_{\rm env}=\frac{\kappa R^2 B F^2}{4}dR\geq0.
\]

Thus, with the **same Mc normalization** in both transforms,

\[
g=\frac{1}{M_c}\int\frac{\kappa R^2BF^2}{4}dR,
\qquad |W(K)-W_{\rm env}(K)|\leq g
\quad\hbox{for every real }K.
\]

This follows from |j0|<=1 and positivity of the omitted measure. For the
retained preparation g=3.31911631e-7. In particular Wenv(0)=-g, whereas
W(0)=0. The envelope is not made artificially compensated. Its residual
becomes comparable to the K-squared term near K=.00233635. The all-K
absolute bound therefore does not promise a relative error bound at
arbitrarily small K. The numerical zero of the exact GR measure is about
3e-12 from independent quadratures, recorded rather than forced to zero.

The envelope normalization is derived from the same scalar turning-point
potential; this comparison is not an instantaneous dust interpretation of
the scalar stress. Scalar pressure, lapse and the null contraction still
matter to full relativistic observations. The inequality is an exact
initial source comparison, not a long-time relativistic error bound.

## A computed mass--vacuum--turnover relation

Let Wc be the transform over 0<=R<=Rc with the same Mc normalization.
Define the resolved crossing in .2<K<8 by

\[
\left[\frac{W(K_{1/2})}{W_c(K_{1/2})}\right]^2=\frac{1}{2}.
\]

This is half of the core-only squared response, not half the total cosmic
power. For the fixed initial member,

\[
K_{1/2}=4.0522389138,\qquad
\frac{\Lambda_v}{M_c k_{1/2}^3}=0.2681493031.
\]

The second identity combines the existing solved supply relation
Lambda Rd^3/Mc=17.8427504271 with the newly computed K. Independent
80/160/320-node quadratures and tighter constraint tolerances change K by
less than 9e-12 and the dimensionless ratio by less than 1.8e-12 (twice
the largest relevant differences). These are integration/constraint
comparison ranges for the fixed numerical core shape, not total physical
uncertainties. The envelope's own half-response crossing is4.0522393354.

At fixed zeta and depletion a, the existing large-frequency proof gives a
point core relative to Rd and a normalized reservoir measure

\[
d\mu_{\rm res}=-\frac{s^2[2ab(s)-a^2b(s)^2]}{J(a)}ds,
\qquad b(s)=1024s^5(1-s)^5.
\]

Hence W tends to 1 minus the j0 transform of that positive normalized
weight. Dominated convergence applies on bounded K intervals. For a=1/2,
the limiting crossing is4.1423205290 and
Lambda/(Mc khalf^3) tends to .2578964084. This follows from the original
supply law and its shape, not a new parameter fit. At finite omega=1e9,
the existing quarter/half/three-quarter depletion preparations instead
give .47250817, .26814930 and .20537800. Inheritance has not selected the
depletion: these different numbers remain distinct physical preparations.

Using the **previously supplied** Lambda=1e-52 m^-2 ruler gives Mc
equivalent to6.322e10 nominal solar masses, Rd=.827653 Mpc,
khalf=4.89606 Mpc^-1 and wavelength 2pi/khalf=1.28331 Mpc. No sky data
selected the field parameters or tested these numbers.

## The projected deficit and what a population would add

Define the spatial density descriptor h(R)=M_R/(4pi R^2). Its line-of-sight
projection is a well-defined spatial statistic,

\[
\Sigma(b)=2\int_b^{R_d}\frac{h(R)R\,dR}{\sqrt{R^2-b^2}},
\qquad M_{\rm cyl}(b)=2\pi\int_0^b u\Sigma(u)du.
\]

Every outer negative-contrast annulus produces a negative projected
contrast outside the core. The retained calculation gives positive
mean-minus-local contrast there, and both vanish outside Rd. Negative
surface contrast is not negative tangential shear. The producer compares
the nonsingular line-of-sight integral with independent spherical-shell
fractions in a cylinder. The planar Fourier transform of this projection
equals the three-dimensional transform at zero line-of-sight wave number
by Fubini's theorem. Thus it carries the same compensation scale.

After a justified averaged weak-field lens dictionary, Sigma and
mean-minus-local Sigma are inputs to convergence and tangential shear.
The exact initial real-scalar optical response cannot be obtained merely
by renaming this spatial projection a lens map. The existing full optical
calculation belongs to its separately specified omega2000 history.

For identical translated profiles in a local weak-field population, write
the number-density fluctuation as Delta n. Then

\[
\delta_{\rm red}(k)=\frac{M_c}{\bar\rho_g}W(kR_d)\Delta n(k),
\quad P_{\rm red}(k)=\left(\frac{M_c}{\bar\rho_g}\right)^2
|W(kR_d)|^2P_n(k).
\]

Here rho_g is the mean geometric-mass density in the same dictionary;
Pn is the spectrum of absolute number-density fluctuations. Independent
centres give Pn=nbar, and any bounded Pn gives O(k^4) redistribution power.
This is **not** the full matter spectrum: retain a large-scale field and
its cross spectrum with Delta n. Their correlation can contribute at
different order. A divergent Pn, different profiles, overlap or correlated
preparation needs the corresponding extension. No population was sampled,
and no abundance or dark-matter fraction follows from this single profile.

This gives an observable route: determine core masses, surrounding deficits
and their spatial turnover in a consistently selected population, then
test their joint relation rather than fitting every radius independently.
Formation, age selection and the light-cone dictionary are the remaining
physical inputs, not reasons to tune this computed kernel to a survey.

## The already evolving parent changes the relation

Postprocessing the five **unchanged** saved spherical envelope fields adds
the missing age dependence without another time step. Every final field
hash is checked against the previous evolution record. Physical radius is
R=A x/(omega zeta), and K uses the original Rd, not a newly fitted radius.
The core-only transform uses the same fixed physical cutoff as before.
At slow time600, A=1.27443638723:

| Quantity | Fine value and numerical comparison range |
| --- | ---: |
| khalf Rd,initial | 3.165445 +/- .000031 |
| Lambda/(Mcore,current khalf^3) | .559473 +/- .000024 |
| K-squared coefficient, original mass/radius units | .0998276 +/- .0000072 |

The ranges use twice the maximum fine/medium, time-control/medium or
outer-control/medium difference. The inherited vacuum is unchanged;
the surrounding field and current core mass evolve. This decisively
excludes using .2681493 as an age-independent coefficient of this history.
The numerical comparison does not certify the long-time full-GR reduction.

The fine finite-domain monopole is -1.19446e-6 of the initial envelope
core mass. Its outer-domain comparison range is1.74105e-6, so a nonzero
asymptotic monopole is not resolved. Each retained domain separately obeys
initial total excess minus its **outer** mass flux, to better than3.5e-9
of the core. Core accretion must not be mistaken for total mass creation.
The evolved turnover is nevertheless resolved under those domain controls.
Neither exact compensation at the original Rd nor a strict evolved
k-to-zero fourth-power law has been imposed.

## Evidence and reproduction

The [compact record](../results/summaries/nsc-atlas-compensation-response-v1.json)
contains source hashes, six initial quadrature/preparation cases, five
authenticated evolution postprocessings and the comparison ranges. About
8 CPU seconds suffice; the earlier evolution is reused, not rerun. Tests
cover an analytic compensated sphere, Taylor remainder, the positive
gradient mass, integration by parts, independent projection, physical
wavenumber scaling, source identity and preparation dependence. Grok
review informed the proper-mass, global-spectrum and total-flux distinctions.

```powershell
# First regenerate retained fields if absent, using the existing owner:
python scripts/run_python.py scripts/check_nsc_atlas_hierarchy_evolution.py --output results/summaries/nsc-atlas-hierarchy-evolution-v1.json
python scripts/run_python.py scripts/check_nsc_atlas_compensation_response.py --output results/summaries/nsc-atlas-compensation-response-v1.json --figure
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_compensation_response.py
```

The preceding evolution command is only needed when its ignored raw fields
are absent. Its established compute budget is separate. The new producer
fails with a specific regeneration instruction if those fields are missing;
it never hides a fresh simulation inside spatial postprocessing.
