# Prediction register

This register states what the fixed model predicts, where it applies and how
agreement could be rejected. Each entry distinguishes derived consequences,
numerical sensitivity and observations. The [machine-readable register](../results/prediction-register.json)
binds the evidence bytes and the exact locations of reported values.

## NSC-P003 — Three freely falling clocks in the inherited parent

Registered retrospectively on 21 October 2026, with a synthetic check and
no astronomical observations. The [clock derivation](nsc-atlas-parent-clocks.md)
uses the existing compact scalar perturbation and two unknown density/phase
amplitudes. Initial proper emitter separations are .04, .08 and .12 times
the initial parent curvature radius. The receiver samples at the fixed
proper age corresponding to background N=.8; all clocks freely fall.

For y_i=log[(received/emitted frequency)_i/(background a_e/a_o)_i],
the conditional linear prediction is

\[
y_3=-2.88990465\,y_1+2.94232704\,y_2.
\]

The first two measurements may calibrate only the two state amplitudes;
the third is withheld. The [record](../results/summaries/nsc-atlas-parent-clocks-v1.json)
retains the fixed couplings, profile, clock events, independent geodesic
variation and numerical comparisons. Refinement changes null weights by
1.35e-9; the comparison tolerance is 1e-7. With independent equal clock
errors sigma_y, the null residual has uncertainty 4.243682044 sigma_y.

This is not a preregistration, empirical confirmation, SI scale selection
or a relation for arbitrary incident profiles. A discrepancy beyond the
full measurement and justified model/remainder error would reject the
declared family. Additional source velocities or arbitrary spatial states
must be included as additional inputs rather than hidden retuning.

## NSC-P002 — A charged vacuum transition and its curvature observable

Registered retrospectively on 20 October 2026; no withheld observations.
This is the explicit Einstein/four-form/charged-membrane alternative,
described in the [source and shared-parent derivation](nsc-atlas-flux-transition.md)
and [compact result](../results/summaries/nsc-atlas-flux-transition-v1.json).
It does not change P001's action or unresolved finite-amplitude status.

With the declared parent flux, bare vacuum and shared charge,
Lambda L_ref^2 is 1.125, .5 and .125 in the parent collar, child and
grandchild core. Two signed moving walls and a common cosmological
join realize those regions on 0<=T_middle/L_ref<=.25. The local
free-fall tidal trace equals Lambda in each vacuum region and removes
its mass contribution. Three consecutive equal-charge regions obey
Delta^2 Lambda=kappa e^2; the benchmark second difference is .25/L_ref^2.
The flat dust parent's H^2+2 a_ddot/a must equal the outer vacuum trace.

These are conditional predictions from stated inputs. The physical ruler,
charge, parent flux, bare vacuum and preparation are not inferred from
observations. Independent quadrature/ODE differences and sampled
residuals check the calculation, not experimental uncertainty.
A discrepancy in these relations beyond declared instrument/model
corrections would reject this candidate's tested claim. Establishing
cross-region measurement access, formation and a realistic cosmological
history remains separate.

## NSC-P001 — Perturbative cap and moving wall

**Registered retrospectively on 16 October 2026.** These results had already
been derived and calculated. This entry is neither a preregistration nor a
test against withheld observations. Its status is a **closed conditional
perturbative result**; observational testing and nonlinear finite-amplitude
closure remain open.

The [cap/wall derivation](nsc-atlas-cap-wall-response.md),
[wall evaluation](nsc-atlas-wall-evolution.md) and
[model snapshot](../results/development/nsc-atlas-cap-wall-response-v1/model-specification.json)
fix a canonical real massless conformal scalar, a de Sitter cap, a
Schwarzschild exterior and a pure-tension wall. The chosen inputs are
$h=H\ell=1$, $b=k\ell=2$, $m=M/\ell=0.4$ and $\xi=1/6$; they are not fitted
observations. The source is the original $C^3$ compact characteristic pulse,
with quiet earlier cap and bottom data. The expanding timelike reference
wall has $x\in[0.45,0.55]$ on the finite no-return domain, with the signed
normal directed from cap to exterior.

Write $s=\tau_W/\ell$ and $\alpha=\kappa_4\varepsilon^2/\ell^2$. In the
Jordan frame the predictions are coefficients in

$$
\frac{\mathfrak m_+}{\ell}=m+\alpha\mu_++O(\alpha^2),\qquad
\frac{\mathfrak m_-}{\ell}=\frac{h^2r^3}{2}+\alpha\nu_-+O(\alpha^2),
$$
$$
r=x+\alpha\rho+O(\alpha^2),\qquad
V_W=V_0(1+\alpha\eta)+O(\alpha^2).
$$

The cap vacuum is subtracted at the actual radius $r$. The coefficient
$\eta$ concerns the exterior null embedding; it is not the central
observer's proper clock. The common quiet event fixes the integration
constants. At reference radius $x=0.55$, the fine side of the retained
129/257-field-node comparison gives:

| Dimensionless coefficient | Value at $x=0.55$ | Maximum 129/257 field-grid difference over nine common cuts |
| --- | ---: | ---: |
| $\mu_+$ | -28.88342232555452 | $3.3761187623326805\times10^{-7}$ |
| $\nu_-$ | 17.759078420520662 | $2.088306541025986\times10^{-7}$ |
| $\rho$ | -5.2868390416906355 | $3.6360701116677774\times10^{-8}$ |
| $\eta$ | -0.8178125677354874 | $3.5651870344821646\times10^{-9}$ |

Values are copied from `/solutions/1/wall_samples/8` of the
[wall field-comparison record](../results/development/nsc-atlas-wall-evolution-v1/field/response.json),
with sensitivity values from its `/wall_differences` object. Its
[manifest](../results/development/nsc-atlas-wall-evolution-v1/field/manifest.json)
retains the original source context. Other numerical settings stay fixed
within this comparison. These sampled differences are not error bars, an
endpoint-specific bound or a total error enclosure. No physical amplitude,
SI length scale or finite-amplitude remainder bound is supplied here.

At the same retained order, the
[exact cap/wall record](../results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json)
and [manifest](../results/development/nsc-atlas-cap-wall-response-v1/manifest.json)
give

$$
\delta\overline K=0,\qquad
\delta a_+=\delta a_-=\frac{6m\rho}{bx^4}
-\frac{2(\mu_+-\nu_-)}{bx^3}.
$$

These are conditional identities for the coefficient of $\alpha$, including
both metric and wall displacement. They do not assert vanishing displacement
or all-order closure of the nonlinear trajectory.

An independent identity failure under the same action, preparation, domain,
orientation and perturbative order would reject this conditional result.
An observational or other independent numerical test still needs a defined
observable, physical input selection and a declared uncertainty. Once those
are supplied, a prediction outside that uncertainty rejects agreement. The
field-grid differences above must not be substituted for that future
acceptance band. No observational calibration has been performed.

The later [finite-amplitude onset theorem](nsc-atlas-finite-amplitude.md#finite-cell-continuation-assembly)
constructs the same candidate on a shorter initial source cell. It does not
reach this entry's terminal reference-radius prediction point. Its later
[quadratic amplitude estimate](nsc-atlas-finite-amplitude.md#second-amplitude-derivative-and-quadratic-onset-remainder)
also stays on that shorter cell. Its [completed uniform remainder and
strict fixed-label signs](nsc-atlas-finite-amplitude.md#certified-finite-amplitude-endpoint-signs)
now identify the shorter cell's first coefficients and bound its exact
Einstein-radius and anchored-clock response. They do not extend to P001's
terminal point. P001's finite-amplitude remainder and observational limits
therefore remain as recorded.

### Full-endpoint observable selected before the three-grid results

The current continuation fixes `alpha=0.00012983127142551443` and the same
source, then evolves to affine incoming label `x=0.4815042082576364`
(`V=1`, reference radius `0.55`). Its observable is the **Jordan material-wall
areal radius** on that cross-section. In the runner's dimensionless variables,

$$
R_{J,\mathrm{metres}}
=\ell\sqrt{r_E^2+\frac{\alpha c^2}{6}},\qquad
\ell=10000\ \mathrm{m}.
$$

The scale is a declared fiducial model input, selected before the endpoint
values; it is not an inferred astrophysical length. A spherical array comoving
with the wall could determine the radius from intrinsic area,
$R_J=\sqrt{\mathcal A_J/(4\pi)}$, or equatorial proper circumference,
$R_J=\mathcal C_J/(2\pi)$. The normalized incoming null label selects the
measurement cross-section. Reception by a remote or central detector is a
separate calculation. This fixed-label event must not be identified with the
fixed Jordan wall proper time used by the historical coefficient $\rho$.

The [completed three-grid comparison](../results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json)
now reports the fine-grid radius and the symmetric interval whose half-width
is the largest difference from either other grid. At the declared ruler:

| Cap/exterior cells | Jordan radius (m) | Maximum actual mass gap | Maximum wall-flux gap |
| --- | ---: | ---: | ---: |
| 2/32 | 5492.322036563843 | 0.00232832534127 | 2.24045303347 |
| 4/64 | 5492.322058752734 | 0.00237798074711 | 2.19447899469 |
| 8/128 | 5492.322035897825 | 0.00240631809367 | 2.14212667980 |

The fine value is **5492.322035897825 m**, with observed resolution half-width
**0.000022854908820590936 m** (22.85 micrometres). The full grid spread has
the same magnitude because the fine value is the smallest. Its symmetric
range is `[5492.322013042916, 5492.322058752734] m`. The shift from the quiet
5500 m reference is `-7.677964102174883 m`.

The wall-flux maximum decreases 4.39 percent, while the actual mass-gap
maximum increases 3.35 percent. Successive radius differences are 22.19 and
22.85 micrometres, so the radius pair gap does not contract either. The
registered finite-amplitude status is **unconverged numerical estimate**.
This is an observed resolution envelope, with the failed mass contraction
reported alongside it. It is not a confidence interval or a rigorous continuum bound.
Common source interpolation, radial numerical error, errors outside this
grid family and physical preparation remain separate uncertainties. An
independent relative uncertainty in $\ell$ transfers directly to the radius.
All three original evolutions have completed. The coarse evolution's
auxiliary accumulated-transport summary remains unavailable after its
post-evolution assembly failure. The historical perturbative result remains
unchanged; this extension does not close its full nonlinear endpoint claim.

The later [active-region diagnosis](../results/summaries/nsc-atlas-endpoint-active-resolution-v1.json)
locates the retained incoming layer in one unbisected outer quintic cell on
all three levels. The clipped endpoint cell narrows only 6.70% then 7.18%,
and the internal knots move. The integrated absolute mass residual grows
5.11%, so maximum-sampling density alone does not explain the failed trend.
This identifies an unresolved cap approximation, not a continuum error bar
or the unique cause of the unarchived peak. No additional grids were run.

The saved-data registration command is:

```sh
python scripts/run_python.py scripts/register_nsc_atlas_fixed_endpoint.py --comparison results/summaries/nsc-atlas-fixed-endpoint-comparison-v1.json
```

It authenticates the three completed summaries and the comparison before
adding an optional `finite_amplitude` block to NSC-P001. It preserves the
historical coefficients and their evidence. An unconverged comparison keeps
that status; a failed check leaves the machine register unchanged.

Check the register with:

```sh
python scripts/run_python.py scripts/check_prediction_register.py
```

The checker authenticates the seven evidence files, their manifest links,
the parameter/domain declarations and the exact saved-value projections.
It performs no symbolic derivation, numerical integration or generator run;
cross-platform numerical replay belongs to the separate evidence validators.
