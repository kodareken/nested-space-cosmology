# Prediction register

This register states what the fixed model predicts, where it applies and how
agreement could be rejected. Each entry distinguishes derived consequences,
numerical sensitivity and observations. The [machine-readable register](../results/prediction-register.json)
binds the evidence bytes and the exact locations of reported values.

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

The [three-grid comparison](../scripts/analyze_nsc_atlas_fixed_endpoint.py)
will report the fine-grid radius and the symmetric interval whose half-width
is the largest difference from either other grid. This is an observed
resolution envelope, with the measured mass/flux contraction reported
alongside it. It is not a confidence interval or a rigorous continuum bound.
Common source interpolation, radial numerical error, errors outside this
grid family and physical preparation remain separate uncertainties. An
independent relative uncertainty in $\ell$ transfers directly to the radius.
The endpoint value and resolution range are pending the three completed
evolutions; the existing perturbative record is preserved meanwhile.

Check the register with:

```sh
python scripts/run_python.py scripts/check_prediction_register.py
```

The checker authenticates the seven evidence files, their manifest links,
the parameter/domain declarations and the exact saved-value projections.
It performs no symbolic derivation, numerical integration or generator run;
cross-platform numerical replay belongs to the separate evidence validators.
