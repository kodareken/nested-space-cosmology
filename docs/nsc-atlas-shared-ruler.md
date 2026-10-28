# One shared DESI ruler and the inherited scale question

## Computed result and scoreboard

The [one-table scoreboard and route ladder](nsc-atlas-shared-ruler-scoreboard.md)
is the reader entrypoint. The [shared-ruler record](../results/summaries/nsc-atlas-shared-ruler-v1.json)
uses all 13 DESI entries and the full SN posterior, with no source retuning.
The joint minimum chi-squares are 1203.22418 and 1191.17680 for the parent
and flat control, respectively: Delta=12.04739, with three fitted quantities
in each. The posterior-predictive BAO deviances are 27.24410 and 16.46889;
their difference is 10.77521 and is not assigned a chi-square distribution.
The 48/96-order marginal-deviance changes are below 1.8e-11. The original
SN scores and six AP scores are preserved at checkpoint 485b3b9.

The parent joint fractions are (.76310716,.32320464,-.08631180) in
(Omega_Lambda,Omega_m,Omega_k) order. The posterior 3*Omega_Lambda is
2.288667 +/- .021121, compared with the supplied 2.054 +/- .022: a
standardized difference 7.695 including both widths. This comparison is
not added again to the joint data likelihood.

The [scale-selection owner](nsc-atlas-scale-selection.md) evaluates eleven
identified route classes. None currently selects metres from the fixed
source. That statement is local to this realization; it is not a rejection
of nesting. The numerical target is never used to choose a free integer,
charge, reference length or counterterm.

## Conditional CMB geometry, with its missing acoustic sector explicit

The fitted curvature triggers the requested check. A full physical Planck
likelihood is undefined for this fixed source and the shape-only control:
neither specifies photons, baryons, a recombination surface or r_s(z_star).
The [computed diagnostic](../results/summaries/nsc-atlas-cmb-geometric-diagnostic-v1.json)
therefore supplies z_star=1089.92 as an external probe label and profiles
the unknown acoustic ratio. No such ratio is called inherited or predicted.

The CMB geometry itself can be continued without a new source. Bootstrap
1/2<A^3E/E0<2. For 0<A<=1, the reference polynomial
E0/6-C*A+rho_v*A^3/3 has minimum
E0/6-(2C/3)*sqrt(C/rho_v) greater than .24^2. Thus
H>.24*A^(-3/2). The backward Hubble-variation bound is

$$
V_-\le\frac{2C}{.24}(A_{\min}^{-1/2}-1)
+\frac{4E_0}{3(.24)}(A_{\min}^{-3/2}-1).
$$

The corrected-energy comparison closes this bootstrap at the actual joint
observer's probe A_min=.00134345355. The density envelope is 7.871e-5 and
the H^2 envelope is 9.625e-5. A point-specific two-way causal estimate has
margin .226636 radians inside the original cap. This is a geometry/source
continuation, not a radiation-era or recombination model. No carrier grid,
scalar-cycle evolution or wall history is computed.

Using the same curvature-relaxed reconstruction of
[Chen--Huang--Wang Table 1](https://arxiv.org/html/1808.05724v1) for both
columns, the mean is (R,lA)=(1.7429,301.409), uncertainties (.0051,.091),
and correlation .54. The two-dimensional covariance is marginalized over
unpredicted baryon density and tilt; it is not a subblock of the precision
matrix. These rounded table numbers are a reconstruction, not a new chain
product. The compression's standard-CMB assumptions remain conditional.

Profiling the unpredicted lA gives chi2_min=(R-1.7429)^2/.0051^2. The
joint-shape geometric probes give R=1.55147562 and 1.75179174, hence
conditional diagnostic scores 1408.82 and 3.04. The unaveraged scalar
phase enclosure still gives an optimistic parent lower bound 1407.72.
The fitted diagnostic sound ratios r_s(z_star)/r_d would be .841850 and
.982668. They are additional nuisance calibrations, not calculated acoustic
physics. These scores are not added to the primary SN+DESI total, and
are not labelled a physical CMB prediction or a general test of nesting.

The primary source explicitly includes radiation, baryon loading and a
decoupling prescription. Those ingredients must be produced by a compatible
extension before this geometry diagnostic becomes the requested full CMB
likelihood. Borrowing its compression does not supply those sources.

## Fixed question and fastest calculation

Continue checkpoint 485b3b9. The original charge, flux, curvature, scalar
normalization, P004/P005 and recorded supernova/AP comparisons stay fixed.
The new question uses all 13 DESI DR2 distance entries and their full
covariance, rather than the six amplitude-free AP contrasts. Brown--Teitelboim,
Israel, Hamilton--Lisle and Turner remain cited inputs, not new derivations.

Rethinking decision: a new field simulation, a multi-parameter cosmology
sampler and a one-dimensional posterior integral were compared. The last
has the shortest dependency chain: reuse the existing light-cone map,
analytically eliminate one BAO scale, and integrate the actual supernova
shape likelihood. No new science package, field grid or source parameter is
introduced. Parallel workers evaluate the two models and one quadrature
control. A failed candidate records a missing selector and the next route;
it does not reject nesting.

The NSC shape parameter is N_o=ln(A_o), with uniform prior on [0,ln(20)]
restricted to the proved emitter domain A_e>=1/3. The flat control has a
uniform Omega_m prior on [0,1]. The common supernova magnitude intercept
is integrated/profiled with the original covariance; its normalization is
independent of shape. Integrate the entire finite prior, including tails,
using piecewise Gauss quadrature around the likelihood peak. Record any
supernova posterior mass outside the additional DESI emitter domain before
conditioning on that domain. Repeat at double quadrature order. This is
posterior propagation, not a Gaussian substitute or a best-fit-only check.

Write t(theta) for the 13 dimensionless DV, DM and DH predictions in c/H0
units, y for the released distance vector, and beta=c/(100 r_d h). Then

$$
\chi^2_{\rm BAO}=a\beta^2-2b\beta+c_0,\qquad
a=t^TC^{-1}t,\quad b=t^TC^{-1}y,\quad c_0=y^TC^{-1}y.
$$

At each shape beta_hat=b/a and chi2_profile=c0-b^2/a. The marginal BAO
likelihood integrates exp(-chi2_BAO/2) against the normalized SN shape
posterior. Also report the minimum joint SN+BAO chi-square and its SN/BAO
parts, with three fitted quantities (one shape, one magnitude intercept,
one ruler). Posterior moments integrate the positive-beta Gaussian exactly;
the reported r_d h interval uses its mixture CDF, not a reciprocal mean.
Use a common flat beta prior for that joint posterior; these normalizations
are not reported as a model Bayes factor. No BAO target retunes the inherited
source or the older frozen results.

The initial 24/48-order control found a 2.26e-4 difference in the marginal
BAO deviance, despite agreeing on the joint optimum and fraction moments.
The conditional ruler likelihood is narrower in shape than the SN posterior.
The admitted comparison therefore uses 48/96 nodes per panel with the same
1e-5 numerical threshold; no physical parameter or error threshold changes.

## Scale route ladder and CMB boundary

Audit every scale-selection class identifiable in the current atlas: flux
quantization, wall radius/EOS, recursion fixed points and level choice,
scalar-clock matching, quantum/action normalization, spectral anomaly or
dimensional transmutation, induced gravity, Casimir/compactification and
global constraints. Separate a scale in model units, a derived L/l_P, and
an empirical SI calibration. The target Lambda=1.0904e-52 m^-2 is a
comparison only. Inverting it to choose a length, integer or counterterm
does not produce an inherited prediction.

If the fitted |Omega_k| exceeds .002, evaluate the Planck18 (R,lA) request
against its actual input dependencies. A last-scattering distance needs
an earlier source history; its acoustic length additionally needs baryons,
photons and recombination. Neither a fitted r_d h nor the A>=1/3 theorem
supplies these. Any supplemental thermal completion is labelled with its
added assumptions and cannot be silently counted as the unchanged parent.

The allowance is 180 process-CPU seconds per model/order worker, four useful
workers, and no field evolution. Success is one common likelihood scoreboard,
posterior propagation controls and a computed route ladder. The SI target
is complete only if a surviving route selects metres without using the
target. Likely failures are a hidden ruler/thermal input, insufficient
posterior quadrature, or an incompatible action imported from another
atlas branch. Each has a separate dependency/control check; no failed route
authorizes an unchanged grid campaign.
