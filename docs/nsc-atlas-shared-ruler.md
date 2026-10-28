# One shared DESI ruler and the inherited scale question

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
