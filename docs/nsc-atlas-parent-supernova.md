# Inherited parent against supernova distance shapes

## Decision and preparation, 28 October 2026

The quiet cap's 7.278-sigma contrast and the controlled inherited matter
background answer different questions. The next calculation compares their
actual distance shapes, together with a LambdaCDM control, using the same
pinned Pantheon+ release and covariance. The latest request explicitly
authorizes that control; the earlier instruction excluding such comparisons
does not describe this task.

Three routes were considered: another local wall repair, adding an adjustable
dark component to the cap, and testing the existing inherited parent without
altering its source. The third connects the new source theorem to real data
and reveals whether the missing component is evolution, preparation, or
observer coverage. It introduces one fitted observer epoch, not a new dust
normalization. A magnitude intercept belongs to the data likelihood and is
profiled analytically. Controls introduce their own explicitly labelled
matter fraction. All state ownership stays local to the calculation.

Use `nsc_atlas_inherited_matter.exact_initial_geometry` and its uniform bound;
reuse `nsc_atlas_cosmological_distance.distance_contrast` and the original
publisher cache/verification. Extend the package with a small FLRW forward
map and common covariance likelihood. No new cosmology dependency or field
simulation is needed. Source parameters remain rho_v=9/32, E0=3/8 and the
exact junction curvature; the high-frequency member is omega=10^9.

1. Compare quiet de Sitter, fixed flat LambdaCDM fractions .3, .315 and .334,
   the inherited parent at A_o=2, and separately calibrated one-parameter
   LambdaCDM/parent histories on the existing 596-row selection.
2. Repeat the likelihood on non-calibrator rows with zHD>=.025 through the
   release maximum. Every inherited-parent emitter must have A_e>=1, so
   the fitted observer must satisfy A_o>=1+max(zHD). Do not extrapolate the
   protected future through its initial Cauchy slice.
3. Report both the unchanged one-dimensional contrast and the full selected
   covariance likelihood with a common magnitude intercept. Count fitted
   shape parameters and distinguish calibration from prediction.
4. Propagate the proved high-frequency density envelope to the light cone;
   check quadrature independently. A domain failure or a poor likelihood
   redirects the source/history construction instead of adding another
   component or silently modifying the data.

The calculation has a 120 process-CPU-second allowance, four BLAS threads,
at most 200 objective evaluations per fit, and zero field/grid evolutions.
Only compact summaries, hashes, producer and regeneration commands are saved.

Completion requires an identical data/zero-point/covariance protocol for both
models, a protected-domain check and a stated numerical/model error. A later
independent distance/growth comparison with the calibrated state held fixed
can falsify its physical interpretation; these already-seen data are a
retrospective comparison. Main risks are an inconsistent heliocentric
prefactor, extrapolation before A=1, and treating a calibrated epoch as a
parameter-free prediction. Each is an explicit code/test boundary.

## Remaining connected closures

The absolute SI ruler, inherited material-wall state, unconverged P001
endpoint, baryon sector and baryon asymmetry retain their separate owners.
A successful supernova shape cannot close them. Grok audits the retained
wall and endpoint equations without running new grids; their useful
consequences will be integrated with this observational result.

## First comparison and the required adaptation

The method was committed in b9f5aba before its evaluation. The retained
[future-only comparison](../results/summaries/nsc-atlas-parent-supernova-v1.json)
reproduces the original contrast and adds the missing control: flat
LambdaCDM with matter .30 and .315 gives standardized residuals -1.51735
and -1.25484. On the two bands the calibrated parent has chi-square
604.29927 versus 604.20151 for calibrated flat LambdaCDM, with the same
two fitted quantities (epoch/matter fraction and magnitude intercept).
This is calibration, not independent confirmation.

The 1322-row sample through zHD=2.26137 forces A_o>=3.26137 in the
future-only domain, hence Omega_chi<=.037845825. Its best permitted epoch
is on that boundary and chi-square is 1636.49112 versus 1178.20879 for
the one-shape-parameter flat control. This rejects the proposed use of
that *future-only observer coverage*, not the earlier source law.

The direct adaptation is a finite backward Cauchy continuation with exactly
the same scalar, momentum, curvature and couplings at A=1. Prove it down to
A=1/3 before making the broader comparison again. No source normalization
or flux is changed. The bootstrap, source and two-way causal proof belong
to `nsc_atlas_parent_past.past_bound`; two independent unaveraged clocks
at omega=200 check it. The high-frequency theorem remains analytic.
The full-sample observer stays at A_o>=1 and every emitter at A_e>=1/3.
This new domain is recorded as v2; the failed v1 comparison stays intact.
The finite allowance is 120 CPU seconds and 160000 RHS evaluations for
the two backward controls, followed by the existing likelihood allowance.

### Backward source and causal proof

Keep the original corrected energy I=A^3 E+(3/2)H A^3 chi Pi,
with dI/dt=(3/2) H_dot A^3 chi Pi and |A^3 chi Pi|<=A^3 E/omega.
Bootstrap 1/2<A^3 E/E0<2 on 1/3<=A<=1. The original Friedmann
constraint then gives 1/3<H<3: the lower comparison
rho_v/3+E0/(6A^3)-C/A^2 decreases up to A=1, and is greater than 1/9
there. The upper comparison is below 9 even at A=1/3.
Consequently the total backward Hubble variation is bounded by

$$
V_-\le\frac{C_{\max}}{2H_-}(A_{\min}^{-2}-1)
 +\frac{2E_0}{3H_-}(A_{\min}^{-3}-1)=19.7736,
\quad C_{\max}=.0228,\ H_-=1/3.
$$

For delta=3*3/(2omega), Gamma=3V_-/(2omega(1-delta)), the same positive
corrected-energy argument works in either time direction. It yields
L=(1-Gamma)/(1+delta), U=1/((1-Gamma)(1-delta)). When 1/2<L<1<U<2,
these strict bounds close the bootstrap. Finite positive A,H,E keep the
smooth homogeneous ODE regular; it reaches A=1/3 in finite proper time.
Combining with the existing future bound gives, for omega=10^9,

$$
\left|\frac{A^3E}{E_0}-1\right|<3.416041\,10^{-8}
\quad(A\ge1/3).
$$

The inherited cap's initial angular radius exceeds 2.83. Its entire future
conformal interval is less than 11/21, and its backward interval to A=1/3
is less than .912. A past ray from any parent-pole receiver at A_o>=1
has angle chi=eta_future(receiver)+eta_past(emitter). For that emitter to
be determined by the initial cap it needs chi+eta_past(emitter) below the
initial angular radius. The two-way estimate is

$$
\eta_{\rm future}+2\eta_{\rm past}
 <\frac{11}{21}+2(.912)=2.347810<2.83.
$$

Thus the whole specified ray remains in the smooth Cauchy domain without
constructing the walls' earlier or later interactions. This proves a finite
past continuation; it is not a formation/bounce or primordial preparation
mechanism. The scalar-to-distance magnitude bound is below 3.836e-8 up
to a common offset, provided the ray remains on the verified monotone sine
branch. The full covariance likelihood removes that common offset.

The [rational proof and unaveraged controls](../results/summaries/nsc-atlas-parent-past-v1.json)
retain exact fractions. Two omega=200 clocks reach A=1/3 with endpoint
comoving-energy difference 5.56e-10 and proper-clock Friedmann residual
5.40e-11, using 39544 RHS calls. Their finite numerical agreement checks
the calculation; the high-frequency bound follows from the bootstrap.

Grok independently checked the supplied equations and two-way domain argument
without reading the paper. Its review distinguishes the backward variation
budget from the future: the 19.7736 budget only covers 1/3<=A<=1. The
implementation explicitly combines that result with the already-proved
future `uniform_bound`, so the stated full A>=1/3 envelope uses both.
The observational rays satisfy the conservative angle<1 check, inside the
larger sine-monotonicity domain.

## Supernova result after the controlled continuation

The [v2 comparison](../results/summaries/nsc-atlas-parent-supernova-v2.json)
retains both failed v1 and all original source inputs. The full 1322-light-
curve selection gives chi-square 1179.38727 for the inherited parent and
1178.20879 for flat LambdaCDM. Each has 1320 degrees of freedom after
one shape calibration and one magnitude intercept. The inherited epoch is
A_o=1.3816545682; its earliest emitter has A=.4236424, safely above 1/3.
The source fixes (Omega_chi,Omega_v,Omega_k)=(.3667843,.7255551,-.0923394)
at that calibrated epoch. These are model fractions inferred within this
comparison, not independently predicted cosmic abundances. The pure quiet
cap remains discrepant; it has not been relabelled as the evolving source.

Both scalar normalization and curvature remain their original values. The
analytic per-distance magnitude envelope is 3.83533e-8 up to a common
offset; the contrast envelope is 8.13220e-8 mag. The quadrature indicator
is 1.21e-15 mag. The finite source/geometry model is conditional, and no
observational uncertainty in the original preparation has been marginalized.

The data and method follow [Brout et al., Pantheon+](https://arxiv.org/abs/2202.04077)
and its [public distance/covariance release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR).
Only the explicitly selected low-redshift bands and the stated zHD>=.025
sample are analyzed here, not the complete collaboration inference chain.

## Next independent geometric channel: DESI AP shape

Freeze the supernova-calibrated epoch and flat-control matter fraction before
evaluating DESI. Use all six anisotropic pairs at z=.51,.706,.934,1.321,
1.484,2.33 from [DESI DR2 II](https://arxiv.org/html/2503.14738v3), Section
III.3, and the public Cobaya mean/covariance at commit
bb0c1c9009dc76d1391300e169e8df38fd1096db. The
[DESI source-provenance clarification](https://help.desi.lbl.gov/index.php?qa=182&qa_1=provenance-external-cobaya-likelihood-public-cosmology-chains)
identifies those data with the collaboration likelihood. No new package,
sound-horizon fit, SI ruler, early thermal history or additional shape
parameter is introduced.

For fixed F(z)=D_M(z)H(z)/c, take the six *linear* contrasts
y_DM-F*y_DH and covariance W C W^T. This avoids assuming that an observed
ratio is Gaussian, and retains cross-bin correlations. It equals profiling
the Gaussian likelihood over an arbitrary unknown radial amplitude in each
bin. The isotropic DV row has no AP pair and is excluded. This is a test
conditional on the fixed supernova point estimate; its uncertainty is not
silently omitted from a claimed joint significance. Public BAO data are a
retrospective external comparison, not newly blinded measurements. The
forecast/method are saved before the numerical BAO comparison; no BAO
retuning is allowed. Its finite compute allowance is 10 CPU seconds and
zero field evolutions. A significant discrepancy redirects the inherited
curvature/preparation law rather than the supplied acoustic ruler.

The [forecast](../results/summaries/nsc-atlas-parent-bao-forecast-v1.json)
was pushed in 18cfdf6 before the numerical DESI comparison. The
[conditional external comparison](../results/summaries/nsc-atlas-parent-bao-validation-v1.json)
gives chi-square 7.15813 for the inherited parent and 8.89401 for the
SN-calibrated flat control, each for six contrasts and with no BAO fit.
The scalar AP-relative envelope is 3.53247e-8. The small difference in
statistics is not evidence of model preference; the supernova calibration
uncertainty has not been marginalized. A common acoustic scale across
redshifts, the full BAO likelihood and early acoustic physics are additional
tests. The chosen AP projection deliberately profiles out a separate unknown
amplitude in each pair; its numerical control verifies that equivalence.

## What each of Douglas's five gaps now requires

1. **SI scale.** The background family has a continuous dimensional freedom:
   scaling (rho_v,E0,C,H^2) by lambda^2 leaves every fraction and the
   epoch-independent correlation unchanged while distances scale as
   lambda^-1. In the full classical scalar action this corresponds to
   scaling the *supplied* masses, mixing and fluxes, with the matching
   radii/masses/tensions scaled consistently; it is not a symmetry changing
   fixed physical couplings. The inherited relation and scale-free data do
   not choose those couplings in SI. A scale-sensitive preparation mechanism
   or independently calibrated dimensional input is still needed. G and hbar
   permit Planck units but do not by themselves select the dimensionless
   flux/mass ratios. No absolute Lambda has been generated here.
2. **Inherited wall state.** The stable gravitational-response branch uses
   a massless surface scalar, sigma=p and sigma R^4=constant. Its shift
   charge is Q=4pi R^2 theta_dot, so sigma=Q^2/(32pi^2 R^4). With no
   bulk-to-surface coupling Q is independently conserved on each wall;
   a zero initial Q stays zero. This derives that EOS but does not inherit
   its state. The evolving parent instead has constant charged tension,
   p=-sigma. Importing the stable branch's Q into it changes the action.
   Derive a compatible surface action/preparation or evolve the existing
   membrane; do not transfer the earlier stability label.
3. **P001.** Grok's source audit confirms the unbisected active outer cell,
   the growing strong incoming-sum mismatch and the closed mass reconstruction
   identity. The new [matched-trace test](../results/summaries/nsc-atlas-matched-trace-v1.json)
   checks adapter compatibility independently. The 3.35% production growth
   remains; a trace-compatible update is required before any new grids.
4. **Observed expansion.** The missing identical-data LambdaCDM control is
   now present, and the sourced parent's controlled earlier continuation
   fits the broad sample comparably with the same number of fitted parameters.
   Its fixed supernova-calibrated geometry also passes the conditional DESI
   AP consistency check. This advances the physical observer connection,
   without establishing a unique theory or primordial selection.
5. **Composition and asymmetry.** This action contains one real scalar
   excitation and a four-form, not a baryon current or photon distribution.
   Omega_chi is therefore neither a derived Omega_c/Omega_b split nor a
   baryon-to-photon ratio. Eta_B is undefined in this sector, not calculated
   to be zero. The next physical extension must identify the conserved/
   violated current, production and entropy history before calculating it.

The immediate outward test is a joint distance/growth/acoustic-scale account
with the source preparation held fixed. The physical construction behind it
must select that preparation, SI scale and surface state. Those source tasks
have priority over more digits in the already controlled background or another
P001 grid family with the same trace defect.
