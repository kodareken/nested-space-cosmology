# Inherited SI scale: computed routes and remaining selector

The target is Lambda_obs=1.0904e-52 m^-2. Checkpoint 485b3b9, its couplings,
regional state and P004/P005 remain unchanged. A failed candidate is not
evidence against nesting. Each local failure below gets one line; the
calculation continues to the next candidate.

## What is fixed before an SI ruler is supplied

The existing source fixes Lambda*L_unit^2=9/32, a0=6.626077076099598 and
Lambda*L_ref^2=12.348252398928544, where L_ref=a0*L_unit is the initial
parent curvature radius. The owner normalization is
phi_hat=sqrt(kappa_nat)*phi, q_hat=sqrt(kappa_nat)*L_unit*q,
e_hat=sqrt(kappa_nat)*L_unit*e and omega_hat=L_unit*omega/c.
Arithmetic kappa=1 has absorbed the gravitational coupling; it does not
identify the numerical simulation unit with one Planck length.

Use CODATA 2022 only for unit conversion: lP^2=hbar*G/c^3 and
kappa_nat=8*pi*lP^2. In mixed SI stress conventions kappa_SI=8*pi*G/c^4,
so kappa_nat=kappa_SI*hbar*c. Substituting kappa_SI directly as a squared
length would be a dimensional error. The independent route audit found
that issue in a draft; no draft SI numbers are promoted to the record.

The computed [route ladder](../results/summaries/nsc-atlas-SI-scale-ladder-v1.json)
and its tests keep target-dependent inverse diagnostics separate from
target-independent conditional predictions.

## A surviving microscopic SI relation, distinct from selecting L_ref

The [coupling-ratio calculation](../results/summaries/nsc-atlas-inherited-microanchors-v1.json)
makes an essential distinction: the free parent curvature scale does not
make every SI coupling combination free. The inherited source has
sigma_hat=kappa_nat*L_unit*sigma_phys and
e_hat=sqrt(kappa_nat)*L_unit*e_phys. Therefore

$$
\frac{e_{\rm phys}}{\sigma_{\rm phys}}
=\frac{\hat e}{\hat\sigma}\sqrt{\kappa_{\rm nat}}
=2.0256752716\times10^{-34}\ {\rm m}.
$$

This microscopic canonical charge/tension length is fixed by the supplied
ratios and G,hbar,c; it is not selected using Lambda_obs.
Canonical charge and tension here are expressed in inverse square and
inverse cubic metres after the hbar=c=1 conversion; the ratio is a length
in that convention. An unconverted SI surface energy must first acquire
its hbar*c factor. These definitions are those of the declared action.
Other fixed SI combinations are (omega_phys/c)/e_phys=1.6205402173e-25 m and
sqrt[(omega_phys/c)/sigma_phys]=5.7294748842e-30 m. They all remain fixed
in the explicit dimensional-completion witnesses while the parent Lambda
changes. Thus they survive as conditional coupling relations, not as a
cosmological-radius closure. The needed parent/microscopic hierarchy would
be 1.6612691531e60 if the target were inserted. That hierarchy is unselected.
In the ladder record, `inherited_SI_scale_selected` refers to the parent
reference scale, not to these microscopic combinations. The target goal
remains unachieved.

## Ranked route ladder

1. **Flux/charge quantization:** q/e=3, screened vacuum levels (9/32,1/8,1/32,0), second difference 1/16; fundamental SI charge and L_ref remain free.
2. **Inherited wall EOS/radius:** e_phys/sigma_phys fixes 2.0256753e-34 m, but frozen stationary radii (1.846303,1.155844) give V_second=(-1.28718,-11.32694); no cosmological-radius hierarchy is selected.
3. **Recursive fixed point/level:** vacuum-length ratios are 1.5 then 2 and the flux-map fixed-point residual is -1/2; no constant geometric ratio/root length or full M,sigma,ell closure is selected.
4. **Scalar frequency/epoch:** Lambda*c^2/omega^2=2.8125e-19 and the calibrated equality redshift is .331586; an SI scalar frequency is absent.
5. **Compact four-cycle period:** the formal added S4 premise gives Lambda*lP^2=1/n; n=1 gives 3.8280731e69 m^-2, while no inherited winding or admissible cycle/contour for the excited parent has been derived.
6. **Semiclassical transition weighting:** isolated pure-flux action coefficients (.6956584,3.1684104,70.9872067) still multiply L_unit^2/kappa_nat and are not the interacting parent's rate or a selector for L_unit.
7. **Quantum occupation/entropy:** a hypothetical one-quantum normalization gives Lambda=9.17759e61 m^-2, but the actual real classical scalar supplies neither that number state nor an inherited occupation.
8. **Anomaly/running/transmutation:** prior added classical vacuum trace is zero; a quantum stationary scale needs a compatible renormalized action/state and subtraction condition which the fixed source does not specify.
9. **Induced gravity/Casimir/compactification:** the stable alternative body's Lambda*M^2=.1144906 still requires an SI mass/modulus/boundary state and uses a different surface action.
10. **Infinite hierarchy/global similarity:** ell_n=r^n ell_0 transports a scale without selecting ell_0 or the observer level.
11. **Legacy S-one spectral fixed point:** the old gap threshold 3*pi/2 and alternative coefficient 18*pi retain a supplied length_anchor and no solved physical fixed point; that action is not the fixed parent.

These are all route classes identified in the existing atlas/source inventory,
not an assertion that every future physical mechanism has been exhausted.
The [stationary-wall calculation](../results/summaries/nsc-atlas-wall-scale-stationary-v1.json)
uses the existing constant-tension potential at the actual initial field.
Its instantaneous frozen representation is not a coupled-future stability
result. The separate stable surface-scalar action is not imported as its EOS.

## Full evaluation of the quantization candidate

Brown--Teitelboim and its quantized-flux extension are reused from
[Bousso--Polchinski, Section 2](https://arxiv.org/html/hep-th/0004134).
Their jump condition and absolute/period quantization are distinct.
In the present scalar source the late actual four-form is F_hat=epsilon*q_hat,
while the displacement q_hat has the electric integer q_hat/e_hat=3.
That integer does not fix a fundamental charge in metres.

The strongest explicit added period candidate assumes a compact Euclidean
S4 of radius sqrt(3/Lambda), an admissible real integration contour and
the period e*integral(F4)=2*pi*n. This is a candidate global preparation,
not a constructed continuation of the excited three-region initial atlas.
Its unit-volume factor is 24*pi^2/Lambda_hat^2. With the owner's field
normalization the period is

$$
\hat e\,\epsilon\hat q\,
\frac{24\pi^2}{\hat\Lambda^2}\frac{L_{\rm unit}^2}{\kappa_{\rm nat}}
=2\pi n,
\qquad
\Lambda\ell_P^2=\frac{3\hat e}{n\hat q}=\frac{1}{n}.
$$

No integer is selected from the target. The minimal illustrative winding
n=1 gives L_ref=5.6795331e-35 m and Lambda=3.8280731e69 m^-2,
3.5107054e121 times the target. Matching the target would require
n approximately 3.5107054e121, which the supplied inheritance does not
produce. The electric integer 3 and this four-cycle winding are different
objects. Equating them is an extra condition, and would still give a
Planck-scale Lambda. Lorentzian/Euclidean contour compatibility and the
actual global saddle must be established before this formal family is a
physical NSC spectrum. This candidate does not survive the requested test.

The existing BT instantons are borrowed, not re-derived. Their finite
pure-flux controls are [Hirano's](https://arxiv.org/html/1804.09985) zero-mass
vacuum setting, not P002 or the interacting scalar parent. For illustration,
inserting the target-derived ruler would make the first exponent
2.7330216e119 in hbar units. Its enormous value is not a selected scale.

## Why assigning the target is not the requested result

Inverting the target would give L_unit=5.0787091733e25 m and
L_ref=3.3651918430e26 m. This would merely force Lambda/Lambda_obs=1.
It is stored as an inverse diagnostic, not a prediction. The corresponding
scalar angular frequency would be 5.9029263e-9 s^-1 (period 33.7294 years),
which is an independent possible physical clock to investigate, not a
measured or derived inherited clock here.

Two exact dimensional completions demonstrate the remaining freedom:
L_ref=1 m and L_ref=2 m reproduce the same normalized q/e=3, charge .5,
flux 1.5 and source geometry, but give Lambda=12.3482524 and 3.08706310
m^-2. They are mathematical witnesses, not proposed astronomical values.
The physical quantum-action normalization changes with this dimensional
completion; no quantum selector is supplied by the frozen classical tests.

No audited route yet returns an inherited, uniquely SI-fixed cosmological
scale. The target remains unresolved. This statement concerns the missing
selector in the fixed realization; it is not a rejection of nesting.

## Proposal for Douglas's review: derive and select the physical charge unit

Status: proposal only. No new selection sector, potential or numerical
campaign is implemented or admitted by this note. The existing source,
likelihood results, P004/P005 and scalar/flux couplings remain unchanged.

The unresolved quantity can be expressed as one dimensionless coupling,
g_scale=kappa_nat/L_unit^2, or eta_e=kappa_nat*e_can=ehat*sqrt(g_scale).
For this source Lambda=(9/32)*g_scale/kappa_nat. The target would require
g_scale=2.5453818753e-120 and eta_e=7.9771264803e-61; these are diagnostics,
not values supplied to a selector. The microscopic e/sigma length is
already fixed while this charge magnitude remains free.

The recommended research construction is one physical geometric modulus
of an enclosing parent that controls the canonical charge normalization,
together with its derived dynamics/state. This is new physical content:
the present four-dimensional action treats e,m,sigma as parameters. A
uniform change of units is not a modulus; varying L_unit at fixed hats
changes the physical coupling family. A normalized Gaussian state exists
for many couplings and does not by itself select one.

The alternatives considered are a full-cycle winding/state rule (still
missing its actual cycle, contour and inherited integer), an independently
measured SI clock (a calibration rather than inheritance-alone), and more
classical wall/shape work (retains the scale freedom). The geometric
charge-normalization construction directly addresses the free coupling.
Its credibility depends on deriving the missing sector rather than
postulating a potential with the desired minimum.

Proposed execution after review:

1. Specify a physical parent geometry and one gauge-inequivalent modulus s.
   Derive its kinetic coefficient and the three-form kinetic normalization
   Z(s). Preserve gauge invariance by keeping the topological charge fixed;
   the canonical charge would be e_can=e_UV/sqrt(Z(s)). Derive the Planck
   normalization of e_UV too. Setting e_UV=1 would leave a supplied scale.
2. Derive a selection action Gamma_sel[s] from that same parent and one
   explicit inherited-state prescription. Every geometric, flux, boundary
   and quantum term and every subtraction condition must have an owner.
   If a controlled local potential exists, solve V_eff'(s*)=0 with positive
   kinetic norm and V_eff''(s*)>0. A nonlocal causal action needs its actual
   stability condition instead. Do not assume a minimum or tune coefficients
   using Lambda_obs. State normalization alone is not a selecting equation.
3. Match canonical charge, tension, scalar frequency and the total vacuum
   stress back to the fixed normalized parent. Include the new sector's
   gravitational stress and vacuum energy; do not subtract its minimum by
   hand. A mismatch defines a separate candidate rather than an invisible
   modification of 485b3b9.
4. Only when these equations exist, solve the resulting one-mode problem,
   compute g_scale, L_unit, L_ref and Lambda in SI with uncertainty, then
   compare the target and predict a second observable with calibration held
   fixed. No additional cosmological fit is needed to formulate the selector.

Implementation ownership would remain here and in the atlas owner. The
current scalar/flux modules are frozen matching references. A prospective
new module `nsc_atlas_parent_scale_selector.py` would own only the derived
physical scale sector and its checks; it is not created until the proposed
construction is approved. No extra physics library or field-grid campaign
is justified for this first derivation.

Goldberger--Wise's [bulk/boundary modulus stabilization](https://arxiv.org/html/hep-ph/9907447v2)
is a method precedent in a different geometry. Bousso--Polchinski's
[charge normalization and compactification discussion](https://arxiv.org/html/hep-th/0004134#S2.SS5)
identifies microscopic charge and geometry dependencies. Neither paper
supplies this NSC sector or its numerical parameters.

Completion criteria for the proposed selector: one physical scale degree
with a derived selecting equation; no adjustable SI coupling, root integer,
initial modulus value or finite vacuum offset fixing the answer; an isolated
stable state/value with controlled approximation; and explicit matching to
the frozen parent. A claimed prediction is compared with Lambda_obs only
after it is computed. Its second-observable consequence can falsify that
candidate. The existing dimensionless 3*Omega_Lambda/curvature tensions are
separate: changing the overall scale cannot remove them.

Failure forecast:

- The proposed scale mode is gauge or absent from the parent action: the
  construction needs an explicit geometric/theory extension, not a variation
  of a unit choice.
- Z(s), e_UV, a boundary state or a renormalized constant remains freely
  selectable: the hierarchy has been renamed rather than derived.
- The derived state is unstable, gives the wrong scale or fails normalized
  matching/its second observable: record that candidate's failure and move
  to a separately specified mechanism without retuning it to the target.
