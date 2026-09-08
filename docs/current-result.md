# Current result: the connected operator and its physical closure

Version 0.2.0 presents a 77-record working collection: the original 58 historical records and 19 scoped follow-ups. The [paper](../paper/nested-space-cosmology.pdf) organizes these results around six physical questions, with technical derivations and historical corrections in appendices.

## The strongest completed chain

| Connection | Reproducible result | Physical meaning and domain |
|---|---|---|
| Geometry to boundary response | [Finite radial maps](../results/nsc-3-boundary-response.json) | Oriented first-order maps, direct/Schur agreement and independent ODE checks |
| Smooth geometry to spectrum | [Smooth geometry](../results/nsc-4-smooth-geometry.json) | Seam-free periodic confinement and required Einstein stress; imposed geometry |
| Horizon to fermionic transport | [Dirac tetrad](../results/nsc-4-dirac-tetrad.json), [Lorentzian evolution](../results/nsc-4-lorentzian-transport.json) | Correct lapse/shift/spin connection and conserved norm on the fixed benchmark |
| Geometry to quantum stress | [Covariant vacuum controls](../results/nsc-4-covariant-measure.json), [shape response](../results/nsc-4-shape-response.json) | Absolute cylinder null response and varying-neck spin-structure stress difference |
| State to energy transfer | [Energy ledger](../results/nsc-6-energy-transfer.json) | Same link: zero static-vacuum injection and nonzero prepared-state transport |
| Geometric work to vacuum excitations | [Vacuum-work calculation](../results/nsc-6-vacuum-work.json) | Prescribed radius pulse creates pairs, with work and regional energy conserved |
| Sheet structure to Dirac chirality | [Observable bridge](../results/nsc-7-observable-bridge.json) | Exact invariant sector for a scalar sheet coupling; actual coupling and sector selection open |
| Recursive depth to a defined response | [Tail endpoint](../results/nsc-5-tail-limit.json) | Fixed spatial chain has a quantified endpoint criterion; no physical Omega selected |
| Source to cosmic evolution | [Plateau requirements](../results/nsc-5-plateau-conditions.json) | Internal exchange, pressure and external room supply cannot be conflated |
| Clock to physical history | [Clock/horizon audit](../results/nsc-5-clock-horizon.json) | Infinite coordinate time, affine continuation and mapped duration are distinct |

## Latest computed physical structure

The [full-spinor boundary result](../results/nsc-8-chiral-boundary.json) gives a nonzero, energy-dependent four-component spatial response with direct/Schur/continuum agreement. It preserves physical chirality in the declared domain. Its action kernel occupies the vector/axial Clifford subspace, so the scalar mass mechanism requires a derived additional interaction or domain selection. Opposite normal or frame conventions cannot supply this by themselves.

The [finite metric-response result](../results/nsc-8-finite-terms.json) identifies four independent bulk coefficient channels on the actual smooth profile. Euler and total-divergence terms have zero closed-cell bulk variation. The exact lapse/radius/shape sensitivities now specify which independent ultraviolet matching conditions the common functional must supply before solving self-sourcing.

## The next connection to close

The full spinor boundary map and its transmission domain must determine whether the physical chiral sector is admissible. In parallel, the common covariant functional and state must determine the absolute stress. These two tracks join at

$$
\frac{\delta\Gamma_{\rm one}}{\delta g^{\mu\nu}}=0.
$$

The calculation must retain every independent metric, field, boundary and relative-scale equation. The Einstein term already induced by the adopted construction is not added again with a new weight.

For the prescribed smooth neck, anomaly matching alone leaves finite invariant action terms undetermined. The written flat common-scale average also diverges in the tested finite prescription. A complete compensator, measure, state and ultraviolet matching are required to determine the source rather than fit it.

The actual boundary response must be transported into a common spin frame before identifying its scalar, pseudoscalar or other components. The scalar sheet benchmark supplies exact target identities, but the geometric radial potential preserves massless chirality in the paired angular continuum problem. No value of the current radial gap is labeled an electron mass.

## Physical source accounting

At a common dimensional energy, parent normalization gives

$$
\Gamma_p(x)=K_p(x)-\Omega^{-1}b\,\Gamma_c(x/\Omega)^{-1}b^\dagger.
$$

The retarded kernel alone does not determine an occupation state or energy current. Internal Q_b cancels from total room continuity; external supply J requires reservoir/boundary gravitational accounting. Pressure and perturbation response are needed in addition to the scalar transfer history. The [observational map](nsc-observational-targets.md) gives the required path to H(z), clustering and lensing.

## Historical results that constrain the construction

- [Unit correction](../results/nsc-2-zeta1-unit-closure-check.json): the additive-gap exponent is lambda/zeta + mu². The corrected finite-family sign theorem supersedes its old stationary-scale candidates.
- [Finite regulated calculus](../results/nsc-3-regulated-recursion.json): the proper-time relative Hessian is not removed by raw determinant invariance; a spatial trace and a covariant spacetime determinant differ.
- [Isolated radial spectrum](../results/nsc-3-radial-spectrum.json): the asymptotically widening spatial throat is gapless.
- [Threshold response](../results/nsc-3-threshold-response.json): a positive Dirichlet-to-Neumann jump can coexist with a gapless bulk spectrum.
- [Stellar benchmark](../results/nsc-1-gravitating-bps-observation-link.json): 3.34 solar masses is an imported branch from a nuclear-calibrated Einstein–BPS model, not an NSC-derived universal limit.

Physical scale selection, absolute self-sourcing, constrained metric stability, identified particle/nuclear sectors, detector probabilities and independent cosmological predictions remain the successive completion targets. The six-target map is maintained in the release specification and manuscript; no historical terminal flag promotes these open targets to solved physics.
