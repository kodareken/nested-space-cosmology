# Current result

**The current frontier is the unsolved energy-resolved recursive child tail.**

This page is the public checkpoint for the committed compact JSON chain. It is not a proof of Nested-Space Cosmology, not a final \(\zeta\), not a particle spectrum, not a dark-sector fit, and not an identification of our universe with a black-hole interior. Labels follow [THEORY.md](../THEORY.md). Imported launch surfaces are in [docs/prior-art-and-open-claim.md](prior-art-and-open-claim.md). Regeneration is in [docs/reproducing.md](reproducing.md).

The live compact record for the unsolved tail is `results/nsc-2-zeta1-recursion-map.json` (`NSC-2-ZETA1-RECURSION-MAP`). Scoped follow-ups are `results/nsc-2-zeta1-unit-closure-check.json` and the nsc-3 triple `nsc-3-regulated-recursion.json`, `nsc-3-radial-spectrum.json`, `nsc-3-geometric-chain.json`. The accepted closure plan is not completed. The recursion-map classification is that norm preservation and first-order Dirac scaling fix a minimal parent/child dilation and the energy-resolved recursive tail. Dilation alone does not derive \(\zeta=\Omega^2\). Its `nonclaims` are all `false`:

- `Omega_value_selected`
- `functional_tail_solution_computed`
- `physical_zeta_promoted`
- `scale_root_restored`

Its `gate` likewise records `mode_resolved_tail_solved: false` and `zeta_derived: false`. The recorded next result is to solve the functional tail on the warped mode spectrum and recompute the anomaly-compensated scale derivative.

## What is now derived

**Repository derivation.** The dilation

$$
(U_\Omega\psi)(x)=\Omega^{d/2}\psi(\Omega x),
\qquad
\mathcal{T}_\Omega=\sigma_1 U_\Omega
$$

preserves the \(L^2\) norm, composes as a semigroup, and scales first-order spectral quantities together, so \(E/\Lambda\) and \(\Phi/\Lambda\) are inherited on a transformed mode. It does not by itself derive \(\zeta=\Omega^2\). The quantities \(\zeta=(\Lambda L_\star)^2\) and \(\Omega=\Lambda_{\mathrm{child}}/\Lambda_{\mathrm{parent}}\) remain distinct until the geometry relates them. With \(b=B_{\mathrm{dim}}/\Lambda_{\mathrm{parent}}\), the recursive outside equation must be solved as

$$
\Gamma_p(x)=K_p(x)-\frac{1}{\Omega}\,b\,\Gamma_c(x/\Omega)^{-1}b^{\dagger}.
$$

The unweighted form \(\Gamma(x)=K(x)-B^{\dagger}\Gamma(x/\Omega)^{-1}B\) absorbs the factor into a symmetric link convention. At vanishing argument this recovers the earlier quadratic fixed point \(\Gamma=K-B^2/\Gamma\) on the branch continuous with \(\Gamma\to K\) as \(B\to 0\). An implementation fixture with \(\Omega=2\) then \(\Omega=3\) composes to \(\Omega=6\); that fixture does not select nature’s \(\Omega\).

This map removes an undefined symbol from the one equation. It does not compute the mode-resolved tail.

**Repository derivation.** An explicit finite proper-time matrix calculus supplies noncommuting first and second variations and a finite scale cocycle (`results/nsc-3-regulated-recursion.json`). An ultrastatic frequency factor differs from a spatial heat trace. The trapped coordinate-time Hamiltonian is not elliptic; that general phenomenon is documented by Finster and Röken, arXiv:1512.00761. Along common-scale dilation at fixed cutoff the finite action is strictly monotone, so that family has no isolated stationary scale.

**Repository derivation.** A Weyl sequence proves that the isolated unwarped radial throat has essential spectrum \(\mathbb{R}\) (`results/nsc-3-radial-spectrum.json`). Smooth gluing alone does not supply a mass gap.

**Repository derivation.** Periodic repetition of \(\rho\in[-R,R]\) with \(w=1/\sqrt{1+\rho^2}\) produces a band gap without inserted mass (`results/nsc-3-geometric-chain.json`). Continuum first band edges are \(0.7034881641\), \(0.4552655377\), and \(0.2479700183\) at \(R=2,4,8\), with second-order lattice agreement. The link is a stencil entry. \(R\) and \(\Omega\) are inputs; seams have unproved gravitational stress. This is a controlled project realization, not a Lorentzian nested cosmology, not an electron, and not a new-to-world Floquet discovery.

## What the anomaly invalidated

**Repository derivation.** The previous scale candidates were owned by the regulated joined-minus-disconnected fermionic determinant, including after exact compact \(y\) warp. Compensating the pure cutoff anomaly of the invariant partition function removes those roots.

At the warped determinant candidate \(\zeta=4.748389947082489\):

| Piece | Value |
|---|---:|
| Determinant derivative | \(0.01201522080638507\) |
| Pure cutoff anomaly | \(12.336590431250533\) |
| Physical child-link derivative | \(-12.324575210444092\) |

The compensated child-link derivative stays negative on the scanned domain \(3\pi/2<\zeta\le 200\), from about \(-19.42\) to about \(-0.588\). Compact record: `results/nsc-2-zeta1-anomaly-decomposition.json`. Those historical numbers used the inconsistent conversion \((\lambda_j+\mu^2)/\zeta\).

**Repository derivation.** The declared additive-gap operator requires \(q_j=\lambda_j/\zeta+\mu^2\). At the same warped candidate, the unit-consistent subtracted derivative is about \(-36.58293\). For the inherited cutoff subtraction at frozen geometry, Hermitian interlacing proves that derivative is strictly negative for every \(\zeta>3\pi/2\) in this finite family. A diagnostic root of the unsubtracted determinant near \(\zeta=6.09675\) still has subtracted derivative about \(-14.15\), so it does not solve the adopted subtracted scale equation. Compact record: `results/nsc-2-zeta1-unit-closure-check.json`.

**Therefore those determinant-only \(\zeta\) values are invalidated as stationary points of the one equation.** They remain numerical diagnostics. Refining the same finite family cannot restore a root.

**Repository derivation.** The missing owner is the recursive child tail, not an independently weighted Einstein–Gauss–Bonnet or heat action. Geometry is already the compensating anomaly of the same fermionic determinant. Adding that geometry again is forbidden double-counting. Compact record: `results/nsc-2-zeta1-anomaly-owner-correction.json`, whose `nonclaims` include `recursive_tail_is_known_to_restore_a_scale_root: false`.

## Earlier diagnostics that must not be promoted

These numbers were obtained in this repository and are useful. Their JSON nonclaims forbid reading them as Nested-Space Cosmology, a final \(\zeta\), or an observed particle.

| Compact record | What it is | What it is not |
|---|---|---|
| `nsc-2-zeta1-lowest-mode.json` | Isolated lowest-mode heat minimum \(\zeta\approx 5.096657\) | Complete spectrum or physical \(\zeta\) |
| `nsc-2-zeta1-angular-tower.json` | Angular tower removes the heat-only root | Proof that no scale stationarity exists |
| `nsc-2-zeta1-regulated-determinant.json` | Determinant restores \(\zeta_{\det}\approx 4.990724\) | Final \(\zeta\) or physical gap |
| `nsc-2-zeta1-y-boundary-sensitivity.json` | Factorized \(y\) root tracks an even zero mode | Selection of a zero mode to force a root |
| `nsc-2-zeta1-orbifold-parity.json` | Orbifold parity selects one chiral zero mode; \(\zeta_{\mathrm{orb}}\approx 4.747863\) | Exact warped operator or promoted \(\zeta\) |
| `nsc-2-zeta1-warped-y.json` | Exact compact warp preserves \(\zeta\approx 4.748390\) | Full Dirac operator with lapse and shift |
| `nsc-1-s-one-child-orientation.json` | Child Ricci sign selects source orientation and the conditional gap \(1006/1015\) | Electron \(\Xi_{\mathrm{NSC}}\) or our universe |
| `nsc-1-s-one-child-scale-correction.json` | Exact child vacuum rejects \(\zeta=1\); requires \(\zeta>3\pi/2\) | A predicted \(\zeta\) or electron mass |

The child-scale correction also supersedes reading \(54/503\) as a child cosmological prediction.

## Other compact closures, with their nonclaims

The chain contains the original 58 v0.1.0 JSON records, preserved byte-for-byte, plus scoped follow-ups 59–64. The following are the principal closures and the statements they explicitly do not make.

**Imported bindings, not novelty**

- `nsc-1-exact-black-universe-defocusing.json` — imported regular black-universe geometry has a trapped positive-\(Q\) interval and child expansion. Nonclaims: our observed universe is this solution; phantom microscopic stability; Skyrme already coupled.
- `nsc-1-s-one-particle-wave-identity.json` — imported BPS collective reduction: one charge configuration and its wave share rest energy. Nonclaims: electron solution; spin-half and electric charge; Bell reduction; complete \(\Theta\).
- `nsc-1-gravitating-bps-observation-link.json` — nuclear-fitted Einstein–BPS maximum mass lands in the observed compact-object transition region. Nonclaims: precision observational fit of \(3.34\,M_\odot\); static star is a collapse trajectory; BPS stress already defocuses; parent/child energy ledger closed.
- `nsc-1-s-one-spectral-room-bind.json` — known local particles, forces, and gravity bound as spectral projections of one room operator. Nonclaims: Standard Model mass ratios predicted; finite Dirac spectrum derived; dark/black kernel completed.
- `nsc-1-s-one-outside-black-bind.json` and `nsc-1-s-one-shadow-matter-bind.json` — imported RS2 Weyl and two-wall shadow projections. The \(3/4\) lensing discriminator is imported, not discovered here.
- `nsc-1-s-one-shadow-lensing-observation.json` — constant two-wall shadow cannot dominate the tested kiloparsec gravitating mass. This does not reject recursive outside geometry as a whole.

**Repository identities**

- Nested \(L0\)–\(L6\) self-equality leaves one overall normalization (`nsc-1-s-one-nested-pair-closure.json`, `nsc-1-s-one-relevant-direction-count.json`).
- Local constants are named spectral outputs (`nsc-1-s-one-constant-dictionary.json`). Nonclaims: electron/proton/neutron spectrum solved; numerical constants of nature predicted; inside-black-hole observationally shown.
- Two-sheet carrier, Gauss–Bonnet coupling \(1015/144\), on-room vacuum-form \(108/1015\), and dark/black invariant written three ways.
- Scalar recursive outside kernel with unit-normalized continuum.
- Exponential heat profile as a heat-kernel diagnostic, not the complete physical anomaly; flat Euclidean graviton form factor without extra cutoff-disk poles; background-adjusted massless room graviton with residue \(12\).
- Local identity that one \(\Phi\) is mass gap, Schur self-energy, and relative-metric link.
- Global spectral foliation and computed finite radial Weyl/DtN maps; full Lorentzian throat maps and the physical matching field remain open.
- Unit-consistent additive-gap argument \(q_j=\lambda_j/\zeta+\mu^2\) and a finite, fixed-geometry sign theorem.
- Finite regulated variations and \(1/\Omega\) chain controls; ultrastatic frequency factor; nonelliptic trapped coordinate-time Hamiltonian.
- Isolated radial throat gapless by Weyl sequence.
- Controlled periodic-throat geometric gap with stencil-derived link. \(R\) and \(\Omega\) remain inputs.

**Invalidated truncations**

- Pure critical Einstein–Gauss–Bonnet loses half the local time-principal rank; induced Einstein restores rank but no scanned diagonal normalization restores real characteristics (`nsc-1-s-one-gauss-bonnet-kinetic-rank.json`, `nsc-1-s-one-induced-interface-rank.json`, `nsc-1-s-one-interface-characteristic-scan.json`). Nested-Space theory is not rejected by that scan; the omitted object is the off-diagonal \(\Phi\) block.
- Standalone bosonic heat-trace graviton covariance fails a necessary Stieltjes condition (`nsc-1-s-one-reflection-positivity-obstruction.json`).
- Inverse relative heat Hessian fails the finite Osterwalder–Schrader test (`nsc-1-s-one-relative-reflection-test.json`); the fermionic sheet-resolution observable does not (`nsc-1-s-one-fermionic-relative-observable.json`). Neither is the electron, and \(\sigma_3\) is not the complete graviton.

## Exact current-frontier statement

**The current frontier is the unsolved energy-resolved recursive child tail.**

Operationally: solve

$$
\Gamma_p(x)=K_p(x)-\frac{1}{\Omega}\,b\,\Gamma_c(x/\Omega)^{-1}b^{\dagger}
$$

on the warped parent/child mode spectrum, with consistent mass units and one regulated operator, then recompute the anomaly-compensated scale derivative. Do not restore a root by adding a second geometric weight, retuning \(\Phi\) after seeing the target, inserting a dark function, or refining the finite additive-gap family. The full covariant measure, physical stationarity, metric health, and nuclear or cosmological predictions remain open. The periodic-throat gap does not close those equations.

Until that tail is solved, public documents must keep every earlier \(\zeta\) candidate, \(1006/1015\), \(54/503\), and every particle or cosmological reading in the diagnostic column.

## Computed finite boundary response — local September 8 revision

The [first-order radial boundary calculation](nsc-boundary-response.md) now evaluates oriented parent/child Weyl maps with sparse LU and independent ODE integration. The parent is rho>0 and its throat normal is negative; the child is rho<0 and its normal is positive. Their jump is E(m_child-m_parent). These finite spatial maps supersede the earlier unevaluated declarations only in their stated domain; full Lorentzian horizon maps and the physical matching field remain open.

The [threshold example](nsc-threshold-response.md) independently gives a zero-energy DtN jump tending to 3/2 in inverse throat-length units while the radial bulk spectrum is gapless. Boundary stiffness is therefore not automatically a particle mass. The compact-mode warp factor is a projection with measured mixing, not an invariant-subspace theorem. Full covariant stationarity, physical metric health, nonlinear collapse and observational predictions remain unfulfilled.
