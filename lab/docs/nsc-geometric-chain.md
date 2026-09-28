# A constructive intervention: repeated geometry creates a spectral gap

The isolated spatial throat is gapless, as the [Weyl-sequence proof](nsc-radial-spectrum.md)
shows. This experiment changes one geometrical condition: repeat the finite
segment `rho in [-R,R]`, with `w=1/sqrt(1+rho²)`, instead of continuing into two
asymptotically widening ends. The [calculation](../scripts/check_nsc_geometric_chain.py)
and [record](../results/nsc-3-geometric-chain.json) test this distinct realization.

No independent constant mass Phi is added. R is a declared intervention
parameter, not a derived scale. The repeated sphere-radius profile is bounded;
its derivatives change at the seams. The calculation solves a spatial Dirac
problem, not the gravitational junction equations of a cosmological geometry.

## The gap follows from the same first-order operator

At zero energy the exact transfer through one motif has multipliers

\[
\exp(\pm I),\qquad I=\int_{-R}^{R}w\,d\rho=2\operatorname{asinh}R>0.
\]

Neither multiplier has modulus one. Thus zero is absent from every real
Bloch-phase fiber. Continuity and the compact Bloch-phase interval give an
open spectral gap. This is a periodic spatial confinement result, rather than
the massive spectrum of the elementary constant-Phi symbol.

The continuum calculation integrates the first-order 2-by-2 transfer matrix
and brackets the first crossing of its trace through 2. An independent
node/edge lattice diagonalizes the Bloch operators at 33 phases, using 32, 64
and 128 intervals per motif. The observed errors converge at second order.

| R in throat units | Continuum first band edge | Finest lattice discrepancy |
|---|---:|---:|
| 2 | 0.7034881641 | 0.00000688 |
| 4 | 0.4552655377 | 0.00000620 |
| 8 | 0.2479700183 | 0.00000243 |

The zero-energy gap exclusion is analytic; the reported band-edge values and
the first-crossing identification have numerical resolution limits. Removing
w closes the zero-phase gap in the control. That formal kappa=0 control is not
an allowed spinor-sphere eigenvalue; it isolates the angular geometric term.

## The link is computed from the discretization

For one motif the staggered derivative is
`A u_i=(-1/h+w_i/2)u_i+(1/h+w_i/2)u_(i+1)`.
Grouping its degrees of freedom by motif derives the next-room link:

\[
B_{\mathrm{last\ edge},\mathrm{first\ node}}=1/h+w_{\mathrm{last}}/2.
\]

All other cross-motif entries vanish. This B is a geometric stencil entry,
not a fitted particle mass. The normalized recursion with 1/Omega agrees with
direct multi-motif inversion at depths 2, 4, 8 and 16 for Omega=1 and 1.3.
Only the equal-scale infinite periodic family has the Bloch gap proof; the
unequal-scale cases are finite response controls.

This closes a useful mathematical connection: the same repeated spatial
operator supplies both the band gap and the eliminated outside response.
It does not determine R or Omega, identify a particle species, establish
Lorentzian nesting, or solve stationary gravity. The earlier finite-regulator
common-scale monotonicity is not removed by this positive gap result.

## Originality and next physical question

Floquet theory and geometric confinement are established mathematics; this
experiment applies them to the project's throat profile. No new-to-world
mechanism is claimed. [Teschl's primary mathematical treatment](https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode)
provides the Floquet framework used here.

The next physical question is whether the common action dynamically produces
the required repeated/confining structure and fixes its spacing while retaining
the parent and child causal domains. Until it does, this is a controlled
positive realization with an explicit remaining geometrical input.

Reproduce: `python3 scripts/check_nsc_geometric_chain.py --check`.
