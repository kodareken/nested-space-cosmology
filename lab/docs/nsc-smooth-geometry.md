# Seam-free periodic radius: gap without derivative seams, and the required Einstein residual

The [repeated-throat chain](nsc-geometric-chain.md) opened a spatial band gap by
periodizing `r=sqrt(1+rho^2)` on `[-R,R]`. That motif is continuous, but `r'`
jumps at every seam, so it is not a gravitationally matched periodic geometry.
This calculation changes only that geometric condition. The
[runner](../scripts/check_nsc_smooth_geometry.py) and
[record](../results/nsc-4-smooth-geometry.json) keep `a` and `R` as explicit
inputs and do not edit the published NSC-3 sources.

## Smooth profile

\[
r(x)=\sqrt{a^{2}+\Bigl(\frac{\sin(kx)}{k}\Bigr)^{2}},
\qquad
k=\frac{\pi}{2R},\qquad x\in[-R,R].
\]

For `a=1`, `r_min=r(0)=1` and `r''(0)=1`. The endpoints satisfy `r'(±R)=0` and
`r'''(±R)=0`, so the `2R`-periodic extension is `C^∞`. No independent mass
`Phi` is inserted. The spatial Dirac superpotential remains `w=kappa/r` with
`kappa=1`.

This removes derivative seams. It is not a stationary Einstein solution: the
metric is imposed, not varied to an extremum of `S_one` or of the Einstein-Hilbert
action.

## The gap is still a Floquet fact of the same first-order operator

At zero energy the monodromy multipliers are `exp(±I)` with

\[
I=\int_{-R}^{R}\frac{\kappa}{r}\,dx=\frac{2\kappa}{ak}\,K\Bigl(-\frac{1}{a^{2}k^{2}}\Bigr)>0.
\]

Zero is absent from every real Bloch fiber, so a gap exists. Continuum transfer
matrices bracket the first crossing of the discriminant through `2`. An
independent node/edge staggered lattice, the same structure-preserving stencil
family as the seamed chain, is diagonalized at 33 Bloch phases on 32, 64 and
128 intervals.

| R in throat units | Continuum first band edge | Finest lattice discrepancy | `8πG ∫(ρ+p_r) dx` |
|---|---:|---:|---:|
| 2 | 0.7642320355 | ~1.16e-6 | -0.5839376400 |
| 4 | 0.5142554634 | ~1.26e-7 | -1.3586948086 |
| 8 | 0.2845408508 | ~2.47e-6 | -2.0867274728 |

Lattice errors decrease at second order. Removing `w` closes the zero-phase
gap in the formal `kappa=0` control. Compared with the seamed `sqrt(1+rho^2)`
motif at the same `R`, the smooth profile has a slightly larger gap because its
areal radius stays smaller.

## Required 4D ultrastatic Einstein stress

For `ds²=dt²-dx²-r(x)² dΩ²` in signature `(+,-,-,-)`, with the Riemann
convention

\[
R^{\rho}{}_{\sigma\mu\nu}
=\partial_{\mu}\Gamma^{\rho}{}_{\nu\sigma}
-\partial_{\nu}\Gamma^{\rho}{}_{\mu\sigma}
+\Gamma^{\rho}{}_{\mu\lambda}\Gamma^{\lambda}{}_{\nu\sigma}
-\Gamma^{\rho}{}_{\nu\lambda}\Gamma^{\lambda}{}_{\mu\sigma},
\qquad
R_{\sigma\nu}=R^{\rho}{}_{\sigma\rho\nu},
\]

direct Christoffel/Ricci evaluation gives

\[
8\pi G\rho=\frac{1-r'^{2}-2rr''}{r^{2}},
\qquad
8\pi G p_{r}=\frac{r'^{2}-1}{r^{2}},
\qquad
8\pi G p_{t}=\frac{r''}{r}.
\]

The same `ρ` is the ADM Hamiltonian density `³R/2` of the spatial metric
`dx²+r² dΩ²`, and it is the lapse-normal projection `G_{μν}n^{μ}n^{ν}` after
`g_{tt}` is restored as `N(x)²`. The Euler-Lagrange derivative of
`∫ N r² ⁴R\,dx` equals `-r² ³R`. A metric finite-difference Ricci check on the
actual profile recovers the same stress.

These identities imply

\[
p_r'+2\frac{r'}{r}(p_r-p_t)=0
\]

and, because `r'(±R)=0`,

\[
8\pi G\int_{-R}^{R}(\rho+p_r)\,dx
=-2\int_{-R}^{R}\Bigl(\frac{r'}{r}\Bigr)^{2}dx<0.
\]

A cosmological constant has `p=-ρ`, so its radial-null projection vanishes. It
cannot supply this integral. The throat value is `8πG(ρ+p_r)|_{x=0}=-2/a²`.

This is the **required effective Einstein stress of the imposed 4D metric**. It
is not the derived matter stress of `S_one`, of the fermionic spectral action,
or of the 5D Gauss-Bonnet carrier. No phantom field and no fitted coupling were
added to produce the gap or the null integral: both are simultaneous
consequences of the same `r(x)`.

## Originality and next physical equation

Periodic Dirac gaps and wormhole flare-out stresses are established
mathematics. The project-level result is that the seam obstruction of the
previous chain can be removed while the gap survives, and that the same
smooth profile then forces a strictly negative radial-null Einstein residual
which a cosmological constant cannot cancel.

The next equation must identify that residual with the variational stress of
the same frozen `Θ` one-action / spectral-boundary functional:

\[
G_{\mu\nu}[r]
=8\pi G\,T^{\mathrm{rest}}_{\mu\nu}[\Psi,\mathbb D_{\Theta},\Phi_{\mathrm{throat}}],
\]

or, as a radial-null constraint that does not introduce a new coupling,

\[
\int_{-R}^{R}\ell^{\mu}\ell^{\nu}
\bigl(G_{\mu\nu}[r]-8\pi G\,T^{\mathrm{rest}}_{\mu\nu}\bigr)\,dx=0,
\qquad
8\pi G\int(\rho+p_r)_{\mathrm{one}}\,dx
=-2\int\Bigl(\frac{r'}{r}\Bigr)^{2}dx.
\]

Until that stress is computed, `a` and `R` remain inputs and the geometry is
not a stationary nested cosmology.

Reproduce: `python3 scripts/check_nsc_smooth_geometry.py --check`.

Here T_rest means the metric variation of Gamma_one after isolating the Einstein term already contained in that same action. Its coefficient determines G; this is bookkeeping, not adding an independent Einstein action or counting it twice. The null contraction uses the same null vector in both indices. The old seam diagnostic records the magnitude of the derivative jump, whose right-minus-left sign was negative.
