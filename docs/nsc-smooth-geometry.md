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


## Relative vacuum stress on the same compact axial circle

For the compact axial $S^1$ control, changing periodic to antiperiodic
spin structure changes the vacuum stress on this actual varying radius.
With $ds^2=N^2dt^2-q^2dx^2-r^2d\Omega^2$, half-density normalization
$u=r\sqrt q\,\psi$ gives

$$H_\kappa=-i\sigma_2(c\partial_x+c'/2)+v\sigma_1,
\qquad c=N/q,\quad v=N\kappa/r.$$

For one complex four-component massless Dirac field,
$E_\eta=-4\sum_{\kappa\ge1}\kappa\sum_js_j(A_{\kappa,\eta})$,
including both sphere signs and the particle/antiparticle determinant
count once. Subtract each angular sector before summing. Varying every
$N_i,q_i,r_i$ by the polar-factor Hellmann--Feynman identity, including
edge-coefficient variation, gives the relative stress for
$\Delta E=E_P-E_{AP}$:

$$\Delta\rho=\frac{g_N}{4\pi qr^2h},\quad
\Delta p_x=-\frac{g_q}{4\pi Nr^2h},\quad
\Delta p_\perp=-\frac{g_r}{8\pi Nqrh}.$$

The pointwise Weyl identity $Ng_N+qg_q+rg_r=0$ makes the relative trace
zero; equal local anomalies/counterterms cancel between domains, so this
does not assert a zero absolute trace. Lapse normalization and conservation
are checked independently. At512 points and $\kappa_{\max}=8$:

| R | $\Delta E$ | Throat $\Delta(\rho+p_x)$ | $\int\Delta(\rho+p_x)dx$ | RMS conservation residual |
|---:|---:|---:|---:|---:|
| 2 | 0.1586544847 | 0.0186082412 | 0.0409127521 | $4.44\times10^{-7}$ |
| 4 | 0.0194597880 | 0.0017988910 | 0.0041368012 | $1.53\times10^{-7}$ |

Positive P-minus-AP null differences mean AP lowers the source in the
needed direction. Axial64/128/256/512 controls showed second-order
movement, with the conservation residual decreasing roughly fourfold per
doubling; high-angular subtractive SVD roundoff is not a proved tail sign
or bound. Absolute AP stress, covariant determinant normalization,
counterterms and the same-action metric equations remain required.
No G or subtraction coefficient was fitted to the imposed Einstein stress.
This compact axial circle is distinct from an unwrapped periodic array
and the transverse compact carrier.

## Geometric variation on the actual smooth transmitting surface

The separate shared PG surface at $\rho=0$ has one smooth metric jet and
opposite outward normals. In the future-child convention
$a=\sqrt{\beta^2-1}$ and
$n=(\beta/a)\partial_\tau-a\partial_\rho$, with $n_p=n$, $n_c=-n$.
Only on stationary background scalars is $D_Tf=-a\,df/d\rho$; the full
normal derivative must remain on time-dependent fields. Matching intrinsic
variations gives $\delta a_c=\delta a_p$, $\delta r_c=\delta r_p$,
while outward shear variations obey $\delta\sigma_c=-\delta\sigma_p$.
The Weyl metric covectors $\overline P_a,\overline P_r$ are odd under
normal reversal, and $\Pi_\sigma$ is even. Consequently

$$(\overline P_{a,p}+\overline P_{a,c})\delta a+
(\overline P_{r,p}+\overline P_{r,c})\delta r+
(\Pi_{\sigma,p}-\Pi_{\sigma,c})\delta\sigma=0.$$

The actual reference values are
$a_0=\sqrt{3\pi/2-1}$, $r_0=1$, $D_Ta=3$, $D_Tr=0$,
$D_T^2r=a_0^2$, $F=2$, $D_TF=-6a_0$. For unchanged
$C_W=-0.0008441130342798952$, parent values are
$\overline P_a=0.4150924626502038$,
$\overline P_r=-0.6300629292752326$,
$\Pi_\sigma=0.10900236876297109$ and
$\sigma=1.557021169290531$; the child has the oriented counterparts.
Einstein--GHY momenta and the already included odd Euler/BoxR endpoint
primitives cancel as well. This is geometric seam matching in the shared
smooth-jet domain, with tolerance $3\times10^{-11}$, not a thin-shell
result, the old93.54264532195464 end-node diagnostic or full stationarity.
The Gaussian Schur factor and complementary bulk are already counted
once. Remaining field variations and bulk Einstein--quantum source
matching stay OPEN.
