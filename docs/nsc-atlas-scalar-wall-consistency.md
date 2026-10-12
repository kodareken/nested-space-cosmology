# Conformal probe stress at the existing wall

The [prepared response](nsc-atlas-prepared-response.md) propagates a linear
classical scalar through a fixed spherical wall. Its normalized waveform did
not need an absolute field normalization. Its gravitational stress does.
This note fixes that normalization and checks the first local interface
conditions at quadratic order in the field amplitude.

The result is a compatibility statement: the transmitted conformal probe
adds no singular surface stress and its regular tangential flux has no jump.
The constant material tension therefore needs no extra scalar surface source
for these two conditions. The smooth bulk stress and the metric response to
it remain separate work. This is not a self-consistent gravitating solution,
a small-backreaction estimate, or a stability result.

The [admission](prior-art-and-open-claim.md#conformal-probe-stress-and-shell-admission-13-october-2026)
was recorded before calculation. The [literature card](arxiv-studies.md#conformal-stress-library)
records the imported junction laws and reading boundaries. The analytic
record is [consistency.json](../results/development/nsc-atlas-scalar-wall-consistency-v1/consistency.json);
its [producer](../scripts/record_nsc_atlas_scalar_wall_consistency.py) checks
exact symbolic identities and authenticates the unchanged prerequisites.

## Fixed action, domain and order

Use the same four-dimensional Jordan-frame metric, signature $(-+++)$,
continuous induced metric $\gamma_{ab}$, and unit spacelike normal from the
de Sitter cap to the Schwarzschild side. All jumps use this one normal.
Restore $G$ only to keep the gravitational and scalar normalizations clear:

$$
\kappa_4=8\pi G,\qquad k=4\pi G\sigma=\frac{\kappa_4\sigma}{2},
\qquad [K_{ab}]=-k\gamma_{ab},\qquad [K]=-3k.
$$

Here $K_{ab}=\frac12\mathcal L_n\gamma_{ab}$ and $\sigma>0$ is the fixed
material tension. The recorded geometry uses $G=1$. Its dimensionless
parameters remain $H\ell=1$, $k\ell=2$, $M/\ell=0.4$. The field domain is
the v7/v8 retained strip and its finite central arrival domain, with wall
$R/\ell\in[0.45,0.55]$, exterior floor $R/\ell=0.4$, and the same zero
earlier cap data. No data are supplied by the causally excluded right-parent
observer.

Fix a canonical positive kinetic coefficient, rather than leaving its
overall multiplier arbitrary. In each smooth bulk region take

$$
S_{\rm bulk}=\int\sqrt{-g}\left[
 f(\phi)\mathcal R-\frac12(\nabla\phi)^2+
 \mathcal L_{\rm vac}\right]d^4x,\qquad
f(\phi)=\frac{1}{2\kappa_4}-\frac{\xi}{2}\phi^2,
\quad \xi=\frac16.
$$

The scalar potential is zero, the background vacuum terms are independent
of $\phi$, and the material surface action is
$-\sigma\int_\Sigma\sqrt{-\gamma}\,d^3y$.
The split action includes the $2fK$ boundary terms of the borrowed
variational junction construction. There is no added scalar surface mass,
kinetic term or source. Multiplying the background vacuum term by $f$ would
change the scalar equation and is not done here.

There is no background scalar. Write $\phi=\varepsilon\Phi$ for the same
small-amplitude family used for the prepared pulse. The stress is quadratic
in $\varepsilon$ and is evaluated on the existing geometry; the order
$\varepsilon^2$ metric correction has not been solved. Assume $f>0$ in a
weak-field interpretation. This assumption alone bounds neither derivatives
nor the size of the gravitational response.

The proof uses continuous $\phi$ with sufficiently differentiable one-sided
traces. The original $C^3$ compact preparation supplies the required
piecewise second derivatives on the finite regular domain, including the
regular-centre construction. The broader $C^2$ reduced-field inverse class
does not by itself establish this four-dimensional stress regularity for
every source. No derivative is estimated from sampled waveform data.

## Borrowed stress and a distributional calculation

Put the nonminimal Einstein-tensor term on the scalar side, exactly once:

$$
T^\phi_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
 -\frac12g_{\mu\nu}(\nabla\phi)^2
 +\xi\left[G_{\mu\nu}\phi^2+
 (g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2\right].
$$

This follows from the canonical nonminimal action; equivalently one can
move $\xi G_{\mu\nu}\phi^2$ to the gravitational side. These are two
accounting conventions for the same equations, not two sources to add.

Let $d=[n\phi]$ and
$C_{ab}=[K_{ab}]-\gamma_{ab}[K]=2k\gamma_{ab}$.
The existing scalar equation gives $d=-2\xi\phi[K]=6\xi k\phi$, hence
$d=k\phi$ for the selected conformal coupling. In Gaussian normal
distance, the singular Einstein tensor is $-C_{ab}\delta_\Sigma$ and
$[n(\phi^2)]=2\phi d$. The first-derivative kinetic terms are only
piecewise regular. The full singular scalar coefficient is therefore

$$
S^\phi_{ab}=\xi\{2\phi d\,\gamma_{ab}-\phi^2 C_{ab}\}
            =2\xi\phi(d-k\phi)\gamma_{ab}=0.
$$

It has no normal-normal or normal-tangential delta coefficient. Products
of a delta with a discontinuous normal derivative are never used:
$\phi^2$ is continuous and the second derivatives are taken
distributionally. Vanishing $S^\phi_{ab}$ does not mean the regular scalar
stress vanishes at or away from the wall.

## Independent check using the coupled junction law

Aviles, Maeda and Martinez [1910.07534v2](https://arxiv.org/html/1910.07534v2),
Eqs. (3.43)-(3.46), give the scalar-independent material surface tensor
for a timelike interface:

$$
t_{ab}=-2fC_{ab}+2d f'(\phi)\gamma_{ab},\qquad
d=2f'(\phi)[K].
$$

Substitution of the declared action and background jumps yields

$$
t_{ab}=\left[-\frac{2k}{\kappa_4}
       +2\xi\phi(k\phi-d)\right]\gamma_{ab}
       =-\sigma\gamma_{ab}.
$$

The scalar correction is precisely the negative of the independently
computed $S^\phi_{ab}$. This agreement checks the sign and the stress
accounting without solving for a new trajectory. The earlier
Barcelo-Visser junction laws supply the same coupled structure after
converting their opposite extrinsic-curvature convention; their different
bulk vacuum coupling is not imported.

## Regular tangential energy flux

The surface coefficient alone does not check energy exchange. Let
$p_\pm=(n\phi)_\pm$ and $v_a=D_a\phi$; $v_a$ is the common tangential
derivative. For one-sided regular traces in Gaussian normal coordinates,

$$
\nabla_n\nabla_a(\phi^2)
 =D_a\{n(\phi^2)\}-K_a{}^bD_b(\phi^2).
$$

Thus the full, improved scalar flux is

$$
T^\phi_{na}=(1-2\xi)p v_a-2\xi\phi D_a p
             +2\xi\phi K_a{}^b v_b+\xi\phi^2G_{na}.
$$

Both unperturbed bulk vacuum regions have $G_{na}=0$. With constant $k$,
$[p]=k\phi$, $D_a[p]=k v_a$ and $[K_a{}^b]=-k\delta_a{}^b$, its jump is

$$
[T^\phi_{na}]=(1-6\xi)k\phi v_a=0.
$$

This is compatible with the unchanged constant-tension tangential
conservation condition. The canonical kinetic term alone would instead
give $k\phi v_a$ and suggest an exchange that the full stress cancels.
Zero jump means equality of the two one-sided fluxes, not zero flux on
either side and not zero transmitted energy.

## Checks, limits and the next unlocked calculation

The [tests](../tests/test_nsc_atlas_scalar_wall_consistency.py) compare
exact symbolic expressions for the distributional and variational
calculations, the one-sided flux, and conservation. With a symbolic
coupling, the inherited scalar jump is $d=6\xi k\phi$; the surface and
flux coefficients become $2\xi(6\xi-1)k\phi^2$ and
$4\xi(1-6\xi)k\phi v_a$. They obey the tangential conservation identity
and vanish at the selected coupling. This algebraic control is not an
evolved alternative-coupling model. Altering a jump or omitting an
improvement term is detected. Zero amplitude, a common positive action
normalization, and simultaneous reversal of the common normal and exchange
of the side labels are checked. One-sided normal quantities then negate
and swap, leaving the jump statements and positive material tension
unchanged. Reversing a normal while retaining every old oriented formula
would not be a consistent convention change.

These exact identities have no quadrature or floating-point error. Their
scope depends on the stated regularity, vacuum mixed components, constant
tension and fixed coupling. Symbolic simplification and evidence hashes
check implementation and provenance; the independent analytic review
checks the derivation. None supplies a formal proof of the entire spacetime.

The smooth bulk tensor is still an order-$\varepsilon^2$ source in
$G_{\mu\nu}=\kappa_4(T^{\rm vac}_{\mu\nu}+T^\phi_{\mu\nu})$.
This note neither evaluates that tensor throughout the pulse nor solves
the corresponding metric and wall corrections. Normal-force balance,
initial gravitational constraints, finite-amplitude validity and stability
remain unverified at that order. A vanishing surface coefficient cannot
stand in for those equations.

The useful next connection is a controlled source or backreaction criterion:
the fixed preparation and physical normalization can now be used to ask
what stress the pulse supplies in a smooth bulk region and which geometric
correction that requires. A later admission must identify that quantity,
the regularity needed to evaluate it, and a finite budget. It must not assume
pointwise positivity for nonminimal scalar stress or replace it by a
perfect fluid. Physical source selection, formation, quantum states and
other attachment routes remain open.
