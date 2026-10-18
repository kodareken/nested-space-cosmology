# Cap mass and the displacement of the shared wall

The same leading field that supplies the exterior mass functional also
determines the cap's vacuum-subtracted mass and a linear equation for the
wall's quadratic displacement. Keeping the cap vacuum at the **actual**
radius is essential. The temporal and normal junction conditions then agree
with the angular condition and two-sided mass transport.

This is the next connection in the [atlas programme](nsc-cosmological-atlas.md).
It extends the [exterior mass result](nsc-atlas-wall-mass-closure.md) under the
same action, preparation and finite no-return domain. Its
[calculation admission](prior-art-and-open-claim.md#cap-mass-and-moving-wall-response-admission-15-october-2026)
fixes those inputs. The [implementation](../src/recursive_horizons/nsc_atlas_cap_wall_response.py)
and [independent checks](../tests/test_nsc_atlas_cap_wall_response.py) make the
identities executable. No sampled waveform is differentiated.

## One action, one proper clock, two mass coefficients

Use the original Jordan metric, signature $(-+++)$, a canonical real
massless scalar, and

$$
f(\phi)=\frac{1}{2\kappa_4}-\frac{\phi^2}{12},\qquad
S_{\rm bulk}=\int\sqrt{-g}\left[f(\phi)\mathcal R
-\frac12(\nabla\phi)^2-V_0\right]d^4x.
$$

Here $V_{0,-}=3H^2/\kappa_4$ in the cap and $V_{0,+}=0$ outside;
these are fixed, scalar-independent vacuum terms. The material action is
$-\sigma\int\sqrt{-\gamma}\,d^3y$, with the corresponding $2fK$
boundary terms. Set $k=\kappa_4\sigma/2$, $h=H\ell=1$, $b=k\ell=2$,
and $m=M/\ell=2/5$. The scalar background vanishes and $f>0$.
The [model specification](../src/recursive_horizons/nsc_cosmological_atlas.py)
records the field, preparation, parameter roles and observers separately.

Let $s=\tau/\ell$ be proper wall time, $x(s)$ its reference radius, and
$r(s)$ the corrected radius, both in units of $\ell$:

$$
\alpha=\frac{\kappa_4\varepsilon^2}{\ell^2},\qquad
r=x+\alpha\rho+O(\alpha^2),\qquad v=x_s>0.
$$

On $x\in[0.45,0.55]$ the signed cap-to-exterior normals satisfy

$$
\beta_+=\frac{m}{bx^2}-\frac{b^2+h^2}{2b}x,\qquad
\beta_-=\beta_++bx,
$$
$$
\beta_+^2-v^2=1-\frac{2m}{x},\qquad
\beta_-^2-v^2=1-h^2x^2.
$$

Define geometric mass by $\mathfrak m=R(1-\nabla R\cdot\nabla R)/2$.
The exterior and cap expansions are

$$
\frac{\mathfrak m_+}{\ell}=m+\alpha\mu_++O(\alpha^2),\qquad
\frac{\mathfrak m_-}{\ell}=\frac{h^2r^3}{2}
+\alpha\nu_-+O(\alpha^2).
$$

The vacuum-subtracted scalar has zero background, so its leading coefficient
$\nu_-$ is unchanged by a first-order change of coordinate identification.
The total cap coefficient pulled back to the reference wall is therefore
$\nu_-+3h^2x^2\rho/2$. Replacing the cap background by a constant
would remove that displacement term and change the matching problem.
The shared quiet corner fixes $\mu_+(s_a)=\nu_-(s_a)=\rho(s_a)=0$.

## Cap transport from the leading field

Write $w=\chi_W/\varepsilon$, $j=w_s$, $p=w/x$ and

$$
A=p_s=\frac{j-vp}{x},\qquad
B_\pm=N_\pm p=\frac{j-\beta_\pm p}{x},\qquad
B_+-B_-=bp.
$$

These are one-sided smooth derivatives in the no-return domain. The
original $C^3$ compact preparation and quiet cap are unchanged. The
Hessian below is covariant: $H_{uu}$ and $H_{uN}$ are not derivatives of
stored samples. On the cap the conformally coupled scalar equation and
full improved stress give

$$
H_{NN}=H_{uu}+\frac{2}{x}(vA-\beta_-B_-)+2h^2p,
$$
$$
q_-=t_{uN,-}=\frac23AB_--\frac13pH_{uN},\qquad
P_-=t_{NN,-}=\frac16A^2+\frac12B_-^2
-\frac13pH_{NN}+\frac16h^2p^2.
$$

Here $t_{ab}=\ell^4T^\phi_{ab}/\varepsilon^2$. Subtracting the exact
vacuum contribution to the geometric mass law leaves

$$
\nu_{-,s}=\frac{x^2}{2}(\beta_-q_--vP_-).
$$

It has a first-derivative integral representation:

$$
L_-=\frac{x^2}{2}\left[\beta_-AB_-
-\frac v2(A^2+B_-^2)\right]+\frac{h^2x^2vp^2}{2},
$$
$$
Q_-=\frac{x^2p}{6}(vA-\beta_-B_-)-\frac{h^2x^3p^2}{12}
=\frac{(v-\beta_-)wj}{6}
+\left(\frac{1}{6x}-\frac{h^2x}{4}\right)w^2,
$$
$$
\boxed{\nu_-(s)=\int_{s_a}^s L_-(s')\,ds'+Q_-(s).}
$$

The extra term in $L_-$ also follows by transforming the fixed cap vacuum
action. For $\Omega^2=1-\kappa_4\phi^2/6$,

$$
V_E=\frac{3H^2}{\kappa_4}\Omega^{-4}
=\frac{3H^2}{\kappa_4}+H^2\phi^2+O(\phi^4).
$$

It would be inconsistent to carry the exterior's zero-potential shortcut
into this cap. The vacuum background's conformal radius change also enters
$Q_-$. The frame calculation uses smooth bulk limits and leaves the
original material wall action in the Jordan frame.

There is a direct check without a frame transformation. Use the cap
identities

$$
x_s=v,\quad p_s=A,\quad
v_s=a_-\beta_-+h^2x,\quad \beta_{-,s}=a_-v,
$$
$$
A_s=H_{uu}+a_-B_-,\qquad B_{-,s}=H_{uN}+a_-A.
$$

Differentiation gives $L_-+Q_{-,s}=x^2(\beta_-q_--vP_-)/2$;
every $a_-$ term cancels. Omitting the vacuum term spoils the identity.
For the exterior, use the already checked $L_+,Q_+$ and $\mu_+$ from
the [previous owner](nsc-atlas-wall-mass-closure.md).

## The conformal junction and wall displacement

With $[X]=X_+-X_-$, the scalar-tensor trace junction for pure tension is

$$
-3\sigma=4(f+3f'^2)[K],\qquad
f+3f'^2=\frac1{2\kappa_4}.
$$

Consequently $[K]=-3k$, $[K_{ab}]=-k\gamma_{ab}$ and
$[n\phi]=k\phi$ on the $f>0$ branch. These follow from the same
scalar-tensor junction laws used in the existing
[surface calculation](nsc-atlas-scalar-wall-consistency.md), rather than
assuming that a scalar has no gravitational stress.

Put $D=\mu_+-\nu_-$. Linearizing the angular junction and the two
geometric mass definitions gives

$$
\delta\beta_\pm=\beta_{\pm,x}\rho+\frac{D}{bx^2},
$$
$$
\boxed{\rho_s=v_x\rho+
\frac{\mu_+/x+\beta_+D/(bx^2)}{v}.}
$$

The cap form replaces the numerator by $\nu_-/x+\beta_-D/(bx^2)$;
the two forms agree exactly. Equivalently, since $x_s=v$,

$$
\rho(s)=v(s)\int_{s_a}^{s}
\frac{\mu_+(q)/x(q)+\beta_+(q)D(q)/(b x(q)^2)}{v(q)^2}\,dq.
$$

The homogeneous solution is $C v(s)$, a relative proper-time origin
shift. The shared quiet event sets $C=0$. These formulas divide by $v$,
which is positive on the declared interval, and stay regular through
$\beta_+=0$. They specify a displacement functional; no numerical
displacement curve has yet been evaluated.

## Temporal and normal conditions from the same transport

The previously derived scalar junction identities give $q_+=q_-=q$
and $P_+=P_-=P$ at this order. Subtracting the mass laws yields

$$
D_s=-\frac{bx^3q}{2}.
$$

For a spherical sourced metric, write
$\mathcal T_{ab}=\kappa_4\ell^2T_{ab}$. The mixed Einstein equation
and geometric mass definition imply, on each side,

$$
\beta_{\pm,s}=a_\pm r_s-\frac r2\mathcal T_{uN,\pm},\qquad
a_\pm\beta_\pm=r_{ss}+\frac{\mathfrak m_\pm/\ell}{r^2}
+\frac r2\mathcal T_{NN,\pm}.
$$

Using $D_s$ and the displacement equation gives equal acceleration
corrections on the two sides,

$$
\delta a_+=\delta a_-=
\frac{6m\rho}{bx^4}-\frac{2D}{bx^3}.
$$

The temporal jump is preserved. The Gauss relation then gives

$$
\overline K=\overline a+\frac{2\overline\beta}{r}
=-\frac{[\mathcal T_{NN}]}{2b}
=-\frac{3h^2}{2b}+O(\alpha^2).
$$

Here $\overline K$ is the dimensionless mean trace. The fixed vacuum jump
is $3h^2$ and the quadratic scalar normal-stress jump vanishes. Thus the
combined metric and displacement response satisfies $\delta\overline K=0$.
This does not require $\rho=0$. At $\beta_+=0$, the exterior acceleration
equation's compatibility condition is precisely
$\mu_{+,s}=-x^2vP/2$, its regular mass law.

## What this connection supplies next

The boundary calculation now supplies both mass functionals, the
quiet-anchored displacement functional, and compatibility of the angular,
temporal and normal conditions. It requires no freely chosen extra wall
force or cap mass constant. First scalar derivatives suffice to evaluate
the mass and displacement functionals; Hessians enter the independent
transport proof, not their numerical evaluation.

The next dynamical calculation must reconstruct the smooth sourced bulk
metric on both sides with these wall data and the existing characteristic
initial data. In particular, mass does not alone determine the null lapse,
the wall's null-coordinate embedding, or the central observer's clock.
Those fields are needed to calculate the corrected arrival map and the
next scalar response. A nonzero allowed amplitude and an error bound
require control of this coupled solution's remainder on the finite domain.
Formation, longer-time stability and composition into a second child use
that evolved solution as their input.

The smooth bulk equations can be written in each side's Einstein frame as
$ds_E^2/\ell^2=-C\,dU\,dV+r_E^2d\Omega^2$. Set $h_+=0$ outside and
$h_-=1$ in the cap, and use the already determined leading field $p$.
At the retained order the equations to expand about the reference are

$$
r_{E,UU}-(\log C)_U r_{E,U}=-\frac\alpha2 r_E p_U^2,
\qquad
r_{E,VV}-(\log C)_V r_{E,V}=-\frac\alpha2 r_E p_V^2,
$$
$$
r_{E,UV}=-\frac{r_{E,U}r_{E,V}}{r_E}-\frac{C}{4r_E}
+\frac{Cr_E}{4}(3h_\pm^2+\alpha h_\pm^2p^2),
$$
$$
(\log C)_{UV}=\frac{2r_{E,U}r_{E,V}}{r_E^2}
+\frac{C}{2r_E^2}-\alpha p_Up_V.
$$

Products in these equations are expanded only through order $\alpha$.
For example, direct curvature calculation gives
$G_{UV}=2r_{E,UV}/r_E+2r_{E,U}r_{E,V}/r_E^2+C/(2r_E^2)$ and
$G_{\theta\theta}=-4r_Er_{E,UV}/C-2r_E^2(\log C)_{UV}/C$.
The potential contributes to the radius equation and cancels from the
mixed lapse equation. Before matching the material boundary, transform
back using $r_J=r_E(1+\alpha p^2/12)+O(\alpha^2)$ and
$\log C_J=\log C+\alpha p^2/6+O(\alpha^2)$. Impose the existing
source/corner normalization, quiet initial cap and regular-centre data,
and match the induced metric and wall embedding in the Jordan frame.

For the subsequent order-$\varepsilon\alpha$ scalar response, the
higher-order continuation of the preparation must also be specified.
Continuing $\chi(U_0,V)=\varepsilon g(V)$ exactly in the anchored source
gauge is one explicit choice. Holding $\phi$ instead gives the same leading
probe but a different next-order source. The present quadratic mass and
wall identities do not depend on that choice. In each smooth Jordan region,
the on-shell improved scalar trace is zero, so $\mathcal R=12H_\pm^2$
with fixed vacuum; the next field forcing comes from the changed wave
operator and moving boundary, not a separately adjustable curvature term.

The [record](../results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json)
and [manifest](../results/development/nsc-atlas-cap-wall-response-v1/manifest.json)
bind the equations, implementation, checks and unchanged exterior evidence.
The executable algebra verifies the identities under their stated
assumptions; numerical integration and bulk evolution remain separate
calculations. Source citations are the already used Aviles--Maeda--Martinez
[1910.07534v2](https://arxiv.org/html/1910.07534v2), Hayward
[gr-qc/9710089v2](https://arxiv.org/html/gr-qc/9710089v2), and
Faraoni--Vitagliano [1401.1189v1](https://arxiv.org/html/1401.1189v1).
