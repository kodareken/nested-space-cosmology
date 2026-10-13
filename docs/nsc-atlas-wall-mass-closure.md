# Exterior geometric mass from the existing wall trace

The leading field already available on the wall determines a conditional
quadratic geometric-mass coefficient on its smooth exterior side. Both null
mass fluxes contribute along this timelike path. Their integral and the full
conformal-frame endpoint term close the answer without solving another bulk
field equation or differentiating sampled waveforms. An independent calculation
with the original improved stress gives the same expression.

This closes an input to later wall matching. It does not solve the sourced
wall, cap or observer geometry. The [admission](prior-art-and-open-claim.md#exterior-wall-mass-closure-admission-13-october-2026)
precedes the calculation; the [record](../results/development/nsc-atlas-wall-mass-closure-v1/wall-mass.json)
and [manifest](../results/development/nsc-atlas-wall-mass-closure-v1/manifest.json)
bind its exact checks and prerequisites. No numerical mass curve is evaluated.

## Borrowed laws and the actual missing connection

An--Lim [2005.04090v1](https://arxiv.org/html/2005.04090v1), (2.14)--(2.15),
supplies both minimal-scalar null mass laws. Hayward
[gr-qc/9710089v2](https://arxiv.org/html/gr-qc/9710089v2), (2.9d) and (8.9a),
already supplies general-matter mass transport and its timelike projection.
Faraoni--Vitagliano [1401.1189v1](https://arxiv.org/html/1401.1189v1), (7),
gives the complete conformal transformation of geometric mass. The smooth
vacuum field/frame map is from Aviles--Maeda--Martinez
[1910.07534v2](https://arxiv.org/html/1910.07534v2), (2.21)--(2.24).

Scalar-induced changes in shell-adjacent mass are also prior art:
Danielsson--Lehner--Pretorius [2109.09814v1](https://arxiv.org/html/2109.09814v1),
Section VI.2 and (98), evolves a minimal scalar through a different shell
with an AdS interior and multicomponent material. Neither their dynamics nor
their mass amplification is imported into this conformal, de Sitter-cap,
constant-tension model. The question here is whether this atlas's fixed
traces provide every input to the borrowed mass laws. It is not a new general
transport law or a priority claim about scalar backreaction.

## Fixed data, order and the two scalar derivatives

Retain the canonical real massless scalar at $\xi=1/6$, zero scalar background
and potential, original $C^3$ preparation, zero bottom/cap data, fixed
vacuum and wall parameters, and the original no-return window. The wall has
$R/\ell\in[0.45,0.55]$, $H\ell=1$, $k\ell=2$ and $m=M/\ell=0.4$.
The preparation is in retained white-hole sector IV; the right-parent
observer remains causally excluded. All quantities below are smooth
exterior limits. No field/frame map is carried through the material action.

Use dimensionless proper wall time $s=\tau/\ell$, background radius $x=R/\ell$
and the outward radial normal $\beta=\beta_+$, with

$$
v=\frac{dx}{ds}>0,\qquad
\beta^2-v^2=f=1-\frac{2m}{x},\qquad
\beta=\frac{m}{(k\ell)x^2}-\frac{(k\ell)^2+(H\ell)^2}{2k\ell}x.
$$

Write the geometric mass length as

$$
\mathfrak m=\frac{R}{2}\left(1-\nabla R\cdot\nabla R\right),\qquad
\frac{\mathfrak m_J}{\ell}=m+\alpha\mu_J+O(\alpha^2),\qquad
\alpha=\frac{\kappa_4\varepsilon^2}{\ell^2},\quad \kappa_4=8\pi G.
$$

The subscript $J$ denotes the original metric. This is a formal coefficient;
no nonzero amplitude or bound on the omitted terms is chosen. Assume a
perturbative family with the fixed leading preparation and a regular smooth
exterior. The coefficient is a necessary response of such a family, not an
existence proof for its coupled material boundary problem.

The existing continuum wall trace supplies

$$
w=\frac{\chi_W}{\varepsilon},\qquad
j=\frac{dw}{ds}=\ell\dot V\frac{J}{\varepsilon},\qquad
p=\frac{w}{x}=\frac{\ell\phi_W}{\varepsilon}.
$$

Dots on $V$ denote physical proper-time derivatives. The inherited normal
is $n=(-\dot U,\dot V)$; $\chi_U|_W=0$ implies $n\chi=\dot\chi$.
For the dimensionless unit normal $N=\ell n$, both scalar derivatives are
therefore known:

$$
A=\frac{dp}{ds}=\frac{j-vp}{x},\qquad
B=Np=\frac{j-\beta p}{x}.
$$

In particular $\chi_U=0$ does not mean $\phi_U=0$: differentiating
$\phi=\chi/R$ retains the radius derivative. The original quiet bottom
anchors $\mu_J(s_a)=0$ at $x_a=0.45$, with $w=j=0$ there. No independent
outgoing scalar trace or new mass integration constant is supplied.

## Both null fluxes and the full frame correction

In the smooth exterior, $\bar g=\omega^2g$ with
$\omega^2=1-\kappa_4\phi^2/6=1-\alpha p^2/6$. The leading Einstein-frame
scalar derivatives agree with the original ones. In dimensionless null
coordinates $ds_4^2/\ell^2=-C\,dU\,dV+x^2d\Omega^2$, the two minimal
mass laws pulled back to the wall give

$$
\frac{d\mu_E}{ds}
=-\frac{x^2}{C}\left(x_Vp_U^2U_s+x_Up_V^2V_s\right)
=\frac{x^2}{2}\left[\beta AB-\frac{v}{2}(A^2+B^2)\right]
\equiv L(s).
$$

The normalization is $CU_sV_s=1$; the radial and scalar null components
obey $x_UU_s=(v-\beta)/2$, $x_VV_s=(v+\beta)/2$,
$p_UU_s=(A-B)/2$ and $p_VV_s=(A+B)/2$. No one-sided flux can be
discarded merely because $\chi_U=0$. The reference lies inside $x<2m$,
so $v>|\beta|$; the minimal coefficient is nonincreasing along this
future timelike path. This is not a sign statement about the physical
Jordan coefficient or radiation reaching infinity.

Expand the borrowed geometric transformation

$$
\bar{\mathfrak m}=\omega\mathfrak m
-\frac{R^3}{2\omega}(\nabla\omega)^2
-R^2\nabla\omega\cdot\nabla R.
$$

Its gradient-square term first enters at order $\alpha^2$. Keeping both
first-gradient components at the retained order yields

$$
\mu_J=\mu_E+Q,\qquad
Q=\frac{mp^2}{12}+\frac{x^2p}{6}(vA-\beta B)
=\frac{(v-\beta)wj}{6}
+\left(\frac{1}{6x}-\frac{m}{4x^2}\right)w^2.
$$

Consequently the anchored answer is the finite functional

$$
\boxed{\displaystyle
\mu_J(s)=\int_{s_a}^{s}L(s')\,ds'+Q(s),\qquad Q(s_a)=0.}
$$

Every integrand and endpoint input is fixed by $w,j$ and the reference
geometry. First scalar derivatives suffice for evaluation. The quiet
anchor is essential; a different mass datum would add a constant.

## Independent original-frame stress check

Use a radial orthonormal frame $(u,N)$ along the background wall. For the
dimensionless scalar Hessian write $H_{uu}$, $H_{uN}$ and $H_{NN}$;
these are covariant components, not second coordinate derivatives.
The Schwarzschild reference is Ricci-flat. With
$t_{ab}=\ell^4T^\phi_{ab}/\varepsilon^2$ in the physical orthonormal frame,
the full improved stress gives

$$
t_{uN}=\frac{2AB}{3}-\frac{pH_{uN}}{3},\qquad
t_{NN}=\frac{A^2}{6}+\frac{B^2}{2}-\frac{pH_{NN}}{3}.
$$

The scalar equation supplies
$H_{NN}=H_{uu}+2(vA-\beta B)/x$. Hayward's general-matter law then reads

$$
\frac{d\mu_J}{ds}=\frac{x^2}{2}(\beta t_{uN}-v t_{NN}).
$$

Let $a_n$ be the dimensionless signed acceleration along $N$. On the
Schwarzschild reference,

$$
x_s=v,\quad p_s=A,\quad
v_s=a_n\beta-\frac{m}{x^2},\quad \beta_s=a_nv,\quad
A_s=H_{uu}+a_nB,\quad B_s=H_{uN}+a_nA.
$$

Differentiating $Q$ with these geometric identities cancels every acceleration
term. The remaining expression gives exactly

$$
L+Q_s=\frac{x^2}{2}(\beta t_{uN}-v t_{NN}).
$$

This verifies the endpoint correction using the original stress and mass
law, without evaluating a second derivative of a stored waveform. Dropping
the mixed improved stress, a null flux or the frame gradient spoils the
comparison. The proof uses ordinary one-sided smooth tensors, not products
of distributions at the wall.

## Invariance, regularity and limits of the result

Neither $L$ nor $Q$ divides by $\beta$. Their expressions stay finite at
the reference $\beta=0$ crossing; no altered wall branch is needed. Null
reparametrizations cancel between coordinate derivatives, velocities and
the metric coefficient. The scalar geometric mass is coordinate invariant.
Because the background exterior mass $m$ is constant, an order-$\alpha$
change in the identification of wall cuts does not change $\mu_J$ at this
order. It does not follow that the wall displacement itself vanishes.

At fixed dimensionless geometry and preparation, changing the regional
length unit leaves $w,j,p,A,B,L,Q$ and $\mu_J(s)$ unchanged. The physical
mass-length correction is $\alpha\ell\mu_J$. Geometric mass is not
conformally invariant: the endpoint $Q$ is part of the answer, not an error.
Unlike the original source leg, the wall has a nonzero post-input tail.
No terminal frame agreement, constant mass plateau or late Jordan mass
sign is inferred from the compact support of the exterior input.

The [exact recorder](../scripts/record_nsc_atlas_wall_mass_closure.py) and
[tests](../tests/test_nsc_atlas_wall_mass_closure.py) check both flux routes,
the frame transformation, signs, dimensional and coordinate behavior and
omitted-term controls. The original source-mass and normal-force contexts
authenticate before their checkers are reused. The normal-force result is
related boundary evidence, not a prerequisite for the mass identity itself.
Ordinary checks authenticate stored bytes; explicit algebra replay is
separate. Evidence creation is bounded by 60 process-CPU seconds.

The next decision is whether evaluating this functional supplies the most
useful input to the remaining coupled junction equations. Any numerical
evaluation must independently control the new integral and its field/geometry
errors; the old 257 waveform samples alone do not certify it. A cap calculation
must retain its vacuum stress and radius-dependent background mass; the
Ricci-flat shortcut used here does not apply across that interface. Full
evolution, stability, finite-amplitude validity and the corrected observer
signal remain open. Other atlas routes remain available.
