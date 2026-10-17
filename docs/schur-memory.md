# Schur memory and the local effective Lambda target

For $D=\begin{pmatrix}D_p&\Phi\\\Phi^\dagger&D_c\end{pmatrix}$, eliminating
the exterior gives $\Sigma(z)=\Phi(z-D_c)^{-1}\Phi^\dagger$ and the interior
resolvent $(z-D_p-\Sigma)^{-1}$. The corresponding causal memory carries an
exterior initial-state drive and cross correlations. The intended physical
target is a local effective $\Lambda$ response from this unresolved exterior
structure. The finite identities, spectral controls and geometric candidates
below supply its retained dictionary; the physical identification of
$\Sigma$ with the actual exterior remains open.

## Step 7: averaging the explicit sourced initial slice

The next geometric calculation uses the existing
[cap--wall--bridge slice](nsc-atlas-connection-slice.md), from its regular
centre to the retained right-end sphere $R_{\mathrm{cut}}/\ell=4$. This is
the isolated sourced body, not an unevaluated common-parent connector.
The reference length remains the already supplied $\ell$.

Before evaluation, the finite question is whether proper-volume scalar
averaging of that slice supplies a vacuum-form correction beyond its fixed
cap volume term. The inputs remain $H\ell=1$, $k\ell=2$,
$M/\ell=0.4$, $R_w/\ell=1.5$ and the existing smooth vacuum slicing taper.
No observation selects them. The calculation imports the already checked
spatial geometry and wall jump; it derives its scalar averages directly.
[Buchert, gr-qc/0102049v2](https://arxiv.org/html/gr-qc/0102049v2), Sections
2.1--3.2, equations (5), (7), (10b), (12)--(16), supply a comparison for
the averaging identities. His fluid-comoving perfect-fluid assumptions
are not imposed on the moving material wall.

The budget is 30 CPU seconds and 20,000 evaluations of the existing spatial
fields, with no evolution or parameter scan. Fixed Gauss orders 32 and 64
are compared with adaptive quadrature at this one domain. The independent
controls are the integrated Hamiltonian including the shell delta,
direct wall-stress projection, the homogeneous cap limit and unit scaling.
A nonzero density is called a vacuum contribution only if its pressure
also has the required sign and magnitude. A non-vacuum result redirects
the next calculation toward a physical common-parent averaging domain
and its clocks rather than refining this isolated slice.

With the owned convention $K=-\mathcal L_T h/2$, set
$\theta=-\operatorname{tr}K$ and
$\sigma^2=\tfrac12(K^i{}_jK^j{}_i-(\operatorname{tr}K)^2/3)$.
Use $dV=4\pi R^2dl$, retaining both branches of the bridge. Then

$$
Q_D=\tfrac23(\langle\theta^2\rangle_D-\langle\theta\rangle_D^2)
      -2\langle\sigma^2\rangle_D,\qquad
3H_D=\langle\theta\rangle_D.
$$

The intrinsic scalar curvature includes
$-4[R']/R_w\,\delta(l-l_w)$. The wall's integrated slice energy is
$E_S=\sigma\gamma$ and its spatial stress trace is
$S_S=\sigma(\gamma^2-3)/\gamma$. Dropping that curvature delta would
remove a real term from the averaged constraint.

For a Gaussian-normal continuation at the initial instant, with the
domain transported along $T$, the scalar identities are

$$
3H_D^2=\kappa\langle E\rangle_D-
       \tfrac12(\langle{}^{(3)}R\rangle_D+Q_D),\qquad
3\ddot a_D/a_D=-\tfrac\kappa2\langle E+S\rangle_D+Q_D.
$$

This continuation specifies the clock for the instantaneous second
identity; it does not construct the later Einstein--shell solution.
Putting the curvature and variance into a spatially flat scalar fitting
model gives

$$
\kappa\rho_g=-\tfrac12(Q_D+\langle{}^{(3)}R\rangle_D),\qquad
\kappa p_g=-\tfrac12 Q_D+\tfrac16\langle{}^{(3)}R\rangle_D.
$$

Their trace projection is
$\Lambda_{g,\mathrm{tr}}=(Q_D-\langle{}^{(3)}R\rangle_D)/4$;
the vacuum condition is $3Q_D+\langle{}^{(3)}R\rangle_D=0$.
Scalar fitting is not a tensor averaging prescription or an observer
redshift--distance calculation. The supplied cap vacuum and the material
wall remain explicit in $E$ and $S$.

### Finite result and the next parent input

The [compact calculation](../results/summaries/nsc-atlas-geometric-average-v1.json)
gives $V_D/\ell^3=355.4931547895$ and

$$
H_D\ell=0.08352417231,\quad Q_D\ell^2=0.48101493580,\quad
\langle{}^{(3)}R\rangle_D\ell^2=0.35791137287.
$$

Consequently $\kappa\rho_g\ell^2=-0.41946315433$ and
$\kappa p_g\ell^2=-0.18085557242$. The trace projection is
$\Lambda_{g,\mathrm{tr}}\ell^2=0.030775890732$, but the necessary vacuum
condition has residual $(3Q_D+\langle{}^{(3)}R\rangle_D)\ell^2
=1.80095618027$. This correction is therefore not a vacuum source.
The largest 32/64/adaptive difference in the trace projection is
$1.0\times10^{-10}$; this is a quadrature comparison, not a certified or
physical-model uncertainty. The wall-inclusive Hamiltonian residual is
$1.2\times10^{-16}$ in these units. Seven focused controls and the
eighteen existing slice controls pass. No physical trajectory is run.

An outward extension can be assessed without scanning cutoffs. Define the
compact-body moments
$I_B=\int_B\theta\,dV$, $J_B=\int_B(2\theta^2/3-2\sigma^2)\,dV$,
and $\mathcal R_B=\int_B{}^{(3)}R\,dV$, including the wall.
On the owned time-symmetric right end, $K=0$ and ${}^{(3)}R=0$; enlarging
the isolated domain changes its volume but none of these moments. Hence

$$
\Lambda_{g,\mathrm{tr}}=
 \frac{J_B-\mathcal R_B}{4V_D}-\frac{I_B^2}{6V_D^2}.
$$

The isolated correction vanishes as the exact AE domain grows. This
analytic limit follows from the supplied vacuum end, rather than a
numerical extrapolation to cosmological scales.

For $N$ unchanged disjoint bodies and a hypothetical time-symmetric vacuum
connector, the vacuum condition would require volume per body
$v_* = 2I_B^2/(3J_B+\mathcal R_B)=23.1714109008\,\ell^3$.
The retained domain of just one body already occupies
$355.4931547895\,\ell^3$. Adding positive connector volume cannot meet
that condition. Joining copies through a $K=0$ vacuum connector is therefore
insufficient under this scalar averaging prescription. The existing
general gluing theorem does not assert a time-symmetric connector;
this is a conditional calculation, not a property assigned to its output.

For a vacuum connector with nonzero $K$, its constraint instead gives
$\mathcal R_C=-J_C$. The smallest additional scalar data are its proper
volume and the two moments $I_C,J_C$, evaluated with a physically selected
normal clock. For one retained body the vacuum test becomes

$$
V_D(3J_B+\mathcal R_B+2J_C)-2(I_B+I_C)^2=0.
$$

Choosing these moments to satisfy the equation would fit a desired form;
they must come from a compatible parent geometry and observer prescription.
The existing existence family does not supply them. This identifies the
necessary Step-5 input to Step 7 without requiring a new formation or
stability campaign. A common-parent development and tensor/observer
interpretation remain necessary before a cosmological Lambda is predicted.

For the local interior target there is a further, already established
[observer boundary](nsc-atlas-causal-classification.md#observer-access-preserve-direction-and-the-emission-event).
In the isolated reference an inward wall signal has finite central reception
only for $HR_w<1$; this averaging slice uses $HR_w=1.5$. Its parent volume
is not a causal coarse-graining domain for that central observer. Parent
moments would complete the parent scalar average, not identify a new central
vacuum term. That identification requires a physical reduction tied to the
retained observer and its compatible earlier data or quantum state. Global
causal exclusions from the isolated reference are not assigned to the
unevaluated many-body development.

## Step 7: the current candidate's first vacuum projection

The current conformal-scalar cap–wall candidate gives the exact smooth-bulk
result

$$
\boxed{\Delta\Lambda_{\rm classical}\,L_{\rm ref}^2=0.}
$$

This is the additional trace projection of its classical scalar stress. The
cap already contains the supplied term $\Lambda_{0,-}\ell^2=3(H\ell)^2=3$;
the exterior contains $\Lambda_{0,+}=0$. Those baseline values were inputs
to the geometry. The zero above holds for any positive reference length,
so it requires no choice of physical ruler or observed cosmological value.
For displaying the baseline, $L_{\rm ref}=\ell$ uses the existing reference
units; the P001 value $\ell=10000\,\mathrm m$ remains a fiducial scale.

Reuse the [owned action and trace identity](nsc-atlas-cap-wall-response.md),
with $\Lambda_0=\kappa_4 V_0$, rather than extending the cap computation:

$$
T^\phi_{\mu\nu}=\nabla_\mu\phi\nabla_\nu\phi
-\tfrac12g_{\mu\nu}(\nabla\phi)^2
+\tfrac16\left[G_{\mu\nu}\phi^2+
(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)\phi^2\right],
$$
$$
T^\phi{}^\mu{}_\mu=\phi(\Box-\mathcal R_J/6)\phi=0,
\qquad \mathcal R_J=4\Lambda_0,
\qquad
\Delta\Lambda_{\rm classical}L_{\rm ref}^2
=\tfrac14(\mathcal R_J-4\Lambda_0)L_{\rm ref}^2=0.
$$

This is exact on shell, including nonzero classical pulses. The pulse still
has nonzero stress, mass transport and wall response: a zero trace does not
make its stress tensor zero. In the decomposition
$T^\phi_{\mu\nu}=-(\Delta\Lambda/\kappa_4)g_{\mu\nu}
+T^{\rm TF}_{\mu\nu}$ with $g^{\mu\nu}T^{\rm TF}_{\mu\nu}=0$, the
vacuum projection is zero. This calculation concerns the unaveraged smooth
bulk. It does not assign a four-dimensional vacuum density to the material
wall, perform nonlinear metric averaging, or include quantum stress.

The [executable calculation](../src/recursive_horizons/nsc_atlas_effective_vacuum.py)
contracts the full improved stress and checks the result independently in
the Einstein frame. With $q=\Psi/\sqrt6$, $F=\operatorname{sech}^2q$ and
$\widehat V=\Lambda_0\cosh^4q$, use
$\mathcal R_E=6(\nabla_Eq)^2+4\widehat V$ and
$\Box_Eq=\widehat V_{,q}/6$. Then

$$
\mathcal R_J=F\left[\mathcal R_E+3\Box_E\log F
-\tfrac32(\nabla_E\log F)^2\right]=4\Lambda_0.
$$

The apparent field-dependent Einstein-frame potential therefore supplies
no extra Jordan-frame vacuum curvature. Seven
[checks](../tests/test_nsc_atlas_effective_vacuum.py) also reconstruct this
curvature from the actual nonlinear bulk operator's metric and field jets,
in both regions at zero and nonzero amplitude. The exact result has no
quadrature error; these checks do not estimate errors in an evolved solution.
Reproduce with
`python scripts/run_python.py -m recursive_horizons.nsc_atlas_effective_vacuum`
and `python scripts/run_python.py -m pytest tests/test_nsc_atlas_effective_vacuum.py -q`.
No endpoint replay, new spatial grid or observational fit enters this calculation.

For Step 7 this resolves the classical scalar route to an additional smooth
vacuum term in the current candidate. A nonzero contribution needs an
explicit additional mechanism, such as a specified quantum effective action
or a defined geometric averaging operation. The existing finite Schur
matrices have not been derived from this candidate's exterior. Its physical
response operator, state and metric variation must be supplied before those
matrices can yield an effective stress. In particular, the zero-frequency
coefficient of a scalar self-energy multiplies a field quadratic term;
it is not by itself a metric volume term or a cosmological constant.
The older two-sheet and Gauss–Bonnet results below retain their own actions.

## Step 7: conditional quantum curvature contribution

The next calculation imports one **quantized real massless conformal scalar**
on the quiet geometry, conditional on a Hadamard state with the usual
short-distance singularity and locally covariant renormalization. It uses the
four-dimensional field, including all angular sectors; the classical `l=0`
pulse alone cannot supply its vacuum fluctuations. The calculation applies
inside each smooth region, away from the material wall. It does not quantize
the metric or select a joint quantum state across the attachment.

For the scalar Laplace operator with $E=-\mathcal R/6$, the imported
[heat-kernel coefficient, equation (4.28)](https://arxiv.org/pdf/hep-th/0306138)
contains $(5\mathcal R^2+60\mathcal R E+180E^2)/360=0$ and
$(12\Box\mathcal R+60\Box E)/360=\Box\mathcal R/180$.
Including the $(4\pi)^{-2}$ factor gives the
[conformal trace coefficient, equations (15)–(16)](https://arxiv.org/pdf/hep-th/9308075):

$$
\langle T^\mu{}_{\mu}\rangle
=\frac{\hbar}{2880\pi^2}
  (\mathcal R_{\mu\nu\rho\sigma}^2-\mathcal R_{\mu\nu}^2)
+\hbar b\Box\mathcal R.
$$

The finite curvature-counterterm choice can change $b$. Both quiet metrics
have constant $\mathcal R$, so this ambiguity vanishes here. The renormalized
volume and Newton terms remain supplied action coefficients. The imported
heat-kernel and quantum-stress rules are assumptions of this extension;
the algebraic cancellation and geometric contractions are calculated locally.

Use the existing reference unit $L_{\rm ref}=\ell$, and set
$\zeta=\kappa_4\hbar/\ell^2$, $h=H\ell$, $m=M/\ell$, $r=R/\ell$.
Direct contraction for
$ds^2/\ell^2=-f\,dt^2+dr^2/f+r^2d\Omega^2$ gives

| Quiet region | $f$ | $\ell^4\mathcal R_{\mu\nu\rho\sigma}^2$ | $\ell^4\mathcal R_{\mu\nu}^2$ |
| --- | --- | --- | --- |
| Cap | $1-h^2r^2$ | $24h^4$ | $36h^4$ |
| Exterior | $1-2m/r$ | $48m^2/r^6$ | $0$ |

The local contractions also apply in the Schwarzschild interior; no global
Euclidean continuation or stationary state is inferred from this chart.
Consequently the trace projection $\Delta\Lambda_{\rm tr}=-\kappa_4\langle T^\mu{}_{\mu}\rangle/4$ is

$$
\frac{\Delta\Lambda_{{\rm tr},-}\ell^2}{\zeta}
=\frac{h^4}{960\pi^2},\qquad
\frac{\Delta\Lambda_{{\rm tr},+}\ell^2}{\zeta}
=-\frac{m^2}{240\pi^2r^6}.
$$

At the already supplied $h=1$, $m=2/5$ and quiet wall radius $r=11/20$,
these coefficients are **0.0001055428996274352** and
**−0.002440241780407079**, respectively. They multiply $\zeta$; neither
is an absolute cosmological Lambda. No cosmological data selected the
reference length or the coefficients.

The cap trace is constant, but that alone does not establish a stress
proportional to the metric: the wall can break the symmetries of the state.
The exterior trace projection varies as $r^{-6}$ and requires the rest of
the stress for conservation. A joint state, its traceless stress and boundary
exchange, the gravitational renormalization conditions, and a controlled
backreaction calculation are still needed for a full effective Lambda.
In particular, this calculation does not derive the supplied cap scale $H$.
It isolates a calculable quantum contribution while the classical result
above remains zero. No additional cap-grid evolution enters either result.

For the already declared P001 fiducial scale $\ell=10000\,\mathrm m$,
restore SI units with $\kappa_{\rm SI}=8\pi G/c^4$ and the quantum trace
normalization $\hbar c$. Thus

$$
\zeta=\frac{8\pi G\hbar}{c^3\ell^2}
=6.56537648967\times10^{-77}.
$$

Using [CODATA 2022 inputs](https://physics.nist.gov/cuu/pdf/wall_2022.pdf),
$G=6.67430(15)\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}$,
exact $h_{\rm Planck}=6.62607015\times10^{-34}\,\mathrm{J\,s}$ and
exact $c=299792458\,\mathrm{m\,s^{-1}}$, gives

| Quiet-region trace projection | $\Delta\Lambda_{\rm tr}\ell^2$ | Standard uncertainty from $G$ alone |
| --- | ---: | ---: |
| Cap at $h=1$ | $6.92928871865\times10^{-81}$ | $1.5573\times10^{-85}$ |
| Exterior at $m=2/5$, $r=11/20$ | $-1.60211060142\times10^{-79}$ | $3.6006\times10^{-84}$ |

These are scale illustrations conditional on the fixed geometry and the
fiducial length, not cosmological predictions or total error intervals.
Keeping $h,m,r$ fixed while doubling $\ell$ divides the dimensionless
contribution by four and the dimensional $\Delta\Lambda$ by sixteen.
The supplied cap baseline remains $\Lambda_0\ell^2=3$; this quantum
contribution does not explain that input.

The same [calculation module](../src/recursive_horizons/nsc_atlas_effective_vacuum.py)
derives the conformal cancellation and contracts Christoffel curvature for
both metrics. An independent sectional-curvature contraction checks the
result. The [focused tests](../tests/test_nsc_atlas_effective_vacuum.py)
cover these identities, normalization and flat limits as well as the
classical calculation and SI scaling. The commands above reproduce the exact expressions
and their displayed decimal coefficients; there is no quadrature error in
these symbolic results.

### Does this term generate the cap curvature?

The same quantum extension supplies a necessary test without a new evolution.
Let $X=(H\ell)^2$, $\lambda_0=\Lambda_0\ell^2$ and
$a=\zeta/(960\pi^2)$. The constant-curvature trace equation is

$$
3X=\lambda_0+aX^2.
$$

Set the supplied renormalized volume coefficient to $\lambda_0=0$ explicitly.
The roots are

$$
X=0\quad\hbox{or}\quad X=\frac{2880\pi^2}{\zeta}.
$$

Since $\zeta=8\pi\ell_P^2/\ell^2$, the nonzero root obeys
$H^2\ell_P^2=360\pi$. It lies outside the small-curvature regime
$H^2\ell_P^2\ll1$ in which neglecting quantum-gravitational corrections is
controlled. Thus this one-scalar anomaly does not generate a nonzero,
small-curvature de Sitter cap when the supplied volume term is zero.
This is a consequence of the displayed equation, not an assumption that a
preferred cosmological model must hold.

For a supplied $\lambda_0\geq0$ and $9-4a\lambda_0\geq0$, the lower branch
is $X=2\lambda_0/(3+\sqrt{9-4a\lambda_0})$ and tends to $\lambda_0/3$ as
$\zeta\to0$. It corrects a supplied curvature. Neither branch establishes
the full tensor equation, a joint state or a compatible quantum wall. The
calculation concerns a constant-curvature smooth cap with the selected free
scalar and renormalized gravitational coefficients; a nonconstant geometry
or another physical mechanism is a separate calculation.

The executable `constant_curvature_trace_condition()` checks the roots and
classical limit. A twelfth test derives the nonzero root directly from the
dimensional SI trace equation, independently of the dimensionless conversion.
No extra grid, fitted scale or observational input enters this test.

### Quantum interface: the actual transmission law

The material interface supplies a separate operator contribution. For the
quiet candidate, $[\phi]=0$ and $[n\phi]=k\phi$. These are the transmission
conditions of $D=-(\nabla^2+E)+U\delta_\Sigma$, with $E=-\mathcal R/6$ and
$U=k$. Use the [Gilkey–Kirsten–Vassilevich transmittal coefficients](https://arxiv.org/html/hep-th/0101105v1)
as imported local polynomials; do not square distributional curvature.
Their theorem concerns an elliptic operator on Riemannian manifolds. The
substitution below checks its invariant algebra on the Lorentzian quiet-wall
geometry. It does not construct a global Euclidean continuation, a Lorentzian
quantum state or a renormalized interface stress.

With the project's common normal from cap to exterior, the source's inward
forms are $L^+=-K_+$ and $L^-=K_-$. Thus
$S=L^++L^-=kI$, $D=L^+-L^-=-2\bar K$. The second heat coefficient's
surface density is proportional to $\operatorname{tr}S/3-U=0$.
For unit smearing, the fourth coefficient also contains curvature
contractions; keeping only its extrinsic-curvature terms would give an
incomplete answer.

In the source index convention $\rho_{ij}=R_{ikkj}$, put
$Q_{ac}=R_{abcb}$, with the tangential metric used for contractions. The
Gauss equation on the shared induced geometry gives

$$
Q_+-Q_-=\frac{k}{2}\big[(\operatorname{tr}D)I+D\big],\qquad
\operatorname{tr}(Q_++Q_-)
=-\mathcal R_+-\mathcal R_-+2\rho_+(n,n)+2\rho_-(n,n)
=-6H^2.
$$

The first identity cancels the difference-curvature and difference-shape
terms in Theorem 7.1's $\mathcal A_1$. The normal-curvature terms in
$\mathcal A_2$ cancel because $S=kI$. This uses no static Schwarzschild
observer and remains a contraction identity inside its horizon. The
remaining three groups are

$$
\mathcal A_1=0,\qquad
\mathcal A_2=56k^3-12kH^2,\qquad
\mathcal A_3=-60k^3.
$$

After the factor $1/[360(4\pi)^2]$, the local surface polynomial is

$$
a_4^\Sigma(1)=\int_\Sigma\sqrt{|\gamma|}\,q_\Sigma\,d^3y,
\qquad q_\Sigma\ell^3
=-\frac{b(b^2+3h^2)}{1440\pi^2}
=-\frac7{720\pi^2}\quad(h=1,\ b=2).
$$

The Schwarzschild mass and radius cancel in the full expression. This is
an interface coefficient, not a four-dimensional volume term. Its value
on this background does not establish that its off-shell metric variation
is merely a shift of the material tension. The full counterterm operators,
their renormalized coefficients and the joint-state stress remain needed
before obtaining a physical effective Lambda from the attachment.

There is an independent boundary-curvature check. For outward curvature
$K$, define
$J=\operatorname{tr}K\operatorname{tr}(K^2)
-\frac23\operatorname{tr}(K^3)-\frac13(\operatorname{tr}K)^3$.
The [Euler boundary polynomial](https://arxiv.org/html/hep-th/0208205v3),
equations (6)–(7) after conversion from inward normals, is
$4(J-2\widehat G_{ab}K^{ab})$.
Adding the cap and exterior sides gives
$B_\Sigma=4[J(K_-)-J(K_+)+k\widehat{\mathcal R}]
=4k(k^2+3H^2)$, using the intrinsic Gauss scalar.
The scalar's Euler coefficient $-1/[360(4\pi)^2]$ reproduces $q_\Sigma$.
This checks the local polynomial and its orientation; it introduces no
new gravitational coupling or global topological assertion.
The calculation is implemented by `conformal_interface_heat_coefficients()`
in the same module. Its checks and the independent Euler-boundary test bring
the focused module to fourteen passing tests.

### Connecting the actual exterior to the Schur operation

The operator to eliminate in this candidate is the exterior conformal wave
operator with the transmission domain above. Its boundary variable is
$\varphi=\phi|_\Sigma$. Define $\mathcal N_\pm$ to map that boundary value
to the corresponding **outward** normal derivative, with the initial or outer
data stated separately. The exterior outward normal is $-n$ and the cap's
is $n$. Therefore the homogeneous boundary equation is

$$
(\mathcal N_-+\mathcal N_++U)\varphi=0.
$$

It is exactly $[n\phi]=U\varphi$ in the common-normal convention. Supplied
incoming data add a driving term to this equation. In real time the maps
are causal history operators; the finite preparation does not justify a
stationary zero-frequency limit.

For an auxiliary elliptic problem with specified outer boundary data, write
the exterior field as its Dirichlet solution plus a fluctuation vanishing
on the interface. Integration by parts removes their bulk cross term.
Gaussian elimination then separates the boundary quadratic form
$\frac12\langle\varphi,\mathcal N_+\varphi\rangle_\Sigma$ from
$\frac\hbar2\operatorname{Tr}\log D_{+,D}$, where $D_{+,D}$ has homogeneous
Dirichlet interface data. The boundary measure is the induced geometric
measure, not a freely normalized matrix inner product.

The first term changes field transmission; the second and the metric
dependence of the first contribute to the effective stress. A numerical
coefficient of the first term alone is not a vacuum volume density. Real-time
expectation values additionally require the joint state's symmetric
two-point function, along with the retarded response. The existing prepared
classical pulse specifies a mean field, not that covariance. This identifies
the missing Step 3 input needed by Step 7: a physical quantum preparation
and its covariance, with compatible boundary exchange and renormalization.
The old finite Schur matrices remain useful identities, but are not yet
evaluations of these physical operators.

### Exterior response to a prescribed wall field

The [prepared pulse](nsc-atlas-prepared-response.md) used the particular
condition $\chi_U|_W=0$ selected by its quiet cap and no-return interval.
The exterior operator needed above instead prescribes the wall value
$y(V)=\chi(U_W(V),V)$. Keep the same quiet geometry and finite wedge,
with source value $g(V)$ at $U=U_0$ and bottom value $b(U)$ at $V=V_a$.
Here $\chi=r\phi$ and all coordinates use the existing reference unit.

For a spherical harmonic of degree $l$, direct reduction of the wave
operator gives

$$
\chi_{UV}+K_l\chi=0,\quad
K_l=-\frac{r_{UV}}r+\frac{C\,l(l+1)}{4r^2},\quad
C=\frac{32m^3}r e^{-r/(2m)},
$$
$$
K_l=\frac{16m^4}{r^4}e^{-r/(2m)}
 +\frac{8m^3l(l+1)}{r^3}e^{-r/(2m)}.
$$

The exterior is Ricci flat, so its conformal-curvature coupling contributes
no additional bulk potential. Define the entry time
$v_0(U)=V_a$ for $U\le U_a$ and $v_0(U)=U_W^{-1}(U)$ otherwise, and
$a(U)=b(U)$ or $y(v_0(U))$ on those respective intervals. The Dirichlet
integral equation is

$$
\chi(U,V)=a(U)+g(V)-g(v_0(U))
-\int_{v_0(U)}^Vdq\int_{U_0}^Udp\,K_l(p,q)\chi(p,q).
$$

Every integration point is retained because $U_W(q)\ge U$ when
$q\ge v_0(U)$. Differentiating in $V$ and then $U$ recovers the bulk
equation. The empty integrals recover all three prescribed faces when
$b(U_0)=g(V_a)$ and $y(V_a)=b(U_a)$.

Writing $I_l(V)=\int_{U_0}^{U_W(V)}K_l\chi\,dU$, the first wall jets are

$$
\chi_V|_W=g'-I_l,\qquad
\chi_U|_W=\frac{y'-g'+I_l}{U_W'}.
$$

At the wall-bottom corner a common first jet additionally requires
$y'(V_a)=U_W'(V_a)b'(U_a)+g'(V_a)-I_l(V_a)$.
Further smoothness requires the corresponding higher-jet transport
conditions. Boundary derivatives are supplied or evaluated through these
equations, rather than estimated from stored waveform samples.

For the proper-time tangent $u=(\dot U,\dot V)$, with
$C\dot U\dot V=1$, the exterior outward normal is
$\nu_+=(\dot U,-\dot V)$. Consequently

$$
\mathcal N_+^R\varphi+\text{initial-data drive}
=\dot U\left(\frac{\chi_U}{r}-\frac{yr_U}{r^2}\right)
-\dot V\left(\frac{\chi_V}{r}-\frac{yr_V}{r^2}\right),
\qquad \varphi=y/r.
$$

The radius terms follow from differentiating the physical scalar and are
essential. This determines the actual causal field-response map associated
with the exterior and its induced boundary measure. Its quantum covariance
and metric variation are additional ingredients for effective stress.
For fixed geometry and identical supplied traces, a uniform field error
$e_\chi$ gives a conditional wall-normal error bound
$|\delta(\nu_+\phi)|\le 2\dot V K_{l,\max}(U_W-U_0)e_\chi/r$.
This follows from $|\delta I_l|\le K_{l,\max}(U_W-U_0)e_\chi$;
it supplies no bound on the numerical field error itself.

Both terms in $K_l$ decrease over $0.4\le r\le0.55$. The inherited analytic
area bound therefore gives
$q_l\le q_0[1+l(l+1)/2]$, with $q_0=0.13882099234357653$.
The bounds for $l=0,1,2,3$ are respectively
$0.138821,0.277642,0.555284,0.971747$, all below one. This proves uniqueness
for each of these mode problems in the finite continuous wedge and bounds
the inverse by $1/(1-q_l)$. It does not supply an angular-tail estimate or
a four-dimensional quantum determinant.

The admitted implementation checks this operator with at most 120 CPU
seconds and 30 million kernel values, using existing geometry and product
integration. It checks the zero-kernel limit, the normal and radius factors,
corner compatibility and recovery of the already solved prepared pulse as
Dirichlet input. Numerical product-integration and geometry errors remain
distinct from the analytic contraction. No source selection, new physical
trajectory, parameter scan or state calibration enters those checks.
The [implementation](../src/recursive_horizons/nsc_atlas_exterior_dtn.py)
uses a fixed-U grid with analytic entry values at the moving wall. Each
time segment owns its active nodes, including the endpoint Runge-Kutta
stages, so a new node retains its prescribed value until entry.
Six [controls](../tests/test_nsc_atlas_exterior_dtn.py) pass, including a
constant physical scalar on the curved exterior and preservation of entry
data. That one 129-node control gives maximum outward derivative
$5.993\times10^{-6}$ instead of its exact zero; it is a numerical diagnostic,
not an error enclosure. The existing 129-node prepared pulse is recovered
within the stated numerical tolerances, with no new physical input.

## The transition and vacuum-form coefficient

The frozen positive-$F$ DEV3 transition uses radius $32/25$, expansion
magnitude $225/64$, thickness $1/50$ and endpoint $F=4$. Its coefficients obey

$$
c_{\rm transition}=-\frac{18}{1015},\qquad
v_{\rm room}=-6c_{\rm transition}=\frac{108}{1015},\qquad
\alpha_{\rm GB}/L_*^2=-\frac1{8c_{\rm transition}}=\frac{1015}{144},\qquad
\alpha_{\rm GB}v_{\rm room,physical}=\frac34.
$$

The explicit two-sheet carrier is

$$
ds_5^2=e^{-36y^2/1015}
\left[A(\rho)dt^2-\frac{d\rho^2}{A(\rho)}-(\rho^2+1)d\Omega_2^2-dy^2\right].
$$

Both areal branches are needed: one radius can have distinct parent and
child metric values, joining at the throat. Pure Einstein gravity requires
a negative null source for this carrier; a canonical positive-kinetic scalar
alone cannot provide it. The transition-fixed Gauss–Bonnet coefficient
cancels that local null projection geometrically. On the room, its full
component identity is

$$
[G^A{}_B+\alpha_{\rm GB}H^A{}_B]_{y=0}
=\operatorname{diag}\left(\frac{108}{1015},\frac{108}{1015},
\frac{108}{1015},\frac{108}{1015},
-\frac{R_4+(1015/144)\,GB_4}{2}\right).
$$

Thus the four tangential components have vacuum form and the normal component
contains the room curvature action. The same $\alpha_{\rm GB}=1015/144$ lies
at the five-dimensional EGB critical point: with
$\Lambda_0=-108/1015$, $\alpha_{\rm GB}\Lambda_0=-3/4$ and the vacuum
discriminant vanishes. Kinetic and characteristic conclusions must therefore
use the actual two-sheet background. Bulk-source microphysics and full
stability remain open.

The expanding-child Ricci sign selects positive $\Lambda_{\rm effective}$
in the retained doubled-FLRW convention $R=-2\Lambda_{\rm effective}$.
It selects $\Phi^2/\Lambda^2=1006/1015$. The retarded boundary pole then has
positive weights and subluminal group velocity, and its propagating band
ends at $p^2+\Phi^2=\Lambda^2$. The squared bandwidth is $9/1015$, half the
transition-curvature magnitude. This boundary excitation is not an electron
identification or the full graviton propagator.

The imported RS2 black-universe equation uses
$G^\nu{}_\mu=-\Lambda_4\delta^\nu{}_\mu-\kappa_4^2T^\nu{}_\mu
-\kappa_5^4\Pi^\nu{}_\mu-E^\nu{}_\mu$ and
$\Lambda_4=(\Lambda_5+\kappa_5^4\lambda^2/6)/2$.
The outside-kernel variation projects to $-E^\nu{}_\mu$ in that convention.
It matches the black-universe exterior response dictionary; the bulk metric
and its healthy response spectrum still have to determine that projection.

## Spectral and boundary controls

Room geometry forces originate in the regulated fermionic spectrum and its
scale anomaly. The compact allocation separates one massless Dirac field,
two fields at each $m_n=n\pi/\ell$ for $n\geq1$, and the warped complement.
It supplies a Euclidean subtraction, not a real-time state. With
$J=\ell_{\rm proper}/\ell$, the compact light limit is

$$
\lim_{u\downarrow0}[h_\Lambda(u)-E_1(u/\nu^2)]
=2\log\frac{\Lambda J}{\nu}
+2\sum_{n\geq1}E_1\!\left[\left(\frac{n\pi}{\ell J\Lambda}\right)^2\right].
$$

The light singular slope is $1/J^2$ and its canonical mass stays zero.
The imported spectral bulk/boundary constraint is $c_K=2c_R$, leaving one
coefficient for $c_R\int R+c_K\int K$. This universal ratio differs from
the particular transition coefficient $-18/1015$.

For the exponential heat profile,
$h(z)=\int_0^1e^{-\alpha(1-\alpha)z}d\alpha$ and
$F_2(z)=(1+z/4)h(z)-2$. Its nearest complex zeros lie beyond the unit cutoff
disk. Canceling the stationary vacuum background gives
$K_2(z)=F_2(z)+1$, $K_2(0)=0$, $K'_2(0)=1/12$ and a positive massless
propagator residue $12$, with no additional pole in that disk. These flat
and Euclidean controls leave the curved retarded positivity test open.
Reciprocal room symmetry requires $B(\chi)=B(-\chi)$: for
$B=b_0+b_2\chi^2+b_4\chi^4$, the stationary linear scalar mixing is zero,
while $B''(0)=2b_2$ retains the nonlinear transition channel.

The first-order curvature EFT change of inverse metric is

$$
\delta g^{\mu\nu}=-\frac{2c_W}{A_4}
\left[R^{\mu\nu}-\frac{R}{6}g^{\mu\nu}
+\frac{T^{\mu\nu}-Tg^{\mu\nu}/3}{2A_4}\right]
+\frac{c_R}{A_4}\left(R-\frac{T}{2A_4}\right)g^{\mu\nu}.
$$

Here old inverse metric equals new inverse metric plus this shift. Boundary
data, observables and the quantum Jacobian must follow the same perturbative
map; it is not a fundamental fourth-order theory.

## Finite state and energy controls

The coordinate generator is gravity plus $M\operatorname{Tr}(CH)$ with
$M=4\kappa$ counted once. Nodal force partials must be divided by node spacing
before conversion to stress densities. With $N=rL$, $q=rQ$ and
$V=4\pi qr^2$, the unit normal is $N^{-1}(\partial_t-\beta\partial_x)$,
whose dual is $Ndt$. Coordinate energy
$e=L F_L+\beta F_\beta=NV\rho-\beta qVj$ differs from normal-observer energy.
The latter obeys

$$
\partial_t(V\rho)+\partial_x[V(Nj/q-\beta\rho)]
=-NV(p_rK_r+2p_\perp K_\perp)-V(j/q)N_x.
$$

Independent metric variations control this identity and the chart-work
ledger. Mode-projection leakage is extra exchange; the controls establish
neither a continuum limit nor renewal.

The frozen six-mode incoherent imbalance witness uses the Cesàro projection
onto the commutant, with finite-time remainder
$\|C_T-\bar C\|_F\leq2\|Y\|_F/(\gamma T)$. Its retained mean populations are
$(1.3837895948404104,1.0058995585808204,0.6103108465787691)$;
the channel current is $0.0012511140479781876$. Reprojection moves the
response by $0.004920298239642384$, whereas the frozen shortcut gives
$0.00030421076749148856$. These finite controls do not close the physical
local gate. The smooth-motif Weyl-disk tail also controls a projected
resolvent only; it is not identified with a dark-energy density or a cosmic
fraction. Original calculations are at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).

The retained coordinate Hamiltonian and normal-observer stress API lives in
[nested.energy_exchange](../src/recursive_horizons/nested/energy_exchange.py).
It keeps the generator gravity energy plus `M Tr(C H)` with `M=4 kappa`
counted once, the finite regional controls, and the separate normal-observer
shell energy. Renewal and a continuum limit remain unclaimed.
