# Field, mass and clock on the moving wall

The [rectangle calculation](nsc-atlas-bulk-response.md) supplies the first
smooth exterior metric. The same leading scalar also determines radial
mass and lapse integrals across the full retained wedge. Their wall limits
connect to the [two-sided mass and displacement laws](nsc-atlas-cap-wall-response.md)
without another free boundary source. The
[calculation admission](prior-art-and-open-claim.md#exterior-wedge-and-moving-wall-calculation-admission-15-october-2026)
keeps the original action, preparation, normalization and finite domain.

All quantities below are dimensionless formal coefficients at order
$\alpha=\kappa_4\varepsilon^2/\ell^2$. No finite amplitude is selected.
Use $a=2m=0.8$, $s=\tau_W/\ell$, $x_s=v>0$, and the original signed
cap-to-exterior normals $\beta_\pm$. The source and wall functions are
those of the existing leading solution.

## The field on the retained wedge

Write $\chi/\varepsilon$ as $\chi$ in this section and set $p=\chi/r_0$.
The exterior equation is $\chi_{UV}=-K\chi$. For a fixed $U$, its
initial point is

$$
L(U)=\begin{cases}V_a,&U\le U_a,\\ U_W^{-1}(U),&U>U_a.\end{cases}
$$

The quiet bottom and the transmitted wall law both give $\chi_U=0$ at
that point. Consequently

$$
\chi_U(U,V)=-\int_{L(U)}^V K(U,q)\chi(U,q)\,dq,\qquad
\chi_V(U,V)=g'(V)-\int_{U_0}^{U}K(q,V)\chi(q,V)\,dq.
$$

The existing dense Volterra solution supplies the integrands. In particular,
$\chi_U=0$ at the wall does not set the physical derivative to zero:
$p_U=-p r_{0,U}/x$ there. Differentiating waveform samples would lose this
distinction. Numerical queries must respect $U\le U_W(V)$ before the
exterior kernel is evaluated. An inverse-wall bracket and the omitted
interval, if any, belong in the derivative error account.

## Independent paths to the same wall mass

In the Einstein frame, constant-$V$ transport gives

$$
\mu_E(U,V)=\mu_{E,0}(V)
-\int_{U_0}^{U}\frac{r_0^2r_{0,V}}{C_0}p_U^2\,dq,
$$
$$
\zeta_E(U,V)=\zeta_{E,0}(V)
+\int_{U_0}^{U}\frac{r_0}{2r_{0,U}}p_U^2\,dq.
$$

Here $Z=C/(2r_U)$, $Z_0=-2a/V$, and $\zeta_E$ is the coefficient
of $\log(Z_E/Z_0)$. Source values are the already fixed characteristic
data. The first integrand is nonnegative and the second is nonnegative in
the retained sector. This gives useful orientation checks.

Transforming back at the wall gives

$$
\mu_+^{\mathrm{radial}}=\mu_{E,W}+Q_+,\qquad
\zeta_{J,W}=\zeta_{E,W}+\frac{p^2}{4}.
$$

The second identity includes the transformed radial derivative. Adding only
the conformal correction to $C$ would give the wrong lapse. Independently,
the timelike mass laws give

$$
\mu_+=\int_{s_a}^{s}L_+\,ds'+Q_+,\qquad
\nu_-=\int_{s_a}^{s}L_-\,ds'+Q_-.
$$

The definitions of $L_\pm,Q_\pm$ and their independent improved-stress
derivation remain in the cap/wall owner. The cap term $\nu_-$ subtracts
the vacuum mass at the actual radius. Comparing $\mu_+^{\mathrm{radial}}$
with $\mu_+$ checks different gravitational integration paths with the same
quiet anchor and leading scalar; the equality is not imposed numerically.

## Displacement and the exterior clock

Set $D=\mu_+-\nu_-$ and
$F=\mu_+/x+\beta_+D/(b x^2)$. The angular junction yields

$$
\rho_s=v_x\rho+F/v,\qquad
\frac{\rho(s)}{v(s)}=\int_{s_a}^{s}\frac{F}{v^2}\,ds',\qquad
\delta\beta_+=\beta_{+,x}\rho+\frac{D}{b x^2}.
$$

The cap form of the displacement equation is an independent algebraic
check. All denominators remain regular at the zero of $\beta_+$.

For $V_W=V_0+\alpha V_1$, proper normalization of the wall tangent and
the signed normal give exactly

$$
1=C_J U_{W,s}V_{W,s},\qquad
r_{W,s}-\beta_+=2r_{J,U}U_{W,s},\qquad
V_{W,s}=\frac{1}{Z_J(r_{W,s}-\beta_+)}.
$$

Let $d=v-\beta_+>0$ and $\eta=V_1/V_0$. Since
$V_{0,s}=-V_0/(2ad)$, its linearization is

$$
E=\zeta_{J,W}+\frac{\rho_s-\delta\beta_+}{d},\qquad
\boxed{\eta_s=\frac{E}{2ad}},\qquad \eta(s_a)=0.
$$

This includes the background evaluation term from the changed $V_W$.
Under the residual constant relabeling $U\mapsto cU$, $V\mapsto V/c$,
$\mu_+$, $\nu_-$, $\rho$, $\zeta$ and $\eta$ are unchanged. $V_1$ itself
rescales as $V_0$ does. The coefficient $\eta$ describes the wall's
exterior null embedding; the central observer has a separate clock.

## Bounded evaluation and comparison

The [wedge evaluator](../src/recursive_horizons/nsc_atlas_wedge_field.py)
queries the original dense field solution with paired space/time arguments.
Its derivative quadrature splits at the two source support edges and the
next spatial-node activation. Inverse wall labels use bounded bisection on
the monotone reference spline, choosing the retained side. The reported
omitted-interval estimate uses the analytic kernel and field bounds. It is
conditional on that numerical spline bracket and excludes spline and
floating evaluation error; it is not a total derivative enclosure.

The [wall integrator](../src/recursive_horizons/nsc_atlas_wall_evolution.py)
uses the same field's first derivatives. Radial Gauss rules split at each
active field cell, including its partial last cell. Proper-time mass and
displacement integrals use dense DOP853 evolution split at the field's
activation events and source edges. A separate four-point Gauss rule per
time panel integrates the exterior clock; partial-panel values use the
integrated cubic through those rates. This clock quadrature has its own
time-panel refinement.

The [comparison producer](../scripts/record_nsc_atlas_wall_evolution.py)
varies one setting at a time: field grid, radial quadrature, field derivative
quadrature, clock time panels, geometry interpolation, source quadrature or
wall ODE tolerance. A complete comparison shares a 120 process-CPU-second
and 30-million-scalar-kernel budget across its two wall solutions and field
queries. Compact records retain nine common time cuts, maximum defects also
checked at clock quadrature nodes, and the execution counts. Refinement
differences measure numerical sensitivity rather than a rigorous enclosure.
The [current result](current-result.md) owns the numerical outcomes.

Existing rectangle evidence remains bound to commit `3470152e1766` through
the explicit `check-bulk-v1` and `bulk-v1` historical commands. The new dense
field interface does not replace those records or claim that the historical
record was produced by changed source.

## Next connection

These formulas determine formal mass, displacement and exterior-clock
coefficients once evaluated with controlled field and quadrature resolution.
The smooth cap reconstruction must retain its scalar-induced potential,
regular centre and fixed central clock. It supplies the next observer map.
The higher-order scalar also needs a stated continuation of the physical
source and the moving interface. Finite-amplitude control, formation and
stability use that coupled solution rather than the leading coefficients
alone.
