# Charged light restoration in the common KS action

For the retained magnetic light field, the exact angular continuation gives

$$
Z_q(-1)=\frac1{30}-\frac{q^2}{6},\qquad
Z_{q,\ne0}(0)=-|q|-\frac13,\qquad
\operatorname{FP}_{s=1} Z_q(s)=4\gamma_E-2H_{|q|}.
$$

These quantities fix the charged cylinder normalization and the finite
fourth-order restoration. Their gauge contribution agrees with the existing
Dirac heat coefficient $(2/3)F^2$. The new action reproduces the charged
homogeneous restoration plus the physical-metric LLL geometry under all four
compact metric variations.

The [record](../results/development/nsc-magnetic-light-restoration.json)
closes this light-allocation/action connection. It supplies neither the full
state integral nor a complete incoming stress or constraint solution. The
previous charged-source and neutral-reference records remain unchanged.

## The same spectrum and normalization

For $b=|q|/2$, reuse the already declared nonzero magnetic angular spectrum:

$$
Z_q(s)=4\sum_{\ell=b+1}^{\infty}
\ell(\ell^2-b^2)^{-s}.
$$

Only $|q|=4$ is evaluated as a physical sector. The $q=0$ case checks the
existing [neutral cylinder normalization](nsc-angular-stress.md). It is not
a flux-sector scan. The subtraction scale stays at $\mu=1$.

The differentiated Hurwitz continuation is

$$
\begin{aligned}
Z_q'(-1)=4\bigg[&2\zeta'(-3,b+1)+b^2\zeta(-1,b+1)
-2b^2\zeta'(-1,b+1)\\
&+b^4\left(\frac14+\frac{\psi(b+1)}2\right)
-\sum_{j=3}^{\infty}
\frac{b^{2j}\zeta(2j-3,b+1)}{j(j-1)}\bigg].
\end{aligned}
$$

The $j=2$ pole multiplied by a zero must be retained. Setting $s=-1$
termwise before continuing would lose its finite contribution. The verifier
evaluates the continuation away from the poles and independently takes the
special-point limits and derivative.

The finite cylinder density and fourth-order harmonic are

$$
C_q(R;\mu)=
-\frac{Z_q'(-1)+[\log(\mu^2R^2)+1-\gamma_E]Z_q(-1)}
{16\pi^2R^4},\qquad
h_q=\log\mu+\frac18\operatorname{FP}Z_q(1).
$$

The barred sphere below has $R=1$. At $q=0$ these expressions reproduce
`static_cylinder_density()` and $\gamma_E/2$. At $q=4$, $H_4=25/12$.
There is no new normalization choice or adjustable coefficient.

For $a=b+1$, $\varrho=(b/a)^2$ and first omitted $J\ge3$, the exact rational
tail bound used for the differentiated series is

$$
|\Delta Z_q'(-1)|\le
\frac{4a^3\varrho^J}{J(J-1)(1-\varrho)}
\left(1+\frac{a}{2J-4}\right).
$$

This is a spectral truncation bound. It is kept separate from floating-point
evaluation of the Hurwitz constants. The record retains the term count,
precision, exact rational bound and a second series tolerance comparison.

## Charged trace and LLL allocation

The existing [compact matching](nsc-compact-matching.md) owns the gauge heat
coefficient. The charged light trace is

$$
\mathcal A_{4,q}=\mathcal A_{4,0}
+\frac{1}{16\pi^2}\frac23F^2
=\mathcal A_{4,0}+\frac{q^2}{48\pi^2r^4}.
$$

The identity
$-[Z_q(-1)-Z_0(-1)]/8=q^2/48$ binds its coefficient without another gauge
calculation. Also $Z_{q,\ne0}(0)+|q|=-1/3$: the LLL zero-mode count is
included once in the full angular heat coefficient.

The [homogeneous tensor owner](../src/recursive_horizons/nsc_magnetic_light_reference.py)
returns the **nonzero-angular restoration only**. It adds the charged gauge
Weyl cocycle to the existing neutral four-dimensional conformal map and
subtracts the LLL two-dimensional Weyl cocycle. Consequently its trace is

$$
\mathcal A_{4,q}-\frac{|q|R_{2,\rm physical}}{96\pi^2r^2}.
$$

Adding the existing physical LLL geometric tensor restores exactly
$\mathcal A_{4,q}$. Both the nonzero-angular restoration and that sum satisfy
the existing homogeneous source-conservation equation. The LLL thermal/state
piece remains separate.

## General action in the canonical KS foliation

Set $\bar g_2=g_2/r^2$, $\bar r=1$ and $\sigma=\log r$. The same complete
light restoration is

$$
S_{\rm light,restore}=S_{\rm cyl}+S_{\bar R_2^2}
+S_{\rm LLL,geo}[\bar g_2]+W_{4,q}[\bar g_4,\sigma],
$$

$$
S_{\rm cyl}=-4\pi C_q\int\sqrt{|\bar g_2|},\qquad
S_{\bar R_2^2}=\frac{h_q}{240\pi}
\int\sqrt{|\bar g_2|}\,\bar R_2^2.
$$

The LLL term is the [existing canonical ADM allocation](nsc-lll-geometric-history.md),
applied to $\bar N=N/r$ and $\bar a=a/r$. Its spatial foliation must be the
same as the canonical source. The API accepts explicit KS coordinates or
a reparameterization of KS time; it rejects a silent substitution of PG
coordinates.

With $g_{s,2}=r^{2s}\bar g_2$, $r_s=r^s$, and
$(\alpha,\beta,\gamma)=(-1/20,11/360,-1/30)$, the remaining action is

$$
\begin{aligned}
W_{4,q}=-\frac{4\pi}{16\pi^2}\int\sqrt{|\bar g_2|}\bigg[
&\alpha\sigma\bar C^2
+\beta\sigma\int_0^1r^{4s}E_4[g_s]\,ds\\
&+\frac\gamma{12}(\bar R_4^2-r^4R_4^2)
+\frac23\sigma\bar F^2\bigg],\qquad \bar F^2=q^2/2.
\end{aligned}
$$

The three-point Gauss rule integrates the formal Weyl parameter. That
parameter is not physical time, a matter field or a new history. Analytic
product/power rules supply the transformed metric/radius jets; all curvature
contractions reuse `spherical_invariants`.

[`LightRestorationAction.actions`](../src/recursive_horizons/nsc_light_restoration_action.py)
returns separate action channels for arbitrary supplied spherical metrics
and analytic first/second jets. The action **already includes LLL geometry
once**. The final source must add its LLL state contribution, the nonzero
angular state-minus-reference integral, and the locked compact action in
their existing allocations. Another physical LLL geometric term would count
it twice. This owner introduces no independently weighted action term.

## Focused evidence and replay

The record keeps the following checks separate:

- exact charged special values, zero-mode counting and gauge heat matching;
- neutral cylinder/harmonic normalization and differentiated-series convergence;
- charged trace, LLL-once counting and homogeneous source conservation;
- all four compact action derivatives against the independently evaluated
  homogeneous tensor, including neutral control;
- resolved coordinate quadrature, exact conformal path endpoints and
  three-point versus five-point formal Weyl quadrature.

The action variation uses the old geometry in $t=-\theta$, with the same KS
spatial slices. The integration interval and boundary-flat variations are
numerical controls, not selected physical Cauchy data or a history duration.
General boundary variations still require their own retained endpoint terms.

```sh
python3 scripts/check_nsc_magnetic_light_restoration.py --check
```

This is an all-field replay with exact structure, fractions and hashes, plus
declared float tolerances. It runs the new allocation controls and reads the
frozen source records as provenance pointers. It launches no historical
generator, state/scattering solve, metric evolution or publication step.
