# Noncompact limit of the finite stationary query

This independent assessment explains the saved
[critical-v3 result](../results/development/nsc-discovery-stationary-critical-v3.json)
without another optimization, source eigensolve or trajectory. Its
[consumer](../scripts/assess_nsc_discovery_stationary_limit.py) and
[tests](../tests/test_nsc_discovery_stationary_limit.py) authenticate the old
producer closure against `697efe0155cdbf796c0174468c90a4c368bb7dbe`, its locked
inputs, and its unchanged JSON/NPZ. The
[stationary owner](nsc-discovery-stationary.md) retains the construction.
This note does not establish nonexistence of finite branches elsewhere.

## Resolved small-amplitude lemma

For $S=1+a\cos(kx)$, $Q=sS$, $k=2\pi h/L$, the actual finite radial operator is
$L_Q=-3D^2+Q^2-D(DQ/Q)$. Expanding its lowest critical eigencondition gives

$$
s=|a|k/\sqrt6+o(|a|),\qquad R=1-\tfrac a3\cos(kx)+O(a^2).
$$

The sign in $R$ is negative: the first perturbation of $-D(DQ/Q)$ is
$+ak^2\cos(kx)$. Its second-order ground-energy contribution is
$-a^2k^2/6$, cancelled by $s^2$. The projected auxiliary SPD equation gives

$$
\chi=-12\cos(kx)/a+O(1),\qquad
\mathcal W=\alpha\langle Q^2\chi^2\rangle\longrightarrow12\alpha Lk^2.
$$

The modes needed through this order are $0,h,2h,3h$. For the saved NF128,
$h=8$ and these fit the retained geometry band $|m|\le63$ and fine derivative
band $|m|\le255$. The leading formulas therefore do not require replacing the
finite auxiliary equation by a continuum product rule. The saved best also has
higher coefficients of order $a^2$ or smaller; they leave these leading limits
unchanged. They can still affect local force components.

As $Q\to0$, the six occupied positive AP eigenvalues approach
$(2\pi/L)(\tfrac12,\tfrac12,\tfrac32,\tfrac32,\tfrac52,\tfrac52)$.
With base weights $(.75,.75,.5,.5,.25,.25)$ and $M=4\kappa$ once,

$$
E_0=M\sum n_a\epsilon_a=7\pi M/L=28\kappa\pi/L,
\qquad
\eta_0=\frac{12\alpha Lk^2}{E_0}
=\frac{12\pi\alpha h^2}{7\kappa}.
$$

Thus $E_0=10.99557428756$ and $\eta_0=1.21871531420$ at the saved
$L=8,h=8,\kappa=1$. Both source energy and auxiliary limit scale as $L^{-1}$
at fixed integer $h$, so this ratio is independent of $L$. The AP momentum
pairs do not add another angular copy factor.

The exact normal denominator is
$16\pi A\langle Q^2R^2\rangle\sim16\pi ALs^2$.
At the auxiliary critical envelope, the leading normal remainder is
$\langle Q,\mathrm{rem}\rangle\to-2\mathcal W$;
magnetic and fixed-source mass derivatives vanish at this order. Consequently

$$
r_{\mathrm{amp}}^2\sim\frac{9\alpha}{\pi Aa^2},\qquad
r|a|\to3\sqrt{\alpha/(\pi A)},\qquad
rQ\to k\sqrt{3\alpha/(2\pi A)}.
$$

The saved radius diverges along this limit while $rQ$ and radial proper period
length have finite limits. These are parameter-space asymptotics, not time
evolution or a physical singularity claim.

## Coordinate conditioning and unsatisfied force

Let $F=\dot p_Q=r_{\mathrm{amp}}^2P_R+\mathrm{rem}$ and
$\delta Q_j=s\cos(jkx)-\mu_jQ$. In exact arithmetic the chosen normal radius
has $\langle Q,F\rangle=0$, so

$$
\partial_{a_j}\eta=\langle\delta Q_j,F\rangle/E_{\mathrm{base}}
=s\langle\cos(jkx),F\rangle/E_{\mathrm{base}}.
$$

The conversion back to cosine force coefficients has gain
$2E_{\mathrm{base}}/(sL)$, which diverges as $s\to0$.
A shrinking coordinate gradient is therefore not uniformly equivalent to a
shrinking local force. The actual physical variations hold the selected $\eta$
fixed. Differentiating the source-preparation function $\eta(a)$ inside the
physical force would add an unowned term.

For the saved best, $a_1\simeq.0003685593$, $s\simeq.0009453913$,
$r\simeq1287.02$, $\eta\simeq1.21871369$, and the local retained
$\dot p_Q$ remains about $32.96$. Its normal denominator is about
$1.618\times10^{-5}$, versus $3.2103$ initially. The gain from gradient to
force coefficients is about $2907.68$. This supports a noncompact escape
explanation for the saved search, rather than a finite stationary solution.

## Independent arithmetic assessment

The old strict `--check` fails its absolute $10^{-11}$ gradient identity
comparison. The producer evaluates $\Delta x_g(T^TB)/E$; the checker evaluates
$(\Delta x_gT^T)B/E$. In this saved case the arithmetic reordering gap is
about $2.05\times10^{-10}$. The new assessment records that failure as
`numeric_identity_unresolved`; it does not modify the old checker or call the
old numeric replay successful.

The radius coefficient is independently extracted using the original
geometric $\dot p_Q$ at unit $R$, zero $\chi$, zero field columns and momenta,
then subtracting the known magnetic term $-4\pi C_F\mathrm{flux}^2Q$.
This isolates the Einstein affine coefficient without subtracting two large
auxiliary forces. The zero fields occur only in this algebraic coefficient
call; the actual saved candidate remains unchanged. The consumer verifies this
coefficient against ordinary unit/2R extraction on the well-conditioned saved
initial state, then reports the selected state's ordinary/stable coefficient,
normal-pairing and gradient differences.

At the saved large radius, ordinary affine terms each reach about
$1.92\times10^7$ and their extrapolation differs from the directly evaluated
selected force by about $.255$. Small mean-Q pairings and large $\mu$ amplify
arithmetic defects. The assessment reports these observed defects alongside
$\epsilon\sum|x_iy_i|$ weighted-dot indicators. These indicators are neither
certified bounds nor a replacement acceptance test. Stable extraction is a
diagnostic; no alternate candidate, source, or old outcome is substituted.

The new assessment can authenticate its own record while explicitly retaining
both the old strict replay failure and the unsatisfied full force. Long-time
holding, stability, a continuum solution and finite-branch nonexistence remain
unestablished.

```sh
# Read-only assessment, at most ten aggregate CPU seconds.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_stationary_limit.py

# Root freezes these new assessment bytes before exclusive creation.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_stationary_limit.py \
  --producer-commit HEAD --write

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/assess_nsc_discovery_stationary_limit.py --check

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py -m pytest \
  tests/test_nsc_discovery_stationary_limit.py -q
```

The sole production destination is the new
`lab/results/development/nsc-discovery-stationary-limit-v1.json`; temporary
JSON output may use `/tmp`. Existing records are refused before computation.
No scientific record is overwritten.
