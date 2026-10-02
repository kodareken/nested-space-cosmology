# Static source-selected geometry control

This control asks whether the locked common action admits a static,
inhomogeneous finite candidate with a rank-six occupied Dirac source. It removes
the supplied constant initial $Q_0$ from this query. It does not change the
discovery family or supply another force. Long-time existence, stability,
black-hole bounces and cosmology are outside this finite control. The separate
[source-free control](nsc-discovery-vacuum-control.md) has a different covariance
and does not replace the occupied static query.

The owners are the [module](../src/recursive_horizons/nsc_discovery_stationary.py),
[driver](../scripts/derive_nsc_discovery_stationary.py), and
[focused tests](../tests/test_nsc_discovery_stationary.py). The source and gauge
come from the [common conformal action](nsc-spherical-conformal-gauge.md) and
[retained Dirac operator](nsc-nested-parent-child.md). The sealed
[first attempt](../results/development/nsc-discovery-stationary-v1.json) remains
preserved. Its necessary balance deficit motivates the higher-harmonic selector;
no stationary solution follows from a radial seed or a virial root alone.

## Static radial condition and seed

The existing induced action has

$$
F=-4\pi A r^2+2\alpha\chi,\qquad Z=-24\pi A,\qquad
V=8\pi A r^2-\alpha(4\chi+\chi^2)-2\pi C_F\mathrm{flux}^2.
$$

Its radial Euler derivative is $E_r=-8\pi A r^3R_4$. A static positive
periodic metric $g=(rQ)^2(dt^2-dx^2)-r^2d\Omega^2$ therefore requires

$$
3r_{xx}/r+(\log Q)_{xx}-Q^2=0.
$$

Integrating the equation directly, and after multiplying by $r$, gives

$$
\int Q^2dx=3\int[(\log r)_x]^2dx,\qquad
\int r_x(\log Q)_xdx=-\int Q^2r\,dx<0.
$$

Constant $Q>0$ would require $3r_{xx}=Q^2r$, impossible for positive periodic
$r$ after integration. This excludes that static ansatz; it does not exclude
oscillatory or bounded evolution. These are continuum identities. The code
reports their sampled defects and uses the original summation-by-parts (SBP)
action for every finite residual, without assuming a discrete product rule.

Choose $\widetilde u=b\cos(2\pi h x/8)$ and $Q=s\exp\widetilde u$. Then

$$
L_s=-3\partial_x^2-\widetilde u_{xx}+s^2\exp(2\widetilde u).
$$

At $s=0$, the constant trial function has Rayleigh quotient zero but is not
an eigenfunction because $\widetilde u_{xx}$ is nonconstant. The lowest
continuum eigenvalue is strictly negative. Adding the positive potential makes
it strictly increasing for $s>0$ and positive for large $s$, hence there is one
zero. Its positive ground eigenfunction supplies a radial shape. The finite
odd Fourier seed checks those signs, its positive sampled ground state and
its radial residual explicitly. It normalizes the seed shape $R$ to mean one;
this convention leaves the physical radius amplitude free. The lowest
continuum eigenfunction argument does not certify a finite or continuum full
stationary solution.

The auxiliary continuum equation gives the seed
$\chi=2(\log Q)_{xx}/Q^2-2$. The selector varies a finite nodal geometry,
interpolates $Q$ and $\chi$ through the owned carrier, and reports the resulting
fine auxiliary relation and projection defects. It never substitutes the
continuum expression for the actual fine virial source or for a discrete rate.

## Actual occupied source

Both controls use the existing period-eight periodic geometry and antiperiodic
fermion carrier: the initial joint attempt used $n_f=32,n_g=31,n_q=128$;
the bounded balance selector declares $n_f=128,n_g=127,n_q=512$.
The nested operator container uses explicit orthonormal placeholders because
its separated packet initializer is unresolved at NF32. Its full-rank canonical
frame uses one coarse harmonic and two child details; every geometry coordinate
remains present and free. The placeholders never supply a measured source.
Each seed installs the actual standing spectral columns before evaluating forces.

The source is obtained from
$H_G=U_f^\dagger H_{\mathrm{fine}}U_f$ through the actual
`nsc_nested_parent_child.hamiltonian`. The code checks that this static operator
is real symmetric and occupies its six lowest positive eigenmodes with fixed
weights $(.75,.75,.5,.5,.25,.25)$. The finite covariance
$C=\sum_an_a|v_a\rangle\langle v_a|$ obeys the CAR bounds and commutes with
$H_G$ up to measured arithmetic error. Real standing eigenvectors supply zero
current up to measured arithmetic error. No current is removed. Eigenphase
rotation leaves the covariance and source invariant.

Equal-weight pairs allow rotations within their eigenspaces. Gaps across
unequal weights, including the sixth-to-empty-mode gap, are reported. Weights
are assigned to ordered positive levels; unobserved crossings between samples
are not excluded. Equal-weight pair degeneracies do not veto a measurement.
An unresolved unequal-weight gap makes a smooth spectral branch interpretation
conditional, and prevents the original joint solver's numerical candidate label.

A full-carrier eigenfield would give
$\rho=M\sum_an_a\epsilon_a|v_a|^2/Q$, with $M=4\kappa$ once.
Retained eigenvectors need not be fine-grid eigenvectors. The actual fine
eigen-tail, quadrature $\rho_{\min}$, positive eigen-carrier expression, and
pointwise source-tail difference are retained. Pointwise or continuum positivity
is not certified from retained positive eigenvalues. The representative field
Hamiltonian has no extra copy factor; forces and energy contain $M$ once.

## Exact finite virial identity

At zero geometric momenta, let $\langle\cdot\rangle=\Delta x_q\sum$ and
$C=C_g+\rho$ be the actual pre-gauge lapse constraint. The exact SBP identity is

$$
\Delta=\alpha\langle Q^2\chi^2\rangle
-\langle Q\rho\rangle
-2\pi C_F\mathrm{flux}^2\langle Q^2\rangle
=-\tfrac12\langle r\dot p_r\rangle
-\langle\chi\dot p_\chi\rangle-\langle QC\rangle.
$$

This follows directly from the amplitude variations of the finite static
energy. Write $H_g=H_E+H_\chi+H_{\mathrm{mag}}$, where the Einstein part is
quadratic in $r$, and

$$
H_\chi=\alpha\langle Q^2(4\chi+\chi^2)\rangle
+4\alpha\langle D\chi,DQ/Q\rangle.
$$

Thus $H_E=-\tfrac12\langle r\dot p_r\rangle$ and
$-\langle\chi\dot p_\chi\rangle=H_\chi+
\alpha\langle Q^2\chi^2\rangle$. SBP gives
$\langle QC\rangle=H_E+H_\chi+H_{\mathrm{mag}}+\langle Q\rho\rangle$.
Combining them proves the stated identity without expanding $D(r^2)$ or
$D\log Q$ by a continuum product rule. The exact adjoint pullback also gives
$\Delta x_g r^T\dot p_r=\Delta x_q( A_gr)^T\dot p_{r,\mathrm{fine}}$,
and the corresponding $\chi$ and lapse pairings. The diagnostic records both
pairings and their arithmetic defects.

The source term equals the actual field energy,
$\langle Q\rho\rangle=M\sum_an_a\epsilon_a$, in this conformal standing
source. The diagnostic checks that equality independently. At stationarity,
all three residual pairings vanish and the necessary balance is

$$
\alpha\langle Q^2\chi^2\rangle
=\langle Q\rho\rangle+2\pi C_F\mathrm{flux}^2\langle Q^2\rangle.
$$

Increasing radius amplitude cannot supply this balance after radial criticality
has made $H_E=0$: the Dirac source depends on $Q$, not $r$. The selector therefore
asks whether shape and curvature variation of $Q$ supply the missing balance.
It does not replace the remaining local equations with this integrated identity.

## Preserved first attempt and fine residuals

The sealed NF32 attempt at amplitude $.35$ used 1000 optimizer function
iterations, 77,632 actual residual/Jacobian calls and about 101.68 CPU seconds.
It stopped with normalized maximum residual about $.7295$, raw
$\dot p_Q$ maximum $3.87033$ and lapse maximum $5.78929$.
The independently authenticated saved arrays have virial terms

| Term | Saved best iterate |
|---|---:|
| Auxiliary $\alpha\langle Q^2\chi^2\rangle$ | $0.9998807437$ |
| Source $\langle Q\rho\rangle$ | $12.2548273808$ |
| Magnetic | $1.7043003173$ |
| Gap $\Delta$ | $-12.9592469545$ |

The exact finite decomposition agrees within arithmetic error. Small coarse
radial and auxiliary residuals do not imply pointwise stationarity: the saved
best fine maxima are $\|\dot p_r\|_\infty\simeq0.45172235$ and
$\|\dot p_\chi\|_\infty\simeq0.00296554$. New diagnostics always report
fine momentum residuals, their differences from prolonged coarse rates,
full lapse/shift residuals and the fine auxiliary seed-relation defect.
`--diagnose` authenticates the sealed arrays, requires unchanged numerical
owners, and measures these quantities without selecting another eigensource
or running another optimizer. It can create a separate v2 diagnostic record.
The old v1 JSON and NPZ remain unchanged.

## Bounded higher-harmonic balance selector

The predeclared selector uses NF128 and harmonics $h=6,8,10,12$, with amplitude
endpoints $b=.025,.1,.35,.65$. Each trial selects its mean $Q$ through the radial
lowest-eigenvalue zero, installs its actual positive spectral source and measures
the fine virial gap. It retains every successful, failed and unadmitted endpoint.
Brent amplitude selection is admitted only across adjacent measured endpoints
whose actual gaps have opposite signs; it never brackets across a failed endpoint.
All evaluated root points and all nonbracketing data remain in the record.

One aggregate budget of at most 120 CPU seconds covers grid construction,
endpoints, root selection and optional radius fit. The smaller of five seconds
and ten percent is reserved for completion overhead. Derivative kernels and
exact $(h,b)$ trial results are cached. The first measured trial forecasts the
endpoint and worst-case costs; subsequent admission uses the largest measured
trial cost and remaining allowance. Dense cost can reduce admission and leave
an explicit incomplete bracket. The selector performs at most 64 actual trials
and 32 Brent iterations per bracket, with an additional conservative payload
admission cap. It does not start NF512, a trajectory, a new solver or another
full geometry Jacobian optimization. `--preflight-only` constructs no matrices
and evaluates no trial.

A completed root is labelled `VIRIAL_SELECTED_SHAPE_ONLY`. All remaining coarse
and fine momentum residuals, local lapse/shift constraints, fine field eigen-tail,
source minimum and projection defects are retained at every measured point.
No bracket, budget exhaustion or failed trial proves nonexistence. If several
harmonics supply roots, all are retained; the first in the declared order is the
primary readout. Without a completed root, the smallest measured absolute gap
is reported only as a diagnostic.

## Conditional radius amplitude from the finite action

The optional `--radius-select` keeps $Q,\chi,C$ fixed and sets $r=aR$.
The actual discrete $\dot p_Q$ is affine in $a^2$:

$$
\dot p_Q(a)=a^2P+B,\qquad
P=[\dot p_Q(2R)-\dot p_Q(R)]/3,\qquad B=\dot p_Q(R)-P.
$$

The code obtains these quantities from the original rates, verifies them at
$1.5R$, and checks that the source forces are unchanged. It selects the least
squares scalar $a^2=-P^TB/P^TP$ in the unweighted coarse nodal norm. A positive
value gives a conditional radius scale and its full residual profile; a
nonpositive value is a named radius-branch mismatch. It is not clamped into a
physical radius. No restoring force or dynamic prescription is introduced.
The fitted coarse/fine $\dot p_Q$, other momentum rates, full lapse/shift,
unchanged field eigen-tail and virial are recorded. This auxiliary fit can be
omitted by budget or fail without erasing the primary shape-selection result.

## Seven-coordinate local-critical successor

The preserved [v2 selector](../results/development/nsc-discovery-stationary-balance-seed-v2.json)
found virial roots for harmonics 8, 10 and 12, at amplitudes approximately
$.1641573$, $.2872187$ and $.3467144$. The first selected shape still had a
large local force: selecting only its radius amplitude left a lapse residual
about $72.55$. Necessary integrated balance is therefore an input to the next
query, not its objective.

The v3 query uses the actual NF128 common action and seven shape coordinates,

$$
Q=sS,\qquad S=1+\sum_{j=1}^7a_j\cos(8j\,2\pi x/8).
$$

All seven modes lie inside the retained geometry band. Its starting coefficients
are the truncated exponential coefficients
$2I_j(b)/I_0(b)$ at the v2 harmonic-eight root
$b=.16415733483383863$. This literal seed is part of the pinned producer. The
trigonometric polynomial retains exact cell symmetry without introducing
non-cell alias modes from exponentiating a coarse nodal field. Actual fine
$S,Q,R$ positivity is checked; no failed value is clipped into the branch.

The radial and auxiliary envelopes now use exact finite reductions. With $D$
the owned fine derivative and $U_g$ its geometry interpolation,

$$
L_Q=-3D^2+\operatorname{diag}\{Q^2-D(DQ/Q)\},\qquad
\dot p_r=16\pi A L_Qr.
$$

The operator is pulled back with the actual geometry adjoint. Since
$D(DQ/Q)=D(DS/S)$ is independent of $s$, the lowest eigenvalue zero selects
$s$; its positive ground state $R$ is normalized to mean one. The projected
auxiliary critical equation is the SPD system

$$
\operatorname{pull}[Q^2U_g\chi]
=2\operatorname{pull}[D(DQ/Q)]-2\operatorname{pull}[Q^2].
$$

No continuum $D^2\log Q$ product-rule substitution enters either solve.
Their arithmetic residuals are recorded.

This query declares a new **source preparation branch**. Let
$E_{\mathrm{base}}=M\sum_an_a^{\mathrm{base}}\epsilon_a$ and
$\mathcal W=\alpha\langle Q^2\chi^2\rangle$. Select

$$
\eta=(\mathcal W-H_{\mathrm{mag}})/E_{\mathrm{base}},\qquad
c=\eta c_{\mathrm{base}},\qquad 0<\eta\le4/3.
$$

The resulting six occupations and total trace $3\eta$ are explicit.
The upper bound is the CAR bound for the largest base occupation $.75$.
This selects a different initial covariance when $\eta\ne1$, with all action
coefficients unchanged. It is not a source pump, reset, or new evolution law.
Every physical force and Hellmann--Feynman variation **holds its selected
$\eta$ fixed**. Differentiating the preparation function $\eta(a)$ is a
separate reduced calculation.

Obtain $P_R$ and the remainder from actual $\dot p_Q$ rates at $R$ and $2R$,
with this selected source. The finite normal identity is

$$
\langle Q,P_R\rangle=16\pi A\langle Q^2R^2\rangle>0.
$$

For $\delta S_j=\cos(8j\,2\pi x/8)$, define

$$
\delta Q_{0,j}=s\delta S_j,\qquad
\mu_j=\langle\delta Q_{0,j},P_R\rangle/\langle Q,P_R\rangle,
\qquad \delta Q_j=\delta Q_{0,j}-\mu_jQ.
$$

These variations are tangent to the radial-critical envelope. The auxiliary
critical variation and the fixed-weight spectral Hellmann--Feynman identity
give the analytic reduced gradient and normal radius scale,

$$
\partial_{a_j}\eta=\langle\delta Q_j,\mathrm{remainder}\rangle/E_{\mathrm{base}},
\qquad r^2_{\mathrm{amplitude}}
=-\langle Q,\mathrm{remainder}\rangle/\langle Q,P_R\rangle.
$$

The implementation checks the normal identity, tangent pairings, radius-source
independence, and affinity at an independent $1.5R$. It compares the analytic
$\eta$ gradient with central differences at two step sizes. Reduced NF64 tests
also independently check the spectral energy variation with $\eta$ held fixed.

The optional execution solves the seven analytic gradients using a small
Newton/LM loop. Finite differences form the seven-by-seven gradient Jacobian;
there is no full geometry finite-difference Jacobian. Backtracking rejects and
records nonpositive fine shapes/ground states, nonpositive normal radius squares
and source preparations outside the CAR interval. It never clips those values.
The loop permits at most 40 iterations and 600 actual trials, with the same
120 aggregate CPU-second ceiling, measured-cost admission and completion
reserve. A bounded cache keeps sixteen full trial states. Every trial retains
its coefficients, preparation, gradient and raw residual summaries; the initial
and best full states, all their raw fine/coarse residuals, and source carriers
are saved in the payload. Budget or line-search termination preserves the best
admissible finite result and does not imply nonexistence.

Covariance cell-translation and Dirac-reflection gaps are measured directly.
No group average is applied to the covariance or force. Reflection uses the
antiperiodic spatial reflection together with $\sigma_1$, which preserves the
even static Dirac operator. Even if seven gradients and the normal pairing
vanish, the full retained $\dot p_Q,\dot p_r,\dot p_\chi$, lapse and shift
residuals are measured, as are all fine residuals and the fine field eigen-tail.
The label `CRITICAL_RETAINED_CANDIDATE` requires both the seven-gradient and full
retained numerical residual criteria, with a resolved unequal-occupation branch.
It is a conditional finite-domain label, not holding, stability or a continuum
stationary solution. A restriction to seven cosines alone does not certify
unmeasured force directions.

`--critical-shape` defaults to preflight with zero numerical trials.
`--execute-critical` explicitly enables the bounded query. Its fresh default
stem is `nsc-discovery-stationary-critical-v3`. Earlier v1/v2 evidence and their
immutable historical producer checks remain preserved.

## Commands and provenance

### Executed local-critical query

The [v3 record](../results/development/nsc-discovery-stationary-critical-v3.json)
is produced by `697efe0` and preserves all 600 evaluations, using 31.26 CPU
seconds. Its verdict is `UNSATISFIED_LOCAL_CRITICAL_RELATIONS`. The selected
trial 579 has mean $Q=0.00094539$, mean areal radius 1287.02 and occupation
scale $\eta=1.21871$, within the CAR domain. The shape approaches uniformity
and the normal radius denominator falls from 3.2103 to
$1.6180\times10^{-5}$. Retained $p_Q$ remains 32.96 and the full lapse maximum
35.07. A small reduced shape gradient therefore does not establish local
criticality. The areal-radius growth is physical, rather than merely a clock
normalization. This noncompact tendency ends the present static search;
additional iterations are not justified by these results. It is not a
class-wide obstruction or a dynamical stability test.

Default execution only measures the NF32 radial seed. Root freezes the exact
producer bytes before executing a scientific selector. These commands run from
the repository root with the existing validation interpreter:

```sh
# No optimizer or record write.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py

# Authenticate preserved evidence against its immutable producing commit.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --check \
  --record lab/results/development/nsc-discovery-stationary-v1

# Forecast the new selector without any numerical trials.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --balance-seed --preflight-only

# After root's producer freeze: new diagnostic, no eigensource selection/solve.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --diagnose \
  --record lab/results/development/nsc-discovery-stationary-v1 \
  --producer-commit HEAD --write

# After root's producer freeze: one bounded discriminating selection.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --balance-seed --radius-select \
  --cpu-limit 120 --producer-commit HEAD \
  --write lab/results/development/nsc-discovery-stationary-balance-seed-v2

# New critical query: preflight only, no numerical state evaluation.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --critical-shape --preflight-only

# After root's producer freeze: seven-coordinate local-critical execution.
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --critical-shape --execute-critical \
  --cpu-limit 120 --producer-commit HEAD \
  --write lab/results/development/nsc-discovery-stationary-critical-v3

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py -m pytest \
  tests/test_nsc_discovery_stationary.py -q
```

`--write` creates only new JSON/NPZ pairs named `nsc-discovery-stationary-*`
directly under `lab/results/development`, or a specified `/tmp` stem.
The default new stem is `nsc-discovery-stationary-v2`. An existing JSON or NPZ
is refused before compute. Records and payloads are each capped at 16 MiB;
per-trial payloads retain all raw field, source and geometry measurements while
omitting repeated dense $H$ and $C$ matrices. Their measured finite spectral gaps
are recorded, and the source columns/geometry preserve reconstruction inputs.

A producing commit is accepted only when every pinned producer file matches
it byte for byte. New measurement checks pin producer and locked input hashes
before and after computation. Read-only replay authenticates sealed producer
and input hashes against the recorded immutable commit, so newer code/docs/tests
do not invalidate v1. An unpinned temporary record instead requires exact current
producer bytes. `--check` validates saved geometry, CAR columns, positive levels,
raw residual summaries, virial pairings and selected roots without solving.
It reports which producer authentication was used.

The historical `--solve` option remains available for the original NF32 joint
query, but is not invoked by diagnosis, balance selection or the critical query. Its logarithmic
box remains $[-8,4]$ for $Q,r$, with $\chi\in[-10^4,10^4]$;
its normalized $10^{-7}$ label is a finite optimizer criterion only. The new
balance selector adds no physical precision gate. No dynamical holding claim,
stability result or continuum certificate is assigned to these static seeds.
