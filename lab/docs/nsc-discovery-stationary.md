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

## Commands and provenance

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
query, but is not invoked by diagnosis or balance selection. Its logarithmic
box remains $[-8,4]$ for $Q,r$, with $\chi\in[-10^4,10^4]$;
its normalized $10^{-7}$ label is a finite optimizer criterion only. The new
balance selector adds no physical precision gate. No dynamical holding claim,
stability result or continuum certificate is assigned to these static seeds.
