# Static source-selected geometry control

This control asks whether the same locked common action admits a static,
inhomogeneous finite candidate with a rank-six occupied Dirac source. It removes
the earlier supplied constant initial $Q_0$ from this query. It does not change
the discovery family, run a trajectory, introduce another force, or assert a
black-hole bounce, cosmology, long-time existence, or stability. The separate
[source-free control](nsc-discovery-vacuum-control.md) addresses a different
covariance and does not replace the occupied static query.

The owners are the [module](../src/recursive_horizons/nsc_discovery_stationary.py),
[driver](../scripts/derive_nsc_discovery_stationary.py), and
[focused tests](../tests/test_nsc_discovery_stationary.py). The source and gauge
come from the [common conformal action](nsc-spherical-conformal-gauge.md) and
[retained Dirac operator](nsc-nested-parent-child.md).

## Derived static radial condition

For the existing induced action,

$$
F=-4\pi A r^2+2\alpha\chi,\qquad Z=-24\pi A,\qquad
V=8\pi A r^2-\alpha(4\chi+\chi^2)-2\pi C_F\mathrm{flux}^2.
$$

The covariant radial Euler derivative is $E_r=-8\pi A r^3R_4$. A static
positive periodic metric $g=(rQ)^2(dt^2-dx^2)-r^2d\Omega^2$ therefore requires

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
reports their sampled defects and uses the original summation-by-parts action
for the finite solver, without assuming a discrete product rule.

## Radial seed and occupied source

Choose a nonconstant cosine $\widetilde u$ and $Q=s\exp\widetilde u$. Then

$$
L_s=-3\partial_x^2-\widetilde u_{xx}+s^2\exp(2\widetilde u).
$$

At $s=0$, the constant trial function has Rayleigh quotient zero, but is not
an eigenfunction because $\widetilde u_{xx}$ is nonconstant. The lowest
continuum eigenvalue is strictly negative. The positive added potential makes
it strictly increasing for $s>0$ and positive for sufficiently large $s$;
there is one zero. Its positive ground eigenfunction supplies a radial shape.
The finite odd Fourier seed checks those signs, its positive sampled ground
state and its radial residual explicitly. It normalizes the seed shape to
mean one. That normalization does not fix the physical radius in the solve.
The seed solves only this radial equation, not the remaining stationary
balance. The strict ground-state monotonicity argument belongs to the
continuum operator; finite seed positivity and eigenvalue signs are measured.

The query uses $n_f=32$, $n_g=31$, and $n_q=128$ on the existing period-eight
periodic geometry and antiperiodic fermion carrier. The nested operator container
uses explicit orthonormal placeholder columns because its separated packet
initializer is unresolved at this coarse band. Its full-rank canonical frame
uses one coarse harmonic and two child details to fit the NF32 collar; all
31 geometry coordinates remain present and free. The placeholders never supply
a measured source: every residual evaluation installs the standing spectral
columns before evaluating forces. It obtains
$H_G=U_f^\dagger H_{\mathrm{fine}}U_f$ through the actual
`nsc_nested_parent_child.hamiltonian`, checks that this static operator is real
symmetric, and occupies its six lowest positive eigenmodes with fixed weights
$(.75,.75,.5,.5,.25,.25)$. Its covariance
$C=\sum_an_a|v_a\rangle\langle v_a|$ obeys the finite CAR bounds and commutes
with $H_G$ up to measured arithmetic error. Real standing eigenvectors supply
zero current up to measured arithmetic error. No current is removed. Eigenphase
rotation leaves the covariance and its source invariant.

Equal-weight pairs allow rotations inside their eigenspaces. Gaps across
unequal weights, including the sixth-to-empty-mode gap, are reported. A
numerically unresolved crossing prevents the candidate label and makes a
smooth Hellmann--Feynman branch interpretation conditional. The optimizer
records its minimum sampled unequal-weight gap and unresolved-gap evaluation
count. It assigns weights to ordered positive levels, and does not claim to
exclude crossings between sampled evaluations. The occupied
positive eigenvalues and all retained eigenvalues are saved.

A full-carrier eigenfield would give
$\rho=M\sum_an_a\epsilon_a|v_a|^2/Q$ with $M=4\kappa$ once. Retained
eigenvectors need not be fine-grid eigenvectors. The actual fine eigen-tail,
quadrature $\rho_{\min}$, positive eigen-carrier expression, and their
pointwise difference are saved. No continuum or pointwise positivity
certificate is promoted from retained positive eigenvalues.

## One bounded finite solve

The unknowns are all nodal $\log Q$, $\log r$, and $\chi$. Both physical
means and the radius amplitude are free. Momenta are zero; a single sine
translation anchor selects a phase. The residual consists of the original
projected $\dot p_Q,\dot p_r,\dot p_\chi$ and the original lapse-plus-matter
and shift-plus-current constraints. The conformal lapse chain rule is included
by the owning rate function. The full quadrature constraints are also reported.
No continuum product expansion replaces a discrete residual.

The least-squares scales are frozen RMS values at the radial seed, floored at
$10^{-3}$ to balance unlike units. Every raw residual is retained. The
computational box is $-8\le\log Q,\log r\le4$,
$-10^4\le\chi\le10^4$; it limits this attempt and cannot establish global
nonexistence. The fine interpolated positive chart is checked by the original
rate owner. The numerical candidate target is a normalized maximum residual
$10^{-7}$ with resolved unequal-weight gaps. This is a finite optimizer
criterion, not a physical error certificate. An unsatisfied attempt names its
largest remaining normalized relation and preserves every residual. An
optimizer stop never becomes a nonexistence result.

The representative field Hamiltonian has no extra copy factor. The source
forces and field energy contain $M$ once. A fixed-direction finite difference
of $M\sum_an_a\epsilon_a$ is compared with the owned full conformal source
variation $\sum(F_Q+F_L)\delta Q$ (Hellmann--Feynman). Its branch-gap caveat is
retained. Long-time evolution, dynamical stability and continuum certification
remain unperformed.

## Preview, bounded execution and replay

The default command only prepares and measures the seed; it does not optimize
or write a result. `--solve` explicitly allows one coarse attempt. The default
CPU ceiling is 120 seconds. The solver reserves the smaller of two seconds and
ten percent of the budget for fixed-size final diagnostics; the
actual CPU cost and termination are recorded. No larger grid or refinement
campaign is started by this driver.

```sh
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --solve --cpu-limit 120 \
  --max-nfev 1000 --producer-commit HEAD --write

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py \
  scripts/derive_nsc_discovery_stationary.py --check

OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py -m pytest \
  tests/test_nsc_discovery_stationary.py -q
```

Root must freeze the exact producer bytes before executing the bounded solve.
The optional producing commit is accepted only when every pinned producer
file matches that commit byte for byte. Without that pin, the record explicitly
stores a null producing commit. Producer and locked-input hashes are checked
before and after measurement. `--write` creates only the new
`lab/results/development/nsc-discovery-stationary-v1.json` and `.npz`, or a
specified `/tmp` stem. Existing JSON or NPZ files are refused before compute.
`--check --record /tmp/example` authenticates a temporary result. Replay is
read only: it validates producer/input/payload hashes and the saved geometry,
CAR columns, positive eigenmodes and residual summaries without rerunning the
solver. The focused tests use small seed/finite algebra checks and do not run
a production solve.
