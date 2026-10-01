# Same-action conformal gauge and constraint forcing

The optional `gauge="conformal"` in
[`nsc_spherical_galerkin_coupling.py`](../src/recursive_horizons/nsc_spherical_galerkin_coupling.py)
uses $L=Q$, $\beta=0$ at every stage. The default remains `prescribed` and
replays the calibrated gauge. This is a formulation of the existing
$\Gamma_{\mathrm{one}}$ coupling, with its original coefficients and source;
it introduces no prescribed radius, source pulse, restoring force, filtering,
or constraint reset. It does not establish renewal or a continuum Cauchy
solution. No new physical trajectory is recorded by this note.

## Constraint variation and gauge ordering

First retain lapse and shift variations of the action in
[the first-order chart](nsc-spherical-feedback-action.md) and
[the direct source](nsc-conformal-adm-source.md). They give

$$
C=C_g+\rho,\qquad D=D_g+j,
$$

$$
C_g=\frac{p_Qp_\chi}{2f}+\frac{\Pi^2}{4ZQ}
 +\frac{Zr_x^2}{Q}-QV-2\partial_x(F_x/Q),
\qquad
D_g=p_r r_x+p_\chi\chi_x-Q\partial_xp_Q,
$$

where $f=F_\chi\ne0$, $\Pi=p_r-F_rp_\chi/f$, and the unchanged owned
coefficients obey

$$
F=-4\pi A r^2+f\chi,\quad Z=-24\pi A,
\quad V=8\pi A r^2-\frac f2(4\chi+\chi^2)-2\pi C_F\mathrm{flux}^2.
$$

For canonical half-density columns, let $K$, $S_1$, $P_0$ denote the
weighted bilinears owned by `column_moments`. In continuum notation their
sampled densities are $K_c=K/\Delta x$, $S_c=\operatorname{Re}S_1/\Delta x$,
$P_c=P_0/\Delta x$. With the isotropic copy weight $M=4\kappa$ applied once,

$$
\rho=M(K_c/Q+\kappa S_c),\qquad j=-MP_c.
$$

Only then choose $L=Q$, $\beta=0$, so
$h_{ab}=Q^2\operatorname{diag}(1,-1)$ and the gauge-fixed Hamiltonian is

$$
H_c=\int h_c\,dx,\qquad h_c=QC.
$$

The constraints $C=D=0$ are retained. A variation of $H_c$ differentiates
its $Q$-dependent lapse. Compared with evaluating the original fixed-lapse
rates at $L=Q$, it adds $-C_g$ to $\dot p_Q$ and the matter contribution is
$-(F_Q+F_L)/\Delta x$. Omitting $-C$ is a distinct extension away from the
constraint surface. The extra term vanishes on that surface; it is a
constraint term derived from this same gauge-fixed action, not an additional
matter force. Both retained constraints have exact homogeneous continuum
transport under the full gauge-fixed equations below.

## Equations and continuum transport

Write $a=F_{rr}=-8\pi A$, so $F_r=ar$ and $Z=3a$. The continuum equations are

$$
\dot Q=\frac{Qp_\chi}{2f},\qquad
\dot r=\frac\Pi{2Z},\qquad
\dot\chi=\frac{Qp_Q}{2f}-\frac{F_r\Pi}{2fZ},
$$

$$
\dot p_Q=-\frac{p_Qp_\chi}{2f}+2QV+\frac{2F_{xx}}Q-M\kappa S_c,
$$

$$
\dot p_r=\frac{\Pi a p_\chi}{2Zf}+2Zr_{xx}+Q^2V_r
 +F_r\partial_x(2Q_x/Q),
$$

$$
\dot p_\chi=Q^2V_\chi+f\partial_x(2Q_x/Q),
\qquad
 i\dot\Phi=(\sigma_2P+\kappa Q\sigma_1)\Phi,
\quad P=-i\partial_x.
$$

The field evolves by the representative block; $M$ and the occupation
weights scale its geometric force and energy, not its generator. Its
high-frequency kinetic coefficient is one. The retained-band stability
frequency is bounded by $k_{\max}+\kappa\max Q$; geometric nonlinearity and
the integration error still require their own checks.

For these equations, direct differentiation gives

$$
\partial_t h_c=\partial_xD,\qquad
\partial_tD=\partial_xh_c.\tag{1}
$$

This is an exact continuum identity on smooth periodic fields in the
positive chart. In particular, $C=D=0$ remains zero; positivity allows
$h_c=0$ to imply $C=0$. The proof includes the angular mass term and the
source force. In real and imaginary components
$\Phi=(u+iv,w+iz)^T$, one weighted column has

$$
K_c=-u w_x-v z_x+w u_x+z v_x,\quad
S_c=2(uw+vz),\quad
P_c=u v_x-v u_x+w z_x-z w_x.
$$

The field equations are
$\dot u=-z_x+\kappa Qz$, $\dot v=w_x-\kappa Qw$,
$\dot w=v_x+\kappa Qv$, $\dot z=-u_x-\kappa Qu$.
Differentiating $h_c$ and $D$ with these and the six geometric equations
cancels the mass and geometric source terms in (1). The independent symbolic
test derives all ten rates as Euler derivatives of $h_c$, with the column's
weighted symplectic factor $1/(2M\omega)$, and simplifies both residuals to
zero for free spatial jets. Arbitrary column weights add linearly.

## Actual finite forcing and the energy estimate

The implementation retains the existing summation-by-parts energy and
exact adjoint rate pullback. It uses $\partial_xF$ as `derivative @ F`,
including the $-C_g$ correction as that sampled constraint. It does **not**
replace it by a product-rule expansion such as $2F_{xx}/Q$ on the finite
grid. That continuum simplification need not equal the discrete derivative.
The finite projected model therefore need not satisfy (1) exactly.

`conformal_constraint_transport` differentiates the actual sampled
constraints analytically in the direction of the actual lifted Galerkin
rate. It returns

$$
f_h=\dot h_c-\mathsf D_x D,\qquad
f_D=\dot D-\mathsf D_x h_c.
$$

The same calculation in the unprojected fine rate gives the quadrature
forcing. The difference between the actual and unprojected derivatives is
the rate-projection forcing. Their vector sum is the actual forcing; norms
of those pieces do not generally add as equality. All three are exposed.
Smooth low-band data have small forcing; a top geometry mode produces a
nonzero projection defect. That defect is not called constraint preservation.

Set

$$
R=\left[\Delta x\sum(h_c^2+D^2)\right]^{1/2},\quad
F=\left[\Delta x\sum(f_h^2+f_D^2)\right]^{1/2}.
$$

Because the periodic derivative is skew-adjoint,

$$
\frac d{dt}R^2=2\Delta x\sum(h_c f_h+D f_D),
\qquad |R'|\le F,
$$

with the integrated bound

$$
R(t)\le R(0)+\int_0^t F(s)\,ds.\tag{2}
$$

The inequality at $R=0$ follows by regularizing the norm and taking the
limit. There is no exponential in $k_{\max}$ in (2). The identity is exact
for the finite ODE and sampled norm; applying the continuum version requires
controlling quadrature/product defects and the temporal residual of the
chosen reconstruction. An RK4 endpoint sequence and a sampled trapezoid
integral do not by themselves certify that integral. The API accordingly
returns `continuum_constraint_certified=false` and
`time_integral_certified=false`.

A bound for $R$ controls $C$ by $\|C\|_2\le\|h_c\|_2/\min Q$ in the
positive chart. A total state-error bound or an error bound on a regenerated
local observable needs its additional dynamical stability and measurement
map. Neither follows from (2) alone.

## Initial compatibility, work, and clocks

$F_L$ and $F_\beta$ are independent of $L,\beta$ in this action. Thus the
same initial $Q,r,\chi$, momenta and Gaussian columns have exactly the same
lapse and shift constraints after selecting the conformal gauge. Canonical
momenta and the initial intrinsic and normal-extrinsic data are unchanged;
coordinate derivatives and the coordinate clock change. This allows reuse
of prepared initial data without adding a new source or solving a different
initial-radius model.

The system stored in `grid.fine` remains the calibrated carrier.
`active_fine_system(grid,fine_state)` returns a new stage-local conformal
system. No shared grid is mutated. The rate bundle exposes that system,
the actual lifted `lapse_dot`, and `shift_dot=0`. For the conformal gauge,
$\dot L=\dot Q_{\mathrm{lifted}}$, so the field work is

$$
\dot E_{\mathrm{field}}=
\sum(F_Q\dot Q+F_L\dot L+F_\beta\dot\beta)
=\sum(F_Q+F_L)\dot Q_{\mathrm{lifted}}.
$$

This is the existing `rate.fieldwork_power`. A consumer using only
$F_Q\dot Q$ does not apply to this dynamic lapse. The force is a nodal
energy derivative: this work sum has no additional $\Delta x$ or $M$.
An independent finite difference of the full field energy checks the
identity, including column evolution.

New runs must record `gauge="conformal"` explicitly. Historical requests
without a gauge select `prescribed`; existing records remain in their
original gauge. Comparing responses between gauges requires a stated
observer and clock protocol. The conformal normal observer has
$d\tau=rQ\,dt$ and zero coordinate shift. Coordinate-time endpoints cannot
be equated with the recorded prescribed-gauge endpoints without that
protocol, nor can old fixed-lapse curvature helpers be used unchanged.

## Reproduction

The owned independent checks cover the exact continuum identities, all
Hamiltonian chain-rule derivatives, dynamic work, the constraint directional
derivative, forcing decomposition, initial-data compatibility, rejection of
unsupported gauge/dimensions/chart, and the unchanged prescribed rate.

```sh
OPENBLAS_NUM_THREADS=1 .venv/validation/bin/python scripts/lab.py -m pytest \
  tests/test_nsc_spherical_conformal_gauge.py \
  tests/test_nsc_spherical_galerkin_coupling.py -q
```

No expensive production run or scientific record is produced by these tests.
