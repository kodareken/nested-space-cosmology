# Unresolved vortices can drive a resolved circulation

This entry connects the [flow dictionary](nsc-atlas-flow-dictionary.md) to an
explicit nonlinear momentum-transfer calculation. The source is OpenAI's
[Finite time blowup for Navier--Stokes](https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf),
released with an [official announcement](https://openai.com/index/navier-stokes-solution/)
and [Lean formalization](https://github.com/openai/NavierStokesAndEuler).
These are attributed OpenAI research results, not discoveries made in this chat.

## Preparation and decision

The source's Theorem 1.1 concerns three-dimensional incompressible flow with
positive viscosity, zero initial velocity, a smooth compactly supported force,
bounded kinetic energy and finite-time velocity blowup. This entry does not
independently verify its entire 166-page argument or rebuild its Lean proof.
The statement, core scalings and selected mean-stress mechanisms are read and
reused under their stated hypotheses.

Three routes were considered: rebuilding the complete singularity proof,
continuing only the prior linear gravitational pulse, or extracting and
independently verifying a finite nonlinear transport mechanism. The third
is chosen. It directly addresses how fine circulating structure affects a
larger resolved flow, without requiring an imposed reflector or singular state.
The [earlier gravitational result](nsc-atlas-gravitational-response.md) remains
complete in its separate linear sector.

Executed route:

1. Construct compact C-infinity divergence-free velocities from a cylindrical
   vector potential with angular mode m and radial phase proportional to m log r.
2. Average the exact nonlinear momentum flux. Prove a nonzero initial mean
   tangential acceleration for an unforced Navier--Stokes evolution, even
   though the initial angular mean velocity is zero.
3. Fix the envelopes, phase ratio -2, modes 4/8/16, radii .75/1.25 and viscosity
   .001 before evaluation. Calculate energy, angular/axial momentum flux and
   viscous cost with controlled quadrature, keeping total momentum balance.
4. Link the exact finite calculation to the source's concentrating-core
   scalings and to the boundary-traction interface in the existing atlas.
   Record which further dynamical/source identifications remain open.

The variables are two fixed envelopes, an angular mode and a phase ratio;
the other vector-potential component is fixed at minus one-half of the first.
Density and length/velocity units are supplied. Pressure is not fitted: its
angular mean has no tangential gradient, which removes it from the predicted
initial tangential rate. The subsequent solution and its lifespan are not
supplied by this initial-rate calculation.

Allow 60 CPU seconds and 500000 quadrature/sample points for this finite task,
with useful independent work in up to four workers. No full 3D evolution grid
or search over a blowup parameter is needed. Check solenoidality symbolically,
the mean force by direct angular advection, and energy/dissipation by
independent formulas and quadrature. Preserve only compact coefficients and
plot data.

Completion requires a verified nonzero rate, its conservation/energy budget,
source attribution and an integrated atlas entry. Likely failures are a lost
cylindrical basis term, confusing zero mean with zero full velocity, or
mistaking a fixed-time transfer for a sustained cascade. Each has a separate
check and scope statement. A direct short-time Navier--Stokes evolution that
fails the predicted initial rate under the same preparation would falsify
the new finite calculation.

## Source identity and reading boundary

The public PDF has 166 pages and the author line OpenAI. The theorem and
physical description match the supplied extract. The supplied UTF-8 text has
523447 bytes and SHA256
`46c37997b150309499f6c4383ad0df2d7b0ac91fdd9ba9b52ceb4a3296af36fd`.
The formalization was inspected at
[f9e8bc5b38b6e212696e8a30e3e91517af887bbd](https://github.com/openai/NavierStokesAndEuler/tree/f9e8bc5b38b6e212696e8a30e3e91517af887bbd),
dated 10 September 2026. Its entry files and metadata were read; Lean was not
rebuilt and its dependencies/axioms were not independently audited here.

Direct reading covers Theorem 1.1, Sections 2.1--2.2 and selected 3.2--3.3 stress
relations. Grok independently inspected selected Proposition 7.5/Appendix C
statements and the formalization's source metadata. Profile existence,
all-order residual cancellation and the full singularity proof were not
independently verified in this task. The reference supplies those results;
the finite mechanism below has its own derivation and controls.

The whole-space source problem starts from zero full velocity and uses a
smooth compact force. It is not an unforced zero-data blowup. Our finite
example below instead starts from specified nonzero data and sets the force
to zero after release. Its initial mean velocity, not its full velocity, is zero.

## A compact exact preparation

Use density 1 and dimensionless reference length/speed 1. Let

$$
g(s)=\begin{cases}\exp(1-1/(1-s^2)),&|s|<1,\\0,&|s|\geq1,\end{cases}
\quad a(r)=g(2(r-1)),\quad b(r)=-a(r)/2,\quad \chi(z)=g(z).
\tag{T1}
$$

For an integer m>=1 set kappa=-2 and phi=m theta+kappa m log r. Define

$$
\mathbf A=\frac{\chi(z)\cos\phi}{m}
\left(0,b(r),a(r)\right),\qquad u_0=\nabla\times\mathbf A.
\tag{T2}
$$

The support is .5<=r<=1.5, |z|<=1. All derivatives vanish at its edges;
the support stays away from the cylindrical axis. Extension by zero is smooth
on all of R3 and div(u0)=0 exactly. It is admissible initial data for a local
classical incompressible Navier--Stokes solution at each positive viscosity.

Write u0=Re(W exp(i phi)). Taking the cylindrical curl gives

$$
W_r=-\frac{b\chi'}m+i\frac{a\chi}r,\quad
W_\theta=-\frac{a'\chi}m-i\kappa\frac{a\chi}r,\quad
W_z=\frac{(b'+b/r)\chi}m+i\kappa\frac{b\chi}r.
\tag{T3}
$$

The mean is the rotation-group average, with the cylindrical basis carried
by the rotation. Every component has zero mean. The covariance does not:

$$
R_{ij}=\langle u_{0i}u_{0j}\rangle
=\frac12\operatorname{Re}(W_i\overline{W_j})
=R^{(0)}_{ij}+\frac{R^{(2)}_{ij}}{m^2}.
\tag{T4}
$$

Both terms are positive semidefinite rank-one covariances: R0 comes from the
imaginary part of W, R2 from its real part before division by m. In particular,

$$
R_{r\theta}=\frac{b a'\chi\chi'}{2m^2}
-\frac{\kappa a^2\chi^2}{2r^2},\qquad
R_{z\theta}=-\frac{a'(b'+b/r)\chi^2}{2m^2}
-\frac{\kappa^2ab\chi^2}{2r^2}.
\tag{T5}
$$

As m grows the prepared velocities converge weakly to zero while R0 remains.
This is a concrete nonlinear averaging defect, not a claim that kinetic
energy disappears or that a uniform time interval survives this limit.
For any positive covariance, half its trace is at least
sqrt(R_rtheta^2+R_rz^2). A specified momentum flux therefore has a nonzero
energy cost even when the mean velocity vanishes.

## Initial-transfer theorem for the actual unforced equation

Let U be the rotational mean of the local classical solution of
u_t+div(u tensor u)-nu Delta u+grad p=0, with initial data T2. Since U(0)=0,
the mean viscous term vanishes initially. The mean pressure is axisymmetric,
so its theta gradient is zero. Averaging the theta momentum equation gives

$$
\left.\partial_t U_\theta\right|_0
=-\left(\partial_rR_{r\theta}+\frac{2R_{r\theta}}r
+\partial_zR_{z\theta}\right).
\tag{T6}
$$

The factor2/r includes the cylindrical basis term. This is the initial
derivative of a Navier--Stokes evolution, not a velocity field made to satisfy
the equation by defining a residual force. At z=0, chi=1 and chi'=0, while
R_ztheta is even in z. Substitution into T6 proves, for every m and nu>0,

$$
\left.\partial_tU_\theta(r,0,t)\right|_{t=0}
=\kappa\frac{a a'}{r^2}=-\frac{2a a'}{r^2}.
\tag{T7}
$$

At the two fixed radii the rates are exactly

$$
A_{\rm inner}=-\frac{1024}{81}e^{-2/3}=-6.49060654184413,\qquad
A_{\rm outer}=\frac{1024}{225}e^{-2/3}=2.33661835506389.
\tag{T8}
$$

Thus finer prepared angular structure drives the same nonzero initial
resolved circulation at these locations. The signs produce opposite rotating
bands. With local circulation Gamma=2pi r Utheta, T7 also fixes its initial
rate. No sustained inverse cascade or finite-time blowup follows from this
initial-rate theorem.

### A parent flow and composition of resolutions

The result also gives an exact incremental law around a prepared parent flow.
Let B be any smooth axisymmetric divergence-free background with admissible
decay, and compare local Navier--Stokes evolutions from B and B+epsilon u0
with the same force and viscosity. Rotational averaging kills both cross
advection terms because B is invariant under the averaging group and <u0>=0.
The added viscous mean is zero and both mean pressures have zero theta
gradient. Therefore

$$
\delta\left(\partial_tU_\theta\right)_{t=0,z=0}
=\epsilon^2\left(-\frac{2aa'}{r^2}\right).
\tag{T8a}
$$

For any finite collection of distinct angular modes, all cross covariances
also average to zero. With initial data B+sum_j epsilon_j u_mj the increment
is T8a with epsilon^2 replaced by sum_j epsilon_j^2. The added kinetic energy
and initial dissipation are respectively sum epsilon_j^2 E_mj and
sum epsilon_j^2 D_mj. This is a concrete nonlinear composition law for
prepared resolutions, with its energy and momentum budget retained.
It is not an arbitrary multiplication by the number of levels: since
E_m>=E0, a fixed added-energy budget bounds sum epsilon_j^2 and the
mean-force increment. Extending it to infinitely many modes needs separate
regularity and dissipation summability assumptions.

## Momentum budget and viscous cost

Let I=integral chi^2 dz=.983380812912727. The outward angular-momentum and
axial-momentum fluxes through a cylinder are

$$
\mathcal J_L(r)=2\pi r^2\int R_{r\theta}dz
=-\pi\kappa a^2 I=2\pi a^2 I,
\tag{T9}
$$

$$
\mathcal J_z(r)=2\pi r\int R_{rz}dz
=\frac{\pi\kappa abI}{r}=\frac{\pi a^2I}{r}.
\tag{T10}
$$

Terms containing chi chi' integrate to zero. Both fluxes vanish at the
radial support edges. Height-integrated initial mean momentum changes
are minus their radial derivatives; the global initial torque and axial
momentum-change budgets sum to zero. This is redistribution within the
prepared field, not an external source of total angular momentum.

Direct integration and integration by parts give

$$
E_m=E_0+\frac{E_2}{m^2},\qquad
D_m=\nu\left(D_2m^2+D_0+\frac{D_{-2}}{m^2}\right).
\tag{T11}
$$

Here E=(1/2) integral |u0|^2 dV and D=nu integral |grad u0|^2 dV. The
independently computed curl loss is identical because the field is compact
and solenoidal. Coefficients are
E0=4.69732546551740, E2=12.4674608205504,
D2=54.9521280290796, D0=679.082984142553 and Dminus2=2628.54761911414.
For nu=.001,

| Angular mode m | Kinetic energy | Initial dissipation rate |
| --- | --- | --- |
| 4 | 5.476541766802 | 1.722601258802 |
| 8 | 4.892129540839 | 4.237090234552 |
| 16 | 4.746026484348 | 14.757095523724 |

The unforced energy identity is E'(0)=-D_m. Finer angular structure retains
the initial mean force but pays an increasing viscous cost. A steady supply,
long-lived cascade or m-independent lifespan has not been established.

The [compact calculation record](../results/summaries/nsc-atlas-vortex-transfer-v1.json)
retains 64/128/256 quadrature controls, direct angular advection, independent
curl/gradient loss and source identities. Four workers evaluate 258944 points,
using about .44 solver CPU seconds. There is no time evolution grid.
Fourteen new controls check the construction and conservation laws. An
[additional direct parent/composition check](../tests/test_nsc_atlas_vortex_composition.py)
verifies the nonlinear rate and covariance of three simultaneous modes on
a smooth axisymmetric parent vortex.

For an explicitly supplied laboratory-sized example L=1 m, U=.1 m/s,
rho=1000 kg/m^3 and nu=1e-4m^2/s, the predicted initial accelerations are
-.0649060654m/s^2 and+.0233661836m/s^2. The torque flux at r=L is
61.7876388N m and the axial momentum flux 30.8938194 N. These are unit-scaled
conditional outputs of the stated preparation, not experimental measurements.

## How the OpenAI result joins the map

The source's concentrating core has radial scale tau^(1/2), axial scale
tau^(1/2-h) and azimuthal/axial velocity scale tau^(-1/2-h), with 0<h<.01.
The independently checked power counting gives core energy tau^(1/2-3h),
dominant radial-gradient loss scale tau^(-1/2-3h), and core L3 norm scale
tau^(-4h/3). Energy tends to zero and the loss remains time-integrable,
although pointwise velocity and the critical norm grow. Those scalings are
consistent with the stated energy identity; they do not replace its global
construction or prove that the source's force has the required regularity.

The reusable connecting equation is the nonlinear momentum correlation
R=<u tensor u>-U tensor U. The atlas's coarse state must retain it if that
flow sector is to predict momentum transfer. The tensor supplies a candidate
material traction for the existing Brown--York/Codazzi balance once its
relativistic stress and geometry have been specified. Our earlier linear
gravitational response cannot itself generate this nonlinear mean transfer.

There are two separate gravitational joints. In Bredberg et al.'s construction,
a fluid with p spatial dimensions has a p+2-dimensional vacuum-gravity lift:
the source's 3D fluid corresponds to 5D gravity in that particular map. The
previous 4D Rindler example used a 2D fluid. Also the singular endpoint leaves
any fixed small-velocity/long-wavelength regime, so its hydrodynamic lift
cannot simply be extrapolated to that endpoint. A compressible extension
would be needed before interpreting this incompressible preparation as an
acoustic horizon. These are specified connections to build, rather than a
reason to exclude the useful fluid blocks from the atlas.

The next task is a short retarded evolution of this fixed initial preparation,
checking T7 and its energy transfer while retaining the exact source law.
Coupling the resulting stress to a selected gravitational branch is a later
calculation; no dark-component abundance or local vacuum term is assigned
to the covariance by its name alone.

```powershell
python scripts/run_python.py scripts/check_nsc_atlas_vortex_transfer.py
python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_vortex_transfer.py tests/test_nsc_atlas_vortex_composition.py
```
