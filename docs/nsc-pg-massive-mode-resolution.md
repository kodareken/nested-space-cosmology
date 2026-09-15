# Massive modes on the common PG Cauchy surface

The matched horizon preparation now has a spatial mode owner on both sides
of the fixed horizon. It uses the same charged Dirac operator, angular bundles,
affine horizon covariance and inherited incoming occupation. The source is
the horizon/infinity data, with the historical and matched seed generations
retained as numerical controls.

The new owners are
[`nsc_pg_massive_modes.py`](../src/recursive_horizons/nsc_pg_massive_modes.py)
and [`nsc_massive_jost_modes.py`](../src/recursive_horizons/nsc_massive_jost_modes.py).
The [runner](../scripts/derive_nsc_pg_massive_modes.py) evaluates all existing
positive-energy nodes and required angular signs. Its finite spatial samples
are checks of a continuous-energy construction; they are not a finite set
whose Gram matrix may replace the continuum CAR identity.

## 1. Operator, coordinates and source channels

Use the already owned common-time half-density operator

$$
H_D=-i\left(v\partial_\rho+\frac{v'}2\right)+M,
\qquad v=\sigma_2-\beta I,
\qquad M=-m\sigma_1-\frac{\lambda}{r}\sigma_3.
$$

The measure is $d\rho$, after the angular and compact normalizations already
recorded in the paired-basis calculation. The transparent neck remains at
$\rho=0$. The horizon is at the existing $\rho_h$; it is not a reflecting
endpoint. The whole-line self-adjoint closure is imported from the
[Dirac-domain result](nsc-dirac-tetrad.md), with the bounded compact mass
term included. The two bulk support projectors and the T seam are unchanged.

For a signed angular channel, the incoming spectral coordinates are the
exterior outgoing-horizon mode, its interior partner, and the incoming mode
from parent infinity. Their physical rank is

$$
P_{\mathrm{src}}(E)=\mathrm{diag}
\left(1,1,\mathbf{1}_{|E|>m}\right).
$$

The threshold is the parent dispersion threshold $|E|=m$. It is distinct
from the inherited occupation threshold $|E|=\Omega m$. An open but empty
incoming channel still belongs to the CAR space; small transmission does
not close it. Below the mass gap, the third coordinate is a zero column
used for storage, and the physical fiber has rank two.

For $E>0$, the two horizon sewings are

$$
M_{\mathrm{in}}=
\begin{pmatrix}0&1&0\\R_\lambda&0&\sqrt{T_\lambda}\end{pmatrix}
P_{\mathrm{src}},
\qquad
M_{\mathrm{ext}}=
\begin{pmatrix}1&0&0\\R_\lambda&0&\sqrt{T_\lambda}\end{pmatrix}
P_{\mathrm{src}}.
$$

The interior partner is independent of the exterior outgoing-horizon
coordinate. Only the ingoing characteristic is continued across the horizon.
This differs from imposing full-spinor continuity at the horizon or using
two globally continuous columns. At the regular neck, ordinary T transparency
continues to apply.

The source covariance is the existing $C_H(E)\oplus n_{\mathrm{in}}(E)$,
restricted by $P_{\mathrm{src}}$. Its off-diagonal horizon-partner correlation
is kept. Between opposite angular signs, $Z=0$ follows from the already
verified angular identity and sign-preserving radial evolution. Negative
energies use the locked paired $\sigma_3 K$ map and the signed source law.

## 2. Massive propagation and PG frame

The inner frame starts with the versioned `matched_delta_q` argument

$$
z_H=\sqrt{\lambda^2+(m r_h)^2}\sqrt{2\delta_q/\kappa_h}.
$$

It evolves under the existing massive radial generator. No time duration is
chosen: the integration coordinate labels points of the same stationary
geometry. The outer frame uses a radial collar, for which the existing
inverse-$r_h$ conversion remains correct.

The exact clock change is $t=\tau-F(\rho)$ with $F'=\beta/A$ and $F(0)=0$.
Only the geometry primitives of the earlier flow atlas are reused. The LLL
mode functions and LLL covariance are not assigned to massive channels.
Writing $R_0=(I+i\sigma_1)/\sqrt2$, the exterior half-density conversion is

$$
S_e=A^{-1/4}e^{\mathrm{atanh}(\beta)\sigma_2/2}
R_0\,\mathrm{diag}(e^{-i\pi/4},e^{i\pi/4}),
\qquad S_e^\dagger vS_e=\sigma_3.
$$

For the interior, with $a=\sqrt{-A}$,

$$
S_i=a^{-1/2}e^{\mathrm{atanh}(1/\beta)\sigma_2/2}R_0,
\qquad S_i^\dagger vS_i=-I.
$$

Both conversions include their radial connection and $iEF'$ in the
transformed Dirac equation. These terms are checked against the T-basis
operator, not dropped as a phase convention.

## 3. Stable exterior amplitudes and their numerical error

A reflection coefficient alone is insufficient for stable evaluation of an
evanescent exterior mode. Propagating a nearly canceled pair of large basis
columns outward can spoil its amplitude while leaving the net current tiny.
The new owner carries the complex outgoing/decaying Jost amplitude inward
and fixes its outgoing-horizon coefficient to one. Above the mass gap it also
propagates the horizon-ingoing solution. Below the gap it evolves the real
Riccati phase, preserving the decaying mode's zero current.

The ratio equation is the same one used by the old scattering owner. Its
outer data now retain the NSC massive asymptotic spinor series through order
eight. The zero-mass limit agrees with the existing massless ratio series.
The numerical outer radius increases only when the series residual requires
it. A doubled outer radius checks the resulting probe-region mode error.
This is numerical boundary control, not a new physical scale or a fit.

If $t_\mathrm{out}$ is the complex transmitted amplitude in the declared
outer Jost coordinate, the third future row is

$$
\left(t_\mathrm{out},0,
-e^{i\arg t_\mathrm{out}}R_\lambda^*\right).
$$

Together with $M_{\mathrm{in}}$ it satisfies $S^\dagger S=P_{\mathrm{src}}$.
The incoming horizon coefficient $\sqrt T$ is an allowed phase convention
because the incoming occupation is diagonal. The complex outer amplitude is
still calculated and recorded. Its phase origin is the specified outer Jost
coordinate; it does not select a Cauchy history.

The original scattering cache and both seed generations stay immutable.
Changes from their finite outer approximation and 6,000-step seed integrator
are reported as control differences. Their old residuals are not relabeled
as the accuracy of this spatial calculation.

## 4. Whole-line Green function and spectral measure

Import the first-order Green jump and Stone formula. The new calculation
binds them to this domain and these source columns. In the trapped region,
let $\Psi_i^\dagger(-v)\Psi_i=I$. For two trapped points,

$$
R^+(x,y)=i\,\mathbf{1}_{x<y}\Psi_i(x)\Psi_i(y)^\dagger.
$$

For an exterior source, write the outgoing and ingoing exterior fields as
$u_o,u_h$. The coefficient rows are fixed by

$$
\begin{pmatrix}u_o(y)&-u_h(y)\end{pmatrix}
\begin{pmatrix}\alpha(y)\\\beta_h(y)\end{pmatrix}
=i\,v(y)^{-1}.
$$

The field is $u_o(x)\alpha(y)$ above the source and
$u_h(x)\beta_h(y)$ below it; the latter continues into the interior ingoing
column. For an interior source and exterior observation, $R^+=0$.
The lower boundary value is the adjoint with exchanged points. These four
cases include the independent interior partner and the causal asymmetry.
They are not an exterior-only Jost construction.

With the complete source columns $\Phi_E$, the normalization check is

$$
\frac{R^+(x,y)-R^+(y,x)^\dagger}{i}
=\Phi_E(x)\Phi_E(y)^\dagger,
\qquad
\frac{dP}{dE}(x,y)=\frac{\Phi_E(x)\Phi_E(y)^\dagger}{2\pi}.
$$

The runner evaluates interior/interior, exterior/exterior and both mixed
regions, above and below the mass gap. A rejection check removes an open
infinity column and verifies that the spatial spectral identity fails.
At coincident points the symmetric Green value is used; a one-sided value
must not double the spectral density.

## 5. Conditions for the continuum completion

The global identification also needs spectral completeness; a finite grid
or an isometric sewing matrix does not establish it. The relevant application
of the imported spectral theorem has three parts:

1. **No point spectrum.** For a real-energy eigenfunction, the conserved
   current is constant. In the interior, $v$ is negative definite and
   $\beta\sim\sqrt{3\pi}|\rho|$ toward child infinity. Nonzero current gives
   $|\chi|^2\geq |j|/(\beta+1)$, whose integral diverges. Thus an $L^2$
   eigenfunction vanishes throughout the interior. On the exterior horizon
   branch the outgoing amplitude behaves as $|\rho-\rho_h|^{-1/2}$; its
   squared norm is not integrable, including at $E=0$. The remaining ingoing
   branch must match the zero interior. This also excludes threshold atoms.
2. **Continuous boundary values away from thresholds.** The owned exterior
   metric has $A=1-2/\rho+O(\rho^{-3})$. With
   $p=\sqrt{E^2-m^2}$, the outgoing phase derivative is
   $p+(2E^2-m^2)/(p\rho)+O(\rho^{-2})$. The exact PG clock factor and this
   Coulomb phase must accompany the radial spinor expansion; the angular
   $\lambda/\rho$ term is retained in that expansion. After this change of
   basis, its derivative is retained and the remaining Jost coefficients
   are integrable at infinity. The construction is locally uniform in
   energy away from $E=\pm m$. At the horizon the Frobenius coupling decays
   in the tortoise coordinate. Opposite current signatures keep the
   ingoing/outgoing Wronskian nonzero on the real axis away from thresholds.
3. **Use the whole-line jump.** Applying Stone to the four-region Green
   function above gives the spectral measure on each such energy interval.
   The exceptional threshold set is finite and has no atoms by part 1;
   it cannot support a non-atomic singular measure. This is the route to
   the continuum CAR identity, rather than a finite-frequency Gram fill.

The executable record must distinguish these analytic conditions from the
measured mode and outer-boundary residuals. Subsequent spatial covariance
assembly must integrate the complete spectral law with its contact/tail
contribution; the finite arrays saved here do not by themselves perform that
integration or determine a stress tensor.

## 6. Executed result

All **32 massive channel groups** are represented: 3,584 positive-energy
and angular-sign labels, with negative energies fixed by the signed map.
The [record](../results/development/nsc-pg-massive-mode-resolution.json)
and its immutable mode-field artifact contain the per-group checks.

| Check | Maximum residual | Tolerance |
|---|---:|---:|
| PG Dirac-frame conversion | $2.24\times10^{-13}$ | $3\times10^{-11}$ |
| Complex scattering isometry | $1.12\times10^{-12}$ | $3\times10^{-11}$ |
| Physical source covariance law | $0$ | $3\times10^{-11}$ |
| Source CAR violation | $2.66\times10^{-17}$ | $3\times10^{-11}$ |
| Interior current | $5.34\times10^{-12}$ | $3\times10^{-9}$ |
| Exterior current | $1.81\times10^{-10}$ | $3\times10^{-9}$ |
| Whole-line spectral jump, relative | $8.29\times10^{-10}$ | $3\times10^{-9}$ |
| PG fields under outer-radius doubling | $7.79\times10^{-10}$ | $3\times10^{-9}$ |

The Green residual is scaled by $\max(1,\|\Phi(x)\Phi(y)^\dagger\|_F)$;
its absolute maximum is $4.59\times10^{-9}$. The matched seed-control
difference reaches $3.88\times10^{-10}$ with the new spatial integration
and outer Jost approximation. That difference is recorded explicitly;
the old 6,000-step seed generation is not replaced.

**Mode/Green construction: PASS. Spatial covariance integration: OPEN.**
The next calculation is the packet energy integral with its full contact and
tail account. The PDF remains v0.25.0. No metric history, stress or transmitting
CTP endpoint jet is assigned by this record.

```sh
# New global mode-field evaluation; old seed/stress generators are not called.
python3 scripts/derive_nsc_pg_massive_modes.py --workers 4

# Focused verification reuses the authenticated field artifact.
python3 scripts/check_nsc_pg_massive_modes.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_pg_massive_modes.py
```
