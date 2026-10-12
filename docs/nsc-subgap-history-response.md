# Compact subgap contribution to the nonlinear CTP variation

This owner connects the existing group-13 subgap state to the spatial
metric jets of the same nonlinear supplied history. It computes a finite
spectral-panel contribution to the Gaussian action variation. It does not
assign an absolute stress or select a metric history.

The inherited panel is $1<E<m$, $m=\pi/2$, with zero angular label. The
spatial response directions have support inside the trapped region, with
a causal buffer to the computational boundary. On this initial causal patch
the exterior horizon column $u$ vanishes. The remaining analytic columns
are ordered $(v,w)$ and use the existing matched-$\delta_q$ preparation.
The physical subgap modes are $(Rv,w,0)$, where $R(E)$ is the already owned
whole-line reflection coefficient. No stored seed covariance is an input.

For each of the four raw KS amplitude directions, the existing history owner
evolves the spatial fields and their variations together. Write

$$
B_{ab}(z)=i\int d\rho\,[U_g\phi_a(\bar z)]^\dagger
                         \delta U_g\phi_b(z).
$$

The conjugate-frequency solution supplies the bra. This is an analytic
bilinear; a complex frequency is never passed to a Gaussian state constructor.
For real energy, $B=iU_g^\dagger\delta U_g$ is Hermitian. With the inherited
horizon functions $d=f-1/2$ and $s$, the centered positive-energy trace is

$$
d(B_{vv}-B_{ww})+2\operatorname{Re}[-isR B_{wv}].
$$

The known signed-frequency map in this zero-angular sector gives the two
signed panels' Gaussian derivative

$$
\delta\Gamma_{13,\mathrm{subgap}}=
-\frac1\pi\left[\int_1^m d(B_{vv}-B_{ww})\,dE
+2\operatorname{Re}\int_{\cal C}(-isR B_{wv})\,dz\right].
$$

The minus sign is fixed by the existing $-i\operatorname{Tr}(C U^\dagger
\delta U)$ convention. The half-identity subtraction here is exactly the
algebraic signed-pair combination of that bare action, not a new vacuum
subtraction. Angular degeneracy, copy counts and compact weights are not
inserted into this per-reduced-family result.

Only the reflection term uses the existing pole-free upper contour. The
smooth unsewn response is interpolated on nested real Chebyshev nodes and
checked against an independently propagated conjugate-frequency pair. The
winding reflection itself is never interpolated. Existing numerical
reflection samples may be reused only from the exact source/configuration
cache key; the new immutable payload stores every sample it actually uses.

The record separates spatial, time, interpolation and contour errors. A
real-frequency control compares the new unsewn basis and contraction to the
authenticated earlier history-jet payload. No old generator is rerun.

The history and four perturbation amplitudes remain the published response
controls. This calculation supplies neither physical time endpoints nor
the full spectral integral. The other energy windows/families and the
reference/local allocation on the spatially varying history must be joined
before an absolute source or extended stationarity claim. No independent
$\Gamma_{\rm rest}$ term is introduced. Metric evolution remains closed.

```sh
python3 scripts/derive_nsc_subgap_history_response.py --prepare
python3 scripts/derive_nsc_subgap_history_response.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_subgap_history_response.py
```
