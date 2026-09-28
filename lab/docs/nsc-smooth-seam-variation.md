# Geometric matching on the actual smooth transmitting surface

The local geometric boundary variation is now evaluated on the shared
$\rho=0$ surface of the already declared global metric. It cancels between
the two regions under the metric-jet matching built into the smooth domain.
The old $93.54264532195464$ history-node diagnostic is preserved and is not
used as either side of this interface calculation.

## One surface and its two outward normals

Z1 defines parent and child by the two spatial supports of one field on one
smooth PG geometry. It permits smooth compactly supported coefficient
variations of that geometry. The boundary data therefore come from a common
metric jet, with opposite outward normals; they are not two independently
chosen endpoint geometries.

In the existing future-child convention, $a=\sqrt{\beta^2-1}$ and the unit
future normal in PG coordinates is

$$
n=\frac{\beta}{a}\partial_\tau-a\partial_\rho.
$$

On the stationary background scalars $f(\rho)$ its action reduces to
$D_Tf=-a\,df/d\rho$. That abbreviated derivative must not be used for
time-dependent spinor fields. The parent region is on the earlier side of
the $\rho=0$ surface, so $n_p=n$ and $n_c=-n$. The intrinsic metric
variations agree, while the outward-normal shear variations oppose:

$$
\delta a_c=\delta a_p,\qquad
\delta r_c=\delta r_p,\qquad
\delta\sigma_c=-\delta\sigma_p,
\qquad \sigma=H_a-H_r.
$$

For the boundary-jet covector derived from the existing Weyl density,
$\overline P_a$ and $\overline P_r$ are odd under normal reversal;
$\Pi_\sigma$ is even. The common-interface variation consequently contains

$$
(\overline P_{a,p}+\overline P_{a,c})\delta a
+(\overline P_{r,p}+\overline P_{r,c})\delta r
+(\Pi_{\sigma,p}-\Pi_{\sigma,c})\delta\sigma.
$$

All three coefficients vanish for the same smooth jet and equal locked
coefficients. This imports the boundary-orientation rule and applies it to
the current NSC domain. It does not derive a new gravity theory or select a
geometry that satisfies the bulk source equations.
Both sides have the same jets in one common future-oriented coordinate;
the odd signs arise from opposite outward normals, not from an assumed
time-reflection symmetry of the background.

## Actual reference values

The existing analytic metric at $\rho=0$ gives

$$
a_0=\sqrt{\frac{3\pi}{2}-1},\quad r_0=1,\quad
D_Ta=3,\quad D_Tr=0,\quad D_T^2r=a_0^2.
$$

Its Weyl combination is $F=2$ and $D_TF=-6a_0$. With the unchanged
$C_W=-0.0008441130342798952$, the parent outward-normal coefficients are

| Coefficient | Value |
|---|---:|
| $\overline P_a$ | $0.4150924626502038$ |
| $\overline P_r$ | $-0.6300629292752326$ |
| $\Pi_\sigma$ | $0.10900236876297109$ |
| $\sigma$ | $1.557021169290531$ |

The child has opposite metric covector entries and opposite shear, with the
same shear-conjugate coefficient. Einstein–GHY momenta cancel under the
same orientation. The already included Euler and $\Box R$ endpoint
primitives are odd and also cancel on matched smooth jets. No geometric
boundary term is introduced or counted a second time.

The [record](../results/development/nsc-smooth-seam-variation.json) contains
the exact reference jets, each side separately and the common-covector
residuals. The independent highest-jet check uses the previous Weyl density;
the declared tolerance is $3\times10^{-11}$.

## The remaining physical calculation

**Geometric seam matching PASS; full stationarity OPEN.** This result closes
the local geometric part of the comparison on the actual smooth surface.
It applies to the shared smooth-jet domain, not to an independently varied
thin shell or arbitrary initial/final time boundaries.

It does not set $\Gamma_{\rm rest}$ to zero. The full Gaussian determinant
already contains its Schur factor and the complementary bulk; the latter
cannot be added again as a new interface force. Quantum/field variations
not yet identified in the complete same-action functional remain unevaluated.
In particular, geometric boundary cancellation does not make the bulk
Einstein–quantum-stress residual vanish.

The original time-node coefficients, physical scales, prepared state and
earlier certificates are unchanged. No new stress or metric timestep is
assigned. Work now belongs to the complete causal source functional and its
physical variational domain, rather than forcing a counterterm to cancel
the old time-node number.

```sh
python3 scripts/derive_nsc_smooth_seam_variation.py --check
python3 -m pytest -q tests/test_nsc_smooth_seam_variation.py
```
