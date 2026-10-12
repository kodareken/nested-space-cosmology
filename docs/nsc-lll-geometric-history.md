# The existing LLL geometric allocation on the common KS history

The massless magnetic sector already carries a state-independent conformal
allocation, separately from its occupation moments and the compact induced
coefficients. The new owner evaluates that same allocation on the supplied
spatially varying history. Its coefficient is the retained multiplicity
$c=|q_{\rm mag}|=4$, counted once. This fills a missing contribution to the
common source; it introduces no additional field or independently weighted
functional.

## Allocation and canonical reference

The [serialized-state owner](nsc-mode-resolved-cauchy-state.md) separates the
LLL occupation moments from `conformal_stress(A, A_prime, A_second, 0, 0, c)`.
The [horizon source](nsc-horizon-source.md) defines its physical radial metric,
anomaly normalization and sphere projection. Its imported conformal stress
and Schwarzian formulas are reused here. The
[reference-band owner](nsc-reference-band-action.md) has a constant chiral
LLL projector on each nonzero-momentum sign chart. We use the corresponding
canonical spatial normal ordering, with no change of initial state.

Write the actual KS metric as

$$
ds_2^2=N^2dT^2-a^2(dz+\beta\,dT)^2,\qquad
L=\log a,\qquad H_a=\frac{a_T-(a\beta)_z}{Na}.
$$

Here $H_a$ is the local axial expansion of a supplied history, not a solved
cosmological expansion. For its null maps $u,v$, let
$\sigma=\log a-\tfrac12\log|u_zv_z|$. Subtraction of the canonical chiral
moments from the covariant conformal tensor eliminates the null maps and
both state functions. In detail, use

$$
j_+=u_z^2t_u-\frac{c}{24\pi}\{u,z\},\qquad
j_-=-v_z^2t_v+\frac{c}{24\pi}\{v,z\},
$$

$$
\rho_{\rm can}=p_{\parallel,\rm can}=\frac{j_+-j_-}{a^2},\qquad
T_{\hat0\hat1,\rm can}=-\frac{j_++j_-}{a^2}.
$$

The minus sign in the covariant orthonormal momentum component agrees with
the existing child-frame source. The required difference is

$$
\Delta\rho=\frac{c}{24\pi}
\left[\frac{2L_{zz}-L_z^2}{a^2}-H_a^2\right],\qquad
\Delta T_{\hat0\hat1}=\frac{c}{12\pi a}\partial_zH_a,
$$

$$
\Delta p_\parallel=\Delta\rho-\frac{c}{24\pi}R_2,
\qquad \Delta p_\perp=0.
$$

For an exact cancellation check set $P=u_z\sigma_u$, $Q=v_z\sigma_v$,
$U=\log|u_z|$, $V=\log|v_z|$. Null transport gives
$P=(L_z-aH_a-U_z)/2$, $Q=(L_z+aH_a-V_z)/2$. The reproducer substitutes
these identities into the inherited stress formulas; every state/ray-map
term cancels. No ray map, transport duration or physical history is selected.

The difference alone is reference dependent, not a separately covariant
tensor. Covariance and conservation belong to its sum with the consistently
normal-ordered state source. The physical radial LLL action is independent
of the sphere radius, so its angular pressure vanishes in this sector and
the radial four-dimensional allocation is $\Delta t_{ab}/(4\pi r^2)$.
There is no change to a radius-rescaled conformal frame.

## Action and all four metric directions

Up to boundary terms fixed by the existing compact, smooth variation domain,
the same allocation has the local action representative

$$
S_{\rm LLL,geo}=\frac{c}{24\pi}\int dT\,dz
\left[-NaH_a^2+\frac{2N_zL_z}{a}-\frac{NL_z^2}{a}\right].
$$

Its metric gradient densities are

$$
\left(\frac{\delta S}{\delta N},\frac{\delta S}{\delta\beta},
\frac{\delta S}{\delta a},\frac{\delta S}{\delta r}\right)
=(-a\Delta\rho,-a^2\Delta T_{\hat0\hat1},N\Delta p_\parallel,0).
$$

The CTP contribution is $S[g_+]-S[g_-]$ with these relative-branch
derivatives. The existing source convention is force = minus action gradient.
We evaluate both the action derivative and the independent stress pairing.
Their equality is tested on the same supplied compact KS/PG history as the
[nonlinear Gaussian calculation](nsc-nonlinear-ks-source.md).
The fixed Jacobian is $dT\,dz=d\tau\,d\rho/a_0$, with
$\partial_T=(\beta_0/a_0)\partial_\tau-a_0\partial_\rho$ and
$\partial_z=\partial_\tau$. Exact metric jets and $R_2$ are reused from the
[spherical local-history owner](nsc-spherical-local-history.md).

On the old homogeneous KS geometry this reduces to

$$
\Delta\rho=-\frac{cA'^2}{96\pi(-A)},\qquad
\Delta p_\parallel=\Delta\rho-\frac{cA''}{24\pi},\qquad
\Delta T_{\hat0\hat1}=0,
$$

exactly its previously counted geometric source. It replaces that source's
stationary evaluation; adding the old `conformal_stress(..., 0, 0)` tensor
afterward would count the allocation twice. It contains neither the thermal
state constants nor an extra compact multiplicity.

## Computation and remaining composition

The [record](../results/development/nsc-lll-geometric-history.json) compares
the action and stress derivatives in all four raw-KS directions, records
separate time/radial quadrature refinements, checks the stationary source
and the trace, and authenticates all inputs. Its reproducer is:

```sh
python3 scripts/derive_nsc_lll_geometric_history.py --check
```

This is the geometric LLL allocation, not the full renormalized stress.
The physical state, fourth-order subtraction and band-connection terms
still need their common spectral trace and convergence calculation.
No exact equality to the entire raw heat-vacuum functional is imposed as
an additional gate. The geometric cancellation at the fixed smooth seam is
preserved; $93.54264532195464$ remains its historical time-node diagnostic.
No extra $\Gamma_{\rm rest}$, stationary solution, updated null/constraint
claim, metric step or PDF revision is produced.
