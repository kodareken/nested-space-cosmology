# The recorded Dirac source in the PG constraint frame

The already computed canonical parent-state tensor is now expressed in
the frame of the [ADM source equations](nsc-adm-source-constraints.md).
This is a transformation of authenticated results, not another angular
sum or vacuum calculation. Its Ward energy flux reproduces the previously
recorded outward parent power.

The [input tensor](../results/development/unruh-state.json) uses the future
child time T and its orthonormal spatial direction z. Inside the horizon
write f=-A=beta²-1>0. At the neck, r=N=q=1 and beta²=3 pi/2. The coframes are

\[
\theta^0_{\rm PG}=d\tau,\quad\theta^1_{\rm PG}=d\rho+\beta d\tau,
\qquad
\theta^0_c=dT=-d\rho/\sqrt f,\quad
\theta^1_c=\sqrt f\,dz=\sqrt f\,d\tau+\beta d\rho/\sqrt f.
\]

Thus

\[
\theta_c=L\theta_{\rm PG},\qquad
L=\frac1{\sqrt f}\begin{pmatrix}\beta&-1\\-1&\beta\end{pmatrix},
\qquad T_{\rm PG}=L^T T_c L.
\]

Using the full stored tensor, including its mixed component, gives:

| Canonical massless contribution | Child frame | PG normal frame |
|---|---:|---:|
| rho | -0.00330689624 | -0.01752638060 |
| T01 | -0.00000304859 | 0.03087075875 |
| radial pressure | -0.04949459659 | -0.06371408094 |
| sphere pressure | -0.04299577129 | -0.04299577129 |

The units and canonical renormalization mu=1 are unchanged. At this neck
the four source-force densities are

\[
(F_N,F_\beta,F_q,F_r)_{\rm canonical}
=(-0.22024299417,\ 0.38793339565,\ 0.80065475448,\ 1.08060159380).
\]

They enter different equations and must not be replaced by a single null
contraction. The ADM Ward flux is

\[
\mathcal J=-\beta F_N-(\beta^2+1)F_\beta+\beta F_q
=-4\pi f\,T_{\hat T\hat z}
=0.0001422206795425.
\]

It agrees with the stored parent Killing power to about 4e-17. This fixes
the orientation and normalization of the connection to the constraint
equations. It does not identify this power with a cosmological density
injection Q.

The required two-derivative geometric tensor must be transformed too.
Keeping a_EH symbolic, its PG components are

\[
(\rho,T_{01},p_r,p_\perp)_{\rm required}
=a_{\rm EH}(-2,\ 4\beta,\ -2(1+3\pi),\ 2(1-3\pi)).
\]

The source-side force target is consequently

\[
F_{\rm required}
=a_{\rm EH}(-8\pi,\ 16\pi\beta,\ 8\pi(1+3\pi),\ 16\pi(3\pi-1)).
\]

The Einstein action contributes minus this target to the total Euler
force. The remaining source equation is therefore
F_canonical+F_rest=F_required(a_EH), with F_rest derived from the same
action. This equation defines the residual; it does not assign a fitted
value to F_rest or a_EH. The earlier child-frame comparison and this
full-tensor comparison are equivalent.

The [new record](../results/development/adm-neck-source-map.json) preserves
all transformed components, the invariant trace, the power and the
symbolic geometric target. Run

```sh
python -B scripts/map_nsc_adm_neck_source.py --check
```

Only the frame/source map is evaluated. The original canonical tensor,
its uncertainties and its remaining physical assumptions are reused.
The full vacuum-conversion and other-sector forces, and the coupled
geometry solution, are still required for self-sourcing.
