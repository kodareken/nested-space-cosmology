# Directed covariance transport for an original subgap source mode

This owner bounds the unperturbed source-preparation covariance at rho=1 for one original
group14/low16_1 energy and angular sign. It keeps the original columns and
source covariance bytes. It does not yet cover the energy panels, the opposite
angular partner, numerical continuation to rho_up, or the complete preparation
component of the incoming gate. The changed geometry's incoming state is
evaluated later by the existing retarded evolution. The gate remains OPEN.

Owner: `src/recursive_horizons/nsc_subgap_source_covariance.py`.
It consumes the [enclosed horizon reflection](nsc-metric-horizon-frame.md).

## Matching the original phase conventions

Let q=pi/2+atan(rho), u=q-q_h and r_h=csc(q_h). The direct frame has leading
columns u^(i s E/(2 kappa_geometry)) e_s. The original radial convention has
the additional diagonal phase D(c), with

$$c=\frac{E\log r_h^2}{2\kappa_{\rm geometry}}
-\frac12\operatorname{atan2}(m r_h,\ell).$$

The interior convention also has D_pi=diag(exp(-i pi/4),exp(i pi/4)). These
follow from the original rotations and rho-h=r_h^2 u+O(u^2). Uniqueness of
the normalized Frobenius frame propagates the constant changes of basis:
F_ext,legacy=F_ext,direct D(c) and
F_in,legacy=F_in,direct D(c) D_pi. The reflection changes as
R_legacy=exp(2ic) R_direct.

The existing sewing identity

$$D(c)S(e^{2ic}R,T)=S(R,T)\operatorname{diag}(e^{ic},e^{ic},e^{-ic})$$

leaves the full original three-port covariance invariant: its horizon pair
gets one common phase, and the infinity port has no cross-correlation with
that pair. The known D_pi must still be retained. No relative coherence is
discarded. In the subgap T=0 and the true reflection has unit modulus by
zero exterior current and the frame's conserved signed current.

Using f=1/(1+exp(2 pi E/kappa_source)) and
c_h=1/(2 cosh(pi E/kappa_source)), the covariance in the direct frame is

$$C_h=\begin{pmatrix}1-f&c_h\overline R\\c_hR&f\end{pmatrix}.$$

The stored source kappa remains unchanged. The geometric slope belongs to
the mode equation and is not substituted into the source occupation.
The reflection enclosure is represented by its enclosed argument; the true
unit-circle point is preserved. No stored field or covariance is normalized.

The source has trace one because T=0 and |R|=1. Its coherence obeys
c_h^2=f(1-f). The exact interior frame is unitary: its generator is
anti-Hermitian and its normalized horizon limit is unitary.

## Complete covariance at each energy

Write Q=(I+n dot sigma)/2. This is the complete 2x2 covariance at this energy,
not a closure of the energy continuum or the evolved inhomogeneous field.
For delta=q_h-q>0 and y=log(delta), define H(delta)=-W(q_h-delta)/(-delta),
where W is the original exact compact metric. The spin Hamiltonian is

$$\mathbf h=
\left(-m\,\csc(q_h-\delta)\sqrt{\delta/H},
\ell\sqrt{\delta/H},-E/H\right).$$

The matrix equation Q_y=-i[h dot sigma,Q] becomes

$$\mathbf n_y=2\mathbf h\times\mathbf n.$$

Its real generator is skew-symmetric, so its propagator preserves the
Euclidean norm of an error. Common wave phases do not enter Q. The original
PG-to-KS restriction differs from the working frame by a common phase, which
also cancels in its covariance. Both claims are specific to this homogeneous
source preparation; the later state still needs its full retarded evolution.

This last cancellation follows from the actual two production maps. With
$B_\xi=\cosh(\xi/2)I+\sinh(\xi/2)\sigma_2$, the interior
`MassivePGModeResolution.frame` is $f=B_\xi U/\sqrt a$, whereas
`TransmittingDiracSeamDomain.trace_map` at unit quadrature weight is
$T=B_\xi U/(r\sqrt a)$. Thus $f=rT$. Writing the sewing matrix as
$\mathsf S$, the resolved PG columns are
$e^{iE\,\mathrm{clock}}fV\mathsf S$, and `restrict_resolved_modes` applies
$e^{iES}T^{-1}/r$. Their canonical columns are consequently the working
columns multiplied by the scalar phase $e^{iE(\mathrm{clock}+S)}$.
That phase cancels in the covariance at each energy. Here $V$ denotes the
working-frame fundamental solution and $S$ the chart shift. This comparison changes
neither the original columns nor their three-port covariance.

## Exact metric jets and continuous defect

Let C=1.5 sin(2q_h)+0.5 cos(2q_h) and
S=1.5 cos(2q_h)-0.5 sin(2q_h). The stable analytic metric has
H(0)=3-2S and, for n>=1,

$$H_n=-\frac{2^{n+1}}{(n+1)!}
\left[C\cos\frac{(n+1)\pi}{2}+S\sin\frac{(n+1)\pi}{2}\right].$$

Because |C|,|S|<=2, |H_n|<=4*2^n/(n+1)!. The omitted series and all required
derivatives are bounded by positive exponential tails on |u|<0.4. Their
derivative bounds are composed into u(y), including its higher coefficients.
They are not copied from one coordinate to another. Positivity of H is
checked on each validation cell.

A DOP853 curve supplies a numerical witness. Each cell is reconstructed from
its actual binary endpoint values and six dense correction rows, using the
exact difference of its stored times. For x in [0,1] and h=y1-y0, the normalized
defect is

$$d(x)=p_x-2h\,\mathbf h(y_0+hx)\times p(x).$$

Arb Taylor coefficients at the midpoint, together with a coefficient bound
over the full cell, enclose its supremum norm. Thus

$$\|n_{\rm true}-p\|_{\rm end}
\le\epsilon_{\rm initial}+\sum_{\rm cells}\sup_{[0,1]}\|d\|
+\epsilon_{\rm endpoint}.$$

There is no extra cell-width factor: h is already in the normalized defect.
The endpoint bridge encloses the difference between the numerical log endpoint
and the exact target rho. The original source's initial uncertainty comes
from the horizon-frame and reflection enclosures and source arithmetic.

## Comparison with the unchanged original data

The stored Q_num=A_num C_src,num A_num^dagger is evaluated in ball arithmetic
from the original binary 2x3 columns and Hermitian 3x3 source covariance.
Let t_num=Tr Q_num and n_num denote its Bloch vector. The operator error is
bounded by

$$\|Q_{\rm true}-Q_{\rm num}\|_{\rm op}
\le\tfrac12\left(|1-t_{\rm num}|+
\|n_{\rm end}-n_{\rm num}\|+\epsilon_{\rm Bloch}\right).$$

This comparison includes numerical source occupation/coherence arithmetic
and original modal preparation error at the selected energy. It does not
pretend that a reconstructed packet's diagnostic accuracy is a certificate.
The saved source data remain the reference inputs for later field evaluation.

The record binds the original inventory, horizon result, source code and small
trajectory payload. Missing cells, non-Hermitian data, a non-subgap channel,
or an unresolved metric are rejected. A regression detects multiplying the
normalized defect by the width twice. Dropping source coherence produces a
large covariance discrepancy and cannot pass as the original source.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_subgap_source_covariance.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_subgap_source_covariance.py --check
```
