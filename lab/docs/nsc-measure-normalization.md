# Common-action normalization: exact profiles versus the physical measure

The previous calculation supplies a covariant four-dimensional proper-time
modulus and an explicit ultraviolet subtraction on the smooth static cell.
Its finite-matrix predecessor still contains an additive rank term

\[
N\bigl[\log(M/\Lambda)+\gamma_E/2\bigr].
\]

That term cannot be promoted to an unrestricted frequency/angular integral of
\(\operatorname{Tr}1\). This note compares the exact finite spectral profiles
that people actually write for that missing normalization, computes their
metric response with the current covariant frequency operator, and keeps the
unknown finite/compensator data explicit. No coefficient is chosen to hold a
throat open. The identity is not inserted as a completed \(\Gamma_{\rm one}\).

The record is [`results/nsc-10-measure-normalization.json`](../results/nsc-10-measure-normalization.json).
The owner is `src/recursive_horizons/nsc_measure_normalization.py`. It reads
the existing operator; it does not replace it.

## Prescriptions, not one measure

Write \(x=d^2/\Lambda^2\) for a Dirac eigenvalue \(d\). The three smooth
profiles compared here are

\[
g_{\rm PT}(x)=\tfrac12 E_1(x),
\]

\[
g_{\rm log}(x)=-\tfrac12 e^{-x}\log(d^2/M^2),
\]

\[
h(x)=-\tfrac12\bigl[e^{-x}(\log x+\gamma_E)+E_1(x)\bigr].
\]

\(g_{\rm PT}\) is the calculated determinant modulus of
[the covariant source](nsc-covariant-source.md), not a completed
\(\Gamma_{\rm one}\): one eigenvalue contributes \(\tfrac12 E_1(D^2/\Lambda^2)\),
and the four-dimensional counting is the existing frequency integral. \(M\)
is a normalization mass, not a fermion mass.

\(g_{\rm log}\) is a **smooth weighted-log** profile. It is not the literal
hard projector used as the primary regularization in
[Andrianov, Kurkov and Lizzi, arXiv:1106.3263](https://arxiv.org/html/1106.3263v1).
That paper cuts the spectrum with

\[
P_N=\Theta\bigl(1-D^2/\Lambda^2\bigr),\qquad
\log Z=\sum_{|\lambda_n|\le\Lambda}\log(\lambda_n/\mu),
\]

and states that the split into Weyl-invariant and anomalous factors is not
unique. A smooth \(\chi(D/\Lambda)\) is mentioned there only as a possible
cutoff function; the exponential weight on \(\log(d^2/M^2)\) is a different
object. The comparison below therefore keeps a hard-log diagnostic

\[
g_{\rm hard}(x)=-\tfrac12\,\Theta(1-x)\log(d^2/M^2)
\]

separate from \(g_{\rm log}\).

The old finite-rank convention is equivalent, on a matrix of size \(N\), to
adding \(N[\log(M/\Lambda)+\gamma_E/2]\) to \(\sum g_{\rm PT}\). Replacing
\(N\) by the heat trace \(\operatorname{Tr}e^{-D^2/\Lambda^2}\) is a
covariantization of that term, not a theorem that the continuum measure
equals \(g_{\rm log}\).

## Exact identity

With \(\log(M/\Lambda)\) written as a single ratio,

\[
\log(d^2/M^2)=\log x-2\log(M/\Lambda).
\]

Direct substitution gives the identity

\[
\boxed{
g_{\rm log}=g_{\rm PT}+\bigl[\log(M/\Lambda)+\gamma_E/2\bigr]e^{-x}+h(x).
}
\]

The heat-covariantized rank is the middle term. The remainder \(h\) is not a
multiple of the identity and is not a multiple of the heat kernel. Therefore
simply covariantizing the old rank term to a heat trace **changes the finite
action** relative to \(g_{\rm log}\) by \(\operatorname{Tr}h(D^2/\Lambda^2)\).
It also changes the finite-matrix theory relative to the old rank term by

\[
\bigl[\log(M/\Lambda)+\gamma_E/2\bigr]\operatorname{Tr}\bigl(e^{-D^2/\Lambda^2}-1\bigr).
\]

The calculated determinant modulus is \(g_{\rm PT}\) without either addition.
That is not a completed \(\Gamma_{\rm one}\). The identity is used as a
comparison of profiles, not as a replacement of the already computed source.

## Low-\(x\) limit and Seeley weights

As \(x\to 0^+\),

\[
E_1(x)=-\gamma_E-\log x+x-\frac{x^2}{4}+\cdots,
\]

so \(h(x)\to 0\). The explicit series used for \(x<0.05\) is

\[
h(x)=\sum_{n=1}^\infty\frac{(-1)^{n+1}}{2\,n!}\,x^n\Bigl(\log x+\gamma_E-\frac1n\Bigr).
\]

The leading term is \(\tfrac12 x(\log x+\gamma_E-1)\). Consequently \(g_{\rm log}\)
and \(g_{\rm PT}\) plus the heat-rank term agree for fixed modes as the cutoff
is removed, and for modes below the cutoff; their difference lives at finite
\(x=\mathrm{O}(1)\). This is not the ultraviolet of each mode.

The moments of \(h\) are exact:

\[
\int_0^\infty h(x)\,dx=-\frac12,\qquad
\int_0^\infty x\,h(x)\,dx=-\frac34,\qquad
h(0^+)=0.
\]

In the heat-kernel counting where the \(a_4\) coefficient multiplies
\(\chi(0)\), \(h\) does **not** change the logarithmic conformal density.
It does change the power-divergent volume and Einstein moments. Anomaly
matching of \(a_4\) therefore cannot decide between heat-covariantized rank
and weighted-log; that is the same non-uniqueness already stated in
Andrianov §§3–4 and in AGENTS 20.28–20.29.

Independent quadrature of \(E_1\) as \(\int_1^\infty dt\,e^{-xt}/t\) and of
those two moments is stored in the record.

## Metric, Weyl, and normalization derivatives

For a spectral function \(f(\lambda)\) on an invertible fiber,

\[
\delta\operatorname{Tr}f=\operatorname{Tr}\bigl(f'(D)\,\delta D\bigr).
\]

The metric variation of \(D_\omega\) is the existing covariant frequency
operator. Uniform Weyl rescaling of \((N,q,r)\) is the physical metric Weyl
used here; a finite Fourier sandwich need not equal that operator before
cutoff weighting.

Ordinary \(x\)-derivatives of the smooth profiles are

\[
\begin{aligned}
\partial_x g_{\rm PT}&=-\tfrac12 e^{-x}/x,\\
\partial_x\bigl(c\,e^{-x}\bigr)&=-c\,e^{-x},\\
\partial_x h&=\tfrac12 e^{-x}(\log x+\gamma_E),\\
\partial_x g_{\rm log}&=\tfrac12 e^{-x}\bigl(\log x-2\log(M/\Lambda)-1/x\bigr),
\end{aligned}
\]

with \(c=\log(M/\Lambda)+\gamma_E/2\). They satisfy the differentiated
identity. The corresponding \(\lambda\)-derivatives enter the Hellmann–Feynman
formula on each fiber.

Normalization derivatives at fixed \(D\) and \(\Lambda\) are

\[
\frac{\partial g_{\rm PT}}{\partial\log M}=0,\qquad
\frac{\partial h}{\partial\log M}=0,\qquad
\frac{\partial g_{\rm log}}{\partial\log M}=e^{-x}=\frac{\partial}{\partial\log M}\bigl(c\,e^{-x}\bigr).
\]

The hard projector instead yields the cutoff rank,

\[
\frac{\partial g_{\rm hard}}{\partial\log M}=\Theta(1-x).
\]

Under \(\mu\to\gamma\mu\), Andrianov’s partition function changes by
\(e^{-(\log\gamma)\operatorname{Tr}P_N}\). That is \(\operatorname{Tr}\Theta\),
not the heat trace. The smooth weighted-log and the primary projector
therefore do not even share the same owner of \(d\Gamma/d\log M\).

A further obstruction is metric variation of \(g_{\rm hard}\): eigenvalues
that cross the wall \(x=1\) contribute a spectral boundary term omitted by
the interior derivative \(-1/\lambda\) on \(x<1\). Independent finite
differences detect that wall. The smooth profiles have no such jump.

Uniform clock scaling obeys \(E_f[cN]=cE_f[N]\) for each spectral profile at
fixed proper scales \(M,\Lambda\). Indeed \(D_\omega[cN]=D_{\omega/c}[N]\),
so the change of frequency variable \(\omega=c\nu\) supplies the factor
\(c\). The smooth profiles share this identity. The hard profile does too
when its frequency cutoff wall is integrated correctly; coarse quadrature
of a discontinuous integrand can obscure it.

## What the covariant operator actually does

Both a constant cylinder and one nonconstant \((N,q,r)\) metric are evaluated
with the existing frequency operator, the same \(4\kappa/\pi\) counting, and
independent controls:

- analytic cylinder eigenvalues versus the finite Fourier matrix;
- \(\operatorname{Tr}e^{-D^2/\Lambda^2}\) from `expm` versus the eigendecomposition;
- the finite-matrix `RegulatedOperator` action versus \(g_{\rm PT}\) plus old rank;
- adaptive \(\omega\) quadrature of the analytic cylinder fibers versus Gauss–Legendre;
- independent Hellmann–Feynman checks on a general-metric fiber;
- centered metric finite differences of the integrated smooth energies;
- the existing `cutoff_response` energy and Weyl trace as a binding control
  on \(g_{\rm PT}\) only.

The identity-trace integral grows linearly with the frequency endpoint. The
heat trace saturates. This is the concrete obstruction to integrating the old
rank term as \(\operatorname{Tr}1\).

The remainder energy \(\operatorname{Tr}h\) is not a rounding error on either
geometry. Heat-covariantizing the old rank is therefore not an exact rewrite
of \(g_{\rm log}\), and neither object is the calculated determinant modulus
or a completed \(\Gamma_{\rm one}\). The sharp projector requires an explicit
cutoff wall; a Gauss–Legendre integral of the discontinuous integrand on a
wide frequency interval is not that formula.

## Computed response on the covariant frequency operator

All smooth-profile numbers below use \(\Lambda=2\), \(M=1\), angular maximum 8
and 32 positive-frequency nodes, on the existing 16-point AP grids. The
proper-time column reproduces `cutoff_response` to \(2\times 10^{-15}\).
Hard-log Gauss values are preserved only as coarse diagnostics.

Cylinder \(L=4\), \(a=1\):

| Profile | Energy | Neck \(\rho+p_x\) |
|---|---:|---:|
| \(g_{\rm PT}\) (calculated determinant modulus) | 4.603085860503167 | -0.003804143338 |
| heat-covariantized rank | -7.888676203741531 | \(1.13\times 10^{-6}\) |
| remainder \(h\) | -14.864141722456047 | -0.001236166910 |
| weighted-log | -18.149732065694400 | -0.005039179925 |
| hard projector log, 32-node Gauss (coarse) | -2.841515236464148 | diagnostic only |

The exact AP-cylinder hard energy, for each \(\kappa,p\) with
\(E=\sqrt{p^2+\kappa^2/a^2}<\Lambda\) and \(W=\sqrt{\Lambda^2-E^2}\), is

\[
E_{\rm hard}=\frac{4\kappa}{\pi}\bigl[(2-\ln(\Lambda^2/M^2))W-2E\arctan(W/E)\bigr].
\]

The uniform-lapse interior piece is \((4\kappa/\pi)\,2[W-E\arctan(W/E)]\); the
cutoff wall is \((4\kappa/\pi)[-\ln(\Lambda^2/M^2)W]\); they sum to
\(E_{\rm hard}\). The exact normalization derivative is
\((8\kappa/\pi)W\). On this grid the two retained modes give

| Quantity | Value |
|---|---:|
| exact \(E_{\rm hard}\) | -3.2978447242754125 |
| interior | 2.1518382062568135 |
| wall | -5.4496829305322265 |
| \(dE_{\rm hard}/d\log M\) | 7.862230538296076 |

An independent adaptive \(\omega\) integral of \(-\log((\omega^2+E^2)/\Lambda^2)\)
on \([0,W]\) plus the explicit shell \(-\ln(\Lambda^2/M^2)W\) reproduces those
values. Gauss–Legendre on the operator's \([0,16]\) map does not:

| Frequency nodes | Hard energy | Error versus exact |
|---:|---:|---:|
| 32 | -2.841515236464143 | 0.456329487811270 |
| 64 | -3.526168496206960 | -0.228323771931548 |
| 128 | -3.467784650665305 | -0.169939926389892 |

Convergence of that discontinuous integrand is not assumed. The 32-node
operator value is kept so the error against the closed form remains visible.

General smooth \((N,q,r)\):

| Profile | Energy | Neck \(\rho+p_x\) |
|---|---:|---:|
| \(g_{\rm PT}\) | 8.527482803351840 | 0.003978152861 |
| heat-covariantized rank | -14.323647014747939 | -0.012003049523 |
| remainder \(h\) | -26.803003706342647 | -0.018598459963 |
| weighted-log | -32.599167917738725 | -0.026623356625 |
| hard projector log (coarse Gauss / interior derivative) | -8.328495023761750 | not a physical stress |

The identity \(g_{\rm log}=g_{\rm PT}+\text{heat rank}+h\) holds on both
integrated spectra to floating-point noise. \(\operatorname{Tr}h\) is larger
than the proper-time modulus itself, so the finite action really does change.

A mean-zero Weyl probe on the cylinder is consistent with zero for every
smooth profile (homogeneous background). On the general metric the same probe
gives physical metric Weyl variations

| Profile | \(\delta_\sigma E\) |
|---|---:|
| \(g_{\rm PT}\) | -0.906949047943905 |
| heat-covariantized rank | 1.520034779856396 |
| remainder \(h\) | 2.842186870717168 |
| weighted-log | 3.455272602629659 |
| hard projector log | coarse interior-derivative diagnostic, not physical stress |

Independent centered metric differences of the mixed \((N,q,r)\) probe match
the smooth spectral gradients to \(5\times 10^{-8}\) or better. Existing
general-metric hard-log energies, neck values, Weyl numbers and interior
Hellmann–Feynman gradients remain unresolved coarse-quadrature / missing-wall
diagnostics. They are not physical stress comparisons. The exact wall term is
owned by the product cylinder formula above.

Normalization derivatives of the smooth profiles on the general metric are
\(dE_{\rm PT}/d\log M=dE_h/d\log M=0\) and
\(dE_{\rm log}/d\log M=dE_{\rm heat}/d\log M=35.40730235933=\operatorname{Tr}e^{-D^2/\Lambda^2}\).
The exact cylinder hard derivative is \((8\kappa/\pi)W\), not the coarse
projected-rank Gauss value. Uniform clock scaling holds for all four smooth
profiles to \(1.5\times10^{-14}\) in the stored controls.

With angular maximum 4, the identity-trace integral is linear in the
frequency endpoint (ratio \(407.436654315252\) at extents 16, 24 and 32)
while the heat trace stays at \(19.394183\). That is the obstruction to
integrating the old rank as \(\operatorname{Tr}1\).

On one cylinder fiber the finite-matrix `RegulatedOperator` action equals
\(g_{\rm PT}\) plus old rank exactly; old rank minus heat-rank is
\(-12.35099411955739\). Analytic cylinder eigenvalues match the Fourier
matrix to \(10^{-14}\); `expm` matches the heat trace to \(10^{-15}\).
Adaptive \(\omega\) quadrature of the \(\kappa=1\) analytic fibers reproduces
the Gauss–Legendre \(\kappa=1\) proper-time energy
\(2.60518485999\) to the reported quadrature error \(2.6\times 10^{-8}\).

The exact identity residual over the stored \(x\)-grid is at most
\(1.11\times 10^{-16}\). Independent \(E_1\) integrals and the moments of
\(h\) agree with \(-\tfrac12\) and \(-\tfrac34\).

AGENTS 20.28–20.29 keep the regulated determinant and its induced heat as one
object: the bosonic spectral action is the scale anomaly of the same fermionic
operator, not a second weighted action. Replacing the calculated determinant
modulus \(g_{\rm PT}\) by \(g_{\rm log}\) would mix a different finite
prescription into that object without computing the compensator and would not
complete \(\Gamma_{\rm one}\). The finite local basis
\(M^4,M^2\mathcal R,C^2,\mathcal R^2,E_4,\Box\mathcal R\) and the common
compensator remain unknown. Neck-null values of the smooth profiles are
stored as diagnostics; they are not solved for a stationary throat.

The real-time Gaussian influence functional is a separate owner. This
calculation is Euclidean-static.

## Reproduction

```sh
python3 scripts/check_nsc_measure_normalization.py --check
PYTHONPATH=src:. python3 -m unittest discover -s tests -p test_nsc_measure_normalization.py -v
```

`--output` and `--check` are mutually exclusive. Existing records cannot be
overwritten. The all-field comparison includes every sample, control, scope
flag and source hash. Float tolerances are atol \(2\times 10^{-7}\) and
rtol \(2\times 10^{-8}\); primitive identity and fiber gates are tighter.

The next owner remains the invariant completion and coupled state of the
**same** calculated determinant modulus already used for the subtracted
source, not a new finite profile chosen to manufacture a root.
