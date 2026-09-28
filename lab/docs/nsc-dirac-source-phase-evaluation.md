# Numerical owner of the formal source-phase coefficient

The [transport identity](nsc-dirac-source-phase-transport.md) supplies the
proven formal coefficient

\[
\delta n_{2,s}(z)=s\,\partial_z f_s(z),\qquad
f_s(z)=\frac{\ell^2}{2}\int_1^{\rho_u}
\bigl(r_g^{-2}-r_{\rm ref}^{-2}\bigr)
(\rho,z-s D(\rho))\,d\rho,
\]

with \(D(\rho)=\int_1^\rho a(\rho')^{-2}\,d\rho'\) on the original chart.
This module evaluates \(f_s\), \(f_{s,z}\) and the full amplitude tangents
\(\delta f_s\), \(\delta f_{s,z}\) for both original source signs. It is
geometry integration of that boxed formula. The uniform remainder and
ordered limit are established by the separate
[source-cutoff bridge](nsc-source-cutoff-bridge.md); numerical tail accuracy
is still required. This geometry owner does not evolve source fields.

## Reused geometry

The radius family is the existing
`LocalIncomingFamily` or `CompatibleIncomingMetric` pure-radius history.
The original chart supplies \(T(\rho)\), \(a(\rho)\) and
\(r_{\rm ref}=\sqrt{1+\rho^2}\). The normal coordinate is \(s=T(\rho)-T(1)\).
The characteristic length \(D\) is the owned continuum speed distance of
that same chart. Because \(T\) and \(a\) are frozen, there is no
\(\delta D\). Normal windows and the actual \(w,U\) profiles are retained;
zero-amplitude directions stay in the tangent.

Caller data are \(\rho_u\ge 1.03\), the original angular \(\ell\), and the
target \(z\). The analytic profile identity includes every direction and
normal window. The coefficient binding adds \(\rho_u\), \(\ell\), target
\(z\), source-sign order \((+1,-1)\) and the Gauss--Legendre count. No
value cache is used.

## One quadrature path

On each panel the four integrands are accumulated together:

\[
\begin{aligned}
f_z&=-\ell^2\int r_z/r^3,\\
\delta f&=-\ell^2\int \delta r/r^3,\\
\delta f_z&=-\ell^2\int\bigl(\delta r_z/r^3-3 r_z\,\delta r/r^4\bigr).
\end{aligned}
\]

Axial derivatives of \(w\) and \(U\) are analytic. The value integrand uses
the stable difference of reciprocals
\(r^{-2}-r_{\rm ref}^{-2}=-\delta r(2 r_{\rm ref}+\delta r)/(r^2 r_{\rm ref}^2)\).
Known normal-window inner/outer radii are inverted through \(T(\rho)\) and
split the \(\rho\) interval; each panel uses a configurable Gauss--Legendre
count. Absolute panel contributions are reported separately; they measure
neither convergence nor quadrature error. A convergence indicator requires
a separate resolution comparison. The perturbation is retained before its
addition to the reference radius, including below floating-point spacing.
The owned analytic coefficient bound must keep radius positive throughout
the slab, and the history and all tangent directions must be flat at the
fixed upstream slice. Nodes and panel ends are also checked.
\(\ell=0\) returns exact zeros without quadrature. A zero-amplitude history
has \(f=f_z=0\) but keeps a nonzero tangent when the profiles are nonzero.

Returned missing coefficient, remainder and physical error bounds are
`None`. The physical local gate stays `OPEN`. No term is added to \(N\) or
\(\beta\) assembly.

## Shapes

Axis 0 is `SOURCE_SIGNS=(+1,-1)`. Amplitude axis 1 follows the bound
metric (`w` then `U` for a `LocalIncomingFamily`). Axis `-1` is the caller
\(z\) order.

## Focused geometry checks

The tests finite-difference both \(w\) and \(U\) on a nonzero eight-mode
local family, finite-difference \(f\) in \(z\), cover both characteristics,
check the leading odd reflection \(\alpha\to-\alpha\) rather than false
finite-\(\alpha\) evenness, distinguish analytic profiles, and freeze the
same upstream \(D(\rho_u)\). They do not run a source generator.

Seven focused checks also cover a perturbation of amplitude1e-20 without
subtractive loss and rejection of a nonpositive radius bound or nonflat
upstream history. Missing coefficient/remainder/physical bounds remain
`None`. The [joined evaluator](nsc-ks-cutoff-bridge-evaluator.md) separately
checks the full N,beta derivative after the same-action contraction.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_dirac_source_phase.py
```
