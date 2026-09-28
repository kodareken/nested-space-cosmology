# A matched warp-conversion source on the actual neck

A concrete nonzero contribution to the remaining source is now evaluated
on the actual black-universe geometry. It comes from the already matched
vacuum and Weyl-squared terms of the warp/covariant conversion J. No
coefficient is fitted to the required stress, and the existing canonical
Dirac tensor is reused.

This is an identified summand of the full functional. The nonlocal
remainder is retained explicitly and has no established small bound on
this geometry. The result is not an approximation to the full source or
a self-sourced solution.

The [record](../results/development/warp-local-neck-source.json) is generated
by [derive_nsc_warp_local_source.py](../scripts/derive_nsc_warp_local_source.py).
It uses the [stored coefficient allocation](nsc-canonical-spectral-bridge.md),
the existing [general metric variational density](nsc-finite-terms.md), and
the [PG source-frame map](nsc-adm-neck-source-map.md).

## Keep the whole functional while evaluating one part

Write the exact subtraction identity

\[
\mathcal J_\Lambda=\mathcal J_{\le4}+\mathcal J_{\rm remainder}.
\]

The inherited local terms, in the original metric convention, are

\[
S_{J,\le4}=-\int\sqrt{-g}\,
[V_J+A_J R_L+C_{W,J}C^2+C_{E,J}E_4+C_{\Box,J}\Box R_L].
\]

The free Dirac contribution supplies no independent R_L² coefficient in
this combination. This does not set other sectors' finite R_L² terms to
zero. Euler and total-derivative terms give no bulk source at the interior
neck; their actual boundary variations remain part of the full action.

For the stored Lambda=2, ell=2 and warp s=1,

\[
V_J=-0.00528370024139,\qquad
A_J=-0.00011105299768,\qquad
C_{W,J}=0.00000420159048.
\]

Their dimensions are respectively length^-4, length^-2 and one. They are
read from the prior record, not estimated from the source to be matched.
The coefficient A_J belongs once in the complete Einstein coefficient
a_EH. Its -2 A_J G tensor is therefore excluded from the source-side sum
reported below. The order-reduced curvature EFT contact is not added to
this calculation in the original metric.

## Vary the full metric before imposing the profile

Use the existing general static density N q r² C² and differentiate N,
q and r independently, including their first and second radial derivatives.
Only afterwards impose

\[
r(\rho)=\sqrt{1+\rho^2},\qquad
A(\rho)=1+3\rho+3(1+\rho^2)(\arctan\rho-\pi/2),
\qquad N=\sqrt A,\quad q=1/\sqrt A.
\]

At the interior neck A<0. The rational covariant variation continues to
that chart with Nq=1; its intermediate static time/radial labels must be
converted to the real child frame. In particular rho_child=-p_static
and p_parallel,child=-rho_static. No physical imaginary stress is inferred.

For unit coefficient of C² this gives the exact child-frame source

\[
\boxed{
(\rho,T_{01},p_\parallel,p_\perp)_{C^2}
=\left(-48,\ 0,\ \frac{112}{3}+16\pi,
                 -\frac{128}{3}-8\pi\right).
}
\]

Its trace is zero. This is a metric variation of the complete local
invariant, not a derivative along the imposed one-parameter geometry.
The vacuum term contributes (V_J,0,-V_J,-V_J).

Transporting the Weyl-squared source to the PG normal frame gives

\[
(\rho,T_{01},p_r,p_\perp)_{C^2,\rm PG}
=\left(-\frac{112}{3},-\frac{32\beta}{3},
                    48+16\pi,-\frac{128}{3}-8\pi\right),
\qquad\beta^2=3\pi/2.
\]

The actual contribution is C_W,J times this tensor. Both frames refer
to the same geometry and the same source normalization.

## The new contribution to the four source equations

With V=V_J and c=C_W,J, the counted source-side force vector is

\[
\boxed{
F_{J,\rm counted}=
\begin{pmatrix}
4\pi(V-112c/3)\\
-128\pi\beta c/3\\
4\pi[V-(48+16\pi)c]\\
8\pi[V+(128/3+8\pi)c]
\end{pmatrix}.
}
\]

This supplies nonzero lapse, momentum and pressure contributions with no
new scalar field. The term proportional to A_J is kept separately as an
Einstein-coefficient increment, so it is not counted twice.

Its parent Killing power vanishes:

\[
-\beta F_N-(\beta^2+1)F_\beta+\beta F_q=0.
\]

The nonzero mixed component in the PG frame is therefore not a new
injection of Killing energy. Adding this static local contribution leaves
the previously computed parent power unchanged. It cannot by itself
balance that flux in an exactly static solution.

The remaining equation is now

\[
F_{\rm canonical}+F_{J,\rm counted}+F_{\rm unknown}
=F_{\rm required}(a_{\rm EH}).
\]

The record contains this known subtotal in the same PG frame. Unknown
still includes J_remainder, the remaining vacuum conversion, other
canonical compact sectors and the applicable gauge/interface terms.
No remainder is assigned the desired residual by definition, and no
a_EH is selected from one component. The local contribution cannot be
used as a controlled truncation of the complete action at this neck.

## Reproduction and scope

```sh
python -B scripts/derive_nsc_warp_local_source.py --check
```

This differentiates the previously available general density, applies the
actual profile, and combines authenticated coefficients and tensor data.
It checks the source trace and frame/flux accounting and compares all
record fields. It does not run an old spectrum, state, angular-stress or
geometry generator. The nonzero source above is the completed calculation;
full physical closure still requires the remaining forces and their joint
geometry solution.
