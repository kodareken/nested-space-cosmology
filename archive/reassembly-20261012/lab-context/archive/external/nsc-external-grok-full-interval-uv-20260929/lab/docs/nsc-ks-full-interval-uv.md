# Full-interval UV remainder attempt

This is one attempt to bound the history-minus-reference \(L_0 A_4\) integrals
on the whole interval \(I=S(1)+[0.12,0.18]\) of the live profile `0b0e4ced…`.
It does not produce a numerical \(C_M\) or an induced \(N/\beta\) remainder.
No exclusive result file is published. The local gate stays OPEN.

Owner: `src/recursive_horizons/nsc_ks_full_interval_uv.py`.
Driver: `scripts/derive_nsc_ks_full_interval_uv.py`. `--record` refuses to
write a file while \(C_M\) is null.

## What was evaluated

The upstream column is the live phase-fixed vacuum spinor at
\(\rho_{\mathrm{up}}\), order 5, 16 Darboux cells, group 1, both angular
signs, cutoff 160. The chart makes the major component real. It is not the
zero solution of a vanishing imaginary transport right-hand side.

For the positive square root of 5 the order-2 and order-4 major coefficients
are about \(-0.204558058653\) and \(+0.376262327369\). Neither interval
contains 0. The normalization \(A_{2,\mathrm{major}}+|A_{1,\mathrm{minor}}|^2/2\)
contains 0. The opposite angular sign is also strictly negative at order 2.
The degree-4 pointwise tail upper from this order-5 call is about
\(0.16447\). That is a valid upper for \(E\ge 160\), not a sharp order-8
constant, and its axial derivatives in the envelope are zero because the
slice is homogeneous. The carrier derivative is not included there.
Thermal and coherent Fermi occupations stay out of this column.

The existing jet call is then made on a midpoint-radius ball containing the
closed interval \(I\). A proper subinterval of the flat plateau, shrunk by
\(2^{-26}\) and kept off the center, does return a finite jet. On \(\Sigma\)
its \(b_z\) contains 0, as required when the normal coordinate vanishes.

Owned commutator integrals on this history are about
\(B_z=1.2063840260691774\times 10^{-3}\) and
\(B_{zz}=0.9555571601910049\). They do not certify \(L_0 A_4\). Supplying them
to the conditional \(M=4\) majorant with the \(L_0 A_4\) integrals absent
leaves \(C_M\) null. The manufactured omission control still fails to enclose
when the residual is dropped. The 60-family moment fold is unchanged and
still has no \(C_M\).

## Domain failure

Two obstructions stop a single characteristic jet of all of \(I\).

The axial inner radius is the binary value \(0.03\). Any Arb ball whose
radius is that value is rounded outward. The ball used here has radius about
\(0.03000000002794\), so `abs(z-center) <= inner` is false. The ball also
meets both sides of the plateau center, and
`local_axial_interval_series` raises

```text
SubdivisionNeeded: z ball straddles the plateau center; |z-center| is not analytic
```

Separately, the characteristic displacement from \(\rho=1\) to
\(\rho_{\mathrm{up}}\) is about \(0.04328030266673\). The flat-to-support
margin is \(0.03\). The image of the closed interval \(I\) under either sign
therefore leaves the axial support. One analytic plateau chart cannot cover
that preimage.

## The term that is still unconstructed

On shell the truncation defect is

\[
\left(L_0-\frac{2ie}{a^2}\Pi_{-s}\right)\sum_{j=0}^{4}e^{-j}A_j
=e^{-4}L_0 A_4,
\]

with

\[
\Pi_{-s}L_0 A_4=\frac{2i}{a^2}\Pi_{-s}A_5,\qquad
A_{5,-s}=\frac{a^2}{2i}T_{-s}A_{4,-s}-\frac{a}{2}V_{-s,s}A_{4,s}.
\]

\(C_M\) needs the \(H^2\) integrals of that spinor and of its first two axial
derivatives along both characteristics through the preimage of all of \(I\),
after history-minus-reference subtraction. Those integrals are not returned
by the jet call above. The triangular \((Q_z,Q_{zz},J_2,J_4,J_M)\) transport
is the \(e^{-3}\) current, not this defect. \(B_z\) and \(B_{zz}\) multiply
the propagated field; they are not the defect.

Upstream invariance cancels a common scalar \(A_2\) in the cubic \(J_3\) and
\(N\) difference on \(\Sigma\). The same coefficient reenters minor \(A_3\).
The owned minor-\(A_3\) identity still depends on major \(A_2\). The massless
statement \(\operatorname{Im}(T_s A_2)=0\) does not set the upstream value to
zero. At \(\rho=1.02\), on the ray through the center of \(I\), the live
history and the zero history have \(\Delta b\) about \(+1.71914\times 10^{-7}\),
separated from zero, and \(R_{\mathrm{up}}\Delta(1/r)\) about
\(-3.79492\times 10^{-8}\), also separated from zero. Replacing
\(R_{\mathrm{up}}\) by zero deletes that term. This ray is not the integral
over \(I\). The crude kernel majorant times \(|R_{\mathrm{up}}|\) is about
\(2.35323\times 10^{-5}\). That is an upper bound on one algebraic coupling,
not \(|\delta A_3|\) and not \(C_M\). Zero deformation makes the kernel and
the minor-\(A_3\) difference exact zeros without setting \(R_{\mathrm{up}}\)
to zero.

## Carrier weight, not a gate comparison

For a \(z\)-independent envelope error the restored axial field satisfies
\(F_z=EF\). An \(E^{-4}\) field error therefore contributes an \(E^{-3}\)
current. With the order-5 tail constant \(C\),

\[
\int_{160}^{\infty}\frac{C}{E^3}\frac{dE}{2\pi}
=\frac{C}{4\pi\cdot 160^2}
=\frac{C}{2\pi}\cdot\frac{1}{2\cdot 160^2}.
\]

The second factor is the same cutoff weight as the owned leading moment of
one family at \(E_c=160\). The numerical upper is about \(5.11265\times 10^{-7}\).
The provisional allocation \(5\times 10^{-12}\) is smaller. That comparison
does not say the difference remainder exceeds the allocation, and it does
not say the remainder is physically absent. The two histories share this
upstream tail. The history-minus-reference piece is the unconstructed
\(L_0 A_4\) integral. A loose absolute-jet sum is not used in its place.

## Cost

One run of the attempt used about \(0.55\) CPU seconds and \(143\) MiB.
No 129-node march and no second heavy process were started.

```sh
OPENBLAS_NUM_THREADS=1 python3 scripts/lab.py -m pytest -q tests/test_nsc_ks_full_interval_uv.py
OPENBLAS_NUM_THREADS=1 python3 scripts/lab.py scripts/derive_nsc_ks_full_interval_uv.py --check
```
