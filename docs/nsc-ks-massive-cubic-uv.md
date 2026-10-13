# Massive cubic coefficient on the directed integrator

This is the generic-mass extension of the massless
[reduced current](nsc-ks-cubic-uv-current.md) and its
[directed enclosure](nsc-ks-cubic-uv-enclosure.md). The action, original
source sectors and history `0b0e4ced` stay fixed. The massless owners are
not edited. The result is formal $e^{-3}$ mathematics plus one directed
point. It is not a UV tail, not coverage of $I$, and not a closed gate.

Owner: `src/recursive_horizons/nsc_ks_massive_cubic_uv.py`. Record:
`results/development/nsc-ks-massive-cubic-uv-v1.json`.

## Parent objects from the owned operator

For $s=\pm1$, $a=a(\rho)$, $b=a\ell/(2r)$ and $d=am/2$,

$$
l=d-isb,\qquad
L=T_s=\partial_\rho-\frac{s}{a^2}\partial_z,\qquad
M=T_{-s}.
$$

The occupied column has major $A_1=iq$ and minor $A_1=l$. Actual $L_0$
gives the minor coefficients

$$
B=i\Bigl(ql-\frac{a^2}{2}M(l)\Bigr),\qquad
C=-i\frac{a^2}{2}M(B)+l(R+ih),
$$

with $R=(sq_z+k-q^2-d^2-b^2)/2$ and arbitrary common transported $k$.
Major transport is $L(q)=-m^2/2-\ell^2/(2r^2)$. For $m\neq0$,

$$
L(h)=-\frac{m\ell}{4}\frac{sa^2 r_\rho+r_z}{r^2},
$$

which is not the massless equation $L(h)=0$. That equation is not
substituted. Continuity at $e^{-4}$ still reads

$$
L(n_4)=-\frac{s}{a^2}\partial_z\bigl(2|B|^2+4\operatorname{Re}(\overline l\,C)\bigr),
$$

and $L(h_3)=-2\operatorname{Re}(\overline l\,C)/a^2$.

## Cancellation without setting $h$ to zero

The raw interference depends on $h$ through $q_z$ and on $h_z$ through
$q$. Both contributions are real, so

$$
J_3=h_{3,z}+Rq_z-qR_z+\operatorname{Im}(\overline l\,B_z+\overline B\,l_z)-sn_4
$$

is independent of $h$ and $h_z$. The same $h$ multiplies $l$ inside
$C$, and $\operatorname{Re}(\overline l\,C)$ cancels it before any
transport of $h$. Combining $L(h_{3,z}-sn_4)$ with this current gives,
for both signs,

$$
L(q_z)=-\frac{4bb_z}{a^2},\qquad
L(q_{zz})=-\frac{4(b_z^2+bb_{zz})}{a^2},
$$

$$
L(J_3)=L(J_3)_{\mathrm{massless}}+m^2(sq_{zz}-4bb_z).
$$

The massless right-hand side is the existing seven-term formula. On
$\Sigma$, where $a$ and $r$ match the reference and $b_z=b_{zz}=0$,

$$
\begin{aligned}
\Delta N_{\mathrm{bracket},3}
&=\frac{s}{a}\Delta J_3+\frac{a^3}{2}\Delta(b_\rho^2)
-\frac{4sb^2}{a}\Delta q_z\\
&\quad+sab\,\Delta b_{\rho z}-sam^2\Delta q_z.
\end{aligned}
$$

Common $k$ cancels. The action remains
$N=-\mu/(2\pi)N_{\mathrm{bracket}}$ and $\beta=+\mu/(2\pi)J_3$.
Every residual in `massive_cubic_identities` is the exact string `0`.

## Three-power factorization

Write $b=\ell b_0$, $q_z=\ell^2 Q_z$ and $q_{zz}=\ell^2 Q_{zz}$. The
triangular system and the Sigma coefficient split as

$$
J_3=\ell^2 J_2+\ell^4 J_4+m^2\ell^2 J_M,
$$

with the same powers in $N_{\mathrm{bracket}}$,

$$
\begin{aligned}
Q_z'&=-4b_0 b_{0z}/a^2,\\
Q_{zz}'&=-4(b_{0z}^2+b_0 b_{0zz})/a^2,\\
J_2'&=g_2,\\
J_4'&=-16b_0^3 b_{0z}/a^2+(4sb_0^2/a^2)Q_{zz},\\
J_M'&=-4b_0 b_{0z}+sQ_{zz}.
\end{aligned}
$$

Here $g_2$ is the quadratic geometric part of the massless current.
Upstream differences remain the homogeneous zero representative, so the
splitting survives integration. Both characteristic signs are still
required. This can later replace sixty separate $(m,\ell)$ integrations
for one history. It has not been used as a family campaign, and it does
not create a complete $C_M$.

## Reuse of the directed integrator

`MassiveCubicGeometry` is a subclass of the historical `CubicGeometry`.
It calls the existing forcing jets and adds $m^2 a(\rho)^2$ times the
$q_z$ forcing to the current, and $s m^2$ to the coupling. The same
additions are made on the uniform cell forcing. `enclose_characteristic`
is called unchanged. The returned surface bracket is then corrected by
$-s a m^2 q_z$. Exact Arb inputs enclose the original $\pi/2$ and
$\sqrt5$ labels. At $m=0$ the additions vanish and the historical
owner is recovered.

## One group-14 point

The recorded pilot is original group 14 only: positive compact mass
$\pi/2$, $|\ell|=\sqrt5$, ledger multiplicity $\mu=12$, and
quadrature cutoff 160 on both angular signs. Both characteristic signs
are enclosed at the history's incoming center, with 1024 cells and Taylor
order 8. The angular pair contributes the usual factor of two because the
coefficient is even in $\ell$. The formal integral of that coefficient
from the cutoff, divided by $2E_c^2$, is stored only as a labeled
non-tail. No physical budget entry is filled.

## What remains open

The higher-order remainder $C_M$, uniform $C_4$ on $I$, the physical
source column, every other original family, and the local gate stay
unbounded. One point and the formal $e^{-3}$ term are not the UV tail.
The massless evaluator still rejects $m\neq0$.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_massive_cubic_uv.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_massive_cubic_uv.py --check
```


## Minor-A3 primitive and transported upstream normalization

For the actual $L_0=\partial_\rho-S_3a^{-2}\partial_z-iV/a$, the minor
recurrence is

$$A_{j+1,-s}=\frac{a^2}{2i}T_{-s}A_{j,-s}
-\frac a2V_{-s,s}A_{j,s},\qquad V_{-s,s}=-m+is\ell/r.$$

It gives

$$A_{2,-s}=\frac{s\ell a}{4r^2}
(-a_\rho ar+a^2r_\rho+sr_z+2qr)
-\frac{ima}{4}(aa_\rho-2q),$$

$$A_{3,-s}=\frac{a^2}{2i}T_{-s}A_{2,-s}
-\frac a2V_{-s,s}(R+ih).$$

The Stage-2 primitive on `0b0e4ced` bounded second geometry/phase jets,
minor $A_2$, and the jet part of minor $A_3$ for original massless group1,
$|\ell|=\sqrt5$. It uses the same normal-window constants
$|\chi'|\le8/\mathrm{width}$ and
$|\chi''|\le105/\mathrm{width}^2$, with $s_\rho=-1/a$ on
$\rho\in[1,33/32]$. The shared unowned upstream $R_{\rm up}+ih_{\rm up}$
still multiplies $\delta(1/r)$ in the whole-slab minor difference. Massless
$T_sh=0$ does not set that datum to zero. Complete minor $A_3$, numerical
$C_4,C_M$ and the vacuum tail remain null. The massless difference-bound
API rejects massive channels; the generic recurrence above retains mass.

Continuity gives
$T_sn_2=sT_sq_z=s\ell^2r_z/r^3$, so
$k=n_2-sq_z$ is transported upstream data. Both histories share it;
it is not zeroed. The lowest raw coefficients are

| Bilinear | $e^1$ | $e^0$ | $e^{-1}$ |
|---|---:|---:|---:|
| $J_I$ | $-s$ | $0$ | $-sk$ |
| $A^\dagger VA+J_{S3}/a$ | $-1/a$ | $0$ | $-aQ/2-k/a$ |

All three history differences vanish on $\Sigma$, where the intrinsic
radii and axial scales agree. A common spatially homogeneous upstream
factor $f(e)=1+c_2e^{-2}+c_3e^{-3}+c_4e^{-4}+O(e^{-5})$ consistently
transforms every transported coefficient. For either bilinear,

$$J_3'-J_3=2\operatorname{Re}(c_2)J_1
+2\operatorname{Re}(c_3)J_0
+(|c_2|^2+2\operatorname{Re}(c_4))J_{-1}.$$

Its history difference therefore vanishes on $\Sigma$. Varying $R$ alone
would miss the associated $h_{3,z}$ and $n_4$ shifts: for a real common
$c$ their contributions are $cq_z+cq_z-2cq_z=0$. This is a formal
occupied-vacuum coefficient invariance, not physical source matching or a
numerical higher-order remainder. Spatially varying factors or unmatched
incoming radii lie outside it.

The retired cubic-residual demonstration applied the shared production
channel basis `(N2,N4,NM,J2,J4,JM)` only at completed exact saved value
nodes. Missing nodes stayed empty; a center value was never broadcast.
No field/source evolution was performed, and $C_M$, uniform $C_4$ on
$I$ and the complete UV tail stayed null. A changed finite-cutoff nodal
residual does not close the gate.
