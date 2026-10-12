# Generated massless cubic N,beta coefficient

This is the next Stage-2 UV construction after the accepted
[surface upstream invariance](nsc-ks-uv-upstream-invariance.md). Common
spacetime-constant scalar upstream coefficients drop out of the formal
occupied-vacuum C4 **difference on Sigma** when higher jets are transformed
together. They do not replace the generated characteristic contribution to
$A_3$ and $n_4$. This page constructs that contribution for the original
massless group-1 pair and contracts it with the actual Pauli vertices.

The action, preparation, source inventory and KS subtraction are unchanged.
The local gate stays OPEN. Historical UV remainder and invariance records are
not rewritten.

Owner: `src/recursive_horizons/nsc_ks_cubic_uv_current.py`. Record:
`results/development/nsc-ks-cubic-uv-current-v1.json`.

## Massless reductions from the owned L0 recurrence

On $L_0=\partial_\rho-S_3 a^{-2}\partial_z-iV/a$ with $V=\ell S_2/r$ and
$T_s=\partial_\rho-s a^{-2}\partial_z$, the massless original pair uses

$$
b=\frac{a\ell}{2r},\qquad
l=-isb,\qquad
B=sbq-\frac{s a^2}{2}b_\rho-\frac12 b_z,
$$

$$
n_2=s q_z+k,\qquad
R=\frac{s q_z+k-q^2-b^2}{2}.
$$

$k=n_2-s q_z$ is common transported data. It is not assigned to zero as a
physical statement; $k=0$ is only a representative of the invariance class.
$R$ is never zeroed independently of $h_3$ and $n_4$.

Minor $C=A_{3,-s}$ is the vacant recurrence, not a local metric jet:

$$
C=sbh-i\Bigl(\frac{a^2}{2}T_{-s}B+sbR\Bigr).
$$

All linked terms are kept. The major-$A_3$ and density transports are

$$
T_s h_3=-sb\,T_{-s}B-\frac{2b^2 R}{a^2},\qquad
T_s n_4=-\frac{s}{a^2}\partial_z\bigl(2B^2+4\operatorname{Re}(\overline l C)\bigr).
$$

On-shell $T_s q=-\ell^2/(2r^2)$ and $[T_s,\partial_z]=0$. Massless
$T_s h=0$, so $h=h_{\mathrm{up}}$ is a common constant; the identity
current $J_I$ at $e^{-3}$ does not see it. Representative ICs at
$\rho_{\mathrm{up}}$ are the owned vanishing major $A_1$,
$h_{3,z}=0$ and $n_4=0$.

## Actual vertices, not the old v3 shorthand

The raw KS vertices and `source_column_matter` fix the action signs:

$$
N=-\frac{\mu}{2\pi}\Bigl(A^\dagger V A+\frac{J_{S_3}}{a}\Bigr),\qquad
\beta=+\frac{\mu}{2\pi}J_I.
$$

The raw beta momentum vertex is $-I$ and the Gaussian action minus is
applied once, so $\beta=+J_I$. The older remainder shorthand
$\beta_3=-\mu/(2\pi)J_3$ is not this identity. Columns already contain
$\sqrt{dE/(2\pi)}$; a second measure is not inserted. Direct Pauli
substitution gives

$$
J_{I,3}=h_{3,z}+R q_z-q R_z+s(b B_z-B b_z)-s n_4.
$$

A common real scalar shift $R\to R+c$ also sends $h_{3,z}\to h_{3,z}+c q_z$
and $n_4\to n_4+2cs q_z$. Their contributions cancel. Holding $h_3$ and
$n_4$ fixed leaves the spurious $c q_z$.

On $\Sigma$, the declared radius history vanishes, so $r_g=r_{\mathrm{ref}}$
and $a_g=a_{\mathrm{ref}}$. The three lower current differences vanish there.
The cubic difference is the generated bulk transport, contracted on
$I=S(1)+[0.12,0.18]$, summed over both original angular signs with equal
$\mu=6$.

## What the diagnostic quadrature measures

The evaluator marches $q$, $q_z$, $q_{zz}$, $h_{3,z}$ and $n_4$ along
$T_s$ characteristics by RK4, using algebraic $z$-jets of $B$, $C$ and
the $n_4$ source (no dummy $L_0 A_j$). History-minus-reference is formed
on $I$. Zero deformation gives exact zero without setting $R=0$. The
$e^{-3}$ coefficient is even in $\ell$; the odd remainder of the pair is
zero at working precision.

The recorded numbers are a **measured** coefficient. They are not a
remainder-validated $C_4$. The 24-versus-16 step discrepancy is a numerical
control, not a directed quadrature remainder. No $L_0 A_4$ $H_2$ integrals
are owned, so $C_M$ stays `null`. The physical source column and the local
gate stay OPEN.

## Reduced system for the next directed quadrature

Combining the same major transport and continuity equations eliminates the
separate $h_{3,z}$ and $n_4$ variables from the current evolution. With
$J_3=J_{I,3}$ and the same $b=a\ell/(2r)$, six exact symbolic controls verify
the following reduction for both $s=\pm1$:

$$
T_s q_z=-4bb_z/a^2,\qquad
T_s q_{zz}=-4(b_z^2+bb_{zz})/a^2,
$$

$$
\begin{aligned}
T_s J_3={}&-a^2 b b_{\rho\rho z}+a^2 b_\rho b_{\rho z}
-2aa_\rho b b_{\rho z}\\
&-s b b_{\rho zz}+s b_z b_{\rho z}
-16b^3b_z/a^2+4s b^2q_{zz}/a^2.
\end{aligned}
$$

Thus only $q_z,q_{zz},J_3$ need transport. The system is affine and triangular:
geometry supplies the first two right-hand sides, and the third depends
linearly on $q_{zz}$. It can be written as single and nested Volterra integrals.
It still contains retarded history dependence; it is not an instantaneous
metric approximation. Upstream values are zero in the same z-homogeneous
representative already used above.

On the incoming slice $a_g=a_{\rm ref}$ and $b_g=b_{\rm ref}=b$, while all pure
z derivatives of that common radius vanish. Direct Pauli contraction gives

$$
\Delta N_{\rm bracket,3}
=\frac{s}{a}\Delta J_3
+\frac{a^3}{2}\Delta(b_\rho^2)
-\frac{4s b^2}{a}\Delta q_z
+sab\,\Delta b_{\rho z}.
$$

Multiply by $-\mu/(2\pi)$ for the action's N contribution; beta is
$+\mu\Delta J_3/(2\pi)$. The common $k$ cancels explicitly here. This identity
holds for a common constant imaginary major-A2 datum as well; that datum is
retained in the symbolic contraction. The existing recorded numerical pilot
still uses the full $h_3/n_4$ transport as an independent construction. A
directed quadrature and higher-order UV remainder have not yet been supplied.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_cubic_uv_current.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_cubic_uv_current.py --check
```

The direct-coordinate mixed jets are checked against finite differences, and
the transported phase and its z derivative are compared with the existing
independent split-Gauss phase owner for both characteristics. These controls
check signs and implementation, without promoting quadrature to a certificate.
The normal-window derivative uses separately evaluated logistic factors, so
rounding the window to one cannot erase a nonzero ramp derivative.
