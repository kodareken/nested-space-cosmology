# Actual Weyl defect for the finite transmitting history

This owner connects the [formal spatial projector](nsc-spatial-reference-symbol.md)
to the actual Weyl-operator commutator. The projector recurrence is reused
from [Panati–Spohn–Teufel](https://arxiv.org/abs/math-ph/0201055), with the
project's existing sign convention. No adiabatic theorem or new physical
state is inferred from its formal order.

The evolved state remains `C_g=U_g C_up U_g†`. Any sufficiently regular
auxiliary operator P gives the exact decomposition

\[
D=C_g-P,\qquad R_\rho=P_\rho-[G_\rho,P],\qquad
D(\rho)=U(\rho,u)D(u)U(\rho,u)^\dagger
 -\int_u^\rho U(\rho,s)R_\rho(s)U(\rho,s)^\dagger ds.
\]

This add-and-subtract identity does not identify P with the source covariance
or add P to the action. If P is instead used to define a subtraction
functional, that functional's existing band-action identity remains a
separate requirement; this helper does not establish it.

## Exact kinetic and potential products

In the owned canonical frame,

\[
G_\rho=iH_\rho,\qquad H_\rho={k\over a^2}\sigma_3
 +{1\over a}(-m\sigma_1+\ell r_g^{-1}\sigma_2),
\]

because `dT/d rho=-1/a`. Thus `R_rho=P_rho-i[H_rho,P]`; in T coordinates,
`R_T=P_T+i[H_T,P]` and `R_rho=-R_T/a`.

Let `h=C*k+V(z)`, with C independent of z. On a numerical period L,
`V_l=L^{-1} integral_cell exp(-i omega_l z)V(z) dz` and
`omega_l=2*pi*l/L`. The Weyl kernel uses `exp(i*k*(z-z')) dk/(2*pi)`.
For Fourier coefficient p_n(k), direct operator multiplication gives

\[
[kC,p]_\#{}_n=k[C,p_n]+{\omega_n\over2}\{C,p_n\},
\]
\[
[V,p]_\#{}_n=\sum_l\left[
 V_l p_{n-l}(k-\omega_l/2)-p_{n-l}(k+\omega_l/2)V_l\right].
\]

There are no higher kinetic star terms. Potential shifts, rather than a
finite formal series, give the actual potential commutator. The formulas
are checked against independent Bloch-mode matrix multiplication, including
both boundary values of the Bloch parameter. A finite supplied Fourier
potential has an exact finite sum; a missing potential tail requires its own
bound. The numerical period is not a physical boundary condition.

## Quantitative Taylor remainder

Taylor truncation through M uses actual kth derivatives divided by `j!`.
For `s=omega_l/2`, the one-mode Frobenius remainder is bounded by

\[
 \|V_l\|\,{|s|^{M+1}\over(M+1)!}
 \left(\sup_{t\in[0,1]}\|\partial_k^{M+1}p(k-ts)\|_F
      +\sup_{t\in[0,1]}\|\partial_k^{M+1}p(k+ts)\|_F\right).
\]

This follows from the integral Taylor remainder. The two closed momentum
segments are required; a derivative at k alone is insufficient.
For omitted large-transfer modes, supply moments
`A_j=sum_omitted ||V_l||*|omega_l/2|^j`. If B_minus and B_plus bound the
shifted symbols there, and B_j bounds the jth derivative at k, triangle
inequalities give

\[
 \|\mathcal R_M\|_F\le A_0(B_-+B_+)
       +2\sum_{j=0}^M {A_j B_j\over j!}.
\]

The transfer split is numerical. It neither truncates the physical source
nor equates source energy with outgoing momentum. Both transfer signs must
be included in the supplied tail moments. Directed arithmetic retains
missing inputs as missing.

## Connection to the fourth-order recurrence

For `P=sum_{j=0}^J P_j`, the exact formal transport identities through J
imply, in T coordinates,

\[
iP_T-[h,P]_\#=i(P_J)_T-[kC,P_J]_{\#,1}
 -\sum_{j=0}^J \mathcal R_{J-j}(V,P_j),
\]

where `R_M` is the exact potential commutator minus its Taylor expansion
through M, including the zeroth term, and the kinetic subscript1 denotes its
first star term. Numerical errors in the formal identities must be added if
one uses numerical coefficient values. The expression is not obtained by
assigning a numeric error to the formal symbol `O(epsilon^5)`.

For J=4, a fourth axial derivative of the defect needs mixed geometry
derivatives through total order9: the time derivative of P4 adds one order,
and four axial derivatives add four more. The current SymbolJet owner trims
P4 to its value. It cannot supply those missing derivatives by relabeling.
The spatial profile extension, temporal derivatives, explicit momentum-decay
constants, nuclear-norm propagation and initial `C_up-P_up` bounds remain
separate inputs. No complete finite-history UV or physical-gate certificate
is produced here.
