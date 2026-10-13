# Nested qualities: a finite-window construction and its local response

## Result and scope

This note connects four existing ingredients into one explicit construction:
coupled regions, inherited law, normalized spectra and an effective local
response. The mathematics below is a finite-operator theorem, valid for every
finite depth under the stated hypotheses. It is not a certificate of the
separate gravitational incoming gate. Regional states need not be identical.

The construction reuses the [normalized recursion](nsc-regulated-recursion.md),
[local-observer correspondence](nsc-local-observer-correspondence.md) and
[energy-transfer identities](https://doi.org/10.5281/zenodo.23090015). Their historical owners
and records are unchanged. A successor is needed because the correspondence
v1 record hashes its original document; editing that document would invalidate
historical replay.

## 1. A concrete meaning of the inheritance transformation

Fix a finite dimension d, a Hermitian matrix H, a nonzero d-by-d link B,
and a positive scale ratio Omega. They form the shared dimensionless/model
data Theta. Work in orthonormal regional frames and use hbar=1 in evolution.
For an integer starting label m and depth N>=1, define the Hermitian operator
J_(m,N) on the direct sum of N copies of C^d by

$$
(J_{m,N})_{jj}=\Omega^{m+j}H,\qquad
(J_{m,N})_{j,j+1}=\Omega^{m+j}B,\qquad
(J_{m,N})_{j+1,j}=\Omega^{m+j}B^\dagger.
$$

Indices j run from zero; all other blocks vanish. Nonzero B makes the
regions coupled. A regional basis change is allowed if applied consistently
to the diagonal blocks, links and state.

Let U_m relabel the equal-sized window beginning at m as the one beginning at
m+1. It is unitary because it preserves the regional orthonormal frames. Then

$$
J_{m+1,N}=\Omega U_mJ_{m,N}U_m^\dagger,
\qquad
\mathcal T_\Theta^*J_{m+1,N}
:=\Omega^{-1}U_m^\dagger J_{m+1,N}U_m=J_{m,N}.
$$

**Proof.** Each diagonal and link block gains exactly one factor Omega;
unitary relabeling changes only its regional label. The equality follows
block by block. Iteration gives the same relation for every finite number
of label shifts. This realizes the organizing condition
`T_Theta^* mathbb_D_Theta = mathbb_D_Theta` as invariance of a family/law
under a specified normalized pullback. It is not the assertion that a
single fixed finite matrix is invariant under unscaled dilation.

For fixed N, spectra obey

$$
\sigma(J_{m,N}/\Omega^m)=\sigma(J_{0,N}).
$$

The isolated regional block likewise has eigenvalues
`omega_n,k = Lambda_n * hat_omega_k(Theta)` in the chosen frequency units.
The eigenvectors map by U_m. Thus dimensionless spectral relationships
are inherited without fixing a numerical value of Omega from observation.
Changing the number N of coupled regions changes the boundary problem and
may change the dressed spectrum; fixed-window covariance must not be
mistaken for independence from the surrounding system.

## 2. Nested reduction and composition

Choose z with positive imaginary part. J_(m,N) and all its principal
subblocks are Hermitian, so the required resolvents exist. Eliminate the
last block and proceed inward using the standard Schur complement. If
S_j(z) is the inverse response after eliminating levels deeper than j,

$$
S_{N-1}(z)=zI-\Omega^{m+N-1}H,
$$
$$
S_j(z)=zI-\Omega^{m+j}H
-\Omega^{2(m+j)}B S_{j+1}(z)^{-1}B^\dagger.
$$

**Proof.** Apply the two-block inverse identity once to the terminal
region, then induct on the remaining depth. Consequently
`S_0(z)^(-1)` is exactly the first regional block of `(zI-J_(m,N))^(-1)`.
Eliminating an already-grouped exterior gives the same result. No
approximation, change of local law or fitted external force is used.

For three blocks K=[[A,B,0],[B^dagger,C,E],[0,E^dagger,F]], the result is
`K_local=A-B(C-E F^(-1) E^dagger)^(-1) B^dagger` whenever the indicated
inverses exist. This displays the influence of a region's own surroundings
inside the response seen one level further in.

Normalizing S_j by its regional scale and writing x=z/Omega^(m+j) gives

$$
\Gamma_j(x)=xI-H-\Omega^{-1}B\Gamma_{j+1}(x/\Omega)^{-1}B^\dagger.
$$

This is the already-owned recurrence, including its necessary 1/Omega.
It proves coherent nesting at arbitrary finite depth. An infinite operator
or geometric spacetime would require its own stated domain/limit; neither
is a requirement or a result of this finite theorem.

## 3. Shared law, distinct regional states

Choose any Hermitian one-particle fermionic covariance C with 0<=C<=I.
It can have unequal diagonal regional blocks and nonzero cross-region
correlations. The common finite Hamiltonian evolves it by

$$
C(t)=e^{-iJt}C(0)e^{iJt}.
$$

Unitary evolution preserves admissibility. It imposes no equality between
regional occupations. For example C_stat=I/2 is stationary for the same J
for which a nonuniform diagonal C_var can have `[J,C_var] != 0`.
This establishes a distinction between the inherited operator law and the
chosen state; it does not identify either example with stellar collapse.

Split the full operator into a retained region A and its complement E,
with coupling V. Eliminating the exterior from the same evolution gives

$$
i\dot\psi_A=J_{AA}\psi_A
-i\int_0^t V e^{-iJ_{EE}(t-s)}V^\dagger\psi_A(s)\,ds
+V e^{-iJ_{EE}t}\psi_E(0).
$$

This follows by solving the exterior linear equation with its initial
condition and substituting it into the retained equation. The retarded
operator depends on J and V; the initial-noise covariance depends on C_EE,
and cross correlations on C_AE. Equal laws therefore permit different
local histories and responses through different regional state data.

For static J the full energy obeys `d Tr(CJ)/dt=0` by cyclicity of trace.
Regional energy can still be exchanged. This is the existing finite
energy-transfer result, not a new source inventory.

## 4. What has been established

| Nested quality | Result |
|---|---|
| Coupled regional hierarchy | Explicit nonzero-link Hermitian finite windows and associative elimination |
| Inherited common law | Exact normalized window-shift covariance with a specified transformation |
| Normalized spectrum | Unitary scale covariance for homologous fixed-depth windows and regional blocks |
| Exterior influence retained locally | Exact Schur response, retarded memory and initial-state/correlation terms |

The theorem supplies a compatible mathematical realization of nested
qualities using standard finite quantum dynamics. Schur algebra and unitary
spectral covariance are imported mathematics, not new physical laws.
It does not determine c or Lambda, derive cosmological abundance fractions,
identify charge conjugation with spatial exchange or prove an eternal engine.
Those are not acceptance conditions for this result. A claim about a
particular curved Dirac realization still inherits that realization's domain
and source assumptions; no old local-gate verdict is changed.

## 5. Bounded verification

`scripts/check_nsc_nested_qualities.py --check` authenticates and replays the
small successor record `results/development/nsc-nested-qualities-v1.json`.
It checks exact rational complex matrices at three coupled levels, shifted
window covariance, the nested/direct resolvent, the 1/Omega recurrence and
two different admissible states. Mutation controls remove the scale factor,
reverse a noncommuting coupling or omit state information. These are
implementation controls of the displayed proof, not cosmological simulations.

The supporting v1 correspondence and normalized-recursion records are reused
by hash. No older scientific record or its hashed source document is edited.

## 6. The massless reference and its finite-cutoff allocation

On the child tail, $u=-1/\rho$ gives the smooth conformal metric

$$
\bar g=g/r^2=\frac{du^2}{(1+u^2)^2W(u)}-W(u)dz^2-d\Omega_2^2,
\qquad W(u)=3(\pi-\arctan u)+\frac{3u-u^2}{1+u^2}.
$$

The auxiliary continuation agrees with this metric above $u=-0.05$, joins
smoothly to $W=3\pi$ on $[-0.1,-0.05]$, and stays constant below $-0.1$.
Its ultrastatic negative-energy projector propagates by the conformal Dirac
operator, with $\psi=r^{-3/2}\bar\psi$. Imported Hadamard propagation and
the anomalous conformal stress law supply this massless reference; they do
not select the physical parent state. The unwrapped angular spectrum has
$\omega^2=k^2/(3\pi)+\kappa^2$, $|\kappa|\geq1$.

In the retained proper-time convention,
$T^\mu{}_\mu=(-C^2/20+11E_4/360-\Box R/30)/(16\pi^2)$. The asymptotic
reference has $H_{\rm child}=\sqrt{3\pi}$ and

$$
\rho_{\rm ref}=\frac{33}{320}-\frac{3u^2}{32}+O(u^3),\qquad
p_{\parallel,\rm ref}=p_{\perp,\rm ref}=-\frac{33}{320}+\frac{u^2}{32}+O(u^3),
\qquad \rho_{\rm ref}+p_{\parallel,\rm ref}=-\frac{u^2}{16}+O(u^3).
$$

This tail expansion cannot be applied at the neck. Its finite $R^2$
ambiguity must be matched by the complement. Define
$\mathcal C_{\nu,\mu}=\Gamma_{\rm light,\nu}-\Gamma_{\rm light,ren}[\mu]$;
the exact allocation is
$\Gamma_5=\Gamma_{\rm light,ren}+(\Gamma_{H,\nu}+\mathcal C_{\nu,\mu})$.
At curvature radius $a=0.325735$ and $\nu=1$, the reference projection is
$0.103125000$, the original light projection $2.8628\times10^{-16}$, and
the converted complement $-0.102079088$. Their unchanged total is
$0.00104591235$. Reallocation does not create a new source.

## 7. Normalization, state conversion and the Weyl cocycle

For $x=d^2/\Lambda^2$, the exact profile comparison is

$$
g_{\rm PT}=\tfrac12E_1(x),\qquad
g_{\rm log}=-\tfrac12e^{-x}\log(d^2/M^2),\qquad
h=-\tfrac12[e^{-x}(\log x+\gamma_E)+E_1(x)],
$$
$$
g_{\rm log}=g_{\rm PT}+[\log(M/\Lambda)+\gamma_E/2]e^{-x}+h.
$$

$M$ is a normalization mass. The remainder has $h(0^+)=0$,
$\int_0^\infty h\,dx=-1/2$ and $\int_0^\infty xh\,dx=-3/4$.
A heat-weighted rank and a sharp spectral projector therefore differ in
their finite volume and Einstein terms; anomaly agreement does not choose
the measure. No profile is fitted to hold a throat open.

For the existing remainder, $dE_{\rm sub}/d\log M=A_4/(16\pi^2)$ with
$a_4=(-18C^2+11E_4-12\Box R)/360$. At fixed reference mass, changing
$M\to sM$ is canceled by $\Delta c_i=\beta_i\log s$, where the basis is
$(M_{\rm ref}^4,M_{\rm ref}^2R,C^2,R^2,E_4,\Box R)$ and

$$
\beta_i=\left(0,0,\frac1{320\pi^2},0,-\frac{11}{5760\pi^2},
\frac1{480\pi^2}\right).
$$

The continuum closed-cell bulk shift is $C^2$ alone. The discrete Euler
column is kept for exact finite-code cancellation; at $s=2$ its leftover
falls from about $1.6\times10^{-10}$ at 32 points to below $10^{-12}$ at
64 and 128. This subtraction ratio $s$ is distinct from recursive
$\Omega$, and it supplies no physical Newton or cosmological running.

Finite proper-time conversion also depends on the state. For one massless
Dirac field on the flat thermal circle,

$$
\Delta F_{\mathcal C}=-\frac{4T^4}{\pi^2}
\sum_{n\geq1}\frac{(-1)^n}{n^4}(1+x_n)e^{-x_n},
\qquad x_n=\frac{n^2\nu^2}{4T^2}.
$$

At $T/\nu=0.5$, canonical thermal energy is $0.0719658654\nu^4$,
while the finite-endpoint thermal determinant gives $-0.0012266736\nu^4$.
The unit pole residue does not equate these state prescriptions. The common
CTP conversion must retain the canonical state, complement and contour
correction together; a state-independent local vacuum term cannot remove
the thermal difference.

Under common metric Weyl scaling,
$D_\omega[g_\phi]=e^{-\phi/2}D_\omega[g]e^{-\phi/2}$. The same modulus gives

$$
\mathcal A_\Lambda[g,\phi]=E_\Lambda[g_\phi]-E_\Lambda[g]
=\int_0^1dt\,\operatorname{Tr}[\phi e^{-D[g_t]^2/\Lambda^2}],
\qquad g_t=e^{2t\phi}g.
$$

Independent heat-path and determinant evaluations agree. The local volume
coefficient is $1/4$, whereas the sharp-projector coefficient is $1/8$.
For $\phi=0.05+0.04\cos(2\pi x/L)$, the retained local anomaly integral is
approximately $-0.002804148$. The cocycle composes as
$\mathcal A[g,\phi+\sigma]=\mathcal A[g,\phi]+\mathcal A[g_\phi,\sigma]$.
It supplies the compensating integrand once; it does not fix the complete
measure or select the physical relative-sheet scale.

## 8. The charged throat dictionary and relational vacuum law

For $\int\sqrt{-g}(AR-CF^2-V)$, imported MMP normalization gives
$G_N=1/(16\pi A)$, $g^2=1/(4C)$, $r_e^2=q^2C/(4A)$ and
$\ell=16r_e^3/(|q|G_N)$. Integer magnetic flux and a closed antiperiodic
return path are part of that domain. At the first matching point $\nu=1$,
$\Xi_{\rm partial}/q^2=V_DC_D/(8A_D^2)=3.4459938566757384$ exceeds the
charged seed bound $1/4$. The full $V_{\rm full}$ is the missing
common-action input in that partial embedding.

The relational normalization route imposes, after the full assembly,

$$
\Gamma_{\rm rel}=(1-\mathcal P_0)\Gamma_{\rm one},\qquad
\mathcal P_0(V,A,C)=(V,0,0),\qquad
\left.\frac{\delta\Gamma_{\rm rel}}{\delta g_{\mu\nu}}
\right|_{\rm homogeneous,unlinked}=0.
$$

It removes the local $a_0$ tadpole while preserving Einstein, gauge,
curvature, finite Casimir and link response. Under inherited units,
$V\to\Omega^4V$ and $\mathrm{Vol}\to\Omega^{-4}\mathrm{Vol}$ keep its
CTP subtraction invariant. This route sets $V_{\rm full}^{\rm rel}=0$,
$\lambda_4=0$ and $\Xi=0$; the physical return state and complete induced
coefficients remain separate conditions.

Keeping the charged light field explicit at
$\mu_B=\sqrt{|q|}/r_e$ gives
$A=[Q_1[h_\Lambda]-\mu_B^2]/[6(4\pi)^2]$ and
$C=H_{\Lambda,\mu_B}(0)/[3(4\pi)^2]$. With $r_e=L_*=1$,
$\Omega=\Lambda r_e$ and $\zeta=\Omega^2$, the radius condition is

$$
F_q(\Omega)=\frac{q^2C(\Omega,\sqrt{|q|})}{4A(\Omega,\sqrt{|q|})}-1=0.
$$

The first checked cutoff-resolved branch is $q=4$:
$\Omega=3.973074368754324$, $\zeta=15.78531993965257$,
$AL_*^2=0.04501936182826088$, $C=0.011254840457065235$ and
$V_{\rm full}L_*^4=0$. Its radius residual is $1.33\times10^{-15}$,
$H_{\rm child}/\Lambda=0.7726963653$, and 26-to-34 compact modes change
$\Omega$ by $8.64\times10^{-9}$. The $q=2,3$ roots fail the same cutoff
condition. This finite candidate does not select flux, establish the
large-$q$ regime, or supply the physical charged boundary state.
Original merged records are at
[Zenodo](https://doi.org/10.5281/zenodo.23090015).
