# Nested qualities: a finite-window construction and its local response

## Result and scope

This note connects four existing ingredients into one explicit construction:
coupled regions, inherited law, normalized spectra and an effective local
response. The mathematics below is a finite-operator theorem, valid for every
finite depth under the stated hypotheses. It is not a certificate of the
separate gravitational incoming gate. Regional states need not be identical.

The construction reuses the [normalized recursion](nsc-regulated-recursion.md),
[local-observer correspondence](nsc-local-observer-correspondence.md) and
[energy-transfer identities](nsc-energy-transfer.md). Their historical owners
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
