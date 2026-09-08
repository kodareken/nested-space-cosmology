# Evaluated full-spinor throat response and the chiral boundary obstruction

The smooth massless spatial throat has a nonzero energy-dependent parent/child
boundary link, but that link preserves physical chirality in a compatible
domain. Its common-frame action kernel has vector and axial Clifford content;
the scalar, pseudoscalar and tensor coefficients vanish to numerical precision.
It therefore does not supply the scalar sheet mass proposed in the
[kinematic observable bridge](nsc-observable-bridge.md). This identifies the
specific local realization that needs an additional, action-derived interaction
or domain mechanism. It does not reject the recursive architecture.

The calculation is stored as `NSC-8-CHIRAL-BOUNDARY-v1` in
[`results/nsc-8-chiral-boundary.json`](../results/nsc-8-chiral-boundary.json).
The [new module](../src/recursive_horizons/nsc_chiral_boundary.py) evaluates the
full four-component paired-angular Green function and reduced port response.
It preserves the historical staggered calculation as a separate regulator.
The PG current and transmission domain are derived below; the numerical maps
in this record are finite **spatial** maps, not the full Lorentzian problem.

## The spinor and the domain

Use the same angular convention as the observable bridge,

\[
\psi_\kappa=r^{-1}(F\Omega_\kappa,iG\Omega_{-\kappa})^T,
\qquad \sigma\cdot\hat r\,\Omega_\kappa=-\Omega_{-\kappa}.
\]

Physical chirality sends `(F,G,kappa)` to `(iG,-iF,-kappa)`. In tensor order
paired angular sign `eta`, then radial component `rho`, this is

\[
\Gamma_5=-\eta_1\otimes\rho_2,\quad
\alpha=I\otimes\rho_2,\quad
\beta_D=I\otimes\rho_3,\quad
V(x)=\frac{\kappa}{\sqrt{1+x^2}}\eta_3\otimes\rho_1.
\]

Here `eta` is the paired angular label, not a geometric sheet. The flat radial
measure follows the original rescaling. In units `hbar=c=L_star=1`,

\[
H_0=-i\alpha\partial_x+V(x),\quad
[\alpha,\Gamma_5]=[V,\Gamma_5]=0,\quad
\{\beta_D,\Gamma_5\}=0.
\]

The Green identity is

\[
\langle\psi,H_0\chi\rangle-\langle H_0\psi,\chi\rangle
=-i[\psi^\dagger\alpha\chi]_{-R}^{R}.
\]

An operator symmetry requires an invariant domain. The finite spatial
benchmark uses the two-dimensional `Q=+1` subspace at each exterior endpoint,
where `Q=eta3 tensor rho3`. Since `{Q,alpha}=0`, these are maximal isotropic
wall spaces for the spatial current; since `[Q,Gamma5]=0`, chirality preserves
the domain. This wall is an explicit finite-box control choice, diagonal in
the paired angular modes. It is not derived as a cosmological exterior, and
there is no assertion of a local angular wall law before restoring all modes.
At the throat the spinor is continuous in the common coordinate frame.

The historical condition `F=0` at both exterior endpoints instead selects
`beta_D=-1`. It is a valid reflecting spatial wall but is not invariant under
physical chirality. Applying the same formal bulk commutator to that domain
would omit its boundary breaking. The experiment evaluates it separately as
an artificial wall control.

## The PG current requires transmission

The [tetrad derivation](nsc-dirac-tetrad.md) and
[Lorentzian propagation calculation](nsc-lorentzian-transport.md) give, after
pairing angular modes,

\[
H_{\rm PG}=-i A(x)\partial_x+\frac{i}{2}\beta_s'(x)I+V(x),
\qquad A=I_\eta\otimes(\rho_2-\beta_s I).
\]

The shift `beta_s` and the Dirac mass matrix `beta_D` are distinct. Since
`A'=-beta_s' I`, the derivative term is `-i(A d+A'/2)`. Its Green identity
contains `-i psi^dagger A chi`, so the number current is

\[
j=\psi^\dagger[I_\eta\otimes(\rho_2-\beta_s I)]\psi.
\]

This current, its characteristic projectors, and the zero incoming-data
domain all commute with physical chirality. Parent is `x>0`, child is `x<0`.
The parent outward throat normal is `-d_x`; the child outward throat normal
is `+d_x`. With identity transmission and a common trace, their oriented
boundary forms cancel. Outward-normal sign is not an extra spinor reflection.

At the throat `beta_s=sqrt(3 pi/2)=2.17080376367`. Both characteristic speeds,
`1-beta_s=-1.17080376367` and `-1-beta_s=-3.17080376367`, point from parent
toward child. There are four incoming child-side channels and four outgoing
parent-side channels when the angular pair is included. At the truncation
`x=-4` every channel is outflow; at `x=4` two channels are incoming and their
data vanish. The horizon is `x=1.90069160547`.

The spatial reflecting wall is invalid as a PG outflow boundary: its allowed
space has restricted current `-beta_s I2`, which is nonzero. At a trapped
cut the full current is definite, so no nonzero reflecting isotropic space
exists. The present numerical spatial map cannot be promoted to a PG map by
inserting a shift into a Hermitian reflecting matrix. A stationary transfer
equation that inverts `A` also meets the characteristic degeneracy at the
horizon. This record does not make that continuation.

## Why elimination cannot manufacture this mass

Let `K=H-z`, and suppose the closed operator domain is invariant under
`Gamma5`. In a compatible finite representation let the retained and
eliminated projections commute with chirality. Every block then intertwines
its source and target chirality operators. If the eliminated block is
invertible,

\[
[K_{ee}^{-1},\Gamma_{5,e}]=0,\qquad
S=K_{rr}-K_{re}K_{ee}^{-1}K_{er},\qquad
[S,\Gamma_{5,r}]=0.
\]

The same conclusion follows for traces of the resolvent when trace
evaluation and source insertion intertwine chirality. For a first-order
continuum operator, a sharp room cut also requires trace and transmission
data; writing a block matrix alone does not settle its domain. Here both the
domain and these data are explicit. The numerical Schur step is an exact
finite algebra operation on the evaluated port kernel.

In a common parent/child frame the nonzero off-diagonal block `B` therefore
obeys `[B,Gamma5]=0`. If it also obeyed `{B,Gamma5}=0`, their sum would give
`2 B Gamma5=0`, hence `B=0`. This is incompatible with the nonzero link.
The kinematic projector

\[
\Pi_{\rm sheet}=\frac{I-\tau_3\otimes\Gamma_5}{2}
\]

requires chirality-preserving diagonal blocks and a chirality-odd link to
define its invariant parent-left/child-right sector. The computed port
response fails that restriction. This is a concrete test of the proposed
identification on the boundary traces, not a derivation of a complete
two-room Hamiltonian or of a particle pole.

## Controlled regulator and the actual map

Every endpoint of every cell carries all four spinor components. On the
cell `[x_j,x_{j+1}]`, let `h=x_{j+1}-x_j` and

\[
F_j=i\alpha[zI-V(x_{j+1/2})],\qquad
T_j=(I-hF_j/2)^{-1}(I+hF_j/2).
\]

The midpoint equation is `psi_{j+1}=T_j psi_j`. It exactly respects
`[T_j,Gamma5]=0`. At real energy it satisfies `T_j^dagger alpha T_j=alpha`.
At complex energy, with `M_j=(I+T_j)/2`, it satisfies

\[
T_j^\dagger\alpha T_j-\alpha
=-2\,\operatorname{Im}(z)\,h\,M_j^\dagger M_j.
\]

This is a mimetic first-order transfer regulator with an exact cell current
identity. It is a new controlled representation, not the old staggered
node/edge matrix with a four-by-four matrix naively attached. Its finite
response is validated against the continuum. No spectral doubling theorem
or ultraviolet spectral completeness is claimed for it.
The assembled Green map also satisfies conjugate-energy Schwarz reflection
and has positive imaginary part in the upper energy half-plane, checking the
source-jump orientation and the spatial resolvent convention independently.

Point sources impose the exact integrated jump

\[
G(y^+,y)-G(y^-,y)=i\alpha.
\]

At coincident ports the Green trace is the arithmetic mean of the two sides.
This convention is part of the kernel definition, including the comparison
with the continuum. Unit-integrated radial point sources make this Green
kernel dimensionless. Its inverse is a port response; it must not be assigned
the units or locality of a bulk mass coefficient without further trace
normalization and localization data.

The primary sparse LU solve simultaneously imposes all midpoint cell
equations, all source jumps and both exterior walls. It evaluates the 12 by
12 Green matrix at ports `(parent x=1, throat x=0, child x=-1)`. Inverting this
small matrix gives the three-port kernel. Schur elimination of the throat
port gives the full eight by eight parent/child response `S(z)`. A separate
joined solve with only the two retained source ports gives the same result.

Two independent comparisons are evaluated: multiplication of four by four
cell transfer matrices with the explicit endpoint wall solve, and DOP853
integration of the continuum first-order equation with the same domain.
The record contains the actual complex Green matrices, full responses and
link entries at each energy. All are recomputed by `--check`.

## Clifford content, values and controls

Choose an oriented algebraic Clifford completion of the physical radial
`alpha`, `beta_D` and `Gamma5`:

\[
\alpha_1=I\rho_2,\quad \alpha_2=\eta_3\rho_1,\quad
\alpha_3=\eta_2\rho_1,\quad
\gamma^0=\beta_D,\quad \gamma^a=\beta_D\alpha_a.
\]

It has signature `+---` and `i gamma0 gamma1 gamma2 gamma3=Gamma5`. The
two additional axes complete the matrix algebra; their identification with
a local angular tetrad has not been separately derived. The vanishing of
the chirality-even action sector is independent of this completion. The
individual vector/axial component labels below refer to this chosen basis.
The
action-level port kernel is `K_action=beta_D B`; it is decomposed in the 16
orthogonal matrices `I, Gamma5, gamma_mu, gamma_mu Gamma5, sigma_mu_nu` using
`c_C=Tr(C^dagger K_action)/4`. A partial-wave nonlocal kernel is not yet a
local Lorentz-covariant interaction merely because this decomposition exists.

For the baseline `[B,Gamma5]=0` and `{beta_D,Gamma5}=0` imply
`{K_action,Gamma5}=0`. Scalar, pseudoscalar and tensor components commute with
chirality and therefore vanish; vector and axial components are allowed.
At `R=4`, `kappa=±1` and 160 cells (`h=.05`), the evaluated results are:

| Complex energy `z` | Link Frobenius norm | Relative continuum response error | `||[S,Pi_sheet]|| / ||S||` |
|---|---:|---:|---:|
| `0.2+0.25i` | 0.974498 | 0.000094703 | 0.224717 |
| `0.6+0.35i` | 1.395915 | 0.000110174 | 0.249895 |
| `1.1+0.30i` | 2.453605 | 0.000364517 | 0.383806 |

The full response chiral commutator is below `4e-15` relative, and direct
versus Schur error is below `2e-15`. The three-port Green condition numbers
are approximately `4.51`, `9.67` and `14.60`, so these probes avoid an
ill-conditioned inversion. At the first energy the nonzero action
coefficients are approximately

| Clifford component | Complex coefficient |
|---|---:|
| `gamma0` | `0.0205425+0.184406i` |
| `gamma1` | `-0.0388172+0.337328i` |
| `gamma2` | `-0.283364-0.0636012i` |
| `gamma3 Gamma5` | `-0.0538344+0.0211764i` |

These are evaluated energy-dependent finite-domain coefficients. None was
fit to a desired mass or observation. Scalar, pseudoscalar and tensor
coefficients are at rounding level. The 40, 80, 160 and 320-cell family at
the first energy has response convergence orders `2.0023`, `2.0006`,
`2.0001`. The one-percent gate is a loose numerical agreement bound, not a
claimed physical accuracy of the cosmological model.

Three controls expose distinct false positives:

1. **Inserted scalar mass.** The local link `B=phi beta_D`, with prescribed
   `phi=.35`, exactly preserves the candidate `Pi_sheet` and gives the usual
   squared Dirac dispersion. It is the kinematic control. Separately adding
   `m=.35 beta_D` to the spatial bulk operator produces a boundary link with
   chiral-odd norm approximately `0.328475`. This verifies sensitivity to
   a scalar interaction; neither inserted parameter is throat-derived.
2. **Domain breaking.** Retaining the massless differential expression and
   switching the exterior walls to historical `F=0` produces a chiral-odd
   link of norm approximately `0.115699`. That signal is an imposed wall
   effect and cannot be called a generated interior mass.
3. **Frame transport.** The normal-reflection lift
   `U=eta1 tensor rho3` gives `Gamma5_child=U Gamma5 U^dagger=-Gamma5`.
   A mixed-frame link `B_local=B_common U^dagger` then appears to
   anticommute with the *untransported* matrix. The correct equation is
   `Gamma5_parent B_local=B_local Gamma5_child`. Restoring the common frame
   restores the nonzero anticommutator norm ratio `2`. Changing coordinates
   or spin frame has not generated a scalar interaction or selected a
   chiral sheet sector. This is passive frame transport; the wall operator
   must transform too, as `U Q U^dagger=-Q`. The same `Q=+1` wall choice at
   both outer endpoints is not invariant under active reflection. This
   control makes no finite-wall parity claim.

## Reproduction and next intervention

```sh
python3 -m unittest tests.test_nsc_chiral_boundary -v
python3 scripts/check_nsc_chiral_boundary.py --check
```

`--check` recomputes every numerical array, convention, flag and source hash.
`--output PATH` opens a new file exclusively and refuses to overwrite an
existing record. Sources include this document and the numerical test and
runner, together with the imported convention owners. Historical result
bytes are not changed.

The decisive result is an evaluated massless boundary response with a
domain-aware obstruction to the scalar identification. The next constructive
intervention must identify a chirality-odd term or physical boundary/state
mechanism from the common action, then transport frames and rerun the same
test. Choosing an arbitrary coefficient or an incompatible wall would only
repeat the controls. The physical PG boundary map, local particle pole,
sector selection, stationary common-action stress and observational
prediction remain separate unresolved dependencies.
