# The compact carrier fixes a matrix of interaction overlaps

Applying the published torsion contact to the [actual free compact
modes](nsc-compact-mass-map.md) produces nonzero couplings between levels.
A single massive level is not an exact interacting truncation of this
candidate bulk action. The geometric warp cancels from free propagation,
but remains in the interaction after canonical field normalization. The
calculated overlaps specify this dependence without choosing a state,
an absolute gravitational strength, or an extra scalar field.

This is a development application of an imported interaction. It does not
adopt a new independently weighted gravitational action. The five-dimensional
torsion coefficient and derivative/boundary completion must still be matched
to the common NSC functional.

## Reuse the published contact with the same conventions

[Castillo-Felisola et al., equations13 and18–20](https://arxiv.org/html/1405.0397v1)
provide the higher-dimensional Einstein–Cartan contact, its Clifford
decomposition, and a reduction with chosen zero-mode profiles. We reuse
the contact and substitute the NSC profiles rather than importing their
Randall–Sundrum geometry or zero-mode approximation.

Use signature +----, Gamma^Y=i gamma5 and bar(Psi)=Psi-dagger gamma0.
Relative to that paper's -++++ convention, gamma_paper=i Gamma and
bar(Psi)_paper=-bar(Psi). Their contracted triple-current square keeps its
sign under this conversion. Write the local candidate as

\[
S_{4f}=\frac{\kappa_5^2}{32s_T}
\int d^5x\sqrt{|g_5|}\,
(\bar\Psi\Gamma^{ABC}\Psi)(\bar\Psi\Gamma_{ABC}\Psi).
\]

Here kappa5 squared=8 pi G5 denotes the coefficient of the already matched
two-derivative bulk Einstein term. The ratio s_T is the five-dimensional
torsion stiffness relative to Einstein–Cartan; s_T=1 is the imported EC
case. **Its NSC value is left symbolic.** The previous four-dimensional
factor-nine comparison does not by itself fix this five-dimensional
functional, its finite terms or its derivative contributions.

The triple-current contraction separates into 6 times an axial-vector
contraction and 3 times an ordered axial-tensor contraction, with
gamma_ab=gamma_[a gamma_b] and both a,b orders included. In four dimensions
the axial coefficient at s_T=1 is 3 kappa4 squared/16=3 pi G4/2, matching
[Lucat–Prokopec's equation9](https://arxiv.org/html/1512.06074v1) at xi=1.
These are action-coefficient conventions; they do not establish a stress sign.

## Normalize the interaction and gravity together

For g5=exp(2 sigma)(g4-dY²) and chi=exp(2 sigma)Psi, four spinor factors
and the volume element leave exp(-3 sigma) in the contact. The coefficient
of R4 in the two-derivative bulk gravitational action has weight exp(3 sigma):

\[
I_3=\int_{-L_\star}^{L_\star}e^{3\sigma}\,dY,
\qquad
\kappa_{4,\mathrm{bulk}}^2=\frac{\kappa_5^2}{I_3}.
\]

This identifies the bulk contribution, not the total measured G4 when
boundary curvature terms or further quantum matching are present. The
compact size and warp are held fixed in this reduction. A varying compact
size requires its additional metric terms.

Let F_n=f_L,n P_L+f_R,n P_R be the normalized mode embedding and
bar(F_n)=beta F_n-dagger beta. The complete reduced vertex is

\[
V_{ij;kl}^{ab;cd}=I_3\int dY\,e^{-3\sigma}
 (\bar F_i\Gamma^{ABC}F_j)_{ab}
 (\bar F_k\Gamma_{ABC}F_l)_{cd},
\]

\[
S_{4f}=\frac{\kappa_{4,\mathrm{bulk}}^2}{32s_T}
\int d^4x\sqrt{|g_4|}\sum_{ijkl}
V_{ij;kl}^{ab;cd}\bar\psi_{i,a}\psi_{j,b}\bar\psi_{k,c}\psi_{l,d}.
\]

The mode sum retains the Weyl zero mode and all massive Dirac levels. This
formula is a bulk contact approximation, requiring the neglected torsion
derivative term to be controlled. It is not an integration over a dynamical
torsion propagator or an already renormalized in-in source.

## Evaluated coefficients on the NSC warp

Define dimensionless profile overlaps

\[
W_{ijkl}^{stuv}=I_3\int dY\,e^{-3\sigma}
f_{i,s}f_{j,t}f_{k,u}f_{l,v}.
\]

For sigma=-18(Y/L_star)²/1015, I3/L_star=1.96509101243726. The common
factor kappa4,bulk squared/(32 s_T) and the Clifford factors above are not
included in these W values:

| Modes and chiralities | Unwarped control | NSC warp |
|---|---:|---:|
| 0000, LLLL | 1 | 1.00025162265549 |
| 1111, LLLL | 3/2 | 1.52339608755034 |
| 1111, RRRR | 3/2 | 1.48012862049561 |
| 1111, LLRR | 1/2 | 0.49874089128800 |
| 0011, LLLL | 1 | 1.01106848941917 |
| 0011, LLRR | 1 | 0.98943475589180 |
| 3111, LLLL | 1/2 | 0.52112422157451 |
| 3111, RRRR | -1/2 | -0.48743692131945 |

These ratios are fixed by the declared geometry and domain. They are not
particle masses, decay rates, observed coupling constants or new adjustable
parameters. The zero-mode ratio is I3 I_minus3/(4 L_star²), which is at least
1 by Cauchy–Schwarz; it equals 1 for a constant warp.

## Why neither the other levels nor the tensor product can be dropped

The implementation evaluates the full spinor vertex and antisymmetrizes its
barred and unbarred fermion legs. The nonzero 3111 coupling remains after
this projection: the left-current component with Weyl matrix indices
(a,b,c,d)=(2,0,3,1) is -6.25349065889412 in the stated vertex normalization.
The corresponding unwarped value is -6. Thus a generic level-1 field sources
the equation for level3 through the interaction. A small value of the overall
coupling can justify an approximation; it does not make a one-level sector
exactly invariant. These component signs are not energy-density signs.

There are genuine selection rules. The profiles obey
f_L,n(-Y)=(-1)^n f_L,n(Y) and
f_R,n(-Y)=(-1)^(n+1) f_R,n(Y). With this even bulk warp, vertices with
i+j+k+l odd vanish. The 0111 vertex vanishes, whereas 3111, 0011 and 0101
do not. Boundary interactions may change this symmetry. Nonzero vertices
permit processes; actual rates additionally require states and kinematics.

For the massive level1, f_L,1 f_R,1 is odd, but its square is even and
W1111,LLRR is nonzero. The tensor part of the antisymmetrized 1111 vertex
also remains nonzero. Dropping an odd *mean bilinear* therefore does not
justify dropping its contact product. The pure Weyl zero mode, in contrast,
has no tensor bilinear. This is why the imported zero-mode truncation does
not establish the same simplification for the massive tower. Fierz
rearrangements can change the displayed channel basis; the full fermionic
vertex, rather than the name of a channel, is the retained interaction.

## Source interface and reproduction

The reusable microscopic source must carry a mode-matrix propagator
S_nm(x,x') and this vertex in its two Wick contractions. A scalar mass and a
single axial coupling do not specify that interacting tower. The remaining
physical inputs are its common-functional stiffness/derivative matching,
compact boundary action, physical scale and initial state. Absolute stress,
energy transfer and self-sourced geometry must then follow from that same
functional; none is inferred from the overlaps alone.

The [development record](../results/development/compact-interaction.json)
authenticates its code and mass-map input. The
[runner](../scripts/check_nsc_compact_interaction.py) checks every field with
atol=rtol=3e-12, with exact source hashes, keys and types. It compares direct
five-dimensional triple matrices against the axial/tensor decomposition,
checks exact unwarped integrals and uses independent adaptive integration.
The largest fixed-versus-adaptive overlap difference is below 2e-14; changes
from 32 through 96 quadrature points are also recorded. These are numerical
checks of the finite overlap calculation, not a continuum quantum remainder
bound or a self-sourcing test.

Run `python3 scripts/check_nsc_compact_interaction.py --check`.
The record is subsequent development evidence; the v0.3.0 paper and its
81-record released collection remain immutable.
