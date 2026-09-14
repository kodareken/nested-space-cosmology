# Step U: the complete seam traces do not identify a Hamiltonian link

This calculation uses the authenticated T artifact and the existing seed
Hamiltonian. It tests the proposed conversion into the independent canonical
parent/child blocks of `CausalCommonFunctional`. No B is supplied or adopted.

## The literal independent-trace embedding fails

T gives two oriented boundary traces of one field. After applying the stored
trace inverse, use the most favorable normalized redundant embedding

$$
E=\frac1{\sqrt2}\begin{pmatrix}I\\I\end{pmatrix},\qquad b=Ea.
$$

The imported CAR identity then gives

$$
\{b,b^\dagger\}=EE^\dagger,
\qquad \mathrm{rank}(EE^\dagger)=d<2d,
\qquad \|EE^\dagger-I_{2d}\|_2\ge1.
$$

For the locked retained channels, $d=3808$ and the proposed two complete
trace blocks have 7616 slots. The record evaluates this operator-CAR mismatch
on every stored block, with target tolerance $3\times10^{-11}$.

This statement is about the **operator anticommutator**, not the occupation
covariance. The matrix $ECE^\dagger$ can pass the existing covariance
positivity check while the redundant variables still fail to be independent
canonical room variables. Treating them as independent would add the missing
canonical modes, even if those modes were initially empty.

This rejects the literal use of T's two complete traces as the two independent
Hamiltonian rooms. It does not reject an orthogonal spatial partition of the
bulk Dirac Hilbert space, and it is not a non-existence proof for the full
extended theory.

## A constrained redundant lift does not fix B

If the redundant pair is retained only as a representation of the physical
graph, define its orthogonal complement by

$$
F=\frac1{\sqrt2}\begin{pmatrix}I\\-I\end{pmatrix}.
$$

Every Hermitian extension in the following family has the same compression
and evolution on the physical graph:

$$
H_J=EhE^\dagger+FJF^\dagger,
\qquad H_JE=Eh,
\qquad E^\dagger H_JE=h.
$$

In these coordinates the blocks are

$$
H_{pp}=H_{cc}=\frac{h+J}{2},\qquad B=\frac{h-J}{2}.
$$

The two algebraic witnesses $J=0$ and $J=h$ therefore give $B=h/2$ and
$B=0$, respectively, while preserving the same physical $h$. The record
checks this using the actual locked mode Hamiltonians and their owned four
metric vertices. Their B derivatives also differ while the compressed
physical derivatives agree.

The fixed object is the physical single-field $h$. The redundant diagonal
blocks also change between these two witnesses. T has not supplied separately
fixed, independent parent and child bulk Hamiltonians, so those are not
additional conditions imposed by this certificate.

Neither witness is adopted. The complementary block is not counted in a
physical spectral trace or promoted to a new field or interaction. These
are representation counterexamples demonstrating why a projected lift alone
does not determine the requested physical coupling.

## Decision

**U OPEN, with a rejected direct representation and a verified
non-identifiability certificate.** No physical B or transmitting jet is
assigned. The exact missing rule is a common-time decomposition into genuine
independent bulk subspaces, with compatible embeddings and operator domains.
Only with that rule does

$$
B=P_pH_D[g]P_c
$$

define a physical block and its metric/embedding derivatives. T's duplicate
traces are not those orthogonal bulk projectors. The remaining boundary
action and physical Weyl mismatch consequently remain unevaluated.

The physical U minimum is not met, so this local result does not trigger a
new PDF/version. All locked coefficients and certificates remain unchanged;
no optimizer, source integration or metric timestep is performed.

```sh
python3 scripts/check_nsc_hamiltonian_trace_representation.py --check
```
