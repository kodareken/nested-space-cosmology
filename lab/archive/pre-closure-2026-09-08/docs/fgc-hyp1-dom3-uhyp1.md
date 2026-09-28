# FGC-1-HYP1-DOM3-UHYP1: compact radial branch-graph hyperbolicity

**FGC-1-HYP1-DOM3-UHYP1** is an exact-rational certificate for the radial
principal system of the unredefined REF1 modified-harmonic formulation. It
works on the compact graph of the already interval-certified implicit
coordinate-time acceleration branch through the activated COMP1 datum at
coordinate radius \(r=4\). It proves a real, complete, uniformly invertible
12-column radial first-order eigenframe and the bounded radial symmetrizer
\(H=V^{-T}V^{-1}\).

The frozen widths are \(d_z=2^{-300}\), \(d_a=2^{-290}\), root-bracket
half-width \(2^{-200}\), hat-pivot displacement \(2^{-140}\), and regulator
Krawczyk displacement \(2^{-68}\). All are strictly positive. Exact
`Fraction` interval endpoints and interval forward AD are used throughout;
there is no floating-point or sampled sign, inverse, bracket, or frame test.

## Constructive sector proof

The physical, tilde, and hat metric-value cone quadratics each have two
rational root brackets. Every bracket passes endpoint-sign, strict
monotonicity, and positive-discriminant gates. Tilde modes use the conditional
universal pure-gauge right identity. Hat modes use selected rows
\((E_{tr},E_\theta,E_\phi)\), pivot columns \((h_{tt},h_{tr},h_{rr})\),
three-variable Krawczyk inclusions, scalar identities, and the conditional
radial Bianchi row-completion identity. The \(\phi\)-free hat column is
scaled by \(2^{-22}\) only for frame conditioning.

The physical modes are the decoupled \(e_\chi\) modes at metric-null roots.
The regulator is not derived by dividing an interval determinant: for each
sign, rows 0--4 of \(P(c)v=0\) are solved by a normalized five-variable
Krawczyk map with \(\phi=1,\chi=0\), seeded from the exact COMP1 quotient
chart. The sixth row follows from the chi identity. Its speed enclosure must
stay sign-fixed and disjoint from the corresponding physical root bracket.

The twelve second-order modes lift as \((-cv,v)\). A fixed exact center
inverse and a strict Neumann frame-defect bound certify every enclosed frame
is invertible. This yields explicit Euclidean coercivity bounds for
\(H=V^{-T}V^{-1}\) and the radial symmetry identity
\(HA=V^{-T}\operatorname{diag}(c_i)V^{-1}=(HA)^T\).

The public JSON stores every decisive exact interval, inclusion margin, mode
route, speed enclosure, Neumann bound, and symmetrizer bound. Very large
intermediate matrices are represented by SHA-256 commitments rather than a
multi-megabyte decimal dump. The reproducer constructs the complete exact
payload first, hashes its canonical serialization, and only then projects the
review-sized record; rerunning it therefore checks the omitted matrices as
well as the published bounds.

## Exact scope

This is a compact nonflat **radial** strong-hyperbolicity theorem on the
implicit REF1 branch graph. The stored record hashes the MODE1, QIFT1, and
COMP1 predecessor configurations/results and all relevant source files.

It does **not** prove multidirectional strong hyperbolicity, constraint
propagation, compatible IBVP data, retained-EFT validity, an evolution,
collapse solution, null affine defocusing, or singularity resolution. No PDE
evolution authorization follows from this local radial-principal certificate.
