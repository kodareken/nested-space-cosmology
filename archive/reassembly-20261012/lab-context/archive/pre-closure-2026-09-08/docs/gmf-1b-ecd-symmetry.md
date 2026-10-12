# GMF-1B-ECD-SYM1: a spherical Einstein--Cartan--Dirac axial-current certificate

GMF-1B-ECD-SYM1 fixes a candidate **axial-current-sector** preflight for the next
GMF-1B gate. It does not solve a bounce, a radial spacetime, a junction, or a
global topology problem. Its purpose is narrower: determine whether a
two-spinor source can remain rotationally symmetric while retaining a nonzero
minimal Einstein--Cartan axial contact invariant. It does not yet prove that
the complete ECD source, including its nonlinear spinor term and effective
stress, is SO(3)-equivariant.

The candidate is the massless left-chiral two-spinor pair of Ventrella and
Choptuik, used there to form a rotationally symmetric Einstein--Dirac source.
Their calculation is a torsion-free Einstein--Dirac result; this document
freezes the additional conventions and tests required before promoting that
pair to an Einstein--Cartan--Dirac (ECD) source. [Ventrella and
Choptuik (2003)](https://arxiv.org/abs/gr-qc/0304007)

## Frozen minimal-ECD convention and its explicit bridge

Version 0.10 reported the correct real-current fixture but did not explicitly
show the factor of `-i` and the global signature bridge connecting the
published Ventrella--Choptuik matrices to the action convention. Version 0.11
corrects that presentation; the earlier numerical values are unchanged.

The fundamental variables are a tetrad `e`, an independent metric-compatible
spin connection `omega`, and two independent classical complex spinor fields
`psi_+`, `psi_-`. Define the minimal ECD action first in a canonical
`(+,-,-,-)` layer with Dirac matrices `Gamma`,

\[
\{\Gamma^{\hat a},\Gamma^{\hat b}\}=2\eta_+^{\hat a\hat b},
\qquad
\Gamma^5=i\Gamma^{\hat0}\Gamma^{\hat1}\Gamma^{\hat2}\Gamma^{\hat3},
\qquad \kappa=8\pi G,
\]

and total real axial current

\[
J_+^{\hat a}=\sum_{I=\pm}\bar\psi_I
\Gamma^{\hat a}\Gamma^5\psi_I.
\]

Combining the torsion-dependent Dirac term with the gravitational contortion
term gives the complete torsion-eliminated interaction

\[
S_{\rm eff}^{(+)}=S_{\rm EH}^{(+)}+S_{\rm D,LC}^{(+)}
-\int d^4x\,e\,\frac{3\kappa}{16}J_{+\hat a}J_+^{\hat a}.
\]

This is not obtained from the intermediate torsion-substituted Dirac term
alone. The source chain and nonlinear equation are given in Cabral, Lobo, and
Rubiera-Garcia, with their notation `kappa^2=8 pi G` translated here to the
project's `kappa=8 pi G`. [Cabral, Lobo, and Rubiera-Garcia
(2019)](https://arxiv.org/abs/1902.02222)

The VC fixture instead uses imaginary `(-,+,+,+)` matrices. The complete
matter-sector translation is

\[
\Gamma^{\hat a}=-i\gamma_{\rm VC}^{\hat a},\qquad
\bar\psi=\psi^\dagger\Gamma^{\hat0}
=\psi^\dagger(-i\gamma_{\rm VC}^{\hat0}),\qquad
\Gamma^5=\gamma_\star,
\]

\[
A_-^{\hat a}=J_+^{\hat a}
=-i\sum_{I=\pm}\bar\psi_I\gamma_{\rm VC}^{\hat a}
\gamma_\star\psi_I,qquad J_+^2=-A_-^2.
\]

Therefore the same interaction in VC variables is

\[
S_{\rm int}^{(-)}=+\int d^4x\,e\,
\frac{3\kappa}{16}A_{-\hat a}A_-^{\hat a},
\]

and its massless equation contains

\[
\gamma^a_{\rm VC}D_a\psi
-i\frac{3\kappa}{8}A_{-a}\gamma^a_{\rm VC}\gamma_\star\psi=0.
\]

This is a global convention translation, not two physical metrics. Tetrad,
epsilon/orientation, index, curvature/action, and stress signs must transform
together. GMF-1B-ECD-INT1 makes the Clifford, adjoint, gamma-five, norm,
contact-action, and nonlinear-equation residuals executable. The algebraic
spin--torsion coupling is the Hehl--Datta mechanism; ECD is the
gauge-theoretic gravity framework introduced by Kibble. [Kibble
(1961)](https://doi.org/10.1063/1.1703702); [Hehl and Datta
(1971)](https://doi.org/10.1063/1.1665738)

The contact term is lower order in the spinor evolution. It changes the
nonlinear source and stress tensor, but not the first-order principal symbol
of the minimally coupled Dirac operator; the candidate therefore does **not**
derive a new metric-null cone, variable locally measured `c`, or a relative
photon/tensor cone.

## Exact Ventrella--Choptuik fixture

In polar-areal variables,

\[
ds^2=-\alpha(t,r)^2dt^2+a(t,r)^2dr^2+r^2d\Omega_2^2,
\]

take the displayed left-chiral pair from the source paper

\[
\psi_+=\frac{e^{i\phi/2}}{2\sqrt{\pi}\,r\sqrt a}
\begin{pmatrix}
F\sin(\theta/2)\\G\cos(\theta/2)\\
F\sin(\theta/2)\\G\cos(\theta/2)
\end{pmatrix},
\qquad
\psi_-=\frac{e^{-i\phi/2}}{2\sqrt{\pi}\,r\sqrt a}
\begin{pmatrix}
F\cos(\theta/2)\\-G\sin(\theta/2)\\
F\cos(\theta/2)\\-G\sin(\theta/2)
\end{pmatrix},
\]

where `F(t,r)` and `G(t,r)` are complex radial functions. The fixture includes
the source paper's regular-centre starting conditions `F(t,0)=G(t,0)=0`; an
ECD implementation must additionally prove smoothness of every torsion
bilinear in a regular Cartesian tetrad.

Write

\[
\mathcal C=\frac{1}{2\sqrt{\pi}\,r\sqrt a}.
\]

With the frozen convention, the dependency-free direct gamma-matrix
certificate verifies

\[
A^{\hat0}=2\mathcal C^2(|F|^2+|G|^2),
\qquad
A^{\hat r}=2\mathcal C^2(|F|^2-|G|^2),
\qquad
A^{\hat\theta}=A^{\hat\phi}=0,
\]

and hence

\[
A_{\hat a}A^{\hat a}=-16\mathcal C^4|F|^2|G|^2.
\]

The implementation evaluates the displayed spinors directly in the source
paper's imaginary gamma representation and records the binary64 residual from
these formulas. It also asserts these three controlled cases:

```text
F = 0, G != 0          -> A != 0 is null, torsion source != 0, A^2 = 0
G = 0, F != 0          -> A != 0 is null, torsion source != 0, A^2 = 0
F != 0, G != 0         -> A != 0, torsion source != 0, A^2 < 0
F = 0, G = 0           -> A = 0, torsion source = 0, A^2 = 0
```

The `A^2<0` statement is an axial-contact fact, not an energy-condition,
bounce, or defocusing result. Equally important, `A^2=0` in either
one-amplitude control does **not** mean that the axial current or algebraic
Cartan torsion vanishes: the current is nonzero and null. It only makes the
scalar contact-action density zero at that field value; it does not by itself
show that the variation of the contact action, or the Hehl--Datta nonlinear
equation, vanishes.

For an exact algebraic regression, set the normalization to one, `kappa=1`,
`F=3/4`, and `G=1/2`. Then

\[
A^{\hat0}=\frac{13}{8},\qquad
A^{\hat r}=\frac{5}{8},\qquad
A^2=-\frac94,\qquad
\frac{3\kappa}{16}A^2=-\frac{27}{64}.
\]

The result stores this separately from the geometric fixture
`r=1/(2 sqrt(pi))`, `a=1`, for which `C` equals one to binary64
tolerance. It does not interpret that regression radius as a physical
transition scale.

## Why the pair is SO(3), not necessarily O(3)

The individual spinors are not rotationally invariant. Their spin-angular
momenta cancel only in the summed source. Ventrella--Choptuik directly derive
the spherical summed torsion-free stress tensor; the ECD extension must check
the same property for the total axial current, contorsion, nonlinear Dirac
term, and effective stress tensor.

A scalar `A_hat0(t,r)` and a radial axial vector
`A_hatr(t,r) rhat` are SO(3)-equivariant. Vanishing angular components are
therefore compatible with spherical metric data. A nonzero radial axial vector
need not be invariant under spatial parity. GMF-1B-ECD-SYM1 consequently
certifies SO(3) compatibility of the total axial-current/bilinear sector,
**not** SO(3) equivariance of the full ECD source and not O(3) or parity
invariance.

The symmetry proof must establish

\[
T_{t\theta}=T_{t\phi}=T_{r\theta}=T_{r\phi}=T_{\theta\phi}=0,
\qquad
T_{\phi\phi}=\sin^2\theta\,T_{\theta\theta},
\]

for the **full effective ECD** stress tensor, as well as

\[
A^{\hat\theta}=A^{\hat\phi}=0
\quad\hbox{and}\quad
\partial_\theta A^{\hat0}=\partial_\phi A^{\hat0}
=\partial_\theta A^{\hat r}=\partial_\phi A^{\hat r}=0.
\]

INT1 now reconstructs the local Hehl--Datta interaction pointwise in the
declared opposite-chirality retained basis at three nondegenerate tested angles
with a small numerical residual. This is not a complete spin-weighted harmonic
projection or a proof that every omitted `(ell,m)` coefficient vanishes. The
complete symmetry proof still requires the kinetic and
spin-connection terms, their torsion-free Dirac tetrad variation, and the
complete effective-stress identities. KIN1 has since closed only the algebraic
contact-sector variation. A spherical metric ansatz or interaction-only
projection is not the full symmetry proof.

## Full-Dirac cancellation is a different comparator

The static full-Dirac singlet construction used by Finster, Smoller, and Yau
is a useful comparator, not this fixture. In the equal-profile,
opposite-spin full-Dirac comparator, the total axial current cancels
pointwise, giving `A^a_total=0` and therefore no minimal-EC contact term. That
does **not** imply cancellation for the Ventrella--Choptuik left-chiral pair:
the latter has a nonzero total temporal/radial axial current whenever both
radial amplitudes are populated.

No single-spinor Fierz identity implies

\[
\left(A_+^{\hat a}+A_-^{\hat a}\right)
\left(A_{+\hat a}+A_{-\hat a}\right)=0.
\]

Single-Weyl currents are null, but the sum of two independent nonparallel null
currents can be timelike. The executable control evaluates the actual
equal-profile Pauli doublet with `gamma^5=i gamma^0 gamma^1 gamma^2 gamma^3`,
records a nonzero individual-spinor axial current, and then records the summed
cancellation residual; it does not merely return zero by declaration. The
source-level fixture, rather than an analogy to the full-Dirac comparator, is
the decision procedure. [Finster, Smoller, and
Yau (1998)](https://arxiv.org/abs/gr-qc/9801079)

## Classical model boundary

`F` and `G` are ordinary complex c-number fields in this proposed numerical
model. This is a classical multi-spinor ECD effective theory. It is not a
derivation of a quantum two-fermion state, a Fermi gas, a Weyssenhoff ensemble,
or a renormalized expectation value. In a quantum calculation one would need
to specify a state and distinguish

\[
\langle\hat A_{\hat a}\hat A^{\hat a}\rangle
\ne
\langle\hat A_{\hat a}\rangle
\langle\hat A^{\hat a}\rangle
\]

in general. No quantum closure is assumed here. Treating the fields as
Grassmann-valued variables would also not be a direct replacement for a
numerical c-number evolution.

## Support and boundary discipline

Compact annular initial support is allowed by a finite-speed first-order
hyperbolic Dirac evolution, but it is **not** a persistent material boundary.
The implemented profile uses

\[
B(r)=
\begin{cases}
\exp\!\left(4-\dfrac{1}{s(1-s)}\right),
&s=\dfrac{r-r_{\rm in}}{r_{\rm out}-r_{\rm in}}\in(0,1),\\[4pt]
0,&\text{otherwise},
\end{cases}
\]

with `F=F_0 B` and `G=G_0 B`. It is exactly zero at and outside
the two edges and exactly one at the midpoint. The regression uses
`r_in=1`, `r_out=2`, `F_0=3/4`, `G_0=1/2`, and 257 sample points. The profile
therefore provides a smooth spinor taper and exact initial spinor-zero buffers.
It does not establish that those regions solve the gravitational constraints
or are vacuum spacetime regions. In a
later finite-time evolution, a timelike matching worldtube would have to be
chosen inside a sufficiently wide constraint-compatible vacuum buffer, and the implementation would
have to demonstrate that the matter domain of dependence does not reach that
worldtube during the declared interval.

Consequently this source certificate does not establish a shell-free junction,
fixed exterior mass, or a permanent compact object. Those are later matching
conditions, not consequences of compact initial support.

## Computed content and exact nonclaims

GMF-1B-ECD-SYM1 computes only a convention-locked contact scalar, two direct
gamma-matrix axial certificates, the double-null Dirac principal polynomial,
and a finite smooth annular profile audit. The machine record is regenerated
with

```bash
python3 scripts/reproduce_gmf1b_ecd_symmetry.py \
  --output results/gmf-1b-ecd-symmetry.json
```

It makes **no** claim of:

- a solved Einstein--Cartan radial evolution or constraint solution;
- a regular bounce, metric-null defocusing, or trapped-to-anti-trapped path;
- finite curvature, finite torsion, a regular centre, or EFT control;
- SO(3) equivariance of the full nonlinear ECD source or effective stress;
- a constraint-compatible vacuum spacetime region merely because the spinor
  profile is zero there;
- a shell-free finite-mass exterior, global matching, external flux, or `Q`;
- radial, nonradial, inner-horizon, or mass-inflation stability;
- a topology change, a complete `S^3` child, or global geodesic completeness;
- dark-matter abundance, dark-energy origin, a varying local `c`, or a
  relative-cone signal; or
- observational confirmation.

In particular, `A^2<0` does not prove a negative effective energy density or
the Raychaudhuri condition required for metric-null defocusing. Those signs
must be calculated from the full effective stress tensor.

## Next gates: GMF-1B-ECD-INT1, then GMF-1B-ECD-ID1

INT1 now supplies the explicit action/signature bridge, direct Hehl--Datta
local retained-basis reconstruction gate documented in
[gmf-1b-ecd-identity.md](gmf-1b-ecd-identity.md). It does not retroactively
turn SYM1 into a full-source calculation.

GMF-1B-ECD-ID1 is an initial-data and constraint step, not yet an evolution.
Starting from the frozen ECD action and a source certificate that passes the
tests above, it must construct a regular-centre, finite-Misner--Sharp-mass
spherical initial slice with:

1. the full Hamiltonian, momentum, tetrad, and algebraic-torsion constraints
   satisfied to declared residual tolerance;
2. a smooth compact/tapered spinor region plus a nonzero spinor-zero annulus,
   with the constraints separately proving whether it is a vacuum region;
3. a declared outer vacuum worldtube and no distributional surface source on
   the initial slice;
4. finite initial curvature, torsion, stress, and spinor bilinears; and
5. a declared SO(3), parity, classical-field, and asymptotic/exterior scope.

Failure to produce such data rejects this candidate before a `1+1` evolution
solver is interpreted as GMF-1B progress.
