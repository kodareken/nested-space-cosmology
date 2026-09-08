# Matching a quantum flow to the operator's actual fields

The finite datum in the [UV match](nsc-torsion-uv-map.md) cannot be fixed by
repeating its leading heat coefficients. The next possible owner is the
common quantum flow. Before computing that flow, this note identifies what
can be imported and the field-space condition that has to be respected.

## The current quantum action is more than one bare spectral trace

The laboratory's one-operator formulation first writes a spectral trace and
a fermion bilinear. Its later quantum definition combines the regulated
determinant with its compensating scale response and a scale-field integral.
The [normalization comparison](nsc-measure-normalization.md) already shows
that the computed proper-time modulus is not that completed functional.

Consequently, the classical relation between the R and K² terms of one a2
coefficient is not by itself a renormalization condition on the full quantum
theory. Imposing c_T=9c_R throughout would require justification from its
measure, constraints and flow. The earlier anomaly and normalization checks
are reused here; they are not rerun and no finite coefficient is set to zero.

## A published flow is available, with a different field space

[Pagani–Percacci, arXiv:1506.02882v3](https://arxiv.org/html/1506.02882v3),
sectionsIII–VI, gives one-loop flows for a generalized Palatini action and
Dirac matter. Its parity-even formulas hold in dimension d. The connection
has quadratic mass terms and is not given a regulating kernel; the metric
flow uses a specified gauge and cutoff. The authors explicitly limit the
UV interpretation of that truncation because connection derivative terms
would then be needed. Its displayed numerical fixed points are in d=4.

We reuse their flow formulas for a compatibility test. They are not the
flow of our five-dimensional proper-time functional, its compact domain or
its coupled KK state. Their cosmological parameter Lambda is also different
from the NSC spectral cutoff with that name.

## Project the connection before assigning its couplings

The source's independent contortion beta_abc is antisymmetric in its last
two indices. Define

\[
I_1=\beta_{abc}\beta^{abc},\qquad
I_2=\beta_{abc}\beta^{bac},\qquad K_{abc}=\beta_{[abc]}.
\]

The current candidate restricts the connection to the skew three-form K.
For the symmetric geometric Dirac operator, vectorial torsion is excluded
and the Cartan-type component is invisible; these are distinct statements.
[Pfäffle–Stephan, Corollary4.6 and Lemma4.7](https://arxiv.org/html/1101.1424v3).
The Hermitian fermion action in Pagani–Percacci also produces a flow term
depending only on the totally antisymmetric component. The exact projection
identity gives

\[
\|K\|^2=\frac{I_1-2I_2}{3}.
\]

Thus C||K||² embeds into this larger ambient action as
g1=C/3, g2=-2C/3 and g3=0. The defining relation is u=2g1+g2=0.
It is not enough to compare the actions only on an already antisymmetric
background: their quantum configuration spaces must also agree.

## The published metric contribution does not preserve this subspace

Let g1,g2 now denote the dimensionless couplings gi/k^(d-2), and call the
paper's gravitational loop factor kappa_loop to distinguish it from the
dimensionful gravitational coupling. Equations41–42 give

\[
\beta_1=-(d-2)g_1+\kappa_{\rm loop}
\frac{d^2-7d-12}{4}g_1-C_f,
\]
\[
\beta_2=-(d-2)g_2+\kappa_{\rm loop}
\frac{(d-4)(d+1)}{4}g_2+2C_f.
\]

For one Dirac field in d=5, C_f=1/(24 pi³) in that prescription. The
projection is then

\[
\left.\beta_u\right|_{u=0}
=-2(d+2)\kappa_{\rm loop}g_1,
\qquad
\left.\beta_u\right|_{u=0,d=5}=-14\kappa_{\rm loop}g_1.
\]

The fermion offsets cancel in this combination. The quoted metric term
generically does not. This is an exact consequence of the imported formulas,
not a derived NSC beta function or a claim that an NSC symmetry is broken.
It means that copying the larger theory's fixed point does not determine
the renormalization of the current Dirac-only connection sector.

## The concrete next formulation step

Use the actual metric, skew three-form, physical fermions and room links as
the independent variables before constructing the flow Hessian. Account
for the quotient of invisible Cartan directions, the restriction on vectorial
torsion, and the associated measure. These are not automatically the same
operation. If a larger field space is adopted
instead, its additional operators and physical meaning must be stated.

For an embedding iota from the selected variables into an ambient set,
even the ordinary off-shell Hessian obeys

\[
H^{\rm selected}_{ab}
=\iota^I{}_{,a}H^{\rm ambient}_{IJ}\iota^J{}_{,b}
+\Gamma_{,I}\iota^I{}_{,ab}.
\]

A field-dependent embedding cannot in general be implemented by dropping
rows and columns after taking the quantum trace. The gauge/constraint
measure and ghosts also belong to that definition. The proper-time matching
window and the finite connection derivative term already computed must be
included where their approximations apply.

The existing data do not yet supply a closed joint flow or its recursive
fixed point. This comparison determines the next calculation's variables
and excludes borrowing a value from a different ambient theory. It does
not determine c_T-9c_R, a quantum state or a self-sourced geometry.

The [exact development record](../results/development/flow-compatibility.json)
checks the generic five-dimensional contortion projection and the algebraic
flow consequence. Its [runner](../scripts/check_nsc_flow_compatibility.py)
authenticates its sources and compares every exact field. Run
`python3 scripts/check_nsc_flow_compatibility.py --check`.
No flow trajectory, spectrum scan or prior calculation is regenerated by
this compatibility check.
