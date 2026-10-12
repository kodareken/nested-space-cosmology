# FGC-1-EFT0-LED1: declared local retained-EFT ledger

**FGC-1-EFT0-LED1** is a fail-closed bookkeeping gate for the FGC-QR action.
It evaluates a few exact local power-counting quantities on a declared
COMP1-centered scalar component box. It does **not** prove EFT validity.

The frozen action is, in natural units and the declared Jordan frame,

\[
S=\int d^4x\sqrt{-g}\left[
\frac{M_{\rm Pl}^2+\beta\phi^2}{2}R
-\frac12(\nabla\phi)^2-\frac12(\nabla\chi)^2
-\frac{\mu^2}{2}\phi^2-\frac{g_4}{4}\phi^4
+\frac{\eta}{8}\phi^2\mathcal G\right].
\]

The model fixture remains FGC-QR with
\(M_{\rm Pl}=2\), \(\mu=3\), \(g_4=1/2\), \(\beta=-1/4\), and
\(\eta=1/2\). Those are code-unit
values; they neither determine nor imply a physical cutoff.
The reproducer binds those exact ledger values to the frozen FGC-1-ACT1
configuration and its passed declared-action gate.  This is source provenance,
not a promotion of ACT1 into an EFT-validity result.

## Declared dimensional ledger

In four dimensions with \(c=\hbar=1\),

\[
[\phi]=[\chi]=[M_{\rm Pl}]=[\mu]=[\Lambda]=1,
\quad [R]=2,
\quad [\mathcal G]=4,
\]
\[
[\beta]=[g_4]=0,
\qquad [\eta]=-2.
\]

The displayed scalar kinetic terms are canonical terms in the declared
Jordan-frame action. They do not by themselves diagonalize the metric--scalar
kinetic sector when \(F(\phi)R\) is active; that is a later principal-symbol
and symmetrizer burden.

The cutoff \(\Lambda=16\) in the frozen config is an **external assumption**
solely for this conditional ledger. It is not inferred from the Planck mass,
the FGC-QR fixture couplings, a local jet width, or a pointwise characteristic
test.

## Operator scope

The config lists all operators retained by the displayed action and a small
representative omitted set:

\[
X_\phi^2/\Lambda^4,\quad X_\chi^2/\Lambda^4,
\quad X_\phi X_\chi/\Lambda^4,\quad \phi^6/\Lambda^2.
\]

Their order-one Wilson bounds and the \(\phi\mapsto-\phi\), \(\chi\) shift
assumptions are inputs, not deductions. The list is explicitly
**representative only**. It is not a complete Wilsonian basis, does not settle
the treatment of all curvature operators under field redefinitions, and does
not predict Wilson coefficients or radiative stability.

This boundary follows the relevant EFT cautions in
[Kovács--Reall, arXiv:2003.08398](https://arxiv.org/abs/2003.08398), especially
their four-derivative EFT action and weak-coupling discussion, and
[East--Ripley, arXiv:2011.03547](https://arxiv.org/abs/2011.03547), whose
four-derivative scalar--tensor form is stated only up to field redefinitions,
total derivatives, and conformal rescalings. Neither source supplies a
universal cutoff or a numerical ACT1 threshold.

## Exact local component evaluation

The declared local component box is centered on the activated COMP1 datum at
\(r=4\):

\[
\phi=\partial_r\phi=\frac1{131072},
\qquad
|\delta\phi|,|\delta(\partial_r\phi)|\leq\frac1{262144}.
\]

Its upper component bounds are strictly below QIFT1's parameter half-width
\(1/65536\), and its lower component bounds remain strictly positive. This
only verifies that these two declared components form a strict subbox of the
QIFT1 parameter interval; it does not turn the component box into QIFT1's
full thirty-dimensional jet box.

Only quantities determined by this scalar component interval are evaluated
exactly. With \(|\phi|\le3/262144\), the ledger checks

\[
F_{\min}=M_{\rm Pl}^2-|\beta||\phi|_{\max}^2,
\qquad
\frac{F_{\min}}{M_{\rm Pl}^2}>\frac{99}{100},
\qquad
\epsilon_F=\frac{|\beta||\phi|_{\max}^2}{M_{\rm Pl}^2},
\]
\[
\epsilon_{\eta\phi^2}=\eta|\phi|_{\max}^2,
\qquad
\epsilon_{\phi}=|\phi|_{\max}/\Lambda,
\]

against strict declared bounds. The dimensional quantity \(F_{\min}\) is
serialized only for provenance; the gate applies its explicit strict **lower**
limit to the dimensionless ratio \(F_{\min}/M_{\rm Pl}^2\). Every
\(\epsilon\) instead has an explicit strict **upper** limit of \(1/100\).
Every number is a canonical rational string;
floats, missing cutoff/Wilson assumptions, inconsistent dimensions, a promoted
claim, or any control at or above its strict bound reject the ledger.

The radial derivative interval is preserved for provenance but is **not**
relabelled as a covariant derivative norm: a coordinate component at one local
jet does not provide the physical frame norm required for a general derivative
expansion estimate.

## Explicit nonclaims

The certificate serializes the full retained table, the representative omitted
table, each Wilson-bound assumption, and the action-parameter references. The
representative omitted remainders are declared but **not locally evaluated**:
the component-only box lacks covariant derivative, curvature, and frequency
data needed for a remainder bound.

The certificate always keeps `retained_eft_validity=false`. In particular it
does not establish a UV completion, a complete operator basis, Wilson
coefficient prediction, radiative stability, frequency or Fourier-support
control, global or open-domain EFT validity, compatible initial data,
evolution, collapse, affine defocusing, or singularity resolution. A finite
local component/jet box cannot derive a cutoff or bound proper frequencies; it
can only verify the listed interval inequalities under the declared external
assumptions.

Run the focused checks with:

```sh
python3 -m unittest tests/test_fgc_eft_ledger.py
```
