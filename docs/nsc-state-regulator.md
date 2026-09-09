# State dependence of the finite proper-time conversion

The adopted finite Euclidean proper-time factor cannot be identified,
without a state-dependent conversion, with the thermal source of the
canonical Dirac field used in the parent-state calculation. This is an
explicit compatibility result, not a new thermal-field theorem or a
thermal approximation to the expanding child.

The distinction already matters for the simplest equilibrium test. At
T/nu=0.5, a single massless four-dimensional Dirac field has canonical
thermal energy density 0.0719658654*nu^4. Applying the lower proper-time
endpoint 1/nu^2 to the **entire** thermal determinant gives -0.0012266736*nu^4
after subtracting its own zero-temperature energy. Those are different
state/source prescriptions. No local vacuum coefficient can change their
thermal difference.

## Reuse and scope

Use [Gusev–Zelnikov, equations 83–89](https://arxiv.org/html/hep-th/9807038):
the fermionic determinant, antiperiodic thermal images, and canonical
thermal free energy. Apply their heat-image representation to the
already adopted finite lower endpoint; do not rederive finite-temperature
QFT. The [canonical influence functional](nsc-influence.md) is the existing
Klich/Fock trace owner. The [parent-state record](nsc-unruh-state.md) and
[compact matching](nsc-compact-matching.md) are read and authenticated,
without rerunning their generators.

The earlier [spatial-circle control](nsc-covariant-measure.md) already
demonstrated the finite omitted winding contribution. That remains a
spatial spin-structure result. Here the genuine thermal circle adds the
canonical-state energy comparison; a spatial AP condition is not relabeled
as a temperature.

The control domain is flat R3 times an antiperiodic Euclidean thermal
circle, with one ordinary massless Dirac field and hbar=c=k_B=1.
nu is the light-sector matching cutoff. It is neither a measured thermal
ceiling nor the carrier cutoff Lambda. This calculation does not assign
a local temperature to the spacelike interior Killing generator.

## The finite conversion is explicitly state dependent

Write F for free-energy density. The stored proper-time light factor is

\[
\Gamma_{\mathrm{light},\nu}=\tfrac12\operatorname{Tr}
E_1(D_E^2/\nu^2).
\]

For the thermal images, its vacuum-subtracted free energy is

\[
\Delta F_\nu(T)=\frac{4T^4}{\pi^2}
\sum_{n=1}^\infty\frac{(-1)^n}{n^4}
\left[1-(1+x_n)e^{-x_n}\right],\qquad
x_n=\frac{n^2\nu^2}{4T^2}.
\]

The canonical value is Delta F_can=-7*pi^2*T^4/180. Therefore the
thermal part of the existing conversion C_nu,mu=Gamma_light,nu-Gamma_light,ren
is fixed in this control:

\[
\boxed{\Delta F_{\mathcal C}(T)
=-\frac{4T^4}{\pi^2}
\sum_{n=1}^\infty\frac{(-1)^n}{n^4}(1+x_n)e^{-x_n}.}
\]

It is nonzero at finite nu. Its metric/thermal derivatives must also be
retained. In equilibrium the lapse variation gives rho=F-T*dF/dT,
pressure=-F and entropy=-dF/dT. For the endpoint source,

\[
\rho_\nu-3p_\nu
=\frac{\nu^4}{2\pi^2}
\sum_{n=1}^\infty(-1)^n e^{-x_n}.
\]

The canonical thermal trace is zero. A finite endpoint has introduced a
state-dependent trace term. The same free energy and density are obtained
independently from the direct Matsubara heat integrals, with
omega_j=(2j+1)*pi*T. This compares the actual finite prescription; it
does not replace the proper-time endpoint with a spatial energy filter.

| T/nu | Canonical thermal rho/nu^4 | Finite-endpoint thermal rho/nu^4 |
|---:|---:|---:|
| 0.10 | 0.0001151454 | 0.0001151454 |
| 0.25 | 0.0044978666 | 0.0031350475 |
| 0.50 | 0.0719658654 | -0.0012266736 |
| 1.00 | 1.1514538468 | -0.0063244282 |

At very small T/nu the discrepancy is exponentially small, but it is not
an exact identity. The record also evaluates a flat control at the stored
parent Hawking temperature, retaining the logarithm of its tiny mismatch.
That number is not an estimate for the actual curved child source.

## Why a pole-residue argument is insufficient

The regulated Euclidean inverse contains e^(-z/nu^2)/z, whose residue at
z=0 is one. A homogeneous Lorentzian bisolution W satisfying D*W=0 is
formally unchanged by a regulator function with value one on shell,
provided that function is defined on W. This fact does not equate a
regulated imaginary-time determinant to a canonical finite-temperature
trace. Their finite contour and state contributions also matter.

The checker makes this concrete for a pair of energies +/-E. The
canonical thermal derivative of free energy is 2/(exp(E/T)+1). The
derivative of the finite Matsubara determinant differs even though the
pole residue is unchanged. It agrees with the thermal-image expression.
Thus the missing contribution cannot be discarded using the residue alone.

[Barvinsky, equations 2.25–2.29](https://arxiv.org/html/1408.6112)
gives an established Euclidean-to-in-in rule under specified vacuum and
boundary assumptions. Those assumptions, and the regulator's continuation,
must be supplied here. The Unruh characteristic state is not identified
with the theorem's past Poincare-invariant vacuum by notation.

## A source-budget constraint, not another free coefficient

For a canonical Gibbs state, each +/-E pair's energy above the ground
state is 2E/(exp(E/T)+1)>=0. Hence the negative vacuum-subtracted thermal
energy above cannot be that canonical field's expectation value. This
does not invalidate proper-time regularization: removing the endpoint
recovers the standard result. It forbids treating the raw finite factor
as the entire physical state functional with the unchanged canonical field.

The same issue cannot be removed merely by reallocating the light/heavy
split. For any positive finite raw thermal modulus F_raw(T) for which
F_raw and T*dF_raw/dT tend to zero at large T, its thermal excitation
energy tends to -F_raw(0). Fermionic Matsubara frequencies have no zero
mode. The adopted positive compact heat weight has an exponentially
suppressed large-frequency tail, so the free flat carrier has this limit.
The already matched vacuum moment gives F_raw(0)=Q2[h]/(8*pi^2),
approximately 0.17887344 in the stored Lambda=2, ell=2 normalization.
This is a high-temperature compatibility limit of that free raw modulus,
not a computed finite-temperature stress on the child or a claim that
the current low-energy realization remains valid at arbitrary temperature.

On the same flat background, every state-independent local vacuum term
cancels in the thermal-minus-vacuum difference; all curvature invariants
vanish. On a fixed curved background local counterterms likewise cancel
between states, as in [DHP Definition 4.1(1)](https://arxiv.org/html/0904.0612).
They cannot supply this correction. The owning term is
the **state/contour completion**, or an explicitly changed microscopic
operator and state construction with its own physical tests.

Preserve the exact allocation

\[
\Gamma_5^{\rm CTP}
=\Gamma_{\rm light,ren}^{\rm CTP}[C]
 +\Gamma_H^{\rm CTP}[C]
 +\mathcal C_{\nu,\mu}^{\rm CTP}[C].
\]

If the raw endpoint prescription is retained on a thermal contour,
Delta C is the nonzero expression above. If the intended physical field
retains its canonical state functional, the completion must satisfy

\[
\Delta\Gamma_{\rm completion}^{\rm CTP}
=\Delta\Gamma_{\rm canonical}^{\rm CTP}
-\Delta\Gamma_{\rm raw}^{\rm CTP}
\]

for that sector. This is a matching condition, not permission to insert
a fitted source. The required thermal correction is calculated here;
the full nonthermal curved conversion has not been evaluated.

In particular, the previously computed negative neck null stress is a
canonical massless-state result. It remains valid in its stated scope.
Neither its value nor the parent power can yet be substituted for the
complete finite-cutoff source. The next owner is the common contour and
initial-state prescription, including the massive, boundary and recursive
sectors. Another angular refinement or a finite local coefficient scan
does not provide it. A ratio of two normalized influence functionals is
only a relative functional; positivity of the complete influence functional
must still follow from its underlying state and evolution.

## Evidence

The [immutable record](../results/development/state-regulator.json) contains
the imported normalization, independently evaluated thermal source,
metric derivatives, pole/state comparison, finite conversion and scope.

```sh
python -B scripts/check_nsc_state_regulator.py --check
```

Only this new compatibility application runs. Old spectra, source tensors
and canonical influence tests remain authenticated stored inputs.
