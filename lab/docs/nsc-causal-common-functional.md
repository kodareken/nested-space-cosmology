# The canonical causal common functional

Real-time physics is now assigned to the canonical Lorentzian Dirac theory.
The spectral operator supplies the same field content, compact spectrum,
normalization and state-independent induced coefficients. This yields

\[
\Gamma_{\rm one}^{\rm CTP}
=S_{\rm induced}[g_+,A_+]-S_{\rm induced}[g_-,A_-]
-i\log\det(I-C_0+C_0U_-^\dagger U_+).
\]

The complete parent–child Hamiltonian, including gauge and link histories,
generates \(U_\pm\). The initial covariance \(C_0\) retains parent, child and
cross-boundary correlations. The previously continued bare Euclidean vertex
is not used as a retarded correlator.

## Executable interface

`CausalCommonFunctional.evaluate` accepts metric, gauge and link histories,
the initial covariance, the induced spectral source and the declared boundary
domain. It returns the normalized CTP action, all four metric forces, their
stress projection, boundary power, retarded response, noise, discrete Ward
residuals and any unresolved terms.

The finite owner evolves the covariance unitarily. At each history step the
Hamiltonian change performs work

\[
\Delta W=\operatorname{Tr}[C\,(H_{n+1}-H_n)],
\]

so that the recorded discrete energy/work identity is exact up to numerical
matrix error. Metric sources use the CTP sign

\[
F_A^{\rm matter}=-\operatorname{Tr}(C\,\partial_AH),
\qquad
F_A=F_A^{\rm matter}+F_A^{\rm induced}.
\]

The induced contribution is supplied once by one named coefficient owner.

## Binding result

A finite two-mode control verifies equal-history normalization, unitary
covariance evolution, causal response, energy/work balance and one-counted
induced forces. The magnitudes of the stored local neck source are used only
to test source accounting; the two-mode matrices are not a physical throat.

The existing equilibrium calculation is reused rather than rerun. Its direct
Euclidean frequency loop and canonical retarded response at imaginary
frequency differ by at most \(5.6\times10^{-16}\) on the recorded finite grid.
The finite proper-time kernel remains distinct from the canonical Kubo kernel,
as already established.

## Next physical input

The interface is executable, but the physical charged parent–child history is
not yet assembled. Its first dependencies are the complete common-action
vacuum coefficient \(V_{\rm full}\), then the charged Hamiltonian, covariance
and closed magnetic return history on the MMP-compatible domain. Those inputs
must come from the embedding gate; they are not filled by the finite control.

## Reproduction

```sh
python3 scripts/derive_nsc_causal_common.py --check
```

The command reads authenticated stored records and performs only the new
finite CTP accounting check. It does not launch the angular, vacuum-stress,
MMP gravity or historical result generators.
