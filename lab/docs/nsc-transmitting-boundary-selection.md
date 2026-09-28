# Transmitting action: the kernel selection is OPEN

The full extended gate binds the existing boundary declaration to the actual
33-channel, 1,904-block Cauchy state. The declared balance is

$$
\pi_p^{ab}+\pi_c^{ab}+\frac{\delta\Gamma_{\rm rest}}{\delta h_{ab}}=0.
$$

The [selection record](../results/development/nsc-transmitting-boundary-selection.json)
authenticates the relevant implementations and records their arguments and
domains. `CausalCommonFunctional` consumes `LinkHistory.values` and the metric
vertices. `GaussianBoundaryState` consumes the link and child propagators.
The spatial throat owner computes transmission on its fixed spatial interval;
the Landau propagator evolves a supplied history. None of these existing
inputs specifies the transmitting link as a functional of the tilted embedding
and metric, or supplies its associated variations.

The compact heat-boundary record explicitly states that it has not derived
the transmitting action. Its reflecting compact endpoints cannot supply the
missing radial/tilted boundary kernel. These are source-domain findings, not
a new theorem about boundary actions in general.

## Apply the remaining freedom to the recorded state

Reuse the [weighted-interface family](nsc-tilted-landau-interface.md),

$$
K_c=G_{L,c}^{-1/2}V_cG_{0,c}^{1/2},\qquad V_c\in U(2n_c).
$$

Its real dimension is 447,488 before any quotient by gauge or state
stabilizers. The new calculation tests whether this freedom can change the
locked covariance in its fixed recorded basis. In each channel it uses one
infinitesimal Hermitian generator coupling the first and last retained
frequency nodes, normalized by $\|X_c\|_F=1$. At the identity control,

$$
\delta V_c=iX_c,\qquad
\delta(V_c^\dagger V_c)=0,\qquad
\delta C_c=i[X_c,C_{c,0}].
$$

The record evaluates both derivatives from the authenticated covariance.
Nonzero $\delta C_c$ means the norm condition alone does not determine the
transported state in that fixed basis. The probe is an identifiability test;
it selects no physical map and evaluates no stress. A simultaneous basis
change of both state and observables is a different operation. The dimension
after quotienting all physical redundancies is not asserted here.

## Selection result

**OPEN:** the physical selection residual is unevaluated (`null`, tolerance
$3\times10^{-11}$). The exact missing data are the mode-resolved transmitting
link/embedding functional in the retained bases and its channel and metric
variations. Spatial identity matching does not set the propagator between
distinct spacetime hypersurfaces to the identity. Equal-history CTP
normalization also does not replace the difference-history variation.

The executable `ModeResolvedTransmittingBoundaryHistoryAction` binds and
audits these dependencies. It does not implement an absent physical
$\Gamma_{\rm rest}$ by choosing an element of the admissible family.

Continue with the [endpoint binding](nsc-weyl-endpoint-match.md) and
[full extended decision](nsc-full-extended-history-gate.md).

```sh
python3 scripts/derive_nsc_transmitting_boundary_binding.py --check
```

This focused verifier reads existing records and the small covariance payload;
it invokes no historical generator or metric evolution.
