# From the owned Dirac Hamiltonian to endpoint jets

Step J now differentiates the existing midpoint KS propagator and passes the
result through `EndpointBranchJets` and the already verified CTP contraction.
The [record](../results/development/nsc-dirac-endpoint-jets.json) separates
this computed bulk contribution from the absent physical transmitting jets.

The Hamiltonian is exactly the one already used by
`nsc_landau_cauchy_isometry.propagate_blockwise`,

$$
H=N\big(h_{x,0}\sigma_1+h_{y,0}r_0\sigma_2/r+k\sigma_3/a\big)
-\beta kI,\qquad k=h_{z,0}a_0.
$$

Its raw-field vertices for $(N,\beta,a,r)$ are evaluated on the same midpoint
nodes. SciPy's matrix-exponential Fréchet derivative and the ordered-product
chain rule propagate these vertices to $dU/dg_{A,e}$. The CTP difference
coordinate then uses

$$
\frac{\partial U_+}{\partial g_{\Delta,A,e}}
=\frac12\frac{\partial U}{\partial g_{A,e}},\qquad
\frac{\partial U_-}{\partial g_{\Delta,A,e}}
=-\frac12\frac{\partial U}{\partial g_{A,e}}.
$$

This imports the derivative machinery and computes its application to the NSC
mode Hamiltonian. It adds no physical interaction or selected history.

## Focused new check

One fixed middle quadrature block is read from each of the 33 authenticated
channels. The 33-node frozen control is the same diagnostic used by the local
endpoint object. It is used only to check the new derivative; the homogeneous
no-interface class remains excluded as a physical solution. Its duration is
neither changed nor promoted to physical data.

All four fields at both endpoints are differentiated. The verifier compares
the resulting jets with endpoint perturbations passed through the original
propagator, and compares the existing CTP contraction with a finite difference
of the original Gaussian action. Unitarity and tangent residuals are also
retained. The content-addressed NPZ stores the control coordinates, sample
indices, bulk unitaries, endpoint derivatives and Gaussian control covectors.

**Bulk jet component PASS; physical Step J OPEN.** The computed component is
frequency diagonal. It does not contain $dB/dg_\Delta$, the tilted surface's
embedding/normal derivative, or their frequency-mixing contribution. The
initial state and every locked coefficient are unchanged. No total physical
jet, stress or updated constraint is assigned.

The missing data are the transmitting Hamiltonian/link as a functional of
metric and embedding, and the corresponding Cauchy restriction/frame map.
The existing `CausalCommonFunctional` accepts complete matrix vertices when
supplied; its small demonstration vertices do not determine this NSC domain.

```sh
python3 scripts/derive_nsc_dirac_endpoint_jets.py --check
```
