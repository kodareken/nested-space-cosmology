# Full extended gate: OPEN at the transmitting action

The [full extended record](../results/development/nsc-full-extended-history-gate.json)
executes Douglas's A→D sequence using the three completed component records
and two new binding records.

| Step | Executed result | Status |
|---|---|---|
| A — select $V_c$ | Audited existing action inputs; applied remaining kernel freedom to the serialized covariance | OPEN: no physical transmitting selector supplied |
| B — match the endpoint | Extracted all eight Weyl end-node coefficients, with zero extraction residual; defined the common variation basis | OPEN: transmitting derivative and metric pullback absent |
| C — compose | Bound reference, node-wise local action and tilted interface to explicit evaluated/missing residual channels | OPEN composition |
| D — stationarity | Evaluated dependency gate before an optimizer | **OPEN** |

The new contribution is the application to the actual Cauchy covariance and
the explicit endpoint binding. Existing Bloch/reference, local-gradient,
unitarity/CAR and homogeneous non-existence checks are imported from their
authenticated records; their generators are not rerun.

The residual surface has four metric slots,

$$
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}
{\delta g^A_{\Delta}(\tau)},
\qquad A\in\{N,\beta,q_{\rm ADM},r\},
$$

plus the transmitting selector and endpoint matching slots. The latter two
remain unevaluated. Consequently no numerical full-history variation is
assembled and no existence or non-existence is claimed for the extended
class. A successful numerical endpoint sum would still require the remaining
history equations; the code never promotes that subcheck into existence.

The old homogeneous no-interface certificate remains a trusted, passing
non-existence input with its recorded shift residual
$-0.0135927088490713$ (tolerance $3\times10^{-11}$). It is not applied as a
non-existence proof for the tilted class.

The exact missing dependency remains the physical content of
`ModeResolvedTransmittingBoundaryHistoryAction`: the link/embedding functional
that selects channel kernels, and the derivative of that same action on the
declared boundary variation space. This run implements its **binding audit**,
not an arbitrary replacement action. The compact reflecting boundary and
static spatial curvature owners retain their original domains.

$A,q,\Omega,\zeta,V_{\rm full}$ and the local coefficient allocation are
unchanged. Finite stress, null components and updated constraints remain
unevaluated; metric stepping and coupled evolution remain closed.

```sh
python3 scripts/derive_nsc_transmitting_boundary_binding.py --check
python3 -m pytest -q tests/test_nsc_transmitting_boundary_history.py
```
