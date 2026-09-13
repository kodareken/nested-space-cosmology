# The remaining boundary derivative retains its own OPEN slot

Step K uses the same eight-coefficient endpoint object as v0.19.0. The
[record](../results/development/nsc-transmitting-boundary-remainder.json)
preserves its exact extraction residual zero and the diagnostic Weyl maximum
$93.54264532195464$.

The new [Hamiltonian jet calculation](nsc-dirac-endpoint-jets.md) contracts
the known bulk derivative through the existing Gaussian CTP machinery on a
frozen control. That control covector is not an allocated transmitting
boundary contribution. It therefore cannot be subtracted from the Weyl
vector to manufacture a two-sided match.

For the actual transmitting surface, the remaining action still must provide
its induced-metric and normal-jet variations on the common endpoint basis.
The compact reflecting heat-boundary term has a different domain. The
spatial curvature EFT does not supply the tilted temporal boundary
completion. These imported domain distinctions remain intact.

**Step K OPEN:** the transmitting boundary remainder, total physical jets and
two-sided mismatch are unevaluated (`null`); the target matching tolerance is
$3\times10^{-11}$. The independent normal-jet slots remain unavailable. No
existing coefficient is changed and no new boundary term is introduced.
