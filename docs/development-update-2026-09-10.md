# Development update: common source and gravitational backreaction

This update publishes the committed laboratory work through
`984a7a68e559fffc367d2884ac60f23fd490beba`. It contains sixteen derivation
notes, twelve compact result records, eleven generators, eight implementation
modules and four focused tests. The 52 scientific files are imported
byte for byte from that Git checkpoint.

The 100-record manifest describes the preceding preprint checkpoint.
This subsequent development snapshot adds the source derivations below.
The September 10 editorial revision presents the paper and overview pages
around these connected equations; the original result records are preserved.

## Read the connected derivation

Start with [the common-source equations](nsc-common-source-derivation.md).
For the declared finite quadratic fermion system, the same parent–child
coherence supplies link energy and regional transfer:

\[
z_B=\operatorname{Tr}_p(B C_{cp}),\qquad
E_{\rm link}=2\Re z_B,\qquad \dot N_p=2\Im z_B.
\]

The variation of that Hamiltonian supplies the geometric force. The linked
results specify what has been constructed and the physical domain of each
statement:

| Connection | Derivation and implementation |
|---|---|
| Boundary response with occupied, empty and initial-correlation kernels | [Finite Gaussian boundary state](nsc-boundary-state.md) |
| Canonical tower, cutoff conversion and warped determinant in one account | [Canonical–spectral bridge](nsc-canonical-spectral-bridge.md) |
| Imported doubled-Dirac anomaly with the actual compact domain | [Compact anomaly application](nsc-compact-anomaly-bridge.md) |
| Exact finite Gaussian representation of the adopted cutoff measure | [Gaussian measure](nsc-gaussian-cutoff-measure.md) |
| Vacuum matching through a quadratic geometric history phase | [Real-time development prescription](nsc-vacuum-matched-ctp.md) |
| Lapse, shift, radius and radial-metric variations, with Ward identities | [ADM source construction](nsc-adm-source-constraints.md) |
| Recorded child-frame tensor transported into the source-equation frame | [Neck source map](nsc-adm-neck-source-map.md) |
| An identified local contribution on the actual neck | [Local warp source](nsc-warp-local-neck-source.md) |
| Recorded outgoing flux and the necessary evolving-parent mass balance | [Backreaction condition](nsc-parent-backreaction-gate.md) |
| Charged field, state and coefficient dictionary into the imported MMP throat | [NSC–MMP embedding](nsc-mmp-embedding.md) |
| Canonical Lorentzian evolution plus one-counted induced source terms | [Causal common functional](nsc-causal-common-functional.md) |
| Identification and quantitative bound for the unresolved absolute vacuum coefficient | [MMP vacuum owner](nsc-mmp-vacuum-owner.md) |
| Recursive zero-tadpole law fixing the homogeneous unlinked vacuum source | [Relational vacuum normalization](nsc-relational-vacuum-normalization.md) |

The known channel has outward power
\(P_0=0.00014222067954246644\) in the recorded throat units.
Under the stated asymptotic Einstein normalization, its contribution obeys
\(dM_B/du=-P_0\). This is a necessary initial mass balance, not a complete
geometry trajectory or a prediction of cosmological energy transfer.

## Exact remaining definition

Canonical CTP now owns causal response. Its finite interface returns all four
metric forces, stress, power, retarded response, noise and Ward residuals.
The embedding gate passes the field, parity, AP-state and normalization maps.
The retained partial coefficients fail the charged seed condition, and the
owner audit finds no existing contribution that determines the required
absolute $V_{\rm full}$. The subsequent recursive law projects only the local
$a_0$ tadpole and gives $V_{\rm full}=\lambda_4=\Xi=0$ while preserving
$A,C$ and finite topology response. No coefficient, incoming flux or mass
history is fitted. The physical return state and complete $A,C$ normalization
are the next embedding conditions.

## Provenance and inspection

[The development snapshot index](../results/development-snapshot.json)
records every imported path, its SHA-256 and byte count, the full source
commit, and each record's generator and original comparison policy.
The original 100-record [release manifest](../results/manifest.json) is
unchanged; the twelve later records are indexed separately.

`make check` authenticates both collections, the stored dependencies,
publication paths and links. This update reuses the completed scientific
verification recorded with the derivations. It does not replay the angular
sum, radiation calculation, historical chain or PDF build.

The listed generators retain their `--check` entrypoints. As with the
earlier records, direct laboratory reproduction needs its authenticated
configuration rather than the curated package metadata. To prepare an
isolated checkout for an explicitly requested replay, copy
[the existing laboratory configuration fixture](../tests/fixtures/source-checkpoints/161028d52ef206a4bc99a25bf97ea467c30e7f16/pyproject.toml)
to that checkout's `pyproject.toml`, then invoke only the desired generator
from the snapshot index. Preserve the public checkout and original JSON.
