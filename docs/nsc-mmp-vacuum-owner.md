# The MMP embedding stops at the absolute vacuum coefficient

The NSC–MMP dictionary passes for the charged field content, compact parity,
antiperiodic state structure and Einstein/Maxwell normalization. The retained
partial coefficients do not satisfy the charged seed condition. This note
performs the plan's one owner-specific follow-up without rerunning any throat
or spectrum calculation.

At matching cutoff \(\nu=1\),

\[
(V_D,A_D,C_D)
=(0.17254086329081558,\;0.005502284451619701,\;0.004837250697696715).
\]

The seed condition is

\[
V_{\rm full}\leq\frac{2A_{\rm full}^2}{q^2C_{\rm full}}.
\]

Holding the already retained \(A_D,C_D\) only to display the size of the
known mismatch, the \(|q|=1\) upper value is

\[
V_{\rm full}\leq0.01251749643695399
=0.072548010936145\,V_D.
\]

Thus at least 92.7452% of the retained partial volume coefficient would need
an oppositely signed same-action contribution even for the weaker charged
RN–de Sitter seed bound. The imported asymptotically flat MMP branch requires
the complete effective bulk coefficient to vanish.

## Existing terms cannot be reassigned to that job

- Compact heat-boundary terms are supported on the boundary rather than the
  homogeneous four-volume.
- Weyl, Euler, box-R and higher-curvature terms vanish or remain
  curvature-dependent in the flat asymptotic region.
- The finite Casimir term depends on topology and return length.
- Equal-history normalization of the canonical CTP state fixes causal state
  response, not the state-independent local vacuum coefficient.
- The recursive tail could contribute only after its physical return domain
  and stationary solution exist; using it here would assume the downstream
  answer.

The selected canonical-CTP realization is therefore the plan's single
controlled revision: it resolves causal ownership and keeps induced spectral
forces once. It leaves \(V_{\rm full}\) without a same-action ultraviolet
matching condition.

## Gate decision

The current realization does not provide a parameter-free MMP embedding.
No coefficient is fitted, no new field is added, and the recursive scale solve
is not started with an undefined local action. The next admissible input is an
explicit microscopic spectrum/measure completing the supertrace or a derived
recursive normalization equation that determines \(V_{\rm full}\).

```sh
python3 scripts/derive_nsc_vacuum_owner.py --check
```

This reproduces the owner classification and numerical bound from stored
records only.
