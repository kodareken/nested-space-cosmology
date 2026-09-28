# The fourth-order reference now follows a supplied general-KS history

`GeneralKSFourthOrderReferenceHistory` consumes the authenticated
`ModeResolvedCauchyState` payload and a supplied homogeneous
Kantowski--Sachs history.  It uses the retained canonical axes

$$
h_{jn}(k;\tau)=\left(-m_j,\frac{\lambda_n}{r(\tau)},
                         \frac{k}{a_\parallel(\tau)}\right)
$$

and constructs the negative-energy Bloch reference through normal-derivative
order four.  The derivative in every recursion step is

$$
D_\perp=N^{-1}\frac{d}{d\tau}.
$$

For $b^{(0)}=-h/|h|$, the already adopted recursion is evaluated node by node,

$$
b^{(n)}=-\frac{h\times D_\perp b^{(n-1)}}{2|h|^2}
-\frac12 b^{(0)}\sum_{i=1}^{n-1}b^{(i)}\cdot b^{(n-i)},
\qquad n=1,\ldots,4.
$$

This is the project-specific application of the existing reference rule; no
heat-kernel or generic ADM result is rederived.

## Four action variations

The owner returns the reference contribution that composes with the canonical
force as

$$
-\mathrm{Tr}(C\,\delta_AH)
+\mathrm{Tr}(C_{\rm ref}\,\delta_AH)
=-\mathrm{Tr}[(C-C_{\rm ref})\delta_AH]
$$

for all four fields $A=(N,\beta,q_{\rm ADM},r)$, with
$q_{\rm ADM}=a_\parallel$.  The per-channel action measure comes from the
serialized `copy_count` and `degeneracy`; the older channel `factor` values mix
several stress-frame conventions and are not promoted into a new physical
stress normalization.

The stored quadrature contains one axial-frequency side.  The implementation
therefore constructs both $k$ and $-k$ references and averages their action
vertices before integration.  In particular,

$$
\frac12\left[-k\,\mathrm{Tr}C_{\rm ref}(k)
+k\,\mathrm{Tr}C_{\rm ref}(-k)\right]=0.
$$

This is the required parity-completed momentum subtraction.  A one-sided
positive-frequency payload is never reported as vacuum momentum.

## Scope and gate

The executable construction covers all 33 serialized channels and 1,904 mode
blocks on any validated supplied homogeneous history.  A manufactured local
stencil checks the seed Hamiltonian axes, the order-by-order Bloch identities,
the parity cancellation and refinement of all four nodal force arrays.  That
stencil is a numerical regression only; it is not a physical duration or a
candidate metric history.

This owner can pass while the extended same-action gate remains open.  It does
not select $g_\star$, transport through a tilted interface, mix frequency
blocks, provide the local induced history, assign a finite stress or null
sign, or update the constraints.  The two remaining composable owners are
`GeneralKSLocalInducedHistory` and the transmitting tilted-Landau interface
with its frequency-mixing Cauchy map.  Coupled evolution remains closed.

The focused verifier is

```sh
python3 scripts/derive_nsc_general_ks_reference.py --check
```

It authenticates and reads the existing mode-state payload.  It does not run
the historical state/source generators or alter $A$, magnetic $q$, $\Omega$,
$\zeta$, or $V_{\rm full}$.
