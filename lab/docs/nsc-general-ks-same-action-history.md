# The declared homogeneous KS history class has no stationary solution

`GeneralKSSameActionHistoryFunctional` now owns the already declared
variation

$$
\frac{\delta\Gamma_{\rm one}^{\rm CTP}[g,C[g]]}
     {\delta g_\Delta^A(\tau)}=0,
\qquad A\in\{N,\beta,q_{\rm ADM},r\},
$$

for smooth homogeneous Kantowski--Sachs histories.  The blockwise state owner,
Einstein coefficient, magnetic gauge term, relational vacuum and local
allocation are entered in one ledger.  The lapse and shift constraints are
evaluated before an interior collocation because every solution must satisfy
them on the limiting seed slice.

## The shift variation excludes the entire class

For a homogeneous KS history, the geometric momentum-constraint side has no
spatial gradient.  Every already declared local metric invariant and its
fourth-order reference is axial-reflection even.  The static magnetic field has
no normal--axial Poynting component, the compact local contribution has
$T_{01}=0$, and $V_{\rm full}=0$.  No boundary momentum shell is declared.

The complete seed state instead has

$$
T_{01}=0.00122387015580509.
$$

With the locked $A=0.04501936182826115$, the shift variation is

$$
\boxed{
\mathcal E_\beta(0)
=-\frac{T_{01}}{2A}
=-0.0135927088490713.
}
$$

Unitary finite-mode evolution and smooth metric vertices make $T_{01}(\tau)$
continuous.  It therefore remains nonzero on a neighborhood of the seed,
while the homogeneous shift constraint requires it to vanish on every slice.
This contradiction is independent of history duration, radial profile and the
diagonal fourth-order/local terms.

Consequently no smooth stationary history exists in the declared homogeneous
frequency-diagonal, no-interface KS class.  The remaining $q_{\rm ADM}$ and
$r$ interior equations need not be searched after an independent constraint
has failed.

Extending the certificate to the tilted, frequency-mixing interface class
still requires three already named implementations:

1. `GeneralKSFourthOrderReferenceHistory`, including all four ADM variations
   and the parity-completed momentum subtraction;
2. `GeneralKSLocalInducedHistory`, returning nodal action forces rather than
   constant endpoint forces;
3. the transmitting tilted-Landau interface variation and Cauchy map.

Evading the result requires a changed physical scope: a derived same-action
counterflow, spatially inhomogeneous spherical data, a declared momentum shell,
or another initial quantum state.  None is inserted here.  No finite stress or
metric trajectory is fabricated, and coupled evolution remains closed.

The focused verifier is

```sh
python3 scripts/derive_nsc_general_ks_same_action_history.py --check
```

It reads the authenticated state/action ledgers and performs no source,
metric, MMP or historical generator.
