# The retained Gaussian mode state is now serialized

The previous source records retained integrated stresses but discarded the
mode matrices needed by unitary evolution.  This owner performs one authorized
base extraction from the unchanged kernels and stores every retained physical
Gaussian covariance in a deterministic, content-addressed NPZ.

The state contains 33 channels and 1,904 physical $2\times2$ blocks:

| Sector | Channels | Nodes per channel | Covariance blocks |
|---|---:|---:|---:|
| LLL, $(j,n)=(0,0)$ | 1 | 48 | 48 |
| Massive angular, $j=0$, $n=1,\ldots,12$ | 12 | 48 | 576 |
| Positive compact, $j=1,2$, $n=0,\ldots,9$ | 20 | 64 | 1,280 |

Identical compact copies and angular degeneracies remain multiplicities; their
matrices are not duplicated.  The canonical basis is

$$
\chi=r\sqrt{a_\parallel}\,\psi,
$$

with channel, frequency and then upper/lower spinor ordering.  In this basis

$$
H_{jn}(k;T)
=-m_j\sigma_1+\frac{\lambda_n}{r(T)}\sigma_2
 +\frac{k}{a_\parallel(T)}\sigma_3.
$$

The LLL covariance is reconstructed on the same 48-node quadrature from its
existing transparent occupation law.  Its two analytic stress constants are
reproduced within $2.8\times10^{-13}$.  The state-independent conformal term
remains a separate, one-counted local allocation.

For massive angular channels, the payload stores the physical covariance and
the existing $E_2+E_4/P_2+P_4$ subtraction arrays.  For positive compact
channels it stores both the physical covariance and the fourth-order
superadiabatic seed reference.  The locked Wilsonian compact tensor remains a
separate local contribution added exactly once.

Reloading the payload with `allow_pickle=False` verifies Hermiticity, the CAR
eigenvalue interval, channel order, offsets, dtypes, array hashes and the full
payload digest.  Reducing the stored blocks and the declared local allocation
reconstructs the old unit-radius tensor within its inherited tolerance.

## The Landau endpoint is not yet a Cauchy map

The recorded stress rapidity fixes the local current metric

$$
B_L=e^{-\eta\sigma_2},
\qquad B_L^{1/2}=e^{-\eta\sigma_2/2},
$$

but $B_L^{1/2}$ is not unitary in the old Cauchy norm.  Applying it directly as
$C\mapsto B_L^{1/2}CB_L^{1/2}$ can move a covariance eigenvalue above one.
The CAR-preserving map must instead be the full Dirac Cauchy isometry

$$
\boxed{
U_{L0}=J_L\,\mathrm{Res}_{\Sigma_L}E_g
       \mathrm{Res}_{\Sigma_0}^{-1}J_0^{-1},
\qquad U_{L0}^\dagger U_{L0}=I.
}
$$

It needs the Landau hypersurface embedding and the intervening metric, gauge,
Dirac and boundary history.  Endpoint $(v,\eta,r_\star)$ does not determine
that history or its frequency-mixing kernel.

The source-selected finite stress also needs the existing fourth-order
reference and local induced allocation generalized from the Bronnikov profile
to the same Kantowski--Sachs history.  The serialized seed reference alone is
not an instantaneous replacement on $r_\star$.

Consequently this record passes state serialization, CAR, and old-surface
moment reconstruction.  It retains the already verified nonzero direction

$$
\left.\partial_{\log r}\rho\right|_C=-0.460687999769215,
$$

but does not label the linearized displacement as the finite $r_\star$ stress.
The full user gate remains open for the physical Cauchy isometry and the
general-history reference/local providers.  Metric evolution remains stopped.

The fast verifier is

```sh
python3 scripts/derive_nsc_mode_resolved_cauchy_state.py --check
```

It authenticates and reduces the stored payload without rerunning a source.
The original extraction was one base run that produced the state artifact and
legacy reductions from the same in-memory vectors; no standalone historical
generator was launched.
