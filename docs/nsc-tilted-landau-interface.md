# The tilted transmitting interface admits frequency mixing but is not selected

This owner applies the existing spatial boundary domain and the serialized
`ModeResolvedCauchyState` to the source-selected Landau endpoint.  It adds no
bulk term, boundary shell, duration, metric history or state refit.

The imported throat domain identifies parent and child spinors in their common
coordinate frame.  The parent normal is opposite to the child normal, so the
two boundary forms cancel.  The existing boundary calculation evaluates this
on a fixed spatial interval and explicitly does not provide a stationary PG
map across the horizon.  The imported curvature EFT likewise leaves Euler,
$\Box R$ and induced boundary variations with their existing separate owners.

## Weighted interface variation

For one retained channel with quadrature weights $w_i$, the seed and tilted
coordinate Cauchy metrics are

$$
G_0=\mathrm{diag}(w_i)\otimes I_2,
\qquad
G_L=\mathrm{diag}(w_i)\otimes B_L,
\qquad
B_L=e^{-\eta\sigma_2}.
$$

The transmitting interface variation is the finite operator residual

$$
\boxed{R_I[K]=K^\dagger G_LK-G_0.}
$$

Since both metrics are positive, every solution and only every solution is

$$
\boxed{K_c=G_{L,c}^{-1/2}V_cG_{0,c}^{1/2},\qquad
       V_c\in U(2n_c).}
$$

The canonical half-density Cauchy map is therefore

$$
U_c=G_{L,c}^{1/2}K_cG_{0,c}^{-1/2}=V_c.
$$

Thus $B_L^{1/2}$ enters the conversion between coordinate traces and canonical
data; it is not used as $U$.  Arbitrary dense $V_c$ may mix the retained
frequency nodes inside a channel.  Different angular and compact channels are
not mixed because no existing boundary owner changes those symmetry labels.

## Complete retained family and witnesses

The payload has 13 channel blocks of dimension 96 and 20 of dimension 128.
The complete channel-block map has dimension and rank 3,808.  The weighted
isometry equations have 894,976 real kernel variables, Jacobian rank 447,488,
and linearized nullity 447,488.  The same number is the real dimension of
$\prod_c U(2n_c)$, so the endpoint equation cannot select a unique kernel.

The verifier evaluates two members in every one of the 33 blocks:

1. $V_c=I$, the frequency-diagonal control;
2. $V_c=F_{n_c}\otimes I_2$, a dense discrete-Fourier witness that mixes all
   retained frequency nodes while preserving the channel label.

For each member it checks $R_I$, the reconstructed canonical map, unitarity,
full rank, covariance Hermiticity, the CAR eigenvalue interval and transport
of the complete covariance spectrum.  The content-addressed artifact stores
every dense mixed $U_c$ and its coordinate kernel $K_c$; it does not store a
claim that either witness is physical.

## Gate status

The construction and residual gate pass, including a true k-block-mixing CAR
witness.  Physical selection remains **OPEN**.  The fixed spatial boundary
response and the Landau endpoint provide neither the intervening hypersurface
embedding nor a mode-resolved scattering kernel.  The exact next missing owner
is the same-action mode-resolved transmitting boundary scattering/embedding
operator that selects one $V_c$ in every retained channel.

No finite stress, interface $T_{01}$, null signs or updated constraints are
assigned.  The extended existence/non-existence gate must compose this OPEN
selector status with the fourth-order-reference and node-wise-local owners.
Coupled metric evolution remains closed.

The focused verifier is

```sh
python3 scripts/derive_nsc_tilted_landau_interface.py --check
```

It authenticates the existing mode-state, constraint, boundary and curvature
records, reconstructs every channel block, and compares all record fields and
artifact arrays.  It performs no source generator, duration solve or metric
timestep.
