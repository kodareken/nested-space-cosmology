# The paired unweighted reference-band exchange has zero bulk Euler force

For the retained canonical subtraction, through formal order four,

$$
\boxed{\left[\delta S_{\rm band}
-\tfrac12\operatorname{Tr}(P\star\delta H+\delta H\star P)\right]_{\rm bulk}=0.}
$$

This result uses the complete **unweighted** $k\in\mathbb R$ trace, the
actual equal angular-sign weights, and independent compact variations of
$N,\beta,a,r$ on
smooth finite jets with positive nondegenerate $N,a,r$. It closes the bulk
exchange in the [existing reference action](nsc-reference-band-action.md).
It does not delete the finite momentum-cutoff terms or complete physical
time/spatial endpoint variations.

## Reuse the owned action identity

Let $D=\delta U\star U^\dagger$ and $P=U^\dagger\star\Pi\star U$.
The [common source owner](../src/recursive_horizons/nsc_common_subtracted_ks_source.py)
already separates the band variation from the symmetric projector vertex.
Using its existing formal identity, that difference is

$$
i\varepsilon\partial_T\operatorname{tr}(\Pi D)
+\operatorname{tr}[\Pi U,\delta H\star U^\dagger]_\star
+\operatorname{tr}[\Pi D,h]_\star
+\tfrac12\operatorname{tr}[\delta H,P]_\star.
$$

No new reference state or counterfunctional is introduced. Under compact
test variations the time derivative is a boundary term. For the spatial
Weyl product, integration by parts gives

$$
\int dz\,\operatorname{tr}[A,B]_\star
=2\sum_{\ell\ {\rm odd}}\frac{(i\varepsilon/2)^\ell}{\ell!}
\partial_k^\ell\int dz\,
\operatorname{tr}[(\partial_z^\ell A)B].
$$

Here the sum is over odd $\ell$; the executable certificate checks the
Leibniz coefficients for every order one through four. The momentum-boundary
primitive contains $\ell-1$ momentum derivatives. Spatial endpoints are
retained separately when variations are not compact.

## The only nondecaying primitive cancels in the actual angular pair

Write $K=|k|$, $s=\operatorname{sgn}k$, $u=ma$, $v=\lambda a/r$ and
$c=N/a$. For $m>0$ choose $\Pi=(I+S_1)/2$. Its overlap with the lower
normal-energy projector is at least $1/2$ for every $k$. The owned polar
frame expands as

$$
U_0=\frac{I+isS_2}{\sqrt2}
+\frac{\tfrac u2(I-isS_2)-ivS_3}{\sqrt2K}+O(K^{-2}),
\qquad
\operatorname{tr}(\Pi D_0)=\frac{i\delta v}{2k}+O(K^{-2}).
$$

The two possibly nonzero constants in the first-Moyal $k$ primitive are

$$
\frac12(s\delta c+\delta\beta)\partial_zv,
\qquad
\frac12(sc+\beta)\partial_z\delta v.
$$

They are nonzero in general for an individual angular sign. Both are odd
under $\lambda\mapsto-\lambda$. The retained inventory supplies 30 equal
angular pairs and three zero-angular families, giving 63 signed families
from 33 groups. Each sign receives `copy_count * degeneracy / n_signs`,
as in the existing source assembly. Thus the constants cancel exactly in
the paired trace; no physical occupation symmetry is assumed.

For $m>0,\lambda=0$, both constants vanish individually. For
$m=0,\lambda\ne0$, choose $\Pi=(I-\operatorname{sgn}(\lambda)S_2)/2$.
The leading frame is a planar rotation with zero projected connection.
The exact chiral LLL has $U=I$ in each $k$-sign chart and $\delta U=0$.
For this sector the full-line trace means the sum of the two owned fixed
nonzero-$k$ chiral charts. Its zero exchange is a chartwise statement;
this certificate does not differentiate a single discontinuous
$\operatorname{sgn}(k)$ projector through the gap-zero point. The separately
owned $c=4$ geometric allocation remains required and is counted once.
No gap floor or change of the separately owned light geometric allocation
is used. The remaining symmetric-projector commutator has no constant
primitive; its leading growing term is a spatial total derivative.

## Higher formal orders decay, including the shift

The actual Hamiltonian and its variation are at most linear in $k$:

$$
H=-\beta kI+N(-mS_1+\lambda S_2/r+kS_3/a).
$$

The existing projector recursion gives
$F_j=O(K^{-j})$, $G_j=O(K^{-j-2})$ and $P_j=O(K^{-j-1})$ for
$j=1,\ldots,4$. Its bounded polar normalization preserves
$U_j,\delta U_j,D_j=O(K^{-j-1})$, while $h_j=O(K^{-j})$.
Coordinate derivatives and smooth metric variations preserve these powers;
momentum derivatives improve them. The scalar shift cancels in the ordinary
spin commutator and remains included in every Moyal power bound.

Consequently the momentum-boundary bounds at formal orders one through four
are $O(1),O(K^{-1}),O(K^{-2}),O(K^{-3})$. The sole constant term was
evaluated and paired above. The [new owner](../src/recursive_horizons/nsc_reference_band_bulk.py)
checks the leading formulas by exact Pauli/SymPy coefficient operations and
enumerates the integer order bounds of the reused recursion. There is no
large-momentum scan or rerun of a physical generator.

The [record](../results/development/nsc-reference-band-bulk.json) reports
exact residuals and the scoped zero bulk gradient. Finite-$k$ endpoint
objects, physical time/spatial endpoint objects, spectral convergence and
stationarity remain separate. An extra momentum weight would leave its
derivatives in the trace and is outside this result. This is not a
raw-heat equivalence condition and assigns no value to $\Gamma_{\rm rest}$.

The LLL chart scope was made explicit after the initial certificate was
authenticated. Commit `e22b03903d58af44043b5d9b5a34901615cd5c48` preserves
the original record bytes. The clarified record identifies that predecessor
and its SHA-256; all exact coefficient results and zero residuals are unchanged.

```sh
python3 scripts/derive_nsc_reference_band_bulk.py --check
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_reference_band_bulk.py
```
