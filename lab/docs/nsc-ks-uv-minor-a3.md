# Minor A3 jet bound

This is the Stage-2 UV primitive named by the current-history remainder v3:
minor \(A_3=a^2/(2i)\,\Pi_{-s}L_0 A_2\) on the live profile `0b0e4ced…` and
the smallest original angular pair (group 1, massless, \(|\ell|=\sqrt{5}\)).
The action, preparation, source inventory and KS subtraction are unchanged.
The local gate stays OPEN. This bound does not contradict \(\Lambda\)CDM and
does not prove eternity. Historical UV remainder v1/v2/v3 records are not
rewritten.

Owner: `src/recursive_horizons/nsc_ks_uv_minor_a3.py`. Record:
`results/development/nsc-ks-uv-minor-a3-v1.json`.

## Parent formulae

On the actual Dirac operator \(L_0=\partial_\rho-S_3 a^{-2}\partial_z-iV/a\)
with \(V=-m S_1+\ell S_2/r\) and \(T_s=\partial_\rho-s a^{-2}\partial_z\),
the vacant recurrence is

\[
A_{j+1,-s}=\frac{a^2}{2i}T_{-s}A_{j,-s}-\frac{a}{2}V_{-s,s}A_{j,s},
\qquad V_{-s,s}=-m+is\ell/r.
\]

Independent substitution of the constructed jets gives the parent minor \(A_2\)

\[
A_{2,-s}=\frac{s\ell a}{4r^2}\bigl(-a_\rho a r+a^2 r_\rho+s r_z+2qr\bigr)
-\frac{im a}{4}(a a_\rho-2q)
\]

and

\[
A_{3,-s}=\frac{a^2}{2i}T_{-s}A_{2,-s}-\frac{a}{2}V_{-s,s}(R+ih).
\]

The first (jet) piece depends on second geometry jets of \(r\) and on first
jets of the transported phase \(q\). It does not depend on major \(A_2\).
The second piece is algebraic in the transported major \(A_2=R+ih\). Dummy
\(L_0 A_j\) symbols are not used. Minor \(A_3\) is not a local metric jet.

For the representative massless pair,

\[
\operatorname{Re}A_{3,-s}=\frac{s a\ell h}{2r},\qquad
\operatorname{Im}A_{3,-s}=-\frac{a^2}{2}T_{-s}A_{2,-s}-\frac{s a\ell R}{2r}.
\]

On-shell \(T_s q=-(m^2+\ell^2/r^2)/2\) and \([T_s,\partial_z]=0\) give
\(T_s q_z=\ell^2 r_z/r^3\). Upstream major \(A_1=0\) is \(z\)-homogeneous, so
the reference \(q_z\) vanishes and \(|q_z|\) is at most duration times that
right-hand side.

## Interval majorants

Second \(r\) jets reuse the owned Chebyshev \(z\)-bounds of \((w,U)\) and the
normal-window constants \(|\chi'|\le 8/\mathrm{width}\),
\(|\chi''|\le 105/\mathrm{width}^2\). Axial \(a_{\rho\rho}\) is bounded on
\([1,33/32]\) from the same \(u=1-\rho\theta\) estimates as \(|a_\rho|\).
Conversion of \(s\)-jets to \((\rho,z)\) uses \(s_\rho=-1/a\). Evaluation
roundoff is excluded. Samples and fitted decay are not proofs.

The jet piece \(\frac{a^2}{2i}T_{-s}(\mathrm{minor}\,A_2)\) and minor \(A_2\)
itself have interval majorants on each history and for the difference. The
generated real-major-\(A_2\) increment \(I=\int T_s R\) is already owned;
its algebraic coupling is bounded. Leading paired \(E^{-2}\) cancellation
is preserved and is not this coefficient. History-minus-reference is kept
distinct from the full envelope.

The difference-coupling bound is restricted to massless channels. Massive
channels have additional mass and transported-imaginary-major terms; the
current API rejects them instead of applying the massless formula. The symbolic
recurrence identities and geometry jet inequalities retain their stated generality.

## Missing input for a complete \(|\delta A_{3,-s}|\)

Both histories share the unowned upstream major \(A_2\) and the same
\(T_s\)-characteristics (they depend only on \(a(\rho)\)). For \(m=0\),
\(T_s h=0\), so \(\delta h=0\) and \(h=h_{\mathrm{up}}\). Then

\[
\delta\Bigl(\frac{R+ih}{r}\Bigr)
=\frac{\delta R}{r_g}-\frac{\delta r\,I_{\mathrm{ref}}}{r_g r_{\mathrm{ref}}}
-\frac{\delta r}{r_g r_{\mathrm{ref}}}(R_{\mathrm{up}}+ih_{\mathrm{up}}).
\]

The first two terms are bounded. The last term is the shared upstream
major \(A_2\) times \(\delta(1/r)\). Massless vanishing of
\(\operatorname{Im}(T_s\) major \(A_2)\) does not set that datum to zero.
The report stores the proved kernel of that term and leaves the complete
history-minus-reference minor \(A_3\) as `null`. Zero deformation makes the
kernel vanish, so that difference is exactly zero without inventing the
datum.

Numerical \(C_4\), \(C_M\) and the vacuum N/beta tail remain `null`. Gate 2
stays OPEN.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_uv_minor_a3.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_uv_minor_a3.py --check
```
