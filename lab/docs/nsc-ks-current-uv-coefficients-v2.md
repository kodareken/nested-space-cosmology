# Finite-history vacuum-envelope coefficients (current proof: v4)

The current reproducible proof is
`scripts/derive_nsc_ks_current_uv_coefficient_control_v4.py`, producing
`results/development/nsc-ks-current-uv-coefficients-v4.json`.
The v2/v3 development records remain historical. Their recorded
working-tree HEAD did not contain the then-untracked source, and the
old controller's `--check` also allowed an absent record. The v4 owner
reconstructs the equations and pilot from their scientific inputs,
binds the implementation dependencies, and rejects a missing proof record.

v2 constructed the envelope jets and the geometric \(h_{s,z}\) integrand.
Its hash-bound record
`results/development/nsc-ks-current-uv-coefficients-v2.json` is frozen
(`sha256=de8e67e4cd89eb0331aff4e18d577220b83ddca94c0c6bac6d0c5a0fe32585df`)
and is not rewritten. The current construction keeps those
jets, fixes the interval-point comparison, and reduces the leading
coincidence \(E^{-2}\) current without solving all of \(A_3\).

The historical context commit is laboratory `HEAD` at the v3 recording,
together with that v2 SHA. No existing value/runtime owner is changed.
Numerical \(C_M\), the \(e^{-3}\) paired constant \(C_4\), the thermal
remainder, upstream error and the local incoming gate stay `OPEN`/`None`.

## Interval jets versus a binary64 evaluation

A finite-precision `LocalAxialFunction` sample is not the analytic
polynomial. It need not lie in an 80-bit ball of that polynomial. The v2
helper also omitted the axial cutoff \(\chi\). v3 evaluates
`local_axial_interval_series` (cutoff times Chebyshev, with the corrected
`316fca2d` domain guard) and records a directed reconstruction error of the
binary64 value to that enclosure. Containment of the float is diagnostic
only. Tests cover interior, ramp and exterior; the polynomial-only jet is
kept as a contrast on the ramp and outside the cutoff.

## Recurrence, continuity and the Sigma \(E^{-2}\) current

Write \(\ell_{\mathrm{min}}=A_{1,-s}\), \(b=A_{2,-s}\),
\(A_{1,s}=iq\), \(A_{2,s}=R+ih\), and \(T_s=\partial_\rho-s a^{-2}\partial_z\).
Independent substitution of the owned recurrence gives

\[
\operatorname{Re}(\overline{\ell_{\mathrm{min}}}\,b)
=\frac{a^2}{2}\operatorname{Im}(\overline{\ell_{\mathrm{min}}}\,T_{-s}\ell_{\mathrm{min}})
=-\frac{a^2}{2}T_s h.
\]

Exact point-density continuity of the occupied envelope is
\(\partial_\rho n-a^{-2}\partial_z j=0\) with \(j=s(n-2|\mathrm{minor}|^2)\).
The \(e^{-3}\) coefficient of that identity is
\(T_s n_3+(4s/a^2)\partial_z\operatorname{Re}(\overline{\ell_{\mathrm{min}}}\,b)\).
Substituting the recurrence yields \(T_s n_3=2s\,T_s h_z\). A normalized,
\(z\)-homogeneous upstream condition \(n_3=h_z=0\) then implies
\(n_3=2s h_z\) on the slab, including at \(\Sigma\), and the same for
\(g\) minus the reference.

At \(\Sigma\), \(r_z=0\), so \(\operatorname{Im}(\overline{\ell_{\mathrm{min}}}\,\partial_z\ell_{\mathrm{min}})=0\).
The identity-momentum coincidence coefficient is then
\(h_z+\operatorname{Im}(\overline{\ell_{\mathrm{min}}}\,\partial_z\ell_{\mathrm{min}})-s n_3=-h_z\).
The \(S_3\)-momentum coefficient is \(-s h_z+4\operatorname{Re}(\overline{\ell_{\mathrm{min}}}\,b)\).
The \(q\)-channel in \(\operatorname{tr}(V\rho_2)\) vanishes. The remaining
density term cancels \(4\operatorname{Re}(\overline{\ell_{\mathrm{min}}}\,b)/a\)
against the \(S_3\) vertex, including the local
\(\delta r_\rho\) piece
\(d=s a^3\ell\,\delta(r_\rho)/(4r^2)\) with
\(\delta\operatorname{Re}(\overline{\ell_{\mathrm{min}}}\,b)=am d/2\).
The owned vertices, checked against `raw_ks_vertex_coefficients`, are
\(V=-m S_1+\ell S_2/r\) and \(S_3/a\) for \(N\), and \(0\), \(-I\) for
\(\beta\), with the action minus sign once and measure \(de/(2\pi)\).
The action coefficients are therefore

\[
N_2=\frac{\mu s\,\delta(h_z)}{2\pi a},\qquad
\beta_2=-\frac{\mu\,\delta(h_z)}{2\pi}.
\]

\(A_3\) is not solved componentwise. It enters \(n_3\), which continuity
and the upstream condition replace. Coincidence \(J\) already includes
\(-se\,n\), so that density/carrier factor is the \(n_3\) term just
eliminated. The expansion is the occupied vacuum envelope. Thermal
occupations are not set to one; their exponentially small UV remainder
stays unbounded.

## Angular pairs

\(T_s h\) and therefore \(h_z\) are odd in \(\ell\). The original inventory
has 30 nonzero-angular pairs with equal multiplicity and the same energy
grid and weights. With those pairs the leading \(1/\Lambda\) vacuum
coincidence coefficient cancels. A negative control with unequal pair
weights does not cancel. This removes only that algebraic leading
coefficient. It does not bound \(C_4\), \(C_M\), the thermal remainder or
upstream error.

## Geometry indicator

The v2 three-node \(h_{s,z}\) quadrature on `0b0e4ced…` is reused as an
indicator, not a tail bound. Certified quadrature error remains `None`.

```sh
.venv/validation/bin/python -m pytest -q tests/test_nsc_ks_finite_history_uv_coefficients.py
python3 scripts/derive_nsc_ks_current_uv_coefficient_control_v4.py --check
```
