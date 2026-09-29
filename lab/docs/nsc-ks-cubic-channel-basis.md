# Shared cubic channel basis

This transports one massless unit-angular basis per characteristic sign and
assembles a local effective-source cubic coefficient for a later mass and
angular label. The factorization is already proved by
[`massive_cubic_identities`](nsc-ks-massive-cubic-uv.md) and `_factor_blocks`
in `src/recursive_horizons/nsc_ks_massive_cubic_uv.py`. All 74 residuals in
that owner are the exact string `0`. They are not recomputed here.

The integrator, coordinate jets and directed cell rules are the existing
[massless enclosure](nsc-ks-cubic-uv-enclosure.md). This file does not edit
those owners, their records, or a physical residual assembly.

Owner: `src/recursive_horizons/nsc_ks_cubic_channel_basis.py`.

## What is integrated

For each rho box, incoming $z$ and characteristic sign $s=\pm1$, construct
one massless `CubicGeometry` forcing jet at $\ell=1$. That label is an
algebraic basis, not a new source species. $f$ is the $q_z$ forcing, $g$
the $q_{zz}$ forcing, $c$ the coupling and $h$ the massless current.
$a$ is the owned background axial series at the same rho box.

$$
B=s a^2 c f,\qquad A=h-B.
$$

$A$ has to come from that massless jet. The massive current is
$h+m^2 a^2 f$. Putting that extra piece into $A$ counts $m^2 a^2 f$ in
$J_2$ and again in $J_M$.

The initial differences are zero. On each decreasing cell

$$
\begin{aligned}
Q_z'&=f,\\
Q_{zz}'&=g,\\
J_2'&=A,\\
J_4'&=B+c Q_{zz},\\
J_M'&=a^2 f+s Q_{zz}.
\end{aligned}
$$

Positive Taylor order calls `taylor_volterra_cell` three times on the **same**
pre-cell $(Q_z,Q_{zz})$. Order 0 calls the owner's uniform `volterra_cell`
three times in the same way, because the Taylor cell starts at order 1.
Updated $Q_{zz}$ is not fed from one channel into another inside the cell.

At $\rho=1$, with the unit-angular jets $b_0$,

$$
\begin{aligned}
N_2&=s J_2/a+a^3/2\,(b_{0\rho}^2-b_{0\rho}^{\mathrm{ref}\,2})+s a b_0 b_{0\rho z},\\
N_4&=s J_4/a-4 s b_0^2 Q_z/a,\\
N_M&=s J_M/a-s a Q_z.
\end{aligned}
$$

The reference $Q_z$ difference stays the homogeneous zero representative, so
no extra reference $Q_z$ is subtracted. These surface formulas were checked
with the forcing identities for both signs; this module does not reopen that
proof.

## Assembly

For a finite mass $m$ and a finite nonzero angular label $\ell$,

$$
\begin{aligned}
q_z&=\ell^2 Q_z,\\
q_{zz}&=\ell^2 Q_{zz},\\
J_3&=\ell^2 J_2+\ell^4 J_4+m^2\ell^2 J_M,\\
N_{\mathrm{bracket},3}&=\ell^2 N_2+\ell^4 N_4+m^2\ell^2 N_M.
\end{aligned}
$$

The same powers rebuild the massive forcing jet:
$q_z$ and $q_{zz}$ scale by $\ell^2$, the current is
$\ell^2 A+\ell^4 B+m^2\ell^2 a^2 f$, and the coupling is
$\ell^2 c+s m^2$. Squares of $\ell$ and $m$ make the assembly even in
those signs. Evenness does **not** collapse $s=+1$ with $s=-1$. Both
characteristics are integrated.

The directed API keeps ball enclosures, the full path from the exact flat
upstream $\rho$, and the old rule that an unresolved subdivision emits no
enclosure. A supplied $z$ ball is one box when the geometry accepts it.
`entire_incoming_interval` and `z_box_is_full_incoming_interval` stay false.
`C_M`, the higher UV remainder and uniform $C_4$ on $I$ stay null. The
physical local gate stays OPEN. Flint precision and series cap are restored
on every public return, including failures.

If the CPU deadline passes before the path is finished,
`CubicChannelBudgetExceeded` is raised and no partial coefficient is returned.

## What this run checks

Manufactured scalar and polynomial cells match their integrals, including both
characteristic signs and a non-constant axial series. The group-14 forcing jet
at $\ell=\sqrt5$, $m=\pi/2$ is rebuilt from the unit basis for both
characteristic signs and both angular signs. A current contaminated by
$m^2 a^2 f$ is not that jet.

One full-path control uses the original history, the original upstream rho,
the incoming center, 1024 cells and Taylor order 8, for both characteristic
signs. It is compared with the stored group-14 coefficient enclosure. Overlap
is only a numerical consistency check. The factorization proof remains the
74 exact residuals. No full incoming interval and no other family is run.

## What remains open

The next missing connection is to evaluate this same basis across the incoming
interval, assemble each nonzero-angular family algebraically, then fold the assembled
coefficients with the existing channel ledger. That aggregation is not in
this slice. The ledger's quadrature measure is 320 for groups 10, 11, 12, 31
and 32, and 160 for the other nonzero-angular families. Multiplicity is per
signed family and is not uniform. The energy folding factor is 1. After that,
the higher remainder $C_M$ is still unconstructed. One center box is not a
UV tail, and this coefficient does not move an infinity or a cosmological
constant.

```sh
.venv/validation/bin/python scripts/lab.py -m pytest -q tests/test_nsc_ks_cubic_channel_basis.py
```

Parent integration added explicit rejection of missing and boolean scalar
forcings, channel states and assembly entries. All 11 tests pass. The full-path
center control retains the stored group14 widths; no full-I or higher-remainder
claim is added. Independent implementation review passed all 11 tests and found no actionable defect.

The reviewed center control took about 24.85 CPU seconds for both signs.
Its coefficient widths were 0.999 to 1.000 times the stored group14 widths.
This measures one shared march; no all-family runtime or full-I result is claimed.
