# Current finite-history UV remainder

The v4 coefficient owner proves cancellation of the paired `E^-2` vacuum
coefficient. The v1 remainder successor preserves the owned recurrence through
M=2 and tests the formal `L0` telescope through M=4 and the angular involution
for both N and beta. The first paired order not forced to cancel is `E^-3`.
The formal telescope does not construct the production current-history A3 or
A4 coefficients.

The order-four truncated envelope has defect

\[
(L_0-2iE a^{-2}\Pi_{-s})\sum_{j=0}^4E^{-j}A_j
=E^{-4}L_0A_4.
\]

`A3` and the real part of the major-`A2` transport depend on characteristic
integrals; they are not replaced by local metric jets. For the representative
massless pair, the local right-hand side of the imaginary major-`A2` transport
vanishes, but its upstream datum is not owned, so neither the coefficient nor
its characteristic value is recorded as zero. A conditional H2 majorant is
executable once two inputs exist: the upstream affine remainder above order
four and the current-history integrals of `L0 A4`. The existing radius owner
already supplies the required `B_z` and `B_zz` commutator integrals.

For the smallest original angular pair, the Fermi thermal tail above 160 is
less than about `4.52e-460` in N and `3.78e-460` in beta; above 320 it is of
order `1e-920`. These directed decimal uppers are not rounded to zero. The
vacuum coefficient and remainder still have `C4=None`, `C_M=None`, so the UV
budget component remains OPEN. The `E^-2` cancellation and the `E^-3` first
allowed paired order are not that tail bound.

The v1 pilot JSON
`results/development/nsc-ks-current-uv-remainder-pilot.json` and schema
`NSC-KS-CURRENT-UV-REMAINDER-PILOT-v1` are immutable historical bytes.

## v2 successor: represented quantities and authenticated bindings

The current executable owner is
`scripts/derive_nsc_ks_current_uv_remainder_v2.py`, producing
`results/development/nsc-ks-current-uv-remainder-v2.json`. It extends the
same M=4 remainder module without rewriting the v1 pilot.

Each of `A2`, `A3`, `n4`, `C4`, `C_M`, the current-history `L0 A4` H2
integrals, the upstream higher-order H2 remainder, the transported majors and
the vacuum N/beta tail is named with its exact formula or repository input
and bound to

- state law `C_\Sigma[g]=F[g]C_{\rm src}F[g]^\dagger`, `dC_{\rm src}=0`
  (`nsc_evolved_incoming_state.py`)
- the live profile `0b0e4ced…`
- horizon preparation (`nsc_paired_horizon_preparation.py` and the restart
  source scales)
- the original source inventory and cutoff-bridge ledger
- the owned KS subtraction vertices (`nsc_common_subtracted_ks_source.py`)

The small authenticated control computes only what those bindings already
make rigorous: the conditional formal M=4 telescope, equal-mu `E^-2`
cancellation, `E^-3` as the first allowed paired order, massless vanishing of
the local `Im(T_s` major `A2)` transport right-hand side, owned `B_z`/`B_zz`
integrals, Fermi occupation tails, and a direct finite-energy-window check of
those occupation tails. The manufactured H2 and even-parity examples are
diagnostic fixtures only and cannot supply a current-history `L0 A4`, `C4`, or
`C_M`. The production definition connecting the full-envelope cubic term to
the required history-minus-reference `C4` is still OPEN. Unavailable values
stay `null` and `OPEN`. No frozen-background tail certificate is copied onto
the changed history. No null is written as zero.

```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 /Users/admin/Projects/BlackHoles-Infinity/.venv/validation/bin/python -m pytest -q tests/test_nsc_ks_finite_history_uv_remainder.py tests/test_nsc_ks_current_uv_remainder_v2.py
python3 scripts/derive_nsc_ks_current_uv_remainder_v2.py --check
```
