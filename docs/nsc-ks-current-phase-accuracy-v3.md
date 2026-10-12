# CURRENT-HISTORY nodal phase-value enclosure v3

This is the successor of the [v2 current-history phase record](nsc-ks-current-phase-accuracy-v2.md) for the live local history `0b0e4ced…` (`results/development/nsc-ks-gate-history-lm-broyden.json`). The history is n=32 in `w` and `U`. The collocation is the current-family Chebyshev set on `I` (129 or 257). The upstream slice is the original archive binary64 `rho_up=1.0300000000000002`.

It reuses the v2 owner as a library: helpers, context, true-gradient contraction and the production unit-angular 192/384 path. It does not copy that 750-line script and does not assign v2 module globals. The enclosure call is a private wrapper around `directed_source_phase_z_enclosure` with a declared tighter radius

\[
\texttt{target\_radius}=\frac{\texttt{unit\_angular\_target}()}{16}=\frac{2}{10^{16}\cdot 16}.
\]

The divisor `16` and the packed upper of that ball are bound in the schema. Production 192 and 384 Gauss gradients are compared to the **true-gradient ball**, including midpoint displacement, not only the ball radius. The true N,beta factor uses the exact integer 107952 and the enclosed background `a(1)`.

The component target remains `1e-12` per N and beta at each requested node. Between-node phase remainder, source/field error, finite tail and the physical local gate stay `None` / `OPEN`. No source or field evolution is performed. A positive nodal bound is this component only.

The inherited v2 3-node sample left N at about `1.20e-12` because the directed ball radius, not Gauss bias, sat above `1e-12` at the endpoints of `I`. 192 and 384 both sat at the true-ball midpoint to about `1e-16`. This owner tightens that radius; it does not change quadrature order or the history.

The recorded numerical enclosure lives in the immutable result register and its exclusive node files. This note does not restate those measured uppers. Status names the actual failed-node count, not the size of the requested sample. `all_declared_nodes_covered` is true only when every declared collocation node is present and complete.

v2 embeds full node rows in the top-level JSON and has no payload hash descriptor. v3 keeps exclusive node proofs and binds a recursive descriptor: each node file and the binding file as `{path, sha256, bytes}`, and a directory digest that is the hash of those child hashes. `--check` authenticates that tree, replays every requested node, and ignores CPU fields. Checkpoints also bind schema, bits, gauss counts, axial scale, the tighter radius and complete quadrature.

Default `--selection pilot` is nodes `0`, middle and last, with a 60 CPU-second cap. All 129/257 nodes are refused until that pilot is recorded (`--allow-all`). A full 129-node run is justified only when the 3-node actual N,beta errors are both `<=1e-12`, and then uses a 300 CPU-second cap.

```sh
.venv/validation/bin/python scripts/derive_nsc_ks_current_phase_accuracy_v3.py \
  --history results/development/nsc-ks-gate-history-lm-broyden.json \
  --nodes 129 --selection pilot --cpu-budget 60 --output-prefix pilot --record
.venv/validation/bin/python scripts/derive_nsc_ks_current_phase_accuracy_v3.py \
  --history results/development/nsc-ks-gate-history-lm-broyden.json \
  --nodes 129 --selection all --allow-all --cpu-budget 300 --output-prefix n129 --record
.venv/validation/bin/python -m pytest -q tests/test_nsc_current_phase_accuracy_v3.py
```
