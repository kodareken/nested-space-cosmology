# Streaming capture of the KS difference trajectory

The in-memory capture in [nsc-ks-trajectory.md](nsc-ks-trajectory.md) retains
every accepted reconstruction polynomial. This owner streams the same
observations: each accepted DOP853 step becomes one `TrajectorySegment` and
is handed to a callback. The streamer keeps only the previous endpoint for
continuity. It does not grow a segment list.

The evolution code object of `evolve_ks_difference_envelope` still runs in a
private namespace whose solver factory is the recording DOP853 subclass. The
process-global factory is not patched. Endpoint fields, dense coefficients
and the prepared-state identity are those of the existing capture.

A chunk writer may persist the stream. Each NPZ file contains only
`rho_nodes`, `state_nodes` and `dense_corrections` for at most four complete
segments. Shared endpoints of adjoining chunks compare bitwise. Bytes are
deterministic; each file is hashed with SHA256; the manifest is replaced
atomically after every complete chunk. An existing run directory is refused.
An interrupted stream may leave those complete chunks, but the manifest
cannot be labeled complete.

Source and history bindings, and the final prepared arrays, are not stored
here. The caller that already owns those artifacts remains responsible for
them. Replay loads one verified chunk, or iterates chunks, without assembling
the full trajectory. The reader checks declared hashes, shapes, numbering
and endpoint continuity, and rejects an incomplete run unless that check is
explicitly relaxed.

This is a memory-bound recording and replay adapter. It is not a residual
enclosure, a source-energy certificate, or a local physical gate. No clipping
is applied to the field or to the inherited source covariance and weights.

The intended later control, not executed by this owner, is N=2048 with exact
dyadic numerical length `205/512`. Physical `I` and the usual axial support
are unchanged; `rho_up` is taken from the existing source; the solver options
are `rtol=5e-14`, `atol=5e-19`, `max_step=1/8192`, `tangents='zero'`. Parent
runs that one fine control after review.

```sh
PYTHONPATH=src python3 -m pytest -q tests/test_nsc_ks_streaming_trajectory.py
python3 scripts/check_nsc_ks_streaming_trajectory.py --check
```
