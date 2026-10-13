# Captured continuous reconstruction of the KS solve

The capture owner reuses the exact code object of
`evolve_ks_difference_envelope` in a private namespace whose DOP853 factory
records each accepted step. It does not patch a process global or modify the
authenticated evolution owner. A focused test verifies identical endpoint
fields and preparation identity with capture enabled.

The local interpolation coefficients come from
[SciPy's DOP853 implementation](https://github.com/scipy/scipy/blob/v1.17.1/scipy/integrate/_ivp/rk.py).
For step fraction x, the reconstruction stores both binary endpoints and six
correction coefficients:

\[
\widehat y(x)=(1-x)y_0+xy_1+
 \sum_{j=1}^{6}F_j x^{\lfloor(j+2)/2\rfloor}
 (1-x)^{\lfloor(j+1)/2\rfloor}.
\]

The endpoints anchor the mathematical polynomial exactly, so adjoining
segments are continuous. This avoids treating a rounded stored difference
`y1-y0` as an exact endpoint identity. The polynomial's derivative uses the
actual signed rho step. Its coefficients are numerical inputs to validation;
the polynomial is not declared an exact PDE solution.

The residual sampler reconstructs the Fourier envelope and its first three
axial derivatives, then evaluates the actual Dirac equation and the first
two differentiated residuals. Radius derivatives come from the same w,U
functions. It retains `B_z X` and `2 B_z X_z+B_zz X`, the homogeneous
reference equation, the source-energy generator, and the original column
weights. Source carriers are not FFT-differentiated.

Sampled periodic L2 norms are diagnostic values. The saved polynomials make
continuous residual validation possible, but sampled maxima or sums are
not error bounds and are never passed as certified inputs to the
[residual-to-matter estimate](nsc-ks-residual-error.md).


## Streaming capture and replay

The streaming adapter sends each accepted `TrajectorySegment` to a callback
and retains only the preceding endpoint. The existing evolution code object
runs with a private recording DOP853 factory; endpoint fields, dense
coefficients and preparation identity agree with in-memory capture.
Deterministic NPZ chunks store `rho_nodes`, `state_nodes` and
`dense_corrections` for at most four complete segments. Shared endpoints
match bitwise, chunks carry SHA256 hashes, and the manifest updates atomically
after each complete chunk. Existing run directories are refused, and an
interrupted manifest remains incomplete.

Replay verifies numbering, hashes, shapes and continuity one chunk at a
time; incomplete data are rejected unless explicitly requested. Source,
history and final prepared-array bindings remain the caller's responsibility.
The fine value-only protocol uses N=2048, period $205/512$, `rtol=5e-14`,
`atol=5e-19`, `max_step=1/8192`, and `tangents='zero'`. This is bounded
storage/capture, with unchanged physical interval and source, rather than
a residual enclosure or gate certificate.
