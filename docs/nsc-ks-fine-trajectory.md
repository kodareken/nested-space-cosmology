# Fine continuous trajectory from the same upstream source

One N=2048 value-only trajectory reuses the original `14_1` source
covariance, weights, energy labels and canonical upstream columns. The
history remains the amplitude .001, U=0 control. The unchanged preparation
digest is checked after evolution. The full retarded-variation owner remains
available; omitting tangent storage in this value-error capture does not
freeze the incoming state.

Use `rtol=5e-14`, `atol=5e-19`, and maximum rho step `1/8192`. The finer
representation is motivated by the too-large continuous residual majorant
of the 64-node control. It is one prescribed run, capped at 180 CPU seconds.

The numerical period is `205/512=0.400390625`, keeping the same centre,
physical I and compact axial support. For the declared power-of-two grids
this produces exactly uniform binary64 nodes; the nominal .4 period had
resolution-dependent rounded spacing. This is computational padding, not
a physical period, duration or source change. Its Fourier profile enclosure
must use the same numerical period before residual validation.

Every accepted dense segment is streamed to deterministic chunks of at
most four segments. Source/history bindings and final prepared fields are
stored separately. Hashes, declared endpoints and bitwise shared states
are verified one chunk at a time. A failed or interrupted run remains
incomplete; the producer refuses an automatic repeat at an existing path.

Analytic profile identities include amplitudes, all basis coefficients and
binary coordinate maps. The full history identity additionally includes
normal windows. It is independent of the numerical grid, so a changed
profile cannot reuse an amplitude-only key.

Capture completion proves data availability and replay, not a numerical
field error or physical local gate. Continuous residual validation, source
accuracy and both incoming constraint residuals remain required.
