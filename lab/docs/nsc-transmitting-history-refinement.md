# Spatial refinement of the conditional transmitting history

The [first spatial-history record](nsc-transmitting-history-modes.md) retains
its OPEN numerical result. Its largest 401-to-801-point response difference
was $1.71\times10^{-7}$, above the declared $3\times10^{-8}$ tolerance.
Temporal refinement and the algebraic controls passed. The tolerance,
metric control, state preparation and physical coefficients are unchanged.

Only the 30 signed families that failed the spatial comparison receive
finer grids. Previously passing families are imported from that record.
The next comparison uses 801 and 1,601 points; 3,201 points are used only
when that new comparison still fails. Temporal refinement is checked on
the accepted spatial grid as well.

The calculation reuses the [locked radial signed map](nsc-massive-signed-preparation.md)
for the opposite-angular family. The source columns are first checked against
that map. In the characteristic frame, with real metric coefficients,

$$
L_{-\lambda}=\sigma_3 L_\lambda^*\sigma_3,
\qquad
w_{-\lambda,E}=\sigma_3 w_{\lambda,-E}^*.
$$

The second equality follows by uniqueness of the same linear Cauchy problem
with the mapped, authenticated input. It does not choose an angular partner
state from an arbitrary covariance completion. The old, independently
computed nonlinear mode projections also check the map before it is used to
avoid a duplicate negative-angular solve. The projected response is conjugated
with both energy indices reversed.

This refinement resolves a spatial numerical error in a conditional field
propagator. It does not turn the supplied pulse into a stationary metric,
establish the complete spectral covariance integral, or populate physical
endpoint jets. The earlier failed record remains intact, with its original
hashes and comparison policy. No new source sector or action term is added.
