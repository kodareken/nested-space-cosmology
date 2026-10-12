# First continuous residual-cell certificate

This record encloses the actual reconstructed KS residual on one entire
accepted rho step and the entire numerical spatial period. It uses the
saved nonzero-history fields, not a new field run or a prescribed matter
dictionary. The residual and its first axial derivative are both bounded.

The ingredients are the anchored degree-seven trajectory polynomial,
degree-eight time coefficient enclosures, the certified Fourier balls of
the compact axial profiles, and the reciprocal-radius remainder. Polynomial
products are formed by ball-arithmetic Fourier convolution. Bernstein
convexity bounds every point of the time step. Summed coefficient absolute
values bound the spatial supremum; the source-column l1 majorant is a
conservative upper bound for the row l2 norm. The omitted profile bands,
time Taylor remainders and reciprocal-radius remainder are added explicitly.

For `S=X_xi-h L_g X`, a uniform row bound on S or S_z also bounds its
norm integral over xi in [0,1], equal to the corresponding absolute-rho
residual integral. No residual sample maximum is used as an upper bound.

The selected saved segment is index 30 of 61. This is a continuous
residual certificate for that cell, not a completed history-error bound,
constraint root, source-energy UV certificate or local physical gate. The
other time cells, preparation errors and source coverage remain required.
`--check` verifies bindings and the directed error composition without
repeating the convolution; `--recompute` reconstructs the bound from the
saved fields and coefficient balls. The single certificate preparation is
capped at 90 CPU seconds.
