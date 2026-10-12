# Retained source columns for the KS history evaluator

`RetainedSourceInventory` reads the authenticated incoming panel inventory
and applies its existing `selected_pieces` masks. Fine panels replace their
coarse regions. The reader preserves original quadrature weights, source
coherences, compact masses and signed angular channels.

For low real panels it uses the accepted `mode_at_one` fields from the
spectral-mode recovery payload. For middle panels it calls
`archived_middle_boundary_modes`, the original order-16 boundary recipe
used by `incoming_high_energy_source`. It therefore uses the same middle
approximation as the incoming baseline, including its zero first column.
It does not replace that convention with a packet inverse. The compact
subgap entries are the eight inherited physical real-energy rows at rho=1;
their unresolved accuracy is retained explicitly.

The existing `restrict_resolved_modes` converts each PG field to canonical
KS amplitudes. The corresponding `ReferenceSourcePanel` can be partitioned
into complete three-column energy fibers, preserving coherence and weights.
For a positive-energy panel with angular label ell, its negative-frequency
partner has angular label -ell, amplitudes `S3 conjugate(A)`, and the owned
negative-frequency source law `I-conjugate(C_src)`. This reuses
`PairedHorizonSeedMap.at_radius`; it applies no additional factor of two.

`fixed_upstream` continues these reference amplitudes from rho=1 to a fixed
rho_up>=1.03 through the existing homogeneous `ks_generator`. This operation
is performed once per input batch, before varying a history. Subsequent
histories receive the same immutable upstream preparation and evolve it
forward. The saved rho=1 reference is never imposed on their incoming state.
The continuation has its own numerical error requirement.

The reader supplies finite source inputs, not a source certificate. The
three ell=0 groups (0,13,23) have exactly zero radius-history response in
this family; their baseline allocations remain present. Group0 has its
analytic owner and no synthetic mode panel is created. Infinite tails,
the remaining low/subgap accuracy, changed-history subtraction and
continuum errors require their own bounds. Contour bilinears and integrated
moments are not accepted as physical Gaussian mode columns. No claim of
full spectral coverage or physical gate closure follows from loading these
panels.
