# Source insertion with the complete middle-energy window

The v2 successor keeps the exact-response norm method of
[nsc-ks-source-operator-majorant.md](nsc-ks-source-operator-majorant.md).
It replaces the old single middle-energy pair with the 48 positive rows and
48 actual signed partners from the resumable middle-window calculation.
The high-energy rows and low pair are reused. Thus 866 signed rows are covered
and 302 remain in this one family; other families and source quadrature remain
outside this result.

The builder replays every saved middle witness before forming the weighted
source-error moments. It rejects incomplete windows, duplicated rows, altered
signed energies, mismatched preparation digests, and missing signed bounds.
The original single middle pair is removed before the window is added.
Every checkpoint JSON and NPZ is bound by its lab-relative path and hash.
The older v1 scientific record is left unchanged.

This bounds the source-preparation contribution to the exact continuous N/beta
response. It does not change the candidate residual, supply numerical field
error, or complete the physical upstream-budget component. The latter remains
null and the local gate remains OPEN.

Owners: `scripts/derive_nsc_ks_source_operator_majorant_v2.py` and
`tests/test_nsc_ks_source_operator_majorant_v2.py`. Use `--record` once and
`--check` for a full saved-witness replay. Runtime does not enter the scientific
record. No new source evolution is launched by the aggregate or its replay.
