# Family14_1 source insertion v7

v7 joins the accepted v6 aggregate with the 108 disjoint subgap rows. It does not add the v5 total to the v6 total. v6 already contains the v4 rows plus the 18 threshold rows. v5 contains the v4 rows plus these same 108 subgap rows. The two added sets meet only on the shared 1042 rows.

Family `(14, 1)` now has 1168 of 1168 signed rows. Missing keys: none. No new source capture and no new solver. The saved subgap proof replay was not repeated. Its tree hash is unchanged, `c310e9023e945fad5985651cc0811e0f737bdff6d8eab352e091b52ac05ee161`. Original weights are applied once. The local gate stays OPEN.

## Partial N and beta

The continuous insertion of these 1168 rows is

$$
N\le 3646674054457366029024452496769420151764188105426118876609\times 2^{-229},
$$

$$
\beta\le 5793343250846686183515013928067141454707331929528967339969\times 2^{-230}.
$$

Read as decimals, those bounds are about \(4.226958072607968\times 10^{-12}\) and \(3.357610065481873\times 10^{-12}\). They are larger than the accepted v6 bounds because the 108 subgap rows are now inside the sum. They are not a bound on every source family, on quadrature, on numerical field error, or on value-assembly arithmetic.

The stored weighted moments are mantissa `7786760115494097` at exponent `-86` and mantissa `5541269878732561` at exponent `-83`.

## Signed-zero identity

NumPy 2.5.1 drops the negative zero at `C[0, 1].real` of each negative fiber while assembling the negative partner. The numeric arrays remain exactly equal: nonzero differences 0. Restoring that one sign bit on each fiber recovers the negative source and preparation digests stored in v6. The low32 batch's two raw source digests are still `641b187a3c37aa47a5841b559163eba5f78ec3ca2c79b535dc21dae8db657bdf` (signed-zero spelling) and `63e52bb2986c13bd72d91f3965a6f48b3d11380a6d29ae8326658c4fb2bc8d2c` (cleared spelling, which is the live digest). The same corner restores every other negative batch digest used by v6 and by the subgap rows. 584 negative rows are therefore recorded under the live cleared digest. Positive digests were already the live digests. Energies were unchanged.

The historical propagator object checked by the threshold closure is git blob `31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227`, not HEAD. The closure checked 67 source hashes and 64 archive digests at `eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9`.

## Budget

Completing family `(14, 1)` removes the remaining signed-row gap inside this one positive-angular family only. `physical_upstream_budget_component` stays null. Quadrature error, numerical field error, value-assembly arithmetic, and every other source family stay out. This record does not close the physical local gate.

## Files and checks

The aggregate owner is `scripts/derive_nsc_ks_source_operator_majorant_v7.py`. The result is `results/development/nsc-ks-source-operator-majorant-v7.json`, sha256 `c9a98578bf53ee02903403f4b63550ddb1d67ecf870d656f7b1f2ee74de44c12`. Accepted v6 remains sha256 `15615d1ee0206ed0fb6b8c3a10df8d860c581b2ebc206bef2242492f0154d517`. v4 and v5 bytes were not rewritten.

From the repository root, after the validation environment is ready:

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_ks_source_operator_majorant_v7.py --record
```

`--record` checks the witness tree, payload hashes, threshold closure, signed-zero spelling, disjoint census, mutated preparation, negative epsilon, and duplicate rows, then writes the file once. It does not call the 54-row proof replay. The rebuilt validation environment has no pytest, so the separate test module was not launched.
