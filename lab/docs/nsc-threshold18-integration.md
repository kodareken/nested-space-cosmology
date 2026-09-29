# Threshold18 rows in source insertion v6

Imported tree `results/development/nsc-threshold18-import-v1` sha256 `bc7c5c2424810650acbeae20c4bedfc18548737f112725c63702eb3f9920a4c7` is unchanged.

Negative batch raw digests differ only by a negative zero at each fiber's covariance `C[0, 1].real`: 16 zeros in batch rows 32:48, and one in each of the nine imported witnesses. Local raw source digest `641b187a3c37aa47a5841b559163eba5f78ec3ca2c79b535dc21dae8db657bdf`. Original raw source digest `63e52bb2986c13bd72d91f3965a6f48b3d11380a6d29ae8326658c4fb2bc8d2c`, recovered by clearing those zeros and hashing the same weights and energies. Local raw preparation digest `91a903bab1a082cb8a01e8f28c03f745d6f3240b882442d55df3a87db4d1eccf`. Original raw preparation digest `af92bd35174e8d76f368c0aa3f031d630cff06aa128fa28c1e06c4392624ed81`. Positive source and preparation digests match. `source_error_moments` uses the local raw digests. Inserted errors are the imported per-row dyadics. Weights are applied once. Nonzero numeric array differences: 0.

Baseline `src/recursive_horizons/nsc_ks_energy_propagator.py` sha256 `5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94` is git blob `31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227`, not HEAD.

Closure checked: 67 source hashes and 64 archive digests at `eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9`.

v6 covers 1060 of 1168 signed rows, 108 remaining. `N` mantissa `5310277635095228096979754608917466329654725356268274396327` exponent `-230` (`3.0776428853051678e-12`). `beta` mantissa `4218374114141829394308763292324159580719549930561481200407` exponent `-230` (`2.444815501574293e-12`). Record sha256 `15615d1ee0206ed0fb6b8c3a10df8d860c581b2ebc206bef2242492f0154d517`. No new ODE. Gate open.
