# Initial Cauchy data for a nonzero shift current

`initial_cauchy_data(grid, state)` builds the initial point for columns already prepared on a Galerkin grid. It can also take `phi0` and `phi1` on the fermion nodes. The action, column source, radius Newton, and geometric rates stay with their owners. The returned state is not integrated. The saved preparation is the v5 diagnostic, [derive_nsc_spherical_galerkin_refinement_v5.py](../scripts/derive_nsc_spherical_galerkin_refinement_v5.py). The episode driver `scripts/derive_nsc_spherical_feedback_episode.py` reuses that preparation. Its saved \(T=0.05\) record is `results/development/nsc-spherical-feedback-episode-v1.json`. This note is the constructor, not that trajectory. The `1e-8` line below is this constructor's solved-residual flag. It is not an automatic veto of the later v5 diagnostic, whose full Hamilton residual lies above that flag and was still measured.

## Reduction checked on the owners

`shift_constraint` and `hamilton_constraint` are

\[
\mathcal D=p_r r_x+p_\chi\chi_x-Q\partial_x p_Q,
\]
\[
\mathcal C=\frac{p_Q p_\chi}{2F_\chi}+\frac{\Pi^2}{4ZQ}+\frac{Z r_x^2}{Q}-QV-2\partial_x(F_x/Q),
\]

with \(\Pi=p_r-(F_r/F_\chi)p_\chi\). The source owner sets \(j=\mathrm{force}_\beta/\mathrm{d}x\) and \(\mathrm{force}_\beta=-M\,\mathrm{Pmom}\), where \(M=4\kappa\) is applied once.

On \(\chi=0\), \(p_r=0\), \(p_\chi=0\) and the constant gauge ratio \(Q=b_0/a_0>0\), \(\Pi=0\), so \(\mathcal C\) does not depend on \(p_Q\) and

\[
\mathcal D+j=-Q\partial_x p_Q+j.
\]

A periodic \(p_Q\) matches \(j/Q\) on the mean-free part. Its additive mean does not enter \(\mathcal C\) or \(\mathcal D\); that mean is declared \(0\). The same constant does enter the initial auxiliary rate from `geometric_rates`,

\[
\dot\chi=L p_Q/(2F_\chi),
\]

while \(\chi\) itself stays \(0\). Because \(\Pi=0\), \(\dot r=\beta r_x\) and the chart proper radius velocity is \(0\).

\(p_Q\) is solved in the geometry subspace. The radius is the existing projected Newton solve for \(\rho=\mathrm{force}_L/\mathrm{d}x\). On these sources \(\rho\) does not depend on \(r\). A nonzero mean of \(j\), and any current outside the geometry subspace, stay in the shift residual. They are recorded. The current array is not edited.

The conditional gaps are part of the same check. For \(j=0.3+\sin(2\pi x/L)\) the shift residual is the constant \(0.3\), and \(p_Q\) matches \(-\cos(kx)/(kQ)\) to `3.375077994860476e-14`. Adding `0.2` to \(p_Q\) leaves \(\mathcal D\) and, at \(p_\chi=0\), \(\mathcal C\) unchanged. At \(p_\chi=0.01\) that constant changes \(\mathcal C\) by `0.14141021696312883`. A pure constant current \(j=1\) produces \(p_Q=0\) and residual \(1\).

## Recorded slices

Tolerance for a solved residual is the existing `1e-8`. Evidence is `lab/results/development/nsc-spherical-cauchy-data-v1.json`.

| Slice | Current mean | Current max | \(p_Q\) max | Projected \(\mathcal D+j\) | Full \(\mathcal D+j\) | Projected \(\mathcal C+\rho\) | Full \(\mathcal C+\rho\) | Chart proper |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| One positive carrier, \(nf=32\) | `-0.7363107781851074` | `0.7363107781851865` | `1.08e-14` | `0.7363107781851075` | `0.7363107781851805` | `1.16e-11` | `6.75e-10` | `0` |
| Rank-6 packet, \(nf=64\) | `1.24e-17` | `5.42e-14` | `7.76e-15` | `1.24e-17` | `4.89e-14` | `4.73e-11` | `2.3751081949407538` | `0` |
| Finite-window lobes, \(nf=64\), \(nq=256\) | `0.002391744843402927` | `6.48259859941834` | `3.6949927676682854` | `0.002391744843430757` | `0.7438151949709526` | `2.62e-11` | `5.426001843815605` | `0` |

The carrier current equals \(-M n_0 k/L\) to `3.3e-16`, with \(n_0=0.75\) and \(k=2\pi\cdot 2.5/L\). Its shift residual is that constant. Radius stays positive (`3.3619479854663585`).

The rank-6 packet is the existing near-zero-current preparation. Saved refinement v5, not reloaded, reports current max `3.2542857299044724e-12` and chart proper velocity `0` at \(nq=1024\). The new constructor leaves \(p_Q\) at roundoff. The projected Hamilton residual is `4.730299485444789e-11` and is solved. The full Hamilton residual is `2.3751081949407538`. The held-out maximum is `2.375108194930538`. Those are separate recorded fields; `solves_full_hamilton` is false. \(r\) runs from `4.209679963964945` to `4.8757670637835355`.

The lobe slice uses the completed embedding profiles, sampled with the existing \(\sqrt{\mathrm{d}x}\) spin frame and pulled back by \(U_f^\dagger\). Band retention is `0.9986430084314399`. Occupations \((0.75,0.75)\), \((0.5,0.5)\), \((0.25,0.25)\) are equal inside each regional pair. On the classical piecewise derivative the occupation-weighted momentum integral is `6.661338147750939e-16`. The spectral current used by the source owner does not have that mean: `0.002391744843402927`, integral `0.019133958747223415`. The geometry-band antiderivative removes every nonconstant projected mode to `2.68e-14`. The projected shift residual equals the mean. The full shift residual `0.7438151949709526` is the unresolved current. The same momentum diagnostic at \(nq=512\), without a second radius solve, leaves mean `0.0027421102411685203` and full residual `0.7319317569372197`.

\(Q=b_0/a_0\) stays `0.251324933269135`. Embedding \(r\) runs from `4.2042616346503445` to `4.842654377102552`. The Hamilton residual is unchanged by this \(p_Q\). The initial \(\dot\chi\) reaches `84.14095068867056` and matches \(L p_Q/(2F_\chi)\) with gap `0`. \(\chi\), \(p_r\) and \(p_\chi\) remain `0`. The chart proper velocity is `0`; the lifted proper velocity is `0.00012467329606885563`.

`build_record` calculates that record and returns it. `write_record` stores a payload only by exclusive atomic creation at an explicit new path. An existing file, a symlink, and the sealed development record are refused. The module command requires `--record` or `--output` and has no default destination:

```sh
python scripts/lab.py -m recursive_horizons.nsc_spherical_cauchy_data --record path/to/new.json
```

The sealed file remains `lab/results/development/nsc-spherical-cauchy-data-v1.json`, SHA-256 `2c6a9ede084e3a56ae7798e44ac017d1086e591f4a2165d26e75a71ff98150f7`. Its Galerkin source pin is `85d1f3dbbd82a38fe14d2405ad9782af6de68095052ff522ac6602e5d50a14ab`. The working-tree file `lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py` is `aa8a64f593e8f0d931eeee407d944b72568fdf4f2c53436048142c3473ec4911`.

`scripts/derive_nsc_spherical_feedback_episode.py` `verify_saved` compares the episode's frozen hashes with that working-tree pin and stops on `galerkin`. The producer, `results/development/nsc-spherical-feedback-episode-v1.json`, and `results/development/nsc-spherical-feedback-episode-v1.npz` stay the saved bytes. `scripts/check_nsc_spherical_feedback_episode_historical.py` reads them. It recovers the Galerkin source from commit `5f10ecd365843d1616e50eb16a20d7acd8377e2c` through `provenance.resolve_pinned_source_bytes`, checks the episode's full frozen dependency hashes, and checks the saved arrays against the v5 payload and the stored replay series. It parses the pinned import graph before any load. The historical Galerkin module imports the pinned coupling and conformal-source modules; those modules import `nsc_adm_source`, `nsc_covariant_operator`, `nsc_influence`, and `nsc_nested_qualities`, which are outside the frozen set. The historical module stays unloaded, and `run_episode` is not called. The working-tree module also imports `nsc_spherical_feedback_action` directly and is a different pin, so a current-code step is not a replay of the historical producer.

```sh
python scripts/lab.py -m pytest tests/test_nsc_spherical_cauchy_data.py tests/test_nsc_spherical_feedback_episode_historical.py tests/test_nsc_spherical_feedback_episode.py -q
```

Eleven tests passed in 4.99 s. Those runs left the sealed Cauchy JSON, the feedback JSON and NPZ, and the manuscript files at their recorded bytes.
