# Phase -1 cleanup baseline — 2026-08-30

This manifest records the repository boundary immediately before the approved
context and repository-hygiene consolidation. It is an inventory, not a new
scientific result.

## Git and verification boundary

- Baseline commit: `f155dc17fb17ac414760d011cc135db39b802b67`
- Branch: `main`
- Tracked files: `1294`
- Baseline tracked-name manifest (`git ls-tree -r --name-only -z HEAD`):
  `71fa1aba10a422dd3bf3f32e7520a77cb2657e059df5115dc0fe8b75b7d65324`
- Baseline tracked-tree manifest (complete `git ls-tree -r --full-tree HEAD`
  records, including mode, object type, blob identity, and path):
  `c4d4c2dbeef0e2f84f14ec8fbee75d50e46e0c66c563e3915ccd1ffedbbd53fa`
- Canonical analytic partition at the last uninterrupted baseline run: `1680`
  tests
- Canonical evolution/restart partition: `352` tests
- `verify-fgc-sf1-foundation`: passed
- `make verify`: passed, including paper, eleven required/available data
  artifacts, and the final repository audit
- No scientific runner, live binder, writer, test partition, or verification
  process was active when the cleanup began.

Untracked host/user state at the boundary consisted of the zero-byte
`.qdrant-initialized` marker plus the user-authored `AGENTS.md` and a 992-byte
macOS bookmark named `PLAN.md`. The bookmark was moved to recoverable Trash
and replaced by the approved plan as real Markdown.

The ignored inventory consisted of scientific/local inputs (`runs/`, the two
Planck maps, and local harness state) plus generated residue covered by the
narrow existing ignore rules. At the start of cleanup the generated Python
residue comprised six `__pycache__` trees and 1,402 `.pyc` files, alongside
`.pytest_cache/`, `.ruff_cache/`, `build/`,
`src/recursive_horizons.egg-info/`, `paper/.build/`, empty `snapshots/`, and
the Finder files listed by the cleanup plan. Integration tests performed
before trashing regenerated one additional modular-checker cache; the final
pre-trash refresh contained seven cache trees and 1,413 `.pyc` files. This
expected generated drift is kept separate from the immutable `runs/` and
Planck identities below.

## Directory-size inventory

Sizes are the pre-trash `du -sk` observations in KiB. They are routing and
cleanup evidence, not scientific measurements.

| Directory | KiB |
|---|---:|
| `.git/` | 21508 |
| `archive/` | 804 |
| `collected-data/` | 2162836 |
| `configs/` | 1692 |
| `docs/` | 1632 |
| `paper/` | 4184 |
| `results/` | 20676 |
| `runs/` | 146608 |
| `scripts/` | 15276 |
| `src/` | 17076 |
| `tests/` | 13472 |

The generated-path inventory immediately before residue removal contained
1,777 ignored paths. Its NUL-delimited
`git ls-files --others --ignored --exclude-standard -z` SHA-256 was
`356bd42c166dabb73bcb3299b56fc817c303f7ed3d023c27c9b8457217a0473d`.
That digest is intentionally expected to change when the generated residue is
moved to Trash; the scientific run-tree and Planck hashes are the invariants.

## Canonical tracked and external identities

| Artifact | Bytes | SHA-256 |
|---|---:|---|
| `AGENTS.md` (user-authored creative orientation) | 36002 | `179967eee937f98ba224a122c44fbc1bb25825480ffb54ddd677f7511b13bda8` |
| `PLAN.md` (accepted implementation/finish-line plan) | 21058 | `6007167cf5d4e738104186ded4e2ab981aaed3089a79746035ed1ac1d082d1b7` |
| `paper/recursive-horizons.pdf` | 381076 | `88b072e4abae6b7122b2b35aae11181af305ebb048a3029f710b148f54b52ad2` |
| `results/fgc-1-tdg10-qa2-pref1.json` | 253881 | `08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6` |
| `collected-data/smica_2048.fits` | 2013312960 | `60952c645eb33d151905ddf5837477e15ca02a3261feca4ceca3c3feece4f9ac` |
| `collected-data/mask_common.fits` | 201335040 | `23f49a3479073408b83a7f1bb2cca490d9be7f934b2e0e891524b6a1d06a5fc4` |

The two Planck files are intentionally retained locally.

## Ignored scientific run inventory

The complete `runs/` tree contained `297` regular leaves and `149556321`
bytes. The cleanup baseline digest is:

```text
CLEANUP-BASELINE-TREE-v1
SHA-256 866082d52e8cfd7d2d10540c2b9f018385969d084eda096139f3270a0b1c9ace
```

The digest is computed over the lexically sorted stream
`repository-relative-path NUL byte-count NUL leaf-SHA-256 LF` prefixed by the
domain line above.

| Run root | Leaves | Bytes | Baseline digest |
|---|---:|---:|---|
| `proto5` | 4 | 951889 | `38b528a30f4cfbebe8f161dad6279a0aaab130c4192ae86774d5982e8bd5eff9` |
| `proto6` | 4 | 552695 | `2854181da693ff36da60267ea05fed8a88d4c8b5922a8aeeedf6c8cb5e8cea69` |
| `proto7` | 6 | 2726698 | `256d7b319dfb0f2f153775759a69b4200480c6b6985d7934dbc384121d2cba55` |
| `proto8` | 4 | 1421599 | `07b4f56640a0b9dafc133849523de497710cbf24dd5299055544485eeaadae9a` |
| `proto10` | 4 | 1718903 | `172905b134df494c362e614fa6ec87b9621010cecb102842e1fb8fc4b0c9efbe` |
| `proto11` | 4 | 3787872 | `3def5d535b593b4e5ae0e9cae294462ce52bbda8d17be35213a39aa6385e6da6` |
| `proto12` | 4 | 25588939 | `d86e3d9a2e71d3fd035d93195ede10a0459b335143ef1f2fdc7da0770c97443a` |
| `proto13` | 4 | 7832771 | `e8d1a7ca9bfe8ae4b594fea0e330ad3d8a4fa77c238b527a95bee4b319218219` |
| `proto14` | 4 | 4294921 | `b19f845cad7addbf64f27dc7a89bc986b836f74e5eb30721356e17100d15c1b5` |
| `proto17` | 56 | 6054860 | `8e0b1f3ba70a4de7061d322ed041ff0328ae4e6ca7ed78da79421883f47ad777` |
| `proto19` | 59 | 6342040 | `1c5dabc1bd4d9470a82f7b549a35db4ba43064b95089625a3ba70796efeec57a` |
| `rsp1` | 4 | 4904985 | `4896ea4ca2de04df1da701fd72d73b32ecc700be8bc45105241431dc1f6a27be` |
| `rsp2` | 4 | 1486466 | `c1e3f46fe16f40e74faf409eea7b305c2d158773b3282d3bf3231f5277b332b1` |
| `src2` | 4 | 3352779 | `5084405e5048be9c5d1c1b3f5478c4fc314107839f383d97b3c82df542ffa9ae` |
| `tdg8-rcv3` | 116 | 25737875 | `67a6deb68c15d59d80078932a9f7559e24e5db9993097fec3ad835bcc49cce97` |
| `tdg9-ar1` | 2 | 410422 | `adb685aefa1eed52147f5d25c35bcca63766297a0f195636bceedba92e655df1` |
| `tdg9-loc2` | 2 | 51362638 | `8485f340056aa57d34b72136f1715109815dae842a47e43f9105617488b317e5` |
| `tdg9-ti2` | 2 | 665654 | `7aa8456e9aad0e7da517f07a5dcfaf26984be74453348067f3b6f8d07ee82b1a` |
| `tdg9-ac1` | 2 | 29042 | `8459c03fd0cb117914f45ceb0c33bc169783cbed5fe85e3145d5c2a425796a37` |
| `tdg9-ur1` | 2 | 10852 | `d512f648a6883703c44b67fb51eb832e00427641b661cf0775709310e47876f5` |
| `tdg10-qa1` | 4 | 294272 | `88bffce4785b68706db368c750c17c1383d46b4ad6d015c5e1e5dd93bab95e12` |
| `tdg10-qa2` | 2 | 28149 | `d1811945f8f193d1cd22bc9b869cb42e20764a98dcffc8e439d20d370f77ed82` |

Every run root is retained. Terminal locks, repeated payloads, and apparently
duplicate store projections remain evidence and must not be hardlinked or
deduplicated in place.

## Cleanup exclusions

The cleanup must not run `git clean`, `make clean`, Git pruning, history
rewrites, broad lock/tmp ignores, or recursive removal over `runs/`,
`collected-data/`, tracked results, the canonical manuscript/PDF, or the
immutable archive.

## Phase -1 closure observation

After consolidation and before the two cleanup commits:

- `make verify-current-development` passed all six latest compact
  certificates, 46 focused Phase -1 tests, catalog regeneration, and the QA2
  compact repository gate.
- `make verify-fgc-sf1-foundation` passed its complete sealed certificate,
  lifecycle, analytic, evolution, and repository-audit chain.
- Canonical `make verify` selected 1,731 analytic and 352 evolution/restart
  cases; the routed suites completed with only their declared sealed/live
  skips, rebuilt the paper, verified all eleven available data artifacts, and
  passed the final repository audit.
- `python scripts/verify_data.py --require` independently verified both Planck
  maps and every tracked data artifact.
- The rebuilt 49-page `paper/recursive-horizons.pdf` is 335910 bytes with
  SHA-256
  `9915f0ebd8561d25de6058ef20d7c967ec6cb00d8bb3638de8d88d259e9b516a`.
- The QA2-PREF1 compact result remains
  `08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6`.
- All 297 `runs/` leaves, 149556321 bytes, every per-root digest above, and the
  global `866082d5...` tree digest remained byte-identical.
- Seven regenerated `__pycache__` trees containing 1,446 `.pyc` files,
  `.pytest_cache/`, `.ruff_cache/`, `build/`, egg-info, `paper/.build/`, empty
  `snapshots/`, and the resolved Finder files were moved to recoverable Trash.
  No scientific run, source data, compact result, owner, or historical runner
  was removed.
