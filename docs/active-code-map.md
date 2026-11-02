# Active code map

| Current work | Implementation |
| --- | --- |
| Offset history, distance kernel and winding geometry | `src/recursive_horizons/nsc_atlas_offset_generator.py` |
| Independent geometric and symbolic controls | `tests/test_nsc_atlas_offset_generator.py` |
| Existing finite-coverage calculation | `scripts/check_nsc_atlas_offset_generator.py` |
| Isolated scientific runtime | `scripts/run_python.py` |
| Derivation and explicit assumptions | [Offset construction](nsc-atlas-offset-generator.md) |

Run `python scripts/run_python.py -m pytest -q tests/test_nsc_atlas_offset_generator.py`.
The geometry-only reference calculation is
`python scripts/run_python.py scripts/check_nsc_atlas_offset_generator.py --winding-only`;
it prints a small result and writes no data files.
The old full-coverage calculation imports the conditional parent only when
requested; that parent is not the foundation's input.

The [previous code map](active-code-map-history-20261102.md) preserves
the other implementations. They remain executable at their original paths.
Their old next-step instructions do not restart a campaign.

Default ripgrep discovery follows `.rgignore` and shows the active desk.
For an explicit historical investigation, use `rg --no-ignore 'pattern'`
with the particular source or owner path. Test discovery and imports are
unaffected by these search exclusions.
