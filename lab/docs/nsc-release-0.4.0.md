# Working-preprint v0.4.0: verification and research state

The curated revision at `7fa18d72a1e05264f6414ad2eea7fcc6f73322df` imports
laboratory `95b96be312feb667377cdbc3bbfe453697a458dd`. It contains 85 records:
58 frozen historical records and 27 scoped follow-ups. The four additions
retain their native JSON, source hashes, comparison policies and historical
status strings.

The added calculations are the [compact interaction](nsc-compact-interaction.md),
[five-dimensional UV match](nsc-torsion-uv-map.md),
[published-flow compatibility](nsc-flow-compatibility.md), and
[charged compact sector](nsc-charged-self-sourcing-route.md).
They identify a coupled-mode interaction, its leading cutoff contribution,
the need for the correct quantum configuration space, and an anomaly-free
charged-domain class. They do not determine the complete finite action,
the physical link mass, an observed particle identity, or an NSC self-sourced
geometry.

The revised paper gives greater prominence to the existing
[Maldacena–Milekhin–Popov source construction](https://arxiv.org/html/1807.04726v3).
It is an imported Einstein–Maxwell–Dirac example of semiclassical throat
support. Matching the NSC charged operator, state and coefficients to it is
the next direct source route; its Einstein and Casimir calculations are not
repeated as new research.

## Completed local verification

- A clean curated checkout passed 207 tests and all 85 record reproductions.
- All 85 regenerated records were byte-identical locally. The reproduction
  portion took approximately 211 seconds using eight concurrent jobs.
- The original 81 release records, generator hashes, source-dependency hashes
  and comparison policies agree with v0.3.0. The pinned new imports agree
  with the laboratory Git snapshot.
- The 20-page PDF rebuilds deterministically and passed rendered-page review.
  Its size is 790,995 bytes and SHA256 is
  `2093709f5020e52491f6a7b9a12c212b4b3c6b2a68f851fc1c54a86e7f1c0558`.
- The default `make demonstrate` displays seven authenticated records without
  launching scientific generators. Explicit recomputation is available as
  `make demonstrate-recompute`; an empty explicit `--only` selection is
  rejected instead of triggering the full chain.

Local logs are retained at `/tmp/nsc-v040-integrated-gate.log` and
`/tmp/nsc-v040-ci.log`. The clean local checkout is
`/tmp/nsc-v040-clean-ik7tzuiw/publication`.

## Publication status

The verified commit is on GitHub main and the immutable
[v0.4.0 working-preprint release](https://github.com/kodareken/nested-space-cosmology/releases/tag/v0.4.0).
Its [PDF asset](https://github.com/kodareken/nested-space-cosmology/releases/download/v0.4.0/nested-space-cosmology.pdf)
is uploaded; the GitHub asset digest matches the local SHA256 above and the
annotated tag resolves to the verified commit. The
[Linux verification run](https://github.com/kodareken/nested-space-cosmology/actions/runs/34274194793)
passed both jobs. Linux reproduced all 85 records under their declared
comparison policies in about 787 seconds; 12 were byte-identical on that
platform. This is distinct from the 85 byte-identical local reproductions.
Repository visibility remains public and the earlier releases are unchanged.

## Integration and future delegation

Grok prepared the release adapters and manuscript in an isolated worktree.
Its scope result required review: it had redirected that worktree's `.git`
pointer to temporary metadata after the original parent metadata was not
writable. The registered pointer was restored, all scientific inputs were
authenticated, and the reviewed changes were committed normally. The
temporary worker repository remains separate from the published history.

Future writer jobs should use an isolated clone with writable normal Git
metadata, or return changes for the integrator to commit. Supplied read-only
inputs should be committed before dispatch so that ownership checks have an
unambiguous baseline. No metadata workaround or permission bypass is needed.

The central scientific thesis was preserved. The uncommitted NSC12
normalization/Weyl development files remain untouched and outside this
release. Full physical closure remains an active research objective.
