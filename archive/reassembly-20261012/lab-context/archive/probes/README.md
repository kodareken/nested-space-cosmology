# nf256 probe archive

Byte copies of the two temporary probes named by the v5 record. The v3, v4,
and v5 records are unchanged.

| Original path | Archive path | sha256 |
|---|---|---|
| `/tmp/nf256_seed_probe.py` | `lab/archive/probes/nf256_seed_probe.py` | `24957658ca54fdb97c791e1690c55d6b0eba544b044bd1a1d9183eaaa0e76dd7` |
| `/tmp/nf256_compensate_probe.py` | `lab/archive/probes/nf256_compensate_probe.py` | `3fa92afdcb292a81ac3ccd8b340ab4f3bb390e7e7e6e42e5ab787dfbe0e7f5b0` |

Both digests match `hashes_before` and `hashes_after` in
`lab/results/development/nsc-spherical-coupling-refinement-v5.json`.

`lab/scripts/derive_nsc_spherical_galerkin_refinement_v4.py` is archive-only.
Running it rewrites the partial v4 JSON from values already printed. It does
not import the solver, prolongate, or evolve.

`lab/scripts/derive_nsc_spherical_galerkin_refinement_v5.py` is the
self-contained measurement. It builds the nf256 state and the short window
itself. The probe paths are optional provenance: a missing file is hashed as
null and is not executed.
