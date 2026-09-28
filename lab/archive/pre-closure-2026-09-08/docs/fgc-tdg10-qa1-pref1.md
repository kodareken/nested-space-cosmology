# FGC-1-TDG10-QA1-PREF1: independent QA1 terminal binder

`FGC-1-TDG10-QA1-PREF1` independently binds the completed four-leaf QA1 RUN1
retry-3 SSPRK3-on-inherited-SBP4 exact complete-C all-18-channel terminal
produced under immutable authority commit
`aeaf0498fd37a51d0e2c7efcf695147f7a9822e7`, whose single parent is
`49c514ac283ae3ec8091a6638442da6fdaf58db0`. The compact result is
`results/fgc-1-tdg10-qa1-pref1.json` with SHA-256
`3a107bcede479df4326aac5e194beefa1042860a8be498d7dc5d0519fdbc8603`. Its
artifact class is
`independently_bound_retry3_ssprk3_sbp4_exact_complete_C_all_channel_pass_terminal`
against target protocol `FGC-2-SF1-PROTO18`.

The one-time construction path uses no-following, race-checked reads for the
four raw leaves, verifies their exact SHA-256 identities, authenticates the
full 20-blob authority delta, and checks the scientific campaign store before
and after against its sealed 115-leaf identity
`5917d70de3080dcc6b7abeb7fdb7cc3370c6140cbeb37bde0c142d6a2d36a445`. Four raw
hashes remain unchanged. The store remains unchanged.

## Independent replay

The binder imports no QA1 runner serialize/validate/reduce/run/status code,
the QA1 authority module, the exact complete-C runtime, or the exact admission
adapter, and it does not reconstruct the seven SSPRK3 proposals from those
owners. Trusting the raw channel summaries without that replay would be a
common-mode failure with RUN1. The independent replay:

1. authenticates authority `aeaf0498…` and parent `49c514ac…`;
2. no-follow, race-checks the four raw leaves against the hashes below;
3. authenticates the 115-leaf store before and after against the sealed snapshot;
4. independently reconstructs 7 shadow paths / 7 proposals / 28 SSPRK3
   stage+endpoint records;
5. re-encloses every D01/D12 complete-C interval with both
   `localize_absolute_maximum` and `localize_absolute_maximum_independently_v2`;
6. reclassifies all 18 channels through design-only `classify_tdg6_channel`;
7. requires complete admission if and only if `failed_channels` is empty if and
   only if every channel admission flag is true;
8. requires the raw terminal class if and only if that independent reduction.

Route disagreement, malformed input, or a resource ceiling is fail-closed and
cannot be rounded into a pass.

## Bounded outcome

All 18 channels independently classify `resolved_order_pass`.
`failed_channels=[]`. Both exact rational localizer routes agree. The
unchanged discriminator is `8 U_12^2 <= L_01^2`, including equality, which is
the sufficient observed-order condition `p >= 3/2`. Candidate ceilings remain
per channel: 16,352 for D01 and 32,704 for D12. Primary refinement depth is
160. No absolute tolerance, signal normalization, debit cancellation, or
post-outcome threshold change is permitted.

The tightest channel is `u:R`. The exact ratio `8*U12^2/L01^2` is approximately
`0.9931888298284983`, leaving about `0.006811170171501636` below the pass
boundary.

The four raw leaf hashes are:

| Leaf | SHA-256 |
|---|---|
| `manifest.json` | `cdea5f1148486b6d0c2c5a3fee222a224ede758680e481efb8d6a78672c52b36` |
| `terminal.json` | `c0aea714469dd79874220f12011fb6c46a8468c7d5866fcffa6daffe06427c98` |
| `fine-endpoint/69478ba14677cec5dfcc7a2f802f2e5a5fe729907a225a42b560d489c55e2c3e.json` | `2c4a560c87d08f43dec5644d5af65c44663fbd77cb07476742f58a9a0e448532` |
| `payloads/ea961ef96a03ccfd572c53e5bfac56cd9898cbe464a5f10d928022c9a7359a9e.npz` | `740b371a9720b53e28164d67d22364172df17a2a30a738043c12033a722eed4f` |

The diagnostic endpoint remains unaccepted. Its content-addressed descriptor
is `69478ba14677cec5dfcc7a2f802f2e5a5fe729907a225a42b560d489c55e2c3e` and the
payload semantic hash is
`ea961ef96a03ccfd572c53e5bfac56cd9898cbe464a5f10d928022c9a7359a9e`. Retry 3
restores generation 9 / journal sequence 10 only. Sequence 11 remains the
historical rejection. Generation 10 is not an executable retry-3 predecessor.

## Compact and live routes

Ordinary verification consumes only the tracked compact certificate and is
raw-, store-, shadow-, and Git-blind. There is no implicit `--live`.

```bash
make fgc-tdg10-qa1-pref1
make verify-fgc-tdg10-qa1-pref1
python scripts/check_repo.py --only-tdg10-qa1-pref1
python3 scripts/reproduce_fgc_tdg10_qa1_pref1.py
```

Raw-blind `make verify-fgc-tdg10-qa1-pref1` passes 21/21 tests and the focused
repository audit.

`--live` is explicit one-time construction history, not ordinary verification:

```bash
python3 scripts/reproduce_fgc_tdg10_qa1_pref1.py --live
```

## Nonclaims

At its construction boundary, the compact result recorded only the conditional
statement that any later state adoption would require a separate prospective
freeze. It did not authorize RA1 or accept the diagnostic endpoint.

The later QA2 design supersedes the earlier idea of an RA1 endpoint transplant:
the diagnostic used SSPRK3 time stepping on the inherited SBP4 spatial
operator, so it may not be relabeled as an RK4 member or adopted as an old
campaign state. Current authorization is limited to temporal method/admission
design. RA1 does not exist and is not authorized.

This does not claim:

- an accepted diagnostic endpoint or campaign state;
- production SSPRK3 comparator status or two-method agreement;
- a common event or GR-0 calibration;
- SGB-L or FGC-QR execution;
- candidate execution, mechanism, retained-EFT, transition, or physics;
- RA1 authorization.

`FGC-1-TDG10-QA1-FRZ1` remains the historical compact freeze: at freeze it was
unexecuted, measured zero shadows, and authorized only the later diagnostic
RUN1. The freeze certificate is not rewritten by the gitignored raw terminal
or by this independent binder.
