# Historical tests excluded from active discovery

This directory contains immutable, non-canonical test history. Files here are
not imported by the active package, discovered by the canonical test runner,
or included in `verify-current-development`.

## SID3 real-store preflight

`test_fgc_pro19_sid3_real_store_preflight.py` was a transient, live-fixture
preflight proving that the then-current generation-eight SID3 projection was
not mutated by its read-only preflight. The live store subsequently advanced,
and the test already self-skipped once that exact fixture ceased to be the
store tip. Its historical source bytes are preserved unchanged at SHA-256:

```text
09d3515947530384063ad83c11faaea25ef7d43d527f23b15b50597067e63b3e
```

The moved file intentionally retains its original source, including its old
repository-relative `ROOT` expression. It is evidence to read, not a runnable
test in this location. Current SID3/PREF28 and later compact artifacts own the
post-attempt evidence.
