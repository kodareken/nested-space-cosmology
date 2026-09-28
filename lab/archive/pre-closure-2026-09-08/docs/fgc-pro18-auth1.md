# FGC-1-PRO18-AUTH1: calibration-only genesis authority content

**Status:** completed content construction only.  AUTH1 binds the compact
`FGC-1-PRO18-PREF26` evidence preserved at
`ab3dd2613f70028388e0658c873c33b4225233c7`, pins the declared production
source closure and evolution environment, and derives exactly one legal
PROTO17 generation-zero `GenesisSpec` for the GR-0 amplitude-3 calibration.
The payload is **calibration-only**: it names neither an opened holdout nor a
candidate branch.

It is deliberately *not* a self-authenticating launch record.  The result
does not contain an authorization commit, result path, result-file digest, or
any claim that its own bytes are authoritative. The later HLT15 launch received
an external authority tuple and independently verified committed Git/blob
bytes, imports, raw reimport, and only then the output root; this did not alter
AUTH1's own non-self-authenticating claim boundary.

## Bound construction

The source of the six restart/state descriptors is compact PREF26 evidence;
AUTH1 does not reopen either raw checkpoint or event journal.  It binds:

- the PROTO17 protocol config and freeze result;
- the synthetic HLT13 durability result;
- the actual PRO14/HLT12 calibration plan;
- the direct `runtime_module`, `adapter_module`, and `runner` source hashes;
- every project-local transitive source dependency reached by those sources;
- one complete observed evolution-runtime contract: CPython 3.14.3, NumPy
  2.5.1, Accelerate BLAS/LAPACK, Darwin arm64, and the canonical full metadata
  digest.

The campaign identifier is not copied from an earlier calibration.  It is
deterministically derived from the target protocol, GR-0/amplitude-3 branch,
calibration namespace, PREF26 evidence digest, and the frozen PRO14 run-plan
digest.  The result contains only the calibration path as a string required by
the pure GenesisSpec; it does not test whether that path—or the holdout path—
exists.

## Exact boundary

True claims are limited to canonical commit-independent content construction,
PREF26 compact-evidence binding, complete calibration-only GenesisSpec
derivation, and pinning of the production source closure (including the
standard-library-only HLT15 bootstrap/authority module), exact evolution
environment, raw-source config binding, and deterministic fresh campaign ID.
False claims include external Git/blob/live-import authentication, HLT15
authorization or execution, future-root observation, namespace or
generation-zero materialization, state advance, calibration, candidate
execution, and physical transition. In particular, candidate execution remains
unauthorized.

The focused verifier reconstructs the complete content artifact and attacks
every `GENESIS_INPUT_FIELDS` category, the holdout namespace string, a missing
closure member, source-hash drift, environment drift, and claim promotion.
