# Recursive Horizons build and verification router.
#
# Keep this first hop compact: immutable certificate recipes live in the
# historical include, active developer routes in the current include, and
# completed state-changing entrypoints fail before any Python process starts.

.DEFAULT_GOAL := test

include mk/historical-certificates.mk
include mk/current-foundation.mk
include mk/closed-live-targets.mk
