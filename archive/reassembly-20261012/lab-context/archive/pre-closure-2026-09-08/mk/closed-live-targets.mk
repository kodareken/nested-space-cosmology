# Closed historical mutation entrypoints.
#
# The underlying scripts remain byte-for-byte historical evidence. These Make
# names are retained only to refuse accidental reruns deterministically before
# a Python interpreter, Git preflight, or campaign store can be opened.

CLOSED_HISTORICAL_STATE_TARGETS := \
	run-fgc-pro13-calibration \
	run-fgc-pro14-calibration \
	recover-fgc-pro19-sid1 \
	resume-fgc-pro19-sid3-event1 \
	run-fgc-pro19-event1 \
	resume-fgc-pro19-event1 \
	run-fgc-tdg8-successor-event1 \
	recover-fgc-tdg8-rcv1 \
	resume-fgc-tdg8-rcv1 \
	run-fgc-tdg8-rcv2 \
	run-fgc-tdg8-rcv3 \
	run-fgc-tdg8-rcv3-rec1 \
	run-fgc-tdg9-ar1 \
	run-fgc-tdg9-loc1 \
	run-fgc-tdg9-loc2 \
	run-fgc-tdg9-ti1 \
	run-fgc-tdg9-ti2 \
	run-fgc-tdg9-ac1 \
	run-fgc-tdg9-ur1 \
	run-fgc-tdg10-qa1 \
	run-fgc-tdg10-qa2 \
	run-fgc-tdg11-msel1 \
	run-fgc-hlt17-srcq1-rec1 \
	run-fgc-pro20-ev1

.PHONY: $(CLOSED_HISTORICAL_STATE_TARGETS)

$(CLOSED_HISTORICAL_STATE_TARGETS):
	@printf '%s\n' \
		"CLOSED historical operation '$@' completed or terminal;" \
		"state-changing rerun is not authorized;" \
		"use the named compact verifier" >&2
	@exit 2
