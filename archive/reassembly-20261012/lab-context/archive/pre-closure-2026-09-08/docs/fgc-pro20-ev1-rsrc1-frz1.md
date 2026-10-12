# FGC-1-PRO20-EV1-RSRC1-FRZ1 — one fresh common-event authority

This corrected exact freeze binds implementation parent `f1f6dedf96ebe14a54b1c2e7ff41f8a4caddfc8e` and authorizes one fresh six-member GR-0 common-event attempt from the sealed HLT15 `t=23/16` origin to `t=3/2`. The earlier authority commit `a99dcb3...` stopped read-only because the parent checked the wrong namespace field; it constructed no source and created no namespace.

The output namespace is `runs/fgc-2-sf1/pro20-rsrc1-event1/event`. The closed `pro20-event1` store is immutable historical evidence and is neither resumed nor copied as accepted state. Each seed and scheduled attempt runs in a fresh child under the unchanged 4-GiB ceiling; the parent alone authenticates and publishes.

The authority grants no calibration, SGB-L, FGC-QR, holdout, defocusing, mechanism, or physical claim. A typed terminal is a result to bind independently. Abnormal child exit or ambiguous publication remains forensic uncertainty.
