# FGC-1-PRO20-EV1-RSRC2-FRZ1 — high-capacity common-event rerun

This freeze authorizes exactly one fresh six-member GR-0 common-event attempt
from the sealed HLT15 time `t=23/16` to `t=3/2` in
`runs/fgc-2-sf1/pro20-rsrc2-event1/event`.

The earlier RSRC2 freeze `0f936b7...` was not consumed scientifically: its
first launch stopped before source construction because the process profile had
not replaced the runtime's separately imported container constant. Commit
`1a32046...` binds that identity seam; this refreeze is its exact successor.

RSRC1 ended at generation 30 with the accepted state preserved because the
long-lived parent exceeded its 1-GiB current-RSS ceiling while repeatedly
materializing the complete generation history. The terminal spent no
scientific retry and made no physical classification.

RSRC2 changes only the measured resource owner:

- generation-chain validation streams immutable generations and retains only
  the latest authenticated state;
- seed and child handoff objects are released after publication;
- terminal full-chain authentication runs in a fresh process;
- the parent ceiling remains 1 GiB;
- the isolated child ceiling is 8 GiB on the 24-GiB M4 host, with a 2-GiB host
  reserve and the same equations, grids, methods, thresholds, and retry caps.

On the sealed 31-item RSRC1 chain, the streaming route reproduced checkpoint
`2c404555d842067a226889ce70138cbb357a9fd4aec612b6223b806fcf3cabab`
and journal
`e0de9a9ceaccc085fda59c6122fad35fb22c4c36b29165dbb13bcda23c5242ae`
while reducing fresh-process maximum RSS from 397,443,072 bytes to 197,083,136
bytes.

The run begins from the original synchronized source and origin, not from the
partial RSRC1 endpoint. It cannot resume, take over, transplant an endpoint,
change a threshold, or open calibration/candidate work. A valid common-event
terminal opens the trapped-GR search. Any other terminal identifies the next
technical or formulation intervention and is not a rejection of Nested-Space
Cosmology.
