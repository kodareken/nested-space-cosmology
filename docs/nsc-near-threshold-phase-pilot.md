# Near-threshold source control

The original subgap8/14_1 row 7 has E=1.5594631247710442, below mass pi/2.
The same source, Jost trajectory settings, horizon frame and upstream Bloch
transport were used in all controls. Only the enclosure subdivision and tube
changed. No source column or physical parameter was changed.

| Defect subdivisions | Phase tube | Outcome |
|---:|---:|---|
| 4 | 1e-5 | Tube invariance failed at cell 1,138 of 1,544 |
| 16 | 1e-5 | Tube invariance failed at cell 1,248 of 1,544 |
| 16 | 1e-3 | Completed all 1,544 exterior cells and 670 Bloch cells |

The successful control returned phase error 4.983422330112534e-5, positive
covariance error 1.4096431640903264e-11 and negative covariance error
1.4098991103763073e-11. Total CPU time was about 44.59 seconds on this Mac.
Values here are diagnostic decimal reports; no scientific record consumes
them. A production record must retain the directed endpoints and witnesses.

The wider tube is an explicitly checked error domain. It does not accept a
failed enclosure or alter the covariance tolerance: the comparison bound must
stay strictly inside that tube on every cell, and the resulting phase error
is propagated through the frame, source coherence and Bloch validation. The
resulting signed covariance errors, not the tube width alone, are what must
later enter the source-error aggregate.

This resolves a method pilot only. No saved-witness production record was
emitted, and the verified 866-row coverage and gate status remain unchanged.
The tighter-tube failures are not evidence against the physical class.

```sh
.venv/validation/bin/python scripts/lab.py scripts/pilot_nsc_near_threshold_source.py --run
```

This is a single selected-row control, not a full campaign. It reuses the
[independently reviewed signed transport method](nsc-massive-jost-mixed-transport.md).
