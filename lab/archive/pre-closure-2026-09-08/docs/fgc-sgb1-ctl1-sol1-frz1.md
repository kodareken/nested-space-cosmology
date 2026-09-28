# FGC-1-SGB1-CTL1-SOL1-FRZ1 — compact orbit-continuation freeze

This prospective compact freeze binds the committed
[SOL1 instrument](fgc-sgb1-ctl1-sol1.md) on implementation image `2be811b`.
It freezes the exact `A_chi=3` family, authenticated `(C,k)` prefix,
16-step remaining-support grid, `rho=1/8` tube, depth zero, 16384-bit cap,
four predecessor hashes, and the already measured scoped outcome.

The nominal contract SHA-256 is
`5fc44f7023eb82632309c6c15120f2e62e5e6483911e18589579576b36654f4b`.
Its exact classification is:

```text
interval_inconclusive / prefix_endpoint_not_inside_declared_tube
```

The authenticated prefix endpoint enclosure is wider than the prospectively
declared tube. There are zero unique tubes. This is not an on-orbit
compactness nonpass, not a Jacobian nonpass, not trapping, and not SGB-L
rejection. The tube, amplitude, support, threshold, chart, and resource cap
are not retuned.

The named `A_chi=1/8` control remains a control only. Its contract SHA-256 is
`661f06f9e6d6fa91d64f1767cf8b8b810c6ce17a4e9bd7ccf8a293b72ed946b9`.
It never replaces the nominal family.

Ordinary verification is compact-only:

```bash
python3.14 -B scripts/reproduce_fgc_sgb1_ctl1_sol1_frz1.py --check
make fgc-sgb1-ctl1-sol1-frz1
```

That route is Git/source/raw/`runs/`/reconstruction blind and launches no
runner. Explicit initial construction may directly re-execute the committed
SOL1 controls, but this freeze does not claim independent reconstruction.
The separate [FGC-1-SGB1-CTL1-SOL1-PREF1](fgc-sgb1-ctl1-sol1-pref1.md)
subsequently authenticates this freeze and independently reconstructs the
unchanged nominal and control outcomes.

All aggregate flags remain false: `SGBL_branch_owned_and_healthy`, holdout,
execution, trajectory trapping, actual-orbit trapping, model rejection,
spacetime/IBVP theorem, mechanism, and physics. A later matched-data redesign
must be separately frozen before another continuation attempt.
