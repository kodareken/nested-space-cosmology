# FGC-1-SGB1-CTL1-SOL1-PREF1 — independent scoped outcome binder

This compact binder independently authenticates SOL1-FRZ1 commit
`a3fdbf5176d7edb5ab90ebe49c136780985e0f3e` and parent
`2be811bd7868904e00e2acddac75186cde169408`, challenges the tracked freeze
config/result hashes, and reconstructs the SOL1 policy and outcomes directly
from the committed core without importing the FRZ1 certificate.

The independently bound nominal branch remains exactly:

```text
interval_inconclusive / prefix_endpoint_not_inside_declared_tube
```

The authenticated `A_chi=3` endpoint enclosure is wider than the frozen
`rho=1/8` tube. There are zero unique tubes. This is wrapping, not on-orbit
compactness nonpass, Jacobian nonpass, trapping, or SGB-L rejection. No
policy or physical-family parameter is retuned. The `A_chi=1/8` path remains
a control only.

Ordinary verification is compact-only and Git/source/raw/`runs`/reconstruction
blind:

```bash
python3.14 -B scripts/reproduce_fgc_sgb1_ctl1_sol1_pref1.py --check
make fgc-sgb1-ctl1-sol1-pref1
```

Explicit initial reconstruction may inspect the frozen Git commit and the
low-level SOL1 owner, but the binder writes no state. A later matched-data
redesign is separately owned and must not mutate this result.

`SGBL_branch_owned_and_healthy`, holdout, execution, trajectory trapping,
actual-orbit trapping, model rejection, spacetime/IBVP theorem, mechanism,
and physics all remain false.
