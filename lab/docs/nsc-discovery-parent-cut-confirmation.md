# One NF256 confirmation of the physical exterior-population response

The owning driver is `lab/scripts/derive_nsc_discovery_parent_cut_confirmation.py`.
It confirms the completed NF128 cut experiment at the same absolute central
proper clock, 2.8568786296681035, selected by the frozen7004a791 forecast's T2.25
baseline. This is one confirmation, with no NF512 ladder or continuum claim.

The baseline and +5% original-column1 occupation arms hold column0's occupation,
declared k=0.45388971484879426, initial physical Q/r and both spinor columns fixed.
The original NF128 source and reference frames extend by the canonical AP
isometry. Periodic physical geometry and nodal momenta interpolate to the new
odd geometry grid before encoding in the existing full identity frame. Columns
are neither reselected nor normalized. Actual Gram, CAR and before-correction
constraint differences are retained.

Each arm then solves its own finite C/D momenta with the frozen cut adapter and
hard seed anchor. Only momenta may change; positivity, sign or constraint failure
leaves preparation OPEN. This initial Cauchy correction is recorded separately
from the evolved response. There is no hidden geometry relaxation.

The two serial full leading RK4 arms use the existing admission and RK4 clock
increment with cap0.0005. A partial-step clock bracket reaches the absolute target
without extrapolation. The total600-CPU allowance includes authentication,
preparation and both arms. Budget/chart stops retain endpoint arrays and report
comparison null unless both arms reached the same proper clock.

Primary comparisons are original-column0 child probability, total child content
and child proper length. The NF128 measured gains were respectively
0.0004292078493, 0.0028420021143 and0.0047672394727. Holding c0 fixed makes the first
readout informative about the original body's response; tagged columns are not
independent particles or energies. Incoming normal radius/current channels are
reported through the original cut observer. The endpoint also reports finite
metric-rate projection defects. Those are not continuum Euler error bounds.

Input prepare/prediction JSON and NPZ and their dependencies authenticate against
7004a7915a9375c2d5a420853158c4cf19bcdca1. Executing source owners must retain their
frozen bytes. New confirmation producer hashes are pinned at root's new commit.
Output is a creation-only JSON/NPZ pair with the existing64-MiB bound. Numerical
movement against NF128 is an indicator, not a genuine continuum error bound.

```sh
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_confirmation.py
.venv/validation/bin/python scripts/lab.py scripts/derive_nsc_discovery_parent_cut_confirmation.py --run --producer-commit FULL_NEW_FROZEN_SHA --cpu-budget 600 --output results/development/nsc-discovery-parent-cut-confirmation-v1
.venv/validation/bin/python scripts/lab.py -m pytest tests/test_nsc_discovery_parent_cut_confirmation.py -q
```

The default is read-only preview. Tests use manufactured bands/constraints and
at most0.003 coordinate time. Root alone runs the physical confirmation.
