# Physical-metric subspace step

`scripts/derive_nsc_ks_gate_lm.py` builds one Levenberg–Marquardt step from measured nodal residuals. It does not certify the gate.

The clock is the recorded CF4 setting: length 1.6, grid 256, max step 0.0005. At that clock the last step-halving moved N by `5.53e-7`. A later acceptance must lower the merit by more than ten times that movement. The prediction stored in `results/development/nsc-ks-gate-lm-subspace.json` is a linear model, not a residual.

The metric is the Chebyshev Gram on the axial support: H^3 for `w`, H^1 for `U`, scaled by the squared columns of the certified principal matrix. The trust radius is the physical size of the iterate6 repair in that metric. The first Jacobian uses eight coordinates, `U0`–`U3` and `w0`–`w3`.
