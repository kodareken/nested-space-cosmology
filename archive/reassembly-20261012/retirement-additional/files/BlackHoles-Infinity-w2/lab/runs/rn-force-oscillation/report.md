Worktree: `C:/Users/swedi/Projects/Science/nester-space-cosmology/BlackHoles-Infinity-w2`, branch `agent/blackholes-infinity/work2`. All five scientific deliverables remain there as new files; no existing repository file was changed and no commit was made. This report and the requested subagent status are administrative outputs at the requested lab paths.

# RN recorded-pressure oscillation audit

## Scope and replay

Finite read-only analysis of all 177 `station-[0-9.]+.npz` files in the six `nsc-rn-neutral-family-v1` cases. `prepare.npz` is excluded. Radius runs from r_inner = mu to 32 r_m, with r_m = 1 in the recorded units. Grids have 3921, 1750, and 869 points for mu = 1.002, 1.010, and 1.040; not every case has 869 points. These are the stored signed radial and angular pressures, not a newly derived pressure-gradient force. No evolution, smoothing, producer import, or power-law fit was performed. Conservation and constraint damping cannot be certified by this crossing-only consumer.

Owners: module `lab/src/recursive_horizons/nsc_rn_force_oscillation.py`, driver `lab/scripts/derive_nsc_rn_force_oscillation.py`, tests `lab/tests/test_nsc_rn_force_oscillation.py`, creation-only record `lab/results/development/nsc-rn-force-oscillation-v1.json`.

Run from the assigned worktree using the lab interpreter:
```text
C:/Users/swedi/Projects/Science/nester-space-cosmology/BlackHoles-Infinity/.venv/validation/Scripts/python.exe scripts/lab.py scripts/derive_nsc_rn_force_oscillation.py
C:/Users/swedi/Projects/Science/nester-space-cosmology/BlackHoles-Infinity/.venv/validation/Scripts/python.exe scripts/lab.py -m pytest tests/test_nsc_rn_force_oscillation.py -q
```
The driver succeeded and deliberately refuses subsequent overwrite. Replay creation in a clean destination checkout. All 10 tests passed. Both numerical kernels are float64, parallel cached Numba `CPUDispatcher`s. Independent numpy sign comparisons verified all 354 pressure crossing counts and all 177 source SHA-256 hashes after creation.

## Definitions and error statement

Opposite nonzero signs count as one crossing; adjacent anchors use linear interpolation. Exact-zero runs use their midpoint only if bounding signs differ. Same-sign touches and boundary zeros do not count. Half-wavelength is the mean of successive crossing separations, not domain width divided by crossing count. Fewer than two crossings gives no wavelength (—). Isolated reversals and separated clusters are not demonstrated periodic waves. Long gaps inflate whole-profile means, so regional medians and interval classes are also retained.

Bins are [r_inner,r_inner+1), [r_inner+1,4), [4,8), [8,16), [16,32]. Crossings belong to their radial bin; intervals belong to the bin containing their midpoint, including gaps spanning bins. Interval and crossing counts therefore differ. Diagnostic classes: <=3 cells grid-scale; >3 through 8 marginal; >8 well sampled spatially, not continuum-converged. The thresholds are disclosed conventions, not physical transitions.

Grid uniformity error is recorded per station. Subcell interpolation is not a subcell accuracy certificate: assuming one simple continuum zero per sign bracket, each crossing is uncertain by at most one cell and an interval by two cells (wider for zero plateaux). Without that assumption, hidden zeros and aliasing cannot be bounded. No continuum truncation estimate is possible from one grid per mu. Raw signs include tiny tails. Relative anchor floors 1e-12, 1e-8, and 1e-4 of each field's maximum are additionally recorded. Suppression by a floor establishes small amplitude, not numerical origin alone.

## Interpretation

Rapid ringing is predominantly in r >= 16, with later radial-pressure ringing entering 8 <= r < 16 in the mu=1.040 cases. It is not primarily an inner-boundary train. Isolated inner grid-scale reversals also occur. Many inner/mid-domain angular crossing separations are tens to hundreds of cells and well sampled, but this alone does not establish physical causation.

At t=3, eps=0.001, outer-bin angular-pressure median spacings are 2.22951, 2.44326, and 2.69934 cells for increasing mu. The imaginary second spinor component has corresponding medians 3.18426, 2.97485, and 2.93088 cells. These rapid scales track cell size rather than a common wavelength in r_m. Mu and dr change together, so this is not a controlled convergence test. At t=11.5, eps=0.001, mu=1.010 radial pressure has 86 crossings, all at r>=16, mean 3.54416 cells; mu=1.040 has 32, mean 3.53304 cells.

Time dependence is strong and nonmonotone; early radial pressure often never changes sign. At fixed mu, changing eps from 0.001 to 0.010 barely affects early patterns but alters late near-zero tail counts. At mu=1.010, t=12 the global radial mean is 16.4830 cells for eps=0.001 versus 3.29747 for eps=0.010: a distant inner crossing inflates the former. This is not evidence of a fivefold change in the local oscillatory scale. Case endpoints differ; all available stations are included below.

## Spinor association

`coupled_phi` is complex128 with shape (2,n,1). It has no intrinsic sign. Four real/imaginary parts are analyzed separately in the stored phase convention. A common phase rotation changes these sign crossings without necessarily changing pressure bilinears; the correlation is not gauge-invariant evidence of physical forcing.

Each station records crossing-edge binary Pearson correlation, same-edge overlap, pressure-crossing fraction within two edges of a spinor crossing, and the spinor's dilated edge coverage as a uniform-edge null baseline. The baseline is descriptive, not a p-value: localization and boundaries violate a spatially uniform independent-null assumption. The strongest consistent association is angular pressure with phi1_imag. Radial pressure has little such association; other spinor parts are generally weaker, with occasional isolated coincidences. Shared numerical spinor ringing can drive pressure ringing too. Correlation does not distinguish a physical drive from a numerical one.

The summary tables pool over saved times. Totals are repeated measurements, not unique spatial zeros; near fractions are weighted by pressure crossing count. Pooled regional means are weighted by interval count and can still include long gaps. Full crossing positions, per-station/per-region means and medians, all four spinor associations, and floor sensitivities are in the record.

## Case coverage

| Case | n | dr/r_m | saved t/r_m | radial counts total (min–max) | angular counts total (min–max) |
|---|---:|---:|---|---:|---:|
| mu-1.002-eps-0.001 | 3921 | 0.00790765 | 0–8 | 288 (0–83) | 1119 (16–142) |
| mu-1.002-eps-0.010 | 3921 | 0.00790765 | 0–8 | 288 (0–81) | 1118 (16–150) |
| mu-1.010-eps-0.001 | 1750 | 0.0177187 | 0–15.5 | 861 (0–86) | 861 (4–80) |
| mu-1.010-eps-0.010 | 1750 | 0.0177187 | 0–16 | 868 (0–84) | 874 (4–80) |
| mu-1.040-eps-0.001 | 869 | 0.0356682 | 0–18.5 | 506 (0–36) | 678 (2–51) |
| mu-1.040-eps-0.010 | 869 | 0.0356682 | 0–19.5 | 578 (0–42) | 680 (3–51) |

## Radial distribution

| Case | Field | crossing counts in five bins | pooled mean cells in bins | outer intervals: <=3 / marginal / >8 |
|---|---|---|---|---|
| mu-1.002-eps-0.001 | radial_pressure | 10, 0, 7, 1, 270 | 1.9208, 582.242, 327.919, 2224.76, 13.6932 | 60 / 136 / 65 |
| mu-1.002-eps-0.001 | angular_pressure | 6, 8, 17, 59, 1029 | 206.937, 242.927, 409.193, 271.788, 33.4145 | 650 / 233 / 136 |
| mu-1.002-eps-0.010 | radial_pressure | 10, 0, 7, 1, 270 | 1.81392, 580.757, 327.531, 2224.97, 13.7369 | 60 / 134 / 67 |
| mu-1.002-eps-0.010 | angular_pressure | 6, 8, 17, 59, 1028 | 207.109, 243.035, 409.303, 271.797, 33.46 | 645 / 236 / 137 |
| mu-1.010-eps-0.001 | radial_pressure | 19, 0, 11, 15, 816 | 1.06499, 220.752, 122.782, 554.102, 8.5625 | 256 / 495 / 48 |
| mu-1.010-eps-0.001 | angular_pressure | 7, 11, 42, 108, 693 | 92.0358, 97.7968, 149.811, 116.137, 41.2609 | 395 / 135 / 149 |
| mu-1.010-eps-0.010 | radial_pressure | 10, 0, 11, 29, 818 | 0.834007, 241.199, 93.6771, 256.672, 9.5691 | 260 / 489 / 52 |
| mu-1.010-eps-0.010 | angular_pressure | 7, 13, 43, 109, 702 | 92.0834, 89.9298, 148.766, 116.587, 44.3319 | 398 / 140 / 149 |
| mu-1.040-eps-0.001 | radial_pressure | 4, 1, 12, 103, 386 | 0.0127483, 139.464, 22.8236, 30.2007, 12.4132 | 136 / 204 / 30 |
| mu-1.040-eps-0.001 | angular_pressure | 4, 11, 53, 106, 504 | —, 48.2836, 63.4199, 58.837, 37.0718 | 257 / 83 / 149 |
| mu-1.040-eps-0.010 | radial_pressure | 2, 1, 19, 164, 392 | —, 149.851, 17.1346, 21.0472, 14.5137 | 145 / 198 / 33 |
| mu-1.040-eps-0.010 | angular_pressure | 4, 15, 54, 105, 502 | —, 45.5426, 62.6396, 58.8504, 41.8394 | 254 / 81 / 152 |

## Correlation and amplitude sensitivity

| Case | angular near phi1_imag / null | radial near phi1_imag / null | angular Pearson range | counts surviving 1e-4 floor: radial/angular |
|---|---:|---:|---:|---:|
| mu-1.002-eps-0.001 | 0.6881 / 0.0846 | 0.0208 / 0.0185 | -0.0058–0.7449 | 18/84 |
| mu-1.002-eps-0.010 | 0.7013 / 0.0860 | 0.0208 / 0.0186 | -0.0058–0.6942 | 18/84 |
| mu-1.010-eps-0.001 | 0.5273 / 0.0932 | 0.0012 / 0.0238 | -0.0117–0.6612 | 39/104 |
| mu-1.010-eps-0.010 | 0.5217 / 0.0920 | 0.0000 / 0.0235 | -0.0096–0.6649 | 40/108 |
| mu-1.040-eps-0.001 | 0.4617 / 0.1200 | 0.0020 / 0.0374 | -0.0170–0.6160 | 98/116 |
| mu-1.040-eps-0.010 | 0.4721 / 0.1195 | 0.0069 / 0.0366 | -0.0170–0.6130 | 160/122 |

## Result

**Mixed evidence: predominantly unresolved grid-scale/marginal ringing in the outer region, coexisting with well-sampled broad sign changes. Resolved physical forcing is not established.** At t=3 outer angular median spacings are 2.23–2.70 cells (eps=0.001); correlated spinor spacings are about 3 cells. The strong angular/spinor association does not demonstrate resolved physics. Radial ringing near t=11.5 is about 3.53–3.54 cells for mu=1.010 and 1.040, with many <=3-cell subintervals. The isolated mu=1.002 inner pair at t=3 is 1.9208 cells (0.0151890 r_m) for eps=0.001 and 1.81392 cells (0.0143439 r_m) for eps=0.010. No single power law or wavelength describes these fields.

Below are ranges of whole-profile mean half-wavelength across saved stations, NOT local resolution certificates. Counts are summed over stations, followed by per-station min–max. Individual station counts and wavelengths are in the owning note and creation-only record.

| Case | radial counts total (range) | radial mean half-wave range | angular counts total (range) | angular mean half-wave range |
|---|---:|---|---:|---|
| mu-1.002-eps-0.001 | 288 (0–83) | 1.9208–155.59 cells (0.015189–1.23035 r_m) | 1119 (16–142) | 27.0998–249.93 cells (0.214296–1.97636 r_m) |
| mu-1.002-eps-0.010 | 288 (0–81) | 1.81392–141.887 cells (0.0143439–1.122 r_m) | 1118 (16–150) | 25.6491–249.93 cells (0.202824–1.97636 r_m) |
| mu-1.010-eps-0.001 | 861 (0–86) | 1.54789–338.971 cells (0.0274265–6.00613 r_m) | 861 (4–80) | 21.594–392.854 cells (0.382618–6.96087 r_m) |
| mu-1.010-eps-0.010 | 868 (0–84) | 1.5277–338.902 cells (0.0270689–6.00491 r_m) | 874 (4–80) | 21.5976–404.723 cells (0.382681–7.17117 r_m) |
| mu-1.040-eps-0.001 | 506 (0–36) | 3.12143–158.964 cells (0.111336–5.66995 r_m) | 678 (2–51) | 16.4712–385.205 cells (0.587499–13.7396 r_m) |
| mu-1.040-eps-0.010 | 578 (0–42) | 3.12244–158.9 cells (0.111372–5.66769 r_m) | 680 (3–51) | 16.4745–392.301 cells (0.587615–13.9927 r_m) |

A same-mu finer-grid replay, phase-invariant spinor/bilinear diagnostics, and boundary-treatment comparison are needed to identify the mechanism. These controls were not run; no continuum error bar or physical causation is claimed.
