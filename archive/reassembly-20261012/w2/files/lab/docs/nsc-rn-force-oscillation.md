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

## Every recorded station

N_R/N_A: raw pressure crossing counts. H: mean crossing-to-crossing half-wavelength.

### mu-1.002-eps-0.001

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | — | — | 16 | 249.93 | 1.97636 |
| 0.5 | 0 | — | — | 41 | 94.0852 | 0.743993 |
| 1 | 0 | — | — | 66 | 58.947 | 0.466132 |
| 1.5 | 0 | — | — | 91 | 43.4166 | 0.343323 |
| 2 | 1 | — | — | 113 | 32.3674 | 0.25595 |
| 2.5 | 1 | — | — | 138 | 27.1691 | 0.214844 |
| 3 | 2 | 1.9208 | 0.015189 | 142 | 27.0998 | 0.214296 |
| 3.5 | 2 | 2.25476 | 0.0178299 | 101 | 39.1205 | 0.309351 |
| 4 | 2 | 4.01812 | 0.0317739 | 81 | 47.3761 | 0.374634 |
| 4.5 | 4 | 8.25318 | 0.0652633 | 68 | 58.202 | 0.460242 |
| 5 | 8 | 8.29369 | 0.0655836 | 56 | 52.187 | 0.412677 |
| 5.5 | 14 | 7.68412 | 0.0607633 | 47 | 63.6962 | 0.503688 |
| 6 | 21 | 155.59 | 1.23035 | 40 | 76.4508 | 0.604547 |
| 6.5 | 35 | 90.9183 | 0.71895 | 35 | 77.9121 | 0.616102 |
| 7 | 49 | 64.108 | 0.506944 | 31 | 90.4722 | 0.715423 |
| 7.5 | 65 | 47.9552 | 0.379213 | 28 | 123.248 | 0.974605 |
| 8 | 83 | 37.3982 | 0.295732 | 25 | 144.056 | 1.13914 |

### mu-1.002-eps-0.010

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | — | — | 16 | 249.93 | 1.97636 |
| 0.5 | 0 | — | — | 41 | 94.087 | 0.744007 |
| 1 | 0 | — | — | 66 | 58.9493 | 0.466151 |
| 1.5 | 0 | — | — | 91 | 43.4186 | 0.343339 |
| 2 | 1 | — | — | 113 | 32.3709 | 0.255977 |
| 2.5 | 1 | — | — | 139 | 26.9858 | 0.213394 |
| 3 | 2 | 1.81392 | 0.0143439 | 150 | 25.6491 | 0.202824 |
| 3.5 | 2 | 2.13222 | 0.0168609 | 103 | 38.3618 | 0.303352 |
| 4 | 2 | 4.04814 | 0.0320113 | 83 | 46.2327 | 0.365592 |
| 4.5 | 4 | 8.25276 | 0.06526 | 66 | 60.0152 | 0.474579 |
| 5 | 8 | 8.28081 | 0.0654818 | 56 | 52.1955 | 0.412744 |
| 5.5 | 14 | 7.67686 | 0.060706 | 47 | 63.7101 | 0.503797 |
| 6 | 23 | 141.887 | 1.122 | 41 | 74.58 | 0.589753 |
| 6.5 | 35 | 90.8917 | 0.71874 | 30 | 91.4229 | 0.722941 |
| 7 | 49 | 64.0919 | 0.506817 | 27 | 104.436 | 0.825847 |
| 7.5 | 65 | 47.9463 | 0.379142 | 24 | 144.786 | 1.14492 |
| 8 | 81 | 38.3296 | 0.303097 | 25 | 144.108 | 1.13955 |

### mu-1.010-eps-0.001

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | — | — | 16 | 111.64 | 1.97811 |
| 0.5 | 0 | — | — | 27 | 64.6439 | 1.14541 |
| 1 | 0 | — | — | 38 | 46.2492 | 0.819476 |
| 1.5 | 0 | — | — | 49 | 36.3447 | 0.643981 |
| 2 | 0 | — | — | 58 | 28.4023 | 0.503252 |
| 2.5 | 0 | — | — | 69 | 24.4459 | 0.43315 |
| 3 | 0 | — | — | 80 | 21.594 | 0.382618 |
| 3.5 | 0 | — | — | 75 | 23.5921 | 0.418021 |
| 4 | 0 | — | — | 58 | 29.7043 | 0.526322 |
| 4.5 | 0 | — | — | 49 | 36.2992 | 0.643174 |
| 5 | 0 | — | — | 41 | 32.0412 | 0.567728 |
| 5.5 | 2 | 1.54789 | 0.0274265 | 34 | 39.6419 | 0.702402 |
| 6 | 5 | 338.971 | 6.00613 | 29 | 47.5703 | 0.842883 |
| 6.5 | 9 | 166.642 | 2.95267 | 25 | 49.2853 | 0.873271 |
| 7 | 15 | 95.1315 | 1.68561 | 22 | 57.6852 | 1.02211 |
| 7.5 | 21 | 66.3697 | 1.17598 | 23 | 67.4714 | 1.19551 |
| 8 | 29 | 47.4879 | 0.841424 | 20 | 81.1153 | 1.43726 |
| 8.5 | 39 | 35.1983 | 0.623668 | 18 | 92.8403 | 1.64501 |
| 9 | 49 | 27.9925 | 0.495991 | 19 | 88.9476 | 1.57604 |
| 9.5 | 60 | 23.0108 | 0.407722 | 14 | 106.336 | 1.88413 |
| 10 | 70 | 19.7012 | 0.349079 | 13 | 120.017 | 2.12655 |
| 10.5 | 80 | 4.48471 | 0.0794633 | 13 | 124.548 | 2.20682 |
| 11 | 84 | 3.87772 | 0.0687082 | 12 | 139.235 | 2.46707 |
| 11.5 | 86 | 3.54416 | 0.0627979 | 11 | 156.163 | 2.76701 |
| 12 | 86 | 16.483 | 0.292058 | 10 | 176.197 | 3.12199 |
| 12.5 | 84 | 17.003 | 0.30127 | 9 | 175.163 | 3.10366 |
| 13 | 50 | 28.9927 | 0.513712 | 5 | 362.164 | 6.41707 |
| 13.5 | 34 | 43.4722 | 0.770271 | 6 | 296.918 | 5.261 |
| 14 | 20 | 10.6435 | 0.18859 | 5 | 379.06 | 6.71644 |
| 14.5 | 14 | 60.0744 | 1.06444 | 4 | 116.649 | 2.06687 |
| 15 | 10 | 93.7389 | 1.66093 | 5 | 392.854 | 6.96087 |
| 15.5 | 14 | 69.8462 | 1.23758 | 4 | 103.445 | 1.8329 |

### mu-1.010-eps-0.010

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | — | — | 16 | 111.64 | 1.97811 |
| 0.5 | 0 | — | — | 27 | 64.6451 | 1.14543 |
| 1 | 0 | — | — | 38 | 46.2501 | 0.819491 |
| 1.5 | 0 | — | — | 49 | 36.3471 | 0.644023 |
| 2 | 0 | — | — | 58 | 28.4058 | 0.503313 |
| 2.5 | 0 | — | — | 69 | 24.4501 | 0.433223 |
| 3 | 0 | — | — | 80 | 21.5976 | 0.382681 |
| 3.5 | 0 | — | — | 75 | 23.5963 | 0.418096 |
| 4 | 0 | — | — | 56 | 30.7917 | 0.545589 |
| 4.5 | 0 | — | — | 47 | 37.8918 | 0.671393 |
| 5 | 0 | — | — | 41 | 32.0497 | 0.567878 |
| 5.5 | 2 | 1.5277 | 0.0270689 | 36 | 37.3853 | 0.662419 |
| 6 | 5 | 338.902 | 6.00491 | 34 | 40.4074 | 0.715966 |
| 6.5 | 9 | 166.581 | 2.9516 | 25 | 49.306 | 0.873639 |
| 7 | 15 | 95.0967 | 1.68499 | 24 | 52.6909 | 0.933614 |
| 7.5 | 21 | 66.3479 | 1.1756 | 23 | 67.504 | 1.19608 |
| 8 | 29 | 47.4754 | 0.841202 | 21 | 77.1564 | 1.36711 |
| 8.5 | 37 | 36.9523 | 0.654746 | 20 | 83.0979 | 1.47239 |
| 9 | 49 | 27.9819 | 0.495803 | 17 | 100.074 | 1.77319 |
| 9.5 | 57 | 24.0843 | 0.426743 | 14 | 106.373 | 1.88479 |
| 10 | 70 | 19.6927 | 0.348929 | 14 | 110.942 | 1.96575 |
| 10.5 | 82 | 4.44956 | 0.0788404 | 13 | 124.608 | 2.20789 |
| 11 | 84 | 3.87564 | 0.0686713 | 12 | 139.311 | 2.46841 |
| 11.5 | 84 | 3.57078 | 0.0632695 | 11 | 156.233 | 2.76825 |
| 12 | 84 | 3.29747 | 0.0584268 | 11 | 158.78 | 2.81338 |
| 12.5 | 74 | 3.56543 | 0.0631747 | 8 | 200.445 | 3.55163 |
| 13 | 48 | 5.02788 | 0.0890874 | 5 | 362.375 | 6.4208 |
| 13.5 | 32 | 7.13299 | 0.126387 | 6 | 296.952 | 5.2616 |
| 14 | 20 | 10.9286 | 0.193641 | 5 | 379.279 | 6.72033 |
| 14.5 | 14 | 60.095 | 1.06481 | 5 | 386.581 | 6.84971 |
| 15 | 12 | 76.7212 | 1.3594 | 4 | 112.008 | 1.98463 |
| 15.5 | 14 | 69.8743 | 1.23808 | 5 | 399.334 | 7.07567 |
| 16 | 26 | 38.861 | 0.688567 | 5 | 404.723 | 7.17117 |

### mu-1.040-eps-0.001

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | — | — | 16 | 55.6562 | 1.98516 |
| 0.5 | 0 | — | — | 21 | 41.8211 | 1.49168 |
| 1 | 0 | — | — | 27 | 32.8057 | 1.17012 |
| 1.5 | 0 | — | — | 31 | 26.1154 | 0.93149 |
| 2 | 0 | — | — | 35 | 23.6727 | 0.844362 |
| 2.5 | 0 | — | — | 41 | 20.6782 | 0.737556 |
| 3 | 0 | — | — | 46 | 18.8535 | 0.672469 |
| 3.5 | 0 | — | — | 51 | 16.4712 | 0.587499 |
| 4 | 0 | — | — | 44 | 19.5923 | 0.698823 |
| 4.5 | 0 | — | — | 35 | 18.3478 | 0.654431 |
| 5 | 0 | — | — | 31 | 21.2656 | 0.758505 |
| 5.5 | 0 | — | — | 28 | 24.1137 | 0.860092 |
| 6 | 0 | — | — | 26 | 26.5977 | 0.948693 |
| 6.5 | 2 | 15.3597 | 0.547852 | 21 | 29.337 | 1.0464 |
| 7 | 5 | 158.964 | 5.66995 | 21 | 30.0968 | 1.0735 |
| 7.5 | 9 | 80.3496 | 2.86593 | 21 | 36.819 | 1.31327 |
| 8 | 12 | 57.9652 | 2.06751 | 18 | 44.904 | 1.60164 |
| 8.5 | 14 | 41.6841 | 1.4868 | 18 | 46.008 | 1.64102 |
| 9 | 20 | 29.813 | 1.06337 | 15 | 47.5727 | 1.69683 |
| 9.5 | 22 | 5.0456 | 0.179968 | 13 | 57.2944 | 2.04359 |
| 10 | 28 | 4.74347 | 0.169191 | 13 | 59.4953 | 2.12209 |
| 10.5 | 32 | 4.27553 | 0.1525 | 12 | 67.1466 | 2.395 |
| 11 | 32 | 3.80956 | 0.13588 | 12 | 68.9459 | 2.45918 |
| 11.5 | 32 | 3.53304 | 0.126017 | 11 | 77.3147 | 2.75768 |
| 12 | 32 | 3.27426 | 0.116787 | 11 | 78.5934 | 2.80329 |
| 12.5 | 30 | 3.12143 | 0.111336 | 7 | 115.684 | 4.12622 |
| 13 | 24 | 3.7455 | 0.133595 | 6 | 143.117 | 5.10474 |
| 13.5 | 14 | 5.89178 | 0.210149 | 6 | 147.06 | 5.24537 |
| 14 | 10 | 8.12842 | 0.289926 | 5 | 187.881 | 6.70139 |
| 14.5 | 6 | 13.4081 | 0.478242 | 5 | 191.551 | 6.83229 |
| 15 | 6 | 79.5831 | 2.83859 | 4 | 54.8868 | 1.95771 |
| 15.5 | 12 | 39.2201 | 1.39891 | 5 | 197.98 | 7.06159 |
| 16 | 16 | 30.9191 | 1.10283 | 5 | 200.485 | 7.15095 |
| 16.5 | 22 | 23.5735 | 0.840826 | 3 | 363.925 | 12.9805 |
| 17 | 26 | 21.0239 | 0.749885 | 2 | 43.3031 | 1.54454 |
| 17.5 | 30 | 19.1404 | 0.682705 | 4 | 252.364 | 9.00138 |
| 18 | 34 | 17.6969 | 0.631215 | 3 | 385.205 | 13.7396 |
| 18.5 | 36 | 17.4047 | 0.620796 | 5 | 195.936 | 6.9887 |

### mu-1.040-eps-0.010

| t/r_m | N_R | H_R cells | H_R/r_m | N_A | H_A cells | H_A/r_m |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | — | — | 16 | 55.6562 | 1.98516 |
| 0.5 | 0 | — | — | 21 | 41.8217 | 1.4917 |
| 1 | 0 | — | — | 27 | 32.8072 | 1.17017 |
| 1.5 | 0 | — | — | 32 | 25.3177 | 0.903036 |
| 2 | 0 | — | — | 35 | 23.6754 | 0.844461 |
| 2.5 | 0 | — | — | 41 | 20.6817 | 0.73768 |
| 3 | 0 | — | — | 46 | 18.8567 | 0.672585 |
| 3.5 | 0 | — | — | 51 | 16.4745 | 0.587615 |
| 4 | 0 | — | — | 43 | 20.0948 | 0.716745 |
| 4.5 | 0 | — | — | 35 | 18.3521 | 0.654585 |
| 5 | 0 | — | — | 31 | 21.2718 | 0.758725 |
| 5.5 | 0 | — | — | 28 | 24.1194 | 0.860295 |
| 6 | 0 | — | — | 24 | 28.9193 | 1.0315 |
| 6.5 | 2 | 15.403 | 0.549398 | 22 | 28.0104 | 0.99908 |
| 7 | 5 | 158.9 | 5.66769 | 19 | 33.4539 | 1.19324 |
| 7.5 | 9 | 80.3578 | 2.86622 | 21 | 36.8427 | 1.31411 |
| 8 | 10 | 57.7718 | 2.06062 | 18 | 44.9237 | 1.60235 |
| 8.5 | 14 | 41.7043 | 1.48752 | 18 | 46.028 | 1.64174 |
| 9 | 20 | 29.8407 | 1.06436 | 14 | 51.3536 | 1.83169 |
| 9.5 | 22 | 5.05247 | 0.180213 | 13 | 57.3279 | 2.04478 |
| 10 | 28 | 4.74197 | 0.169137 | 13 | 59.5401 | 2.12369 |
| 10.5 | 32 | 4.27895 | 0.152623 | 12 | 67.1738 | 2.39597 |
| 11 | 32 | 3.81668 | 0.136134 | 12 | 68.9826 | 2.46048 |
| 11.5 | 32 | 3.53561 | 0.126109 | 12 | 70.4147 | 2.51157 |
| 12 | 32 | 3.26958 | 0.11662 | 9 | 98.277 | 3.50536 |
| 12.5 | 30 | 3.12244 | 0.111372 | 5 | 173.611 | 6.19241 |
| 13 | 20 | 4.52128 | 0.161266 | 6 | 143.217 | 5.10829 |
| 13.5 | 14 | 5.9033 | 0.21056 | 6 | 147.138 | 5.24813 |
| 14 | 8 | 10.5013 | 0.374563 | 5 | 187.994 | 6.70541 |
| 14.5 | 6 | 13.4074 | 0.478218 | 5 | 191.696 | 6.83746 |
| 15 | 6 | 79.6221 | 2.83998 | 5 | 195.113 | 6.95933 |
| 15.5 | 12 | 39.2291 | 1.39923 | 4 | 51.4724 | 1.83593 |
| 16 | 14 | 35.6859 | 1.27285 | 5 | 200.711 | 7.15899 |
| 16.5 | 22 | 23.5902 | 0.84142 | 3 | 364.23 | 12.9914 |
| 17 | 28 | 19.4689 | 0.69442 | 4 | 247.351 | 8.82258 |
| 17.5 | 30 | 19.1364 | 0.68256 | 4 | 252.578 | 9.00899 |
| 18 | 32 | 18.7407 | 0.668447 | 3 | 385.448 | 13.7482 |
| 18.5 | 38 | 16.4772 | 0.58771 | 3 | 392.301 | 13.9927 |
| 19 | 38 | 17.2327 | 0.614658 | 5 | 199.361 | 7.11085 |
| 19.5 | 42 | 15.5275 | 0.553836 | 4 | 269.587 | 9.61567 |

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
