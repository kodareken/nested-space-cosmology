# Tests

Thirteen reproducible scripts. Each maps to a specific claim in the corpus. Run them in the order below; the de Sitter benchmark certifies the pipeline, the audit closes the instruments, and the forward test is the pre-registered experiment.

| # | Script | Claim it tests | Result | Dependencies |
|---|--------|----------------|--------|--------------|
| 1 | `consistency_audit.py` | All instruments mutually consistent (15-digit identities) + the data checks | **8/9 close** (the spectrum named open) | numpy only |
| 2 | `dS_benchmark.py` | The mode pipeline reproduces the exact de Sitter answer | **CERTIFIED**: n_s = 1.0002 vs 1.000000 | numpy only |
| 3 | `substitution_audit.py` | The closures survive non-Planck inputs (five DESI DR2 datasets) | **PASS**: 0.75–2.0σ (the known dataset tension) | numpy only |
| 4 | `junction_prediction.py` | The pre-registered quadrupolar-tilt prediction (frozen before the map test) | Prediction: A_Q/(1−n_s) = a*², axis in the anomaly family | numpy only |
| 5 | `quadrupole_test.py` | The executed forward test on the Planck SMICA map | **NOT CONFIRMED** at current sensitivity (ratio 0.154; marginal ~2σ) | numpy + healpy + the maps (see `collected-data/README.md`) |
| 6 | `junction_ns.py` | The anisotropic crossing mechanism | n_s = 1 − 12.24(σ/θ)²; observed 0.9649 ⟺ σ/θ ≈ 5.4% | numpy only |
| 7 | `spectrum_benchmarks.py` | The stiff-bounce spectrum + the setting of c | n_s = 2.94 (blue, own patch); braking 72.8 e-folds | numpy only |
| 8 | `ks_background.py` | The KS+radiation+spin background (Einstein constraint) | conserved to 10⁻¹²–10⁻¹⁴; closed oscillator, pancake string | numpy only |
| 9 | `ks_matter.py` | The KS+matter+spin background (bounce condition) | bounce requires ε_S > ε_m/2; satisfied | numpy only |
| 10 | `ks_pk.py` | The KS direction-resolved test-scalar spectra | n_s = 2.77 (R), 2.81 (Ω) — blue, converged | numpy only |
| 11 | `p1_scan.py` | The fertility scan (12⁴ points over α, m_e/m_p, G, Λ) | **UNSUPPORTED as computed**: p = 0.64, no local maximum | numpy only |
| 12 | `spin_ladder.py` | The spin-per-mass ladder, the implosion-settle temperature ladder, the spin floor | parent a* = 0.39 ± 0.04; T_GH/T_H = 2.4096 | numpy only |
| 13 | `g9_duality.py` | The deflation–inflation duality and the chain numbers | two-clock locking 2.68 = O(1); descendant stack 5.88 R_H | numpy only |

## Reproducing the forward test (#5)

```bash
python3 -m venv cmbvenv && cmbvenv/bin/pip install healpy numpy
cmbvenv/bin/python quadrupole_test.py
```

The maps: download into `collected-data/` per `collected-data/README.md` (2.0 GB SMICA + 201 MB common mask, both from the public IRSA mirror).

## Honesty rules

- The pre-registered prediction (#4) was frozen before the map test (#5) ran — the verdict is what it is.
- Scripts that were superseded during development (the self-made HEALPix geometry, the debugging iterations) are in `../archive/`, not here.
- Every number cited in the corpus is reproducible from these scripts alone (plus the public data files).
