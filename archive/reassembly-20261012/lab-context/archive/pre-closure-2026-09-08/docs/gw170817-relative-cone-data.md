# GW170817/GRB 170817A: compact relative-cone input record

This repository stores four small public metadata products for a narrow,
auditable check of the **published** GW170817/GRB 170817A timing inputs.  The
machine record is [`published-constraints.json`](../collected-data/published-constraints.json);
the dependency-free validator is
[`observations.py`](../src/recursive_horizons/observations.py), and the generated
artifact is [`gw170817-timing.json`](../results/gw170817-timing.json).

It is deliberately **not** a raw strain or gamma-ray light-curve reanalysis.
It neither measures a raw dimensional local `c` nor compares our domain with a
parent domain.  It checks a stated low-redshift tensor/photon propagation
subclass only after the source-emission assumptions of the collaboration paper
are declared.

## Frozen sources

| Local file | Format and role | Source |
| --- | --- | --- |
| `gwosc-event-v1.json` | GWOSC event release metadata: event identity, GPS time, detectors, release DOI | [GWOSC GW170817 JSON API](https://gwosc.org/api/v2/event-versions/GW170817-v1?format=json) |
| `glg_tcat_all_bn170817529_v03.fit` | Fermi/GBM primary FITS header: `TRIGTIME`, `OBJECT`, `DATE-OBS` | [HEASARC GBM trigger file](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2017/bn170817529/current/glg_tcat_all_bn170817529_v03.fit) |
| `G298048.lvc` | LVC notice UTC trigger time | [GCN/LVC notice](https://gcn.gsfc.nasa.gov/notices_l/G298048.lvc) |
| `GCN-21520.txt` | Fermi/GBM trigger identifier and rounded UTC | [GCN Circular 21520 text](https://gcn.nasa.gov/circulars/21520.txt) |

The JSON records SHA-256 values for those exact local files.  A matching hash
shows only that the local input is the declared input; it does not validate an
astrophysical model.

The primary source for the published timing and speed interval is the joint
LIGO/Virgo/Fermi-GBM/INTEGRAL paper,
[arXiv:1710.05834v2](https://arxiv.org/abs/1710.05834v2),
DOI [10.3847/2041-8213/aa920c](https://doi.org/10.3847/2041-8213/aa920c).
GWOSC event data are released under [CC BY 4.0](https://gwosc.org/events/GW170817/).
Fermi states that GBM data were never proprietary; retain archive filenames and
version suffixes when using its public products.

## Timing reconstruction

The validator checks the following declared collaboration inputs:

```text
GW geocentric merger time     = 1187008882.430 +/- 0.002 GPS s
                              = 2017-08-17T12:41:04.430000Z
GBM trigger                   = 2017-08-17T12:41:06.474598Z
                              = 524666471.474598 Fermi MET
gamma main-pulse onset        = 0.310 +/- 0.048 s before the GBM trigger
Fermi arrival correction      = 0.003176 s before geocenter
```

Thus the reconstructed *geocentric onset* difference is

```text
Delta t = t_gamma,onset,geo - t_GW,merger,geo
        = 1.737774 s,
```

which is consistent with the reported `+1.74 +/- 0.05 s`.  Combining the
declared `0.002 s` GW and `0.048 s` onset uncertainties in quadrature gives
`0.04804... s`, which rounds to `0.05 s`.  This is a reconstruction of the
published input arithmetic.  It is not a refit of the GW waveform, a GBM onset
estimator, or a new uncertainty analysis.

The LVC notice's low-latency UTC time, `12:41:04.445710`, is intentionally kept
as a provenance cross-check; it is not substituted for the paper's geocentric
merger-time posterior.

## Conditional speed interval

For the paper's defined quantity

```text
delta_v = (v_GW - v_EM) / v_EM,
```

and a propagation-only `d = t_gamma - t_GW` over distance `D`, the compact
arithmetic uses

```text
delta_v = d / (D/c - d).
```

The collaboration used the conservative distance `D=26 Mpc`.

- Taking first photons and the GW peak to have been emitted simultaneously sets
  `d=+1.74 s` and gives `+6.50e-16`, rounded to `+7e-16`.
- Allowing the GRB signal to have been emitted `10 s` after the GW signal sets
  `d=1.74-10 s` and gives `-3.09e-15`, rounded to `-3e-15`.

This reproduces the stated interval

```text
-3e-15 <= (v_GW-v_EM)/v_EM <= +7e-16.
```

The interval is conditional on that source-lag bracket.  In general,

```text
observed arrival difference = intrinsic source emission lag + propagation delay.
```

No compact public metadata check can remove this degeneracy.  A model with a
redshift-, direction-, frequency-, or mode-dependent cone requires a separate
propagation calculation; an inaccessible parent speed is not an observable in
this data set.

## Run

```bash
python3 scripts/reproduce_observations.py
python3 -m unittest discover -s tests -p 'test_observations.py' -v
```

The test verifies source hashes, the relevant FITS header cards, LVC/GCN
trigger metadata, timing/sign conventions, the published bound rounding, and
malformed-input rejection.
