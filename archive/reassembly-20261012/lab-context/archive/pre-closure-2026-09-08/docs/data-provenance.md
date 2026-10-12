# Data provenance

[`collected-data/manifest.csv`](../collected-data/manifest.csv) is the machine-readable authority for local third-party artifacts. Every row records its repository-relative filename, upstream release/source, bytes, SHA-256, retrieval date, citation, redistribution note, and `required_in_clone` status.

```bash
python3 scripts/verify_data.py
python3 scripts/verify_data.py --require
```

Default verification requires every small tracked artifact and permits the two large ignored Planck products to be absent. `--require` also requires those optional Planck files. Every present entry receives exact byte-size and SHA-256 verification; only appropriate format/content identity checks are then applied.

## DESI DR2 BAO

`collected-data/desi-dr2-bao/` contains the 13-element DR2 BAO mean vector, matching 13 by 13 covariance, and small official `base_w/desi-bao-all` reference outputs. The likelihood inputs are pinned to [CobayaSampler/bao_data commit `b7b8a36e9bccb063081f811f323cada21ab5fbdd`](https://github.com/CobayaSampler/bao_data/commit/b7b8a36e9bccb063081f811f323cada21ab5fbdd), which the official [DESI DR2 BAO results README](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html) names as the public likelihood source. The three release outputs are checked against the official [v1.0 SHA-256 manifest](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum).

This is a Gaussian likelihood of measured **BAO distances**, not a Gaussian posterior summary of `w`. It is BAO-only, not a full-shape likelihood, and fitting a phenomenological scaling cannot by itself establish an external boundary, variable local `c`, or a recursive-domain mechanism. Cite [DESI DR2 Results II](https://arxiv.org/abs/2503.14738). DESI data are [CC BY 4.0](https://data.desi.lbl.gov/doc/acknowledgments/) with citation, modification-notice, and acknowledgment requirements.

## GW170817 provenance

`collected-data/gw170817/` holds compact provenance records: [GWOSC event metadata](https://gwosc.org/api/v2/event-versions/GW170817-v1?format=json), [GCN/LVC G298048 notices](https://gcn.gsfc.nasa.gov/notices_l/G298048.lvc), [GCN Circular 21520 text](https://gcn.nasa.gov/circulars/21520.txt), and the [Fermi GBM trigger-entry FITS](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2017/bn170817529/current/glg_tcat_all_bn170817529_v03.fit). Its Fermi signature requires `SIMPLE=T`, `FILETYPE=TRIGGER ENTRY`, `TELESCOP=GLAST`, `INSTRUME=GBM`, `OBJECT=GRB170817529`, `TRIGTIME=524666471.474598`, and the manifest’s 5,760 bytes.

These files establish event association and support reconstruction of the published timing-input arithmetic, not an independent propagation inference. No raw GW strain, calibrated gamma-ray time series, detector response, source-emission model, or line-of-sight likelihood is included. They cannot remove the intrinsic emission-lag assumption behind the GW170817 speed interval; see [Abbott et al. 2017](https://arxiv.org/abs/1710.05834). GWOSC data are CC BY 4.0; cite [the GW170817 dataset DOI](https://doi.org/10.7935/K5B8566F). The public Fermi/GCN records are attributed to their sources; no separate Creative Commons license is asserted here.

## Deliberate data boundary

Raw DESI spectra/clustering products and full posterior chains are omitted: this repository does not reproduce target selection, reconstruction, BAO fitting, full-shape analysis, or a complete posterior pipeline. The retained vector/covariance are sufficient for a specified BAO-distance test, while the many-tens-of-megabytes official chains are comparison benchmarks rather than likelihood inputs. The raw GW strain is omitted for the analogous reason: event notices cannot replace calibrated detector analysis.

The ignored Planck products remain Planck 2018 **PR3**, not PR4/NPIPE. Their hash and Planck-specific HEALPix extension signatures are checked only if locally present; no active repository claim depends on them.

When present, the verifier requires the documented HEALPix binary-table identity: `PIXTYPE=HEALPIX`, `COORDSYS=GALACTIC`, `ORDERING=NESTED`, `NSIDE=2048`, and `NAXIS2=50,331,648`. The SMICA extension additionally requires `COMP-MAP`, ten 40-byte float32 fields, `I_STOKES`, `K_CMB`, and data offset 8,640 bytes; the common mask requires `MASK-INT`, `TMASK`, one four-byte field, and data offset 5,760 bytes. These checks establish local file identity/layout, not the validity of any cosmological estimator.

The previously archived patch-based quadrupole script remains invalid as a scientific result: its pixel reorder was wrong, its patch power-law slope was not a validated estimator of primordial `n_s`, and it lacked transfer, covariance, and scan-calibrated null simulations. Any replacement needs a release-matched map/likelihood, beam/pixel/mask/noise treatment, repeated optimization in every null realization, and a frozen analysis choice before unblinding.
