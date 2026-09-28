# Collected data

This directory is a small, reproducible input bundle—not a substitute for a complete DESI, Planck, Fermi, or gravitational-wave analysis. Every artifact has an exact local name, upstream URL, byte count, SHA-256, source citation, and redistribution note in [manifest.csv](manifest.csv).

## Tracked small artifacts

| Directory | Contents | Why retained |
|---|---|---|
| `desi-dr2-bao/` | 13-point DR2 BAO mean vector, 13 by 13 covariance, and small official constant-`w` configuration/summary/MAP outputs | Supports a real BAO-distance likelihood and records the released reference setup without vendoring posterior chains. |
| `gw170817/` | GWOSC event metadata, GCN/LVC notices, the Fermi GBM GCN circular, and one 5,760-byte GBM trigger-entry FITS | Records public event-time and detector-product provenance for the relative-cone benchmark. It is not strain or gamma-ray light-curve analysis input. |

The DESI likelihood vector and covariance are pinned to commit `b7b8a36e9bccb063081f811f323cada21ab5fbdd`, which the official [DESI DR2 BAO cosmology-results README](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/README.html) identifies as the public likelihood source. The three DESI release outputs are authenticated by the official [v1.0 SHA-256 manifest](https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum).

## Optional Planck downloads

The two large Planck 2018 PR3 FITS files are ignored by Git and optional in a fresh clone:

| File | Purpose | Bytes |
|---|---|---:|
| `smica_2048.fits` | Planck 2018 PR3 SMICA I/Q/U map | 2,013,312,960 |
| `mask_common.fits` | Planck 2018 PR3 common intensity mask | 201,335,040 |

Download them only for a release-matched, validated CMB pipeline:

```bash
curl --fail --location --continue-at - --output smica_2048.fits \
  "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/COM_CMB_IQU-smica_2048_R3.00_full.fits"
curl --fail --location --continue-at - --output mask_common.fits \
  "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/masks/COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"
```

## Verify

From the repository root:

```bash
python3 scripts/verify_data.py
```

By default, all small artifacts marked `required_in_clone=yes` must be present and validate. Missing ignored Planck maps are reported as optional skips, so a clean clone passes. To require every manifest entry, including those two local Planck files:

```bash
python3 scripts/verify_data.py --require
```

The verifier performs byte-size and SHA-256 checks for every present manifest file, then performs only the format-specific identity checks that are meaningful for that file: Planck extension headers, the Fermi trigger-entry primary header, DESI vector/covariance/configuration markers, and GWOSC/GCN record markers. A pass proves the local artifact matches its recorded source identity; it does not validate a cosmological inference.

## Deliberate omissions

- **Raw GW strain:** omitted because an arrival-time/reconstruction analysis needs calibrated strain, detector data-quality choices, and a dedicated GW parameter-estimation pipeline. Event metadata and notices cannot replace that work.
- **Full DESI chains:** omitted because the official four-chain posterior releases are tens of megabytes per data combination and are posterior benchmarks, not new likelihood inputs for a boundary-scaling model. The retained vector/covariance are the minimal testable BAO data.
- **Raw DESI spectra and clustering catalogs:** omitted because this repository does not reproduce DESI target selection, reconstruction, BAO fitting, or full-shape analysis. The DR2 BAO likelihood is the appropriate public data-product boundary here.

## Terms and citations

DESI data are [CC BY 4.0](https://data.desi.lbl.gov/doc/acknowledgments/): cite the applicable DESI release/paper, indicate modifications on redistribution, and include DESI’s required acknowledgment. GWOSC data are CC BY 4.0; cite the [GWOSC GW170817 dataset](https://doi.org/10.7935/K5B8566F). The Fermi GBM product and GCN materials are kept with source attribution; this repository does not assert a Creative Commons license for them. Consult each recorded source and citation before redistributing or publishing derived work.
