#!/usr/bin/env python3
"""Verify collected-data hashes, sizes, and format-specific identity markers."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import BinaryIO

REPOSITORY = Path(__file__).resolve().parents[1]
MANIFEST = REPOSITORY / "collected-data" / "manifest.csv"

PLANCK_FITS_SIGNATURES = {
    "smica_2048.fits": {"XTENSION": "BINTABLE", "PIXTYPE": "HEALPIX", "ORDERING": "NESTED", "COORDSYS": "GALACTIC", "NSIDE": 2048, "NAXIS2": 50_331_648, "TFIELDS": 10, "NAXIS1": 40, "EXT-NAME": "COMP-MAP", "TTYPE1": "I_STOKES", "TUNIT1": "K_CMB", "data_offset": 8_640},
    "mask_common.fits": {"XTENSION": "BINTABLE", "PIXTYPE": "HEALPIX", "ORDERING": "NESTED", "COORDSYS": "GALACTIC", "NSIDE": 2048, "NAXIS2": 50_331_648, "TFIELDS": 1, "NAXIS1": 4, "EXT-NAME": "MASK-INT", "TTYPE1": "TMASK", "data_offset": 5_760},
}
FERMI_TRIGGER_SIGNATURE = {"SIMPLE": True, "FILETYPE": "TRIGGER ENTRY", "TELESCOP": "GLAST", "INSTRUME": "GBM", "OBJECT": "GRB170817529", "TRIGTIME": 524_666_471.474598}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=REPOSITORY / "collected-data", help="directory containing the manifest filenames")
    parser.add_argument("--require", action="store_true", help="also require optional manifest files such as the ignored Planck maps")
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _parse_value(card: str) -> object:
    value = card[10:80].split("/", 1)[0].strip()
    if value.startswith("'"):
        return value.strip().strip("'").strip()
    if value == "T":
        return True
    if value == "F":
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value.replace("D", "E"))
        except ValueError:
            return value


def read_header(handle: BinaryIO, offset: int) -> tuple[dict[str, object], int]:
    handle.seek(offset)
    cards: list[str] = []
    while True:
        block = handle.read(2_880)
        if len(block) != 2_880:
            raise ValueError(f"incomplete FITS header block at byte {handle.tell()}")
        block_cards = [block.decode("ascii")[index:index + 80] for index in range(0, 2_880, 80)]
        cards.extend(block_cards)
        if any(card.startswith("END") for card in block_cards):
            break
    header: dict[str, object] = {}
    for card in cards:
        keyword = card[:8].strip()
        if keyword == "END":
            break
        if keyword and card[8:10] == "= ":
            header[keyword] = _parse_value(card)
    return header, len(cards) * 80


def padded_data_size(header: dict[str, object]) -> int:
    naxis = int(header.get("NAXIS", 0))
    if naxis == 0:
        return 0
    elements = 1
    for axis in range(1, naxis + 1):
        elements *= int(header[f"NAXIS{axis}"])
    raw = elements * (abs(int(header["BITPIX"])) // 8) + int(header.get("PCOUNT", 0))
    raw *= int(header.get("GCOUNT", 1))
    return int(math.ceil(raw / 2_880.0) * 2_880)


def fits_extension(path: Path) -> tuple[dict[str, object], int]:
    with path.open("rb") as handle:
        primary, primary_header_size = read_header(handle, 0)
        if primary.get("SIMPLE") is not True:
            raise ValueError("primary HDU does not declare SIMPLE = T")
        extension_offset = primary_header_size + padded_data_size(primary)
        extension, extension_header_size = read_header(handle, extension_offset)
    return extension, extension_offset + extension_header_size


def _check_signature(actual: dict[str, object], expected: dict[str, object], label: str) -> list[str]:
    return [f"{label}: {key} is {actual.get(key)!r}, expected {value!r}" for key, value in expected.items() if actual.get(key) != value]


def verify_fits(path: Path) -> tuple[str | None, list[str]]:
    if path.name in PLANCK_FITS_SIGNATURES:
        try:
            extension, data_offset = fits_extension(path)
        except (OSError, UnicodeDecodeError, ValueError, KeyError) as exc:
            return "Planck FITS extension", [f"{path}: cannot parse FITS headers: {exc}"]
        extension["data_offset"] = data_offset
        return "Planck FITS extension", _check_signature(extension, PLANCK_FITS_SIGNATURES[path.name], str(path))
    if path.name == "glg_tcat_all_bn170817529_v03.fit":
        try:
            with path.open("rb") as handle:
                primary, _ = read_header(handle, 0)
        except (OSError, UnicodeDecodeError, ValueError) as exc:
            return "Fermi GBM trigger-entry FITS", [f"{path}: cannot parse primary FITS header: {exc}"]
        return "Fermi GBM trigger-entry FITS", _check_signature(primary, FERMI_TRIGGER_SIGNATURE, str(path))
    return None, []


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def verify_semantics(path: Path) -> tuple[str | None, list[str]]:
    """Return narrow format/content checks for tracked small artifacts."""
    try:
        if path.name == "desi_gaussian_bao_ALL_GCcomb_mean.txt":
            rows = [line.split() for line in _text(path).splitlines() if not line.startswith("#")]
            labels = ["DV_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DM_over_rs", "DH_over_rs", "DH_over_rs", "DM_over_rs"]
            return "DESI DR2 BAO 13-element mean vector", [] if len(rows) == 13 and [row[2] for row in rows] == labels else [f"{path}: expected 13 ordered DESI BAO observables"]
        if path.name == "desi_gaussian_bao_ALL_GCcomb_cov.txt":
            rows = [[float(value) for value in line.split()] for line in _text(path).splitlines()]
            valid = len(rows) == 13 and all(len(row) == 13 for row in rows) and all(math.isfinite(value) for row in rows for value in row) and all(abs(rows[i][j] - rows[j][i]) < 1e-12 for i in range(13) for j in range(13))
            return "DESI DR2 BAO 13 by 13 covariance", [] if valid else [f"{path}: expected finite symmetric 13 by 13 covariance"]
        if path.name == "chain.updated.yaml":
            text = _text(path)
            markers = ("camb:", "version: 1.5.4", "dark_energy_model: ppf", "desi_bao_all:", "desi_gaussian_bao_ALL_GCcomb_mean.txt", "desi_gaussian_bao_ALL_GCcomb_cov.txt", "wa:\n    latex:", "value: 0.0")
            return "DESI constant-w Cobaya configuration", [] if all(marker in text for marker in markers) else [f"{path}: expected DESI constant-w CAMB/Cobaya markers"]
        if path.name == "chain.margestats":
            text = _text(path)
            return "DESI constant-w posterior summary", [] if all(marker in text for marker in ("parameter", "w", "hrdrag", "omegam")) else [f"{path}: expected DESI posterior-summary columns"]
        if path.name == "bestfit.minimum.txt":
            text = _text(path)
            return "DESI constant-w MAP summary", [] if all(marker in text for marker in ("minuslogpost", "chi2__BAO", "w", "hrdrag")) else [f"{path}: expected DESI MAP columns"]
        if path.name == "gwosc-event-v1.json":
            record = json.loads(_text(path))
            valid = record.get("name") == "GW170817" and record.get("grace_id") == "G298048" and record.get("gps") == 1_187_008_882.4 and record.get("detectors") == ["H1", "L1", "V1"]
            return "GWOSC GW170817 event metadata", [] if valid else [f"{path}: expected GW170817/G298048 metadata"]
        if path.name == "G298048.lvc":
            text = _text(path)
            return "GCN/LVC GW170817 notices", [] if all(marker in text for marker in ("GCN/LVC NOTICE", "TRIGGER_NUM:      G298048", "SEQUENCE_NUM:     2")) else [f"{path}: expected G298048 initial/update LVC notices"]
        if path.name == "GCN-21520.txt":
            text = _text(path)
            return "GCN Circular 21520", [] if all(marker in text for marker in ("NUMBER:  21520", "GRB 170817A", "trigger 524666471")) else [f"{path}: expected GCN 21520 Fermi GBM notice"]
    except (OSError, UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        return None, [f"{path}: cannot parse semantic content: {exc}"]
    return None, []


def _is_required(row: dict[str, str], require_all: bool) -> bool:
    return require_all or row.get("required_in_clone", "no").strip().lower() == "yes"


def verify_manifest(data_dir: Path, require_all: bool = False, manifest: Path = MANIFEST) -> tuple[list[str], list[str], int, int]:
    """Verify a manifest and return ``(passes, failures, verified, skipped)``."""
    with manifest.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required_columns = {"local_filename", "bytes", "sha256", "source_url", "required_in_clone"}
    if not rows or not required_columns.issubset(rows[0]):
        raise ValueError("manifest lacks required provenance columns")
    passes: list[str] = []
    failures: list[str] = []
    verified = skipped = 0
    for row in rows:
        relative = Path(row["local_filename"])
        if relative.is_absolute() or ".." in relative.parts:
            failures.append(f"{row['local_filename']}: unsafe manifest path")
            continue
        path = data_dir / relative
        if not path.exists():
            skipped += 1
            qualifier = "required" if _is_required(row, require_all) else "optional"
            passes.append(f"SKIP {row['local_filename']}: absent ({qualifier}); download from {row['source_url']}")
            if _is_required(row, require_all):
                failures.append(f"{row['local_filename']}: required file is absent")
            continue
        actual_size = path.stat().st_size
        if actual_size != int(row["bytes"]):
            failures.append(f"{row['local_filename']}: size {actual_size}, expected {row['bytes']}")
            continue
        actual_hash = sha256_file(path)
        if actual_hash != row["sha256"]:
            failures.append(f"{row['local_filename']}: SHA-256 {actual_hash}, expected {row['sha256']}")
            continue
        fits_label, fits_failures = verify_fits(path)
        semantic_label, semantic_failures = verify_semantics(path)
        failures.extend(fits_failures)
        failures.extend(semantic_failures)
        if fits_failures or semantic_failures:
            continue
        labels = [label for label in (fits_label, semantic_label) if label]
        check = "; ".join(labels) if labels else "size and SHA-256"
        passes.append(f"PASS {row['local_filename']}: {actual_size} bytes, SHA-256 {actual_hash}; {check}")
        verified += 1
    return passes, failures, verified, skipped


def main() -> int:
    args = parse_args()
    try:
        passes, failures, verified, skipped = verify_manifest(args.data_dir, args.require)
    except (OSError, ValueError, csv.Error) as exc:
        print(f"FAIL manifest: {exc}", file=sys.stderr)
        return 1
    for message in passes:
        print(message)
    for failure in failures:
        print(f"FAIL {failure}", file=sys.stderr)
    if failures:
        return 1
    print(f"data verification complete: {verified} verified, {skipped} absent/skipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
