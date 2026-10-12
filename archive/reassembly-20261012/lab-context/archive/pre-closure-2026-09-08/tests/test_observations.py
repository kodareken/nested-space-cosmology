from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.observations import (  # noqa: E402
    ObservationalInputError,
    combined_timing_uncertainty_seconds,
    file_sha256,
    load_published_constraints,
    parse_fits_primary_header,
    parse_gcn_21520_trigger,
    parse_lvc_notice_trigger_utc,
    propagation_delta_from_gamma_minus_gw_delay,
    reproduce_gw170817_timing,
    reconstructed_geocentric_delay_seconds,
    validate_gw170817_inputs,
)


DATA_DIRECTORY = REPOSITORY / "collected-data" / "gw170817"
CONSTRAINTS = REPOSITORY / "collected-data" / "published-constraints.json"


class PublishedConstraintTests(unittest.TestCase):
    def test_hashes_and_primary_source_fields_are_frozen(self) -> None:
        payload = load_published_constraints(CONSTRAINTS)
        record = payload["GW170817_GRB170817A"]
        self.assertEqual(record["primary_paper"]["arxiv"], "1710.05834v2")
        self.assertEqual(record["primary_paper"]["doi"], "10.3847/2041-8213/aa920c")
        for filename, metadata in record["local_sources"].items():
            self.assertEqual(file_sha256(DATA_DIRECTORY / filename), metadata["sha256"])

    def test_fits_notice_and_gcn_metadata_parse(self) -> None:
        header = parse_fits_primary_header(DATA_DIRECTORY / "glg_tcat_all_bn170817529_v03.fit")
        self.assertEqual(header["OBJECT"], "GRB170817529")
        self.assertEqual(header["DATE-OBS"], "2017-08-17T12:38:54")
        self.assertAlmostEqual(header["TRIGTIME"], 524666471.474598)
        self.assertEqual(
            parse_lvc_notice_trigger_utc((DATA_DIRECTORY / "G298048.lvc").read_text()),
            "12:41:04.445710",
        )
        self.assertEqual(
            parse_gcn_21520_trigger((DATA_DIRECTORY / "GCN-21520.txt").read_text()),
            {"utc": "12:41:06.47", "met_integer": "524666471", "trigger_name": "170817529"},
        )

    def test_published_input_arithmetic_and_conservative_rounding(self) -> None:
        result = validate_gw170817_inputs(DATA_DIRECTORY, CONSTRAINTS)
        self.assertAlmostEqual(result["reconstructed_geocentric_delay_seconds"], 1.737774)
        self.assertAlmostEqual(result["combined_input_uncertainty_seconds"], math.sqrt(0.002**2 + 0.048**2))
        self.assertAlmostEqual(round(result["combined_input_uncertainty_seconds"], 2), 0.05)
        self.assertEqual(result["published_rounded_delta_v_over_v_em_interval"], [-3e-15, 7e-16])
        self.assertLess(result["reconstructed_delta_v_over_v_em"]["lower"], -3e-15)
        self.assertGreater(result["reconstructed_delta_v_over_v_em"]["upper"], 6e-16)
        json.dumps(result, allow_nan=False)

    def test_versioned_timing_artifact_is_json_safe(self) -> None:
        artifact = reproduce_gw170817_timing(DATA_DIRECTORY, CONSTRAINTS)
        self.assertEqual(artifact["project_version"], "0.11.0")
        self.assertEqual(artifact["schema_version"], 1)
        self.assertEqual(
            artifact["classification"],
            "published_timing_input_reconstruction_not_raw_reanalysis",
        )
        self.assertEqual(
            artifact["validation"]["published_rounded_delta_v_over_v_em_interval"],
            [-3e-15, 7e-16],
        )
        json.dumps(artifact, allow_nan=False)

    def test_propagation_conversion_requires_an_explicit_lag_model(self) -> None:
        self.assertTrue(
            math.isclose(
                propagation_delta_from_gamma_minus_gw_delay(1.74, 26.0),
                6.501986418812904e-16,
                rel_tol=1.0e-12,
            )
        )
        self.assertTrue(
            math.isclose(
                propagation_delta_from_gamma_minus_gw_delay(1.74 - 10.0, 26.0),
                -3.08657516203416e-15,
                rel_tol=1.0e-12,
            )
        )


class MalformedInputTests(unittest.TestCase):
    def test_malformed_fits_lvc_gcn_and_times_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.fit"
            path.write_bytes(b"SIMPLE  =                    T".ljust(80, b" "))
            with self.assertRaises(ObservationalInputError):
                parse_fits_primary_header(path)
        with self.assertRaises(ObservationalInputError):
            parse_lvc_notice_trigger_utc(
                "TRIGGER_TIME: 1 SOD {12:00:00.000001} UT\n"
                "TRIGGER_TIME: 2 SOD {12:00:01.000001} UT"
            )
        with self.assertRaises(ObservationalInputError):
            parse_gcn_21520_trigger("At 12:41:06.47 UT")
        with self.assertRaises(ObservationalInputError):
            reconstructed_geocentric_delay_seconds("not-time", "2017-08-17T12:41:04Z", 0.1, 0.1)

    def test_nonphysical_numeric_inputs_are_rejected(self) -> None:
        for bad in (0.0, -1.0, True, math.inf, math.nan):
            with self.subTest(bad=repr(bad)):
                with self.assertRaises(ObservationalInputError):
                    combined_timing_uncertainty_seconds(bad, 0.1)
                with self.assertRaises(ObservationalInputError):
                    propagation_delta_from_gamma_minus_gw_delay(1.0, bad)
