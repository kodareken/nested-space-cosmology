"""Focused synthetic and sealed-family controls for the LOC2 runner."""

from __future__ import annotations

import ast
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]

from scripts import run_fgc_tdg9_loc2 as runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_loc2_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import build_static_gr0_shells  # noqa: E402


@dataclass(frozen=True)
class _SealedFailure:
    classification: str


class LOC2RunnerTests(unittest.TestCase):
    def test_no_forbidden_runner_or_candidate_imports(self) -> None:
        source = (ROOT / "scripts/run_fgc_tdg9_loc2.py").read_text()
        imported = {
            alias.name
            for node in ast.walk(ast.parse(source))
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertFalse(any("run_fgc_tdg9_ar1" in name for name in imported))
        self.assertFalse(any("tdg8_persisted_retry_replay" in name for name in imported))
        self.assertFalse(any("tdg9_ar1_pref1_binder" in name for name in imported))
        self.assertNotIn("acquire_writer", source)
        self.assertNotIn("publish_checkpoint", source)

    def test_digest_cancellation_is_rejected(self) -> None:
        evidence = SimpleNamespace(
            polynomial_count=2,
            candidate_count=4,
            classification="nonunique_or_interval_inconclusive",
            candidates=(),
            global_absolute_lower=Fraction(0),
            global_absolute_upper=Fraction(0),
            stationary_count_stream_sha256="b" * 64,
        )
        with self.assertRaisesRegex(runner.LOC2RunnerError, "independent_disagreement"):
            runner._require_route_agreement(evidence, evidence, "a" * 64, "digest-cancellation")  # type: ignore[attr-defined]

    def test_candidate_key_and_interval_mutations_are_rejected(self) -> None:
        cubics = [
            runner.primary.LocalCubic((0, Fraction(1, 10), Fraction(-1, 2), Fraction(1, 3)), {"case": "two-root"})
        ]
        independent_cubics = [
            runner.independent.IndependentLocalCubicV2(item.coefficients, item.metadata)
            for item in cubics
        ]
        first = runner.primary.localize_absolute_maximum(
            cubics, maximum_candidates=4, refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH
        )
        second = runner.independent.localize_absolute_maximum_independently_v2(
            independent_cubics, maximum_candidates=4
        )
        digest = runner._primary_stationary_count_digest(cubics)  # type: ignore[attr-defined]
        runner._require_route_agreement(first, second, digest, "control")  # type: ignore[attr-defined]
        key_changed = replace(
            second,
            candidates=(replace(second.candidates[0], location_ordinal=7), *second.candidates[1:]),
        )
        with self.assertRaisesRegex(runner.LOC2RunnerError, "independent_disagreement"):
            runner._require_route_agreement(first, key_changed, digest, "key")  # type: ignore[attr-defined]
        interval_changed = replace(
            second,
            candidates=(
                replace(second.candidates[0], parameter_lower=Fraction(2), parameter_upper=Fraction(2)),
                *second.candidates[1:],
            ),
        )
        with self.assertRaisesRegex(runner.LOC2RunnerError, "independent_disagreement"):
            runner._require_route_agreement(first, interval_changed, digest, "interval")  # type: ignore[attr-defined]

    def test_typed_root_inconclusive_is_invalid_and_not_published(self) -> None:
        sealed = _SealedFailure(classification="sufficient_contraction_failure")
        sealed_hash = runner.sha256(runner._canonical(runner._json_exact(sealed))).hexdigest()  # type: ignore[attr-defined]
        cubic = runner.primary.LocalCubic((0, 1, 0, 0), {})
        independent_cubic = runner.independent.IndependentLocalCubicV2((0, 1, 0, 0), {})
        with (
            patch.object(runner.loc1, "_row_hash", return_value="row"),
            patch.object(
                runner.loc1.pref1_exact,
                "assess_exact_temporal_refinement",
                return_value=sealed,
            ),
            patch.object(
                runner.loc1.pref1_exact_independent,
                "assess_exact_temporal_refinement_independently",
                return_value=sealed,
            ),
            patch.object(
                runner,
                "_cubics",
                return_value=([cubic, cubic], [independent_cubic, independent_cubic]),
            ),
            patch.object(
                runner.independent,
                "localize_absolute_maximum_independently_v2",
                side_effect=runner.independent.RootIsolationInconclusive(
                    "root_location_inconclusive", "synthetic"
                ),
            ),
            patch.object(runner, "_publish") as publish,
            patch.object(authority, "OWNED_ROW_COUNT", 1),
            self.assertRaisesRegex(
                runner.LOC2RunnerError, "root_location_inconclusive"
            ),
        ):
            runner._localize_occurrence(  # type: ignore[attr-defined]
                (),
                3,
                "u:alpha",
                {
                    "row_stream_sha256": "row",
                    "primary_evidence_sha256": sealed_hash,
                    "independent_evidence_sha256": sealed_hash,
                },
            )
        publish.assert_not_called()

    def test_prepare_error_is_wrapped_under_loc2_replay_owner(self) -> None:
        original = runner.loc1.LOC1RunnerError("replay", "synthetic", "detail")
        with patch.object(runner.loc1, "_prepare", side_effect=original):
            with self.assertRaisesRegex(
                runner.LOC2RunnerError, "sealed_LOC1_reconstruction_failed"
            ):
                runner._prepare_replay(ROOT, object(), {}, {})  # type: ignore[arg-type,attr-defined]

    def test_publication_rejects_symlink_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary, tempfile.TemporaryDirectory() as foreign:
            root = Path(temporary)
            (root / "runs").mkdir()
            (root / "runs" / "fgc-2-sf1").symlink_to(foreign)
            with self.assertRaises(OSError):
                runner._publish(root, {"kind": "manifest"}, {"kind": "terminal"})  # type: ignore[attr-defined]
            self.assertEqual(tuple(Path(foreign).iterdir()), ())

    def test_retry3_value_D01_hash_bound_family_matches(self) -> None:
        store = HLT16CampaignStore(ROOT / ar1.PREF2_STORE_PATH)
        shells = build_static_gr0_shells(ROOT)
        replay = next(item for item in ar1.REPLAYS if item["retry"] == 3)
        prepared, _historical = runner.loc1._prepare(ROOT, store, shells, replay)
        rows = tuple(runner.loc1._rows(runner.loc1._surface(prepared), "u:alpha"))
        self.assertEqual(
            runner.loc1._row_hash(rows),
            "6909fcf0d9b3b860ccb0a6cf9c535344f28c30c106db3667e08c9bcef77e389a",
        )
        first_cubics, second_cubics = runner._cubics(rows, "D01", "value_V")  # type: ignore[attr-defined]
        first = runner.primary.localize_absolute_maximum(
            first_cubics,
            maximum_candidates=4 * len(first_cubics),
            refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        )
        second = runner.independent.localize_absolute_maximum_independently_v2(
            second_cubics,
            maximum_candidates=4 * len(second_cubics),
        )
        digest = runner._primary_stationary_count_digest(first_cubics)  # type: ignore[attr-defined]
        runner._require_route_agreement(first, second, digest, "retry3-family")  # type: ignore[attr-defined]
        self.assertEqual((first.polynomial_count, first.candidate_count), (4088, 8176))
        self.assertEqual((second.polynomial_count, second.candidate_count), (4088, 8176))
        self.assertEqual(len(first.candidates), 2)
        self.assertEqual(len(second.candidates), 2)
        expected = Fraction(7421011090394622581, 10384593717069655257060992658440192)
        self.assertEqual(first.global_absolute_lower, expected)
        self.assertEqual(second.global_absolute_upper, expected)


if __name__ == "__main__":
    unittest.main()
