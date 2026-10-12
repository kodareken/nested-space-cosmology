"""Phase -1 public authority-boundary invariants."""

from __future__ import annotations

import inspect
from pathlib import Path
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]

from scripts.repo_checks import public  # noqa: E402
from scripts.repo_checks.public import (  # noqa: E402
    HISTORICAL_CHRONOLOGY_TOKENS,
    PHASE_MINUS1_REQUIRED_PHRASES,
    SLIMMED_PUBLIC_SURFACES,
    check_fgc_alignment,
)


class PhaseMinusOnePublicAuthorityTests(unittest.TestCase):
    def test_live_surfaces_satisfy_the_compact_authority_boundary(self) -> None:
        self.assertEqual(check_fgc_alignment(), [])

    def test_alignment_check_is_compact_and_not_a_chronology_matrix(self) -> None:
        source = inspect.getsource(check_fgc_alignment)
        self.assertLessEqual(len(source.splitlines()), 20)
        public_source = (ROOT / "scripts/repo_checks/public.py").read_text(
            encoding="utf-8"
        )
        start = public_source.index("SLIMMED_PUBLIC_SURFACES")
        end = public_source.index("def check_git_data_boundary()")
        self.assertLessEqual(public_source[start:end].count("\n"), 200)
        self.assertNotIn("missing PRO18 public-boundary phrase", public_source)
        self.assertNotIn("missing HLT16/PRO19 public-boundary phrase", public_source)
        self.assertNotIn("missing CFL1 current-boundary phrase", public_source)
        self.assertNotIn("missing framework contract phrase", public_source)

    def test_slimmed_surfaces_do_not_require_historical_chronology(self) -> None:
        self.assertEqual(
            SLIMMED_PUBLIC_SURFACES,
            {
                "README.md",
                "docs/finite-gradient-closure.md",
                "docs/epistemic-status.md",
                "docs/review-summary.md",
                "paper/recursive-horizons.md",
                "paper/fgc-local-defocusing/README.md",
            },
        )
        for relative in SLIMMED_PUBLIC_SURFACES:
            joined = "\n".join(PHASE_MINUS1_REQUIRED_PHRASES[relative])
            for token in HISTORICAL_CHRONOLOGY_TOKENS:
                with self.subTest(relative=relative, token=token):
                    self.assertNotIn(token, joined)

    def test_claim_ledger_owns_qa_hashes_and_nonpromotion(self) -> None:
        owned = "\n".join(PHASE_MINUS1_REQUIRED_PHRASES["docs/claim-ledger.md"])
        self.assertIn("| C231 | **FGC-1-TDG10-QA1-PREF1.**", owned)
        self.assertIn("| C233 | **FGC-1-TDG10-QA2-PREF1.**", owned)
        self.assertIn(
            "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6",
            owned,
        )
        self.assertIn("successor_remedy_selected=false", owned)
        self.assertIn("state_advance_authorized=false", owned)
        self.assertIn("nonpromotion controls", owned)

    def test_roadmap_owns_frontier_stop_and_finish_line(self) -> None:
        owned = "\n".join(PHASE_MINUS1_REQUIRED_PHRASES["docs/research-roadmap.md"])
        self.assertIn("This is **FGC-2-SF1**", owned)
        self.assertIn("No production temporal successor has been selected", owned)
        self.assertIn("FGC-1-TDG11-MSEL1-FRZ1", owned)
        self.assertIn("No push, publication", owned)

    def test_readme_owns_current_qa2_routing_without_state_advance(self) -> None:
        owned = "\n".join(PHASE_MINUS1_REQUIRED_PHRASES["README.md"])
        self.assertIn("FGC-1-TDG10-QA2-PREF1", owned)
        self.assertIn("No production temporal successor is selected", owned)
        self.assertIn("or physical claim is", owned)
        self.assertIn("mutate a terminal store, or run without a new committed", owned)

    def test_plan_owns_scientific_authority_order(self) -> None:
        owned = PHASE_MINUS1_REQUIRED_PHRASES["PLAN.md"]
        self.assertEqual(
            owned[0],
            "1. `AGENTS.md` — creative research orientation.",
        )
        self.assertEqual(
            owned[-1],
            "7. `README.md` — short project entrypoint and routing page only.",
        )
        self.assertTrue(
            all(row.startswith(f"{index}. ") for index, row in enumerate(owned, start=1))
        )

    def test_agents_hypothesis_boundary_does_not_pin_user_metadata(self) -> None:
        self.assertFalse(any(
            phrase.startswith("status:")
            for phrase in PHASE_MINUS1_REQUIRED_PHRASES["AGENTS.md"]
        ))
        for status in (
            'status: "Working hypothesis and exploration map"',
            'status: "Work in progress hypothesis. No new concepts claimed."',
        ):
            with self.subTest(status=status):
                self.assertEqual(self._failures_after(
                    "AGENTS.md",
                    lambda text: "\n".join(
                        status if line.startswith("status:") else line
                        for line in text.splitlines()
                    ),
                ), [])

    def test_missing_agents_hypothesis_or_nonadvocacy_boundary_fails(self) -> None:
        for phrase in PHASE_MINUS1_REQUIRED_PHRASES["AGENTS.md"]:
            with self.subTest(phrase=phrase):
                failures = self._failures_after("AGENTS.md", lambda text: text.replace(phrase, ""))
                self.assertTrue(any("AGENTS.md" in item and phrase in item for item in failures))

    def test_missing_qa2_hash_fails_before_other_surfaces(self) -> None:
        failures = self._failures_after(
            "docs/claim-ledger.md",
            lambda text: text.replace(
                "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6",
                "0" * 64,
            ),
        )
        self.assertTrue(failures)
        self.assertTrue(all("docs/claim-ledger.md" in item for item in failures))
        self.assertTrue(
            any(
                "08f88f02966ea5936e417a3a6a2c583e11783341a3d76f4ef5dc821e992691f6"
                in item
                for item in failures
            )
        )

    def test_missing_ledger_nonpromotion_fails(self) -> None:
        failures = self._failures_after(
            "docs/claim-ledger.md",
            lambda text: text.replace("successor_remedy_selected=false", "successor_remedy_selected=true"),
        )
        self.assertTrue(
            any("successor_remedy_selected=false" in item for item in failures)
        )

    def test_missing_roadmap_finish_line_fails(self) -> None:
        failures = self._failures_after(
            "docs/research-roadmap.md",
            lambda text: text.replace("This is **FGC-2-SF1**", "This is a later local test"),
        )
        self.assertTrue(any("FGC-2-SF1" in item for item in failures))

    def test_missing_readme_qa2_routing_fails(self) -> None:
        failures = self._failures_after(
            "README.md",
            lambda text: text.replace("FGC-1-TDG10-QA2-PREF1", "FGC-1-TDG10-QA1-PREF1"),
        )
        self.assertTrue(any("README.md" in item and "QA2-PREF1" in item for item in failures))

    def test_missing_fgc_law_fails(self) -> None:
        failures = self._failures_after(
            "docs/finite-gradient-closure.md",
            lambda text: text.replace("Finite Gradient Closure", "Finite Placeholder Closure"),
        )
        self.assertTrue(
            any("docs/finite-gradient-closure.md" in item for item in failures)
        )

    def test_missing_ngs_nonclaim_fails(self) -> None:
        failures = self._failures_after(
            "docs/nested-gradient-spectrum.md",
            lambda text: text.replace(
                "The golden ratio is a possible eigenvalue, not a premise",
                "The golden ratio is a required premise",
            ),
        )
        self.assertTrue(
            any("possible eigenvalue, not a premise" in item for item in failures)
        )

    def test_missing_epistemic_label_fails(self) -> None:
        failures = self._failures_after(
            "docs/epistemic-status.md",
            lambda text: text.replace("### Q — Open question", "### Q — Quiet note"),
        )
        self.assertTrue(any("Open question" in item for item in failures))

    def test_missing_review_archive_stub_fails(self) -> None:
        failures = self._failures_after(
            "docs/review-summary.md",
            lambda text: text.replace(
                "archive/reviews/review-summary-2026-08-21.md",
                "docs/review-summary.md",
            ),
        )
        self.assertTrue(
            any("review-summary-2026-08-21.md" in item for item in failures)
        )

    def test_missing_paper_live_routing_fails(self) -> None:
        failures = self._failures_after(
            "paper/recursive-horizons.md",
            lambda text: text.replace("**Live development routing.**", "**Status log.**"),
        )
        self.assertTrue(any("Live development routing" in item for item in failures))

    def test_missing_workspace_no_successor_fails(self) -> None:
        failures = self._failures_after(
            "paper/fgc-local-defocusing/README.md",
            lambda text: text.replace("It selects no successor", "It selects a successor"),
        )
        self.assertTrue(any("selects no successor" in item for item in failures))

    def test_missing_authority_order_fails(self) -> None:
        failures = self._failures_after(
            "PLAN.md",
            lambda text: text.replace(
                "1. `AGENTS.md` — creative research orientation.",
                "1. `configs/fgc/artifact-catalog.json` — scientific authority.",
            ),
        )
        self.assertTrue(any("AGENTS.md" in item for item in failures))

    def test_absent_surface_is_reported(self) -> None:
        def reader(relative: str) -> str:
            if relative == "docs/review-summary.md":
                raise FileNotFoundError("review stub missing")
            return (ROOT / relative).read_text(encoding="utf-8")

        with patch.object(public, "_read_public_surface", side_effect=reader):
            failures = check_fgc_alignment()
        self.assertEqual(len(failures), 1)
        self.assertIn("cannot read framework surface docs/review-summary.md", failures[0])

    def test_removing_historical_ids_from_readme_still_passes(self) -> None:
        def transform(text: str) -> str:
            mutated = text
            for token in HISTORICAL_CHRONOLOGY_TOKENS:
                mutated = mutated.replace(token, "")
            return mutated

        self.assertEqual(self._failures_after("README.md", transform), [])

    def _failures_after(self, relative: str, transform) -> list[str]:
        original = (ROOT / relative).read_text(encoding="utf-8")
        mutated = transform(original)

        def reader(name: str) -> str:
            if name == relative:
                return mutated
            return (ROOT / name).read_text(encoding="utf-8")

        with patch.object(public, "_read_public_surface", side_effect=reader):
            return check_fgc_alignment()


if __name__ == "__main__":
    unittest.main()
