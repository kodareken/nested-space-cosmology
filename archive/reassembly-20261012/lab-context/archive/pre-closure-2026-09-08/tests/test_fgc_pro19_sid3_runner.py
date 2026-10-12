from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace
from unittest.mock import patch
import unittest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)

from recursive_horizons.fgc.evolution.proto17_pure_construction import (  # noqa: E402
    MEMBER_KEYS,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    STATIC_INPUT_PATHS,
)


def _load_runner():
    spec = importlib.util.spec_from_file_location(
        "fgc_pro19_sid3_runner_test_module",
        ROOT / "scripts" / "run_fgc_pro19_event1.py",
    )
    if spec is None or spec.loader is None:
        raise AssertionError("event-one runner cannot be imported")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


runner = _load_runner()


def _member(time: float, *, mode: str = "FRESH_READY") -> SimpleNamespace:
    return SimpleNamespace(
        cursor={
            "accepted_boundary_time": {
                "rational": str(time),
                "binary64_hex": time.hex(),
            },
            "mode": mode,
        },
        pending_owner=None,
        pending_cap_hex=None,
        cfl_current=0,
        cfl_total=0,
        ledger={"cumulative_temporal_retry_count": 0},
    )


def _checkpoint() -> SimpleNamespace:
    return SimpleNamespace(
        generation=9,
        sha256="9" * 64,
        journal_sequence=9,
        journal_tip_sha256="8" * 64,
        authorization_commit="1" * 40,
        plan_sha256="2" * 64,
        campaign_id="campaign",
        disposition="nonterminal",
        terminal=None,
        event=23,
        target={"rational": "3/2", "binary64_hex": (3 / 2).hex()},
        members={key: _member(23 / 16) for key in MEMBER_KEYS},
    )


def _authority() -> SimpleNamespace:
    return SimpleNamespace(
        store_path="runs/fgc-2-sf1/proto17/calibration",
        progression_plan=object(),
        static_input_bytes={
            **{path: path.encode("ascii") for path in STATIC_INPUT_PATHS},
            "extra": b"bound compact evidence",
        },
        implementation_commit="3" * 40,
        authority_commit="4" * 40,
        closure_sha256="5" * 64,
        config_sha256="6" * 64,
        result_sha256="7" * 64,
        checkpoint_generation=8,
        checkpoint_sha256="a" * 64,
    )


class Proto19SID3RunnerTests(unittest.TestCase):
    def test_template_gate_rejects_non_gr0_operator_or_amplitude(self) -> None:
        class Operator:
            pass

        class Member:
            def __init__(self, key, *, amplitude="3", operator=None):
                self.key = key
                self.amplitude = amplitude
                self.operator = Operator() if operator is None else operator

        templates = {key: Member(key) for key in MEMBER_KEYS}
        with (
            patch.object(runner, "Proto14RunMember", Member),
            patch.object(runner, "Proto12GR0EvolutionOperator", Operator),
        ):
            runner._require_sid3_gr0_templates(templates)
            wrong = dict(templates)
            wrong[MEMBER_KEYS[0]] = Member(MEMBER_KEYS[0], amplitude="4")
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError, "template identity"
            ):
                runner._require_sid3_gr0_templates(wrong)
            wrong = dict(templates)
            wrong[MEMBER_KEYS[0]] = Member(MEMBER_KEYS[0], operator=object())
            with self.assertRaisesRegex(
                runner.Proto19Event1RunnerError, "template identity"
            ):
                runner._require_sid3_gr0_templates(wrong)

    def test_importing_persisted_runtime_and_runner_does_not_open_gen0_factory(
        self,
    ) -> None:
        script = f"""
import importlib, importlib.util, pathlib, sys, types
root = pathlib.Path({str(ROOT)!r})
sys.path[:0] = [str(root / 'src')]
for name, path in (
    ('recursive_horizons', root / 'src' / 'recursive_horizons'),
    ('recursive_horizons.fgc', root / 'src' / 'recursive_horizons' / 'fgc'),
    ('recursive_horizons.fgc.evolution', root / 'src' / 'recursive_horizons' / 'fgc' / 'evolution'),
    ('scripts', root / 'scripts'),
):
    package = types.ModuleType(name)
    package.__package__ = name
    package.__path__ = [str(path)]
    sys.modules[name] = package
importlib.import_module('recursive_horizons.fgc.evolution.hlt16_campaign_runtime')
spec = importlib.util.spec_from_file_location('sid3_import_probe', root / 'scripts' / 'run_fgc_pro19_event1.py')
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
for forbidden in (
    'recursive_horizons.fgc.evolution.proto19_progression_inputs',
    'recursive_horizons.fgc.evolution.static_initial_admission',
    'scripts.run_fgc_gr0_calibration_v13',
):
    assert forbidden not in sys.modules, forbidden
"""
        completed = subprocess.run(
            [sys.executable, "-c", script],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_read_only_preflight_restores_all_six_without_writer_or_advance(
        self,
    ) -> None:
        authority = _authority()
        checkpoint = _checkpoint()
        snapshot = SimpleNamespace(
            checkpoint=checkpoint,
            terminal_lock_present=False,
            suffix_classification="clean",
            suffix_records=(),
            staging_paths=(),
            orphan_descriptor_sha256=None,
            orphan_payload_semantic_sha256=None,
        )
        status = SimpleNamespace(
            terminal=False,
            active_write=False,
            writer_state=None,
            state="clean_checkpoint",
        )
        store = object()
        templates = {key: SimpleNamespace(key=key) for key in MEMBER_KEYS}

        def restore(_store, _checkpoint, template, *, key):
            self.assertIs(_store, store)
            self.assertIs(_checkpoint, checkpoint)
            self.assertEqual(template.key, key)
            template.time = 23 / 16
            return template

        with (
            patch.object(runner, "_validated_sid3_authority", return_value=authority),
            patch.object(
                runner, "HLT16CampaignStore", return_value=store
            ) as store_type,
            patch.object(
                runner,
                "_require_sid3_generation8_ancestor",
                return_value=(checkpoint, snapshot, status),
            ),
            patch.object(
                runner, "build_static_gr0_shells", return_value=templates
            ) as factory,
            patch.object(runner, "_require_sid3_gr0_templates") as template_guard,
            patch.object(
                runner.runtime, "restore_member_with_overlay", side_effect=restore
            ) as restore_member,
            patch.object(runner, "_member_to_advance", return_value=MEMBER_KEYS[0]),
            patch.object(runner, "_acquire_writer") as acquire,
            patch.object(runner, "_advance_authenticated_event") as advance,
        ):
            result = runner.sid3_resume_preflight(
                ROOT,
                authority=authority,
                store_root=ROOT / authority.store_path,
            )

        store_type.assert_called_once_with((ROOT / authority.store_path).resolve())
        factory.assert_called_once_with(
            static_input_bytes={
                path: authority.static_input_bytes[path] for path in STATIC_INPUT_PATHS
            }
        )
        template_guard.assert_called_once_with(templates)
        self.assertEqual(restore_member.call_count, len(MEMBER_KEYS))
        acquire.assert_not_called()
        advance.assert_not_called()
        self.assertTrue(result["six_member_restore_passed"])
        self.assertTrue(result["safe_to_resume_trajectory"])
        self.assertFalse(result["writer_lease_acquired"])
        self.assertFalse(result["state_advanced"])

    def test_resume_consumes_the_same_preflight_templates_and_original_plan(
        self,
    ) -> None:
        authority = _authority()
        checkpoint = _checkpoint()
        templates = {key: object() for key in MEMBER_KEYS}
        context = runner.SID3ResumePreflightContext(
            authority=authority,
            store=object(),
            checkpoint=checkpoint,
            templates=templates,
            suffix_classification="clean",
            recovery_state="clean_checkpoint",
        )

        def advance(plan, store, observed, supplied):
            self.assertIs(plan, authority.progression_plan)
            self.assertIs(store, context.store)
            self.assertIs(observed, checkpoint)
            self.assertIs(supplied, templates)
            return {"schema": "event", "state": "event_one_complete"}

        with (
            patch.object(
                runner, "_sid3_resume_preflight_context", return_value=context
            ) as preflight,
            patch.object(
                runner, "_advance_authenticated_event", side_effect=advance
            ) as execute,
            patch.object(runner, "_require_sid3_gr0_templates") as template_guard,
            patch.object(runner, "build_static_gr0_shells") as rebuild,
        ):
            result = runner.sid3_resume(
                ROOT,
                authority=authority,
                store_root=ROOT / authority.store_path,
            )
        preflight.assert_called_once_with(
            ROOT,
            authority=authority,
            store_root=ROOT / authority.store_path,
            permit_recoverable_suffix=True,
        )
        execute.assert_called_once()
        template_guard.assert_called_once_with(templates)
        rebuild.assert_not_called()
        self.assertEqual(result["state"], "event_one_complete")
        self.assertFalse(result["candidate_branch_opened"])
        self.assertFalse(result["calibration_result_earned"])
        self.assertFalse(result["physical_result_earned"])


if __name__ == "__main__":
    unittest.main()
