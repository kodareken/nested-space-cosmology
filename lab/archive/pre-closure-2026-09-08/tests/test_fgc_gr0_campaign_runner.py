from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration import (  # noqa: E402
    NormalFlowTracers,
    _campaign_checkpoint_metadata,
    _event_log_identity,
    _recover_event_log_from_checkpoint,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
)


def _flat_state() -> tuple[EvolutionState, np.ndarray]:
    coordinates = np.linspace(0.0, 4.0, 17)
    u = np.zeros((coordinates.size, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = coordinates
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q), coordinates


class FGCGR0CampaignRunnerTests(unittest.TestCase):
    def test_tracer_preview_is_atomic_until_commit(self) -> None:
        state, coordinates = _flat_state()
        tracers = NormalFlowTracers.create(
            minimum=0.5,
            maximum=3.5,
            spacing=0.5,
            state=state,
            coordinates=coordinates,
            cutoff=16.0,
            outer_radius=4.0,
        )
        positions_before = tracers.positions.copy()
        proper_before = tracers.proper_times.copy()
        positions, proper = tracers.preview_advance(
            old_state=state,
            new_state=state,
            coordinates=coordinates,
            step_size=0.125,
        )
        np.testing.assert_array_equal(tracers.positions, positions_before)
        np.testing.assert_array_equal(tracers.proper_times, proper_before)
        np.testing.assert_array_equal(positions, positions_before)
        np.testing.assert_allclose(proper, 0.125)
        tracers.commit_advance(positions, proper)
        np.testing.assert_array_equal(tracers.positions, positions)
        np.testing.assert_array_equal(tracers.proper_times, proper)

    def test_resume_archives_only_a_complete_uncommitted_tail(self) -> None:
        checkpoint_log = b'{"event_index":0}\n'
        uncommitted = b'{"event_index":1}\n'
        with tempfile.TemporaryDirectory() as directory:
            event_log = Path(directory) / "events.jsonl"
            event_log.write_bytes(checkpoint_log + uncommitted)
            recovery = _recover_event_log_from_checkpoint(
                event_log,
                checkpoint_log,
            )
            self.assertEqual(event_log.read_bytes(), checkpoint_log)
            self.assertEqual(
                recovery["action"],
                "archived_uncommitted_tail_and_restored_checkpoint",
            )
            orphan = event_log.with_name(recovery["orphaned_tail_path"])
            self.assertEqual(orphan.read_bytes(), uncommitted)
            self.assertEqual(
                recovery["orphaned_tail"],
                _event_log_identity(uncommitted),
            )

    def test_resume_rejects_divergent_event_history(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            event_log = Path(directory) / "events.jsonl"
            event_log.write_bytes(b'{"event_index":99}\n')
            with self.assertRaisesRegex(ValueError, "diverges"):
                _recover_event_log_from_checkpoint(
                    event_log,
                    b'{"event_index":0}\n',
                )

    def test_checkpoint_metadata_carries_prior_amplitude_and_log_identity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            event_log = Path(directory) / "events.jsonl"
            payload = b'{"amplitude":"5/2","event_index":0}\n'
            event_log.write_bytes(payload)
            prior = [{"amplitude": "5/2", "eligible": False}]
            metadata = _campaign_checkpoint_metadata(
                manifest={"campaign_id": "campaign"},
                amplitude="3",
                amplitude_index=1,
                event_index=0,
                consecutive_qualified_events=0,
                amplitude_records=prior,
                event_log_path=event_log,
                members={},
                terminal=False,
            )
            self.assertEqual(metadata["completed_amplitude_records"], prior)
            self.assertEqual(metadata["event_log"], _event_log_identity(payload))


if __name__ == "__main__":
    unittest.main()
