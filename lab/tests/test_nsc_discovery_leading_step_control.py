"""Full local scaling, immutable handoffs, spawn propagation and bounded local replay."""
import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_leading_step_control as step


def endpoint(case="nf128_coherent_dt0.0005"):
    record, arrays = step.leading.episode.load_checkpoint(step.SOURCE, case)
    pair = step.leading.backend.make_fft_pair(step.leading.episode.pair_from_arrays(arrays, record))
    return pair, step.leading.state_from_arrays(arrays), record


def test_full28_row_norm_includes_every_field_link_and_moving_scale_triangle():
    with threadpool_limits(limits=1), step.leading.backend.fft_thread_limit(1):
        pair, state, record = endpoint()
        rate, bundle = step.leading.rates(pair, state, return_bundle=True)
        fine, system = bundle["fine_state"], bundle["fine_system"]
        rows, scales, _frequency, _kg, _kf = step.scaled_local_reaction(pair.grid, fine, system)
        # Independent unscaled Jacobian at one actual point, transformed by
        # the declared coordinate scales (including 24 real field slots).
        point = pair.grid.nq//4
        one = step.leading.State(*(getattr(fine, name)[point:point+1].copy() for name in step.leading.FIELDS))
        one_system = step.local_system(system)
        one_system.length_density = one.Q
        one_system.shift = np.zeros(1)
        unscaled = np.zeros((28, 28))
        for column in range(28):
            t = step.leading.State(*(np.zeros_like(getattr(one, name)) for name in step.leading.FIELDS))
            if column < 4:
                setattr(t, step.leading.REAL_FIELDS[column], np.ones(1))
            else:
                index = column-4
                field = t.phi0 if index<12 else t.phi1
                index %= 12
                field[0, index%6] = 1 if index<6 else 1j
            image = step.leading.fine_jvp(one_system, one, t)
            unscaled[:, column] = np.concatenate(
                [np.array([getattr(image, name)[0] for name in step.leading.REAL_FIELDS]),
                 image.phi0.real[0], image.phi0.imag[0], image.phi1.real[0], image.phi1.imag[0]])
        physical = np.r_[scales[point], np.full(24, np.sqrt(pair.grid.dx_q))]
        weighted = unscaled*physical[None, :]/physical[:, None]
        np.testing.assert_allclose(rows[point], np.sum(abs(weighted), axis=1), rtol=1e-12, atol=1e-10)
        dt, metadata = step.step_restriction(pair, state, record["step_cap"])
        assert metadata["local_real_variables"] == 28
        assert metadata["absolute_moving_scale_rate"] > 0
        assert metadata["reaction_majorant"] == pytest.approx(
            metadata["frozen_scaled_reaction_norm"]+metadata["absolute_moving_scale_rate"])
        assert not metadata["physical_growth_removed"] and not metadata["stability_certificate"]
        assert 0 < dt <= record["step_cap"]


def test_exact_payload_and_prefix_copy_completed_arms_remain_terminal(tmp_path):
    with threadpool_limits(limits=1):
        for case in ("nf128_coherent_dt0.0005", "nf128_coherent_frozen_geometry_dt0.0005"):
            old, arrays = step.leading.episode.load_checkpoint(step.SOURCE, case)
            info = step.copy_case(step.SOURCE, tmp_path, case, producer_commit="unit-test")
            copied, values = step.leading.episode.load_checkpoint(tmp_path, case, 0)
            assert copied["array_sha256"] == old["array_sha256"]
            assert copied["arrays_sha256"] == old["arrays_sha256"]
            for name in arrays:
                np.testing.assert_array_equal(values[name], arrays[name])
            original_prefix = (step.SOURCE/(case+"-observations.jsonl")).read_bytes()
            assert (tmp_path/(case+"-observations.jsonl")).read_bytes() == original_prefix
            assert info["prefix_sha256"] == hashlib.sha256(original_prefix).hexdigest()
            if "frozen" in case:
                assert copied["status"] == "station_reached" and copied["already_complete"]
                with step.numerical_adapter(), step.leading.episode_adapter():
                    result = step.leading.episode._case_worker({
                        "directory": str(tmp_path), "case_id": case, "cpu_budget_seconds": 165.,
                        "forecast_factor": 1.5})
                assert result["child_cpu_seconds"] == 0 and result["chunk"] is None


def test_paired_new_steps_preserve_rhs_source_hashes_and_small_local_error():
    report = step.paired_step_check(case_ids=("nf128_coherent_dt0.001", "nf128_coherent_dt0.0005"))
    assert report["frozen_source_hashes_unchanged"] and report["bounded_local_steps_only"]
    for row in report["checks"]:
        assert row["RHS_bit_identical"] and row["new_dt"] > 10*row["old_dt"]
        assert row["one_step_two_half_geometric_scaled_max_gap"] < 1e-6
        assert row["one_step_two_half_weighted_field_relative_gap"] < 1e-7


def test_schema_adapter_restores_and_legacy_cli_rejects_successor(tmp_path):
    before = step.leading.SCHEMA, step.leading.step_restriction, step.leading.pool, step.leading.rates
    with step.numerical_adapter():
        assert step.leading.SCHEMA == step.SCHEMA
        assert step.leading.step_restriction is step.step_restriction
        assert step.leading.pool is step.pool
        assert step.leading.rates is before[3]
    assert (step.leading.SCHEMA, step.leading.step_restriction, step.leading.pool, step.leading.rates) == before
    (tmp_path/"manifest.json").write_text(json.dumps({"schema": step.SCHEMA}))
    with pytest.raises(ValueError, match="not a leading"):
        step.leading.run(tmp_path)
    with pytest.raises(ValueError, match="165"):
        step.prepare(cpu_budget=166.)


def test_spawn_safe_initializer_installs_new_schema_and_step_without_a_trajectory():
    with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn"),
                             initializer=step.worker_initializer) as executor:
        result = executor.submit(step.worker_identity).result(timeout=30)
    assert result["leading_schema"] == result["episode_schema"] == step.SCHEMA
    assert result["scaled_step_installed"]
    assert result["original_RHS_module"].endswith("nsc_discovery_leading_einstein")
    assert not result["trajectory_executed"]
