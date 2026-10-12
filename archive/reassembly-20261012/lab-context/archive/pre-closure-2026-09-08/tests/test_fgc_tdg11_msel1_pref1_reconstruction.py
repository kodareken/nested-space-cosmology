"""Outcome-blind independent recorded-RK controls; no campaign inputs."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import unittest

import numpy as np

from recursive_horizons.fgc.evolution.numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    propose_step,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg11_compensated_rk import (
    ARITHMETIC_ID as COMPENSATED_ARITHMETIC_ID,
    propose_compensated_step,
)
from recursive_horizons.fgc.evolution.tdg11_msel1_pref1_reconstruction import (
    ORIGINAL_ARITHMETIC_ID,
    IndependentChannelReconstruction,
    IndependentRecordedFamily,
    IndependentRecordedFamilyError,
    reconstruct_channel_independently,
    validate_recorded_family_independently,
)


Q = Fraction
ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT / "src/recursive_horizons/fgc/evolution/tdg11_msel1_pref1_reconstruction.py"
)
OWNED = slice(1, -4)
BLOCKS = ("u", "p", "q")
DERIVATIVES = ("du", "dp", "dq")
FIELDS = ("alpha", "v", "lambda", "R", "phi", "chi")
POINT_COUNT = 9
FIELD_COUNT = 6
OWNED_ROWS = POINT_COUNT - 5
LARGE_BASE = 2.0**52
MACRO_WIDTH = 0.25
HASH_DOMAIN = b"TDG11-MSEL1-RECORDED-FAMILY-v1\n"
FORBIDDEN_IMPORT_MARKERS = (
    "tdg11_msel1_reconstruction",
    "tdg11_rational_complete_c",
    "tdg11_msel1_runtime",
    "tdg11_msel1_authority",
    "tdg11_msel1_contract",
    "tdg5_",
    "tdg6_temporal_admission_runtime",
    "tdg6_temporal_admission_theorem",
    "tdg7_",
    "tdg8_",
    "tdg9_ac1",
    "tdg9_ar1",
    "tdg9_loc1",
    "tdg9_loc2",
    "tdg9_ti1",
    "tdg9_ti2",
    "tdg9_ur1",
    "tdg10_",
    "hlt16",
    "binder",
    "scripts",
    "proto",
    "campaign",
    "store",
    "propose_step",
    "propose_compensated_step",
    "validate_recorded_family",
    "reconstruct_channel",
    "assess_rational_complete_c",
)
_PROPOSERS = (
    (propose_step, ORIGINAL_ARITHMETIC_ID),
    (propose_compensated_step, COMPENSATED_ARITHMETIC_ID),
)
_RK4_WEIGHTS = (Q(1, 6), Q(1, 3), Q(1, 3), Q(1, 6))
_SSP_WEIGHTS = (Q(1, 6), Q(1, 6), Q(2, 3))


def _imported_names(path: Path) -> set[str]:
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                names.add(node.module)
            names.update(alias.name for alias in node.names)
    return names


def _state(*, u=1.0, p=0.0, q=0.0, points=POINT_COUNT, fields=FIELD_COUNT):
    return EvolutionState(
        np.full((points, fields), u, dtype=np.float64),
        np.full((points, fields), p, dtype=np.float64),
        np.full((points, fields), q, dtype=np.float64),
    )


def _rhs_from_functions(du, dp, dq):
    def rhs(_time, state):
        return EvolutionRHS(du(state), dp(state), dq(state), {})

    return rhs


def _constant_rhs(du=0.0, dp=0.0, dq=0.0):
    return _rhs_from_functions(
        lambda state: np.full_like(state.u, du),
        lambda state: np.full_like(state.p, dp),
        lambda state: np.full_like(state.q, dq),
    )


def _unit_u_rhs():
    return _constant_rhs(du=1.0)


def _linear_rhs():
    return _rhs_from_functions(
        lambda state: np.array(state.u, copy=True),
        lambda state: np.array(state.p, copy=True),
        lambda state: np.array(state.q, copy=True),
    )


def _zero_rhs():
    return _constant_rhs()


def _fingerprint_state(state):
    return tuple(getattr(state, name).tobytes() for name in BLOCKS)


def _fingerprint_rhs(rhs):
    return (rhs.du.tobytes(), rhs.dp.tobytes(), rhs.dq.tobytes())


def _fingerprint_paths(paths):
    parts = []
    for path in paths:
        for proposal in path:
            parts.append(_fingerprint_state(proposal.initial_state))
            parts.append(_fingerprint_state(proposal.candidate_state))
            for stage in proposal.stages:
                parts.append(_fingerprint_state(stage.state))
                parts.append(_fingerprint_rhs(stage.rhs))
    return tuple(parts)


def _path(proposer, method, state, rhs, start, width, steps):
    proposals = []
    current_state, current_time, step = state, start, width / steps
    for _ in range(steps):
        proposal = proposer(
            method=method,
            time=current_time,
            step_size=step,
            state=current_state,
            rhs=rhs,
        )
        proposals.append(proposal)
        current_state = proposal.candidate_state
        current_time = proposal.final_time
    return tuple(proposals)


def _paths(proposer, method, state, rhs, *, width=MACRO_WIDTH, start=0.0):
    return (
        _path(proposer, method, state, rhs, start, width, 1),
        _path(proposer, method, state, rhs, start, width, 2),
        _path(proposer, method, state, rhs, start, width, 4),
    )


def _validate(paths, method, arithmetic_id, rows=OWNED_ROWS):
    return validate_recorded_family_independently(
        paths,
        method=method,
        arithmetic_id=arithmetic_id,
        owned_row_count=rows,
    )


def _family(proposer, arithmetic_id, method, state, rhs, **kwargs):
    paths = _paths(proposer, method, state, rhs, **kwargs)
    return paths, _validate(paths, method, arithmetic_id)


def _replace_stage(proposal, index, **changes):
    stages = list(proposal.stages)
    stages[index] = replace(stages[index], **changes)
    return replace(proposal, stages=tuple(stages))


def _mutate_owned_state(state, *, block="u", amount=1.0):
    data = {name: np.array(getattr(state, name), copy=True) for name in BLOCKS}
    data[block][1, 0] += amount
    return EvolutionState(data["u"], data["p"], data["q"])


def _mutate_owned_rhs(rhs, *, amount=1.0):
    du = np.array(rhs.du, copy=True)
    du[1, 0] += amount
    return EvolutionRHS(du, rhs.dp, rhs.dq, dict(rhs.diagnostics))


def _path_deltas(raw_segments, corrected_segments):
    deltas = []
    for raw, corrected in zip(raw_segments, corrected_segments, strict=True):
        deltas.append((raw[2] - raw[0]) - (corrected[2] - corrected[0]))
    return tuple(deltas)


def _feed_array(hasher, array):
    hasher.update(np.asarray(array.shape, dtype=np.int64).tobytes())
    hasher.update(np.ascontiguousarray(array, dtype=np.float64).tobytes())


def _bits(value):
    if type(value) not in (float, np.float64):
        raise TypeError("times must be binary64")
    return np.float64(float(value)).tobytes()


def _independent_family_digest(paths, method, arithmetic_id, owned_row_count):
    hasher = sha256(HASH_DOMAIN)
    hasher.update(f"{method}\n{arithmetic_id}\n{owned_row_count}\n".encode("ascii"))
    for path in paths:
        for proposal in path:
            hasher.update(proposal.method.encode("ascii") + b"\n")
            hasher.update(_bits(proposal.initial_time))
            hasher.update(_bits(proposal.final_time))
            for state in (proposal.initial_state, proposal.candidate_state):
                for block in BLOCKS:
                    _feed_array(hasher, getattr(state, block))
            for stage in proposal.stages:
                hasher.update(stage.stage_name.encode("ascii") + b"\n")
                hasher.update(_bits(stage.time))
                for block in BLOCKS:
                    _feed_array(hasher, getattr(stage.state, block))
                for derivative in DERIVATIVES:
                    _feed_array(hasher, getattr(stage.rhs, derivative))
    return hasher.hexdigest()


class ImportBoundaryTests(unittest.TestCase):
    def test_module_does_not_import_production_reconstruction_or_runner(self):
        imported = _imported_names(MODULE_PATH)
        source = MODULE_PATH.read_text()
        self.assertIn("numerical_engine", imported)
        self.assertIn("tdg11_compensated_rk", imported)
        self.assertIn("compensated_update", imported)
        self.assertIn("CompensatedArithmeticStop", imported)
        self.assertNotIn("tdg6_temporal_admission_design", imported)
        self.assertNotIn("TDG6_COMPLETE_STATE_CHANNELS", imported)
        self.assertIn("TDG11-MSEL1-RECORDED-FAMILY-v1", source)
        self.assertIn(
            "never emits a scientific pass label",
            " ".join(source.split()),
        )
        self.assertIn("slice(1, -4)", source)
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertFalse(
                any(marker in name for name in imported),
                msg=f"forbidden import marker {marker!r} in {imported}",
            )
        self.assertNotIn("from .tdg11_msel1_reconstruction", source)
        self.assertNotIn("import tdg11_msel1_reconstruction", source)
        self.assertNotIn("from .tdg11_rational_complete_c", source)
        self.assertNotIn("from .tdg11_msel1_runtime", source)
        self.assertNotIn("propose_step", imported)
        self.assertNotIn("propose_compensated_step", imported)
        self.assertNotIn("import scripts", source)
        self.assertNotIn("from scripts", source)
        test_imported = _imported_names(Path(__file__))
        self.assertNotIn(
            "recursive_horizons.fgc.evolution.tdg11_msel1_reconstruction",
            test_imported,
        )
        self.assertNotIn(
            "recursive_horizons.fgc.evolution.tdg11_rational_complete_c",
            test_imported,
        )
        self.assertNotIn(
            "recursive_horizons.fgc.evolution.tdg11_msel1_runtime",
            test_imported,
        )
        self.assertFalse(
            any("binder" in name or "scripts" in name for name in test_imported)
        )


class RecordedFamilyTests(unittest.TestCase):
    def test_both_proposers_methods_and_all_channels_validate(self):
        self.assertEqual(ORIGINAL_ARITHMETIC_ID, "tdg11_recorded_legacy_binary64_rk_v1")
        self.assertEqual(
            COMPENSATED_ARITHMETIC_ID,
            "tdg11_binary64_twofold_dot2_coherent_rk_v1",
        )
        self.assertEqual(len(TDG6_COMPLETE_STATE_CHANNELS), 18)
        state = _state()
        rhs = _unit_u_rhs()
        for proposer, arithmetic_id in _PROPOSERS:
            for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
                paths, family = _family(proposer, arithmetic_id, method, state, rhs)
                self.assertIsInstance(family, IndependentRecordedFamily)
                self.assertEqual(tuple(len(path) for path in family.paths), (1, 2, 4))
                self.assertEqual(family.method, method)
                self.assertEqual(family.arithmetic_id, arithmetic_id)
                self.assertEqual(family.owned_row_count, OWNED_ROWS)
                self.assertEqual(len(family.family_sha256), 64)
                self.assertEqual(
                    family.family_sha256,
                    _independent_family_digest(
                        paths, method, arithmetic_id, OWNED_ROWS
                    ),
                )
                self.assertEqual(
                    family.family_sha256,
                    _validate(paths, method, arithmetic_id).family_sha256,
                )
                for channel in TDG6_COMPLETE_STATE_CHANNELS:
                    rebuilt = reconstruct_channel_independently(family, channel)
                    self.assertIsInstance(rebuilt, IndependentChannelReconstruction)
                    self.assertEqual(rebuilt.channel, channel)
                    self.assertEqual(rebuilt.row_count, OWNED_ROWS)
                    self.assertEqual(len(rebuilt.raw_rows), OWNED_ROWS)
                    self.assertEqual(len(rebuilt.corrected_rows), OWNED_ROWS)
                    self.assertEqual(len(rebuilt.accumulation_bounds), 3)
                    self.assertEqual(len(rebuilt.embedded_defect_bounds), 3)
                    self.assertNotIn("classification", rebuilt.__dataclass_fields__)
                    self.assertNotIn("pass", rebuilt.__dataclass_fields__)
                    for row in rebuilt.raw_rows:
                        self.assertEqual(len(row), 3)
                        self.assertEqual(len(row[0]), 5)
                        self.assertEqual(len(row[1]), 2)
                        self.assertEqual(len(row[2]), 4)
                        for item in row[0]:
                            self.assertIsInstance(item, Fraction)

    def test_inputs_remain_immutable(self):
        state = _state(u=-0.0)
        rhs = _zero_rhs()
        before_state = _fingerprint_state(state)
        paths = _paths(propose_step, PRIMARY_METHOD, state, rhs)
        before_paths = _fingerprint_paths(paths)
        family = _validate(paths, PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID)
        reconstruct_channel_independently(family, "u:alpha")
        self.assertEqual(_fingerprint_state(state), before_state)
        self.assertEqual(_fingerprint_paths(paths), before_paths)
        self.assertIsInstance(family.paths, tuple)
        with self.assertRaises(FrozenInstanceError):
            family.method = COMPARATOR_METHOD  # type: ignore[misc]
        with self.assertRaises(TypeError):
            family.paths[0] = ()  # type: ignore[index]

    def test_signed_zero_owned_stage_tampering_is_rejected(self):
        paths = _paths(propose_step, PRIMARY_METHOD, _state(u=0.0), _zero_rhs())
        stage = paths[0][0].stages[1]
        tampered = np.array(stage.state.u, copy=True)
        tampered[1, 0] = -0.0
        if tampered.tobytes() == stage.state.u.tobytes():
            self.skipTest("owned byte already encodes minus zero")
        broken = _replace_stage(
            paths[0][0],
            1,
            state=EvolutionState(tampered, stage.state.p, stage.state.q),
        )
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                ((broken,), paths[1], paths[2]), PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID
            )

    def test_manual_hashed_constructor_cannot_skip_arithmetic_validation(self):
        paths, family = _family(
            propose_step,
            ORIGINAL_ARITHMETIC_ID,
            PRIMARY_METHOD,
            _state(),
            _unit_u_rhs(),
        )
        moved = _mutate_owned_state(paths[0][0].candidate_state)
        altered = replace(
            paths[0][0],
            candidate_state=moved,
            stages=(
                *paths[0][0].stages[:-1],
                replace(paths[0][0].stages[-1], state=moved),
            ),
        )
        forged = ((altered,), paths[1], paths[2])
        digest = _independent_family_digest(
            forged,
            family.method,
            family.arithmetic_id,
            family.owned_row_count,
        )
        with self.assertRaises(IndependentRecordedFamilyError):
            IndependentRecordedFamily(
                paths=forged,
                method=family.method,
                arithmetic_id=family.arithmetic_id,
                owned_row_count=family.owned_row_count,
                family_sha256=digest,
            )

    def test_public_constructor_rejects_hash_relabel_and_wrong_digest(self):
        paths, family = _family(
            propose_step,
            ORIGINAL_ARITHMETIC_ID,
            PRIMARY_METHOD,
            _state(),
            _unit_u_rhs(),
        )
        with self.assertRaises(IndependentRecordedFamilyError):
            replace(family, arithmetic_id=COMPENSATED_ARITHMETIC_ID)
        with self.assertRaises(IndependentRecordedFamilyError):
            replace(family, method=COMPARATOR_METHOD)
        with self.assertRaises(IndependentRecordedFamilyError):
            replace(family, family_sha256="a" * 64)
        with self.assertRaises(IndependentRecordedFamilyError):
            IndependentRecordedFamily(
                paths=paths,
                method=family.method,
                arithmetic_id=family.arithmetic_id,
                owned_row_count=family.owned_row_count,
                family_sha256="A" * 64,
            )

    def test_postvalidation_numpy_mutation_is_rejected(self):
        _, family = _family(
            propose_step,
            ORIGINAL_ARITHMETIC_ID,
            PRIMARY_METHOD,
            _state(),
            _unit_u_rhs(),
        )
        array = family.paths[0][0].candidate_state.u
        array.setflags(write=True)
        array[1, 0] += 0.5
        array.setflags(write=False)
        with self.assertRaisesRegex(
            IndependentRecordedFamilyError, "changed after validation"
        ):
            reconstruct_channel_independently(family, "u:alpha")


class RejectionTests(unittest.TestCase):
    def setUp(self):
        self.state = _state()
        self.rhs = _unit_u_rhs()
        self.paths = _paths(propose_step, PRIMARY_METHOD, self.state, self.rhs)

    def test_bad_field_count_name_time_method_arithmetic_and_row_count(self):
        five = _state(fields=5)
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                _paths(propose_step, PRIMARY_METHOD, five, _unit_u_rhs()),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(self.paths, PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID, rows=3)
        with self.assertRaises(TypeError):
            _validate(self.paths, PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID, rows=True)
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(self.paths, "RK4", ORIGINAL_ARITHMETIC_ID)
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(self.paths, PRIMARY_METHOD, "not-an-arithmetic")
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(self.paths, COMPARATOR_METHOD, ORIGINAL_ARITHMETIC_ID)
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                (self.paths[0], self.paths[1], self.paths[2][:3]),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )
        renamed = _replace_stage(self.paths[0][0], 1, stage_name="rk4_kX")
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                ((renamed,), self.paths[1], self.paths[2]),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )
        mistimed = _replace_stage(self.paths[0][0], 1, time=0.2)
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                ((mistimed,), self.paths[1], self.paths[2]),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )
        shortened = replace(self.paths[0][0], stages=self.paths[0][0].stages[:-1])
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                ((shortened,), self.paths[1], self.paths[2]),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )
        for non_binary64 in (False, 0, "0.0", Q(0), np.float32(0.0)):
            malformed = _replace_stage(self.paths[0][0], 0, time=non_binary64)
            with self.assertRaisesRegex(IndependentRecordedFamilyError, "binary64"):
                _validate(
                    ((malformed,), self.paths[1], self.paths[2]),
                    PRIMARY_METHOD,
                    ORIGINAL_ARITHMETIC_ID,
                )

    def test_state_rhs_continuity_layout_and_owned_replay_tampering(self):
        outer, medium, fine = self.paths
        mutated_state = _mutate_owned_state(medium[1].initial_state)
        first = replace(medium[1].stages[0], state=mutated_state)
        broken_state = replace(
            medium[1],
            initial_state=mutated_state,
            stages=(first, *medium[1].stages[1:]),
        )
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "contiguous"):
            _validate(
                (outer, (medium[0], broken_state), fine),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )

        mutated_rhs = _mutate_owned_rhs(medium[1].stages[0].rhs)
        broken_rhs = _replace_stage(medium[1], 0, rhs=mutated_rhs)
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "RHS"):
            _validate(
                (outer, (medium[0], broken_rhs), fine),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )

        gapped = replace(medium[1], initial_time=medium[1].initial_time + 0.01)
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "contiguous"):
            _validate(
                (outer, (medium[0], gapped), fine),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )

        stage = outer[0].stages[2]
        tampered_stage = _replace_stage(
            outer[0],
            2,
            state=_mutate_owned_state(stage.state),
        )
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "owned"):
            _validate(
                ((tampered_stage,), medium, fine),
                PRIMARY_METHOD,
                ORIGINAL_ARITHMETIC_ID,
            )

        endpoint_only = replace(
            outer[0], candidate_state=_mutate_owned_state(outer[0].candidate_state)
        )
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "endpoint"):
            _validate(
                ((endpoint_only,), medium, fine), PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID
            )

        moved = _mutate_owned_state(outer[0].candidate_state)
        last = replace(outer[0].stages[-1], state=moved)
        perturbed = replace(
            outer[0],
            candidate_state=moved,
            stages=(*outer[0].stages[:-1], last),
        )
        with self.assertRaises(IndependentRecordedFamilyError):
            _validate(
                ((perturbed,), medium, fine), PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID
            )

        layout_paths = _paths(propose_step, PRIMARY_METHOD, self.state, self.rhs)
        fortran = np.asfortranarray(
            np.array(layout_paths[0][0].initial_state.u, copy=True)
        )
        self.assertFalse(fortran.flags.c_contiguous)
        object.__setattr__(layout_paths[0][0].initial_state, "u", fortran)
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "binary64"):
            _validate(layout_paths, PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID)

        dtype_paths = _paths(propose_step, PRIMARY_METHOD, self.state, self.rhs)
        narrow = np.array(dtype_paths[0][0].initial_state.u, dtype=np.float32)
        object.__setattr__(dtype_paths[0][0].initial_state, "u", narrow)
        with self.assertRaisesRegex(IndependentRecordedFamilyError, "binary64"):
            _validate(dtype_paths, PRIMARY_METHOD, ORIGINAL_ARITHMETIC_ID)

    def test_unknown_channel_and_wrong_family_type_fail(self):
        _, family = _family(
            propose_step, ORIGINAL_ARITHMETIC_ID, PRIMARY_METHOD, self.state, self.rhs
        )
        with self.assertRaises(IndependentRecordedFamilyError):
            reconstruct_channel_independently(family, "u:shift")
        with self.assertRaises(TypeError):
            reconstruct_channel_independently(
                object(),
                "u:alpha",  # type: ignore[arg-type]
            )


class ChannelMappingTests(unittest.TestCase):
    def test_owned_row_and_eighteen_channel_mapping(self):
        u = np.zeros((POINT_COUNT, FIELD_COUNT), dtype=np.float64)
        p = np.zeros_like(u)
        q = np.zeros_like(u)
        du = np.zeros_like(u)
        dp = np.zeros_like(u)
        dq = np.zeros_like(u)
        u[0] = 999.0
        p[0] = 999.0
        q[0] = 999.0
        u[-4:] = 888.0
        p[-4:] = 888.0
        q[-4:] = 888.0
        for row in range(1, 5):
            for field in range(FIELD_COUNT):
                u[row, field] = 10 * row + field + 0.125
                p[row, field] = 100 * row + field + 0.25
                q[row, field] = 1000 * row + field + 0.5
                du[row, field] = 0.01 * (row + 1) * (field + 1)
                dp[row, field] = 0.02 * (row + 1) * (field + 1)
                dq[row, field] = 0.03 * (row + 1) * (field + 1)

        def rhs(_time, state):
            return EvolutionRHS(
                np.array(du, copy=True),
                np.array(dp, copy=True),
                np.array(dq, copy=True),
                {},
            )

        state = EvolutionState(u, p, q)
        before = _fingerprint_state(state)
        for proposer, arithmetic_id in _PROPOSERS:
            for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
                paths, family = _family(proposer, arithmetic_id, method, state, rhs)
                before_paths = _fingerprint_paths(paths)
                for channel in TDG6_COMPLETE_STATE_CHANNELS:
                    block, field_name = channel.split(":")
                    field = FIELDS.index(field_name)
                    rebuilt = reconstruct_channel_independently(family, channel)
                    source = {"u": u, "p": p, "q": q}[block]
                    slope = {"u": du, "p": dp, "q": dq}[block]
                    for row, raw in enumerate(rebuilt.raw_rows):
                        grid_row = row + 1
                        self.assertEqual(
                            raw[0][0], Q.from_float(float(source[grid_row, field]))
                        )
                        self.assertEqual(
                            raw[0][1], Q.from_float(float(slope[grid_row, field]))
                        )
                        self.assertNotEqual(raw[0][0], Q.from_float(999.0))
                        self.assertNotEqual(raw[0][0], Q.from_float(888.0))
                self.assertEqual(_fingerprint_state(state), before)
                self.assertEqual(_fingerprint_paths(paths), before_paths)


class AlgebraTests(unittest.TestCase):
    def test_delta_telescoping_and_noncancelling_bound(self):
        def signed_rhs(time, state):
            value = 1.0 if time < 0.125 else -1.0
            return EvolutionRHS(
                np.full_like(state.u, value),
                np.zeros_like(state.p),
                np.zeros_like(state.q),
                {},
            )

        state = _state(u=LARGE_BASE)
        paths, family = _family(
            propose_step, ORIGINAL_ARITHMETIC_ID, PRIMARY_METHOD, state, signed_rhs
        )
        before_paths = _fingerprint_paths(paths)
        rebuilt = reconstruct_channel_independently(family, "u:alpha")
        mixed = False
        for row, (raw, corrected) in enumerate(
            zip(rebuilt.raw_rows, rebuilt.corrected_rows, strict=True)
        ):
            for level, raw_segments, corr_segments, bound in (
                (0, (raw[0],), (corrected[0],), rebuilt.accumulation_bounds[0]),
                (1, raw[1], corrected[1], rebuilt.accumulation_bounds[1]),
                (2, raw[2], corrected[2], rebuilt.accumulation_bounds[2]),
            ):
                deltas = _path_deltas(raw_segments, corr_segments)
                y = [raw_segments[0][0]] + [item[2] for item in raw_segments]
                z = [corr_segments[0][0]] + [item[2] for item in corr_segments]
                self.assertEqual(y[0], z[0])
                total = Q(0)
                signed = Q(0)
                for index, delta in enumerate(deltas, start=1):
                    total += abs(delta)
                    signed += delta
                    self.assertEqual(y[index] - z[index], signed)
                    self.assertEqual(
                        corr_segments[index - 1][1], raw_segments[index - 1][1]
                    )
                    self.assertEqual(
                        corr_segments[index - 1][3], raw_segments[index - 1][3]
                    )
                self.assertLessEqual(abs(sum(deltas)), bound)
                self.assertEqual(total, sum(abs(item) for item in deltas))
                if any(item > 0 for item in deltas) and any(
                    item < 0 for item in deltas
                ):
                    mixed = True
                    self.assertGreater(
                        sum(abs(item) for item in deltas), abs(sum(deltas))
                    )
                    self.assertGreater(bound, abs(sum(deltas)))
                if row == 0:
                    _ = level
        self.assertTrue(mixed)
        fine_abs = [
            sum(abs(item) for item in _path_deltas(raw[2], corr[2]))
            for raw, corr in zip(rebuilt.raw_rows, rebuilt.corrected_rows, strict=True)
        ]
        self.assertEqual(rebuilt.accumulation_bounds[2], max(fine_abs))
        self.assertEqual(len(rebuilt.raw_rows[0][2]), 4)
        self.assertEqual(_fingerprint_paths(paths), before_paths)

    def test_accumulation_bound_is_row_maximum_of_absolute_sums(self):
        u = np.ones((POINT_COUNT, FIELD_COUNT), dtype=np.float64)
        u[2] = LARGE_BASE
        state = EvolutionState(
            u,
            np.zeros_like(u),
            np.zeros_like(u),
        )
        _, family = _family(
            propose_step,
            ORIGINAL_ARITHMETIC_ID,
            PRIMARY_METHOD,
            state,
            _unit_u_rhs(),
        )
        rebuilt = reconstruct_channel_independently(family, "u:alpha")
        per_row = [
            sum(abs(item) for item in _path_deltas((raw[0],), (corr[0],)))
            for raw, corr in zip(rebuilt.raw_rows, rebuilt.corrected_rows, strict=True)
        ]
        self.assertEqual(rebuilt.accumulation_bounds[0], max(per_row))
        self.assertEqual(per_row[0], 0)
        self.assertEqual(per_row[1], Q(1, 4))
        self.assertGreater(max(per_row), per_row[0])

    def test_embedded_linear_and_tableau_controls(self):
        self.assertEqual(sum(_RK4_WEIGHTS), Q(1))
        self.assertEqual(sum(_SSP_WEIGHTS), Q(1))
        rk4_c = (Q(0), Q(1, 2), Q(1, 2), Q(1))
        self.assertEqual(
            sum(b * c for b, c in zip(_RK4_WEIGHTS, rk4_c, strict=True)),
            Q(1, 2),
        )
        ssp_c = (Q(0), Q(1), Q(1, 2))
        self.assertEqual(
            sum(b * c for b, c in zip(_SSP_WEIGHTS, ssp_c, strict=True)),
            Q(1, 2),
        )
        state = _state(u=1.0, p=2.0, q=4.0)
        for proposer, arithmetic_id in _PROPOSERS:
            for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
                _, constant = _family(
                    proposer,
                    arithmetic_id,
                    method,
                    state,
                    _constant_rhs(du=1.0, dp=1.0, dq=1.0),
                )
                rebuilt = reconstruct_channel_independently(constant, "u:alpha")
                self.assertEqual(rebuilt.embedded_defect_bounds, (0, 0, 0))
                _, linear = _family(
                    proposer, arithmetic_id, method, state, _linear_rhs()
                )
                observed = reconstruct_channel_independently(linear, "u:alpha")
                expected = []
                weights = _RK4_WEIGHTS if method == PRIMARY_METHOD else _SSP_WEIGHTS
                for path in linear.paths:
                    total = Q(0)
                    for proposal in path:
                        h = Q.from_float(proposal.final_time - proposal.initial_time)
                        ks = [
                            Q.from_float(float(stage.rhs.du[OWNED][0, 0]))
                            for stage in proposal.stages[: len(weights)]
                        ]
                        end = Q.from_float(
                            float(proposal.stages[-1].rhs.du[OWNED][0, 0])
                        )
                        if method == PRIMARY_METHOD:
                            total += abs(h * (ks[3] - end) / 6)
                        else:
                            total += abs(h * (-ks[0] - ks[1] + 2 * ks[2]) / 3)
                    expected.append(total)
                self.assertEqual(observed.embedded_defect_bounds, tuple(expected))


class LargeBaseControlTests(unittest.TestCase):
    def test_reconstruction_agreement_does_not_erase_actual_error(self):
        self.assertEqual(LARGE_BASE + MACRO_WIDTH, LARGE_BASE)
        state = _state(u=LARGE_BASE, p=0.0, q=0.0)
        before = _fingerprint_state(state)
        for proposer, arithmetic_id in _PROPOSERS:
            for method in (PRIMARY_METHOD, COMPARATOR_METHOD):
                paths, family = _family(
                    proposer, arithmetic_id, method, state, _unit_u_rhs()
                )
                self.assertEqual(_fingerprint_state(state), before)
                for proposal in paths[0] + paths[1] + paths[2]:
                    np.testing.assert_array_equal(proposal.candidate_state.u, state.u)
                rebuilt = reconstruct_channel_independently(family, "u:alpha")
                base = Q.from_float(LARGE_BASE)
                for row in rebuilt.corrected_rows:
                    self.assertEqual(
                        row[0], (base, Q(1), base + Q(1, 4), Q(1), Q(1, 4))
                    )
                    self.assertEqual(
                        row[1],
                        (
                            (base, Q(1), base + Q(1, 8), Q(1), Q(1, 8)),
                            (base + Q(1, 8), Q(1), base + Q(1, 4), Q(1), Q(1, 8)),
                        ),
                    )
                    self.assertEqual(
                        [item[0] for item in row[2]] + [row[2][-1][2]],
                        [base + Q(index, 16) for index in range(5)],
                    )
                    self.assertEqual(len(row[2]), 4)
                self.assertEqual(
                    rebuilt.accumulation_bounds, (Q(1, 4), Q(1, 4), Q(1, 4))
                )
                self.assertEqual(rebuilt.accumulation_bounds[2], Q(1, 4))
                self.assertEqual(rebuilt.raw_rows[0][0][2], base)
                self.assertEqual(rebuilt.corrected_rows[0][0][2], base + Q(1, 4))
                self.assertEqual(rebuilt.corrected_rows[0][2][0][2], base + Q(1, 16))
                self.assertNotEqual(
                    rebuilt.corrected_rows[0][2][0][2], rebuilt.raw_rows[0][2][0][2]
                )
                rest = reconstruct_channel_independently(family, "p:alpha")
                self.assertEqual(rest.accumulation_bounds, (0, 0, 0))
                self.assertEqual(rest.embedded_defect_bounds, (0, 0, 0))
                self.assertEqual(rest.raw_rows[0][0][0], 0)
                self.assertEqual(rest.corrected_rows[0][0][2], 0)


if __name__ == "__main__":
    unittest.main()
