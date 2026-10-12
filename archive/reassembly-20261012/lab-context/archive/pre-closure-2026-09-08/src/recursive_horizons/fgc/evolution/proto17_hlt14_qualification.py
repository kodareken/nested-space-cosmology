"""Deterministic temporary-store qualification for HLT14 synthetic machinery."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from typing import Any, Callable

import numpy as np

from .proto17_hlt14_inputs import (
    HLT14AuthorityError, HLT14BundleError, HLT14InputError, HLT14NamespaceForbiddenError,
    LaunchAuthorityTuple, LaunchAuthorityVerifier, VerifiedBundleAdapter,
    canonical_bytes, canonical_sha256, shape_framed_array_content_sha256,
)
from .proto17_hlt14_runtime import Proto17HLT14Error, Proto17HLT14Runtime
from .proto17_pure_construction import (
    build_genesis, canonical, checkpoint, construct_common_event, cursor,
    digest, synthetic_genesis_fixture,
)


MEMBERS = ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385")
PHASES = (
    "before_write", "after_write", "before_flush", "after_flush", "before_fsync",
    "after_fsync", "before_replace", "after_replace", "before_dir_fsync", "after_dir_fsync",
)
SEALED_HASHES = {
    "generation_zero_checkpoint_sha256": "01e8214dcbbef025fe7883670043659d3ec66bb47e9b3e9dd4b8089b8ba5f8de",
    "common_event_receipt_sha256": "17f925e8489341392ac5df4fba1cc47b66a27153e64eca4455c3f701f37cbc17",
    "successor_checkpoint_sha256": "c4c31f53a6a97632addafd54d6b3ec08bdcffe4aa4321269f8c8bca957e5380b",
}


def _git(root: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=root, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=True).stdout.strip()


def _full_synthetic_genesis(root: Path) -> dict[str, object]:
    """Create one complete six-member real-NPZ GenesisSpec under ``root``."""
    root.mkdir(parents=True, exist_ok=False)
    fixture = deepcopy(synthetic_genesis_fixture()); fixture["physical_input_manifest"]["relative_input_root"] = "inputs"; bundles = root / "inputs" / "bundles"; bundles.mkdir(parents=True)
    for index, key in enumerate(MEMBERS):
        descriptor = fixture["member_descriptors"][index]; entry = fixture["physical_input_manifest"]["bundle_members"][index]
        state = descriptor["state_object"]; arrays: dict[str, np.ndarray] = {}
        for array_index, item in enumerate(state["physical_arrays"]):
            array = np.full(tuple(item["shape"]), index + array_index / 16.0, dtype="<f8")
            arrays[item["logical_name"]] = array
            item["little_endian_c_bytes_sha256"] = sha256(array.tobytes(order="C")).hexdigest()
        state["physical_state_sha256"] = shape_framed_array_content_sha256(arrays["u"], arrays["p"], arrays["q"])
        state["restart_payload_sha256"] = shape_framed_array_content_sha256(*(arrays[name] for name in ("u", "p", "q", "tracer_positions", "tracer_proper_times", "event_proper_times", "event_fields")))
        descriptor["state_object_canonical_sha256"] = digest(state)
        metadata = {name: descriptor[name] if name in descriptor else state[name] for name in (
            "member_key", "method", "point_count", "coordinate_time", "accepted_boundary_time", "step_index", "transaction_serial", "input_hash", "source_retry_count", "CFL_retry_count", "runtime_monitor_state", "causal_state", "tracer_state", "event_history_state", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256")}
        path = bundles / f"{key}.npz"; saved = {item["npz_storage_key"]: arrays[item["logical_name"]] for item in state["physical_arrays"]}
        np.savez(path, **saved, metadata_utf8=np.frombuffer(canonical_bytes(metadata), dtype=np.uint8))
        entry["relative_input_path"] = f"bundles/{key}.npz"; entry["raw_file_sha256"] = sha256(path.read_bytes()).hexdigest()
    package = root / "fixturepkg"; package.mkdir(); (package / "__init__.py").write_text("")
    roles = {"runtime_module": "runtime.py", "adapter_module": "adapter.py", "runner": "runner.py"}
    for filename in roles.values(): (package / filename).write_text("VALUE = 1\n")
    for role, filename in roles.items(): fixture[role + "_sha256"] = sha256((package / filename).read_bytes()).hexdigest()
    fixture["physical_input_manifest_sha256"] = digest(fixture["physical_input_manifest"])
    for field in ("member_descriptor_set_sha256", "genesis_checkpoint_sha256", "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256"):
        fixture[field] = "0" * 64
    return build_genesis(fixture, _validate_derived=False)["genesis_spec"]


def _authority_case(root: Path) -> tuple[dict[str, object], LaunchAuthorityTuple]:
    genesis = _full_synthetic_genesis(root)
    _git(root, "init", "-q"); _git(root, "config", "user.email", "hlt14@example.invalid"); _git(root, "config", "user.name", "HLT14")
    roles = {"runtime_module": "runtime.py", "adapter_module": "adapter.py", "runner": "runner.py"}
    pinned = [{"path": f"fixturepkg/{filename}", "sha256": genesis[role + "_sha256"], "module": f"fixturepkg.{filename[:-3]}", "role": role} for role, filename in roles.items()]
    result = {"artifact_id": "FGC-1-HLT14-MON14", "genesis_spec": genesis, "pinned_sources": pinned}
    authority = root / "authority.json"; authority.write_bytes(canonical_bytes(result)); _git(root, "add", "."); _git(root, "commit", "-qm", "synthetic authority"); commit = _git(root, "rev-parse", "HEAD")
    for module in ("fixturepkg", "fixturepkg.runtime", "fixturepkg.adapter", "fixturepkg.runner"): sys.modules.pop(module, None)
    sys.path.insert(0, str(root))
    try:
        launch = LaunchAuthorityTuple(commit, "authority.json", sha256(authority.read_bytes()).hexdigest(), canonical_sha256(genesis))
        verified = LaunchAuthorityVerifier(root).verify(launch)
        bundles = VerifiedBundleAdapter.from_manifest_root(root, verified.genesis_spec).verify_genesis(verified.genesis_spec)
        if set(bundles) != set(MEMBERS): raise AssertionError("all-of six-bundle admission differs")
        return dict(verified.genesis_spec), launch
    finally:
        sys.path.remove(str(root))
        for module in ("fixturepkg", "fixturepkg.runtime", "fixturepkg.adapter", "fixturepkg.runner"): sys.modules.pop(module, None)


def _runtime_nominal(root: Path, spec: dict[str, object]) -> Proto17HLT14Runtime:
    runtime = Proto17HLT14Runtime(spec)
    genesis = runtime.materialize_generation_zero(root)
    if genesis != runtime.genesis:
        raise AssertionError("adapter genesis differs from constructor")
    runtime.replay_six_semantic_accepted_records(root)
    successor = runtime.commit_common_event(root)
    if runtime.open(root) != successor or successor["committed_common_event_index"] != 24:
        raise AssertionError("common event did not complete")
    return runtime


def _sealed_hash_vector(root: Path) -> dict[str, str]:
    """Reproduce the three PREF25 hashes on its unchanged abstract fixture.

    The real-NPZ authority path necessarily has different input-dependent state
    addresses.  Keeping this as a separate vector prevents those two kinds of
    evidence from being conflated.
    """
    runtime = Proto17HLT14Runtime(synthetic_genesis_fixture())
    genesis = runtime.materialize_generation_zero(root)
    runtime.replay_six_semantic_accepted_records(root)
    successor = runtime.commit_common_event(root)
    receipt = runtime._records(root)[-1]
    observed = {
        "generation_zero_checkpoint_sha256": genesis["checkpoint_sha256"],
        "common_event_receipt_sha256": receipt["record_sha256"],
        "successor_checkpoint_sha256": successor["checkpoint_sha256"],
    }
    if observed != SEALED_HASHES:
        raise AssertionError("HLT14 runtime differs from the sealed PREF25 hash vector")
    return observed


def _lone_recovery(root: Path, spec: dict[str, object]) -> None:
    runtime = Proto17HLT14Runtime(spec); runtime.materialize_generation_zero(root); runtime.replay_six_semantic_accepted_records(root)
    def fail(phase: str, target: str) -> None:
        if target == "common_event_receipt" and phase == "after_dir_fsync":
            raise RuntimeError("synthetic interruption")
    try:
        runtime.commit_common_event(root, fault_hook=fail)
    except RuntimeError:
        recovered = runtime.recover(root)
        if recovered["committed_common_event_index"] != 24:
            raise AssertionError("lone receipt did not recover")
    else:
        raise AssertionError("receipt fault did not interrupt")


def _fault_cases(spec: dict[str, object]) -> list[str]:
    names: list[str] = []
    for target in ("common_event_receipt", "checkpoint"):
        for phase in PHASES:
            with TemporaryDirectory() as temporary:
                root = Path(temporary) / "store"; runtime = Proto17HLT14Runtime(spec)
                runtime.materialize_generation_zero(root)
                predecessor = runtime.replay_six_semantic_accepted_records(root)
                prior_records = runtime._records(root)
                expected_successor = construct_common_event(
                    predecessor, prior_journal_records=prior_records,
                )["checkpoint"]
                def fail(observed_phase: str, observed_target: str, *, phase=phase, target=target) -> None:
                    if observed_phase == phase and observed_target == target:
                        raise RuntimeError("synthetic fault")
                try:
                    runtime.commit_common_event(root, fault_hook=fail)
                except RuntimeError:
                    recovered = runtime.recover(root)
                    receipt_is_durable = (
                        target == "checkpoint"
                        or phase in {"after_replace", "before_dir_fsync", "after_dir_fsync"}
                    )
                    expected = expected_successor if receipt_is_durable else predecessor
                    expected_journals = 7 if receipt_is_durable else 6
                    expected_checkpoints = 3 if receipt_is_durable else 2
                    if (
                        recovered != expected
                        or runtime.open(root) != expected
                        or len(list((root / "journal").glob("*.journal"))) != expected_journals
                        or len(list((root / "checkpoints").glob("*.json"))) != expected_checkpoints
                    ):
                        raise AssertionError("fault recovery differs")
                else:
                    raise AssertionError("fault hook did not fire")
                disposition = "successor" if receipt_is_durable else "predecessor"
                case_target = "common_event_checkpoint" if target == "checkpoint" else target
                names.append(f"fault:{case_target}:{phase}:{disposition}")
    return names


def _accepted_prefix_fault_cases(spec: dict[str, object]) -> list[str]:
    """Classify every accepted-prefix publication window exactly.

    An accepted journal record without its complete six-member checkpoint is
    deliberately not recoverable.  Only a fault before journal replacement
    leaves an untouched genesis; only a fully replaced generation-one
    checkpoint makes the complete prefix openable.
    """
    names: list[str] = []
    for target in ("accepted_journal", "checkpoint"):
        for phase in PHASES:
            with TemporaryDirectory() as temporary:
                root = Path(temporary) / "store"
                runtime = Proto17HLT14Runtime(spec)
                genesis = runtime.materialize_generation_zero(root)

                def fail(observed_phase: str, observed_target: str, *, phase=phase, target=target) -> None:
                    if observed_phase == phase and observed_target == target:
                        raise RuntimeError("synthetic accepted-prefix fault")

                try:
                    runtime.replay_six_semantic_accepted_records(root, fault_hook=fail)
                except RuntimeError:
                    replacement_visible = phase in {
                        "after_replace", "before_dir_fsync", "after_dir_fsync",
                    }
                    if target == "accepted_journal" and not replacement_visible:
                        if runtime.open(root) != genesis or runtime.recover(root) != genesis:
                            raise AssertionError("pre-replace accepted fault did not preserve genesis")
                        disposition = "genesis"
                    elif target == "checkpoint" and replacement_visible:
                        opened = runtime.open(root)
                        if (
                            opened["campaign_generation"] != 1
                            or opened["committed_common_event_index"] != 23
                            or runtime.recover(root) != opened
                        ):
                            raise AssertionError("replaced accepted checkpoint did not open exactly")
                        disposition = "accepted_prefix"
                    else:
                        for operation in (runtime.open, runtime.recover):
                            try:
                                operation(root)
                            except Proto17HLT14Error:
                                pass
                            else:
                                raise AssertionError("partial accepted prefix was recovered")
                        disposition = "invalid_prefix"
                else:
                    raise AssertionError("accepted-prefix fault hook did not fire")
                case_target = "accepted_checkpoint" if target == "checkpoint" else target
                names.append(f"fault:{case_target}:{phase}:{disposition}")
    return names


def _authority_mutations(root: Path, launch: LaunchAuthorityTuple) -> list[str]:
    names: list[str] = []; verifier = LaunchAuthorityVerifier(root)
    try:
        verifier.verify(LaunchAuthorityTuple("f" * 40, launch.authorization_result_path, launch.authorization_result_sha256, launch.genesis_spec_sha256))
    except HLT14AuthorityError: names.append("mutation:authority_foreign_commit")
    else: raise AssertionError("foreign authority accepted")
    raw = json.loads((root / "authority.json").read_text()); raw["genesis_spec"].pop("member_descriptors")
    (root / "authority.json").write_bytes(canonical_bytes(raw)); _git(root, "add", "authority.json"); _git(root, "commit", "-qm", "malformed authority")
    malformed = LaunchAuthorityTuple(_git(root, "rev-parse", "HEAD"), "authority.json", sha256((root / "authority.json").read_bytes()).hexdigest(), canonical_sha256(raw["genesis_spec"]))
    try: verifier.verify(malformed)
    except HLT14AuthorityError: names.append("mutation:authority_malformed_spec")
    else: raise AssertionError("malformed authority accepted")
    (root / "fixturepkg/runtime.py").write_text("VALUE = 2\n")
    try: verifier.verify(launch)
    except HLT14AuthorityError: names.append("mutation:authority_pinned_source_drift")
    else: raise AssertionError("pinned source drift accepted")
    return names


def _raw_mutations(root: Path, spec: dict[str, object]) -> list[str]:
    names: list[str] = []; adapter = VerifiedBundleAdapter.from_manifest_root(root, spec); descriptor = deepcopy(spec["member_descriptors"][0]); entry = deepcopy(spec["physical_input_manifest"]["bundle_members"][0])
    raw = deepcopy(entry); raw["raw_file_sha256"] = "0" * 64
    try: adapter.verify(descriptor, raw)
    except HLT14BundleError: names.append("mutation:raw_file_hash")
    else: raise AssertionError("raw-file hash mutation was accepted")
    shaped = deepcopy(descriptor); shaped["state_object"]["physical_arrays"][0]["shape"] = [1, 1]
    try: adapter.verify(shaped, entry)
    except HLT14BundleError: names.append("mutation:array_shape")
    else: raise AssertionError("array-shape mutation was accepted")
    hashed = deepcopy(descriptor); hashed["state_object"]["physical_state_sha256"] = "0" * 64
    try: adapter.verify(hashed, entry)
    except HLT14BundleError: names.append("mutation:aggregate_hash")
    else: raise AssertionError("aggregate-hash mutation was accepted")
    manifest_root = root / spec["physical_input_manifest"]["relative_input_root"]
    linked = manifest_root / "bundles" / "linked.npz"; linked.symlink_to(manifest_root / entry["relative_input_path"])
    symlink = deepcopy(entry); symlink["relative_input_path"] = "bundles/linked.npz"
    try: adapter.verify(descriptor, symlink)
    except (HLT14InputError, HLT14NamespaceForbiddenError): names.append("mutation:bundle_symlink")
    else: raise AssertionError("symlinked raw bundle was accepted")
    source = manifest_root / entry["relative_input_path"]
    with np.load(source, allow_pickle=False) as archive:
        payload = {
            key: archive[key].copy()
            for key in archive.files
            if key != "metadata_utf8"
        }
        original_metadata = json.loads(
            archive["metadata_utf8"].tobytes().decode("utf-8")
        )

    def reject_metadata(filename: str, raw_metadata: bytes, case: str) -> None:
        path = manifest_root / "bundles" / filename
        np.savez(
            path, **payload,
            metadata_utf8=np.frombuffer(raw_metadata, dtype=np.uint8),
        )
        altered = deepcopy(entry)
        altered["relative_input_path"] = f"bundles/{filename}"
        altered["raw_file_sha256"] = sha256(path.read_bytes()).hexdigest()
        try:
            adapter.verify(descriptor, altered)
        except HLT14BundleError:
            names.append(case)
        else:
            raise AssertionError(f"{case} was accepted")

    changed = deepcopy(original_metadata)
    changed["step_index"] += 1
    reject_metadata(
        "metadata-semantic.npz", canonical_bytes(changed),
        "mutation:metadata_semantic_binding",
    )
    reject_metadata(
        "metadata-noncanonical.npz",
        json.dumps(original_metadata, sort_keys=True, indent=1).encode("utf-8"),
        "mutation:metadata_noncanonical_json",
    )
    reject_metadata(
        "metadata-duplicate.npz",
        b'{"member_key":"first","member_key":"second"}',
        "mutation:metadata_duplicate_json",
    )
    return names


def _fully_rehashed_semantic_forgeries(root: Path, spec: dict[str, object]) -> list[str]:
    names: list[str] = []
    # Forgery one: rehash a changed executed plan and every dependent journal/checkpoint digest.
    runtime = Proto17HLT14Runtime(spec); one = root / "forgery-one"; runtime.materialize_generation_zero(one); predecessor = runtime.replay_six_semantic_accepted_records(one); records = runtime._records(one); forged = deepcopy(records); forged[0]["payload"]["executed_plan"] = {"synthetic": False}; parent = "0" * 64
    journal = one / "journal"
    for number, record in enumerate(forged, 1):
        record["journal_parent_sha256"] = parent; bare = {key: value for key, value in record.items() if key != "record_sha256"}; record["record_sha256"] = digest(bare)
        old = records[number - 1]; (journal / f"{number:020d}-{old['record_sha256']}.journal").unlink(); encoded = canonical(record); (journal / f"{number:020d}-{record['record_sha256']}.journal").write_bytes(f"{len(encoded):08x} ".encode() + encoded + b"\n"); parent = record["record_sha256"]
    altered = deepcopy(predecessor); altered["journal_tip_sha256"] = parent; bare = {key: value for key, value in altered.items() if key != "checkpoint_sha256"}; altered["checkpoint_sha256"] = digest(bare); old_checkpoint = one / "checkpoints" / f"00000000000000000001-{predecessor['checkpoint_sha256']}.json"; old_checkpoint.unlink(); (one / "checkpoints" / f"00000000000000000001-{altered['checkpoint_sha256']}.json").write_bytes(canonical(altered))
    try: runtime.open(one)
    except Proto17HLT14Error: names.append("mutation:fully_rehashed_semantic_plan_forgery")
    else: raise AssertionError("fully rehashed plan forgery accepted")
    # Forgery two: rehash a high-water/cursor mutation into a canonical checkpoint.
    runtime = Proto17HLT14Runtime(spec); two = root / "forgery-two"; runtime.materialize_generation_zero(two); predecessor = runtime.replay_six_semantic_accepted_records(two); altered = deepcopy(predecessor); old = altered["cursors"]["RK4-2049"]; old["attempt_serial"] = 2; old["attempt_id"] = f"{altered['campaign_id']}:RK4-2049:23:2"; old.pop("cursor_chain_sha256"); altered["cursors"]["RK4-2049"] = cursor(old); altered["attempt_high_water"]["RK4-2049"] = 2; altered["cursor_generation_high_water"]["RK4-2049"] = 1
    forged_checkpoint = checkpoint(protocol_artifact_id=altered["protocol_artifact_id"], campaign_id=altered["campaign_id"], campaign_generation=1, committed_common_event_index=23, active_event_target_time=altered["active_event_target_time"], cursors=altered["cursors"], ledgers=altered["ledgers"], states=altered["states"], journal_tip_sha256=altered["journal_tip_sha256"], attempt_high_water=altered["attempt_high_water"], cursor_generation_high_water=altered["cursor_generation_high_water"], parent_checkpoint_sha256=altered["parent_checkpoint_sha256"])
    old_path = two / "checkpoints" / f"00000000000000000001-{predecessor['checkpoint_sha256']}.json"; old_path.unlink(); (two / "checkpoints" / f"00000000000000000001-{forged_checkpoint['checkpoint_sha256']}.json").write_bytes(canonical(forged_checkpoint))
    try: runtime.open(two)
    except Proto17HLT14Error: names.append("mutation:fully_rehashed_high_water_forgery")
    else: raise AssertionError("fully rehashed high-water forgery accepted")
    return names


def _mutation_cases(root: Path, spec: dict[str, object]) -> list[str]:
    names: list[str] = []
    # Extra accepted suffix after a committed event is invalid.
    runtime = Proto17HLT14Runtime(spec); store = root / "suffix"; runtime.materialize_generation_zero(store); runtime.replay_six_semantic_accepted_records(store); runtime.commit_common_event(store)
    last = sorted((store / "journal").glob("*.journal"))[-1]; (store / "journal" / "00000000000000000008-extra.journal").write_bytes(last.read_bytes())
    try: runtime.open(store)
    except (Proto17HLT14Error, ValueError): names.append("mutation:invalid_accepted_suffix")
    else: raise AssertionError("invalid suffix accepted")
    incomplete = root / "incomplete"; runtime.materialize_generation_zero(incomplete)
    try: runtime.commit_common_event(incomplete)
    except Proto17HLT14Error: names.append("mutation:incomplete_accepted_prefix")
    else: raise AssertionError("incomplete accepted prefix opened common event")
    existing = root / "preexisting"; existing.mkdir(); (existing / "user-file").write_text("preserve")
    try: runtime.materialize_generation_zero(existing)
    except Proto17HLT14Error: names.append("mutation:preexisting_root")
    else: raise AssertionError("preexisting root overwritten")
    # A store generated from a different injected authority cannot be adopted.
    foreign = root / "foreign"; other = deepcopy(spec); other["campaign_id"] = "foreign"; other = build_genesis(other, _validate_derived=False)["genesis_spec"]; foreign_runtime = Proto17HLT14Runtime(other); foreign_runtime.materialize_generation_zero(foreign)
    try: runtime.open(foreign)
    except Proto17HLT14Error: names.append("mutation:foreign_store")
    else: raise AssertionError("foreign store accepted")
    try: runtime.materialize_generation_zero(Path("runs/fgc-2-sf1/proto17/calibration"))
    except Proto17HLT14Error as error:
        if "production PROTO17 namespace" not in str(error):
            raise AssertionError("production namespace was rejected for the wrong reason") from error
        names.append("mutation:production_namespace_guard")
    else: raise AssertionError("production namespace accepted")
    semantic = root / "semantic"; runtime.materialize_generation_zero(semantic); runtime.replay_six_semantic_accepted_records(semantic)
    record = sorted((semantic / "journal").glob("*.journal"))[0]; record.write_bytes(record.read_bytes().replace(b'"synthetic":true', b'"synthetic":false'))
    try: runtime.open(semantic)
    except Proto17HLT14Error: names.append("mutation:semantic_six_record_replay")
    else: raise AssertionError("forged replay accepted")
    return names


def _generation_zero_staging_faults(spec: dict[str, object]) -> list[str]:
    names: list[str] = []
    for target, phase in (
        ("state_object", "before_write"),
        ("genesis_checkpoint", "after_fsync"),
        ("namespace", "before_namespace_rename"),
        ("namespace", "after_namespace_rename"),
    ):
        with TemporaryDirectory() as temporary:
            root = Path(temporary) / "staging"; runtime = Proto17HLT14Runtime(spec)
            def fail(observed_phase: str, observed_target: str, *, target=target, phase=phase) -> None:
                if observed_target == target and observed_phase == phase:
                    raise RuntimeError("generation zero staging fault")
            try: runtime.materialize_generation_zero(root, fault_hook=fail)
            except RuntimeError:
                if root.exists(): raise AssertionError("staging fault adopted partial root")
                names.append(f"fault:generation_zero:{target}:{phase}")
            else: raise AssertionError("generation zero fault hook did not fire")
    return names


def synthetic_qualification() -> dict[str, object]:
    """Execute every synthetic HLT14 integration category in temp directories."""
    names: list[str] = []
    with TemporaryDirectory() as temporary:
        base = Path(temporary)
        sealed_hashes = _sealed_hash_vector(base / "sealed-vector")
        names.append("nominal:sealed_PREF25_three_hash_vector")
        spec, launch = _authority_case(base / "authority"); names.append("nominal:temporary_git_authority_full_six_member_genesis")
        names.extend(f"bundle:{key}" for key in MEMBERS)
        _runtime_nominal(base / "nominal", spec); names.append("nominal:semantic_six_record_common_event")
        _lone_recovery(base / "recovery", spec); names.append("nominal:lone_receipt_recovery")
        names.extend(_authority_mutations(base / "authority", launch))
        names.extend(_raw_mutations(base / "authority", spec))
        names.extend(_mutation_cases(base / "mutations", spec))
        names.extend(_fully_rehashed_semantic_forgeries(base / "forgeries", spec))
    names.extend(_fault_cases(spec))
    names.extend(_accepted_prefix_fault_cases(spec))
    names.extend(_generation_zero_staging_faults(spec))
    ordered = sorted(names)
    return {
        "case_names": ordered,
        "case_count": len(ordered),
        "case_digest": sha256(canonical_bytes(ordered)).hexdigest(),
        "all_cases_passed": True,
        "sealed_PREF25_hash_vector": sealed_hashes,
        "sealed_PREF25_hash_vector_reproduced": sealed_hashes == SEALED_HASHES,
        "real_NPZ_authority_path_has_input_dependent_hashes": True,
        "temporary_directories_only": True,
        "in_process_exception_faults_not_power_loss_proof": True,
        "real_production_instance_requirements_completed": False,
    }


__all__ = ["PHASES", "MEMBERS", "SEALED_HASHES", "synthetic_qualification"]
